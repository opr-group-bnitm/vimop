// This file includes or is derived from code developed by Oxford Nanopore Technologies Plc.
// Copyright (c) 2020-, Oxford Nanopore Technologies Plc. All rights reserved.
// 
// This file is subject to the terms described in the LICENSE file in the root of this repository.


enum InputType {
    SingleFile,
    TopLevelDir,
    DirWithSubDirs,
}


/**
* Check if a file ends with one of the target extensions.
*
* @param file: path to the file in question
* @param extensions: list of valid file extensions
* @return: boolean whether the file has one of the provided extensions
*/
def is_target_file(file, extensions) {
    extensions.any { ext -> file.name.endsWith(ext) }
}


/**
 * Take a channel of the shape `[meta, reads, path-to-stats-dir | null]` and extract the
 * run IDs from the `run_ids` file in the stats directory into the metamap. If the path
 * to the stats dir is `null`, add an empty list.
 *
 * @param ch: input channel of shape `[meta, reads, path-to-stats-dir | null]`
 * @return: channel with a list of run IDs added to the metamap
 */
def add_run_IDs_to_meta(ch) {
    // HashSet for all observed run_ids
    Set<String> ingressed_run_ids = new HashSet<String>()

    // extract run_ids from fastcat stats / bamstats results and add to metadata as well
    // as `ingressed_run_ids`
    ch = ch.map { meta, reads, stats ->
        ArrayList run_ids = []
        if (stats) {
            run_ids = stats.resolve("run_ids").splitText().collect { line -> line.strip() }
            ingressed_run_ids += run_ids
        }
        // `meta + [...]` returns a new map which is handy to avoid any
        // modifying-maps-in-closures weirdness
        // See https://github.com/nextflow-io/nextflow/issues/2660
        [meta + [run_ids: run_ids], reads, stats]
    }
    // put run_ids somewhere global for trivial access later
    // bit grim but decouples ingress metadata from workflow main.nf
    // additionally no need to use CWUtil as we're not overriding any user params
    ch.subscribe(onComplete: {
        params.wf["ingress.run_ids"] = ingressed_run_ids
    })
    return ch
}


/**
 * Take a map of input arguments, find valid FASTQ inputs, and return a channel
 * with elements of `[metamap, seqs.fastq.gz | null, path-to-fastcat-stats | null]`.
 * The second item is `null` for sample sheet entries without a matching barcode
 * directory. The last item is `null` if `fastcat` was not run (it is only run on
 * directories containing more than one FASTQ file or when `stats: true`).
 *
 * @param arguments: map with arguments containing
 *  - "input": path to either: (i) input FASTQ file, (ii) top-level directory containing
 *     FASTQ files, (iii) directory containing sub-directories which contain FASTQ
 *     files
 *  - "sample": string to name single sample
 *  - "sample_sheet": path to CSV sample sheet
 *  - "analyse_unclassified": boolean whether to keep unclassified reads
 *  - "stats": boolean whether to write the `fastcat` stats
 *  - "fastcat_extra_args": string with extra arguments to pass to `fastcat`
 *  - "required_sample_types": list of required sample types in the sample sheet
 *  - "watch_path": boolean whether to use `watchPath` and run in streaming mode
 * @return: channel of `[Map(alias, barcode, type, ...), Path|null, Path|null]`.
 *  The first element is a map with metadata, the second is the path to the
 *  `.fastq.gz` file with the (potentially concatenated) sequences and the third is
 *  the path to the directory with the `fastcat` statistics. The second element is
 *  `null` for sample sheet entries for which no corresponding barcode directory was
 *  found. The third element is `null` if `fastcat` was not run.
 */
workflow fastq_ingress {
    take:
    arguments

    main:
    // check arguments
    Map margs = parse_arguments(arguments, ["fastcat_extra_args": ""])

    ArrayList fq_extensions = [".fastq", ".fastq.gz", ".fq", ".fq.gz"]

    // `watch_path` will be handled within `get_valid_inputs()`
    def input = get_valid_inputs(margs, fq_extensions)

    def ch_result
    if (margs.stats) {
        // run fastcat regardless of input type
        ch_result = fastcat(input.files.mix(input.dirs), margs["fastcat_extra_args"])
    } else {
        // run `fastcat` only on directories and rename / compress single files
        ch_result = fastcat(input.dirs, margs["fastcat_extra_args"])
        | mix(
            input.files
            | move_or_compress_fq_file
            | map { meta, path -> [meta, path, null] }
        )
    }
    // add sample sheet entries without barcode dirs to the results channel and extract
    // the run IDs into the metamaps before returning
    ch_result = ch_result.mix(input.missing.map { v -> v + [null] })
    def ch_samples = add_run_IDs_to_meta(ch_result)

    emit:
    samples = ch_samples
}


/**
 * Run `watchPath` on the input directory and return a channel of shape [metamap,
 * path-to-target-file]. The meta data is taken from the sample sheet in case one was
 * provided. Otherwise it only contains the `alias` (either `margs["sample"]` or the
 * name of the parent directory of the file).
 *
 * @param input: path to a directory to watch
 * @param margs: Map with parsed input arguments
 * @param extensions: list of valid extensions for the target file type
 * @return: Channel of [metamap, path-to-target-file]
 */
workflow watch_path {
    take:
    input
    margs
    extensions

    main:
    // we have two cases to consider: (i) files being generated in the top-level
    // directory and (ii) files being generated in sub-directories. If we find files of
    // both kinds, throw an error.
    if (input.isFile()) {
        error "Input ($input) must be a directory when using `watch_path`."
    }
    // get existing target files first (look for relevant files in the top-level dir and
    // all sub-dirs)
    ch_existing_input = channel.fromPath(input)
    | concat(channel.fromPath("$input/*", type: 'dir'))
    | map { dir -> get_target_files_in_dir(dir, extensions) }
    | flatten
    // now get channel with files found by `watchPath`
    ch_watched = channel.watchPath("$input/**").until { file_path -> file_path.name.startsWith('STOP') }
    // only keep target files
    | filter { file_path -> is_target_file(file_path, extensions) }
    // merge the channels
    ch_watched = ch_existing_input | concat(ch_watched)
    // check if input is as expected; start by throwing an error when finding files in
    // top-level dir and sub-directories
    String prev_input_type
    ch_watched
    | map { file_path ->
        String input_type = (file_path.parent == input) ? "top-level" : "sub-dir"
        if (prev_input_type && (input_type != prev_input_type)) {
            error "`watchPath` found input files in the top-level directory " +
                "as well as in sub-directories."
        }
        // if file is in a sub-dir, make sure it's not a sub-sub-dir
        if ((input_type == "sub-dir") && (file_path.parent.parent != input)) {
            error "`watchPath` found an input file more than one level of " +
                "sub-directories deep ('$file_path')."
        }
        // we also don't want files in the top-level dir when we got a sample sheet
        if ((input_type == "top-level") && margs["sample_sheet"]) {
            error "`watchPath` found input files in top-level directory even though " +
                "a sample sheet was provided ('${margs["sample_sheet"]}')."
        }
        prev_input_type = input_type
    }
    if (margs.sample_sheet) {
        // add metadata from sample sheet (we can't use join here since it does not work
        // with repeated keys; we therefore need to transform the sample sheet data into
        // a map with the barcodes as keys)
        def ch_sample_sheet = get_sample_sheet(file(margs.sample_sheet), margs.required_sample_types)
        | collect
        | map { rows -> rows.collectEntries { row -> [(row["barcode"]): row] } }
        // now we can use this channel to annotate all files with the corresponding info
        // from the sample sheet
        ch_watched = ch_watched
        | combine(ch_sample_sheet)
        | map { file_path, sample_sheet_map ->
            String barcode = file_path.parent.name
            Map sample_sheet_entry = sample_sheet_map[barcode]
            // throw error if the barcode was not in the sample sheet
            if (!sample_sheet_entry) {
                error "Sub-directory $barcode was not found in the sample sheet."
            }
            [create_metamap(sample_sheet_entry), file_path]
        }
    } else {
        ch_watched = ch_watched
        | map { file_path ->
            // This file could be in the top-level dir or a sub-dir. In the first case
            // check if a sample name was provided. In the second case, the alias is
            // always the name of the sub-dir.
            String alias
            if (file_path.parent == input) {
                // top-level dir
                alias = margs["sample"] ?: file_path.parent.name
            } else {
                // sub-dir
                alias = file_path.parent.name
            }
            [create_metamap([alias: alias]), file_path]
        }
    }

    emit:
    watched = ch_watched
}


process move_or_compress_fq_file {
    label "general"
    cpus 1
    memory "2 GB"
    input:
        // don't stage `input` with a literal because we check the file extension
        tuple val(meta), path(input)
    output:
        tuple val(meta), path("seqs.fastq.gz")
    script:
        String out = "seqs.fastq.gz"
        if (input.name.endsWith('.gz')) {
            // we need to take into account that the file could already be named
            // "seqs.fastq.gz" in which case `mv` would fail
            """
            [ "$input" == "$out" ] || mv "$input" $out
            """
        } else {
            """
            cat "$input" | bgzip -@ $task.cpus > $out
            """
        }
}


process fastcat {
    label "ingress"
    cpus 3
    memory "2 GB"
    input:
        tuple val(meta), path("input")
        val extra_args
    output:
        tuple val(meta),
              path("seqs.fastq.gz"),
              path("fastcat_stats")
    script:
        String out = "seqs.fastq.gz"
        String fastcat_stats_outdir = "fastcat_stats"
        """
        mkdir $fastcat_stats_outdir
        fastcat \
            -s ${meta["alias"]} \
            -r >(bgzip -c > $fastcat_stats_outdir/per-read-stats.tsv.gz) \
            -f $fastcat_stats_outdir/per-file-stats.tsv \
            --histograms histograms \
            $extra_args \
            input \
            | bgzip > $out

        mv histograms/* $fastcat_stats_outdir
        # extract the run IDs from the per-read stats
        csvtk cut -tf runid $fastcat_stats_outdir/per-read-stats.tsv.gz \
        | csvtk del-header | sort | uniq > $fastcat_stats_outdir/run_ids
        """
}


/**
 * Parse input arguments for `fastq_ingress` or `xam_ingress`.
 *
 * @param arguments: map with input arguments (see the corresponding ingress function
 *  for details)
 * @param extra_kwargs: map of extra keyword arguments and their defaults (this allows
 *  the argument-parsing to be tailored to a particular ingress function)
 * @return: map of parsed arguments
 */
def parse_arguments(arguments, extra_kwargs=[:]) {
    ArrayList required_args = ["input"]
    Map default_kwargs = [
        "sample": null,
        "sample_sheet": null,
        "analyse_unclassified": false,
        "stats": true,
        "required_sample_types": [],
        "watch_path": false
    ]
    ArgumentParser parser = new ArgumentParser(
        args: required_args,
        kwargs: default_kwargs + extra_kwargs,
        name: "fastq_ingress")
    return parser.parse_args(arguments)
}


/**
 * Find valid inputs based on the target extensions and return a branched channel with
 * branches `missing`, `files` and `dir`, which are of the shape `[metamap, input_path |
 * null]` (with `input_path` pointing to a target file or a directory containing target
 * files, respectively). `missing` contains sample sheet entries for which no
 * corresponding barcodes were found.
 * Unless `watchPath` was requested, the function checks whether the input is a single
 * target file, a top-level directory with target files, or a directory containing
 * sub-directories (usually barcodes) with target files.
 *
 * @param margs: parsed arguments (see `fastq_ingress` and `xam_ingress` for details)
 * @param extensions: list of valid extensions for the target file type
 * @return: branched channel with branches `missing`, `dir`, and `files`
*/
workflow get_valid_inputs {
    take:
    margs
    extensions

    main:
    log.info "Searching input for $extensions files."
    Path input
    try {
        input = file(margs.input, checkIfExists: true)
    } catch (java.nio.file.NoSuchFileException _e) {
        error "Input path $margs.input does not exist."
    }
    // declare resulting input channel
    def ch_input
    // run `watchPath` if requested
    if (margs["watch_path"]) {
        ch_input = watch_path(input, margs, extensions)
    } else {
        // check which of the allowed input types (single file, top-lvl dir, dir with
        // sub-dirs) we got
        input_type = determine_input_type(
            input, extensions, margs.analyse_unclassified
        )
        // handle case of `input` being a single file
        if (input_type == InputType.SingleFile) {
            ch_input = channel.of(
                [create_metamap([alias: margs["sample"] ?: input.simpleName]), input])
        } else if (input_type == InputType.TopLevelDir) {
            // input is a directory containing target files
            ch_input = channel.of(
                [create_metamap([alias: margs["sample"] ?: input.baseName]), input])
        } else {
            // input is a directory with sub-directories (e.g. barcodes) containing
            // target files --> find these sub-directories
            ArrayList sub_dirs_with_target_files = file(
                input.resolve('*'), type: "dir"
            ).findAll { sub_dir -> get_target_files_in_dir(sub_dir, extensions) }
            // remove directories called 'unclassified' unless otherwise specified
            if (!margs.analyse_unclassified) {
                sub_dirs_with_target_files = sub_dirs_with_target_files.findAll { sub_dir ->
                    sub_dir.baseName != "unclassified"
                }
            }
            // filter based on sample sheet in case one was provided
            if (margs.sample_sheet) {
                // get channel of entries in the sample sheet
                def ch_sample_sheet = get_sample_sheet(
                    file(margs.sample_sheet), margs.required_sample_types
                )
                // get the union of both channels (missing values will be replaced with
                // `null`)
                ch_union = channel.fromPath(sub_dirs_with_target_files).map { sub_dir ->
                    [sub_dir.baseName, sub_dir]
                }.join(ch_sample_sheet.map{ row -> [row.barcode, row]}, remainder: true)
                // after joining the channels, there are three possible cases:
                // (i) valid input path and sample sheet entry are both present
                // (ii) there is a sample sheet entry but no corresponding input dir
                //      --> we'll emit `[metamap-from-sample-sheet-entry, null]`
                // (iii) there is a valid path, but the sample sheet entry is missing
                //      --> drop this entry and print a warning to the log
                ch_input = ch_union.map {barcode, path, sample_sheet_entry ->
                    if (sample_sheet_entry) {
                        [create_metamap(sample_sheet_entry), path]
                    } else {
                        log.warn "Input directory '$barcode' was found, but sample " +
                            "sheet '$margs.sample_sheet' has no such entry."
                    }
                }
            } else {
                // no sample sheet --> simply emit the sub-dirs with the target files
                ch_input = channel.fromPath(sub_dirs_with_target_files).map { sub_dir ->
                    [create_metamap([alias: sub_dir.baseName, barcode: sub_dir.baseName]), sub_dir]
                }
            }
        }
    }
    // finally, we "unwrap" directories containing only a single file and then split the
    // results channel into the three different output types (sample sheet entries
    // without corresponding barcodes -- i.e. with `path == null`, single files, and
    // dirs with multiple files)
    def ch_branched_results = ch_input.map { meta, path ->
        if (path && path.isDirectory()) {
            List fq_files = get_target_files_in_dir(path, extensions)
            if (fq_files.size() == 1) {
                path = fq_files[0]
            }
        }
        [meta, path]
    } .branch { _meta, path ->
        missing: !path
        files: path.isFile()
        dirs: path.isDirectory()
    }

    emit:
    missing = ch_branched_results.missing
    files = ch_branched_results.files
    dirs = ch_branched_results.dirs
}

/**
 * Determine which of the allowed categories (single file, top-level directory, or
 * directory with sub-directory) an input path belongs to.
 *
 * @param margs: parsed arguments (see `fastq_ingress()` or `xam_ingress()` for details)
 * @param extensions: list of valid extensions for the target file type
 * @return: input type represented as an instance of the `InputType` enum
 */
def determine_input_type(
    input, extensions, analyse_unclassified
) {
    if (input.isFile()) {
        if (!is_target_file(input, extensions)) {
            error "Input file is not of required file type."
        }
        return InputType.SingleFile
    } else if (!input.isDirectory()){
        error "Input $input appears to be neither a file nor a directory."
    }
    // `input` is a directory --> we accept two cases: (i) a top-level directory with
    // target files and no sub-directories or (ii) a directory with one layer of
    // sub-directories containing target files. First, check if the directory contains
    // target files and find potential sub-directories (and sub-dirs with target files;
    // note that these lists can be empty)
    boolean dir_has_target_files = get_target_files_in_dir(input, extensions)
    ArrayList sub_dirs = file(input.resolve('*'), type: "dir")
    ArrayList sub_dirs_with_target_files = sub_dirs.findAll { sub_dir ->
        get_target_files_in_dir(sub_dir, extensions)
    }.findAll { sub_dir -> sub_dir.baseName != "unclassified" || analyse_unclassified }

    // define string to re-use in error messages below
    String target_files_str = \
        "target files (ending in ${extensions.collect{ ext -> '\'' + ext + '\''}.join(' / ')})"

    // check for target files in the top-level dir; if there are any, make sure there
    // are no sub-directories containing target files
    if (dir_has_target_files) {
        if (sub_dirs_with_target_files) {
            error "Input directory '$input' cannot contain $target_files_str " +
                "and also sub-directories with such files."
        }
        return InputType.TopLevelDir
    }

    // no target files in the top-level dir --> make sure there were sub-dirs with
    // target files
    if (!sub_dirs_with_target_files) {
        error "Input directory '$input' must contain either $target_files_str " +
            "or sub-directories containing such files (no more than one layer deep)."
    }
    // we don't allow sub-sub-directories with target files
    if (sub_dirs.any { sub_dir ->
        ArrayList subsubdirs = file(sub_dir.resolve('*'), type: "dir")
        subsubdirs.any { subsubdir -> get_target_files_in_dir(subsubdir, extensions) }
    }) {
        error "Input directory '$input' cannot contain more " +
            "than one level of sub-directories with $target_files_str."
    }
    return InputType.DirWithSubDirs
}


/**
 * Create a map that contains at least these keys: `[alias, barcode, type]`.
 * `alias` is required, `barcode` and `type` are filled with default values if
 * missing. Additional entries are allowed.
 *
 * @param arguments: map with input parameters; must contain `alias`
 * @return: map(alias, barcode, type, ...)
 */
def create_metamap(arguments) {
    ArgumentParser parser = new ArgumentParser(
        args: ["alias"],
        kwargs: [
            "barcode": null,
            "type": "test_sample",
            "run_ids": [],
        ],
        name: "create_metamap",
    )
    def metamap = parser.parse_known_args(arguments)
    metamap['alias'] = metamap['alias'].replaceAll(" ","_")
    return metamap
}


/**
 * Get the target files in the directory (non-recursive).
 *
 * @param dir: path to the target directory
 * @param extensions: list of valid extensions for the target file type
 * @return: list of found target files
 */
def get_target_files_in_dir(dir, extensions) {
    file(dir.resolve("*")).findAll { file_path -> is_target_file(file_path, extensions) }
}


/**
 * Check the sample sheet and return a channel with its rows if it is valid.
 *
 * @param sample_sheet: path to the sample sheet CSV
 * @return: channel of maps (with values in sample sheet header as keys)
 */
workflow get_sample_sheet {
    take:
    sample_sheet
    required_sample_types

    main:
    // If `validate_sample_sheet` does not return an error message, we can assume that
    // the sample sheet is valid and parse it. However, because of Nextflow's
    // asynchronous magic, we might emit values from `.splitCSV()` before the
    // error-checking closure finishes. This is no big deal, but undesired nonetheless
    // as the error message might be overwritten by the traces of new nextflow processes
    // in STDOUT. Thus, we use the somewhat clunky construct with `concat` and `last`
    // below. This lets the CSV channel only start to emit once the error checking is
    // done.
    ch_err = validate_sample_sheet(sample_sheet, required_sample_types).map { err_msg ->
        // check if there was an error message
        if (err_msg) error "Invalid sample sheet: ${err_msg}."
        err_msg
    }
    // concat the channel holding the path to the sample sheet to `ch_err` and call
    // `.last()` to make sure that the error-checking closure above executes before
    // emitting values from the CSV
    ch_rows = ch_err.concat(channel.fromPath(sample_sheet)).last().splitCsv(
        header: true, quote: '"'
    )

    emit:
    rows = ch_rows
}


/**
 * Python script for validating a sample sheet. The script will write messages
 * to STDOUT if the sample sheet is invalid. In case there are no issues, no
 * message is emitted.
 *
 * @param: path to sample sheet CSV
 * @param: list of required sample types (optional)
 * @return: string (optional)
 */
process validate_sample_sheet {
    label "general"
    cpus 1
    label "ingress"
    memory "2 GB"
    input:
        path "sample_sheet.csv"
        val required_sample_types
    output: stdout
    script:
    String req_types_arg = required_sample_types ? "--required_sample_types "+required_sample_types.join(" ") : ""
    """
    workflow-glue check_sample_sheet sample_sheet.csv $req_types_arg
    """
}
