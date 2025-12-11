// Copyright (c) 2025 Outbreak Preparedness and Response Group at BNITM
// This file is part of ViMOP and is licensed under the MIT License.
// See the LICENSE file in the root of this repository for full license details.

import groovy.json.JsonBuilder
import DatabaseInput
import SystemRequirements

nextflow.enable.dsl = 2


include { fastq_ingress } from './lib/ingress'

include {
    custom_data_base;
} from './lib/custom_db.nf'

include {
    db_update_get_config;
    update_data_base as update_virus;
    update_data_base as update_centrifuge;
    update_data_base as update_contaminants;
    data_base_transfer;
} from './lib/data_base_update.nf'

include {
    lengths_and_qualities as lengths_and_qualities_trimmed;
    lengths_and_qualities as lengths_and_qualities_cleaned;
    empty_tsv;
    empty_tsv as empty_kraken_style_report;
    empty_fasta;
    trim;
    read_stats;
    minCpus;
    minRAM;
    classify_centrifuge;
    classify_contigs;
    extract_contig_classification;
    no_contig_classification;
    filter_with_centrifuge;
    filter_contaminants;
    filter_virus_target;
    assemble_canu;
    reassemble_canu;
    no_assembly;
    pop_bubbles;
    prepare_blast_search;
    canu_contig_info;
    blast;
    extract_blasthits;
    get_ref_fasta;
    split_custom_ref;
    map_to_ref;
    map_to_sv_consensus;
    calc_coverage;
    subsample_alignments;
    sniffles;
    cutesv;
    structural_variant_consensus;
    medaka_variant_consensus;
    simple_consensus;
    auto_consensus;
    compute_mapping_stats;
    concat_mapping_stats;
    collect_reference_info;
    sample_report;
    get_best_consensus_files;
    simplify_reference_fasta;
    output;
} from './lib/processes'


workflow pipeline {
    take:
        samples
        db_config
    main:
        samplenames = samples
        | map { meta, reads, stats -> meta.alias }

        // trimming
        trimmed = samples
        | map { meta, reads, stats -> [meta, reads] }
        | trim
        | map { meta, reads -> [meta + ["trimmed_reads": reads], reads] }

        if (db_config.doClassify) {
            // metagenomic read classification with centrifuge
            classification = trimmed
            | map { meta, reads -> [meta, reads, db_config.classificationDir, db_config.classificationLibrary] }
            | classify_centrifuge

            kraken_style_reports = classification
            | map { meta, classification, report, kraken, html -> [meta.alias, kraken] }
        } else {
            classification = Channel.empty()

            kraken_style_reports = samplenames
            | empty_kraken_style_report
        }

        if(db_config.doFilterWithCentrifuge) {
            to_clean = classification
            | map { meta, classification, report, kraken, html -> [meta, meta.trimmed_reads, classification, db_config.virusTaxIDFile] }
            | filter_with_centrifuge
        } else {
            to_clean = trimmed
            | read_stats
        }

        // contaminant filtering
        cleaned = to_clean
        | map { meta, reads, stats -> [meta, reads, stats, db_config.contaminationFilterFiles, db_config.contaminationFilters] }
        | filter_contaminants

        // get readstats
        lenquals_trim = trimmed
        | map {meta, reads -> [meta.alias, reads]}
        | lengths_and_qualities_trimmed
    
        lenquals_clean = cleaned.reads
        | map {meta, reads -> [meta.alias, reads]}
        | lengths_and_qualities_cleaned

        // mapping reads to given virus targets to filter them
        mapped_to_virus_target = cleaned.reads
        | combine(Channel.from(db_config.virusTargets))
        | map{ meta, reads, target -> [meta + ["mapping_target": target.target], reads, target.path] }
        | filter_virus_target

        to_assemble_targeted = mapped_to_virus_target
        | map { meta, reads -> [meta, reads, params.nocontigs_enable_read_usage] }

        if (params.assemble_notarget) {
            to_assemble_notarget = cleaned.reads
            | map { meta, reads -> [meta + ["mapping_target": "no-target"], reads, false] }
        } else {
            to_assemble_notarget = Channel.empty()
        }

        if (params.assemble_notarget || db_config.virusTargets.size() > 0) {
            first_assemblies = to_assemble_notarget
            | mix(to_assemble_targeted)
            | assemble_canu
        } else {
            first_assemblies = samples
            | map { meta, reads, stats -> meta + ["mapping_target": "no-assembly"] }
            | no_assembly
        }

        // Re-assemble
        // try to assemble reads, that do not map to any previously created contig
        cleaned_reads = cleaned.reads
        | map { meta, reads -> [meta.alias, reads] }

        re_assemblies = first_assemblies.contigs
        | map { meta, contigs -> [meta.alias, contigs] }
        | groupTuple(by: 0)
        | join(cleaned_reads, by: 0)
        | map { samplename, contigs, reads -> [["alias": samplename, "mapping_target": "re-assembly"], contigs, reads] }
        | reassemble_canu

        contigs = first_assemblies.contigs
        | mix(re_assemblies.contigs)

        if (db_config.doClassify) {
            contig_classification = contigs
            | map { meta, contigs -> [meta, contigs, db_config.classificationDir, db_config.classificationLibrary] }
            | classify_contigs
            | extract_contig_classification
        } else {
            contig_classification = contigs
            | map {meta, contigs -> meta}
            | no_contig_classification
        }

        collected_contig_class_info = contig_classification
        | map { meta, contig_classification -> [meta.alias, contig_classification] }
        | groupTuple(by: 0)

        assembly_stats = first_assemblies.stats
        | mix(re_assemblies.stats)

        // database search for references using blast
        blast_queries = contigs
        | pop_bubbles
        | prepare_blast_search

        collected_contigs_infos = blast_queries
        | canu_contig_info
        | map { meta, contig_info -> [meta.alias, contig_info]}
        | groupTuple(by: 0)

        blast_hits = blast_queries
        | map { meta, contigs -> [meta, contigs, db_config.blastDir, db_config.blastPrefix] }
        | blast
        | extract_blasthits

        // Get mapping targets from BLAST hits stored in CSV files
        sample_ref = blast_hits
        | flatMap { meta, hits -> 
            hits.readLines()
                .drop(1)
                .collect { line -> tuple(meta.alias, meta.mapping_target, line.split(',')[1]) } 
        }
        | map { samplename, target_filter, blast_hit -> [samplename, blast_hit] }
        | unique

        // Get BLAST-derived reference sequences
        ref_seqs = sample_ref
        | map { samplename, ref_id -> ref_id }
        | unique
        | map { ref_id -> [ref_id, db_config.blastDir, db_config.blastPrefix] }
        | get_ref_fasta

        // Split custom_ref_fasta if provided
        custom_ref_fasta = params.custom_ref_fasta ?
                            Channel.fromPath(params.custom_ref_fasta) :
                            Channel.empty()

        custom_ref_files = custom_ref_fasta
        | split_custom_ref
        | flatten
        | map { fasta ->
            def ref_id = fasta.baseName
            [ref_id, fasta]
        }

        // Extract sample names
        samplenames = samples
        | map{ meta, reads, stats -> meta.alias }

        // Extract custom_refs
        custom_refs = custom_ref_files
        | map { ref_id, fasta -> ref_id }

        // combine each sample with every custom reference
        custom_sample_refs = samplenames
        | combine(custom_refs)

        extended_sample_ref = sample_ref
        | mix(custom_sample_refs)

        extended_ref_seqs = ref_seqs
        | mix(custom_ref_files)

        // Match reads with reference sequences for consensus generation
        simplified_ref_seqs = extended_ref_seqs
        | simplify_reference_fasta

        reads_and_ref = extended_sample_ref
        | combine(simplified_ref_seqs)
        | filter { samplename, ref_id_1, ref_id_2, ref_seq -> ref_id_1 == ref_id_2 }
        | map { samplename, ref_id_1, ref_id_2, ref_seq -> [samplename, ref_id_1, ref_seq] }
        | combine(trimmed)
        | filter { samplename, ref_id, ref_seq, meta, reads -> samplename == meta.alias }
        | map { samplename, ref_id, ref_seq, meta, reads -> [meta + ["consensus_target": ref_id], reads, ref_seq] }

        // Map reads against references for final consensus generation
        mapped_to_ref = reads_and_ref
        | map_to_ref
        | map { meta, ref, bam, bai -> [meta + ["ref": ref, "bam_to_ref": bam, "bai_to_ref": bai], ref, bam, bai]}

        coverage = mapped_to_ref
        | map { meta, ref, bam, bai -> [meta, bam, bai]}
        | calc_coverage

        if(params.sv_do_call_structural_variants) {

            subsampled_mapped_to_ref = mapped_to_ref
            | subsample_alignments

            if(params.sv_method == "cutesv"){
                structural_variants = subsampled_mapped_to_ref
                | cutesv
            } else if (params.sv_method == "sniffles") {
                structural_variants = subsampled_mapped_to_ref
                | sniffles
            } else {
                error "${params.sv_method} is not a valid choice for params.sv_method"
            }

            sv_consensus = structural_variants
            | map { meta, sv -> [meta + ["structural_variants": sv], meta.ref, sv] }
            | structural_variant_consensus

            mapped_to_sv_consensus = sv_consensus
            | map { meta, sv_consensus -> [meta, meta.trimmed_reads, meta.structural_variants, sv_consensus, meta.bam_to_ref, meta.bai_to_ref] }
            | map_to_sv_consensus
        } else {
            structural_variants = Channel.empty()

            mapped_to_sv_consensus = mapped_to_ref
        }

        if (params.consensus_method == 'medaka') {
            consensi = mapped_to_sv_consensus
            | map { meta, ref, bam, bai -> [
                meta, ref, bam, bai,
                meta.trimmed_reads] }
            | medaka_variant_consensus
        } else if (params.consensus_method == 'simple') {
            consensi = mapped_to_sv_consensus
            | simple_consensus
        } else if (params.consensus_method == 'auto') {
            consensi = mapped_to_sv_consensus
            | map { meta, ref, bam, bai -> [
                meta, ref, bam, bai,
                meta.trimmed_reads] }
            | auto_consensus
        }

        reference_info = empty_fasta
        | mix(extended_ref_seqs | map {refid, refseq -> refseq})
        | collect
        | collect_reference_info

        // compute the stats for the reads mapped to the different targets and build a big table.
        mapping_stats = mapped_to_ref
        | map {meta, ref, bam, bai -> [meta.alias, meta.consensus_target, meta, ref, bam, bai]}
        | join(consensi | map {meta, cons -> [meta.alias, meta.consensus_target, cons]}, by: [0, 1])
        | join(coverage | map {meta, cov -> [meta.alias, meta.consensus_target, cov]}, by: [0, 1])
        | map {samplename, target_name, meta, ref, bam, bai, cons, cov -> [meta, ref, bam, bai, cons, cov]}
        | compute_mapping_stats
        | map {meta, stats -> [meta.alias, stats]}

        // mix in empty files for cases, where no target was detected 
        // in order to still create a report
        collected_mapping_stats = samplenames
        | empty_tsv
        | mix(mapping_stats)
        | groupTuple(by: 0)
        | concat_mapping_stats

        // Create the report
        collected_assembly_stats = assembly_stats
        | map {meta, stats -> [meta.alias, stats]}
        | groupTuple(by: 0)

        collected_blast_hits = blast_hits
        | map {meta, hits -> [meta.alias, hits]}
        | groupTuple(by: 0)

        sample_results = cleaned.stats
        | join(collected_assembly_stats, by: 0)
        | join(collected_contig_class_info, by: 0)
        | join(collected_blast_hits, by: 0)
        | join(collected_mapping_stats, by: 0)
        | join(collected_contigs_infos, by: 0)
        | join(lenquals_trim, by: 0)
        | join(lenquals_clean, by: 0)
        | join(kraken_style_reports, by:0)
        | combine(Channel.of(db_config.virusConfigFileName))
        | combine(Channel.of(db_config.contaminationConfigFileName))
        | combine(Channel.of(db_config.classificationConfigFileName))
        | combine(reference_info)
        | sample_report

        collected_consensi = consensi
        | map {meta, consensus -> [meta.alias, consensus]}
        | groupTuple(by: 0)

        best_consensi = sample_results.consensus_stats
        | join(collected_consensi, by: 0)
        | get_best_consensus_files

        // define output
        ch_to_publish = Channel.empty()
        | mix(
            // centrifuge classification
            classification | map { meta, classification, report, kraken, html -> [classification, "$meta.alias/classification", "classification_${meta.alias}.tsv"] },
            classification | map { meta, classification, report, kraken, html -> [report, "$meta.alias/classification", "classification_report_${meta.alias}.tsv"] },
            classification | map { meta, classification, report, kraken, html -> [kraken, "$meta.alias/classification", "classification_kraken_${meta.alias}.tsv"] },
            classification | filter { params.output_krona_plot } | map { meta, classification, report, kraken, html -> [html, "$meta.alias/classification", "classification_${meta.alias}.html"] },
            // contigs
            contigs | map { meta, contigs -> [contigs, "$meta.alias/assembly", "${meta.mapping_target}.contigs.fasta"] },
            // consensus
            mapped_to_ref | map { meta, ref, bam, bai -> [ref, "$meta.alias/consensus", "${meta.consensus_target}.reference.fasta"] },
            mapped_to_ref | map { meta, ref, bam, bai -> [bam, "$meta.alias/consensus", "${meta.consensus_target}.reads.bam"] },
            mapped_to_ref | map { meta, ref, bam, bai -> [bai, "$meta.alias/consensus", "${meta.consensus_target}.reads.bam.bai"] },
            structural_variants | map { meta, variants -> [variants, "$meta.alias/consensus", "${meta.consensus_target}.structural_variants.vcf"] },
            consensi | map { meta, consensus -> [consensus, "$meta.alias/consensus", "${meta.consensus_target}.consensus.fasta"] },
            coverage | map { meta, coverage -> [coverage, "$meta.alias/consensus", "${meta.consensus_target}.depth.txt"] },
            // selected consensi
            best_consensi | map { alias, consensus_dir -> [consensus_dir, "$alias", "selected_consensus"] },
            // report and tables
            sample_results.report | map { alias, report -> [report, "$alias", "report_${alias}.html"] },
            sample_results.read_stats | map { alias, read_stats -> [read_stats, "$alias/tables", "reads.tsv"] },
            sample_results.contig_stats | map { alias, contig_stats -> [contig_stats, "$alias/tables", "contigs.tsv"] },
            sample_results.consensus_stats | map { alias, consensus_stats -> [consensus_stats, "$alias/tables", "consensus.tsv"] },
            // advanced output
            cleaned.reads | filter { params.output_cleaned_reads } | map { meta, reads -> [reads, "$meta.alias/clean", null] }
        )

    emit:
        results = ch_to_publish
}


def parseYamlToMap(Path yamlPath) {
    File yamlFile = yamlPath.toFile()
    def yaml = new org.yaml.snakeyaml.Yaml()
    def config = yaml.load(yamlFile.text)
    return config
}


def checkFlag(value) {
    return (value != null && value != 'false' && value != false)
}


params.download_db_all = checkFlag(params.download_db_all)
params.download_db_contamination = checkFlag(params.download_db_contamination)
params.download_db_virus = checkFlag(params.download_db_virus)
params.download_db_centrifuge = checkFlag(params.download_db_centrifuge)
params.download_db_update_existing = checkFlag(params.download_db_update_existing)


workflow db_update {
    main:
        def doUpdateContaminants = params.download_db_all || params.download_db_contamination
        def doUpdateVirus = params.download_db_all || params.download_db_virus
        def doUpdateCentrifuge = params.download_db_all || params.download_db_centrifuge

        download_config = db_update_get_config()

        config_dict = download_config
        | map { yaml -> parseYamlToMap(yaml) }

        db_virus = update_virus(config_dict, 'virus', doUpdateVirus)
        db_centrifuge = update_centrifuge(config_dict, 'centrifuge', doUpdateCentrifuge)
        db_contaminants = update_contaminants(config_dict, 'contaminants', doUpdateContaminants)

        Channel.empty()
        | mix(
            db_virus.database,
            db_centrifuge.database,
            db_contaminants.database
        )
        | toList
        | flatMap
        | data_base_transfer

        download_config
        | map { yaml -> [yaml, '', 'latest.yaml'] }
        | toList
        | flatMap
        | output
}


// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)

def doUpdate = (
    params.download_db_all
    || params.download_db_contamination
    || params.download_db_virus
    || params.download_db_centrifuge
)

def doBuildCustomDB = checkFlag(params.custom_db_do_build)

if (!doBuildCustomDB && !doUpdate && !params.fastq) {
    System.err.println("No fastqs provided! Exit.")
    System.exit(1)
} else if ([doUpdate, !!params.fastq, doBuildCustomDB].count { it } > 1) {
    System.err.println("Either download data base, build custom data base or analyse.")
    System.exit(1)
}

// check the system requirements before starting the workflow
if (doUpdate) {
    new SystemRequirements(true).checkSystemRequirements(
        params.download_db_min_disk_space_work_gb,
        params.download_db_min_disk_space_home_gb,
        params.download_db_min_ram_gb,
        params.download_db_min_cpus,
        params.database_defaults.base,
        session.workDir.toString()
    )
} else if (doBuildCustomDB) {
        new SystemRequirements(true).checkSystemRequirements(
        params.custom_db_min_disk_space_work_gb,
        params.custom_db_min_disk_space_out_gb,
        params.custom_db_min_ram_gb,
        params.custom_db_min_cpus,
        params.out_dir,
        session.workDir.toString()
    )
} else {
    new SystemRequirements(true).checkSystemRequirements(
        params.min_disk_space_work_gb,
        params.min_disk_space_out_gb,
        params.min_ram_gb,
        params.min_cpus,
        params.out_dir,
        session.workDir.toString()
    )
}


workflow {
    if(doUpdate) {
        db_update()
    } else if (doBuildCustomDB) {
        custom_data_base()
    } else {
        samples = fastq_ingress([
            "input": params.fastq,
            "sample_sheet": params.sample_sheet,
            "stats": true
        ])
        db_config = new DatabaseInput(params)
        pipeline(samples, db_config)
        pipeline.out.results
        | toList
        | flatMap
        | output
    }
}


workflow.onComplete {
    File outputFile = new File("${params.out_dir}/params.json")
    def json = new JsonBuilder(params)
    outputFile.withWriter('UTF-8') { writer -> writer.write(json.toPrettyString()) }
}


workflow.onError {
}
