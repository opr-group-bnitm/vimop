// Copyright (c) 2025 Outbreak Preparedness and Response Group at BNITM
// This file is part of ViMOP and is licensed under the MIT License.
// See the LICENSE file in the root of this repository for full license details.

import org.yaml.snakeyaml.Yaml


class DatabaseInput {
    String virusDir
    String virusConfigFileName
    String contaminantsDir
    String contaminationConfigFileName
    String classificationDir
    String classificationConfigFileName
    String blastDir
    String blastPrefix
    List<String> contaminationFilters
    List<String> contaminationFilterFiles
    List<Map<String, String>> virusTargets
    String classificationLibrary
    String virusTaxIDFile
    boolean doClassify
    boolean doFilterWithCentrifuge

    static void exitError(String message) {
        System.err.println(message)
        System.exit(1)
    }

    static String red(String message) {
        return "\u001B[31m${message}\u001B[0m"
    }

    static void assertDir(String dir) {
        def path = new File(dir)
        if (!path.exists() || !path.isDirectory()) {
            def fullMessage = (
                red("ERROR\n")
                + red("The directory '${path}' does not exist.\n")
                + red("Is the database installed? You can easily install it:\n")
                + "\n"
                + red("If you are using EPI2ME Desktop check the respective boxes in the section 'Setup'.\n")
                + "\n"
                + red("If you are running ViMOP from command line download the whole database in one run:\n")
                + red("nextflow run opr-group-bnitm/vimop --download_db_all\n")
                + "\n"
                + red("Exiting.")
            )
            exitError(fullMessage)
        }
    }

    static void assertFile(String fname) {
        def path = new File(fname)
        if (!path.exists() || !path.isFile()) {
            exitError("The file '${path}' does not exist. Exiting.")
        }
    }

    static String getDir(String dir, String defaultDir) {
        def outdir = dir ?: defaultDir
        assertDir(outdir)
        return outdir
    }

    static String getFile(String fname, String defaultFname) {
        def fnameOut = fname ?: defaultFname
        assertFile(fnameOut)
        return fnameOut
    }

    static Map readYamlConfig(String fname) {
        def yamlFile = new File(fname)
        def yamlParser = new Yaml()
        def configData = yamlParser.load(yamlFile.text)
        return configData
    }

    static String getFileFromConfig(Map config, String path, String key) {
        if (!config.containsKey(key)) {
            exitError("Missing required key '$key' in configuration")
        }
        def fname = "${path}/${config[key]}"
        assertFile(fname)
        return fname
    }

    static void checkUniqueKeys(List<Map<String, Object>> maps) {
        def seenKeys = new HashSet<String>()
        maps.each { map ->
            map.each { key, value ->
                def normalizedKey = key.toLowerCase()
                if (seenKeys.contains(normalizedKey)) {
                    exitError("Duplicate key detected: $key")
                }
                seenKeys.add(normalizedKey)
            }
        }
    }

    static Map upperCaseMap(List<Map<String, Object>> maps) {
        checkUniqueKeys(maps)
        def merged = maps.inject([:]) { mergedSoFar, map -> mergedSoFar + map }
        def mapWithUpperCaseKeys = merged.collectEntries {
            name, fasta -> [(name.toUpperCase()): fasta]
        }
        return mapWithUpperCaseKeys
    }

    static Map getVirusFilterPaths(Map yamlConfig) {
        def all = ["ALL": yamlConfig.all.fasta]
        def family_filters = yamlConfig.filters
        def curated = yamlConfig.curated.collectEntries {
            name, entries -> [(name): entries.fasta]
        }
        return upperCaseMap([all, family_filters, curated])
    }

    DatabaseInput(Map dbParams) {

        def baseDir = getDir(dbParams.base_db, dbParams.database_defaults.base)

        // contamination
        def cf = dbParams.contamination_filters?.toLowerCase()
        def doFilterContaminants = !(cf in [null, 'null', 'none', 'false'])

        this.contaminantsDir = getDir(
            dbParams.contaminants_db,
            "${baseDir}/${dbParams.database_defaults.contaminants}"
        )
        this.contaminationConfigFileName = getFile(
            dbParams.contaminants_db_config,
            "${this.contaminantsDir}/${dbParams.database_defaults.contaminants_db_config}"
        )
        def contaminationConfig = upperCaseMap([readYamlConfig(this.contaminationConfigFileName)['filters']])

        if (doFilterContaminants) {
            this.contaminationFilters = dbParams.contamination_filters.tokenize(",")
            this.contaminationFilterFiles = this.contaminationFilters.collect {
                contaminant -> getFileFromConfig(contaminationConfig, this.contaminantsDir, contaminant.toUpperCase())
            }
        } else {
            this.contaminationFilters = []
            this.contaminationFilterFiles = []
        }

        // virus
        def virusTargetNames = (dbParams.targets ?: '').tokenize(',')

        this.virusDir = getDir(
            dbParams.virus_db,
            "${baseDir}/${dbParams.database_defaults.virus}"
        )
        this.virusConfigFileName = getFile(
            dbParams.virus_db_config,
            "${this.virusDir}/${dbParams.database_defaults.virus_db_config}"
        )
        def virusConfig = readYamlConfig(this.virusConfigFileName)
        def filterFilenames = getVirusFilterPaths(virusConfig)
        this.virusTargets = virusTargetNames.collect {
            target -> [
                target: target,
                path: getFileFromConfig(filterFilenames, virusDir, target.toUpperCase())
            ]
        }

        // blast
        this.blastDir = "${this.virusDir}/${virusConfig.all.blast_db}"
        assertDir(this.blastDir)

        this.blastPrefix = virusConfig.all.blast_prefix
        assertFile("${this.blastDir}/${this.blastPrefix}.ndb")

        // classification
        this.doClassify = dbParams.centrifuge_do_classify
        this.doFilterWithCentrifuge = dbParams.centrifuge_do_classify && dbParams.centrifuge_do_filter

        this.classificationDir = getDir(
            dbParams.classification_db,
            "${baseDir}/${dbParams.database_defaults.classification}"
        )
        this.classificationConfigFileName = getFile(
            dbParams.classification_db_config,
            "${this.classificationDir}/${dbParams.database_defaults.classification_db_config}"
        )

        def classificationConfig = readYamlConfig(this.classificationConfigFileName)
        classificationConfig.files.each { fname -> assertFile("${this.classificationDir}/${fname}") }
        this.classificationLibrary = classificationConfig.index_name
        this.virusTaxIDFile = "${this.classificationDir}/${classificationConfig.virus_taxid_file}"
        assertFile(this.virusTaxIDFile)
    }
}
