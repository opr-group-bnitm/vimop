// Copyright (c) 2025 Outbreak Preparedness and Response Group at BNITM
// This file is part of ViMOP and is licensed under the MIT License.
// See the LICENSE file in the root of this repository for full license details.


nextflow.enable.dsl = 2


process db_update_get_config {
    label "general"
    errorStrategy "retry"
    maxRetries 3
    cpus 1
    output:
        path("db.yaml")
    """
    wget -o download_config.log -O db.yaml ${params.download_db_config_url}
    """
}


def checksumDir(String directory) {
    return "find ${directory}/. -type f -exec sha256sum {} \\; | sort | awk '{print \$1}' | sha256sum | awk '{print \$1}'"
}


process check_download_necessary {
    label "general"
    cpus 1 
    input:
        tuple val(do_update), val(expected_checksum), path("db_path")
    output:
        path("dummy.txt"), optional: true
    """
    if [ "${do_update}" == "true" ]
    then
        if [[ -d "db_path" ]]
        then
            if [[ "${params.download_db_update_existing}" == "true" ]]
            then
                checksum_old_data=\$(${checksumDir("db_path")})
                if [ \$checksum_old_data != ${expected_checksum} ]
                then
                    echo "update required because of different checksums (got \$checksum_old_data but expected ${expected_checksum})" > dummy.txt
                fi
            fi
        else
            echo "update required because the directory is missing" > dummy.txt
        fi
    fi
    """
}


process download_parts {
    // Parallelize independent, checksum-verified archive parts. The upstream
    // default is conservative for weak connections; this host has sufficient
    // bandwidth and disk for eight concurrent downloads.
    maxParallel = 8
    label "download"
    errorStrategy "retry"
    maxRetries 3
    cpus 1
    input:
        val(file_config)
    output:
        path("${file_config.name}")
    """
    wget -o dl_${file_config.name}.log -O ${file_config.name} ${file_config.url}

    checksum_new=\$(sha256sum ${file_config.name} | awk '{print \$1}')
    checksum_expected=${file_config.checksum}

    if [[ \$checksum_new != \$checksum_expected ]]
    then
        echo "Checksum mismatch!" >&2
        echo "Expected: \$checksum_expected" >&2
        echo "Computed: \$checksum_new" >&2
        exit 1
    fi
    """
}


process merge_parts_and_extract {
    label "general"
    cpus 1
    input:
        tuple val(db_config), path(parts), val(db_name)
    output:
        path("${db_name}")
    """
    # merge the files
    ls ${parts.join(' ')} | sort | xargs cat > merged.tar.xz

    # check the checksum
    checksum_new=\$(sha256sum merged.tar.xz | awk '{print \$1}')
    checksum_expected=${db_config.checksum_zipped}

    if [[ \$checksum_new != \$checksum_expected ]]
    then
        echo "Checksum mismatch for merged.tar.xz!" >&2
        echo "Expected: \$checksum_expected" >&2
        echo "Computed: \$checksum_new" >&2
        exit 1
    fi
    
    # extract
    tar -xf merged.tar.xz

    # check the directory checksum
    checksum_expected_directory=${db_config.checksum_directory}
    checksum_new_directory=\$(${checksumDir(db_name)})
    if [ "\$checksum_new_directory" != "\$checksum_expected_directory" ]
    then
        echo "Checksum mismatch for ${db_name}!" >&2
        echo "Expected: \$checksum_expected_directory" >&2
        echo "Computed: \$checksum_new_directory" >&2
        exit 1
    fi
    """
}


workflow update_data_base {
    take:
        config_dict
        database_name
        do_update
    main:
        db_config = config_dict
        | map { config -> config.sub_databases[database_name] }

        db_parts = db_config
        | map { db_config -> [do_update, db_config.checksum_directory, "${params.database_defaults.base}/${database_name}"] }
        | check_download_necessary
        | combine(db_config)
        | map { dummy, db_config -> db_config.files }
        | flatten
        | download_parts

        collected_parts = db_parts
        | collect
        | map { parts -> [parts] }

        database = db_config
        | combine(collected_parts)
        | map { db_config, parts -> [db_config, parts, database_name] }
        | merge_parts_and_extract
    emit:
        database = database
}


process data_base_transfer {
    label "general"
    cpus 1
    publishDir (
        params.database_defaults.base,
        mode: "copy",
        saveAs: {
            db
        }
    )
    input:
        path(db)
    output:
        path(db)
    """
    """
}
