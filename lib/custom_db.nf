// Copyright (c) 2025 Outbreak Preparedness and Response Group at BNITM
// This file is part of ViMOP and is licensed under the MIT License.
// See the LICENSE file in the root of this repository for full license details.


nextflow.enable.dsl = 2


def minCpus(int requestedCpus) {
    return Math.min(requestedCpus, params.custom_db_min_cpus as int)
}


def minRAM(int requestesRAM) {
    ram = Math.min(requestesRAM, params.custom_db_min_ram_gb as int)
    return "${ram} GB"
}


process build_contaminants_db {
    label "general"
    cpus 1
    input:
        path("contaminant_dir")
    output:
        path("contaminants")
    script:
    """
    set -euo pipefail
    shopt -s nullglob extglob

    # check fasta file validity
    for fname in contaminant_dir/*.{fna,fna.gz,fasta,fasta.gz}
    do
        if ! seqtk seq "\$fname" > /dev/null 2>&1
        then
            echo "Error: invalid FASTA file (\$fname)" >&2
            exit 1
        fi
    done

    # copy files (only if any exist)
    mkdir -p contaminants
    files=( contaminant_dir/*.{fna,fna.gz,fasta,fasta.gz} )
    if [ "\${#files[@]}" -gt 0 ]
    then
        cp "\${files[@]}" contaminants/
    fi

    cd contaminants

    # gzip any uncompressed files
    for fname in *.{fna,fasta}
    do
        gzip -f "\$fname"
    done

    # write YAML without quotes around keys; normalize keys to safe chars
    {
        echo "version: ${params.custom_db_contaminants_version}"
        echo "description: ${params.custom_db_contaminants_description}"
        echo "filters:"
        for fn in *.{fna.gz,fasta.gz}
        do
            filter_name=\${fn%.@(fasta|fna).gz}
            # normalize: keep A–Z, a–z, 0–9, dot, underscore, dash; others -> '_'
            filter_name=\${filter_name//[^A-Za-z0-9._-]/_}
            echo "  \${filter_name}: \${fn}"
        done
    } > contaminants.yaml
    """
}


process build_virus_db {
    label "general"
    cpus 1
    input:
        tuple path("virus_input_config.yaml"), path("ALL.fasta")
    output:
        path("virus")
    script:
    """
    set -euo pipefail

    mkdir -p virus
    cp ALL.fasta virus/ALL.fasta

    split_virus_db_files.py \\
        --yaml virus_input_config.yaml \\
        --fasta virus/ALL.fasta \\
        --outdir virus/ \\
        --description "${params.custom_db_virus_description}" \\
        --version "${params.custom_db_virus_version}" \\
        --blast-db blast_db \\
        --blast-prefix ALL \\
        --config-out virus.yaml \\
        --strict

    makeblastdb \\
        -dbtype nucl \\
        -in virus/ALL.fasta \\
        -out virus/blast_db/ALL \\
        -parse_seqids \\
        -blastdb_version 5
    """
}


process download_taxdump {
    label "general"
    cpus 1
    errorStrategy "retry"
    maxRetries 3
    output:
        path("taxdump.tar.gz") 
    script:
    """
    set -euo pipefail
    wget ${params.custom_db_centrifuge_taxonomy_url}
    """
}


process seqid_to_taxid {
    label "general"
    cpus { minCpus(4) }
    input:
        tuple path("sequences.fasta"), path("taxdump.tar.gz")
    output:
        tuple path("virus_taxids.txt"), path("seqid2taxid.map"), path("nodes.dmp"), path("names.dmp")
    script:
    """
    set -euo pipefail

    workdir=\$(pwd)
    taxonkit_datadir=\$workdir/taxonkit_data

    mkdir -p \$taxonkit_datadir
    mv taxdump.tar.gz \$taxonkit_datadir/taxdump.tar.gz
    cd \$taxonkit_datadir
    tar -xzf taxdump.tar.gz
    cd \$workdir

    # Get all virus taxids
    taxonkit list --data-dir \$taxonkit_datadir --ids ${params.custom_db_centrifuge_taxid_viruses} | sed 's/.* //' > virus_taxids.txt

    {
        echo "Bacteria\t${params.custom_db_centrifuge_taxid_bacteria}"
        echo "Archaea\t${params.custom_db_centrifuge_taxid_archaea}"
        echo "Eukaryota\t${params.custom_db_centrifuge_taxid_eukaryota}"
        echo "Viruses\t${params.custom_db_centrifuge_taxid_viruses}"
        echo "Viroids\t${params.custom_db_centrifuge_taxid_viroids}"
    } > kingdoms.taxids

    # Get the species and families from the fasta files
    extract_families_and_species.py \\
        --input sequences.fasta \\
        --taxon-table taxontable.tab \\
        --families families.txt \\
        --species species.txt

    taxonkit name2taxid --data-dir \$taxonkit_datadir --threads ${task.cpus} < species.txt > species.taxids
    taxonkit name2taxid --data-dir \$taxonkit_datadir --threads ${task.cpus} < families.txt > families.taxids

    seqids_to_taxids.py \\
        --taxon-table taxontable.tab \\
        --species-taxids species.taxids \\
        --family-taxids families.taxids \\
        --kingdom-taxids kingdoms.taxids \\
        --seqid-to-taxid seqid2taxid.map

    mv \$taxonkit_datadir/names.dmp .
    mv \$taxonkit_datadir/nodes.dmp .
    """
}


process build_centrifuge_index {
    label "centrifuge"
    cpus { minCpus(8) }
    memory { minRAM(30) }
    input:
        tuple path("sequences.fasta"), 
            path("virus_taxids.txt"),
            path("seqid2taxid.map"),
            path("nodes.dmp"),
            path("names.dmp")
    output:
        path("centrifuge")
    script:
    """
    set -euo pipefail

    mkdir -p centrifuge_tmp
    cd centrifuge_tmp

    cp ../virus_taxids.txt virus_taxids.txt

    ktUpdateTaxonomy.sh .

    centrifuge-build \\
        -p ${task.cpus} \\
        --conversion-table ../seqid2taxid.map \\
        --taxonomy-tree ../nodes.dmp \\
        --name-table ../names.dmp \\
        ../sequences.fasta \\
        all

    {
        echo "version: ${params.custom_db_centrifuge_version}"
        echo "description: \\"${params.custom_db_centrifuge_description}\\""
        echo "virus_taxid_file: virus_taxids.txt"
        echo "index_name: all"
        echo "files:"
        for f in all.*.cf
        do
            echo "- \$f"
        done
    } > centrifuge.yaml

    rm images.dmp
    mkdir -p taxonomy
    mv taxonomy.tab taxonomy/taxonomy.tab

    cd ..
    mv centrifuge_tmp centrifuge
    """
}


process custom_data_base_transfer {
    label "general"
    cpus 1
    publishDir (
        params.custom_db_outpath,
        mode: "copy",
        saveAs: {
            db
        }
    )
    input:
        path(db)
    output:
        path(db)
    script:
    """
    """
}


workflow custom_data_base {
    main:
        if(params.custom_db_contaminants_input_path) {
            contaminants_db = build_contaminants_db(params.custom_db_contaminants_input_path)
        } else {
            contaminants_db = Channel.empty()
        }

        if(params.custom_db_virus_yaml && params.custom_db_virus_fasta){
            virus_db = Channel.of([params.custom_db_virus_yaml, params.custom_db_virus_fasta])
            | build_virus_db
        } else {
            virus_db = Channel.empty()
        }

        if(params.custom_db_centrifuge_fasta) {
            centrifuge_db = download_taxdump
            | map { taxdump -> [params.custom_db_centrifuge_fasta, taxdump]}
            | seqid_to_taxid
            | map { virus_taxids, seqid2taxid, nodes, names -> [params.custom_db_centrifuge_fasta, virus_taxids, seqid2taxid, nodes, names] }
            | build_centrifuge_index
        } else {
            centrifuge_db = Channel.empty()
        }

        Channel.empty()
        | mix(
            contaminants_db,
            centrifuge_db,
            virus_db
        )
        | toList
        | flatMap
        | custom_data_base_transfer
}
