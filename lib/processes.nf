// Copyright (c) 2025 Outbreak Preparedness and Response Group at BNITM
// This file is part of ViMOP and is licensed under the MIT License.
// See the LICENSE file in the root of this repository for full license details.


nextflow.enable.dsl = 2


def minCpus(int requestedCpus) {
    return Math.min(requestedCpus, params.min_cpus as int)
}


def minRAM(int requestesRAM) {
    ram = Math.min(requestesRAM, params.min_ram_gb as int)
    return "${ram} GB"
}


process lengths_and_qualities {
    label "general"
    cpus 1
    input:
        tuple val(samplename), path(reads)
    output:
        tuple val(samplename), path("read_len_qual.tsv")
    """
    #!/usr/bin/env python
    from Bio import SeqIO
    from pandas import DataFrame
    lines = [
        (len(record.seq), sum(record.letter_annotations['phred_quality']) / len(record.seq))
        for record in SeqIO.parse('${reads}', 'fastq')
        if len(record.seq) > 0
    ]
    DataFrame(lines, columns=['Length', 'Quality']).to_csv('read_len_qual.tsv', sep='\t')
    """
}


process trim {
    label "general"
    cpus 1
    input:
        tuple val(meta), path("demultiplexed.fastq.gz")
    output:
        tuple val(meta), path("trimmed.fastq")
    """
    if [[ ${params.trim_length} == 0 ]]
    then
        gunzip -c demultiplexed.fastq.gz > trimmed.fastq
    else
        seqtk trimfq -b ${params.trim_length} -e ${params.trim_length} demultiplexed.fastq.gz > trimmed.fastq
    fi
    """
}

process classify_centrifuge {
    label "centrifuge"
    cpus { minCpus(12) }
    memory { minRAM(28) }
    input:
        tuple val(meta), path("seqs.fastq"), path(db_path), val(target_db)
    output:
        tuple val(meta),
            path("classification.tsv"),
            path("classification_report.tsv"),
            path("classification_kraken.tsv"),
            path("classification.html")
    """
    if [[ -s seqs.fastq ]]
    then
        centrifuge \\
            -p ${task.cpus} \\
            --mm \\
            -x ${db_path}/${target_db} \\
            -U seqs.fastq \\
            --report-file classification_report.tsv \\
            -S classification.tsv
        centrifuge-kreport \\
            -x ${db_path}/${target_db} classification.tsv > classification_kraken.tsv
        ktImportTaxonomy \\
            -tax ${db_path}/taxonomy \\
            -m 3 -t 5 \\
            classification_kraken.tsv \\
            -o classification.html
    else
        # write a report even if no sequences were available
        touch classification.tsv
        touch classification_report.tsv
        touch classification_kraken.tsv

        # Create an informative HTML file telling the user that no sequences were found
        cat <<EOF > classification.html
<!DOCTYPE html>
<html>
<head>
    <title>No Sequences Found</title>
    <style>
        body { font-family: Arial, sans-serif; text-align: center; margin-top: 50px; }
        h1 { color: red; }
    </style>
</head>
<body>
    <h1>No sequences were found</h1>
    <p>For this barcode no sequences could be classified by centrifuge.</p>
</body>
</html>
EOF
    fi
    """
}


process classify_contigs {
    label "centrifuge"
    cpus { minCpus(12) }
    memory = { minRAM(28) }
    input:
        tuple val(meta), path("contigs.fasta"), path(db_path), val(target_db)
    output:
        tuple val(meta), 
        path("classification_report.tsv"),
        path("classification.tsv")
    """
    if [[ -s contigs.fasta ]]
    then
        centrifuge \
            -p ${task.cpus} \
            --mm \
            -x ${db_path}/${target_db} \
            -f contigs.fasta \
            --report-file classification_report.tsv \
            -S classification.tsv
        centrifuge-kreport \
            -x ${db_path}/${target_db} classification.tsv > classification_kraken.tsv
    else
        # write a report even if no sequences were available
        echo -e "taxID\tname\ttaxRank" > classification_report.tsv
        echo -e "taxID\treadID" > classification.tsv
    fi
    """
}


process extract_contig_classification {
    label "general"
    cpus 1
    input:
        tuple val(meta), path("classification_report.tsv"), path("classification.tsv")
    output:
        tuple val(meta), path("classification_summary_${meta.mapping_target}.csv")
    """
    #!/usr/bin/env python
    import pandas as pd
    report = pd.read_csv('classification_report.tsv', sep='\t')
    classification = pd.read_csv('classification.tsv', sep='\t')

    if report.empty or classification.empty:
        with open('classification_summary_${meta.mapping_target}.csv', 'w') as f:
            f.write('readID,taxRank,name')
    
    else:
        summary = classification[['taxID', 'readID']].merge(report[['taxID', 'name', 'taxRank']], on='taxID', how='left')
        summary['taxRank'] = summary['taxRank'].fillna('unclassified')
        summary['name'] = summary['name'].fillna('no name')
        summary = summary.drop(columns=['taxID'])
        summary.to_csv('classification_summary_${meta.mapping_target}.csv', index=False)
    """
}


process no_contig_classification {
    label "general"
    cpus 1
    input:
        val(meta)
    output:
        tuple val(meta), path("classification_summary_${meta.mapping_target}.csv")
    """
    #!/usr/bin/env python
    with open('classification_summary_${meta.mapping_target}.csv', 'w') as f:
        f.write('readID,taxRank,name')
    """
}


def seqstatsHeader(String fnameOut) {
    return """
    echo "step\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len" > $fnameOut
    """
}


def isNotEmptyCheck(String fname) {
    if (fname.endsWith(".gz")) {
        return " \$(zcat ${fname} | wc -l) -ne 0 "
    }
    return " -s ${fname} "
}


def seqstatsLine(String rowIdentifier, String fnameIn, String fnameOut) {
    return """
    echo -n "${rowIdentifier}\t" >> ${fnameOut}
    if [[ ${isNotEmptyCheck(fnameIn)} ]]
    then
        seqkit stats -T ${fnameIn} \\
        | tail -n 1 \\
        | awk '{print \$4"\t"\$5"\t"\$6"\t"\$7"\t"\$8}' >> ${fnameOut}
    else
        echo "0\t0\t0\t0.0\t0" >> ${fnameOut}
    fi
    """
}


process filter_with_centrifuge {
    label "general"
    cpus 1
    memory { minRAM(3) }    
    input:
        tuple val(meta), path('seqs.fastq'), path('classification.tsv'), path('virus_taxids.txt')
    output:
        tuple val(meta), path('centrifuge_filtered.fastq'), path('stats.tsv')
    """
    ${seqstatsHeader("stats.tsv")}
    ${seqstatsLine("nofilter", "seqs.fastq", "stats.tsv")}

    filter_reads_by_centrifuge_classification.py \\
        --centrifuge classification.tsv \\
        --fastq seqs.fastq \\
        --out centrifuge_filtered.fastq \\
        --virus-taxids virus_taxids.txt \\
        --min-score ${params.centrifuge_filter_min_score}

    ${seqstatsLine("centrifuge-filtered", "centrifuge_filtered.fastq", "stats.tsv")}
    """
}


process read_stats {
    label "general"
    cpus 1
    input:
        tuple val(meta), path('seqs.fastq')
    output:
        tuple val(meta), path('seqs.fastq'), path('stats.tsv')
    """
    ${seqstatsHeader("stats.tsv")}
    ${seqstatsLine("nofilter", "seqs.fastq", "stats.tsv")}
    """
}


process filter_contaminants {
    label "general"
    cpus { minCpus(8) }
    memory { minRAM(20) }
    input:
        tuple val(meta), path('seqs.fastq'), path('read_stats.tsv'), path(db_paths), val(contaminants)
    output:
        tuple val(meta), path('filtered.fastq'), emit: reads
        tuple val(meta.alias), path("clean_stats.tsv"), emit: stats
    """
    cp read_stats.tsv stats_tmp.tsv

    contaminants=(${contaminants.join(" ")})
    db_paths=(${db_paths.join(" ")})
    fn_input=seqs.fastq

    for i in "\${!contaminants[@]}"
    do
        c=\${contaminants[\$i]}
        db=\${db_paths[\$i]}
        fn_out=filtered_no_\${c}.fastq

        fn_sam=filter_\${c}.sam

        minimap2 \\
            -x map-ont \\
            -a \$db \\
            -t ${task.cpus} \\
            \$fn_input \\
        | samtools fastq \\
            -f 4 \\
            --reference \$db \\
            -1 \${fn_out} -2 \${fn_out} -0 \${fn_out} -s \${fn_out} -n

        ${seqstatsLine("\${c}", "\$fn_out", "stats_tmp.tsv")}

        fn_input=\$fn_out
    done
    cp \$fn_input filtered.fastq
    mv stats_tmp.tsv clean_stats.tsv
    """
}


process filter_virus_target {
    label "general"
    cpus { minCpus(4) }
    memory { minRAM(20) }
    input:
        tuple val(meta), path('seqs.fastq'), path(target)
    output:
        tuple val(meta), path('filtered.fastq')
    """
    ref_size=\$(stat -Lc%s "${target}")

    if [[ "\$ref_size" -gt ${params.target_filter_minimap_index_split_thresh * 1024 * 1024} ]]
    then
        echo "Large reference — using split index"
        minimap2 -I 2G -d map_index.mmi --split-prefix map_index "${target}"
        ref=map_index.mmi
    else
        echo "Small reference - no need to split the index"
        ref="${target}"
    fi
    minimap2 \\
        -ax map-ont \\
        --split-prefix map_index \\
        --secondary=no \\
        "\$ref" \\
        -t ${task.cpus} \\
        seqs.fastq \\
    | samtools fastq \\
        -F 4 > filtered.fastq

    rm -f map_index.mmi
    """
}


process assemble_canu {
    label "canu"
    cpus { minCpus(16) }
    memory { minRAM(24) }
    input:
        tuple val(meta), path("seqs.fastq"), val(get_reads_if_no_contigs)
    output:
        tuple val(meta), path("asm.contigs.fasta"), emit: contigs
        tuple val(meta), path("assembly_stats_${meta.mapping_target}.tsv"), emit: stats
    """
    ${task.ext.conda_init ?: ''}
    ${task.ext.conda_activate ?: ''}

    seqtk seq -L ${params.canu_min_read_length} seqs.fastq > filtered.minlen.fastq

    ${task.ext.conda_deactivate ?: ''}

    outdir=.
    set +e
    canu \\
        -nanopore-raw filtered.minlen.fastq \\
        -fast \\
        -p asm \\
        -d \$outdir \\
        genomeSize=${params.canu_genome_size} \\
        minReadLength=${params.canu_min_read_length} \\
        minOverlapLength=${params.canu_min_overlap_length} \\
        corOutCoverage=${params.canu_cor_out_coverage} \\
        readSamplingBias=${params.canu_read_sampling_bias} \\
        stopOnLowCoverage=${params.canu_stop_on_low_coverage} \\
        minInputCoverage=${params.canu_min_input_coverage} \\
        maxInputCoverage=${params.canu_max_input_coverage} \\
        maxThreads=${task.cpus} \\
        maxMemory=${task.memory.toGiga()}g
    set -e

    ${task.ext.conda_activate ?: ''}

    touch asm.contigs.fasta

    if [[ ! -f asm.correctedReads.fasta.gz ]]
    then
        touch asm.correctedReads.fasta
        gzip asm.correctedReads.fasta
    fi

    ${seqstatsHeader("stats.tsv")}
    ${seqstatsLine("raw_" + meta.mapping_target, "filtered.minlen.fastq", "stats.tsv")}
    ${seqstatsLine("corrected_" + meta.mapping_target, "asm.correctedReads.fasta.gz", "stats.tsv")}
    ${seqstatsLine("contigs_" + meta.mapping_target, "asm.contigs.fasta", "stats.tsv")}

    if [[ "${get_reads_if_no_contigs}" == "true" && ! -s asm.contigs.fasta ]]
    then
        # no contigs, get some reads, prefer corrected reads if they are available
        if [[ ${isNotEmptyCheck("asm.correctedReads.fasta.gz")} ]]
        then
            seqkit sort -l -r asm.correctedReads.fasta.gz \\
            | seqkit head -n ${params.nocontigs_max_reads_precluster} > longest_reads.fasta

            read_type=corrected
        else
            seqkit fq2fa seqs.fastq | seqkit sort -l -r \\
            | seqkit head -n ${params.nocontigs_max_reads_precluster} > longest_reads.fasta

            read_type=raw
        fi

        if [[ -s longest_reads.fasta ]]
        then
            cd-hit-est \\
                -n ${params.nocontigs_cdhit_wordlen} \\
                -c ${params.nocontigs_cdhit_thresh} \\
                -T ${task.cpus} \\
                -M ${task.memory.toMega()} \\
                -i longest_reads.fasta \\
                -o clustered.fasta

            seqkit head -n ${params.nocontigs_nreads} clustered.fasta > selected.fasta

            # rename the reads and add header infos
            rename_seqs.py \\
                --prefix \${read_type}_read_ \\
                --input selected.fasta \\
                --output renamed.fasta

            ${seqstatsLine("longestreads_" + meta.mapping_target, "renamed.fasta", "stats.tsv")}

            # replace the empty contigs-file with the longest reads
            mv renamed.fasta asm.contigs.fasta
        fi
    fi

    mv stats.tsv assembly_stats_${meta.mapping_target}.tsv
    """
}


process reassemble_canu {
    label "canu"
    cpus { minCpus(16) }
    memory { minRAM(24) }
    input:
        tuple val(meta), path("contigs_*.fasta"), path("seqs.fastq")
    output:
        tuple val(meta), path("reassembly.contigs.fasta"), emit: contigs
        tuple val(meta), path("reassembly.stats.tsv"), emit: stats
    """
    ${task.ext.conda_init ?: ''}
    ${task.ext.conda_activate ?: ''}

    ${seqstatsHeader("reassembly.stats.tsv")}
    touch new.contigs.fasta

    seqtk seq -L ${params.canu_min_read_length} seqs.fastq > filtered.fastq

    ${seqstatsLine("input_reassembly", "filtered.fastq", "reassembly.stats.tsv")}

    i=0
    for fn_contig in contigs_*.fasta
    do
        rename_seqs.py \\
            --prefix inputcontigs\${i}_ \\
            --input \$fn_contig \\
            --output contigs_renamed\$i.fasta
        i=\$((i + 1))
    done
    cat contigs_renamed*.fasta > last.contigs.fasta
    rm contigs_renamed*.fasta

    i=0
    while [[ \$i -lt ${params.reassemble_max_iter} ]]
    do
        i=\$((i + 1))

        minimap2 \\
            -ax map-ont \\
            -t ${task.cpus} \\
            --secondary=no \\
            last.contigs.fasta \\
            filtered.fastq \\
        | samtools fastq -f 4 > re.filtered.fastq

        mv re.filtered.fastq filtered.fastq

        if [[ ! -s filtered.fastq ]]
        then
            break
        fi

        ${task.ext.conda_deactivate ?: ''}

        outdir=canu_output_\$i
        set +e
        canu \\
            -nanopore-raw filtered.fastq \\
            -fast \\
            -p asm \\
            -d \$outdir \\
            genomeSize=${params.canu_genome_size} \\
            minReadLength=${params.canu_min_read_length} \\
            minOverlapLength=${params.canu_min_overlap_length} \\
            corOutCoverage=${params.canu_cor_out_coverage} \\
            readSamplingBias=${params.canu_read_sampling_bias} \\
            stopOnLowCoverage=${params.canu_stop_on_low_coverage} \\
            minInputCoverage=${params.canu_min_input_coverage} \\
            maxInputCoverage=${params.canu_max_input_coverage} \\
            maxThreads=${task.cpus} \\
            maxMemory=${task.memory.toGiga()}g
        set -e

        ${task.ext.conda_activate ?: ''}

        if [[ ! -s \$outdir/asm.contigs.fasta ]]
        then
            break
        fi

        cd-hit-est \\
            -n ${params.reassemble_cdhit_wordlen} \\
            -c ${params.reassemble_cdhit_thresh} \\
            -i \$outdir/asm.contigs.fasta \\
            -o clustered.contigs.\$i.fasta

        rename_seqs.py \\
            --prefix newcontigs\${i}_ \\
            --input clustered.contigs.\$i.fasta \\
            --output contigs.\$i.fasta

        cat contigs.\$i.fasta >> new.contigs.fasta

        ${seqstatsLine("raw_reassembly\$i", "filtered.fastq", "reassembly.stats.tsv")}
        ${seqstatsLine("corrected_reassembly\$i", "\$outdir/asm.correctedReads.fasta.gz", "reassembly.stats.tsv")}
        ${seqstatsLine("contigs_reassembly\$i", "contigs.\$i.fasta", "reassembly.stats.tsv")}
    
        mv contigs.\$i.fasta last.contigs.fasta
    done

    mv new.contigs.fasta reassembly.contigs.fasta
    """
}


process no_assembly {
    label "general"
    cpus 1
    input:
        val(meta)
    output:
        tuple val(meta), path("noassembly.contigs.fasta"), emit: contigs
        tuple val(meta), path("noassembly.stats.tsv"), emit: stats
    """
    touch noassembly.contigs.fasta
    ${seqstatsHeader("noassembly.stats.tsv")}
    """
}


process pop_bubbles {
    label "general"
    cpus 1
    input:
        tuple val(meta), path('canu_contigs.fasta')
    output:
        tuple val(meta), path('nobubbles.fasta')
    """
    #!/usr/bin/env python
    from Bio import SeqIO
    with open('nobubbles.fasta', "w") as f_out:
        for record in SeqIO.parse('canu_contigs.fasta', 'fasta'):
            if "bubble=yes" not in record.description:
                SeqIO.write(record, f_out, 'fasta')
    """
}


process prepare_blast_search {
    label "general"
    cpus 1
    input:
        tuple val(meta), path("contigs.fasta")
    output:
        tuple val(meta), path("sorted-contigs.fasta")
    """
    seqkit sort --by-length --reverse contigs.fasta \
    | seqkit rename > sorted-contigs.fasta
    """
}


process canu_contig_info {
    label "general"
    cpus 1
    input:
        tuple val(meta), path("contigs.fasta")
    output:
        tuple val(meta), path("contig-info-${meta.mapping_target}.tsv")
    """
    #!/usr/bin/env python
    import pandas as pd
    with open('contigs.fasta') as f_in:
        fasta_headers = [
            line[1:].strip().split()
            for line in f_in
            if line.startswith('>')
        ]
    info = [
        {
            'WorkflowMappingTarget': '${meta.mapping_target}',
            'Contig': header[0],
            **dict(entry.split('=') for entry in header[1:])
        }
        for header in fasta_headers
    ]
    columns = [
        'WorkflowMappingTarget',
        'Contig',
        'len',
        'reads',
    ]
    df = pd.DataFrame(info, columns=columns)
    df.to_csv('contig-info-${meta.mapping_target}.tsv', sep='\\t', index=False)
    """
}


process blast {
    label "general"
    cpus { minCpus(4) }
    input:
        tuple val(meta), path("sorted-contigs.fasta"), path(db_path), val(target)
    output:
        tuple val(meta), path("blast-results.xml")
    """
    blastn \
        -num_threads ${task.cpus} \
        -db ${db_path}/${target} \
        -query sorted-contigs.fasta \
        -out blast-results.xml \
        -outfmt 5 \
        -max_target_seqs 1
    """
}


process extract_blasthits {
    label "general"
    cpus 1
    input:
        tuple val(meta), path("blast-results.xml")
    output:
        tuple val(meta), path("blast-hits-${meta.mapping_target}.csv")
    """
    #!/usr/bin/env python
    import pandas as pd
    from xml.etree import ElementTree
    column_names = [
        'Query',
        'Reference',
        'Description',
        'Family',
        'Organism',
        'Segment',
        'Orientation',
        'HitLength',
        'Bitscore',
        'QueryFrom',
        'QueryTo',
        'HitFrom',
        'HitTo',
        'IdenticalPositions',
        'AlignmentLength',
        'Gaps',
        'WorkflowMappingTarget',
    ]
    rows = []
    try:
        root = ElementTree.parse('blast-results.xml').getroot()
        for iteration in root.findall(".//Iteration"):
            query_name = iteration.find('Iteration_query-def').text.split()[0]
            hit = iteration.find("./Iteration_hits/Hit[Hit_num='1']")
            if hit is not None:
                try:
                    accession = hit.find('Hit_accession').text
                    fasta_header = hit.find('Hit_def').text.replace(',', '').replace(';', '')
                    # sometimes, the description contains the | symbol, so we have to split around it.
                    id_and_description, family, organism, orientation, segment = fasta_header.rsplit('|', 4)
                    description = id_and_description.split('|', 1)[1]
                    hit_length = int(hit.find('Hit_len').text)
                    hsp = hit.find("./Hit_hsps/Hsp[Hsp_num='1']")
                    if hsp is not None:
                        bit_score = float(hsp.find("Hsp_bit-score").text)
                        query_from = int(hsp.find("Hsp_query-from").text)
                        query_to = int(hsp.find("Hsp_query-to").text)
                        hit_from = int(hsp.find("Hsp_hit-from").text)
                        hit_to = int(hsp.find("Hsp_hit-to").text)
                        identity = int(hsp.find("Hsp_identity").text)
                        alignment_length = int(hsp.find("Hsp_align-len").text)
                        gaps = int(hsp.find("Hsp_gaps").text)
                        rows.append({
                            'Query': query_name.strip(),
                            'Reference': accession.strip(),
                            'Description': description.strip(),
                            'Family': family.strip(),
                            'Organism': organism.strip(),
                            'Segment': segment.strip(),
                            'Orientation': orientation.strip(),
                            'HitLength': hit_length,
                            'Bitscore': bit_score,
                            'QueryFrom': query_from,
                            'QueryTo': query_to,
                            'HitFrom': hit_from,
                            'HitTo': hit_to,
                            'IdenticalPositions': identity,
                            'AlignmentLength': alignment_length,
                            'Gaps': gaps,
                            'WorkflowMappingTarget': '${meta.mapping_target}',
                        })
                except (ValueError, AttributeError):
                    continue
    except ElementTree.ParseError:
        # empty file with no blast hits
        pass
    hits = pd.DataFrame(rows, columns=column_names)
    hits.to_csv('blast-hits-${meta.mapping_target}.csv', header=True, index=False)
    """
}


process get_ref_fasta {
    label "general"
    cpus 1
    input:
        tuple val(ref_id), path(db_path), val(db_name)
    output:
        tuple val(ref_id), path("ref.fasta")
    """
    blastdbcmd -entry ${ref_id} -db ${db_path}/${db_name} -out ref.fasta
    """
}


process split_custom_ref {
    label "general"
    cpus 1
    input:
        path fasta_file
    output:
        path("split_refs/*.fasta")
    """
    #!/usr/bin/env python
    import os
    from Bio import SeqIO

    fasta_file = '${fasta_file}'
    os.makedirs('split_refs', exist_ok=True)
    records = list(SeqIO.parse(fasta_file, 'fasta'))

    if not records:
        print(f'WARNING: No records found in input FASTA {fasta_file}. Skipping.', flush=True)
        exit(0)

    for record in records:
        ref_id = record.id
        file_path = f'split_refs/{ref_id}.fasta'
        SeqIO.write(record, file_path, 'fasta')
    """
}


process map_to_ref {
    label "general"
    cpus { minCpus(8) }
    memory { minRAM(10) }
    input:
        tuple val(meta), path("trimmed.fastq"), path("ref.fasta")
    output:
        tuple val(meta), path("ref.fasta"), path("sorted.bam"), path("sorted.bam.bai")
    """
    minimap2 -ax map-ont ref.fasta trimmed.fastq \
        -t ${task.cpus} --secondary=no \
        -w ${params.map_to_target_minimap_window_size} \
        -k ${params.map_to_target_minimap_kmer_size} \
    | samtools view -F 4 -b -o mapped.bam

    samtools sort --threads ${task.cpus} -o sorted.bam mapped.bam
    samtools index sorted.bam
    """
}

process map_to_sv_consensus {
    label "general"
    cpus { minCpus(8) }
    memory { minRAM(10) }
    input:
        tuple val(meta),
            path("trimmed.fastq"),
            path("structural_variants.vcf"),
            path("sv_consensus.fasta"),
            path("mapped_to_ref.bam"),
            path("mapped_to_ref.bam.bai")
    output:
        tuple val(meta),
            path("sv_consensus.fasta"),
            path("mapped_to_sv_consensus.bam"),
            path("mapped_to_sv_consensus.bam.bai")
    """
    set +e
    variant_count=\$(grep -c -v '^#' structural_variants.vcf)
    set -e

    if [[ "\$variant_count" -eq 0 ]]
    then
        # no need to map again
        mv mapped_to_ref.bam mapped_to_sv_consensus.bam
        mv mapped_to_ref.bam.bai mapped_to_sv_consensus.bam.bai
    else
        minimap2 -ax map-ont sv_consensus.fasta trimmed.fastq \\
            -t ${task.cpus} --secondary=no \\
            -w ${params.map_to_target_minimap_window_size} \\
            -k ${params.map_to_target_minimap_kmer_size} \\
        | samtools view -F 4 -b -o mapped.bam

        samtools sort --threads ${task.cpus} -o mapped_to_sv_consensus.bam mapped.bam
        samtools index mapped_to_sv_consensus.bam
    fi
    """
}


process calc_coverage {
    label "general"
    cpus 1
    input:
        tuple val(meta), path("sorted.bam"), path("sorted.bam.bai")
    output:
        tuple val(meta), path("coverage.txt")
    """
    samtools depth -aa -J sorted.bam > coverage.txt
    """
}


process simplify_reference_fasta {
    label "general"
    cpus 1
    input:
        tuple val(meta), path("ref.fasta")
    output:
        tuple val(meta), path("simple_ref.fasta")
    """
    #!/usr/bin/env python
    from Bio import SeqIO

    ref = SeqIO.read('ref.fasta', 'fasta')
    seq_out = ''.join([
        base if base in 'ATCG' else 'N'
        for base in str(ref.seq).upper()
    ])

    with open('simple_ref.fasta', 'w') as f_out:
        f_out.write(f">{ref.id}\\n{seq_out}\\n") 
    """
}


process subsample_alignments {
    label "general"
    cpus 1
    input:
        tuple val(meta), path("ref.fasta"), path("in.bam"), path("in.bam.bai")
    output:
        tuple val(meta), path("ref.fasta"), path("out.bam"), path("out.bam.bai")
    """
    set -euo pipefail

    # Ensure jvarkit can find java where its wrapper expects it
    if [[ -n "\$CONDA_PREFIX" && ! -x "\$CONDA_PREFIX/bin/java" ]]
    then
        ln -sf "\$(command -v java)" "\$CONDA_PREFIX/bin/java"
    fi

    target_coverage=${params.sv_cap_coverage}
    if [[ \$target_coverage -gt 0 ]]
    then
        samtools dict ref.fasta -o ref.dict

        jvarkit sortsamrefname in.bam -o refname_sorted.bam
        jvarkit biostar154220 -R ref.fasta -n \$target_coverage refname_sorted.bam -o capped.bam

        samtools sort capped.bam -o out.bam
        samtools index out.bam
    else
        mv in.bam out.bam
        mv in.bam.bai out.bam.bai
    fi
    """
}


process sniffles {
    label "structural_variants"
    cpus { minCpus(2) }
    memory { minRAM(10) } 
    input:
        tuple val(meta),
            path("ref.fasta"),
            path("mapped_to_ref.bam"),
            path("mapped_to_ref.bam.bai")
    output:
        tuple val(meta),
            path("sv.filtered.vcf")
    """
    if [[ \$(samtools view mapped_to_ref.bam | wc -l) -ne 0 ]]
    then
        set +e
        sniffles \\
            -t ${task.cpus} \\
            --input mapped_to_ref.bam \\
            --reference ref.fasta \\
            --vcf sv.vcf \\
            --minsupport ${params.sniffles_min_support} \\
            --minsvlen ${params.sniffles_min_sv_len}
        set -e

        bcftools view -i 'INFO/AF>=${params.sniffles_min_variant_allele_fraction}' sv.vcf -Ov -o sv.filtered.vcf
    else
        # create an empty vcf file
        refid=\$(head -n 1 ref.fasta | sed 's/^>//' | awk '{print \$1}')
        reflen=\$(grep -v '^>' ref.fasta | tr -d '\\n' | wc -c)

        # Create the minimal VCF
        {
            echo "##fileformat=VCFv4.2"
            echo "##contig=<ID=\$refid,length=\$reflen>"
            echo -e "#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO"
        } > sv.filtered.vcf
    fi
    """
}


process cutesv {
    label "structural_variants"
    cpus { minCpus(2) }
    memory { minRAM(5) } 
    input:
        tuple val(meta),
            path("ref.fasta"),
            path("mapped_to_ref.bam"),
            path("mapped_to_ref.bam.bai")
    output:
        tuple val(meta),
            path("sv.filtered.vcf")
    """
    if [[ \$(samtools view mapped_to_ref.bam | wc -l) -ne 0 ]]
    then
        mkdir work_sv

        cuteSV mapped_to_ref.bam ref.fasta sv.vcf work_sv \\
            --genotype \\
            --min_support ${params.cutesv_min_support} \\
            --min_mapq ${params.cutesv_min_mapq} \\
            --min_read_len ${params.cutesv_min_read_len} \\
            --min_size ${params.cutesv_min_sv_len}

        # only keep precise variants with sufficiently high allele fraction
        bcftools view \\
            -f PASS \\
            -i 'INFO/PRECISE=1 && INFO/AF>=${params.cutesv_min_variant_allele_fraction}' \\
            sv.vcf -o filtered.vcf

        # only keep variants supported by bcftools consensus
        bcftools view -i '(ALT !~ "<.*>") || (ALT == "<DEL>")' filtered.vcf -o sv.filtered.vcf
    else
        # create an empty vcf file
        refid=\$(head -n 1 ref.fasta | sed 's/^>//' | awk '{print \$1}')
        reflen=\$(grep -v '^>' ref.fasta | tr -d '\\n' | wc -c)

        # Create the minimal VCF
        {
            echo "##fileformat=VCFv4.2"
            echo "##contig=<ID=\$refid,length=\$reflen>"
            echo -e "#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO"
        } > sv.filtered.vcf
    fi
    """
}


process structural_variant_consensus {
    label "structural_variants"
    input:
        tuple val(meta),
            path("ref.fasta"),
            path("sv.filtered.vcf")
    output:
        tuple val(meta),
            path("sv_consensus.fasta")
    """
    if [[ \$(bcftools view -H sv.filtered.vcf | wc -l) -ne 0 ]]
    then
      bcftools sort sv.filtered.vcf -Oz -o structural_variants.vcf.gz
      tabix structural_variants.vcf.gz
      bcftools consensus \\
        --fasta-ref ref.fasta \\
        -o sv_consensus.fasta structural_variants.vcf.gz
    else
      cp ref.fasta sv_consensus.fasta
    fi
    """
}


def medakaGetModel() {
    return """
    model_choice="default"

    if [[ ${params.medaka_consensus_model} == "auto" ]]
    then
        set +e
        model_path=\$(medaka tools resolve_model --auto_model variant trimmed.fastq 2>/dev/null)
        exit_code=\$?
        if [[ \$exit_code -eq 0 && -n "\$model_path" ]]
        then
            model_choice=\$model_path
        fi
        set -e
    else
        model_choice=${params.medaka_consensus_model}:variant
    fi
    """
}


def medakaVariantConsensus(String alias) {
    return """
    \$medaka_cmd

    medaka vcf consensus.hdf ref.fasta variants.vcf --gvcf

    bcftools sort variants.vcf -o sorted.vcf

    medaka tools annotate sorted.vcf ref.fasta sorted.bam annotated.vcf

    consensus_from_medaka_gvcf.py \\
        --ref ref.fasta \\
        --gvcf annotated.vcf \\
        --min_depth ${params.consensus_min_depth} \\
        --min_qual ${params.consensus_medaka_min_qual} \\
        --model \$model_choice \\
        --sample ${alias} \\
        --out consensus.fasta
    """
}


process medaka_variant_consensus {
    label "medaka"
    cpus { minCpus(2) }
    memory { minRAM(24) }
    input:
        tuple val(meta),
            path("ref.fasta"),
            path("sorted.bam"),
            path("sorted.bam.bai"),
            path("trimmed.fastq")
    output:
        tuple val(meta), path("consensus.${meta.consensus_target}.fasta")
    """
    medaka_cmd="medaka inference sorted.bam consensus.hdf --threads ${task.cpus}"

    ${medakaGetModel()}

    if [[ \$model_choice != "default" ]]
    then
        medaka_cmd="\$medaka_cmd --model \${model_choice}"
    fi

    ${medakaVariantConsensus(meta.alias)}

    mv consensus.fasta consensus.${meta.consensus_target}.fasta
    """
}


def simpleConsensus(String alias) {
    return """
    samtools consensus \\
        -a --mark-ins --show-del yes --show-ins yes \\
        --mode simple \\
        --min-depth ${params.consensus_min_depth} \\
        --call-fract ${params.consensus_min_share} \\
        -f fasta sorted.bam \\
        > consensus_draft.fasta

    refid=\$(head -n 1 ref.fasta | awk '{print substr(\$1, 2)}')
    echo ">consensus method=samtools_simple reference=\$refid sample=${alias}" > consensus_draft_2.fasta

    # write a consensus of Ns in case there was nothing mapped at all.
    if [ ! -s consensus_draft.fasta ]
    then
        seq_length=\$(readlength.sh ref.fasta | grep -w '#Bases:' | awk '{print \$2}')
        sequence=\$(printf "%.0sN" \$(seq 1 "\$seq_length") | fold -w 80)
        echo "\$sequence" >> consensus_draft_2.fasta
    else
        seqkit seq -w 80 consensus_draft.fasta | tail -n +2 >> consensus_draft_2.fasta
    fi

    correct_samtools_consensus.py ref.fasta consensus_draft_2.fasta sorted.bam \\
        --call-fract ${params.consensus_min_share} \\
        --min-depth ${params.consensus_min_depth} \\
        --output consensus.fasta
    """
}


process simple_consensus {
    label "general"
    cpus 1
    input:
        tuple val(meta),
            path("ref.fasta"),
            path("sorted.bam"),
            path("sorted.bam.bai")
    output:
        tuple val(meta), path("consensus.${meta.consensus_target}.fasta")
    """
    ${simpleConsensus(meta.alias)}
    mv consensus.fasta consensus.${meta.consensus_target}.fasta
    """
}


process auto_consensus {
    label "medaka"
    cpus { minCpus(2) }
    memory { minRAM(24) }
    input:
        tuple val(meta),
            path("ref.fasta"),
            path("sorted.bam"),
            path("sorted.bam.bai"),
            path("trimmed.fastq")
    output:
        tuple val(meta), path("consensus.${meta.consensus_target}.fasta")
    """
    medaka_cmd="medaka inference sorted.bam consensus.hdf --threads ${task.cpus}"

    ${medakaGetModel()}

    if [[ \$model_choice == "default" ]]
    then
        ${simpleConsensus(meta.alias)}
    else
        medaka_cmd="\$medaka_cmd --model \${model_choice}"
        ${medakaVariantConsensus(meta.alias)}
    fi

    mv consensus.fasta consensus.${meta.consensus_target}.fasta
    """
}


process compute_mapping_stats {
    label "general"
    cpus 1
    input:
        tuple val(meta),
        path('ref.fasta'),
        path('sorted.bam'),
        path('sorted.bam.bai'),
        path('cons.fasta'),
        path('depth.txt')
    output:
        tuple val(meta), path("mapping_stats.csv")
    """
    ref_id=\$(head -n 1 ref.fasta | sed 's/>//' | awk '{print \$1}')
    ref_len=\$(seqtk size ref.fasta | awk '{print \$2}')
    num_mapped_reads=\$(samtools view -F 4 -c sorted.bam)
    n_count=\$(seqtk comp cons.fasta | awk '{print \$9}')
    cons_length=\$(seqtk comp cons.fasta | awk '{print \$2}')
    avg_coverage=\$(awk '{sum+=\$3} END {if (NR > 0) print sum/NR; else print 0}' depth.txt)
    echo "\$ref_id\t\$ref_len\t\$num_mapped_reads\t\$n_count\t\$cons_length\t\$avg_coverage" > mapping_stats.csv
    """
}


process empty_tsv {
    label "general"
    cpus 1
    input:
        val(samplename)
    output:
        tuple val(samplename), path('empty.tsv')
    """
    touch empty.tsv
    """
}


process concat_mapping_stats {
    label "general"
    cpus 1
    input:
        tuple val(samplename), path("collected_stats_*.tsv")
    output:
        tuple val(samplename), path('all_stats.tsv')
    """
    echo "Reference\tReferenceLength\tNumberOfMappedReads\tNCount\tConsensusLength\tAverageCoverage" > all_stats.tsv
    cat collected_stats_*.tsv >> all_stats.tsv
    """
}


process empty_fasta {
    label "general"
    cpus 1
    output:
        path('empty.fasta')
    """
    touch empty.fasta
    """
}


process collect_reference_info {
    label "general"
    cpus 1
    input:
        path('ref_*.fasta')
    output:
        path('reference_info.tsv')
    """
    #!/usr/bin/env python
    import pandas as pd
    import glob
    data = []
    empty = ''
    for fname in glob.glob("ref_*.fasta"):
        with open(fname) as f:
            header = f.readline().lstrip('>').strip()
            if not header:
                continue
            try:
                id, rest = map(str.strip, header.split('|', 1))
                descr, fam, org, orient, seg = map(str.strip, rest.rsplit('|', 4))
                data.append({
                    'Reference': id,
                    'Description': descr,
                    'Family': fam,
                    'Organism': org,
                    'Orientation': orient,
                    'Segment': seg,
                })
            except ValueError:
                data.append({
                    'Reference': header.split('|')[0],
                    'Description': empty,
                    'Family': empty,
                    'Organism': empty,
                    'Orientation': empty,
                    'Segment': empty,
                })
    cols = [
        'Reference',
        'Description',
        'Family',
        'Organism',
        'Orientation',
        'Segment',
    ]
    df = pd.DataFrame(data, columns=cols)
    df.to_csv('reference_info.tsv', sep='\t')
    """
}


process sample_report {
    label "report"
    cpus 1
    input:
        tuple val(samplename),
            path('clean_stats.tsv'),
            path(assembly_stats),
            path(contig_classes),
            path(blast_hits),
            path('mapping_stats.tsv'),
            path(contig_infos),
            path('trimmed_read_stats.tsv'),
            path('cleaned_read_stats.tsv'),
            path('kraken_style_classifications.tsv'),
            path('virus_db_config.yaml'),
            path('contamination_db_config.yaml'),
            path('classification_db_config.yaml'),
            path('reference_info.tsv')
    output:
        tuple val(samplename), path('report.html'), emit: report
        tuple val(samplename), path('stats_reads.tsv'), emit: read_stats
        tuple val(samplename), path('stats_contigs.tsv'), emit: contig_stats
        tuple val(samplename), path('stats_consensus.tsv'), emit: consensus_stats
    """
    export MPLCONFIGDIR=\$(pwd)/.matplotlib
    mkdir -p \$MPLCONFIGDIR

    mergestats_reads.py \\
        --clean-read-stats clean_stats.tsv \\
        --out stats_reads.tsv

    mergestats_contig.py \\
        --contig-classes ${contig_classes.join(" ")} \\
        --blast-hits ${blast_hits.join(" ")} \\
        --contig-info ${contig_infos.join(" ")} \\
        --out stats_contigs.tsv

    mergestats_assembly.py \\
        --stats ${assembly_stats.join(" ")} \\
        --out stats_assembly.tsv

    mergestats_consensus.py \\
        --virus-db-config virus_db_config.yaml \\
        --mapping-stats mapping_stats.tsv \\
        --reference-info reference_info.tsv \\
        --out stats_consensus.tsv

    version_table.py \\
        --virus-db-config virus_db_config.yaml \\
        --contamination-db-config contamination_db_config.yaml \\
        --classification-db-config classification_db_config.yaml \\
        --out db_versions.tsv

    sample_html_report.py \\
        --pipeline-version ${workflow.manifest.version} \\
        --samplename ${samplename} \\
        --virus-db-config virus_db_config.yaml \\
        --reads-stats stats_reads.tsv \\
        --contigs-stats stats_contigs.tsv \\
        --consensus-stats stats_consensus.tsv \\
        --assembly-stats stats_assembly.tsv \\
        --trimmed-read-distribution trimmed_read_stats.tsv \\
        --cleaned-read-distribution cleaned_read_stats.tsv \\
        --db-versions db_versions.tsv \\
        --read-classifications kraken_style_classifications.tsv \\
        --min-classification-frac ${params.centrifuge_plot_min_frac} \\
        --out report.html
    """
}


process get_best_consensus_files {
    label "report"
    cpus 1
    input:
        tuple val(samplename), path('consensus_stats.tsv'), path('cons_*.fasta')
    output:
        tuple val(samplename), path('out')
    """
    get_curated_consensus_genomes.py \
        --consensus-stats consensus_stats.tsv \
        --consensus-files cons_*.fasta \
        --out-dir out
    """
}


process output {
    label "general"
    cpus 1
    publishDir (
        params.out_dir,
        mode: "copy",
        saveAs: {
            dirname ? (fname_out ? "$dirname/$fname_out" : "$dirname/$fname_in") : fname_in
        }
    )
    input:
        tuple path(fname_in), val(dirname), val(fname_out)
    output:
        path(fname_in)
    """
    """
}
