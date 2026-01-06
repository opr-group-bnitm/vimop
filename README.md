# ViMOP

[<img align="left" src="ViMOP_logo.png" width="150"/>](ViMOP_logo.png)

ViMOP is a pipeline for detecting known virus species in untargeted Oxford Nanopore sequencing data and reconstructing their genomes.

Developed by the [Outbreak Preparedness and Response team](https://www.bnitm.de/forschung/forschungsgruppen/pathogen/abt-virologie/laborgruppe-duraffour-pahlmann/team) at the Bernhard Nocht Institute for Tropical Medicine, ViMOP is used to analyze nanopore reads from clinical samples of viruses such as Lassa, Ebola, or Dengue at various sequencing sites.

If you have questions, suggestions, or would like to contribute — or if you require a specific setup (e.g. for licensing) — please feel free to contact us.

## Quickstart

You can run and install the pipeline from the command line using [nextflow](https://www.nextflow.io/) or from the [EPI2ME desktop](https://nanoporetech.com/software/other/epi2me-desktop-application) application from ONT.

Here are tutorials for setup, running the pipeline and interpretation of the output:

- [installation with EPI2ME Desktop](https://github.com/OPR-group-BNITM/vimop/blob/main/tutorials/01a_installation_tutorial_epi2me.md)
- [run an analysis with EPI2ME Desktop](https://github.com/OPR-group-BNITM/vimop/blob/main/tutorials/02a_run_vimop_with_epi2me.md)
- [installation from command line](https://github.com/OPR-group-BNITM/vimop/blob/main/tutorials/01b_installation_tutorial_command_line.md)
- [run an analysis from the command line](https://github.com/OPR-group-BNITM/vimop/blob/main/tutorials/02b_run_vimop_with_commandline.md)
- [output interpretation](https://github.com/OPR-group-BNITM/vimop/blob/main/tutorials/03_example_output_interpretation.md)

## Purpose and limitations

The main purpose of this pipeline is the assembly of virus genomes from human clincical or animal samples.

The pipeline automatically finds well fitting virus genomes and uses them as references to build reference based consensus genomes.
This works well for small and medium size RNA viruses such as Lassa, Dengue, Ebola and many others.
However, for large DNA viruses with extensive repetitive regions (e.g., mpox), assemblies may contain inaccuracies and should be regarded as draft genomes.
In any case, we recommend carefully reviewing your output (e.g. the alignment .bam files).

### Reference database

We have created a reference database with our favourite viruses.
However, you can also easily create your own.
For information on databases read further [down](#database).
If you need assistance for setting up a reference dataset, please contact us.

## Hardware requirements

This pipeline runs best on a powerful laptop or PC.
We recommend at least 30 GB RAM and 16 CPUs.
Depending on your dataset, the pipeline may also work on lower resources.
You can change the parameters `--min_cpus`, `--min_ram`, `--min_disk_space_work_gb` and `--min_disk_space_out_gb` and the pipeline will run with less.
However, this may or may not work, and it may take much longer, as some tools like Canu need quite some resources. 

## Software dependencies

ViMOP uses
- [Nextflow](https://www.nextflow.io/)
- [Docker](https://www.docker.com/)
- [EPI2ME desktop](https://nanoporetech.com/software/other/epi2me-desktop-application) (optional)

For installation and setup see our command line or EPI2ME desktop [tutorials](#quickstart).
If you do not want to use docker, you can also use a conda or apptainer profile as described in the next section.

### Alternative profiles for command line usage: Conda and Apptainer

By default ViMOP uses docker (for setup see the installation tutorials).
This should run on all operating systems and so far has been tested on the following operating systems:
- MacOS 15.6.1 (Intel core i9); Nextflow 25.10.2; Docker 25.0.2
- MacOS 15.5 (Intel core i9); Nextflow 25.10.2; Docker 24.0.5
- Ubuntu 24.04.1; Nextflow 25.10.2; Docker 28.5.1
- Ubuntu 22.04.5; Nextflow 25.10.2; Docker 29.1.2

If you prefer to not use docker and you are using Linux, there are two alternative profiles implemented in ViMOP.

The **conda** profile is activated using the option `-profile conda`.
This has been tested on
- Ubuntu 24.04.1; Nextflow 25.10.2; conda 25.11.0
- Ubuntu 22.04.5; Nextflow 25.10.2; conda 25.7.0

You can also use **mamba** typing `-profile conda,mamba`.
However, some versions of mamba and nextflow may not work together.
We succesfully ran this with
- Ubuntu 24.04.1; Nextflow 25.10.2; conda 25.11.0; mamba 2.4.0
- Ubuntu 22.04.5; Nextflow 25.10.2; conda 25.7.0; mamba 2.3.1

The **apptainer** profile is activated using `-profile apptainer`.
It has been tested on
- Ubuntu 24.04.1; Nextflow 25.10.2; apptainer 1.4.5
- Ubuntu 22.04.5; Nextflow 25.10.2; apptainer 1.4.5

## Workflow
![vimop flowchart](ViMOP_flowchart.png)
In the following, the most important steps of the pipelines are explained with the respective options to set.

### Input

A directory with fastq files is passed to `--fastq`.
The directory can contain subdirectories called barcode01, barcode02, ... 
In this case, the different barcodes are treated as different samples with separate output produced for each.

### Read trimming

At the beginning the pipeline trims the ends of the reads to remove adapter sequences and primers.
`--trim_length` sets the number of bases trimmed from both ends.

### Taxonomic classification and removal of non-viral reads

Centrifuge is used to classify the reads.
This helps to get an overview of how your sample is composed. A Kronaplot of the read classification will appear in the report and optionally you can choose to generate a separate HTML with the Krona plot with --output_krona_plot true.
Centrifuge also classifies the contigs (see later) to get a rough idea about contigs that do not match any entry in the reference database or only partially. 
Use `--centrifuge_do_classify false` to deactivate all centrifuge classifications and save time.

Additionally, the centrifuge classifications is used to remove reads (not contigs) that are non-viral to a given degree of confidence.
Deactivate this with `--centrifuge_do_filter false` and use `--centrifuge_filter_min_score 150` to set the minimum level of confidence (e.g. 150 in our example).
The pipeline will remove all reads that are classified to anything that is not viral with at least this score.
Use this if you have a lot of different contaminations that you want to remove such as bacterial reads from a fecal sample.
But beware that centrifuge is limited in its accuracy and you may also remove some false positives.

### Host and contaminant removal

Host and contaminant reads are removed by mapping them against reference sequences and extracting those that map.
The database config file (contamination.yaml) defines the reference sets and the key values assigned to them.
Use the option `--contamination_filters "reagent,mouse"` for example to remove mouse and reagent reads (the default is human reads).
The following filters are included in our database

| Filter        | Description                                             |
|---------------|---------------------------------------------------------|
| reagent       | Reagent associated sequences                            |
| human_rna     | Human transcriptome                                     |
| human_dna     | Human genome                                            |
| mouse         | Mus musculus genome                                     |
| mastomys      | Mastomys natalensis genome                              |
| aedes_aegypti | Aedes aegypti genome                                    |

### Target virus read enrichment

One can also filter for specific species or groups of viruses.
Only reads that map to the given targets are then used in the following assembly step.
An arbitrary number of filters can be used.
The filters themselves are part of the reference database and the respective names defined in the database configs.
This command `--targets "MARV,EBOV,FILO"` would activate filters for Marburg virus, Ebola and the Filo-virus family.

In our default database there are filters for individual virus species and for virus families.
They are listed in the following:

| Virus                                           | Abbreviation | TaxId     |
| ----------------------------------------------- | ------------ | --------- |
| Emesvirus zinderi                               | MS2          | 329852    |
| Lentivirus humimdef1                            | HIV1         | 3418650   |
| Lentivirus humimdef2                            | HIV2         | 3418651   |
| Mammarenavirus lassense                         | LASV         | 3052310   |
| Mammarenavirus choriomeningitidis               | LCMV         | 305230    |
| Mammarenavirus juninense                        | JUNV         | 2169991   |
| Orthoebolavirus                                 | EBOV         | 3044781   |
| Orthoflavivirus denguei                         | DENV         | 3052464   |
| Orthoflavivirus zikaense                        | ZIKA         | 3048459   |
| Orthomarburgvirus                               | MARV         | 3044783   |
| Orthonairovirus hazaraense                      | HAZV         | 3052519   |
| Severe acute respiratory syndrome coronavirus 2 | COVID        | 2697049   |
| West nile virus                                 | WNV          | 3048448   |
| Yellow fever virus                              | YFV          | 3046277   |

| Virus family            | Abbreviation | TaxId     |
| ----------------------- | ------------ | --------- |
| All                     | ALL          | 10239     |
| Arenaviridae            | ARENA        | 11617     |
| Filoviridae             | FILO         | 11266     |
| Hantaviridae            | HANTA        | 1980413   |
| Nairoviridae            | NAIRO        | 1980415   |

### De novo assembly

An assembly is run for each of the filtered read sets (explained in the previous section) and for the read set with no target filter.
To disable running the no target filter set `--assemble_notarget false`.

After the initial assembly, an iterative re-assembly procedure is run (unless deactivated with `--reassemble_max_iter 0`).
The purpose is to also find virus segments, that are present in very low concentration.
All reads are mapped to the contigs and those that map are removed.
The remaining reads get assembled once more.
This procedure is repeated until a maximum number of cycles is reached, no reads are left after filtering or no contigs are produced by the canu assembly.

For read sets that were filtered to a target (see previous section), there is a special procedure if no contigs were assembled.
In this case, the longest X reads (set with `--nocontigs_max_reads_precluster`) are passed to a clustering by cd-hit-est.
Of these clusters, the longest Y reads (set with `--nocontigs_nreads`) are chosen instead of contigs for the following target search.
Canu-corrected reads are used for this if available, else the raw reads.

A number of parameters for the canu assembler can be set, that determine how many reads are used and corrected for assembly. See options.

### Reference identification

Each contig is used to for a BLAST search in the virus reference database.
The highest scoring hit is then used as a reference genome.

### Reference-guided assembly

Reads are mapped against the reference genome.
The mapping parameters can be changed (see options).
There are multiple options to generate the consensus.
The default option is to choose automatically from samtools consensus and medaka.
If medaka finds a model for you data, medaka is used, else samtools.
Samtools consensus takes the most abundant base and masks all positions with too little coverage or unclear signals below the given threshold.
You can also directly choose medaka (medaka) or samtools (simple).

For medaka you can pass a model name for the option `medaka_consensus_model`.
Choosing auto will let medaka choose the model, which is the default.
However, if medaka does not find a model fitting your data, medaka (if explicitly chosen) will use the medaka default model, which may not be optiomal.

## Options

To get a list of all available options and their default parameters type `nextflow main.nf --help`.

## Output and report

The output is structured like this:

```
output
├── nf-report.html
├── params.json
└── barcode01
    ├── report_barcode01.html
    ├── classification
    │   ├── classification.html
    │   ├── classification.tsv
    │   ├── classification_kraken.tsv
    │   └── classification_report.tsv
    ├── consensus
    │   ├── AB627954.consensus.fasta
    │   ├── AB627954.depth.txt
    │   ├── AB627954.reads.bam
    │   ├── AB627954.reads.bam.bai
    │   ├── AB627954.reference.fasta
    │   ├── AB627954.variants.vcf.gz
    │   ├── CS272305.consensus.fasta
    │   ├── CS272305.depth.txt
    │   ├── CS272305.reads.bam
    │   ├── CS272305.reads.bam.bai
    │   ├── CS272305.reference.fasta
    │   └── CS272305.variants.vcf.gz
    ├── tables
    │   ├── consensus.tsv
    │   ├── contigs.tsv
    │   └── reads.tsv
    ├── assembly
    │   ├── no-target.contigs.fasta
    │   └── re-assembly.contigs.fasta
    └── selected_consensus
        ├── LCMV_L.fasta
        └── LCMV_S.fasta
```

nf-report.html contains technical information about the run and ressource usage.
For each sample there is a directory with results (here barcode01).
The summary of the results is found in report_samplename.html.
Tables contains the information from the html report in .tsv files.
The classification directory contains all files for the centrifuge read classification.
The directory consensus contains consensus genome sequences, alignment files and variants as well as the chosen reference.
The chosen selected genomes for curated virus species are listed in selected_consensus with a separate file for each segment.
Contigs can be found in the fasta files in the assembly directory.

## Database

ViMOP relies on a reference database structure.
In the following, the structure of the database is described.
To use ViMOP to create your own custom database see further [below](#custom-database-creation).

### Database structure

The ViMOP data base is usually placed in your home directory in a folder called `ViMOP_DB`.
I has the following structure:

```
ViMOP_DB/
├── centrifuge/
├── contaminants/
└── virus/
```

The three subdirectories contain files for centrifuge classification, contaminants/host read removal and the virus reference sequences.
Each directory contains a file with a yaml file with the same name prefix (e.g. centrifuge.yaml, contaminants.yaml, virus.yaml).
The configs hold the relevant information about the database parts as well as an entry 'version' with a version number and an entry description with a brief 'description'.

The three database parts are briefly described in the following.
ViMOP also include a module to create your own custom data base.
It is described after the general description of the pipeline parts.

#### centrifuge

The centrifuge config looks like this

```yaml
version: 1.0
description: "Refseq reference genomes plus genbank virus sequences"
index_name: all
files:
- all.1.cf
- all.2.cf
- all.3.cf
- all.4.cf
virus_taxid_file: virus_taxids.txt
```

The index name has to be the prefix of the centrifuge index which are the files listed under files.
These files need to be in the centrifuge DB directory.
The virus_taxid_file contains all virus taxids.
This information is important for the centrifuge based filtering.
In addition to unclassified reads, reads classified to these Tax-IDs will be kept since they are considered to be virus reads.
Version and description are for display in the report.

#### contaminants

This directory holds files with sequences of host or reagents.
The respective config file looks like this

```yaml
filters:
  reagent: "reagent-db.fasta.gz"
  human_rna: "GCF_000001405.39_GRCh38.p13_rna.fna.gz"
  human_dna: "GCF_000001405.39_GRCh38.p13_genomic.fna.gz"
  mouse: "GCA_000001635.8_GRCm38.p6_genomic.fna.gz"
  mastomys: "GCF_008632895.1_UCSF_Mcou_1_genomic.fna.gz"
version: 1.0
description: "Human (GRCh38), mouse (8_GRCm38), mastomys and contaminant filter set"
```

The keys are used to choose the filters using the command `--contamination_filters`.

#### virus

The virus database contains virus reference sequences.
It consists of a config file, a set of sequence files and a blast database.
Let's have a look at a small example.
Our directory content could look like this

```
ViMOP_BUCKET/
├── db.yaml
├── ALL.fasta
├── EBOV.fasta
├── LASV.fasta
├── FILO.fasta
└── blast_db/
    ├── ALL.nhr
    ├── ALL.ndb
    ├── ALL.nin
    ├── ALL.njs
    ├── ALL.nog
    ├── ALL.not
    ├── ALL.nos
    ├── ALL.nsq
    ├── ALL.ntf
    └── ALL.nto
```

And the corresponding config could look like this:

```
all:
  blast_db: blast_db
  blast_prefix: ALL
  fasta: ALL.fasta
curated:
  LASV:
    fasta: LASV.fasta
    name: Mammarenavirus lassense
    organisms:
    - Lassa virus GA391
    - Lassa virus Josiah
    - Mammarenavirus lassaense
    segments:
    - S
    - L
  EBOV:
    fasta: EBOV.fasta
    name: Orthoebolavirus
    organisms:
    - Orthoebolavirus bombaliense
    - Orthoebolavirus bundibugyoense
    - Orthoebolavirus restonense
    - Orthoebolavirus sudanense
    - Orthoebolavirus taiense
    - Orthoebolavirus zairense
    segments:
    - Unsegmented
filters:
  FILO: FILO.fasta
version: 1.0
description: "A database that has Lassa and Ebolavirus curated."
params.fasta_sequences: /data/home/nils.petersen/dev/VirusDatasetCuration/workflow/testset/testset.fasta
params.taxa_config: /data/home/nils.petersen/dev/VirusDatasetCuration/workflow/testset/test_groups_refs_and_organisms.yaml
params.filter_max_n_share: 0.01
params.filter_min_relative_length: 0.8
params.cdhit_threshold: 0.98
```

ALL.fasta contains all sequences.
EBOV.fasta contains ebola virus sequences, LASV Lassa virus and FILO filo virus sequences.
These are the files that are used as mapping filters.
The file ALL.fasta is also used to create the blast database.
Our blast database is build with the blast version used in the pipeline (see `images/general/general.yaml` in this repository).
The curated viruses will be shown in their own sections in the report and have segment wise one fasta file, where the consensus sequence with the highest recovery is chosen.
All other virus targets are added in ALL and 

The headers in the fasta files have to be formatted like in this example

```
>OQ791477.1 |MAG: Lassa mammarenavirus isolate c0212 segment L genomic sequence|Arenaviridae|Mammarenavirus lassaense|forward|L
CAAAATGGGCAACAAGCAAGCCAAGTCAACCAAAGTCGATGAACAACATAGAGCTCAT...
```

Separated with a "|" we have
- genbank ID
- description
- family
- species name
- orientation of the sequence with respect to the original database entry. We re-oriented sequences so that all sequences of a curated dataset have the same orientation. However, this can also simply be set to "Unknown".
- the segment name. Set to "Unknown" for non-curated sequences. For curated sets (e.g. in our example LASV and EBOV) this needs to be assigned. If there is only one segment, use "Unsegmented". The segments also need to be listed in the config file.

### Custom database creation

ViMOP provides a module to create your own custom data base.
You can create a whole data base or only parts of it (`virus`, `contaminants` and/or `centrifuge`).
To run the module in EPI2ME desktop, select the section `Create custom data base` and tick the box `Build a custom data base`.
Then choose paths to you input files as described below.
To run the data base setup from the command line use the option `--custom_db_do_build`.
Use the option `--custom_db_outpath` to define the directory where the data base will be written to.
In the following, the command line arguments are listed but not the EPI2ME fields, as these are self-explanatory.

You can also find an example command on how to create a custom data base in `test/build_custom_db.sh` and example input files in `test/data/custom_db_test`.

#### Build contaminants/host database

Building a host data base takes a directory with fasta files.
A filter will be created for each fasta file.
The contaminant/host filters will have the names in the files (e.g. `human_dna` if a file `humand_dna.fasta.gz` exists in the directory). Fasta files can be compressed (`fasta.gz`) or not (`.fasta`).
Pass this directory via `--custom_db_contaminants_input_path`.
Use `--custom_db_contaminants_version` and `--custom_db_contaminants_description` to add version and description.

#### Build virus database

To create the virus data base pass a fasta file with genomes (`--custom_db_virus_fasta`) and a yaml file with additional configuration information (`--custom_db_virus_yaml`).

The format of the fasta headers needs to be `>SEQID |DESCRIPTION|FAMILY|SPECIES|ORIENTATION|SEGMENT`.
Orientation can be forward, reverse or unknown and is with respect to the data base entry (e.g. if the orentation was swapped for a genbank entry you can add reverse).
It is also ok to always set this column to unknown.
Segment should be Unsegmented or the segment name or Unknown.
Unknown should only be used for non-curated species, whereas highlighted species should always have segment information.
Species must match those used in the config (see below).

An example of a genome entry: `>KU174142.1 |Mutant Zaire ebolavirus isolate Ebola virus/H.sapiens-rec/COD/1976/Yambuku-Mayinga-eGFP-SUDV_GP, complete genome|Filoviridae|Orthoebolavirus zairense|forward|Unsegmented`

The yaml file looks like this

```yaml
curated:
  LASV:
    name: Mammarenavirus lassense
    organisms:
    - Mammarenavirus lassaense
  EBOV:
    name: Orthoebolavirus
    organisms:
    - Orthoebolavirus zairense
    - Orthoebolavirus taiense
    - Orthoebolavirus sudanense
    - Orthoebolavirus restonense
    - Orthoebolavirus bundibugyoense
    - Orthoebolavirus bombaliense
    - Orthoebolavirus
    - Bombali virus
  MARV:
    name: Orthomarburgvirus
    organisms:
    - Orthomarburgvirus marburgense
    - Orthomarburgvirus
filters:
  FILO:
    organisms:
    - Orthomarburgvirus marburgense
    - Orthoebolavirus zairense
    - Orthoebolavirus taiense
    - Orthoebolavirus sudanense
    - Orthoebolavirus restonense
    - Orthoebolavirus bundibugyoense
    - Orthoebolavirus bombaliense
    - Orthomarburgvirus
    - Orthoebolavirus
```

The first part determines the curated virus sets.
Here, each species name is only allowed once.
The species names must match those in the headers of the fasta file.
Results for curated sets will get their own dedicated section in the report and sorted by segment.
Additional filters (e.g. to filter for virus families) can be created using again the species names.
Here, species names can appear in arbitrarily many filters.

Use `--custom_db_virus_version` and `--custom_db_virus_description` to add version and description.

#### Build centrifuge index

To build a centrifuge index pass a single fasta file.
The header must be formatted as `>SEQID |KINGDOM|FAMILY|SPECIES`.
`KINGDOM` must be one of Eukaryota, Archaea, Bacteria or Viruses.

An example: `>KU174142.1  |Viruses|Filoviridae|Orthoebolavirus zairense`

Use `--custom_db_centrifuge_version` and `--custom_db_centrifuge_description` to add version and description.

## Citation

If you are using ViMOP please cite us:

Petersen NP, Le M, Renevey A, Emua E, Ryter S, Annibaldis G, Camara J, Boumbaly S, Erameh C, Laske T, Baumbach J, Lemey P, Günther S, Duraffour S, Kafetzopoulou LE. ViMOP: A user-friendly and field-applicable pipeline for untargeted viral genome nanopore sequencing. *Bioinformatics*. 2025; btaf687. [10.1093/bioinformatics/btaf687](https://doi.org/10.1093/bioinformatics/btaf687)

This repository is linked to [zenodo](https://doi.org/10.5281/zenodo.15592229), where you can find a DOI for the version you are using.

## Acknowledgements

This product includes software developed by Oxford Nanopore Technologies Plc.
