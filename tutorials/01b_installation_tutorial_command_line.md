# Installation of ViMOP with command line
## Table of Contents

- [Prerequisites](#prerequisites)
  - [Docker installation](#docker-installation)
  - [Nextflow installation](#nextflow-installation)
- [Installation and setup of ViMOP](#installation-and-setup-of-vimop)
  - [Run ViMOP demo](#run-vimop)

## Prerequisites

ViMOP uses epi2me, nextflow and docker. All of these dependencies can be installed with and without the usage of the command line depending on the user's preference.

### Docker installation

You can find detailed installation steps for Docker in the [Docker docs](https://docs.docker.com/engine/install/). After following the steps on that web page for your operating system you should make Docker managable as non-root user and make it run on start up. For this, follow the steps on this website: [Post-install instructions](https://docs.docker.com/engine/install/linux-postinstall/).

### Nextflow installation

The detailed installation steps for the Nextflow installation can be found here: [nextflow documentation](https://nextflow.io/docs/latest/install.html#install-page).  

If you have conda installed, you can create a conda environment and install nextflow in there with:

```bash
conda create -n nextflow nextflow
conda activate nextflow
```

## Installation and setup of ViMOP

To run ViMOP you need to install the application and set up the database. Do this in one command with

```bash
nextflow run OPR-group-BNITM/vimop --download_db_all
```

The nextflow run command will automatically install the pipeline while the `--download_db_all` will setup the database. If you want to do this separately, you could run (`nextflow pull OPR-group-BNITM/vimop`) first.

**Resume option**
If the pipeline fails during the process (which may happen due to instable network access), use the `-resume` option to continue your download without having to restart again. This works also for analysis run.

**modular database download**
You can also separate the download of the reference data into three parts by running the pipeline three times using the options `--download_db_virus`, `--download_db_contamination` and `--download_db_centrifuge` in separate runs. We would recommend this especially with bad network connection. 

**Updata an exising database**
If you want to replace an existing database with our latest version, add the option `--download_db_update_existing`.

### Run Vimop

How to run ViMOP on a test data set is described [here](02b_run_vimop_with_commandline.md).
