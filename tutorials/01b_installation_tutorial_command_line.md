# Installation of ViMOP with command line
## Table of Contents

- [Prerequisites](#prerequisites)
  - [Docker installation](#docker-installation)
  - [Nextflow installation](#nextflow-installation)
- [Installation and setup of ViMOP](#installation-and-setup-of-vimop)
  - [Run ViMOP demo](#run-vimop)

## Prerequisites
## Docker installation
Please note that you need administrator rights to install Docker. After completing all installation steps however, you will be able to run Docker without admin rights.  
If you cannot install Docker in your system and are on Linux, you can still use ViMOP using the conda or apptainer profile. Just remember to activate the corresponding profiles everytime you execute ViMOP by using the parameters `-profile conda` or `-profile apptainer`.

### MacOS and Windows
For MacOS and Windows we recommend installing Docker Desktop. For this please refer to [01a_installation_tutorial_epi2me.ipynb](01a_installation_tutorial_epi2me.md#macos-and-windows)

### Linux Distributions
You can find a detailed installation manual for your system in the [Docker docs](https://docs.docker.com/engine/install/).  
For most Linux distributions, docker offers a [convenience script](https://docs.docker.com/engine/install/ubuntu/#install-using-the-convenience-script) for installation. You can execute it with

```
curl -fsSL https://get.docker.com -o get-docker.sh
sudo sh ./get-docker.sh
```

After following the steps on that web page, you should make Docker managable as non-root user. For this, follow the steps on this website: [Manage Docker as non-root user](https://docs.docker.com/engine/install/linux-postinstall/), that tells you to run the following commands:  
```
sudo groupadd docker
sudo usermod -aG docker $USER
newgrp docker
```
You may need to restart your Laptop so that Docker shows up as installed in your system.  

### Nextflow installation

The detailed installation steps for the Nextflow installation can be found here: [nextflow documentation](https://nextflow.io/docs/latest/install.html#install-page).  
ViMOP was developed under nextflow version 24.10. However, newer versions should also work.

If you have conda installed, you can create a conda environment and install nextflow in there with:

```
conda create -n nextflow nextflow=24.10
conda activate nextflow
```

## Installation and setup of ViMOP

To run ViMOP you need to install the application and set up the database. Do this in one command with

```
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
