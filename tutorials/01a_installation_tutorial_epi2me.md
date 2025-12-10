# ViMOP with EPI2ME

This tutorial describes how to install Docker, Nextflow, EPI2ME and then ViMOP via the graphical user interface of EPI2ME.

# Table of Contents

- [ViMOP with EPI2ME](#vimop-with-epi2me)
  - [Prerequisites](#prerequisites)
    - [Docker installation](#docker-installation)
    - [EPI2ME installation](#epi2me-installation)
    - [Nextflow installation via EPI2ME](#nextflow-installation-via-epi2me)
  - [Installation and setup of ViMOP](#installation-and-setup-of-vimop)
    - [Install workflow](#install-workflow)
    - [Set up the data base](#set-up-the-data-base)
    - [Demo run](#demo-run)
  - [Accessing the source code via the command line](#accessing-the-source-code-via-the-command-line)
  
# Prerequisites

## Docker installation
Please note that you need administrator rights to install Docker. After completing all installation steps however, you should be able to run Docker with your regular user account.

### MacOS and Windows
If you do not have the Docker engine already installed in your system, you can install it bundled together with Docker Desktop. In MacOS and Windows the installation of Docker Desktop works completely without the command line. Download the Docker Desktop installer according to your operating system on the bottom of this [Docker docs page](https://docs.docker.com/desktop/) and follow the instructions after opening the installer. Restart your system after you are done with the installation. Remember to start the Docker Desktop App before you procede with the tutorial and whenever you want to use ViMOP.

### Linux distributions
For Linux we recommend the installation of the Docker engine. This requires a bit of command line usage by following the Docker installation steps in [01b_installation_tutorial_command_line.md](01b_installation_tutorial_command_line.md#linux-distributions).

## EPI2ME installation
To use ViMOP in the graphical user interface you need to install Oxford Nanopore's EPI2ME Desktop application. Download the installer according to your operating system from the [EPI2ME website](https://epi2me.nanoporetech.com/downloads/). At the time of this tutorial the most recent version of EPI2ME Desktop is version 5.3. Follow the installation instructions after opening the downloaded file, for Ubuntu/Debian you need to open the file with the software center. EPI2ME also provides an [installation guide](https://epi2me.nanoporetech.com/epi2me-docs/installation/) for all operating systems.  

### Start up EPI2ME
After starting the application, you will be prompted to sign in. However, you can just continue as guest. Sometimes that option is hidden, but you can find it by clicking on the three dots at the bottom.  

![continue_as_guest](01a_installation_tutorial_epi2me_files/continue_as_guest.png)

## Nextflow installation via EPI2ME
Via EPI2ME Desktop you can install nextflow via the EPI2ME interface by going to Settings->Local->System setup->Nextflow->setup.  
This setup will also make sure you have the correct java version installed.

# Installation and setup of ViMOP
As soon as you have EPI2ME, Nextflow and Docker installed in your system you are good to go to! Restart your computer, make sure Docker is running and open EPI2ME.

## Install workflow 
Follow these steps to import ViMOP into EPI2ME
1. Import ViMOP from GitHub  
  
![import_GH](01a_installation_tutorial_epi2me_files/image.png) 
 
2. Paste the URL to our GitHub repository into the interface ```https://github.com/opr-group-bnitm/vimop``` and click Download  
  
![image-2.png](01a_installation_tutorial_epi2me_files/image-2.png)
 
## Set up the data base
Before analyzing your first sequencing run, you need to install the data base. For this follow these steps:  
  
1. Open workflow  
  
![image-3.png](01a_installation_tutorial_epi2me_files/image-3.png)
 
2. Launch the workflow  
  
![image-4.png](01a_installation_tutorial_epi2me_files/image-4.png)
 
3. Open the setup menu  
  
![image-5.png](01a_installation_tutorial_epi2me_files/image-5.png)
 
4. Select to download all databases. If you want to update an already existing database you have downloaded from us, you need to select to overwrite the existing databases as well.  
  
![image-6.png](01a_installation_tutorial_epi2me_files/image-6.png)

In this same menu you could also choose to only download or update each section of the database separately (i.e. centrifuge index, host and contaminant database and virus reference genome database). Which is especially recommended in settings with bad network connection.
 
5. Launch the workflow.  
  
![image-7.png](01a_installation_tutorial_epi2me_files/image-7.png)
 
The download will now take a while depending on your network connection.

## Demo run
To test the functionality of ViMOP you can select to run a demo run with a simulated Lassa virus here:  
  
![image-8.png](01a_installation_tutorial_epi2me_files/image-8.png)  

To run ViMOP with your own data checkout [02a_run_vimop_with_epi2me](02a_run_vimop_with_epi2me.md).

# Accessing the source code via the command line

After installing ViMOP with EPI2ME you can also execute it via the command line by accessing the source code directly in EPI2ME's workflow folder.  
Open the terminal and change the directory to:

```
cd epi2melabs/workflows/opr-group-bnitm/vimop
```

To run the tool from here via the command line, install nextflow as described in [01b_installation_tutorial_command_line.md](01b_installation_tutorial_command_line.md#nextflow-installation) and then refer to our [command line tutorial](https://github.com/OPR-group-BNITM/vimop/blob/ml_tutorials/tutorials/02b_run_vimop_with_commandline.md).