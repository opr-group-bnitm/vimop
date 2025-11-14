# Start an analysis with ViMOP in EPI2ME

## User defined input

If you want to run ViMOP with your own fastqs, you can select your input in the first input menu of the interface. You can select either a fastq or fastq.gz file or a directory that contains subdirectories of the different barcodes that contain the fastq or fastq.gz files.  

For this tutorial, you can also download our [demo data](https://opr.bnitm.de/example_data/vimop-demo/lasv_simulated/barcode01/lasv_l_1000_lasv_s_1000_ms2_800_human_2000_coryne_200.fastq).  

To run the workflow follow these steps:

1. Select the path to the FASTQ file or folder. You can choose either the files themselves or the folder that contains multiple subfolders that contain one or more fastqs.  

![image-9.png](02a_run_vimop_with_epi2me_files/image-9.png)
 
2. Optionally, you can, modify the pipeline parameters in the menu on the left.  

For example, if you want to prioritize one of the curated viruses from our list here you can specify the filters under "Filters and trimming"->"Targets"  

![image-10.png](02a_run_vimop_with_epi2me_files/image-10.png)

Or, under Nextflow configuration you can name your sequencing run.  

![image-11.png](02a_run_vimop_with_epi2me_files/image-11.png)  

Please be aware that you cannot specify the output folder when using epi2me. Read more on that in accessing the output file section.

3. If you have set everything as needed, you can start the analysis by clicking launch workflow.  
You can monitor the progress in the interface.

![image.png](02a_run_vimop_with_epi2me_files/image.png)

## Accessing the output files

### Output files

To access the raw outputfiles, you can find them in the files tab or you can open them in your explorer by going to Options->Open folder.

![image-2.png](02a_run_vimop_with_epi2me_files/image-2.png)

### ViMOP report

ViMOP provides a summarizing report about your analysis.  
Find out more about this report and how to read it in the next tutorial: [03_example_output_interpretation](03_example_output_interpretation.md).

![image-8.png](02a_run_vimop_with_epi2me_files/image-8.png)


