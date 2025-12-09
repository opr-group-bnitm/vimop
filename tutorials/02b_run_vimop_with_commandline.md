# Start an analysis with ViMOP using the command line

Make sure ViMOP is installed as described [here](01b_installation_tutorial_command_line.md).  
If you installed ViMOP without Docker, remember to set the parameter ```-profile``` to ```conda``` or ```apptainer``` at every execution.

## Run ViMOP

To test the pipeline you can get some demo data [here](https://opr.bnitm.de/example_data/lasv_simulated_run.fastq) or download them typing

```bash
wget https://opr.bnitm.de/example_data/vimop-demo.tar.gz
tar -xvzf vimop-demo.tar.gz
```

Subsequently, run ViMOP demo with default settings.

```bash
nextflow run OPR-group-BNITM/vimop --fastq vimop-demo/lasv_simulated --out_dir vimop-demo/output
```

## Options

To show all available options run:

```bash
nextflow run OPR-group-BNITM/vimop --help
```

If you want customize your run and set many options, you may want to use a config file instead of typing all options into your terminal. You can donwload our [nextflow.config](https://github.com/OPR-group-BNITM/vimop/blob/main/nextflow.config) and modify as needed. Then run ViMOP with

```bash
nextflow run OPR-group-BNITM/vimop --fastq /path/to/fastqfiles --out_dir /path/for/your/output -c /path/to/nextflow.config
```
## ViMOP report

ViMOP provides a summarizing report about your analysis.  
Find out more about this report and how to read it in the next tutorial: [03_example_output_interpretation](03_example_output_interpretation.md).

![image-8.png](02a_run_vimop_with_epi2me_files/image-8.png)
