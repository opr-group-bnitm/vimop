# Start an analysis with ViMOP using the command line

Make sure ViMOP, its data base and prerequisites are downloaded and installed as described in [01b_installation_tutorial_command_line.md](01b_installation_tutorial_command_line.md).  

## Run ViMOP

To test the pipeline you can get some demo data [here](https://opr.bnitm.de/example_data/lasv_simulated_run.fastq) or download them by typing

```bash
wget https://opr.bnitm.de/example_data/vimop-demo.tar.gz
tar -xvzf vimop-demo.tar.gz
```

Subsequently, run ViMOP demo with default settings:

```bash
nextflow run opr-group-bnitm/vimop --fastq vimop-demo/lasv_simulated --out_dir vimop-demo/output
```

## Options

To show all available options run:

```bash
nextflow run opr-group-bnitm/vimop --help
```

If you want customize your run and set many options, you may want to use a config file instead of typing all options into your terminal. You can donwload our [nextflow.config](https://github.com/opr-group-bnitm/vimop/blob/main/nextflow.config) and modify as needed. Then run ViMOP with

```bash
nextflow run opr-group-bnitm/vimop --fastq /path/to/fastqfiles --out_dir /path/for/your/output -c /path/to/nextflow.config
```
## ViMOP report

ViMOP provides a summarizing report about your analysis run.  
Find out more about this report and how to read it in the next tutorial: [03_example_output_interpretation](03_example_output_interpretation.md).

![image-8.png](02a_run_vimop_with_epi2me_files/image-8.png)
