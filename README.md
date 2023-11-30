# A CGAT-core pipeline to benchmark tRNA-Seq quantification methods


This pipeline was used for the manuscript: *Benchmarking tRNA-Seq quantification approaches by realistic tRNA-Seq data simulation*. Smith, T., et al.

The pipeline was build using CGAT-Core. CGAT-core is a workflow management system to build scalable data analysis pipelines. CGAT-core includes libraries and helper functions to enable researchers to quickly design and build computational workflows for the analysis of large-scale data-analysis.

<a href="https://github.com/cgat-developers/cgat-core">
  <img src="https://github.com/cgat-developers/cgat-core/blob/master/docs/img/CGAT_logo.png" alt="CGAT-core" width="100">
</a>

[CGAT-core Documentation](https://cgat-core.readthedocs.io/en/latest/ "CGAT-core read the docs")

The pipeline steps include:

- Generates error profiles
- Compares errors with known modifications
- Simulate tRNA-Seq samples
- Align reads
- Compares simulation ground truth and observed alignments
- Tallies tRNA counts
- Compares quantification estimates with ground truths
- Quantifies from real tRNA-Seq data 

### Installation
The same installation instructions apply to run the pipeline locally and on the HPC.
We're using mamba below. See https://cgat-core.readthedocs.io/en/latest/getting_started/Installation.html
for further CGAT-core installation options if needed.


See https://github.com/mamba-org/mamba for instructions on how to install mamba.
Can alternatively only use conda, but mamba is quicker for installation. Even with mamba,
you may find the following steps take a few minutes.

1. Create a new conda environment with the required packages
	```bash
	mamba create --name trna -f <PATH TO THIS REPOSITORY/conda_envs/trna.environment.yml>
	```

2. Create a separate dedicated environment for `mimseq`, called 'mimseq', following the instructions [here](https://github.com/nedialkova-lab/mim-tRNAseq).

3. Activate the trna environment: `conda activate trna`

## Prepare the pipeline direcory

### Download raw data
Create a new directory to run the pipeline and create a subdirectory called `input.dir`. Follow the instructions in the README files [here](download_input/raw_data/) to download the _Homo sapiens_ tRNA-Seq raw data into `input.dir`. Repeat this process for the _Mus musculus_ data in a separate to directory. The pipeline must be run separately for the two species. 

### Copy the reference files and MODOMICS files

The reference files can be found [here](download_input/reference_files/). Copy all files for the species

The MODOMICS files can be found [here](download_input/modomics/). Copy the `index.csv` and the species-specific `sequences.json` file.


### Run the pipeline (locally)

Make sure you are in a directory containing just the raw data in the `input.dir` directory and the reference files and modomics files at the top of the directory. 

The options we use to run the pipeline are:

- `-p1`: One parallel task (Any more may require too much RAM or CPUs to run locally)
- `--checksums = 2`: Keep track of job completion
- `-v10`: Verbose output
- `--local`: Don't submit jobs to a HPC, run them locally.
- `make full`: Run to the task full


```bash
python <PATH TO THIS REPOSITORY/pipelines/pipeline_compare_trnaseq_quant.py> --checksums=1 -p1 -v10 make full  --local
```

The final output of the pipeline are in the `final_results.dir` directory.

If you'd like to run the pipeline faster, change to `-p<NUMBER OF TASKS>`, replacing `NUMBER OF TASKS` with the number of parallel tasks you would like to run. However, be aware CPU/memory usage will quickly exceed your desktop/laptop.... which brings us to...

### Run the pipeline (on the cambridge HPC)

First, log into the cambridge HPC (https://docs.hpc.cam.ac.uk/hpc/) and navigate to a suitable RDS directory. Then follow the same directions above to download the input data and install CGAT-core and the pipeline dependencies with mamba/conda.

Before we run the pipeline, we need to generate a `.cgat.yml` config file in the home directory so that CGAT-core knows how to interact with the HPC workload manager and submit jobs etc. The UoC HPC uses SLURM as the workload manager. There is a template [.cgat.file](https://github.com/MRCToxBioinformatics/Pipeline_examples/blob/main/CGATCore/.cgat.yml), which you can save to your home directory and update to provide the required details. Further details about configuring CGAT-core pipelines to run on clusters are here: https://cgat-core.readthedocs.io/en/latest/getting_started/Cluster_config.html.

After that, it's as simple as running the same commmand as previously, but without the `--local` and allowing more parallel processes.  Below, we allow a maximum of 10 tasks to be run in parallel. Note that this is different to the number of threads used in each individual task, which is parameterised within the pipeline code and/or config file.

We can also increase the number of threads used for each of the hisat tasks to speed them up. This is controlled via the `pipeline.yml` config file. By default, the pipeline uses the config file in the same directory as the pipeline source code. We can also include a config file in the working directory to overwrite the parameters.

> &#x26a0;&#xfe0f; **Check you are in the appropriate working directory in your RDS before you run this command. You should have an input.dir subdirectory with the neccessary input files.**

```bash
nohup python <PATH TO THIS REPOSITORY/pipelines/pipeline_compare_trnaseq_quant.py>  --checksums=1 -p100 -v10 make full &
```

Note that we also include the command `nohup` at the start and `&` at the end of our command. `&` puts the command into the background and `nohup` ensures it will keep running even if we log out or lose connection with the HPC.

Your submitted jobs may not start immediately, but when they do, they will be much quicker and you can run many parallel tasks.  For the example data here, you are unlikely to see any significant run time benefit running on the HPC. However, if you want to analyse many samples, using the HPC becomes essential! This will also allow you to run resource-hungry jobs beyond your desktop/laptop specifications.

### Run the pipeline (on another HPC)
This has not yet been tested. Please get in touch if you'd like help to do so.

## &#x26a0;&#xfe0f; Troubleshooting
On the HPC, you may get the following errors:

------------

`conda Can't locate local/lib.pm in @INC (you may need to install the local::lib module)`

This can be resolved with:
`mamba install -c bioconda perl-local-lib`

------------

`gtfToGenePred: error while loading shared libraries: libssl.so.1.0.0: cannot open shared object file: No such file or directory`

Try the following to resolve
`mamba install -c bioconda samtools=1.13`
`mamba update ucsc-gtftogenepred`

------------



