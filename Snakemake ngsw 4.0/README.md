[![Snakemake](https://img.shields.io/badge/snakemake-≥5.6.0-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)


Greetings,

It is assumed that you have your genome.fa, .fa.ann, .fa.amb, .fa.bwt, .fa.fai, .fa.pac, .fa.sai, .dict. You might need bwa Index and Samtools faidx. Specify the paths in the config file. Also, you have to create a "samples" directory inside the main directory in which you are working, put there your samples and rename them like: a_1.fastq/a_1.fastq ... z_1.fastq/z_2.fastq. In the Snakefile (this is going to change very soon moving to config file), search at the top for the list "SAMPLES", and add as many letters as many samples you have. 

If you have Conda installed on your system, you are pretty ready.

Now, create a working directory "snakemake-yourjob" and the minimal Snakemake environment from the provided file, typing:

		   	conda env create --name Snakemake-yourjob -f environment.yaml

Activate the environment by typing

			conda activate Snakemake-yourjob
	   
With wrappers, all the others packages needed will be provided and installed in your conda environment. 
Here you get a run.py python script, which handles the run of the Snakefile for you. So, type:

	   		./run.py  -h  (for help)

You are able to limit the workflow to your necessities. Choose the workflow between the listed options and specify it via command line in the following way:

			./run.py -w chosen_workflow -c config.json(.yaml) -other_arguments -q X 

In this "beta" phase, we have 2 standard options for "chosen_workflow" : Single (one sample) or Combined (multiple samples).  

To print the DAG of jobs, type:

			./run.py -d | dot -Tsvg > dag.svg

the expected output should be:
	
	-a "mapped_reads" directory with the bam file/s
	-a "genotyped" directory with the vcf file
	-a "sift_dbNSFP" directory with the final fully annotated file
	-a "log" directory wich contains all the error messages you might encounter during the workflow
	-a "stats" directory with some stats to be processed by scripts
	
future improvements:
	
	-update to DeepVariant/GlNexus
	-add new workflow options!! (cnv, somatic and mitochondrial variant calling and annotation)
	-add the SnpSift ExtractFields step based on the self made python script "extractfields.py"
	-add an html report 
