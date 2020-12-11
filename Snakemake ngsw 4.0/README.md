[![Snakemake](https://img.shields.io/badge/snakemake-≥5.6.0-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)


Greetings,

in this version wrappers will be used.

It is assumed that you have your genome.fa, .fa.ann, .fa.amb, .fa.bwt, .fa.fai, .fa.pac, .fa.sai, .dict. You might need bwa Index and Samtools faidx. Specify the path in the config.json file. Also, you have to create a "data" directory inside the main directory in which you are working, create a "samples" directory inside, put your samples there and rename them like: a_R1.fastq/a_R2.fastq ... z_R1.fastq/z_R2.fastq. In the Snakefile (this is going to change very soon), search at the top for the list "SAMPLES", and add as many letters as many samples you have. 

If you have Conda installed on your system, you are pretty ready.

Now, create a working directory "snakemake-yourjob" and the minimal Snakemake environment from the provided file, typing:

		   	conda env create --name Snakemake-yourjob -f environment.yaml

Activate the environment by typing

			conda activate Snakemake-yourjob
	   
With wrappers, all the others packages needed will be provided and installed in your conda environment. 
Here is provided a run.py python script, which handles the run of the Snakefile for you. So, type:

	   		./run.py  -h  (for help)

You are able to limit the workflow to your necessities. Choose the workflow and specify it via command line in the following way:

			./run.py -w chosen_workflow -c config.json(.yaml) -other_arguments -q X 

In this "alpha" phase, we have 4 standar options for "chosen_workflow" : Mapping, Calling, Annotating and Complete.  

For the DAG of jobs, type:

			./run.py -d | dot -Tsvg > dag.svg

the expected output should be:
	
	-a "qc" directory with the .zip and .html report of fastqc
	-a "meqc" directory with the results of the MEQC.R script
	-a "mapped_reads" directory which contains the .bam file (fastq reference genome mapping)
	-a "reheaded_reads" with the reheaded .bam file (reheading)
	-a "dedup" directory which contains the .bam files without duplicates, and a metrics.txt file (remove (pcr) duplicates)
	-a "recal" folder, which contains the base recalibrated .bam file (recalibration)
	-a "sorted_reads" directory which contains the .sorted.bam and .sorted.bam.bai files (sorting)
	-a "calls" directory with .g.vcf file (variant calling)
	-an "annotated" directory which contains annotated .vcf, .html and .csv (annotation)
	-a "combined" directory which contains the combined gvcf from different samples gvcfs
	-a "genotyped" directory which contains the genotyped vcf from the combined gvcf
	-a "log" directory wich contains all the error messages you might encounter
	-a "stats" directory with some stats to be processed by scripts (work in progress)

UPDATES:
	
	-workflow is now dynamic without the need of a workflow.json file
	-iCoverage works
	
future improvements:
	
	-use costumized SNPeff databases
	-upgrade gatk call with deepvariant (the wrapper is currently not working)
	-add a BED file to limit the calling in gatk HaplotypeCaller
