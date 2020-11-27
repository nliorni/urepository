[![Snakemake](https://img.shields.io/badge/snakemake-≥5.6.0-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)


Greetings,

in this version wrappers will be used.

You are provided with a "data" directory, which contains a "samples" directory, and with a file called "environment.yaml".
Put your genome.fa (you have to rename it)(must be indexed with bwa index and samtools faidx) in "data", and your samples in "samples" (you can rename your sample (rename with "a_1(2).fastq")

You need to have Conda installed on your system.

Now, put yourself in the main directory and activate the minimal Snakemake environment from the provided file, typing:

	
	   		conda env create --name Snakemake-yourjob -f environment.yaml
	   
With wrappers, all the others packages needed will be provided and installed in your conda environment. 
Here is provided a run.py python script, which handles the run of the Snakefile for you. So, type:

	   		./run.py  -h  (for help)

You are able to limit the workflow to your necessities. This workflow supports "workflow.json" files, and you should have
several different files. Choose the workflow and specify it via command line in the following way:

			./run.py <workflowfile.json> -c 

Look at how the workflowfile.json is made for a better understanding. 
In this "alpha" phase, we have 4 standar options: Mapping, Calling, Annotating and Complete. 
You can edit one single workflow file, or you may have multiple files, one for each workflow. 
If this, for some reason, does not work, simply hard type:


			snakemake --use-conda target [type of workflow ...] --cores X

For the DAG of jobs, type:

			./run.py <workflowfile.json> -d | dot -Tsvg > dag.svg

the expected output should be:  (some of those are now temporary)
	
	-a "qc" directory with the .zip and .html report of fastqc
	-a "meqc" directory with the results of the MEQC.R script
	-an "icov" directory with the output of iCoverage
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
	-a "sift_annotated" directory with the output of snpsift annotate
	-a "sift_vartyped" directory with the output of snpsift vartype

UPDATES:
	
	-workflow is now dynamic
	-iCoverage works
	
future improvements:
	
	-add SNPSift dbNSFP wrapper rule
	-use costumized SNPeff databases
	-upgrade gatk call with deepvariant (the wrapper is currently not working)
	-add a BED file to limit the calling in gatk HaplotypeCaller

