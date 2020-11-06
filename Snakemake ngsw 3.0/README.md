Greetings,

in this version wrappers will be used.

You are provided with a "data" directory, which contains a "samples" directory, and with a file called "environment.yaml".
Put your genome.fa (you have to rename it)(must be indexed with bwa index and samtools faidx) in "data", and your samples in "samples" (you can rename your sample (rename with "a_1(2).fastq")

You need to have Conda installed on your system.

Now, put yourself in the main directory and activate the minimal Snakemake environment from the provided file, typing:

	
	   		conda env create --name Snakemake-yourjob -f environment.yaml
	   

With the new implementation, it should be necessary to download bedtools, because Snakemake does not support a wrapper for it. Therefore:

			conda install bedtools

(Please, wake up  and insert bedtools in the environment.yaml file!)
With wrappers, all the others packages needed will be provided and installed in your conda environment. 
Here is provided a run.py python script, which handles the run of the Snakefile for you. So, type:

	   		./run.py  -h  (for help)
			./run.py   

If this, for some reason, does not work, simply hard type:

			snakemake --use-conda --cores X
			
For the DAG of jobs, type:

			snakemake --dag | dot -Tsvg > dag.svg

the expected output should be:
	
	-a "qc" directory with the .zip and .html report of fastqc
	-a "mapped_reads" directory which contains the .bam file (fastq reference genome mapping)
	-a "reheaded_reads" with the reheaded .bam file (reheading)
	-a "dedup" directory which contains the .bam files without duplicates, and a metrics.txt file (remove (pcr) duplicates)
	-a "recal" folder, which contains the base recalibrated .bam file (recalibration)
	-a "sorted_reads" directory which contains the .sorted.bam and .sorted.bam.bai files (sorting)
	-a "calls" directory with .g.vcf file (variant calling)
	-an "snpeff" directory which contains annotated .vcf, .html and .csv (annotation)
	-a "log" directory wich contains all the error messages you might encounter
	-a "stats" directory with some stats to be processed by scripts (work in progress)

UPDATES:
	
	-Added a python executable for the workflow
	-config.yaml file now is necessary, working and editable!
	-Workflow is now modularized
	-Added Picard AddOrReplaceReadGroups for reheading the bam file. (needed for BaseRecalibrator)
	-Added gatk BaseRecalibrator and ApplyBQSR.
	-Annotation with SNPEFF is implemented (I still had to manually download the snpeff database, for some problems with the wrapper, I will soon implement the rule itself).
	-Fastqc working, please specify the output in the snakemake call
	
future improvements:
	
	-make costumized scripts work
	-upgrade gatk call with deepvariant (the wrapper is currently not working)
	-generalize the workflow: config.yaml -->> UPGRADE HAS STARTED :)
	-add a BED file to limit the calling in gatk HaplotypeCaller
	-autodownload genome, get files from user interactively
	-implement a params.json file, able to swith path inside the workflow (which rule to execute, and which not), and pass it	
	  to the run.py script
