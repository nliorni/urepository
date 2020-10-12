readme v1.0 Snakemake ngs pipeline

Greetings,

since here are used no wrappers, the version 1.0 works this way: 

	-create and name a conda environment from the provided environment.yaml file
	-activate such environment
	-put yourself in the folder in wich you found this readme file
	 (you can rename it as your current job)
	-put your genome.fa in "data" and your fastq samples in "data/samples"
	-rename your file respectevly "genome.fa" and "[A-Z].fastq"
	-adjust the placeholded config.yaml file as needed 
	-index your genome if not done before: bwa index genome.fa
	-cd to the main directory
	-download gatk4 in your conda environment using conda install gatk4 -y
	-perform a dry run to check if everything is all right: snakemake -n --cores 1
	-perform the run (glhf): snakemake --cores 1

this version supports the following functions:

	-maps your reads to the reference genome using bwa mem, creating a bam file in 
	 the "mapped" directory
	-sorts the alignment creating a sorted.bam file
	-indexes the alignment creating a .bai file
	-calls the variants in the sorted.bam file creating a all.vcf file
	-performs a minimal annotation with gatk
	-creates a plot of the calls qualities

future improvements:

	-use gatk HaplotypeCaller for the variant calling instead of bcftools
	-use VEP or SNPEFF for the annotation instead of gatk
	-use Picard for marking duplicates instead of samtools and make a metrics file
	-use gatk BaseRicalibrator and gatk ApplyBQSR (sites_of_variation.vcf file
	 needed)

	 
	
