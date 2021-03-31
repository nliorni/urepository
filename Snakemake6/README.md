Greetings,

It is assumed that you have your genome.fa, .fa.ann, .fa.amb, .fa.bwt, .fa.fai, .fa.pac, .fa.sai, .dict. You might need bwa Index and Samtools faidx. Specify the paths in the config file. Also, you have to create a "samples" directory inside the main directory in which you are working, put there your samples. After that, open the configuration file and write there the path of your samples.  

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

			./run.py {previous arguments defining workflow} -d | dot -Tsvg > dag.svg

the expected output will be included in a results folder.
	
future improvements:

	-add new workflow options!! (cnv, somatic and mitochondrial variant calling and annotation)
	-add an html report 
	-migrate workflow options to config file and wrap the Snakefile
