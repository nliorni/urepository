if config["fastqc"]==True:

## FastQC ##
    rule FastQC:
        input:
            reads=lambda wildcards: expand(f"{config['samples'][wildcards.sample]}_{{num}}.fastq", num=[1,2])
        output:
            html="results/qc/fastqc/{sample}_{num}.html",
            zip="results/qc/fastqc/{sample}_{num}.zip" 
        params: config["FastQC"]["params"]
        message:
            "Running FastQC. Creating a report in the form of '{output.html}' and '{output.zip}' from the reads '{input.reads}'."
        wrapper:
            "0.67.0/bio/fastqc"
        



