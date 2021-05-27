if config["fastqc"]==True:

## FastQC ##
    rule FastQC:
        input:
            reads=lambda wildcards: expand(f"{config['samples'][wildcards.sample]}_{{unit}}.fastq.gz", unit=["L001_R1_001","L001_R2_001"])
        output:
            html="results/qc/fastqc/{sample}_{unit}.html",
            zip="results/qc/fastqc/{sample}_{unit}.zip" 
        params: config["FastQC"]["params"]
        message:
            "Running FastQC. Creating a report in the form of '{output.html}' and '{output.zip}' from the reads '{input.reads}'."
        wrapper:
            "0.67.0/bio/fastqc"
        



