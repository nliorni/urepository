#FASTQC
rule fastqc:
    input:
        reads=["data/samples/{sample}_1.fastq", "data/samples/{sample}_2.fastq"]
    output:
        html="qc/fastqc/{sample}.html",
        zip="qc/fastqc/{sample}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    log:
        "logs/fastqc/{sample}.log"
    message:
        "Creating a fastqc report in the form of '{output.html}' and '{output.zip}' from the reads '{input.reads}'"
    threads: 1
    wrapper:
        "0.67.0/bio/fastqc"
