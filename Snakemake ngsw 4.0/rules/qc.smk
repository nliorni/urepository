#FASTQC
rule FastQC:
    input:
        reads=["samples/{sample}_1.fastq", "samples/{sample}_2.fastq"]   #edit to {sample}_{number}.fastq
    output:
        html="qc/fastqc/{sample}.html",
        zip="qc/fastqc/{sample}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    message:
        "Creating a fastqc report in the form of '{output.html}' and '{output.zip}' from the reads '{input.reads}' with FastQC"
    threads: 1
    wrapper:
        "0.67.0/bio/fastqc"



