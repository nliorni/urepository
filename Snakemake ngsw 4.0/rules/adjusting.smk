#REHEADING WITH PICARD
rule picard_reheading:
    input:
        "mapped_reads/{sample}.bam"
    output:
        temp("reheaded_reads/{sample}.bam")
    log:
        "logs/picard/replace_rg/{sample}.log"
    message:
        "Reheading '{input}' into '{output}' with Picard AddOrReplaceGroups"
    params:
        "RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=20"
    wrapper:
        "0.67.0/bio/picard/addorreplacereadgroups"

#MARK AND REMOVE DUPLICATES WITH PICARD
rule picard_mark_duplicates:
    input:
        "reheaded_reads/{sample}.bam"
    output:
        bam=temp("dedup/{sample}.bam"),
        metrics="dedup_metrics/{sample}.metrics.txt"
    log:
        "logs/picard/dedup/{sample}_dedup.log"
    message:
        "Marking and removing duplicates from '{input}' into '{output.bam}' and writing metrics into '{output.metrics}' with Picard MarkDuplicates"
    params:
        "REMOVE_DUPLICATES=true"
    wrapper:
        "0.66.0/bio/picard/markduplicates"

#fare in modo che questa regola venga eseguita solo SE MANCA Il file genome.dict... qualcosa come if os.path.isfile("data/genome.dict")==False:
#CREA FILE GENOME.DICT
#rule create_dict:
#    input:
#        genome=config["reference"]
#    output:
#        "data/genome.dict"
#    log:
#        "logs/picard/create_dict.log"
#    message:
#        "Creating reference genome dictionary with Picard CreateSequenceDictionary, because no genome.dict was found"
#    params:
#        extra=""  # optional: extra arguments for picard.
#    wrapper:
#        "0.66.0/bio/picard/createsequencedictionary"

#RECALIBRATION TABLE WITH GATK
rule gatk_baserecalibrator:
    input:
        bam="dedup/{sample}.bam",
        ref=config["reference"],
        dic=config["dict"],
        known=config["snpdb"]  # optional known sites
    output:
        recal_table=temp("recal/{sample}.grp")
    log:
        "logs/gatk/baserecalibrator/{sample}.log"
    message:
        "Recalibration of  {input.bam} based on {input.known} creating {output} with Gatk BaseRecalibrator "
    params:
        extra="",  # optional
        java_opts="", # optional
    wrapper:
        "0.66.0/bio/gatk/baserecalibrator"

#APPLY THE CREATED TABLE FOR BASE RECALIBRATION
rule gatk_applybqsr:
    input:
        bam="dedup/{sample}.bam",
        ref=config["reference"],
        dictio=config["dict"],
        recal_table="recal/{sample}.grp"
    output:
        bam=temp("recal/{sample}.bam")
    log:
        "logs/gatk/gatk_applybqsr/{sample}.log"
    message:   
        "Applying the '{input.recal_table}' to '{input.bam}' outputing '{output}' with Gatk ApplyBQSR"
    params:
        extra="",  # optional
        java_opts="", # optional
    wrapper:
        "0.66.0/bio/gatk/applybqsr"

##local realignment around indels???
