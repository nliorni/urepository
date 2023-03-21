samples = config["cohort"]
## Gatk CombineGvcfs ##
rule gatk_CombineGvcfs:
    input:
        gvcfs = [v for k,v in samples.items()],
        #gvcfs=expand("{sample}.g.vcf", sample=config["cohort"]),
        ref=config["reference"]
    output:
        outputDir+"results/cohort/cohort.g.vcf"
    params:
        extra="",
        java_opts="",
        variant=lambda wildcards, input: ' -V '.join(input.gvcfs)
    message:
        "running gatk combinegvcfs (cohort) for samples {input.gvcfs}"
    shell:
        "gatk --java-options '{params.java_opts}' CombineGVCFs {params.extra} -V {params.variant} -R {input.ref} -O {output}"

## Gatk GenotypeGvcfs ##
rule gatk_GenotypeGvcfs_cohort:
    input:
        gvcf=outputDir+"results/cohort/cohort.g.vcf",  # combined gvcf over multiple samples calls/all.g.vcf
        ref=config["reference"]
    output:
        vcf=temp(outputDir+"results/cohort/cohort.vcf"),
        idx=temp(outputDir+"results/cohort/cohort.vcf.idx")
    message:
        "running gatk genotypegvcfs (cohort)"
    params:
        extra=config["gatk_GenotypeGvcfs"]["params"]["extra"],  
        java_opts=config["gatk_GenotypeGvcfs"]["params"]["java_opts"], 
    shell:
        "gatk --java-options {params.java_opts} GenotypeGVCFs --variant {input.gvcf} --reference {input.ref} --output {output.vcf}"

## Gatk SelectVariants (Snps) ##
rule gatk_snps_SelectVariants_cohort:
    input:
        vcf=outputDir+"results/cohort/cohort.vcf",
        ref=config["reference"],
    output:
        vcf=temp(outputDir+"results/cohort/cohort.snps.vcf"),
        idx=temp(outputDir+"results/cohort/cohort.snps.vcf.idx")
    message:
        "running gatk selectvariants (snps) (cohort)"
    params:
        extra=config["gatk_SelectVariants"]["snps"]["extra"],  # optional filter arguments, see GATK docs
        java_opts=config["gatk_SelectVariants"]["snps"]["java_opts"]
    shell:
        "gatk SelectVariants --variant {input.vcf} --reference {input.ref} {params.extra} --output {output.vcf}"

## Gatk SelectVariants (Indels) ##
rule gatk_indel_SelectVariants_cohort:
    input:
        vcf=outputDir+"results/cohort/cohort.vcf",
        ref=config["reference"],
    output:
        vcf=temp(outputDir+"results/cohort/cohort.indel.vcf"),
        idx=temp(outputDir+"results/cohort/cohort.indel.vcf.idx")
    message:
        "running gatk selectvariants (indels) (cohort)"
    params:
        extra=config["gatk_SelectVariants"]["indel"]["extra"],  # optional filter arguments, see GATK docs
        java_opts=config["gatk_SelectVariants"]["indel"]["java_opts"]
    shell:
        "gatk SelectVariants --variant {input.vcf} --reference {input.ref} {params.extra} --output {output.vcf}"

## Gatk VariantFiltration (Snps) ##   
rule gatk_snps_VariantFiltration_cohort:
    input:
        vcf=outputDir+"results/cohort/cohort.snps.vcf",
        ref=config["reference"],
        vcf_idx = outputDir+"results/cohort/cohort.snps.vcf.idx"
    output:
        vcf=temp(outputDir+"results/cohort/cohort.snps.filtered.vcf"),
        idx=temp(outputDir+"results/cohort/cohort.snps.filtered.vcf.idx")
    message:
        "running gatk variantfiltration (snps) (cohort)"
    params:
        filters=config["gatk_FilterVariants"]["snps"]["filters"],
        extra=config["gatk_FilterVariants"]["snps"]["extra"],  
        java_opts=config["gatk_FilterVariants"]["snps"]["java_opts"]
    shell:
        "gatk VariantFiltration --variant {input.vcf} --reference {input.ref} --cluster-size 3 --cluster-window-size 10 --filter-name '100QUAL' --filter-expression 'QUAL < 100' --filter-name 'SnpWarning' --filter-expression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' --output {output.vcf}"

## Gatk VariantFiltration (Indels) ##
rule gatk_indel_VariantFiltration_cohort:
    input:
        vcf=outputDir+"results/cohort/cohort.indel.vcf",
        ref=config["reference"],
        vcf_idx = outputDir+"results/cohort/cohort.indel.vcf.idx"
    output:
        vcf=temp(outputDir+"results/cohort/cohort.indel.filtered.vcf"),
        idx=temp(outputDir+"results/cohort/cohort.indel.filtered.vcf.idx")
    message:
        "running gatk variantfiltration (indels) (cohort)"
    params:
        filters=config["gatk_FilterVariants"]["indels"]["filters"],
        extra=config["gatk_FilterVariants"]["indels"]["extra"],  
        java_opts=config["gatk_FilterVariants"]["indels"]["java_opts"]
    shell:
        "gatk VariantFiltration --variant {input.vcf} --reference {input.ref} --cluster-size 3 --cluster-window-size 10 --filter-name '100QUAL' --filter-expression 'QUAL < 100' --filter-name 'IndelWarning' --filter-expression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0' --output {output.vcf}"


## Picard MergeVcfs ##
rule picard_MergeVcfs_cohort:
    input:
        snps=outputDir+"results/cohort/cohort.snps.filtered.vcf", 
        indel=outputDir+"results/cohort/cohort.indel.filtered.vcf"
    output:
        vcf=outputDir+"results/cohort/cohort.filter.vcf",
        idx=temp(outputDir+"results/cohort/cohort.filter.vcf.idx")
    message:
        "running picard mergevcfs (cohort)"
    params:
        extra=config["picard_MergeVcfs"]["extra"]
    shell:
        "picard MergeVcfs {params.extra} --INPUT {input.snps} --INPUT {input.indel} --OUTPUT {output.vcf}"

## VEP Annotate ##
rule VEP_annotate_cohort:
    input:
        calls=outputDir+"results/cohort/cohort.filter.vcf",  # .vcf, .vcf.gz or .bcf
        cache=config["vep_annotate"]["cache"], # can be omitted if fasta and gff are specified
        plugins="/software/dataseq/databases/vep/plugins",
    output:
        calls=temp(outputDir+"results/cohort/cohort.vep.vcf"),  # .vcf, .vcf.gz or .bcf
        stats=outputDir+"results/cohort/vep_stats/cohort.html"
    params:
        plugins=config["vep_annotate"]["plugins"],
        extra=config["vep_annotate"]["extra"] # optional: extra arguments
    log:
        outputDir+"logs/vep/cohort.vepanno.log"
    message: "annotating with VEP (cohort)"
    threads: 20
    wrapper:
        "v1.15.0/bio/vep/annotate"

## SnpSift VariantType ##
rule snpsift_vep_VarType_cohort:
    input:
        vcf=outputDir+"results/cohort/cohort.vep.vcf"
    output:
        vcf=temp(outputDir+"results/cohort/cohort.vep.vartype.vcf")
    message:
        "running snpsift vartype (cohort)"
    params:
        java_opts=config["snpsift_VarType"]["params"]["java_opts"],
        extra=config["snpsift_VarType"]["params"]["extra"],
        snpsift_path=config["snpsift_path"]
    log:
        outputDir+"logs/snpsift/vartype/cohort.log"
    shell:
        "java -jar {params.snpsift_path} VarType {input.vcf} > {output.vcf}"


## SnpSift ExtractFields - 1st parsing step ##
rule snpsift_vep_ExtractFields_cohort:
    input:
        outputDir+"results/cohort/cohort.vep.vartype.vcf"
    output:
        temp(outputDir+"results/cohort/cohort.vep.ssef.tsv")
    log:
        outputDir+"logs/snpsift/ExtractFields/cohort.log"
    params:
        fieldsfile=config["vcf_fieldsfile"]
    message:
        "running snpsift extractfields (parsing vcf fields) (cohort)"
    shell:
        "python3 workflow/scripts/extractfields.py --input {input} --fieldsfile {params.fieldsfile} > {output}"    

## VEP ExtractFields - 2nd parsing step ##
rule vep_extractFields_cohort:
    input:
        tsv=outputDir+"results/cohort/cohort.vep.ssef.tsv",
        vcf=outputDir+"results/cohort/cohort.vep.vcf"
    output:
        temp(outputDir+"results/cohort/cohort.vep.tsv")
    params:
        fieldsfile=config["vep_fieldsfile"]
    message:
        "running VEP extractfields (parsing CSQ fields) (cohort)"
    shell:
        "python3 workflow/scripts/vep_extractfields.py --input {input.vcf} --fieldsfile {params.fieldsfile} --tsv {input.tsv}"

## final parsing ##
rule final_parsing_genelist_cohort:
    input:
        outputDir+"results/cohort/cohort.vep.tsv"
    output:
        temp(outputDir+"results/cohort/cohort_parsed.xlsx")
    params:
        genelists = config["genelists"] if "genelists" in config and config["genelists"] else []
    message:
        "final operations for the cohort: formatting columns, checking predictors and adding gene lists"
    shell:
        "python3 workflow/scripts/final_parsing_cohort.py -i {input} -g {params.genelists} -o {output}"

## Add genotype columns ##
rule add_gt_columns:
    input:
        excel = outputDir+"results/cohort/cohort_parsed.xlsx",
        vcf = outputDir+"results/cohort/cohort.vep.vcf"
    output:
        outputDir+"results/cohort/cohort.xlsx"
    message:
        "adding genotype columns"
    shell:
        "python3 workflow/scripts/cohort_genotype.py --vcf {input.vcf} --excel {input.excel} --output {output}"

