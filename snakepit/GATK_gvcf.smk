# If snakemake_lsf is integrated, the creation of index and stats can become a local rule
# If snakemake_lsf is integrated, the log_folder and subfolders will no longer need manual creation

## setting variables
configfile: "config.yaml"
chroms = list(range(1,30))

#Reduced paths
vcf_genot_path = config["resources"]["vcf_genot"]
assembly = config["resources"]["assembly"]
vcf_beagle_path = config["fold_out"]["vcf_beagle"] + '/{family}/'

##Tools
GATK4 = config["tools"]["GATK4"]
LOAD_JAVA = config["tools"]["JAVA"]
BEAGLE = config["tools"]["BEAGLE"]

config['trios']='/cluster/work/pausch/alex/BSW_analysis/relationship.csv'
from pathlib import PurePath

def read_trios(ext='.vcf.gz'):
    if 'trios' not in config:
        return []
    import pandas as pd
    df = pd.read_csv(config['trios'])
    df.fillna('missing',inplace=True)

    targets = []
    for _, row in df.iterrows():
        for chr in range(1,30):
            Path(f'{config["fold_out"]["vcf_beagle"]}/{"_".join(row)}/{chr}').mkdir(parents=True,exist_ok=True)
        targets.append(f'{config["fold_out"]["vcf_beagle"]}/{"_".join(row)}/{"_".join(row)}.vcf.gz')
    return targets

rule all:
    input:
        read_trios()

config['gvcf']='/cluster/work/pausch/temp_scratch/BSW_remap/var_call/trios/haplotype_caller/gvcf'
config['GATK']='/cluster/work/pausch/audald/software/gatk-4.2.0.0/gatk'
config['ref']='/cluster/work/pausch/assemblies/ARS_UCD1.2/ARS_UCD.fasta'
config['out']='/cluster/work/pausch/temp_scratch/BSW_remap/var_call/trios/joint_genotyping'
config['trios']='/cluster/work/pausch/alex/BSW_analysis/relationship.csv'


rule combinegvcf:
    input:
        offspring = config['gvcf'] + '/{chr}/{offspring}_{chr}.g.vcf.gz',
        sire  = lambda wildcards: config['gvcf'] + '/{chr}/{sire}_{chr}.g.vcf.gz' if wildcards.sire != 'missing' else [],
        dam = lambda wildcards: config['gvcf'] + '/{chr}/{dam}_{chr}.g.vcf.gz' if wildcards.dam != 'missing' else []
    output:
        config['out'] + '/{offspring}_{sire}_{dam}/{chr}.gatk4.vcf.gz'
    params:
        dir_ = lambda wildcards, output: PurePath(output[0]).parent,
        gvcf = lambda wildcards,input: f'-V {input.offspring} ' + (f'-V {input.sire} ' if wildcards.sire!='missing' else '') + (f'-V {input.dam}' if wildcards.dam!='missing' else '')
    threads: 1
    resources:
        mem_mb = 6000,
        disk_scratch=30,
        walltime = '60'
    shell:
        '''
        mkdir -p {params.dir_}
        export TILEDB_DISABLE_FILE_LOCKING=1
        cd $TMPDIR
        cp {input} .
        {config[GATK]} \
        GenomicsDBImport \
        --genomicsdb-workspace-path db_{wildcards.chr} \
        -L {wildcards.chr} \
        --reader-threads {threads} \
        -R {config[ref]} \
        {params.gvcf}

        {config[GATK]} \
        GenotypeGVCFs  \
        -R {config[ref]}  \
        -L {wildcards.chr} \
        -O {output} \
        -V gendb://db_{wildcards.chr}
        '''

rule merge_trios:
    input:
        (vcf_beagle_path + f'{chr}/prune_var.vcf.gz' for chr in range(1,30))
    output:
        vcf_beagle_path + '{family}.vcf.gz'
    threads: 2
    resources:
        mem_mb = 2000,
        walltime = '20'
    shell:
        '''
        bcftools concat --threads {threads} -o {output} -Oz {input}
        tabix -p vcf {output}
        '''

rule select_SNP:
    input:
        vcf = vcf_genot_path + "{family}/{chromosome}.gatk4.vcf.gz",
        ref = assembly
    output:
        vcf_beagle_path + "{chromosome}/snp.vcf.gz"
    threads: 1
    resources:
        mem_mb = 3000,
        walltime = '30'
    shell:
        LOAD_JAVA +
        GATK4 +
        " SelectVariants  " +
        " -R {input.ref} " +
        " -V {input.vcf}  " +
        " --select-type-to-include SNP  " +
        " --output {output} "

rule filter_SNP:
    input:
        vcf = rules.select_SNP.output,
        ref = assembly
    output:
        vcf_beagle_path + "{chromosome}/snp_filtered.vcf.gz"
    threads: 1
    resources:
        mem_mb = 3000,
        walltime = '30'
    shell:
        LOAD_JAVA +
        GATK4 +
        "  VariantFiltration " +
        "  -R {input.ref} " +
        "  -V {input.vcf} " +
        " --filter \"QD < 2.0\" --filter-name \"snp_QD2\" "
        " --filter \"QUAL < 30.0\" --filter-name \"snp_QUAL30\" "
        " --filter \"SOR > 3.0\" --filter-name \"snp_SOR3\" "
        " --filter \"FS > 60.0\" --filter-name \"snp_FS60\" "
        " --filter \"MQ < 40.0\" --filter-name \"snp_MQ40\" "
        " --filter \"MQRankSum < -12.5\" --filter-name \"snp_MQRankSum-12.5\" "
        " --filter \"ReadPosRankSum < -8.0\" --filter-name \"snp_ReadPosRankSum-8\" "
        " --missing-values-evaluate-as-failing true "
        " --output {output} "

rule select_indel:
    input:
        vcf = vcf_genot_path + "{family}/{chromosome}.gatk4.vcf.gz",
        ref = assembly
    output:
        vcf_beagle_path + "{chromosome}/indel.vcf.gz"
    threads: 1
    resources:
        mem_mb = 3000,
        walltime = '30'    
    shell:
        LOAD_JAVA +
        GATK4 +
        " SelectVariants  " +
        " -R {input.ref} " +
        " -V {input.vcf}  " +
        " --select-type-to-include INDEL " +
        " --output {output} "


rule filter_indel:
    input:
        vcf = rules.select_indel.output,
        ref = assembly
    output:
        vcf_beagle_path + "{chromosome}/indel_filtered.vcf.gz"
    threads: 1
    resources:
        mem_mb = 3000,
        walltime = '30'
    shell:
        LOAD_JAVA +
        GATK4 +
        "  VariantFiltration " +
        "  -R {input.ref} " +
        "  -V {input.vcf} " +
        " --filter \"QD < 2.0\" --filter-name \"indel_QD2\" "
        " --filter \"QUAL < 30.0\" --filter-name \"indel_QUAL30\" "
        " --filter \"SOR > 10.0\" --filter-name \"indel_SOR10\" "
        " --filter \"FS > 200.0\" --filter-name \"indel_FS200\" "
        " --filter \"ReadPosRankSum < -20.0\" --filter-name \"indel_ReadPosRankSum-20\" "
        " --missing-values-evaluate-as-failing true "
        " --output {output} "

rule merge_variants:
    input:
        sic=rules.filter_SNP.output,
        idc=rules.filter_indel.output
    output:
        vcf_beagle_path + "{chromosome}/merge_var.vcf.gz"
    threads: 1
    resources:
        mem_mb = 4000,
        walltime = '30'
    shell:
        LOAD_JAVA +
        GATK4 +
        " MergeVcfs " +
        " --INPUT {input.sic} " +
        " --INPUT {input.idc} " +
        " --OUTPUT {output} "

rule remove_filtered:
    input:
        var = rules.merge_variants.output
    output:
        vcf_beagle_path + "{chromosome}/prune_var.vcf.gz"
    threads: 1
    resources:
        mem_mb = 3000,
        walltime = '30'
    shell:
        LOAD_JAVA +
        GATK4 +
        " SelectVariants "
        " -V {input.var} "
        " --exclude-filtered "
        "--output {output}"

rule beagle_imputation:
    input:
        rules.remove_filtered.output
    output:
        vcf_beagle_path + "{chromosome}/beagle.vcf.gz"
    params:
        lambda wildcards, output: PurePath(output[0]).with_suffix('').with_suffix('')
    shell:
        LOAD_JAVA +
        " java -Xss25m -Xmx40G -jar " + BEAGLE +
        " gl={input} " +
        " nthreads=18 " +
        " {params}"

rule index_creation:
    input:
        select_SNP = rules.select_SNP.output,
        filter_SNP = rules.filter_SNP.output,
        select_indel = rules.select_indel.output,
        filter_indel = rules.filter_indel.output,
        merge_variants = rules.merge_variants.output,
        remove_filtered = rules.remove_filtered.output,
        beagle_imputation = rules.beagle_imputation.output
    output:
        select_SNP = vcf_beagle_path + "{chromosome}/snp.vcf.gz.tbi",
        filter_SNP = vcf_beagle_path + "{chromosome}/snp_filtered.vcf.gz.tbi",
        select_indel = vcf_beagle_path + "{chromosome}/indel.vcf.gz.tbi",
        filter_indel = vcf_beagle_path + "{chromosome}/indel_filtered.vcf.gz.tbi",
        merge_variants = vcf_beagle_path + "{chromosome}/merge_var.vcf.gz.tbi",
        remove_filtered = vcf_beagle_path + "{chromosome}/prune_var.vcf.gz.tbi",
        beagle_imputation = vcf_beagle_path + "{chromosome}/beagle.vcf.gz.tbi"
    threads: 1
    resources:
        mem_mb = 3000,
        walltime = '30'
    shell:
        "tabix -fp vcf {input.select_SNP} \n" +
        "tabix -fp vcf {input.filter_SNP} \n" +
        "tabix -fp vcf {input.select_indel} \n" +
        "tabix -fp vcf {input.filter_indel} \n" +
        "tabix -fp vcf {input.merge_variants} \n" +
        "tabix -fp vcf {input.remove_filtered} \n" +
        "tabix -fp vcf {input.beagle_imputation} \n"

rule stats_creation:
    input:
        select_SNP = rules.select_SNP.output,
        filter_SNP = rules.filter_SNP.output,
        select_indel = rules.select_indel.output,
        filter_indel = rules.filter_indel.output,
        merge_variants = rules.merge_variants.output,
        remove_filtered = rules.remove_filtered.output,
        beagle_imputation = rules.beagle_imputation.output,
        beagle_index = rules.index_creation.output.beagle_imputation
    output:
        select_SNP = vcf_beagle_path + "{chromosome}/snp.stats",
        filter_SNP = vcf_beagle_path + "{chromosome}/snp_filtered.stats",
        select_indel = vcf_beagle_path + "{chromosome}/indel.stats",
        filter_indel = vcf_beagle_path + "{chromosome}/indel_filtered.stats",
        merge_variants = vcf_beagle_path + "{chromosome}/merge_var.stats",
        remove_filtered = vcf_beagle_path + "{chromosome}/prune_var.stats",
        beagle_imputation = vcf_beagle_path + "{chromosome}/beagle.stats"
    threads: 1
    resources:
        mem_mb = 3000,
        walltime = '30'
    shell:
        "bcftools stats {input.select_SNP} > {output.select_SNP} \n" +
        "bcftools stats {input.filter_SNP} > {output.filter_SNP} \n" +
        "bcftools stats {input.select_indel} > {output.select_indel} \n" +
        "bcftools stats {input.filter_indel} > {output.filter_indel} \n" +
        "bcftools stats {input.merge_variants} > {output.merge_variants} \n" +
        "bcftools stats {input.remove_filtered} > {output.remove_filtered} \n" +
        "bcftools stats {input.beagle_imputation} > {output.beagle_imputation} \n"
