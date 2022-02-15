## Run as follow:
#snakemake --jobs 100 -nrp --cluster-config cluster.json --cluster "bsub -J {cluster.jobname} -n {cluster.ncore} -W {cluster.jobtime} -o {cluster.logi} -R \"rusage[mem={cluster.memo}]\""

## parsing the required files

# Deal with the lsf logfile

from pathlib import Path, PurePath

config['reference'] ='/cluster/work/pausch/inputs/ref/SSC/11.1/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa'
config['bam_path'] = '/nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/SSC/bam'

BAMDIR = config['bam_path']


## Parsing output directory

##Parsing wildcards
chr_list = list(range(1,19))
int_list = " ".join(["-L {}".format(i) for i in chr_list])

config['gvcf']='/nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/alex/DV_GATK_project/GATK_piggs/gvcf'
config['out']='/nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/alex/DV_GATK_project/GATK_piggs/joint_geno'
config['vcf_beagle_path'] = '/nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/alex/DV_GATK_project/GATK_piggs/filtered'
vcf_beagle_path = config['vcf_beagle_path']
OUTDIR = config['out']
config['db'] = '/cluster/work/pausch/adela/cache/pig_db_corr.vcf.gz'
config['chip'] = '/cluster/work/pausch/adela/cache/pig_geno_corr.vcf.gz'

rule all:
    input:
        expand( vcf_beagle_path + "{chromosome}/prune_var.vcf.gz", chromosome = chr_list)

### Recalibration

rule recalibrator_creator:
    input:
        BAMDIR + "/{sample}.bam"
    output:
        OUTDIR + "/bam_recal/{sample}_recalibrator.table"
    threads: 1
    resources:
        mem_mb = 5000
    params:
        chr = int_list
    shell:
        '''
        gatk BaseRecalibrator -I {input} {params.chr} -R {config[reference]} --known-sites {config[db]} --known-sites {config[chip]} --lenient true -O {output}
        '''

## TO DO: print reads each chromosome separately and joined 
## AS SUGGESTED by BULL GENOME PROJECTS
rule base_recalibrator:
    input:
        recal = rules.recalibrator_creator.output,
        sample = BAMDIR + "/{sample}.bam"
    output:
        temp(OUTDIR + "/bam_recal/{sample}_recalibrated.bam")
    threads: 1
    resources:
        mem_mb = 5000
    params:
        chr = int_list
    shell:
        '''
        gatk ApplyBQSR -R {config[reference]} {params.chr} -I {input.sample} --bqsr-recal-file {input.recal} -O {output}
        '''

# Scatter by sample and chromosome

rule Haplotype_caller:
    input:
        rules.base_recalibrator.output
    output:
        config['gvcf'] + "/{chromosome}/{sample}.g.vcf.gz"
    threads: 1
    resources:
        mem_mb = 10000,
        walltime = '24:00'
    shell:
        '''
        gatk --java-options "-Xmx10g -XX:ParallelGCThreads=2" HaplotypeCaller -I {input} -R {config[reference]} -L {wildcards.chromosome} -O {output} --ERC GVCF
        '''

rule combinegvcf:
    input:
        offspring = (config['gvcf'] + f'/{{chr}}/{sample}.g.vcf.gz' for sample in config['animals']),
    output:
        config['out'] + '/{chr}.gatk4.vcf.gz'
    params:
        dir_ = lambda wildcards, output: PurePath(output[0]).parent,
        gvcf = (f'-V {x}' for x in config['animals'])
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
        gatk \
        GenomicsDBImport \
        --genomicsdb-workspace-path db_{wildcards.chr} \
        -L {wildcards.chr} \
        --reader-threads {threads} \
        -R {config[reference]} \
        {params.gvcf}

        gatk \
        GenotypeGVCFs  \
        -R {config[reference]}  \
        -L {wildcards.chr} \
        -O {output} \
        -V gendb://db_{wildcards.chr}
        '''


rule filter_type:
    input:
        vcf = vcf_beagle_path + "{chromosome}/{type}.selected.vcf.gz"
    output:
        vcf_beagle_path + "{chromosome}/{type}_filtered.vcf.gz"
    params:
        filters = lambda wildcards: '--filter "QD < 2.0" --filter-name "snp_QD2" --filter "QUAL < 30.0" --filter-name "snp_QUAL30" --filter "SOR > 3.0" --filter-name "snp_SOR3" --filter "FS > 60.0" --filter-name "snp_FS60" --filter "MQ < 40.0" --filter-name "snp_MQ40" --filter "MQRankSum < -12.5" --filter-name "snp_MQRankSum-12.5" --filter "ReadPosRankSum < -8.0" --filter-name "snp_ReadPosRankSum-8"' if wildcards.type == 'snp' else '--filter "QD < 2.0" --filter-name "indel_QD2" --filter "QUAL < 30.0" --filter-name "indel_QUAL30" --filter "SOR > 10.0" --filter-name "indel_SOR10" --filter "FS > 200.0" --filter-name "indel_FS200" --filter "ReadPosRankSum < -20.0" --filter-name "indel_ReadPosRankSum-20"'
    threads: 1
    resources:
        mem_mb = 3000,
        walltime = '30'
    shell:
        '''
        gatk VariantFiltration -R {config[reference]} -V {input.vcf} \
        {params.filters} --output {output}
        '''

rule select_variants:
    input:
        vcf = config['out'] + "/{chromosome}.gatk4.vcf.gz",
    output:
        vcf_beagle_path + "{chromosome}/{type}.selected.vcf.gz"
    params:
        type = lambda wildcards: 'INDEL' if wildcards.type == 'indel' else 'SNP'
    threads: 1
    resources:
        mem_mb = 3000,
        walltime = '30'
    shell:
        '''
        gatk SelectVariants -R {config[reference]} -V {input.vcf} --select-type-to-include {params.type} --output {output}
        '''

rule remove_filtered:
    input:
        snps = vcf_beagle_path + "{chromosome}/snp_filtered.vcf.gz",
        indels = vcf_beagle_path + "{chromosome}/indel_filtered.vcf.gz"
    output:
        merged = temp(vcf_beagle_path + "{chromosome}/merge_var.vcf.gz"),
        pruned = vcf_beagle_path + "{chromosome}/prune_var.vcf.gz"
    threads: 1
    resources:
        mem_mb = 3000,
        walltime = '30'
    shell:
        '''
        gatk MergeVcfs --INPUT {input.snps} --INPUT {input.indels} --OUTPUT {output.merged}
        gatk SelectVariants -V {output.merged} --exclude-filtered --output {output.pruned}
        '''
