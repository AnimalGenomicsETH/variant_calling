from pathlib import PurePath

rule all:
    input:
        'merfin_filtering.csv',
        expand('merfin/{sample}.{vcf}.merfin.{fitted}.filter.vcf',sample=config['samples'],vcf=config['vcfs'],fitted=config['modes'])

def get_fastq(sample):
    return (str(PurePath(config['fastq']) / f'{sample}_R{N}.fastq.gz') for N in (1,2))

rule meryl_count_reads:
    input:
        lambda wildcards: get_fastq(wildcards.sample)
    output:
        temp(directory('readmers/{sample}.readmers.SR.meryl'))
    threads: 24
    resources:
        mem_mb = 5000
    params:
        mem = lambda wildcards,resources,threads: resources['mem_mb']*threads/1100
    shell:
        '/cluster/work/pausch/alex/software/MERFIN_NEW/build/bin/meryl count k=21 memory={params.mem} threads={threads} output {output} {input}'

rule meryl_count_asm:
    output:
        directory('readmers/ARS.seqmers.meryl')
    threads: 6
    resources:
        mem_mb = 3500
    params:
        mem = lambda wildcards,resources,threads: resources['mem_mb']*threads/1100
    shell:
        '/cluster/work/pausch/alex/software/MERFIN_NEW/build/bin/meryl count k=21 memory={params.mem} threads={threads} output {output} {config[reference]}'

rule meryl_print_histogram:
    input:
        'readmers/{sample}.readmers.SR.meryl'
    output:
        'readmers/{sample}.hist'
    shell:
        '/cluster/work/pausch/alex/software/MERFIN_NEW/build/bin/meryl histogram {input} > {output}'

rule merfin_lookup:
    input:
        'readmers/{sample}.hist'
    output:
        'readmers/{sample}/lookup_table.txt',
        'readmers/{sample}/model.txt'
    params:
        out = lambda wildcards, output: PurePath(output[0]).parent,
        ploidy = 2#lambda wildcards: 2 if wildcards.haplotype == 'asm' else 2
    envmodules:
        'gcc/8.2.0',
        'r/4.0.2'
    shell:
        'Rscript /cluster/work/pausch/alex/software/genomescope2.0/genomescope.R -i {input} -k 21 -o {params.out} --fitted_hist -p {params.ploidy}'

#rule bcftools_extract:
#    input:
#        vcf = lambda wildcards: config['vcfs'][wildcards.vcf]
#    output:
#        temp('

rule merfin_filter:
    input:
        seqmers = 'readmers/ARS.seqmers.meryl',
        readmers = 'readmers/{sample}.readmers.SR.meryl',
        vcf = lambda wildcards: config['vcfs'][wildcards.vcf],
        lookup = 'readmers/{sample}/lookup_table.txt',
        model = 'readmers/{sample}/model.txt'
    output:
        'merfin/{sample}.{vcf}.merfin.{fitted}.filter.vcf'
    params:
        out = lambda wildcards, output: PurePath(output[0]).with_suffix('').with_suffix(''),
        fitted = lambda wildcards, input: f'-prob {input.lookup} -peak {float([line for line in open(input.model)][8].split()[1])}' if wildcards.fitted == 'fitted' else ''
    threads: 12
    resources:
        mem_mb = 7000,
        walltime = '24:00'
    shell:
        '''
        #bcftools view -s {wildcards.sample} | bcftools view 'GT[*]="alt"' 
        /cluster/work/pausch/alex/software/MERFIN_NEW/build/bin/merfin -filter -debug -sequence {config[reference]} -seqmers {input.seqmers} {params.fitted} -readmers {input.readmers} -threads {threads} -vcf <(bcftools view -s {wildcards.sample} {input.vcf} | bcftools view -i 'GT[*]="alt"') -output {params.out}
        '''

rule tabulate_results:
    input:
        expand('merfin/{sample}.{vcf}.merfin.{fitted}.filter.vcf',sample=config['samples'],vcf=config['vcfs'],fitted=config['modes'])
    output:
        'merfin_filtering.csv'
    shell:
        '''
        echo sample,DV,DV_merfin,GATK,GATK_merfin > {output}
        for sample in {params.samples};
        do
          echo $sample,$(grep "records:" $(ls -rt logs/merfin_filter/sample-$sample.vcf-DV.fitted-raw*err | head -n 1) | awk '{{print $2}}'),$(wc -l merfin/$sample.DV.merfin.raw.filter.vcf),$(grep "records:" $(ls -rt logs/merfin_filter/sample-$sample.vcf-GATK.fitted-raw*err | head -n 1) | awk '{{print $2}}'),$(wc -l merfin/$sample.GATK.merfin.raw.filter.vcf) >> {output}

        '''
