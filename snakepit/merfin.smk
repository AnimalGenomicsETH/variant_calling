from pathlib import PurePath
from itertools import product

localrules: tabulate_isec

rule all:
    input:
        'merfin_isec.upset',
        'merfin_filtering.csv',
        expand('merfin/{sample}_merfin.intersections',sample=config['samples']),
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

rule meryl_filter:
    input:
        'readmers/{sample}.readmers.SR.meryl'
    output:
        temp(directory('readmers/{sample}.readmers.gt1.SR.meryl'))
    resources:
        mem_mb = 5000
    shell:
        '''
        /cluster/work/pausch/alex/software/MERFIN_NEW/build/bin/meryl greater-than 1 {input} output {output}
        '''

rule meryl_count_asm:
    output:
        temp(directory('readmers/ARS.seqmers.meryl'))
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
        temp('readmers/{sample}.hist')
    shell:
        '/cluster/work/pausch/alex/software/MERFIN_NEW/build/bin/meryl histogram {input} > {output}'

rule merfin_lookup:
    input:
        'readmers/{sample}.hist'
    output:
        temp('readmers/{sample}/lookup_table.txt'),
        temp('readmers/{sample}/model.txt')
    params:
        out = lambda wildcards, output: PurePath(output[0]).parent
    envmodules:
        'gcc/8.2.0',
        'r/4.1.3'
    shell:
        'Rscript /cluster/work/pausch/alex/software/genomescope2.0/genomescope.R -i {input} -k 21 -o {params.out} --fitted_hist -p 2'

rule bcftools_extract:
    input:
        vcf = lambda wildcards: config['vcfs'][wildcards.vcf]
    output:
        temp('vcfs/{sample}.{vcf}.vcf.gz')
    resources:
        mem_mb = 4000
    shell:
         '''
         bcftools view -s {wildcards.sample} {input.vcf} | bcftools norm -m -any -f {config[reference]} | bcftools view -i 'GT[*]="alt"' -o {output}
         tabix -fp vcf {output}
         '''


def merfin_fitting(wildcards,lookup,model):
    if wildcards.fitted == 'fitted':
        return f'-prob {lookup} -peak {float([line for line in open(model)][7].split()[1])}'
    elif wildcards.fitted == 'coverage':
        return f'-peak {config["samples"][wildcards.sample]}'
    else:
        return ''

rule merfin_filter:
    input:
        seqmers = 'readmers/ARS.seqmers.meryl',
        readmers = 'readmers/{sample}.readmers.gt1.SR.meryl',
        vcf = 'vcfs/{sample}.{vcf}.vcf.gz', #lambda wildcards: config['vcfs'][wildcards.vcf],
        lookup = lambda wildcards: 'readmers/{sample}/lookup_table.txt' if wildcards.fitted == 'fitted' else [],
        model = lambda wildcards: 'readmers/{sample}/model.txt' if wildcards.fitted == 'fitted' else []
    output:
        temp('merfin/{sample}.{vcf}.merfin.{fitted}.filter.vcf')
    params:
        out = lambda wildcards, output: PurePath(output[0]).with_suffix('').with_suffix(''),
        fitted = lambda wildcards, input: merfin_fitting(wildcards,input.lookup,input.model)
    threads: 10
    resources:
        mem_mb = 5000,
        walltime = '4:00'
    shell:
        '''
        /cluster/work/pausch/alex/software/MERFIN_NEW/build/bin/merfin -filter -sequence {config[reference]} -seqmers {input.seqmers} {params.fitted} -readmers {input.readmers} -threads {threads} -vcf {input.vcf} -output {params.out}
        '''
        
rule bcftools_sort:
    input:
        'merfin/{sample}.{vcf}.merfin.{fitted}.filter.vcf'
    output:
        temp('merfin/{sample}.{vcf}.merfin.{fitted}.filter.vcf.gz')
    resources:
        mem_mb = 2000,
        disk_scratch = 2
    shell:
        '''
        bcftools sort -T $TMPDIR -o {output} {input}
        tabix -fp vcf {output}
        '''
        
rule tabulate_results:
    input:
        vcfs = expand('vcfs/{sample}.{vcf}.vcf.gz',sample=config['samples'],vcf=config['vcfs']),
        merfins = expand('merfin/{sample}.{vcf}.merfin.{fitted}.filter.vcf.gz',sample=config['samples'],vcf=config['vcfs'],fitted=config['modes'])
    output:
        filtering = 'merfin_filtering.csv',
        genotypes = 'merfin_genotypes.csv'
    resources:
        mem_mb = 2000
    params:
        samples = list(config['samples'].keys()),
        N_samples = len(config['samples'])-1,
        MAX_V = len(config['vcfs'])-1,
        MAX_M = len(list(product(config['vcfs'],config['modes'])))-1,
        orders = ','.join(list(config['vcfs'])+[f'{V}_{F}' for V,F in product(config['vcfs'],config['modes'])])
    shell:
        '''
        echo sample,{params.orders} > {output.filtering}
        echo sample,GT,count > {output.genotypes}
        vcf=({input.vcfs})
        merfins=({input.merfins})
        samples=({params.samples})
        for s in {{0..{params.N_samples}}}
        do
          echo -n ${{samples[$s]}}, >> {output.filtering}
          for X in {{0..{params.MAX_V}}}
          do
            echo -n $(bcftools index -n ${{vcf[$((s*({params.MAX_V}+1) + X))]}}), >> {output.filtering}
            bcftools query -f '[%GT]\n' ${{vcf[$((s*({params.MAX_V}+1) + X))]}} | sed -E 's/(1\|0|1\/0|1\|0|0\|1)/0\/1/g' | sed -E 's/1\|1/1\/1/g' | sort | uniq -c | awk -v S=${{vcf[$((s*({params.MAX_V}+1) + X))]}} '{{print S","$2","$1}}' >> {output.genotypes}
          done
          for Y in {{0..{params.MAX_M}}}
          do
            echo -n $(bcftools index -n ${{merfins[$((s*({params.MAX_M}+1) + Y))]}}), >> {output.filtering}
            bcftools query -f '[%GT]\n' ${{merfins[$((s*({params.MAX_M}+1) + Y))]}} | sed -E 's/(1\|0|1\/0|1\|0|0\|1)/0\/1/g' | sed -E 's/1\|1/1\/1/g' | sort | uniq -c | awk -v S=${{merfins[$((s*({params.MAX_M}+1) + Y))]}} '{{print S","$2","$1}}' >> {output.genotypes}
          done
          echo >> {output}
        sed -i 's/,$//g' {output.filtering}
        done
        '''

rule bcftools_isec:
    input:
        vcfs = expand('vcfs/{{sample}}.{vcf}.vcf.gz',vcf=config['vcfs']),
        filtered = expand('merfin/{{sample}}.{vcf}.merfin.{fitted}.filter.vcf.gz',vcf=config['vcfs'],fitted=config['modes'])
    output:
        'merfin/{sample}_merfin.intersections'
    threads: 4
    resources:
        mem_mb = 1000,
        walltime = '60'
    shell:
        '''
	    bcftools isec --threads {threads} {input} | awk '{{print $5,(length($3)+length($4))==2}}' | sort -k1,1 | uniq -c | awk '{{print $0" {wildcards.sample}"}}' > {output}
        '''

rule tabulate_isec:
    input:
        expand('merfin/{sample}_merfin.intersections',sample=config['samples'])
    output:
        'merfin_isec.upset'
    params:
        samples = list(config['samples'].keys())
    shell:
        '''
        cat {input} > {output}
        '''
