
rule pbsv_discover:
    input:
        expand(rules.samtools_merge.output,mapper='mm2',allow_missing=True)
    output:
        'variants/SVs/{sample}.svsig.gz'
    conda:
        'pbccs'
    threads: 1
    resources:
        mem_mb_per_cpu = 2500
    shell:
        'pbsv discover --hifi {input[0]} {output}'

rule pbsv_call:
    input:
        signatures = expand(rules.pbsv_discover.output,sample=('BSWCHEF120141990854',))#config['samples'])
    output:
        'variants/SVs/cohort.SV.vcf'
    conda:
        'pbccs'
    threads: 8
    resources:
        mem_mb_per_cpu = 4000
    shell:
        'pbsv call --hifi -j {threads} {config[reference]} {input.signatures} {output}'


CHROMOSOMES = list(map(str,range(1,30)))+['X','Y','MT']

rule sawfish_discover:
    input:
        bam = multiext('alignments/{sample}.pbmm2.bam','','.csi'),
        reference = config['reference']
    output:
        directory('SVs/sawfish/{sample}')
    params:
        regions = ','.join(CHROMOSOMES)
    threads: 12
    resources:
        mem_mb_per_cpu = 7000,
        runtime = '4h'
    shell:
        '''
sawfish discover \
--threads {threads} \
--output-dir {output} \
--ref {input.reference} \
--bam {input.bam[0]} \
--target-region {params.regions} \
--cov-regex "^\d+" \
--expected-cn {input.PAR}
        '''

#TODO: Push samples to main
rule sawfish_joint_call:
    input:
        bam = expand(rules.sawfish_discover.output,sample=metadata.get_column('Interbull ID').unique().to_list())
    output:
        directory('sawfish/SVs')
    params:
        files = lambda wildcards, input: ' '.join(f'--sample {S}' for S in input.bam),
        regions = ','.join(CHROMOSOMES)
    threads: 8
    resources:
        mem_mb_per_cpu = 6000,
        runtime = '4h'
    shell:
        '''
sawfish joint-call \
--threads {threads} \
--output-dir {output} \
{params.files} \
--target-region {params.regions}
        '''
