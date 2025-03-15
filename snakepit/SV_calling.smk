rule pbsv_discover:
    input:
        bam = multiext('alignments/{sample}.pbmm2.bam','','.csi'),
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

CHROMOSOMES = list(map(str,range(1,30))) + ['X','Y','MT']

rule sawfish_discover:
    input:
        bam = multiext('alignments/{sample}.pbmm2.bam','','.csi'),
        reference = config['reference']
    output:
        candidates = directory('SVs/sawfish/{sample}')
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
        candidates = expand(rules.sawfish_discover.output['candidates'],sample=samples)
    output:
        directory('sawfish/SVs')
    params:
        files = lambda wildcards, input: ' '.join(f'--sample {S}' for S in input.candidates),
        regions = ','.join(CHROMOSOMES)
    threads: 8
    resources:
        mem_mb_per_cpu = 6000,
        runtime = '4h'
    shell:
        '''
sawfish joint-call \
--threads {threads} \
--output-dir {output.vcf} \
--target-region {params.regions} \
{params.files}
        '''
