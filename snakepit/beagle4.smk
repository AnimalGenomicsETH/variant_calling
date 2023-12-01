
from pathlib import PurePath

regions = [line.strip().split()[0] for line in open('all_regions.bed')]

rule all:
    input:
        expand('variants/contigs/{region}.beagle4.vcf.gz',region=regions)

rule bcftools_scatter:
    input:
        gvcf = 'variants/contigs.Unrevised.vcf.gz'
    output:
        expand('variants/contigs/{region}.vcf.gz',region=regions,allow_missing=True),
        expand('variants/contigs/{region}.vcf.gz.csi',region=regions,allow_missing=True)
    params:
        _dir = lambda wildcards, output: PurePath(output[0]).with_suffix('').with_suffix('').with_suffix('').with_suffix('')
    threads: 2
    resources:
        mem_mb = 2500,
        walltime = '4h'
    shell:
        '''
        bcftools +scatter {input.gvcf} -o {params._dir} -Oz --threads {threads} --write-index -S <(cut -f 1 all_regions.bed) --no-version
        '''

rule beagle4_imputation:
    input:
        'variants/contigs/{region}.vcf.gz'
    output:
        'variants/contigs/{region}.beagle4.vcf.gz'
    params:
        prefix = lambda wildcards, output: PurePath(output[0]).with_suffix('').with_suffix(''),
        name = lambda wildcards, output: PurePath(output[0]).name
    threads: lambda wildcards, input: 16 if input.size_mb > 10 else (8 if input.size_mb > 1 else 2)
    resources:
        mem_mb = 4000,
        walltime = lambda wildcards, input: '24h' if input.size_mb > 50 else '4h'
    shell:
        '''
        java -jar -Xss25m -Xmx60G /cluster/work/pausch/alex/software/beagle.27Jan18.7e1.jar gl={input} nthreads={threads} out={params.prefix}
        mv {output[0]} $TMPDIR/{params.name}
        bcftools reheader -f {config[reference]}.fai -o {output[0]} $TMPDIR/{params.name}
        tabix -fp vcf {output[0]}
        '''
