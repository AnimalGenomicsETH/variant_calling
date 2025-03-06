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
        vcf = 'variants/contigs/{region}.vcf.gz',
        fai = f"{config['reference']}.fai"
    output:
        multiext('variants/contigs/{region}.beagle4.vcf.gz','','.tbi')
    threads: 6
    resources:
        mem_mb = 4000,
        walltime = lambda wildcards, input: '24h' if input.size_mb > 50 else '4h'
    shell:
        '''
        #add ne for popsize
        java -jar -Xss25m -Xmx60G /cluster/work/pausch/alex/software/beagle.27Jan18.7e1.jar gl={input} nthreads={threads} out=$TMPDIR/imputed.vcf.gz
        bcftools reheader -f {input.fai} -o {output[0]} $TMPDIR/imputed.vcf.gz
        tabix -fp vcf {output[0]}
        '''

#beagle5