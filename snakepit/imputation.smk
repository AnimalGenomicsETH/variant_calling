rule bcftools_scatter_X:
    input:
        gvcf = 'variants/contigs.Unrevised.vcf.gz'
    output:
        expand('variants/contigs/{region}.vcf.gz',region=regions,allow_missing=True),
        expand('variants/contigs/{region}.vcf.gz.csi',region=regions,allow_missing=True)
    params:
        _dir = lambda wildcards, output: PurePath(output[0]).with_suffix('').with_suffix('').with_suffix('').with_suffix('')
    threads: 2
    resources:
        mem_mb_per_cpu = 2500,
        runtime = '4h'
    shell:
        '''
        bcftools +scatter {input.gvcf} -o {params._dir} -Oz --threads {threads} --write-index -S <(cut -f 1 all_regions.bed) --no-version
        '''

#TODO: can we chain filtering steps?
rule bcftools_filter:
    input:
        vcf = lambda wildcards: multiext('{run}/{region}.Unrevised.vcf.gz','','.tbi') if wildcards.filtration == 'pre_filter' else multiext('{run}/{region}.beagle4.vcf.gz','','.tbi') 
    output:
        vcf = multiext('{run}/{region}.{filtration}.vcf.gz','','.tbi')
    shell:
        '''
bcftools view -i 'QUAL>20' -o {output} {input}
        '''    

rule beagle4_imputation:
    input:
        vcf = multiext('{run}/{region}.Unrevised.vcf.gz','','.tbi'),
        fai = f"{config['reference']}.fai"
    output:
        multiext('{run}/{region}.beagle4.vcf.gz','','.tbi')
    threads: 6
    resources:
        mem_mb_per_cpu = 4000,
        runtime = lambda wildcards, input: '24h' if input.size_mb > 50 else '4h'
    shell:
        '''
        #add ne for popsize
        java -jar -Xss25m -Xmx60G /cluster/work/pausch/alex/software/beagle.27Jan18.7e1.jar gl={input} nthreads={threads} out=$TMPDIR/imputed.vcf.gz
        bcftools reheader -f {input.fai} -o {output[0]} $TMPDIR/imputed.vcf.gz
        tabix -fp vcf {output[0]}
        '''

#beagle5