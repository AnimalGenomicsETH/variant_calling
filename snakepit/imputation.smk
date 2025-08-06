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

rule beagle_imputation:
    input:
        vcf = multiext('{run}/{region}.Unrevised.vcf.gz','','.tbi'),
        fai = f"{config['reference']}.fai"
    output:
        vcf = multiext('{run}/{region}.{beagle,beagle4|beagle5.vcf.gz','','.tbi')
    params:
        ne = 100, #need to use?
        imputation_tag = 'gl' if wildcards.beagle == 'beagle4' else 'gt'
    threads: 6
    resources:
        mem_mb_per_cpu = 8000,
        runtime = lambda wildcards, input: '4h' if input.size_mb > 5 else '1h'
    shell:
        '''
        #if there are no variants, beagle errors out ungracefully
        if [[ $(bcftools index -n {input.vcf[0]}) -gt 0 ]]
        then
          {wildcards.beagle} {params.imputation_tag}={input.vcf[0]} nthreads={threads} out=$TMPDIR/imputed
          bcftools reheader -f {input.fai} -o {output.vcf[0]} $TMPDIR/imputed.vcf.gz
        else
          cp {input.vcf[0]} {output.vcf[0]}
        fi
        tabix -fp vcf {output.vcf[0]}
        '''

#bcftools concat --threads 4 --naive-force -o variants_2025/contigs.beagle4.vcf.gz variants_2025/*.beagle4.vcf.gz
