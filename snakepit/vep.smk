rule VEP:
    input:
        vcf = '',
        gtf = config['gtf'],
        reference = config['reference']
    output:
        vcf = multiext('','','.tbi'),
        stats = ''
    threads: 4
    resources:
        mem_mb_per_cpu = 1500,
        runtime = '4h'
    container: '/cluster/work/pausch/alex/software/images/ensembl-vep_release_113.0.sif'
    shell:
        '''
vep \
--fork {threads} \
-i {input.vcf[0]} \
-gtf {input.gtf} \
-fasta {input.reference} \
-o {output.vcf[0]} \
--hgvs --vcf --symbol \
--compress_output bgzip \
--stats_text --stats_file {output.stats}

tabix -p vcf {output.vcf[0]}
'''
