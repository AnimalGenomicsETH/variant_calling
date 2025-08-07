rule VEP:
    input:
        vcf = multiext('{run}/contigs.beagle4.vcf.gz','','.tbi'),
        gtf = config['gtf'],
        reference = config['reference']
    output:
        vcf = multiext('{run}/all.beagle4.vep.vcf.gz','','.tbi'),
        stats = multiext('{run}/all.beagle4.vep','.stats.txt','.warnings')
    params:
        stats = lambda wildcards, output: Path(output.stats[0]).with_suffix('')
    threads: 8
    resources:
        mem_mb_per_cpu = 6000,
        runtime = '24h'
    container: '/cluster/work/pausch/alex/software/images/ensembl-vep_release_114.2.sif'
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
--stats_text --stats_file {params.stats} \
--warning_file {output.stats[1]}

tabix -p vcf {output.vcf[0]}
'''
