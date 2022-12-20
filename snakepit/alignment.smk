from pathlib import Path, PurePath

def get_bam_extensions():
    if config.get('cram',False):
        return ('cram','cram.crai')
    else:
        return ('bam','bam.csi')

rule all:
    input:
        expand('alignments/{sample}.{ext}',sample=config['samples'],ext=get_bam_extensions())

#ASL: -q 15 -u 40 are both defaults, so not needed.
rule fastp_filter:
    input:
        expand('input/{sample}.R{N}.fastq.gz',N=(1,2),allow_missing=True)
    output:
        fastq = expand('fastq/{sample}.R{N}.fastq.gz',N=(1,2),allow_missing=True)
    threads: 4
    resources:
        mem_mb = 2500
    shell:
        '''
        fastp -q 15 -u 40 -g --thread {threads} -i {input[0]} -o {output.fastq[0]} -I {input[1]} -O {output.fastq[1]} --json /dev/null --html /dev/null
        '''

rule bwamem2_index:
    input:
        reference = config['reference']
    output:
        temp(multiext(config['reference'],'.0123','.amb','.ann','.bwt.2bit.64','.pac'))
    threads: 1
    resources:
        mem_mb = 85000
    shell:
        'bwa-mem2 index {input}'

#NOTE: markdup -S flag marks supplementary alignments of duplicates as duplicates.
#NOTE: readgroup as an additional tag is seemingly unused now, so removed.
rule bwamem2_alignment:
    input:
        fastq = rules.fastp_filter.output,
        reference_index = rules.bwamem2_index.output
    output:
        bam = expand('alignments/{sample}.{ext}',ext=get_bam_extensions(),allow_missing=True),
        dedup_stats = 'alignments/{sample}.dedup.stats'
    params:
        bwa_index = lambda wildcards, input: PurePath(input.reference_index[0]).with_suffix('')
    threads: 24 #NOTE: samtools pipes are hardcoded using 4 threads (-@ 4).
    resources:
        mem_mb = 4000,
        scratch = '50G',
        walltime = '24:00'
    shell:
        '''
        bwa-mem2 mem -Y -t {threads} {params.bwa_index} {input.fastq} |\
        samtools collate -u -O -@ 4 - |\
        samtools fixmate -m -u -@ 4 - - |\
        samtools sort -T $TMPDIR -u -@ 4 |\
        samtools markdup -T $TMPDIR -S -@ 4 --write-index -f {output.dedup_stats} - {output.bam[0]}
        '''
