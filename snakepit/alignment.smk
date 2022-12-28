from pathlib import Path, PurePath

def get_bam_extensions():
    if config.get('cram',False):
        return ('cram','cram.crai')
    else:
        return ('bam','bam.csi')

rule all:
    input:
        expand('alignments/{sample}.{caller}{ext}',sample=config['samples'],ext=get_bam_extensions(),caller=('','strobe.'))

#ASL: -q 15 -u 40 are both defaults, so not needed.
#--length_required 15
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
        fastp -q 15 -u 40 -g --length_required 30 --thread {threads} -i {input[0]} -o {output.fastq[0]} -I {input[1]} -O {output.fastq[1]} --json /dev/null --html /dev/null
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
        reference_index = rules.bwamem2_index.output,
        reference = config['reference']
    output:
        bam = expand('alignments/{sample}.{ext}',ext=get_bam_extensions(),allow_missing=True),
        dedup_stats = 'alignments/{sample}.dedup.stats'
    params:
        bwa_index = lambda wildcards, input: PurePath(input.reference_index[0]).with_suffix(''),
        cram_options = '--output-fmt-option version=3.0 --output-fmt-option normal' if config.get('cram',False) else ''
    threads: 18 #NOTE: samtools pipes are hardcoded using 6 threads (-@ 6).
    resources:
        mem_mb = 4000,
        scratch = '50G',
        walltime = '4:00'
    shell:
        '''
        bwa-mem2 mem -Y -t {threads} {params.bwa_index} {input.fastq} |\
        samtools collate -u -O -@ 6 - |\
        samtools fixmate -m -u -@ 6 - - |\
        samtools sort -T $TMPDIR -u -@ 6 |\
        samtools markdup -T $TMPDIR -S -@ 6 --write-index {params.cram_options} -f {output.dedup_stats} --reference {input.reference} - {output.bam[0]}
        '''

rule strobealign_index:
    input:
        reference = config['reference']
    output:
        temp(config['reference'] + '.sti')
    threads: 1
    resources:
        mem_mb = 30000
    shell:
        '''
        strobealign {input} -i -r 150
        '''

rule strobealign_alignment:
    input:
        fastq = rules.fastp_filter.output,
        reference_index = rules.strobealign_index.output,
        reference = config['reference']
    output:
        bam = expand('alignments/{sample}.strobe.{ext}',ext=get_bam_extensions(),allow_missing=True),
        dedup_stats = 'alignments/{sample}.strobe.dedup.stats'
    params:
        strobe_index = lambda wildcards, input: PurePath(input.reference_index[0]).with_suffix('')
    threads: 12
    resources:
        mem_mb = 5000,
        scratch = '50G'
    shell:
        '''
        strobealign {params.strobe_index} {input.fastq} --use-index -t {threads} -r 150 |\
        samtools collate -u -O -@ 6 - |\
        samtools fixmate -m -u -@ 6 - - |\
        samtools sort -T $TMPDIR -u -@ 6 |\
        samtools markdup -T $TMPDIR -S -@ 6 --write-index -f {output.dedup_stats} --reference {input.reference} - {output.bam[0]}
        '''
