from pathlib import Path, PurePath

def get_bam_extensions():
    if config.get('cram',False):
        return ('cram','cram.crai')
    else:
        return ('bam','bam.csi')

rule all:
    input:
        expand('alignments/{sample}.{caller}.{ext}',sample=config['samples'],ext=get_bam_extensions(),caller=config['aligners'])

rule fastp_filter:
    input:
        expand('input/{sample}.R{N}.fastq.gz',N=(1,2),allow_missing=True)
    output:
        fastq = expand('fastq/{sample}.R{N}.fastq.gz',N=(1,2),allow_missing=True)
    params: #Fancy way of handling if there is no fastp values in the configfile
        min_quality = config.get('fastp',{}).get('min_quality',15),
        unqualified = config.get('fastp',{}).get('unqualified',40),
        min_length  = config.get('fastp',{}).get('min_length',15)
    threads: 4
    resources:
        mem_mb = 2500
    shell:
        '''
        fastp -q {params.min_quality} -u {params.unqualified} -g --length_required {params.min_length} --thread {threads} -i {input[0]} -o {output.fastq[0]} -I {input[1]} -O {output.fastq[1]} --json /dev/null --html /dev/null
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

#ASL: fixed read size (-r) at 150, may need changing for different batches
rule strobealign_index:
    input:
        reference = config['reference']
    output:
        temp(config['reference'] + '.sti')
    threads: 1
    resources:
        mem_mb = 30000
    shell:
        'strobealign {input} -i -r 150'

def generate_aligner_command(aligner):
    match aligner:
        case 'strobe':
            return 'strobealign {params.strobe_index} {input.fastq} --use-index -t {threads} -r 150'
        case 'bwa':
            return 'bwa-mem2 mem -Y -t {threads} {params.bwa_index} {input.fastq}'
        case _:
            raise('Unknown aligner')

#NOTE: markdup -S flag marks supplementary alignments of duplicates as duplicates.
#NOTE: readgroup as an additional tag is seemingly unused now, so removed.
#This now selects the alignment code using the "generate_aligner_command" function, and then pipes through samtools.
rule aligner:
    input:
        fastq = rules.fastp_filter.output,
        reference_index = lambda wildcards: rules.strobealign_index.output if wildcards.aligner == 'strobe' else rules.bwamem2_index.output,
        reference = config['reference']
    output:
        bam = expand('alignments/{sample}.{aligner}.{ext}',ext=get_bam_extensions(),allow_missing=True),
        dedup_stats = 'alignments/{sample}.{aligner}.dedup.stats'
    params:
        aligner_command = lambda wildcards: generate_aligner_command(wildcards.aligner),
        index = lambda wildcards, input: PurePath(input.reference_index[0]).with_suffix('')
    threads: 12
    resources:
        mem_mb = 5000,
        scratch = '50G'
    shell:
        '''
        {params.aligner_command} |\
        samtools collate -u -O -@ 6 - |\
        samtools fixmate -m -u -@ 6 - - |\
        samtools sort -T $TMPDIR -u -@ 6 |\
        samtools markdup -T $TMPDIR -S -@ 6 --write-index -f {output.dedup_stats} --reference {input.reference} - {output.bam[0]}
        '''
