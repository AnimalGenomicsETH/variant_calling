from pathlib import Path, PurePath

rule all:
    input:
        expand('alignments/{sample}.{ext}',sample=config['samples'],ext=('bam','bam.csi'))

#ASL: -q 15 -u 40 are both defaults, so not needed.
rule fastp_filter:
    input:
        expand('input/{sample}.R{N}.fastq.gz',N=(1,2),allow_missing=True)
    output:
        fastq = expand('fastq/{sample}.R{N}.fastq.gz',N=(1,2),allow_missing=True)
    shell:
        '''
        fastp -q 15 -u 40 -g -i {input[0]} -o {output.fastq[0]} -I {input[1]} -O {output.fastq[1]} --json /dev/null --html /dev/null
        '''

rule bwamem2_index:
    input:
        reference = config['reference']
    output:
        temp(multiext(config['reference'],'.0123','.amb','.ann','.2bit.64','.pac'))
    threads: 1
    resources:
        mem_mb = 85000
    shell:
        'bwa-mem2 index {input}'
        
rule bwamem2_alignment:
    input:
        fastq = rules.fastp_filter.output,
        reference_index = rules.bwamem2_index.output
    output:
        bam = multiext('alignments/{sample}.bam','','.csi'),
        dedup_stats = 'alignments/{sample}.dedup.stats'
    params:
        bwa_index = lambda wildcards, input: PurePath(input.reference_index[0]).with_suffix(''),
        rg="TODO"#"@RG\\tID:{flowcell}.{lane}\\tCN:TUM\\tLB:{sample}\\tPL:illumina\\tPU:{flowcell}:{lane}\\tSM:{sample}",
    threads: 24
    resources:
        mem_mb = 4000,
        scratch = '50G',
        walltime = '24:00'
    shell:
        '''
        bwa-mem2 -Y -t {threads} {input.reference_index} {input.fastq} |\
        samtools collate -u -O |\
        samtools fixmate -m -u - - |\
        samtools sort -T $TMPDIR |\
        samtools markdup -T $TMPDIR --write-index -f {output.dedup_stats} - {output.bam[0]}
        '''
