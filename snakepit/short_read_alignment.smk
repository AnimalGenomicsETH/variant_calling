from pathlib import Path

import polars as pl

short_read_data = pl.read_csv(config['short_reads']).rows_by_key(key=["sample"]) if 'short_reads' in config else pl.DataFrame()

rule all:
    input:
        expand('alignments/{sample}.{caller}.{ext}',sample=short_read_data,ext='bam',caller=config.get('aligners','bwa'))

rule fastp_filter:
    input:
        fastq = lambda wildcards: short_read_data[wildcards.sample][0]
    output:
        fastq = expand('fastq/{sample}.R{N}.fastq.gz',N=(1,2),allow_missing=True)
    params:
        min_quality = 15,
        unqualified = 40,
        min_length  = 50,
        poly_x = '' if config.get('RNA',False) else '--trim_poly_x'
    threads: 4
    resources:
        mem_mb_per_cpu = 2500,
        runtime = '4h'
    shell:
        '''
fastp \
--dont_eval_duplication \
--qualified_quality_phred {params.min_quality} \
--unqualified_percent_limit {params.unqualified} \
--trim_poly_g  {params.poly_x} \
--length_required {params.min_length} \
--thread {threads} \
--in1 {input.fastq[0]} \
--out1 {output.fastq[0]} \
--in2 {input.fastq[1]} \
--out2 {output.fastq[1]} \
--json /dev/null \
--html /dev/null
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
        'bwa-mem2 index {input.reference}'

#ASL: fixed read size (-r) at 150, may need changing for different batches
rule strobealign_index:
    input:
        reference = config['reference']
    output:
        temp(config['reference'] + '.sti')
    threads: 4
    resources:
        mem_mb = 8000
    shell:
        'strobealign {input} --create-index -t {threads} -r 150'

def generate_aligner_command(aligner,input,threads):
    index = Path(input.reference_index[0]).with_suffix('')
    match aligner:
        case 'strobe':
            return f'strobealign {index} {input.fastq} --use-index -t {threads} -r 150'
        case 'bwa':
            return f'bwa-mem2 mem -Y -t {threads} {index} {input.fastq}'
        case _:
            raise('Unknown aligner')

#NOTE: markdup -S flag marks supplementary alignments of duplicates as duplicates.
#NOTE: readgroup as an additional tag is seemingly unused now, so removed.
#This now selects the alignment code using the "generate_aligner_command" function, and then pipes through samtools.
rule aligner:
    input:
        fastq = rules.fastp_filter.output,
        reference_index = lambda wildcards: rules.strobealign_index.output if wildcards.aligner == 'strobe' else rules.bwamem2_index.output,
        reference = lambda wildcards: [] if wildcards.EXT == 'bam' else config['reference']
    output:
        bam = 'alignments/{sample}.{aligner}.{EXT,cram|bam}',
        dedup_stats = 'alignments/{sample}.{aligner}.{EXT}.dedup.stats'
    params:
        aligner_command = lambda wildcards, input, threads: generate_aligner_command(wildcards.aligner,input,threads),
        cram = lambda wildcards, input: f'--output-fmt-option version=3.1  --reference {input.reference}' if wildcards.EXT == 'cram' else ''
    threads: 16
    resources:
        mem_mb_per_cpu = 4000,
        runtime = '24h'
    shell:
        '''
{params.aligner_command} |\
samtools collate -u -O -@ 6 - |\
samtools fixmate -m -u -@ 6 - - |\
samtools sort -T $TMPDIR -u -@ 6 |\
samtools markdup -T $TMPDIR -S -@ 6 --write-index -f {output.dedup_stats} {params.cram} - {output.bam[0]}
        '''

rule STAR_index:
    input:
        reference = config['reference'],
        GTF = config['GTF']
    output:
        index = directory('STAR/reference')
    threads: 4
    resources:
        mem_mb_per_cpu = 15000,
        runtime = "4h"
    shell:
        '''
STAR \
--runMode genomeGenerate \
--genomeFastaFiles {input.reference} \
--sjdbGTFfile {input.GTF} \
--outTmpDir $TMPDIR/STAR \
--genomeDir {output.index} \
--runThreadN {threads} \
--sjdbOverhang 100
        '''

#--sjdbGTFfile {input.GTF} \
#--sjdbOverhang 100 \
#--outSAMattrRGline {params.rg} \
#--outFileNamePrefix $TMPDIR \
rule STAR_align:
    input:
        fastq = rules.fastp_filter.output['fastq'],
        index = rules.STAR_index.output['index'],
        vcf = config['VCF']
    output:
        bam = multiext('alignments/{sample}.STAR.bam','','.csi'),
        dedup_stats = 'alignments/{sample}.STAR.dedup.stats'
    threads: 12
    resources:
        mem_mb_per_cpu = 7000,
        runtime = "24h"
    shell:
        '''
bcftools view \
--samples {wildcards.sample} \
--types snps \
--genotype het \
--threads {threads} \
{input.vcf} > $TMPDIR/heterozygotes.vcf

STAR \
--runMode alignReads \
--twopassMode Basic \
--genomeDir {input.index} \
--readFilesIn {input.fastq} \
--readFilesCommand zcat \
--varVCFfile $TMPDIR/heterozygotes.vcf \
--outTmpDir $TMPDIR/STAR \
--runThreadN {threads} \
--outSAMmapqUnique 60 \
--waspOutputMode SAMtag \
--outMultimapperOrder Random \
--outFileNamePrefix $TMPDIR \
--outSAMtype SAM \
--outStd SAM |\
samtools collate -u -O -@ 6 - |\
samtools fixmate -m -u -@ 6 - - |\
samtools sort -T $TMPDIR -u -@ 6 |\
samtools markdup -T $TMPDIR -S -@ 6 --write-index -f {output.dedup_stats} - {output.bam[0]}

#samtools view -e '![vW] || [vW]==1' -o {output[0]} --write-index -t {threads} {output.bam}
        '''
