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
        index = temp(multiext(config['reference'],'.0123','.amb','.ann','.bwt.2bit.64','.pac'))
    threads: 1
    resources:
        mem_mb = 85000
    shell:
        'bwa-mem2 index {input.reference}'

#Strobealign is quick to index on the fly and removes read-length dependency
rule strobealign_align:
    input:
        reference = config['reference'],
        fastq = rules.fastp_filter.output['fastq']
    output:
        sam = pipe('alignments/{sample}.strobe.sam')
    shell:
        '''
strobealign {input.reference} {input.fastq} -t {threads} > {output.sam}
        '''

rule bwamem2_align:
    input:
        index = rules.bwamem2_index.output['index'],
        fastq = rules.fastp_filter.output['fastq']
    output:
        sam = pipe('alignments/{sample}.bwa.sam')
    shell:
        '''
bwa-mem2 mem -Y -t {threads} {input.index} {input.fastq} > {output.sam}
        '''

#NOTE: markdup -S flag marks supplementary alignments of duplicates as duplicates.
#NOTE: readgroup as an additional tag is seemingly unused now, so removed.
#This now selects the alignment code using the "generate_aligner_command" function, and then pipes through samtools.
rule samtools_markdup:
    input:
        sam = 'alignments/{sample}.{aligner}.sam'
    output:
        bam = 'alignments/{sample}.{aligner}.{EXT,cram|bam}',
        dedup_stats = 'alignments/{sample}.{aligner}.{EXT}.dedup.stats'
    params:
        cram = lambda wildcards, input: f'--output-fmt-option version=3.1  --reference {input.reference}' if wildcards.EXT == 'cram' else ''
    threads: 2
    resources:
        mem_mb_per_cpu = 1000,
        runtime = '1h'
    shell:
        '''
samtools collate -@ {threads} -O -u {input.sam} |\
samtools fixmate -@ {threads} -m -u /dev/stdin /dev/stdout |\
samtools sort -@ {threads} -T $TMPDIR -u |\
samtools markdup -@ {threads} -T $TMPDIR -S --write-index -f {output.dedup_stats} {params.cram} /dev/stdin {output.bam}
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

#samtools view -e '![vW] || [vW]==1' -o {output[0]} --write-index -t {threads} {output.bam}
rule STAR_align:
    input:
        fastq = rules.fastp_filter.output['fastq'],
        index = rules.STAR_index.output['index'],
        vcf = config['VCF']
    output:
        sam = pipe('alignments/{sample}.STAR.sam')
    threads: 12
    resources:
        mem_mb_per_cpu = 4000,
        runtime = "4h"
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
--outStd SAM > {output.sam}
        '''

rule paftools_gff2bed:
    input:
        GTF = config['GTF']
    output:
        bed = 'mm2.bed' #rename eventually
    localrule: True
    shell:
        '''
paftools.js gff2bed {input.GTF} > {output.bed}
        '''

rule minimap2_align_RNA:
    input:
       fastq = rules.fastp_filter.output['fastq'],
       reference = config['reference'],
       bed = rules.paftools_gff2bed.output['bed']
    output:
        sam = pipe('alignments/{sample}.mm2.sam')
    shell:
        '''
minimap2 -t {threads} -x splice:sr -j {input.bed} --write-junc {input.reference} {input.fastq} > $TMPDIR.bed
minimap2 -t {threads} -ax splice:sr -j {input.bed} --pass1=$TMPDIR.bed {input.reference} {input.fastq} > {output.sam}
        '''
