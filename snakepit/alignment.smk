configfile: "config.yaml"
sample = config["samples"]
reference = config["references"]
chr_list = list(range(1,30))

# Reduced paths
raw_data_path = config["resources"]["raw"]
fastp_output = config["fold_out"]["fastp"]
split_fastq = config["fold_out"]["split_fastq"]
alignment = config["fold_out"]["alignment"]
sorted_alignment = config["fold_out"]["sorted_alignment"]
dedup_alignment = config["fold_out"]["dedup_alignment"]
picard_metrics = config["fold_out"]["picard_metrics"]
reference_path = config["resources"]["reference"]
map_stats = config["fold_out"]["map_stats"]


# Tools
FASTP = config["tools"]["fastp"]
LOAD_PYTHON_374 = config["tools"]["load_python_374"]
FASTQ_SPLITTER = config["tools"]["fastq_splitter"]
LOAD_BWA = config["tools"]["load_bwa"]
LOAD_SAMTOOLS = config["tools"]["load_samtools"]
SAMBLASTER = config["tools"]["samblaster"]
SAMBAMBA = config["tools"]["sambamba"]
LOAD_PICARD_TOOLS = config["tools"]["load_picard_tools"]
LOAD_JAVA = config["tools"]["load_java_jdk"]
LOAD_R = config["tools"]["load_r"]
PICARD_TOOLS = config["tools"]["picard_tools"]
MQ_sam_flags_script = config["tools"]["sam_flags_script"]

rule all:
    input:
        stats_1=expand(sorted_alignment + "{reference}/{sample}/{sample}.stats",sample=sample,reference=reference),
        stats_2=expand(dedup_alignment + "{reference}/{sample}/{sample}.stats",sample=sample,reference=reference),
        index=expand(dedup_alignment + "{reference}/{sample}/{sample}.bam.bai",sample=sample,reference=reference),
        metrics = expand(picard_metrics + "{reference}/{sample}/{sample}.alignment_summary_metrics",sample=sample,reference=reference),
        stats = expand(picard_metrics + "{reference}/{sample}/{sample}_idxstats",sample=sample,reference=reference),
        map_stats = expand(map_stats + "success/{reference}_{sample}_{chrm}.log",reference=reference,sample=sample,chrm=chr_list)

#ASL: -q 15 -u 40 are both defaults, so not needed.
rule fastp_filter:
    input:
        expand('reads_R{N}.fastq.gz',N=(1,2))
    output:
        fastq = expand('reads_R{N}.fastq.gz',N=(1,2))
        json = "{sample}/{sample}_fastp.json"
    shell:
        '''
        fastp -q 15 -u 40 -g -i {input[0]} -o {output.fastq[0]} -I {input[1]} -O {output.fastq[1]} -j {output.json} > /dev/null
        '''

checkpoint split_fastq:
    input:
        R1 = fastp_output + "{sample}/{sample}_R1.fastq.gz",
        R2 = fastp_output + "{sample}/{sample}_R2.fastq.gz"
    output:
        split_fastqs = temp(directory(split_fastq + "{sample}/"))
    params:
        prefix = split_fastq + "{sample}/{sample}_"
    shell:
        LOAD_PYTHON_374 + FASTQ_SPLITTER + " -o {params.prefix} {input.R1} {input.R2}"


rule bwamem2_index:
    input:
        reference = config['references']
    output:
        temp(get_dir('input','{ref}_index.0123'))
    params:
        lambda wildcards, output: PurePath(output[0]).with_suffix('')
    threads: 1
    resources:
        mem_mb = 85000
    shell:
        'bwa-mem2 index -p {params} {input}'
        
rule bwamem2_alignment:
	input:
		R1 = split_fastq + "{sample}/{sample}_{flowcell}_{lane}_R1.fq.gz",
		R2 = split_fastq + "{sample}/{sample}_{flowcell}_{lane}_R2.fq.gz",
		reference_index = rules.bwamem2_index.output
	output:
		bam = multiext('alignments/{sample}.bam','','.csi'),
        dedup_stats = 'alignments/{sample}.dedup.stats'
	params:
		rg="apple"#"@RG\\tID:{flowcell}.{lane}\\tCN:TUM\\tLB:{sample}\\tPL:illumina\\tPU:{flowcell}:{lane}\\tSM:{sample}",
	shell:
        '''
        bwa-mem2 -M -t {threads} {input.reference_index} {input.fastq} |\
        samtools collate -u -O |\
        samtools fixmate -m -u - - |\
        samtools sort -T $TMPDIR |\
        samtools markdup -T $TMPDIR --write-index -f {output.dedup_stats} -o {output.bam[0]}
        '''




## merge
## duplicates
## stats
