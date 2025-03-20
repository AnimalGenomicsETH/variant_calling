rule pbsv_discover:
    input:
        bam = multiext('alignments/{sample}.pbmm2.bam','','.csi'),
    output:
        'SVs/pbsv/{sample}.svsig.gz'
    conda:
        'pbccs'
    threads: 1
    resources:
        mem_mb_per_cpu = 2500
    shell:
        'pbsv discover --hifi {input[0]} {output}'

rule pbsv_call:
    input:
        signatures = expand(rules.pbsv_discover.output,sample=('BSWCHEF120141990854',))#config['samples'])
    output:
        'SVs/pbsv/cohort.SV.vcf'
    conda:
        'pbccs'
    threads: 8
    resources:
        mem_mb_per_cpu = 4000
    shell:
        'pbsv call --hifi -j {threads} {config[reference]} {input.signatures} {output}'

CHROMOSOMES = list(map(str,range(1,30))) + ['X','Y','MT']

from pathlib import Path
## implictly assume indexed bam/cram files, snakemake won't complain but jobs may fail.
def get_sample_bam_path(wildcards):
    bam = alignment_metadata.filter(pl.col('sample ID')==wildcards.sample).get_column('bam path').to_list()[0]
    match Path(bam).suffix:
        case '.bam':
            if Path(f"{bam}.bai").exists():
                return bam, f"{bam}.bai"
            elif Path(f"{bam}.csi").exists():
                return bam, f"{bam}.csi"
        case '.cram':
            if Path(f"{bam}.crai").exists():
                return bam, f"{bam}.crai"
    raise Exception(f"BAM file ({bam}) is not indexed.")

rule sawfish_discover:
    input:
        bam = get_sample_bam_path,
        reference = config['reference'],
        PAR = config['PAR'] if 'PAR' in config else []
    output:
        candidates = directory('SVs/sawfish/{sample}')
    params:
        regions = ','.join(CHROMOSOMES),
        PAR = f"--expected-cn {config['PAR']}" if 'PAR' in config else ''
    threads: 8
    resources:
        mem_mb_per_cpu = 7000,
        runtime = '4h'
    shell:
        '''
sawfish discover \
--threads {threads} \
--output-dir {output} \
--ref {input.reference} \
--bam {input.bam[0]} \
--target-region {params.regions} \
--cov-regex "^\d+" \
{params.PAR}
        '''

rule sawfish_joint_call:
    input:
        candidates = expand(rules.sawfish_discover.output['candidates'],sample=samples)
    output:
        vcf = directory('sawfish/SVs')
    params:
        files = lambda wildcards, input: ' '.join(f'--sample {S}' for S in input.candidates),
        regions = ','.join(CHROMOSOMES)
    threads: 8
    resources:
        mem_mb_per_cpu = 6000,
        runtime = '4h'
    shell:
        '''
sawfish joint-call \
--threads {threads} \
--output-dir {output.vcf} \
--target-region {params.regions} \
{params.files}
        '''

rule pbsv_predict_tandem_repeats:
    input:
        reference = config['reference']
    output:
        bed = 'GCA_002263795.4_ARS-UCD2.0_genomic.trf.bed'
    shell:
        '''
TODO: add pbsv script
        '''

rule sniffles2_call:
    input:
        bam = get_sample_bam_path,
        reference = config['reference'],
        TR = rules.pbsv_predict_tandem_repeats.output['bed'],
        SV_panel = lambda wildcards: rules.sniffles_filter.output if wildcards.mode == 'forced' else []
    output:
        vcf = 'SVs/sniffles2/{sample}.{mode}.vcf.gz',
        snf = 'SVs/sniffles2/{sample}.{mode}.snf'
    params:
        panel = lambda wildcards, input: '--genotype-vcf {input.SV_panel}' if wildcards.mode == 'forced'  else ''
    threads: 4
    resources:
        mem_mb_per_cpu = 2000
    conda:
        'sniffles'
    shell:
        '''
sniffles \
--input {input.bam[0]} \
--reference {input.reference} \
--tandem-repeats {input.TR} \
--sample-id {wildcards.sample} \
--threads {threads} \
--max-del-seq-len 100000 \
{params.panel} \
--vcf {output.vcf} \
--snf {output.snf}
        '''

rule sniffles_merge:
    input:
        snfs = expand(rules.sniffles2_call.output['snf'],sample=samples,allow_missing=True),
        reference = config['reference'],
        TR = rules.pbsv_predict_tandem_repeats.output['bed'],
    output:
        vcf = 'SVs/sniffles2/cohort.denovo.vcf.gz',
        filtered = 'SVs/sniffles2/cohort.denovo.filtered.vcf.gz'
    threads: 12
    resources:
        mem_mb_per_cpu = 3000
    conda:
        'sniffles'
    shell:
        '''
sniffles \
--input {input.snfs} \
--reference {input.reference} \
--tandem-repeats {input.TR} \
--threads {threads} \
--max-del-seq-len 100000 \
--vcf {output.vcf}

bcftools +fill-from-fasta {output.vcf} -- -c REF -f {input.reference} |\
bcftools view --threads {threads} -i 'abs(INFO/SVLEN)<=1000000&&INFO/SVTYPE!="BND"' -o {output.filtered} --write-index
        '''