wildcard_constraints:
    entry = r'|'.join(config['vcf']),
    truth = r'|'.join(config['vcf']),
    query = r'|'.join(config['vcf'])


rule all:
    input:
        expand('F1/{truth}_{query}.tsv',truth='v18',query='small')

checkpoint bcftools_split:
    input:
        vcf = lambda wildcards: config['vcf'][wildcards.entry]
    output:
        _dir = directory('F1/split_{entry}')
    resources:
        mem_mb_per_cpu = 2500
    shell:
        '''
        bcftools +split -i 'GT[*]="alt"' -Ob -o {output._dir} {input.vcf}
        '''

#--stratification {input.regions} \
rule happy:
    input:
        truth = 'F1/split_{truth}/{sample}.bcf',
        query = 'F1/split_{query}/{sample}.bcf',
        reference = config['reference']
    output:
        csv = 'F1/{truth}_{query}.{sample}.summary.csv',
        others = temp(multiext('F1/{truth}_{query}.{sample}','.bcf','.bcf.csi','.extended.csv','.roc.all.csv.gz','.runinfo.json'))
    params:
        _dir = lambda wildcards, output: Path(output.csv).with_suffix('').with_suffix('')
    container: '/cluster/work/pausch/alex/software/images/hap.py_latest.sif'
    threads: 2
    resources:
        mem_mb_per_cpu = 2500
    shell:
        '''
/opt/hap.py/bin/hap.py \
--reference {input.reference} \
--bcf \
--usefiltered-truth \
--no-roc \
--no-json \
--leftshift \
--pass-only \
--scratch-prefix $TMPDIR \
--write-counts \
--threads {threads} \
--report-prefix {params._dir} \
{input.truth} {input.query}
        '''

# ugly solution for https://github.com/snakemake/snakemake/issues/3475
def aggregate_happy(wildcards):
    query_samples = glob_wildcards(Path(checkpoints.bcftools_split.get(entry=wildcards.truth).output['_dir']).joinpath('{sample}.bcf')).sample
    return expand(rules.happy.output['others'][2],sample=query_samples,allow_missing=True)

def agg_hap2(wildcards):
    truth_samples = glob_wildcards(Path(checkpoints.bcftools_split.get(entry=wildcards.query).output['_dir']).joinpath('{sample}.bcf')).sample
    return []

rule gather_happy:
    input:
        aggregate_happy,
        agg_hap2
    output:
        csv = 'F1/{truth}_{query}.tsv'
    localrule: True
    shell:
        '''
        echo -e "variant truth query recall precision truth_TiTv query_TiTv F1_score sample tissue" > {output}
        for i in {input}
        do
          awk -v I=$(basename $i) -F',' '$2=="*"&&$4=="PASS" {{ split(I,a,"."); print $1,$17,$38,$8,$9,$22,$43,$11,a[1],a[2] }}' $i >> {output}
        done
        '''
