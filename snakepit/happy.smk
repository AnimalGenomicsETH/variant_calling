from pathlib import PurePath

wildcard_constraints:
    entry = r'|'.join(config['vcf']),
    truth = r'|'.join(config['vcf']),
    query = r'|'.join(config['vcf'])

rule all:
    input:
        expand('F1/{truth}_{query}.csv',truth='v18',query='small')

checkpoint bcftools_split:
    input:
        vcf = lambda wildcards: config['vcf'][wildcards.entry]
    output:
        _dir = directory('F1/split_{entry}')
    resources:
        mem_mb_per_cpu = 2500
    shell:
        '''
        bcftools +split -i 'GT[*]="alt"' -Oz -o {output._dir} {input.vcf}
        '''

#--stratification {input.regions} \
rule happy:
    input:
        truth = 'F1/split_{truth}/{sample}.vcf.gz',
        query = 'F1/split_{query}/{sample}.vcf.gz',
        reference = config['reference']
    output:
        csv = 'F1/{truth}_{query}.{sample}.summary.csv',
        others = temp(multiext('F1/{truth}_{query}.{sample}','.bcf','.bcf.csi','.extended.csv','.roc.all.csv.gz','.runinfo.json'))
    params:
        _dir = lambda wildcards, output: PurePath(output.csv).with_suffix('').with_suffix('')
    container: '/cluster/work/pausch/alex/software/images/hap.py_latest.sif'
    threads: 1
    resources:
        mem_mb = 5000,
        scratch = '10G'
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
{input.vcf_truth} {input.vcf_query}
        '''

def aggregate_happy(wildcards):
    query_samples = glob_wildcards(PurePath(checkpoints.bcftools_split.get(entry=wildcards.query).output['_dir']).joinpath('{sample}.vcf.gz')).sample
    
    #truth_samples = glob_wildcards(PurePath(checkpoints.bcftools_split.get(entry=wildcards.truth).output['_dir']).joinpath('{sample}.vcf.gz')).sample
    #print(truth_samples,query_samples)
    
    targets = []
    for sample in query_samples:
        targets.extend([f'F1/split_{wildcards.truth}/{sample}.vcf.gz',f'F1/split_{wildcards.query}/{sample}.vcf.gz',f'F1/{wildcards.truth}_{wildcards.query}.{sample}.summary.csv'])
    print(targets)
    return targets

rule gather_happy:
    input:
        aggregate_happy
    output:
        csv = 'F1/{truth}_{query}.csv'
    localrule: True
    shell:
        '''
        echo -e "variant region truth query recall precision truth_TiTv query_TiTv F1_score sample tissue" > {output}
        for i in {input}
        do
          awk -v I=$(basename $i) -F',' '$2=="*"&&($3~/CDS/||$3=="NCE"||$3=="intergenic")&&$4=="PASS" {{ split(I,a,"."); print $1,$3,$17,$38,$8,$9,$22,$43,$11,a[1],a[2] }}' $i >> {output}
        done
        '''
