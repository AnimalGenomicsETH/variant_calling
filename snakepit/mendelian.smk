def read_trios(ext='.{chr}.vcf.gz'):
    if 'trios' not in config:
        return []
    import pandas as pd
    df = pd.read_csv(config['trios'])
    df.fillna('missing',inplace=True)

def get_pedigree_samples(wildcards):
    return ''

#TODO: fix with new code
rule GLnexus_merge_pedigree:
    input:
        expand('{run}/deepvariant/{sample}.{region}.g.vcf.gz',sample=get_pedigree_samples)
    output:
        vcf = '{run}/mendelian/{pedigree}.{region}.vcf.gz'
    params:
        gvcfs = lambda wildcards, input: list('/data/' / PurePath(fpath) for fpath in input),
        out = lambda wildcards, output: '/data' / PurePath(output[0]),
        singularity_call = lambda wildcards: make_singularity_call(wildcards,'-B .:/data', input_bind=False, output_bind=False, work_bind=False),
        mem = lambda wildcards,threads,resources: threads*resources['mem_mb_per_cpu']/1500
    threads: 8
    resources:
        mem_mb_per_cpu = 6000,
        runtime = '4:00'
    shell:
        '''
        ulimit -Sn 4096
        {params.singularity_call} \
        {config[GL_container]} \
        /bin/bash -c " /usr/local/bin/glnexus_cli \
        --dir $TMPDIR \
        --config DeepVariantWGS \
        --threads {threads} \
        --mem-gbytes {params.mem} \
        {params.gvcfs} \
        | bcftools view - | bgzip -@ 2 -c > {params.out}"
        '''

#TODO: build better pedigrees if needed
rule rtg_pedigree:
    output:
        'mendelian/{offspring}_{sire}_{dam}.ped'
    shell:
        '''
        FILE={output}
cat <<EOM >$FILE
#PED format pedigree
#
#fam-id/ind-id/pat-id/mat-id: 0=unknown
#sex: 1=male; 2=female; 0=unknown
#phenotype: -9=missing, 0=missing; 1=unaffected; 2=affected
#
#fam-id ind-id pat-id mat-id sex phen
1 {wildcards.offspring} {wildcards.sire} {wildcards.dam} 1 0
1 {wildcards.sire} 0 0 1 0
1 {wildcards.dam} 0 0 2 0
EOM
        '''

#TODO: fix singularity if needed
rule rtg_format:
    input:
        reference = lambda wildcards: multiext(config['reference'],'','.fai')
    output:
        sdf = 'mendelian/reference.sdf'
    container: ''
    shell:
        '''
        rtg format -o {output.sdf} {input.reference[0]}
        '''

rule rtg_mendelian_concordance:
    input:
        sdf = rules.rtg_format.output['sdf'],
        vcf = rules.GLnexus_merge_pedigree.output['vcf'],
        pedigree = 'mendelian/{offspring}_{sire}_{dam}.ped'
    output:
        temp = temp('mendelian/filled_{offspring}_{sire}_{dam}.{chr}.vcf.gz'),
        results = multiext('mendelian/{offspring}_{sire}_{dam}.{chr}','.inconsistent.vcf.gz','.inconsistent.stats','.mendel.log')
    params:
        vcf_in = lambda wildcards, input, output: '/data' / PurePath(input.vcf) if (wildcards.dam != 'missing' and wildcards.sire != 'missing') else '/data' / PurePath(output.temp),
        vcf_annotated = lambda wildcards, output: '/data' / PurePath(output.results[0]),
    threads: 1
    resources:
        mem_mb_per_cpu = 10000,
        runtime = '30'
    container: ''
    shell:
        '''
        bcftools merge --no-index -o {output.temp} -Oz {input.vcf} {config[missing_template]}
        
        rtg mendelian -i {params.vcf_in} --output-inconsistent {params.vcf_annotated} --pedigree=/data/{input.pedigree} -t /data/{input.sdf} > /data/{output.results[2]}
        
        bcftools stats {output.results[0]} | grep "^SN" > {output.results[1]}
        '''

rule rtg_vcfeval:
    input:
        sdf = rules.rtg_format.output['sdf'],
        vcf = 'mendelian/{offspring}_{sire}_{dam}.{chr}.vcf.gz'
    output:
        logs = 'mendelian/{offspring}_{sire}_{dam}.{chr}.vcfeval.log'
    params:
        samples = lambda wildcards: f'{wildcards.sire if wildcards.sire != "missing" else wildcards.dam},{wildcards.offspring}',
    threads: 4
    resources:
        mem_mb_per_cpu = 5000,
    container: ''
    shell:
        '''
cd $TMPDIR

rtg vcfeval -t /data/{input.sdf} -b /data/{input.vcf} -c /data/{input.vcf} \
--sample {params.samples} -o null_output -T {threads} --no-roc \
--squash-ploidy --output-mode=roc-only > /data/{output}
        '''

rule bcftools_mendelian:
    input:
        vcf = 'mendelian/{offspring}_{sire}_{dam}.{chr}.vcf.gz'
    output:
        logs = 'mendelian/{offspring}_{sire}_{dam}.{chr}.mendelian.log'
    params:
        sample = '{dam},{sire},{offspring}'
    threads: 1
    resources:
        mem_mb_per_cpu = 3000,
        runtime = '30'
    shell:
        '''
bcftools +mendelian {input.vcf} -t {params.sample} -m c -m x > >(bcftools stats - | grep "SN" >> {output}) 2> >(grep -v "#" >> {output})
        '''

rule bcftools_count:
    input:
        'mendelian/{offspring}_{sire}_{dam}.{chr}.vcf.gz'
    output:
        'mendelian/{offspring}_{sire}_{dam}.{chr}.count'
    shell:
        '''
        tabix -fp vcf {input}
        bcftools index -n {input} > {output}
        '''

rule genotype_count:
    input:
        'mendelian/{offspring}_{sire}_{dam}.{chr}.vcf.gz'
    output:
        'mendelian/{offspring}_{sire}_{dam}.{chr}.genotypes'
    shell:
        '''
bcftools annotate -x INFO,^FORMAT/GT {input} | grep -oP "([\\.|\\d]/[/.|\\d])" | sort | uniq -c > {output}
        '''

rule mendel_summary:
    input:
        logs = read_trios('.{chr}.mendelian.log'),
        stats = read_trios('.{chr}.count')
    output:
        'mendelian/{caller}.{chr}.summary.df'
    run:
        import pandas as pd

        rows = []
        #for log_in, stat_in in zip(input.logs,input.stats):
        for log_in in input.logs:
            rows.append({k:v for k,v in zip(('offspring','sire','dam'),PurePath(log_in).with_suffix('').with_suffix('').with_suffix('').name.split('_'))})
            with open(log_in,'r') as fin:
                for i,line in enumerate(fin):
                    parts = line.rstrip().split()
                    if parts[0].isdigit():
                        rows[-1]['consistent'] = int(parts[0])
                        rows[-1]['inconsistent'] = int(parts[1])
                        rows[-1]['uninformative'] = int(parts[2])
                    elif line[0] != '#' and 'number of SNPs:' in line:
                        rows[-1]['SNP'] = int(parts[-1])
                    elif line[0] != '#' and 'number of indels:' in line:
                        rows[-1]['indel'] = int(parts[-1])
            with open(PurePath(log_in).with_suffix('').with_suffix('.count'),'r') as fin:
                rows[-1]['vcf_count'] = int(fin.readline())
        df = pd.DataFrame(rows)
        df['rate'] = df['inconsistent']/(df['consistent']+df['inconsistent'])
        df['duo'] = (df == 'missing').any(axis=1)
        df['caller'] = wildcards.caller
        df.to_csv(output[0],index=False)