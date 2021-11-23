def read_trios(ext='.{chr}.vcf.gz'):
    if 'trios' not in config:
        return []
    import pandas as pd
    df = pd.read_csv(config['trios'])
    df.fillna('missing',inplace=True)

    targets = []
    for _, row in df.iterrows():
        targets.append(get_dir('mendel',f'{"_".join(row)}{ext}'))
    return targets

ruleorder: remove_tigs > GLnexus_merge_families

rule GLnexus_merge_families:
    input:
        offspring = get_dir('output','{offspring}.bwa.{chr}.g.vcf.gz'),
        sire  = lambda wildcards: get_dir('output','{sire}.bwa.{chr}.g.vcf.gz') if wildcards.sire != 'missing' else [],
        dam = lambda wildcards: get_dir('output','{dam}.bwa.{chr}.g.vcf.gz') if wildcards.dam != 'missing' else []
    output:
        get_dir('mendel','{offspring}_{sire}_{dam}.{chr}.vcf.gz')
    params:
        gvcfs = lambda wildcards, input: list('/data/' / PurePath(fpath) for fpath in input),
        out = lambda wildcards, output: '/data' / PurePath(output[0]),
        DB = lambda wildcards, output: f'/tmp/GLnexus.DB',
        singularity_call = lambda wildcards: make_singularity_call(wildcards,'-B .:/data', input_bind=False, output_bind=False, work_bind=False),
        mem = lambda wildcards,threads,resources: threads*resources['mem_mb']/1500
    threads: 8
    resources:
        mem_mb = 6000,
        disk_scratch = 50,
        walltime = '4:00',
        use_singularity = True
    shell:
        '''
        ulimit -Sn 4096
        {params.singularity_call} \
        {config[GL_container]} \
        /bin/bash -c " /usr/local/bin/glnexus_cli \
        --dir {params.DB} \
        --config DeepVariantWGS \
        --threads {threads} \
        --mem-gbytes {params.mem} \
        {params.gvcfs} \
        | bcftools view - | bgzip -@ 2 -c > {params.out}"
        '''

rule remove_tigs:
    input:
        multiext(get_dir('mendel','{offspring}_{sire}_{dam}.all.vcf.gz'),'','.tbi')
    output:
        get_dir('mendel','{offspring}_{sire}_{dam}.autosomes.vcf.gz')
    params:
        regions = ','.join(map(str,range(1,30)))
    threads: 2
    resources:
        mem_mb = 2000,
        walltime = '20'
    shell:
        'bcftools view --threads {threads} -r {params.regions} -O z -o {output} {input[0]}'

rule tabix:
    input:
        '{vcf}'
    output:
        '{vcf}.tbi'
    resources:
        mem_mb = 2000,
        walltime = '20'
    shell:
        'tabix -p vcf {input}'

rule rtg_pedigree:
    output:
        get_dir('mendel','{offspring}_{sire}_{dam}.ped')
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

rule rtg_format:
    input:
        ref = lambda wildcards: multiext(config['reference'],'','.fai')
    output:
        sdf = get_dir('main','ARS.sdf')
    params:
        singularity_call = lambda wildcards,input: make_singularity_call(wildcards,extra_args=f'-B {PurePath(input.ref[0]).parent}:/reference/,.:/data',tmp_bind=False,output_bind=False,work_bind=False),
        ref = lambda wildcards,input: f'/reference/{PurePath(input.ref[0]).name}',
    shell:
        '''
        {params.singularity_call} \
        {config[RTG_container]} \
        rtg format -o /data/{output.sdf} {params.ref}
        '''

rule rtg_mendelian_concordance:
    input:
        sdf = get_dir('main','ARS.sdf'),
        vcf = get_dir('mendel','{offspring}_{sire}_{dam}.{chr}.vcf.gz'),
        pedigree = get_dir('mendel','{offspring}_{sire}_{dam}.ped')
    output:
        temp = temp(get_dir('mendel','filled_{offspring}_{sire}_{dam}.{chr}.vcf.gz')),
        results = multiext(get_dir('mendel','{offspring}_{sire}_{dam}.{chr}'),'.inconsistent.vcf.gz','.inconsistent.stats','.mendel.log')
    params:
        vcf_in = lambda wildcards, input, output: '/data' / PurePath(input.vcf) if (wildcards.dam != 'missing' and wildcards.sire != 'missing') else '/data' / PurePath(output.temp),
        vcf_annotated = lambda wildcards, output: '/data' / PurePath(output.results[0]),
        singularity_call = lambda wildcards: make_singularity_call(wildcards,'-B .:/data',input_bind=False,output_bind=False,work_bind=False)
    threads: 1
    resources:
        mem_mb = 10000,
        walltime = '30'
    shell:
        '''
        #bcftools view -s BSWCHEM120096296001 -Ou RM808_BSWCHEM120096296001_missing.vcf.gz | bcftools reheader --samples <(echo missing) | bcftools merge --no-index -o filled_RM808_BSWCHEM120096296001_missing.vcf.gz -Oz RM808_BSWCHEM120096296001_missing.vcf.gz -
        bcftools merge --no-index -o {output.temp} -Oz {input.vcf} {config[missing_template]}
        {params.singularity_call} \
        {config[RTG_container]} \
        /bin/bash -c "rtg mendelian -i {params.vcf_in} --output-inconsistent {params.vcf_annotated} --pedigree=/data/{input.pedigree} -t /data/{input.sdf} > /data/{output.results[2]}"
        bcftools stats {output.results[0]} | grep "^SN" > {output.results[1]}
        '''

rule rtg_vcfeval:
    input:
        sdf = get_dir('main','ARS.sdf'),
        vcf = get_dir('mendel','{offspring}_{sire}_{dam}.{chr}.vcf.gz')
    output:
        logs = get_dir('mendel','{offspring}_{sire}_{dam}.{chr}.vcfeval.log')
    params:
        samples = lambda wildcards: f'{wildcards.sire if wildcards.sire != "missing" else wildcards.dam},{wildcards.offspring}',
        singularity_call = lambda wildcards: make_singularity_call(wildcards,'-B .:/data',input_bind=False,output_bind=False,work_bind=False)
    threads: 4
    resources:
        mem_mb = 5000,
        disk_scratch = 20
    shell:
        '''
        {params.singularity_call} \
        {config[RTG_container]} \
        /bin/bash -c "cd $TMPDIR; rtg vcfeval -t /data/{input.sdf} -b /data/{input.vcf} -c /data/{input.vcf} --sample {params.samples} -o null_output -T {threads} --no-roc --squash-ploidy --output-mode=roc-only > /data/{output}"
        '''

rule bcftools_mendelian:
    input:
        vcf = get_dir('mendel','{offspring}_{sire}_{dam}.{chr}.vcf.gz')
    output:
        logs = get_dir('mendel','{offspring}_{sire}_{dam}.{chr}.mendelian.log')
    params:
        sample = '{dam},{sire},{offspring}'
    threads: 1
    resources:
        mem_mb = 3000,
        walltime = '30'
    shell:
        '''
        bcftools +mendelian {input.vcf} -t {params.sample} -m c -m x > >(bcftools stats - | grep "SN" >> {output}) 2> >(grep -v "#" >> {output})
        '''

rule bcftools_count:
    input:
        get_dir('mendel','{offspring}_{sire}_{dam}.{chr}.vcf.gz')
    output:
        get_dir('mendel','{offspring}_{sire}_{dam}.{chr}.count')
    shell:
        '''
        tabix -fp vcf {input}
        bcftools index -n {input} > {output}
        '''

rule genotype_count:
    input:
        get_dir('mendel','{offspring}_{sire}_{dam}.{chr}.vcf.gz')
    output:
        get_dir('mendel','{offspring}_{sire}_{dam}.{chr}.genotypes')
    shell:
        '''
        bcftools annotate -x INFO,^FORMAT/GT {input} | grep -oP "([\.|\d]/[/.|\d])" | sort | uniq -c > {output}
        '''

rule mendel_summary:
    input:
        logs = read_trios('.{chr}.mendelian.log'),
        stats = read_trios('.{chr}.count')
    output:
        get_dir('main','mendel.{caller}.{chr}.summary.df')
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

