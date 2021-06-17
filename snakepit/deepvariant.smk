from pathlib import PurePath
from itertools import product

class Default(dict):
    def __missing__(self, key):
        return '{'+key+'}'

def get_dir(base='work',ext='', **kwargs):
    if base == 'input':
        base_dir = config.get('bam_path','')
    elif base == 'output':
        base_dir = 'output'
    elif base == 'work':
        base_dir = get_dir('output','intermediate_results_{animal}',**kwargs)
    elif base == 'mendel':
        base_dir = 'mendel'
    elif base == 'main':
        base_dir = ''
    else:
        raise Exception('Base not found')
    return str(Path(base_dir.format_map(Default(kwargs))) / ext.format_map(Default(kwargs)))

wildcard_constraints:
     haplotype = r'asm|hap1|hap2|parent1|parent2',
     phase = r'unphased|phased',
     model = r'pbmm2|hybrid|bwa|mm2'

def get_model(wildcards,base='/opt/models',ext='model.ckpt'):
    model_location = f'{base}/{{}}/{ext}'
    if wildcards['model'] == 'pbmm2':
        return model_location.format('pacbio')
    elif wildcards['model'] == 'hybrid':
        return model_location.format('hybrid_pacbio_illumina')
    elif wildcards['model'] == 'bwa':
        return model_location.format('wgs')
    elif wildcards['model'] == 'mm2':
        return model_location.format('wgs')


##TODO determine how to pull singularity images
##TODO merge calls with snakemake format

def make_singularity_call(wildcards,extra_args='',tmp_bind='$TMPDIR',input_bind=True, output_bind=True, work_bind=True):
    call = f'singularity exec {extra_args} '
    if input_bind:
        call += f'-B {":/".join([get_dir("input",**wildcards)]*2)} '
    if output_bind:
        call += f'-B {":/".join([get_dir("output",**wildcards)]*2)} '
    if work_bind:
        call += f'-B {":/".join([get_dir("work",**wildcards)]*2)} '
    if tmp_bind:
        call += f'-B {tmp_bind}:/tmp '
    return call

for dir_ in ('input','output','work'):
    for animal,reference in product(config['animals'],config['reference']):
        Path(get_dir(dir_,animal=animal,ref=reference,**config)).mkdir(exist_ok=True)


def capture_logic():
    if trios in config:
        return read_trios()
    else:
        return [get_dir('main','cohort.autosomes.vcf.gz')]

def read_trios(ext='.vcf.gz'):

    import pandas as pd
    df = pd.read_csv(config['trios'])
    df.fillna('missing',inplace=True)

    targets = []
    for _, row in df.iterrows():
        targets.append(get_dir('mendel',f'{"_".join(row)}{ext}'))
    return targets

rule all:
    input:
        capture_logic()

rule deepvariant_make_examples:
    input:
        ref = multiext(config['reference'],'','.fai'),
        bam = multiext(get_dir('input','{animal}.bam'),'','.bai')
    output:
        example = temp(get_dir('work','make_examples.tfrecord-{N}-of-{sharding}.gz')),
        gvcf = temp(get_dir('work','gvcf.tfrecord-{N}-of-{sharding}.gz'))
    params:
        examples = lambda wildcards, output: PurePath(output[0]).with_name(f'make_examples.tfrecord@{config["shards"]}.gz'),
        gvcf = lambda wildcards, output: PurePath(output[1]).with_name(f'gvcf.tfrecord@{config["shards"]}.gz'),
        phase_args = lambda wildcards: '--parse_sam_aux_fields={i} --sort_by_haplotypes={i}'.format(i=('true' if False else 'false')),
        model_args = lambda wildcards: '--add_hp_channel --alt_aligned_pileup diff_channels --realign_reads=false --vsc_min_fraction_indels 0.12' if 'bwa' == 'pbmm2' else '',
        singularity_call = lambda wildcards,input: make_singularity_call(wildcards,extra_args=f'-B {PurePath(input.ref[0]).parent}:/reference/'),
        ref = lambda wildcards,input: f'/reference/{PurePath(input.ref[0]).name}',
        regions = ' '.join(map(str,range(1,30)))
    threads: 1
    resources:
        mem_mb = 10000,
        walltime = '14:00',
        disk_scratch = 1,
        use_singularity = True
    shell:
        '''
        {params.singularity_call} \
        {config[DV_container]} \
        /opt/deepvariant/bin/make_examples \
        --mode calling \
        --ref {params.ref} \
        --reads {input.bam[0]} \
        --examples {params.examples} \
        --gvcf {params.gvcf} \
        {params.model_args} \
        {params.phase_args} \
        --task {wildcards.N}
        #--regions {params.regions}
        '''

rule deepvariant_call_variants:
    input:
        (get_dir('work', f'make_examples.tfrecord-{N:05}-of-{config["shards"]:05}.gz') for N in range(config['shards']))
    output:
        temp(get_dir('work','call_variants_output.tfrecord.gz'))
    params:
        examples = lambda wildcards: (f'make_examples.tfrecord@{config["shards"]}.gz'),
        model = lambda wildcards: get_model({'model':'bwa'}),
        dir_ = lambda wildcards: get_dir('work',**wildcards),
        singularity_call = lambda wildcards, threads: make_singularity_call(wildcards,f'--env OMP_NUM_THREADS={threads}'),
        contain = lambda wildcards: config['DV_container'],
        vino = lambda wildcards: '--use_openvino'
    threads: 32
    resources:
        mem_mb = 4000,
        disk_scratch = 1,
        use_singularity = True,
        walltime = lambda wildcards: '8:00'
    shell:
        '''
        {params.singularity_call} \
        {params.contain} \
        /bin/bash -c "cd {params.dir_}; /opt/deepvariant/bin/call_variants \
        --outfile /{output} \
        --examples {params.examples} \
        --checkpoint {params.model} \
        {params.vino}"
        '''

rule deepvariant_postprocess:
    input:
        ref = multiext(config['reference'],'','.fai'),
        variants = get_dir('work','call_variants_output.tfrecord.gz'),
        gvcf = (get_dir('work', f'gvcf.tfrecord-{N:05}-of-{config["shards"]:05}.gz') for N in range(config['shards']))
    output:
        vcf = get_dir('output','{animal}.bwa.vcf.gz'),
        gvcf = get_dir('output','{animal}.bwa.g.vcf.gz')
    params:
        gvcf = lambda wildcards,input: PurePath(input.gvcf[0]).with_suffix('').with_suffix(f'.tfrecord@{config["shards"]}.gz'),
        singularity_call = lambda wildcards,input: make_singularity_call(wildcards,extra_args=f'-B {PurePath(input.ref[0]).parent}:/reference/'),
        ref = lambda wildcards,input: f'/reference/{PurePath(input.ref[0]).name}',
        contain = lambda wildcards: config['DV_container']
    threads: 1
    resources:
        mem_mb = 80000,
        walltime = '14:00',
        disk_scratch = 1,
        use_singularity = True
    shell:
        '''
        {params.singularity_call} \
        {params.contain} \
        /opt/deepvariant/bin/postprocess_variants \
        --ref {params.ref} \
        --infile {input.variants} \
        --outfile {output.vcf} \
        --gvcf_outfile {output.gvcf} \
        --nonvariant_site_tfrecord_path {params.gvcf} \
        --novcf_stats_report
        '''

rule split_gvcf_chromosomes:
    input:
        get_dir('output','{animal}.bwa.g.vcf.gz')
    output:
        gvcf = temp((get_dir('output','{animal}.bwa.{chr}.g.vcf.gz',chr=CHR) for CHR in range(1,30))),
        tbi = temp((get_dir('output','{animal}.bwa.{chr}.g.vcf.gz.tbi',chr=CHR) for CHR in range(1,30)))
    threads: 1
    resources:
        mem_mb = 3000,
        walltime = '25'
    run:
        for chromosome in range(1,30):
            out_file = output.gvcf[chromosome-1]
            shell(f'tabix -h {{input}} {chromosome} | bgzip -@ {{threads}} -c > {out_file}')
            shell(f'tabix -p vcf {out_file}')
    
rule GLnexus_merge_chrm:
    input:
        vcf = (get_dir('output','{animal}.bwa.{chr}.g.vcf.gz',animal=ANIMAL) for ANIMAL in config['animals']),
        tbi = (get_dir('output','{animal}.bwa.{chr}.g.vcf.gz.tbi',animal=ANIMAL) for ANIMAL in config['animals'])
    output:
        temp(multiext(get_dir('main','cohort.{chr}.vcf.gz'),'','.tbi'))
    params:
        gvcfs = lambda wildcards, input: list('/data/' / PurePath(fpath) for fpath in input.vcf),
        out = lambda wildcards, output: f'/data/{PurePath(output[0]).name}',
        DB = lambda wildcards, output: f'/tmp/GLnexus.DB',
        bed = lambda wildcards: '' if True else '--bed /data/BSW_autosome.bed',
        singularity_call = lambda wildcards: make_singularity_call(wildcards,'-B .:/data', input_bind=False, output_bind=False, work_bind=False),
        mem = lambda wildcards,threads,resources: threads*resources['mem_mb']/1000
    threads: 12 #force using 4 threads for bgziping
    resources:
        mem_mb = 8000,
        disk_scratch = 150,
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
        {params.bed} {params.gvcfs} \
        | bcftools view - | bgzip -@ 4 -c > {params.out}"
        tabix -p vcf {output[0]}
        '''

rule aggregate_merge_chrm:
    input:
        vcf = expand(get_dir('main','cohort.{chr}.vcf.gz'),chr=range(1,30)),
        tbi = expand(get_dir('main','cohort.{chr}.vcf.gz.tbi'),chr=range(1,30)),
    output:
        vcf = get_dir('main','cohort.autosomes.vcf.gz')
    threads: 8
    resources:
        mem_mb = 2000,
        walltime = '30'
    shell:
        '''
        bcftools concat --threads {threads} -o {output.vcf} -Oz {input.vcf}
        tabix -p vcf {output.vcf}
        #zgrep "#" {input[0]} > {output.temp}
        #zgrep -v "#" {input} | cat {output.temp} - | bgzip -@ {threads} -c > {output.vcf}
        #requires indexing, and so not worth it?
        #bcftools view -h {input[0]} > {output.temp}
        #bcftools view -H {input} | cat {output.temp} - | bgzip -@ {threads} -c > {output.vcf}
        '''

rule GLnexus_merge:
    input:
        expand(get_dir('output','{animal}.bwa.g.vcf.gz'),animal=config['animals'])
    output:
        get_dir('main','cohort.vcf.gz')
    params:
        gvcfs = lambda wildcards, input: list('/data/' / PurePath(fpath) for fpath in input),
        out = lambda wildcards, output: f'/data/{PurePath(output[0]).name}',
        DB = lambda wildcards, output: f'/tmp/GLnexus.DB',
        singularity_call = lambda wildcards: make_singularity_call(wildcards,'-B .:/data', input_bind=False, output_bind=False, work_bind=False),
        mem = lambda wildcards,threads,resources: threads*resources['mem_mb']/1500
    threads: 24 #force using 4 threads for bgziping
    resources:
        mem_mb = 7000,
        disk_scratch = 500,
        walltime = "4:00",
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
        | bcftools view - | bgzip -@ 4 -c > {params.out}"
        '''
        
rule GLnexus_merge_families:
    input:
        offspring = 'output/{offspring}.bwa.g.vcf.gz',
        sire  = lambda wildcards: 'output/{sire}.bwa.g.vcf.gz' if wildcards.sire != 'missing' else [],
        dam = lambda wildcards: 'output/{dam}.bwa.g.vcf.gz' if wildcards.dam != 'missing' else []
    output:
        get_dir('mendel','{offspring}_{sire}_{dam}.vcf.gz')
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

rule rtg_pedigree:
    output:
        get_dir('mendel','{offspring}_{sire}_{dam}.ped'
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
        vcf = get_dir('mendel','{offspring}_{sire}_{dam}.vcf.gz'),
        pedigree = get_dir('mendel','{offspring}_{sire}_{dam}.ped')
    output:
        temp = temp(get_dir('mendel','filled_{offspring}_{sire}_{dam}.vcf.gz')),
        results = multiext(get_dir('mendel','{offspring}_{sire}_{dam}'),'.inconsistent.vcf.gz','.inconsistent.stats','.mendel.log')
    params:
        vcf_in = lambda wildcards, input: '/data' / PurePath(input.vcf) if (wildcards.dam != 'missing' and wildcards.sire != 'missing') else '/data/mendel/filled_' + PurePath(input.vcf).name,
        vcf_annotated = lambda wildcards, output: '/data' / PurePath(output.results[0]),
        singularity_call = lambda wildcards: make_singularity_call(wildcards,'-B .:/data',input_bind=False,output_bind=False,work_bind=False)
    threads: 1
    resources:
        mem_mb = 10000,
        walltime = '30'
    shell:
        '''
        bcftools merge --no-index -o {output.temp} -Oz {input.vcf} {config[missing_template]}
        {params.singularity_call} \
        {config[RTG_container]} \
        /bin/bash -c "rtg mendelian -i {params.vcf_in} --output-inconsistent {params.vcf_annotated} --pedigree=/data/{input.pedigree} -t /data/{input.sdf} > /data/{output.results[2]}"
        bcftools stats {output.results[1]} | grep "^SN" > {output.results[1]}
        '''

rule mendel_summary:
    input:
        logs = read_trios('.mendel.log'),
        stats = read_trios('.inconsistent.stats')
    output:
        get_dir('main','mendel.summary.df')
    run:
        import pandas as pd

        rows = []
        for log_in, stat_in in zip(input.logs,input.stats):
            rows.append({k:v for k,v in zip(('offspring','sire','dam'),PurePath(log_in).with_suffix('').with_suffix('').name.split('_'))})
            with open(log_in,'r') as fin:
                for line in fin:
                    if 'violation of Mendelian constraints' in line:
                        violate, total = (int(i) for i in line.split()[0].split('/'))
                        rows[-1]['violate'] = violate
                        rows[-1]['total'] = total
            with open(stat_in,'r') as fin:
                for line in fin:
                    if 'number of SNPs' in line:
                        rows[-1]['SNP'] = int(line.rstrip().split()[-1])
                    elif 'number of indels' in line:
                        rows[-1]['indel'] = int(line.rstrip().split()[-1])
        
        df = pd.DataFrame(rows)
        df['rate']=df['violate']/df['total']
        df['duo']=(df == 'missing').any(axis=1)
        df.to_csv(output[0],index=False)

