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

rule all:
    input:
        get_dir('main','cohort.vcf.gz')

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
        mem_mb = 6000,
        walltime = '4:00',
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
    threads: 36
    resources:
        mem_mb = 4000,
        disk_scratch = 1,
        use_singularity = True,
        walltime = lambda wildcards: '4:00'
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
        mem_mb = 40000,
        walltime = '4:00',
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
        (get_dir('output','{animal}.bwa.{chr}.g.vcf.gz',animal=ANIMAL) for ANIMAL in config['animals'])
    output:
        temp(get_dir('main','cohort.{chr}.vcf'))
    params:
        gvcfs = lambda wildcards, input: list('/data/' / PurePath(fpath) for fpath in input),
        out = lambda wildcards, output: f'/data/{PurePath(output[0]).name}',
        DB = lambda wildcards, output: f'/tmp/GLnexus.DB',
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
        #--bed /data/BSW_autosome.bed \
        --threads {threads} \
        --mem-gbytes {params.mem} \
        {params.gvcfs} \
        | bcftools view - | bgzip -@ 4 -c > {params.out}"
        '''

rule aggregate_merge_chrm:
    input:
        expand(get_dir('main','cohort.{chr}.vcf.gz'),chr=range(1,30))
    output:
        temp = get_dir('main','cohort.temp.vcf'),
        vcf = get_dir('main','cohort.autosomes.vcf.gz')
    threads: 8
    resources:
        mem_mb = 2000,
        walltime = '30'
    shell:
        '''
        zgrep "#" {input[0]} > {output.temp}
        zgrep -v "#" {input} | cat {output.temp} - | bgzip -@ {threads} -c > {output.vcf}
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
        mem = lambda wildcards,threads,resources: threads*resources['mem_mb']/1000
    threads: 24 #force using 4 threads for bgziping
    resources:
        mem_mb = 4000,
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
        
