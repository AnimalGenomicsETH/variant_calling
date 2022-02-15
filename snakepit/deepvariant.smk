from pathlib import Path,PurePath
from itertools import product

class Default(dict):
    def __missing__(self, key):
        return '{'+key+'}'

def get_dir(base='work',ext='', **kwargs):
    if base == 'input':
        base_dir = config.get('bam_path','')
    elif base == 'output':
        base_dir = 'output_DV'
    elif base == 'work':
        base_dir = get_dir('output','intermediate_results_{animal}',**kwargs)
    elif base == 'mendel':
        base_dir = '{caller}_mendel'
    elif base == 'main':
        base_dir = ''
    else:
        raise Exception('Base not found')
    return str(Path(base_dir.format_map(Default(kwargs))) / ext.format_map(Default(kwargs)))

if config.get("per_sample",True):
    ruleorder: split_gvcf_chromosomes > deepvariant_postprocess 
else:
    ruleorder: deepvariant_postprocess > split_gvcf_chromosomes


include: 'mendelian.smk'

wildcard_constraints:
    haplotype = r'asm|hap1|hap2|parent1|parent2',
    phase = r'unphased|phased',
    model = r'pbmm2|hybrid|bwa|mm2',
    caller = r'DV|GATK',
    chr = r'\d*|all|autosomes'
    
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

for dir_ in ('output','work'):
    for animal,reference in product(config['animals'],config['reference']):
        Path(get_dir(dir_,animal=animal,ref=reference,**config)).mkdir(exist_ok=True)


def capture_logic():
    if 'trios' in config:
        return [get_dir('main','mendel.{caller}.autosomes.summary.df',caller=caller) for caller in ('DV','GATK')]
    elif config.get('regions') == 'autosomes':
        return [get_dir('main','cohort.autosomes.vcf.gz')]
    else:
        return [get_dir('main','cohort.all.vcf.gz')]

rule all:
    input:
        capture_logic()

PIG = list(range(1,19))
def get_regions(wildcards):
    if 'regions' not in config:
        return ''
    elif config.get('regions',False) == 'autosomes':
        return f'--regions "{" ".join(map(str,PIG))}"'
    else:
        return f'--regions {wildcards.chromosome}'

rule deepvariant_make_examples:
    input:
        ref = multiext(config['reference'],'','.fai'),
        bam = multiext(get_dir('input','{animal}.bam'),'','.bai')
    output:
        example = temp(get_dir('work','make_examples.{chromosome}.tfrecord-{N}-of-{sharding}.gz')),
        gvcf = temp(get_dir('work','gvcf.{chromosome}.tfrecord-{N}-of-{sharding}.gz'))
    params:
        examples = lambda wildcards, output: PurePath(output['example']).with_suffix('').with_suffix(f'.tfrecord@{config["shards"]}.gz'),
        gvcf = lambda wildcards, output: PurePath(output['gvcf']).with_suffix('').with_suffix(f'.tfrecord@{config["shards"]}.gz'),
        phase_args = lambda wildcards: '--parse_sam_aux_fields={i} --sort_by_haplotypes={i}'.format(i=('true' if False else 'false')),
        model_args = lambda wildcards: '--add_hp_channel --alt_aligned_pileup diff_channels --realign_reads=false --vsc_min_fraction_indels 0.12' if 'bwa' == 'pbmm2' else '',
        singularity_call = lambda wildcards,input: make_singularity_call(wildcards,input_bind=False,extra_args=f'-B {PurePath(input.ref[0]).parent}:/reference/ -B {Path(input.bam[0]).resolve().parent}:{PurePath(input.bam[0]).parent}'),
        ref = lambda wildcards,input: f'/reference/{PurePath(input.ref[0]).name}',
        regions = lambda wildcards: get_regions(wildcards) #f'--regions "{" ".join(map(str,range(1,30)))}"' if config.get('autosomes',False) else ''
    threads: 1
    resources:
        mem_mb = 6000,
        walltime = '4:00',
        disk_scratch = 1,
        use_singularity = True
    shell:
        #--include_med_dp \
        '''
        {params.singularity_call} \
        {config[DV_container]} \
        /opt/deepvariant/bin/make_examples \
        --mode calling \
        --ref {params.ref} \
        --reads {input.bam[0]} \
        --examples {params.examples} \
        --gvcf {params.gvcf} \
        {params.regions} \
        {params.model_args} \
        {params.phase_args} \
        --task {wildcards.N}
        '''

rule deepvariant_call_variants:
    input:
        (get_dir('work', f'make_examples.{{chromosome}}.tfrecord-{N:05}-of-{config["shards"]:05}.gz') for N in range(config['shards']))
    output:
        temp(get_dir('work','call_variants_output.{chromosome}.tfrecord.gz'))
    params:
        examples = lambda wildcards,input: PurePath(input[0]).with_suffix('').with_suffix(f'.tfrecord@{config["shards"]}.gz'),
        model = lambda wildcards: get_model({'model':'bwa'}),
        dir_ = lambda wildcards: get_dir('work',**wildcards),
        singularity_call = lambda wildcards, threads: make_singularity_call(wildcards,f'--env OMP_NUM_THREADS={threads}'),
        contain = lambda wildcards: config['DV_container'],
        vino = lambda wildcards: '--use_openvino'
    threads: 24
    resources:
        mem_mb = 2500,
        disk_scratch = 1,
        use_singularity = True,
        walltime = lambda wildcards: '4:00',
        use_AVX512 = True
    shell:
        '''
        {params.singularity_call} \
        {params.contain} \
        /bin/bash -c "cd /tmp; /opt/deepvariant/bin/call_variants \
        --outfile /{output} \
        --examples /{params.examples} \
        --checkpoint {params.model} \
        {params.vino}"
        '''

rule deepvariant_postprocess:
    input:
        ref = multiext(config['reference'],'','.fai'),
        variants = get_dir('work','call_variants_output.{chromosome}.tfrecord.gz'),
        gvcf = (get_dir('work', f'gvcf.{{chromosome}}.tfrecord-{N:05}-of-{config["shards"]:05}.gz') for N in range(config['shards']))
    output:
        vcf = get_dir('output','{animal}.bwa.{chromosome}.vcf.gz'),
        gvcf = get_dir('output','{animal}.bwa.{chromosome}.g.vcf.gz')
    params:
        gvcf = lambda wildcards,input: PurePath(input.gvcf[0]).with_suffix('').with_suffix(f'.tfrecord@{config["shards"]}.gz'),
        singularity_call = lambda wildcards,input: make_singularity_call(wildcards,extra_args=f'-B {PurePath(input.ref[0]).parent}:/reference/'),
        ref = lambda wildcards,input: f'/reference/{PurePath(input.ref[0]).name}',
        contain = lambda wildcards: config['DV_container']
    threads: 1
    resources:
        mem_mb = 50000,
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
        get_dir('output','{animal}.bwa.{chromosome}.g.vcf.gz',chromosome = config.get('regions','all'))
    output:
        gvcf = temp((get_dir('output','{animal}.bwa.{chr}.g.vcf.gz',chr=CHR) for CHR in PIG)),
        tbi = temp((get_dir('output','{animal}.bwa.{chr}.g.vcf.gz.tbi',chr=CHR) for CHR in PIG))
    threads: 1
    resources:
        mem_mb = 3000,
        walltime = '25'
    run:
        for chromosome in PIG:
            out_file = output.gvcf[chromosome-1]
            shell(f'tabix -h {{input}} {chromosome} | bgzip -@ {{threads}} -c > {out_file}')
            shell(f'tabix -p vcf {out_file}')
    
rule GLnexus_merge_chrm:
    input:
        vcf = (get_dir('output','{animal}.bwa.{chr}.g.vcf.gz',animal=ANIMAL) for ANIMAL in config['animals']),
        #tbi = (get_dir('output','{animal}.bwa.{chr}.g.vcf.gz.tbi',animal=ANIMAL) for ANIMAL in config['animals'])
    output:
        multiext(get_dir('main','cohort.{chr,\d+}.vcf.gz'),'','.tbi')
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

rule aggregate_autosomes:
    input:
        vcf = expand(get_dir('main','cohort.{chr}.vcf.gz'),chr=PIG),
        #tbi = expand(get_dir('main','cohort.{chr}.vcf.gz.tbi'),chr=PIG),
    output:
        multiext(get_dir('main','cohort.autosomes.vcf.gz'),'','.tbi')
    threads: 8
    resources:
        mem_mb = 2000,
        walltime = '60'
    shell:
        '''
        bcftools concat --threads {threads} -o {output[0]} -Oz {input.vcf}
        tabix -p vcf {output[0]}
        '''

rule GLnexus_merge:
    input:
        expand(get_dir('output','{animal}.bwa.all.g.vcf.gz'),animal=config['animals'])
    output:
        multiext(get_dir('main','cohort.all.vcf.gz'),'','.tbi')
    params:
        gvcfs = lambda wildcards, input: list('/data/' / PurePath(fpath) for fpath in input),
        out = lambda wildcards, output: f'/data/{PurePath(output[0]).name}',
        DB = lambda wildcards, output: f'/tmp/GLnexus.DB',
        singularity_call = lambda wildcards: make_singularity_call(wildcards,'-B .:/data', input_bind=False, output_bind=False, work_bind=False),
        mem = lambda wildcards,threads,resources: threads*resources['mem_mb']/1024
    threads: 32 #force using threads for bgziping
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
        | bcftools view - | bgzip -@ {threads} -c > {params.out}"
        tabix -p vcf {output[0]}
        '''
        
