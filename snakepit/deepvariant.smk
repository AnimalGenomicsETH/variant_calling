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
    elif base == 'splits':
        base_dir = 'split_files'
    elif base == 'work':
        base_dir = get_dir('output','intermediate_results_{animal}',**kwargs)
    elif base == 'mendel':
        base_dir = '{caller}_mendel'
    elif base == 'main':
        base_dir = ''
    else:
        raise Exception('Base not found')
    return str(Path(base_dir.format_map(Default(kwargs))) / ext.format_map(Default(kwargs)))

def get_walltime(wildcards, attempt):
    return {1:'4:00',2:'24:00'}[attempt]

if 'regions_bed' in config:
    ruleorder: bcftools_scatter > GLnexus_merge

#if config.get("merging","all") == "all":
#    ruleorder: GLnexus_merge > aggregate_chromosomes
#else:
#    ruleorder: aggregate_chromosomes > GLnexus_merge

#include: 'mendelian.smk'
def determine_run_regions():
    regions = []
    if 'regions_bed' in config:
        with open(config['regions_bed'],'r') as bed:
            for line in bed:
                parts = line.rstrip().split()
                if len(parts) == 1:
                    regions.append(parts[0])
                elif len(parts) == 2:
                    regions.append('{}:0-{}'.format(*parts))
                else:
                    regions.append('{}:{}-{}'.format(*parts[:3]))
    else:
        with open(config['reference']+'.fai','r') as fai:
            regions.append([line.split()[0] for line in fin])
    return regions

CHROMOSOMES=['A','B']

wildcard_constraints:
    haplotype = r'asm|hap1|hap2|parent1|parent2',
    phase = r'unphased|phased',
    model = r'pbmm2|hybrid|bwa|mm2',
    caller = r'DV|GATK',
    #chrom = r'all|chromosomes' + r'|'.join(CHROMOSOMES),
    cohort = r'\w+'
    
def get_model(model,base='/opt/models',ext='model.ckpt'):
    model_location = f'{base}/{{}}/{ext}'
    if model == 'pbmm2':
        return model_location.format('pacbio')
    elif model == 'hybrid':
        return model_location.format('hybrid_pacbio_illumina')
    elif model == 'bwa':
        return model_location.format('wgs')
    elif model == 'mm2':
        return model_location.format('pacbio')


def make_singularity_call(wildcards,extra_args='',tmp_bind='$TMPDIR',input_bind=False, output_bind=False, work_bind=False):
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
        return [get_dir('main','mendel.{caller}.chromosomes.summary.df',caller=caller) for caller in ('DV','GATK')]
    #elif config.get('regions') == 'chromosomes':
    #    return [get_dir('main',f'{config.get("cohort","cohort")}.chromosomes.{config.get("GL_config","WGS")}.vcf.gz')]
    else:
        return [get_dir('main',f'{config.get("cohort","cohort")}.{config.get("regions","all")}.{str(PurePath(config.get("GL_config","WGS")).with_suffix(""))}.vcf.gz')]

REGIONS = determine_run_regions()
rule all:
    input:
        expand(get_dir('main','{cohort}.{regions}.{preset}.vcf.gz{ext}'),ext=('','.tbi'),cohort=config.get("cohort","cohort"),regions=REGIONS,preset=config.get("GL_config","WGS"))
        #"regions" if 'regions_bed' in config else config.get('chromosome',0),preset=config.get("GL_config","WGS"))
        #capture_logic()

def get_regions(wildcards):

    if wildcards.callset == 'regions':
        return f'--regions "{" ".join(map(str,REGIONS))}"'
    else:
        return f'--regions {wildcards.callset}'
    if config.get('regions','all') == 'all':
        return ''
    elif config.get('regions',False) == 'chromosomes':
        #TEMP CODE
        return f'--exclude_regions "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18"'
        #return f'--regions "{" ".join(map(str,CHROMOSOMES))}"'
    else:
        return f'--regions {wildcards.chromosome}'
#how to handle 0 region
#--exclude_regions

## DEEPVARIANT SPECIAL ARG HANDLING
def _is_quoted(value):
  if value.startswith('"') and value.endswith('"'):
    return True
  if value.startswith("'") and value.endswith("'"):
    return True
  return False

def _add_quotes(value):
  if isinstance(value, str) and _is_quoted(value):
    return value
  return '"{}"'.format(value)

def _extra_args_to_dict(extra_args):
  """Parses comma-separated list of flag_name=flag_value to dict."""
  args_dict = {}
  if extra_args is None:
    return args_dict
  for extra_arg in extra_args.split(','):
    (flag_name, flag_value) = extra_arg.split('=')
    flag_name = flag_name.strip('-')
    # Check for boolean values.
    if flag_value.lower() == 'true':
      flag_value = True
    elif flag_value.lower() == 'false':
      flag_value = False
    args_dict[flag_name] = flag_value
  return args_dict

def _extend_command_by_args_dict(command_S, extra_args):
  """Adds `extra_args` to the command string."""
  command = []
  for key in sorted(extra_args):
    value = extra_args[key]
    if value is None:
      continue
    if isinstance(value, bool):
      added_arg = '' if value else 'no'
      added_arg += key
      command.extend(['--' + added_arg])
    else:
      command.extend(['--' + key, _add_quotes(value)])
  return command_S + ' '.join(command)

def get_model_params(model):
    special_args = {}
    if model == 'bwa':
        special_args['channels'] = 'insert_size'
    elif model == 'mm2':
        special_args['add_hp_channel'] = True
        special_args['alt_aligned_pileup'] = 'diff_channels'
        special_args['max_reads_per_partition'] = 600
        special_args['min_mapping_quality'] = 1
        special_args['parse_sam_aux_fields'] = True
        special_args['partition_size'] = 25000
        special_args['phase_reads'] = True
        special_args['pileup_image_width'] = 199
        special_args['realign_reads'] = False
        special_args['sort_by_haplotypes'] = True
        special_args['track_ref_reads'] = True
        special_args['vsc_min_fraction_indels'] = 0.12

    return _extend_command_by_args_dict('',special_args)


## END DEEPVARIANT SPECIAL ARG HANDLING
#error can be caused by corrupted file, check gzip -t -v
rule deepvariant_make_examples:
    input:
        ref = multiext(config['reference'],'','.fai'),
        bam = multiext(get_dir('input','{animal}.bam'),'','.bai')
    output:
        example = temp(get_dir('work','make_examples.{callset}.tfrecord-{N}-of-{sharding}.gz')),
        gvcf = temp(get_dir('work','gvcf.{callset}.tfrecord-{N}-of-{sharding}.gz'))
    params:
        examples = lambda wildcards, output: PurePath(output['example']).with_suffix('').with_suffix(f'.tfrecord@{config["shards"]}.gz'),
        gvcf = lambda wildcards, output: PurePath(output['gvcf']).with_suffix('').with_suffix(f'.tfrecord@{config["shards"]}.gz'),
        model_args = get_model_params(config['model']),
        singularity_call = lambda wildcards, input, output: make_singularity_call(wildcards,input_bind=False,output_bind=True,work_bind=False,extra_args=f'-B {PurePath(output["example"]).parent}:/{PurePath(output["example"]).parent} -B {PurePath(input.ref[0]).parent}:/reference/ -B {Path(input.bam[0]).resolve().parent}:{PurePath(input.bam[0]).parent}'),
        ref = lambda wildcards,input: f'/reference/{PurePath(input.ref[0]).name}',
        contain = lambda wildcards: config['DV_container'],
        regions = lambda wildcards: get_regions(wildcards) 
    threads: 1
    resources:
        mem_mb = config.get('resources',{}).get('make_examples',{}).get('mem_mb',6000),
        walltime = config.get('resources',{}).get('make_examples',{}).get('walltime','4:00'),#get_walltime,
        disk_scratch = 1,
        use_singularity = True
    priority: 50
    shell:
        '''
        {params.singularity_call} \
        {params.contain} \
        /opt/deepvariant/bin/make_examples \
        --mode calling \
        --ref {params.ref} \
        --include_med_dp \
        --reads {input.bam[0]} \
        --examples {params.examples} \
        --gvcf {params.gvcf} \
        --sample_name {wildcards.animal} \
        {params.regions} \
        {params.model_args} \
        --task {wildcards.N}
        '''

rule deepvariant_call_variants:
    input:
        (get_dir('work', f'make_examples.{{callset}}.tfrecord-{N:05}-of-{config["shards"]:05}.gz') for N in range(config['shards']))
    output:
        temp(get_dir('work','call_variants_output.{callset}.tfrecord.gz'))
    params:
        examples = lambda wildcards,input: PurePath(input[0]).with_suffix('').with_suffix(f'.tfrecord@{config["shards"]}.gz'),
        model = get_model(config['model']),
        dir_ = lambda wildcards: get_dir('work',**wildcards),
        singularity_call = lambda wildcards, output, threads: make_singularity_call(wildcards,f'--env OMP_NUM_THREADS={threads} -B {PurePath(output[0]).parent}:/{PurePath(output[0]).parent}'),
        contain = lambda wildcards: config['DV_container'],
        vino = '--use_openvino --openvino_model_dir /tmp' if config.get('openvino',False) else ''
    threads: config.get('resources',{}).get('call_variants',{}).get('threads',12)
    resources:
        mem_mb = config.get('resources',{}).get('call_variants',{}).get('mem_mb',1500),
        disk_scratch = 1,
        use_singularity = True,
        walltime = config.get('resources',{}).get('call_variants',{}).get('walltime','2G4:00'),
        #use_AVX512 = True
    priority: 90
    shell:
        '''
        {params.singularity_call} \
        {params.contain} \
        /opt/deepvariant/bin/call_variants \
        --outfile /{output} \
        --examples /{params.examples} \
        --checkpoint {params.model} \
        {params.vino}
        '''

rule deepvariant_postprocess:
    input:
        ref = multiext(config['reference'],'','.fai'),
        variants = get_dir('work','call_variants_output.{callset}.tfrecord.gz'),
        gvcf = (get_dir('work', f'gvcf.{{callset}}.tfrecord-{N:05}-of-{config["shards"]:05}.gz') for N in range(config['shards']))
    output:
        vcf = get_dir('output','{animal}.bwa.{callset}.vcf.gz'),
        gvcf = get_dir('output','{animal}.bwa.{callset}.g.vcf.gz')
    params:
        gvcf = lambda wildcards,input: PurePath(input.gvcf[0]).with_suffix('').with_suffix(f'.tfrecord@{config["shards"]}.gz'),
        singularity_call = lambda wildcards, input, output: make_singularity_call(wildcards,extra_args=f'-B {PurePath(input.ref[0]).parent}:/reference/ -B {PurePath(output["vcf"]).parent}:/{PurePath(output["vcf"]).parent}'),
        ref = lambda wildcards,input: f'/reference/{PurePath(input.ref[0]).name}',
        contain = lambda wildcards: config['DV_container']
    threads: 1
    resources:
        mem_mb = config.get('resources',{}).get('postprocess',{}).get('mem_mb',40000),
        walltime = config.get('resources',{}).get('postprocess',{}).get('walltime','4:00'),#get_walltime,
        disk_scratch = 1,
        use_singularity = True
    priority: 100
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

rule bcftools_scatter:
    input:
        get_dir('output','{animal}.bwa.regions.g.vcf.gz')
    output:
        vcf = (get_dir('splits','{animal}_{region}.vcf.gz',region=R) for R in REGIONS),
        tbi = (get_dir('splits','{animal}_{region}.vcf.gz.tbi',region=R) for R in REGIONS)
    params:
        prefix = '{animal}_',
        scatter = ','.join(REGIONS)
    threads: 4
    resources:
        mem_mb = 2000,
        walltime = '60'
    shell:
        '''
        bcftools +scatter {input} --threads {threads} -s {params.scatter} -Oz -o {output} -p {params.prefix}
        vcf=({output.vcf})
        parallel -j {threads} -u tabix -fp {{}} ::: ${{vcf[@]}}
        '''

def get_GL_config(preset):
    if preset == 'WGS':
        return 'DeepVariantWGS'
    elif preset == 'Unfiltered':
        return 'DeepVariant_unfiltered'
    elif 'GL_config' in config:
        return '/data/' + config['GL_config'][preset]
    else:
        raise ValueError(f'Unknown config {preset=}')

def get_merging_input(wildcards,ext=''):
    if 'regions_bed' in config: #wildcards.callset == 'regions':
        return expand(get_dir('splits',f'{{animal}}_{{region}}.vcf.gz{ext}'),animal=config['animals'],region=wildcards.callset)
    else:
        return expand(get_dir('output','{animal}.bwa.{callset}.g.vcf.gz'),animal=config['animals'],callset=wildcards.callset)

rule GLnexus_merge:
    input:
        vcf = lambda wildcards: get_merging_input(wildcards),
        tbi = lambda wildcards: get_merging_input(wildcards,'.tbi')
        #vcf = (get_dir('splits','{animal}_{region}.g.vcf.gz',animal=ANIMAL) for ANIMAL in config['animals']),
        #tbi = (get_dir('splits','{animal}_{region}.g.vcf.gz.tbi',animal=ANIMAL) for ANIMAL in config['animals'])
    output:
        multiext(get_dir('main','{cohort}.{callset}.{preset}.vcf.gz'),'','.tbi')
    params:
        gvcfs = lambda wildcards, input: list('/data' / PurePath(fpath) for fpath in input.vcf),
        out = lambda wildcards, output: '/data' / PurePath(output[0]),
        DB = lambda wildcards, output: f'/tmp/GLnexus.DB',
        preset = lambda wildcards: get_GL_config(wildcards.preset),
        singularity_call = lambda wildcards: make_singularity_call(wildcards,'-B .:/data', input_bind=False, output_bind=False, work_bind=False),
        mem = lambda wildcards,threads,resources: threads*resources['mem_mb']/1024,
        container = config['GL_container']
    threads: lambda wildcards: 24 if False == 'all' else 12 #force using 4 threads for bgziping
    resources:
        mem_mb = 8000,
        disk_scratch = 50,
        walltime = '4:00',
        use_singularity = True
    priority: 25
    shell:
        '''
        ulimit -Sn 4096
        {params.singularity_call} \
        {params.container} \
        /bin/bash -c " /usr/local/bin/glnexus_cli \
        --dir {params.DB} \
        --config {params.preset} \
        --threads {threads} \
        --mem-gbytes {params.mem} \
        {params.gvcfs} \
        | bcftools view - | bgzip -@ 4 -c > {params.out}"
        tabix -p vcf {output[0]}
        '''
