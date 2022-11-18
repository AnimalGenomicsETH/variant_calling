from pathlib import Path,PurePath
from itertools import product

if 'regions_bed' in config:
    ruleorder: bcftools_scatter > GLnexus_merge

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


wildcard_constraints:
    haplotype = r'asm|hap1|hap2|parent1|parent2',
    phase = r'unphased|phased',
    model = r'pbmm2|hybrid|bwa|mm2',
    caller = r'DV|GATK',
    cohort = r'\w+',
#    preset = rf'{config.get("GL_config","WGS")}'
    
def get_model(model,base='/opt/models',ext='model.ckpt'):
    model_location = f'{base}/{{}}/{ext}'
    match model:
        case ['pbmm2','mm2']:
            return model_location.format('pacbio')
        case 'hybrid':
            return model_location.format('hybrid_pacbio_illumina')
        case 'bwa':
            return model_location.format('wgs')
        case _:
            return model

for dir_ in ('output','work'):
    for animal,reference in product(config['animals'],config['reference']):
        Path(get_dir(dir_,animal=animal,ref=reference,**config)).mkdir(exist_ok=True)

REGIONS = determine_run_regions()

RUN_NAME = config.get('Run_name','deepvariant')

rule all:
    input:
        expand(get_dir('main','{cohort}.{regions}.{preset}.vcf.gz{ext}'),ext=('','.tbi'),cohort=config.get("cohort","cohort"),regions=REGIONS,preset=config.get("GL_config","WGS"))

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


#NOTE: may need to be updated if deepvariant changes its internal parameters.
def make_custom_example_arguments(model):
    match model:
        case 'RNA':
            return '--channels \'\' --split_skip_reads'
        case 'PACBIO':
            return '--add_hp_channel --alt_aligned_pileup "diff_channels" --max_reads_per_partition "600" --min_mapping_quality "1" --parse_sam_aux_fields --partition_size "25000" --phase_reads --pileup_image_width "199" --norealign_reads --sort_by_haplotypes --track_ref_reads --vsc_min_fraction_indels "0.12"'
        case 'WGS':
            return '--channels "insert_size"'

#error can be caused by corrupted file, check gzip -t -v
rule deepvariant_make_examples:
    input:
        reference = multiext(config['reference'],'','.fai'),
        bam = multiext(config['bam_path']+'{sample}.bam'),'','.bai')
    output:
        examples = temp('deepvariant/intermediate_results_{sample}/make_examples.{callset}.tfrecord-{N}-of-{sharding}.gz'),
        gvcf = temp('deepvariant/intermediate_results_{sample}/gvcf.{callset}.tfrecord-{N}-of-{sharding}.gz')
    params:
        examples = lambda wildcards, output: PurePath(output['examples']).with_suffix('').with_suffix(f'.tfrecord@{config["shards"]}.gz'),
        gvcf = lambda wildcards, output: PurePath(output['gvcf']).with_suffix('').with_suffix(f'.tfrecord@{config["shards"]}.gz'),
        model_args = make_custom_example_arguments(config['model']),
        regions = lambda wildcards: get_regions(wildcards) 
    threads: 1
    resources:
        mem_mb = config.get('resources',{}).get('make_examples',{}).get('mem_mb',6000),
        walltime = config.get('resources',{}).get('make_examples',{}).get('walltime','4:00'),#get_walltime,
        scratch = "1G"
    priority: 50
    container: config.get('containers',{}).get('DV','docker://google/deepvariant:latest')
    shell:
        '''
        /opt/deepvariant/bin/make_examples \
        --mode calling \
        --ref {input.reference[0]} \
        --include_med_dp \
        --reads {input.bam[0]} \
        --examples {params.examples} \
        --gvcf {params.gvcf} \
        --sample_name {wildcards.sample} \
        {params.regions} \
        {params.model_args} \
        --task {wildcards.N}
        '''

rule deepvariant_call_variants:
    input:
        expand('deepvariant/intermediate_results_{{sample}}/make_examples.{{callset}}.tfrecord-{N}-of-{sharding}.gz',sharding=f'{config["shards"]:05}',N=[f'{i:05}' for i in range(config['shards'])])
    output:
        temp('deepvariant/intermediate_results_{sample}/call_variants_output.{callset}.tfrecord.gz')
    params:
        examples = lambda wildcards,input: PurePath(input[0]).with_suffix('').with_suffix(f'.tfrecord@{config["shards"]}.gz'),
        model = get_model(config['model'])
    threads: config.get('resources',{}).get('call_variants',{}).get('threads',12)
    resources:
        mem_mb = config.get('resources',{}).get('call_variants',{}).get('mem_mb',1500),
        walltime = config.get('resources',{}).get('call_variants',{}).get('walltime','24:00'),
        scratch = "1G"
    priority: 90
    container: config.get('containers',{}).get('DV','docker://google/deepvariant:latest')
    shell:
        '''
        /opt/deepvariant/bin/call_variants \
        --outfile /{output} \
        --examples /{params.examples} \
        --checkpoint {params.model}
        '''

rule deepvariant_postprocess:
    input:
        reference = multiext(config['reference'],'','.fai'),
        rules.deepvariant_call_variants.output[0],
        gvcf = expand('deepvariant/intermediate_results_{{sample}}/gvcf.{{callset}}.tfrecord-{N}-of-{sharding}.gz',sharding=f'{config["shards"]:05}',N=[f'{i:05}' for i in range(config['shards'])])
    output:
        vcf = 'deepvariant/{sample}.bwa.{callset}.vcf.gz',
        gvcf = get_dir('output','{sample}.bwa.{callset}.g.vcf.gz')
    params:
        gvcf = lambda wildcards,input: PurePath(input.gvcf[0]).with_suffix('').with_suffix(f'.tfrecord@{config["shards"]}.gz')
    threads: 1
    resources:
        mem_mb = config.get('resources',{}).get('postprocess',{}).get('mem_mb',40000),
        walltime = config.get('resources',{}).get('postprocess',{}).get('walltime','4:00'),
        scratch = "1G"
    priority: 100
    container: config.get('containers',{}).get('DV','docker://google/deepvariant:latest')
    shell:
        '''
        /opt/deepvariant/bin/postprocess_variants \
        --ref {input.reference[0]} \
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
        out = lambda wildcards,output: PurePath(output[0]).parent,
        scatter = ','.join(REGIONS)
    threads: 2
    resources:
        mem_mb = 2000,
        walltime = '60'
    shell:
        '''
        bcftools +scatter {input} --threads {threads} -s {params.scatter} -Oz -o {params.out} -p {params.prefix}
        vcf=({output.vcf})
        parallel -j {threads} -u tabix -fp vcf {{}} ::: ${{vcf[@]}}
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
        preset = lambda wildcards: get_GL_config(wildcards.preset),
        mem = lambda wildcards,threads,resources: threads*resources['mem_mb']/1024
    threads: lambda wildcards: 24 if False == 'all' else 12 #force using 4 threads for bgziping
    resources:
        mem_mb = 8000,
        disk_scratch = 50,
        walltime = '4:00',
        use_singularity = True
    priority: 25
    config.get('containers',{}).get('GLnexus','docker://google/deepvariant:latest')
    shell:
        '''
        /usr/local/bin/glnexus_cli \
        --dir $TMPDIR/GLnexus.DB \
        --config {params.preset} \
        --threads {threads} \
        --mem-gbytes {params.mem} \
        {params.gvcfs} |\
        bcftools view - | bgzip -@ 4 -c > {params.out}

        tabix -p vcf {output[0]}
        '''
