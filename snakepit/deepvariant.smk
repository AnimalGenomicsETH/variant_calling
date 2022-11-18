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
            regions.append([line.split()[0] for line in fai])
    return regions


wildcard_constraints:
    haplotype = r'asm|hap1|hap2|parent1|parent2',
    phase = r'unphased|phased',
    model = r'pbmm2|hybrid|bwa|mm2',
    caller = r'DV|GATK',
    cohort = r'\w+',
    region = r'\w+',
    run = r'\w+'
#    preset = rf'{config.get("GL_config","WGS")}'
    
def get_model(model,base='/opt/models',ext='model.ckpt'):
    model_location = f'{base}/{{}}/{ext}'
    match model:
        case 'PACBIO' | 'mm2':
            return model_location.format('pacbio')
        case 'hybrid':
            return model_location.format('hybrid_pacbio_illumina')
        case 'WGS':
            return model_location.format('wgs')
        case _:
            return model

#for dir_ in ('output','work'):
#    for animal,reference in product(config['animals'],config['reference']):
#        Path(get_dir(dir_,animal=animal,ref=reference,**config)).mkdir(exist_ok=True)

REGIONS = determine_run_regions()

config.setdefault('Run_name', 'deepvariant')

if True:
    print('Binding in relevant paths for singularity')
    workflow.singularity_args = f'-B $TMPDIR -B {config["bam_path"]} -B {PurePath(config["reference"]).parent}'
    for preset in config.get('GL_config',[]):
        if Path(config['GL_config'][preset]).exists():
            workflow.singularity_args += f' -B {PurePath(config["GL_config"][preset]).parent}'
    if Path(config['model']).exists():
        workflow.singularity_args += f' -B {PurePath(config["model"]).parent}'

print(workflow.singularity_args)

def get_regions(region):
    match region:
        case 'bed':
            return '--regions '

        case 'all' | _:
            return ''

rule all:
    input:
        expand('{name}/{region}.{preset}.vcf.gz',name=config['Run_name'],region=config['regions'],preset=config['GL_config'])
        #expand('{cohort}/{chromosome}.{preset}.vcf.gz{ext}',ext=('','.tbi'),cohort=config['Run_name'],chromosome=get_regions(),preset=config.get("GL_config","WGS"))

        #expand(get_dir('main','{cohort}.{regions}.{preset}.vcf.gz{ext}'),ext=('','.tbi'),cohort=config.get("cohort","cohort"),regions=REGIONS,preset=config.get("GL_config","WGS"))

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
        bam = multiext(config['bam_path']+'{sample}/{sample}.bam','','.bai')
    output:
        examples = temp('{run}/deepvariant/intermediate_results_{sample}_{region}/make_examples.tfrecord-{N}-of-{sharding}.gz'),
        gvcf = temp('{run}/deepvariant/intermediate_results_{sample}_{region}/gvcf.tfrecord-{N}-of-{sharding}.gz')
    params:
        examples = lambda wildcards, output: PurePath(output['examples']).with_suffix('').with_suffix(f'.tfrecord@{config["shards"]}.gz'),
        gvcf = lambda wildcards, output: PurePath(output['gvcf']).with_suffix('').with_suffix(f'.tfrecord@{config["shards"]}.gz'),
        model_args = make_custom_example_arguments(config['model']),
        regions = lambda wildcards: get_regions(wildcards.region) 
    threads: 1
    resources:
        mem_mb = config.get('resources',{}).get('make_examples',{}).get('mem_mb',6000),
        walltime = config.get('resources',{}).get('make_examples',{}).get('walltime','4:00'),
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
        expand('{{run}}/deepvariant/intermediate_results_{{sample}}_{{region}}/make_examples.tfrecord-{N}-of-{sharding}.gz',sharding=f'{config["shards"]:05}',N=[f'{i:05}' for i in range(config['shards'])])
    output:
        temp('{run}/deepvariant/intermediate_results_{sample}_{region}/call_variants_output.tfrecord.gz')
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
        --outfile {output} \
        --examples {params.examples} \
        --checkpoint {params.model}
        '''

rule deepvariant_postprocess:
    input:
        reference = multiext(config['reference'],'','.fai'),
        variants = rules.deepvariant_call_variants.output[0],
        gvcf = expand('{{run}}/deepvariant/intermediate_results_{{sample}}_{{region}}/gvcf.tfrecord-{N}-of-{sharding}.gz',sharding=f'{config["shards"]:05}',N=[f'{i:05}' for i in range(config['shards'])])
    output:
        vcf = '{run}/deepvariant/{sample}.{region}.vcf.gz',
        gvcf = '{run}/deepvariant/{sample}.{region}.g.vcf.gz'
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
        '{run}/deepvariant/{sample}.{region}.g.vcf.gz'
        #rules.deepvariant_postprocess.output[1]
    output:
        vcf = '3'#expand('{{run}}/scattered/{{sample}}_{region}.vcf.gz',region=config['regions']),
    params:
        prefix = '{sample}_',
        out = lambda wildcards,output: 3,#PurePath(output[0]).parent,
        scatter = ','.join(['a','b'])
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
    match preset:
        case 'WGS':
            return 'DeepVariantWGS'
        case 'Unfiltered':
            return 'DeepVariant_unfiltered'
        case _:
            return config['GL_config'][preset]

def get_merging_input(wildcards,ext=''):
    #if 'regions_bed' in config: #wildcards.callset == 'regions':
    #    return expand(get_dir('splits',f'{{animal}}_{{region}}.vcf.gz{ext}'),animal=config['animals'],region=wildcards.callset)
    #else:
    return expand('{{run}}/deepvariant/{sample}.{{region}}.g.vcf.gz',sample=config['samples'])
    

rule GLnexus_merge:
    input:
        #vcf = expand('{{run}}/deepvariant/{sample}.{{region}}.g.vcf.gz',sample=config['samples'])
        gvcfs = lambda wildcards: get_merging_input(wildcards),
        tbi = lambda wildcards: get_merging_input(wildcards,'.tbi')
    output:
        multiext('{run}/{region}.{preset}.vcf.gz','','.tbi')
    params:
        preset = lambda wildcards: get_GL_config(wildcards.preset),
        mem = lambda wildcards,threads,resources: threads*resources['mem_mb']/1024
    threads: config.get('resources',{}).get('merge',{}).get('threads',12),
    resources:
        mem_mb = config.get('resources',{}).get('merge',{}).get('mem_mb',6000),
        walltime = config.get('resources',{}).get('merge',{}).get('walltime','4:00'),
        scratch = '50G'
    priority: 25
    container: config.get('containers',{}).get('GLnexus','docker://ghcr.io/dnanexus-rnd/glnexus:v1.4.3')
    shell:
        '''
        /usr/local/bin/glnexus_cli \
        --dir $TMPDIR/GLnexus.DB \
        --config {params.preset} \
        --threads {threads} \
        --mem-gbytes {params.mem} \
        {input.gvcfs} |\
        bcftools view --threads 4 -o {output[0]} - 

        tabix -p vcf {output[0]}
        '''
