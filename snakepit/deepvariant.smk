from pathlib import Path,PurePath

wildcard_constraints:
    region = r'\w+',
    run = r'\w+',
    preset = r'|'.join(config['GL_config'])
    
config.setdefault('Run_name', 'DV_variant_calling')

if 'binding_paths' in config:
    for path in config['binding_paths']:
        workflow._singularity_args += f' -B {path}'

regions = config.get('regionsz',list(map(str,range(1,30))) + ['X','Y','MT','unplaced'])

if config.get('scatter_merge',True):
    ruleorder: bcftools_scatter > deepvariant_postprocess
else:
    ruleorder: deepvariant_postprocess > bcftools_scatter

rule all:
    input:
        expand('{name}/{region}.{preset}.vcf.gz',name=config['Run_name'],region=regions,preset=config['GL_config']),
        expand('{name}/{region}.{preset}.vcf.gz.tbi',name=config['Run_name'],region=regions,preset=config['GL_config'])

cohort_samples = config['samples'] if 'glob' not in config['samples'] else glob_wildcards(config["bam_path"] + config["samples"]["glob"]).sample

#NOTE: may need to be updated if deepvariant changes its internal parameters.
def make_custom_example_arguments(model):
    match model:
        case 'RNA' | 'WES':
            return '--channels \'\' --split_skip_reads'
        case 'PACBIO':
            return '--alt_aligned_pileup "diff_channels" --max_reads_per_partition "600" --min_mapping_quality "1" --parse_sam_aux_fields --partition_size "25000" --phase_reads --pileup_image_width "147" --norealign_reads --sort_by_haplotypes --track_ref_reads --vsc_min_fraction_indels "0.12" --trim_reads_for_pileup'
        case 'WGS':
            return '--channels "insert_size"'
        case 'ONT_R104':
            return '--alt_aligned_pileup "diff_channels" --max_reads_per_partition "600" --min_mapping_quality "5" --parse_sam_aux_fields --partition_size "25000" --phase_reads --pileup_image_width "99" --norealign_reads --sort_by_haplotypes --track_ref_reads --vsc_min_fraction_snps "0.08" --vsc_min_fraction_indels "0.12" --trim_reads_for_pileup'
        case 'MASSEQ':
            return '--alt_aligned_pileup "diff_channels" --max_reads_per_partition 0 --min_mapping_quality 1 --parse_sam_aux_fields --partition_size "25000" --phase_reads --pileup_image_width "199" --norealign_reads --sort_by_haplotypes --track_ref_reads --vsc_min_fraction_indels 0.12 --trim_reads_for_pileup --max_reads_for_dynamic_bases_per_region "1500"'
        case _:
            return '--channels \'\' --split_skip_reads'

def get_regions(region):
    if region == 'all':
        return ''
    if region not in config['regions']:
        return ''
    elif isinstance(config['regions'][region],str) and Path(config['regions'][region]).exists():
        return f' --regions {config["regions"][region]}' #+ '"' + " ".join([l.strip() for l in open(config['regions'][region])]) + '"'
    elif isinstance(config['regions'][region],str):
        return f' --regions "{config["regions"][region]}"'
    else:
        return ' --regions ' + '"' + f'{" ".join(map(str,config["regions"][region]))}' + '"'

BAM_EXT = config.get('bam_index','.bai')

#error can be caused by corrupted file, check gzip -t -v
rule deepvariant_make_examples:
    input:
        reference = multiext(config['reference'],'','.fai'),
        bam = multiext(config['bam_path']+config['bam_name'],'',BAM_EXT)
    output:
        examples = temp('{run}/deepvariant/intermediate_results_{sample}_{region}/make_examples.tfrecord-{N}-of-{sharding}.gz'),
        gvcf = temp('{run}/deepvariant/intermediate_results_{sample}_{region}/gvcf.tfrecord-{N}-of-{sharding}.gz'),
        info = temp('{run}/deepvariant/intermediate_results_{sample}_{region}/make_examples.tfrecord-{N}-of-{sharding}.gz.example_info.json')
    params:
        examples = lambda wildcards, output: PurePath(output['examples']).with_suffix('').with_suffix(f'.tfrecord@{config["shards"]}.gz'),
        gvcf = lambda wildcards, output: PurePath(output['gvcf']).with_suffix('').with_suffix(f'.tfrecord@{config["shards"]}.gz'),
        model_args = make_custom_example_arguments(config['model']),
        regions = lambda wildcards: get_regions(wildcards.region) 
    threads: 1
    resources:
        mem_mb = config.get('resources',{}).get('make_examples',{}).get('mem_mb',6000),
        walltime = config.get('resources',{}).get('make_examples',{}).get('walltime','4h'),
        scratch = "1G",
        storage_load = 1
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

def get_checkpoint(model):
    match model:
        case 'PACBIO' | 'WGS':
            return f'/opt/models/{model.lower()}'
        case 'hybrid':
            return '/opt/models/hybrid_pacbio_illumina'
        case _:
            return config['model']

rule deepvariant_call_variants:
    input:
        examples = expand('{run}/deepvariant/intermediate_results_{sample}_{region}/make_examples.tfrecord-{N}-of-{sharding}.gz',sharding=f'{config["shards"]:05}',N=[f'{i:05}' for i in range(config['shards'])],allow_missing=True),
        json = expand('{run}/deepvariant/intermediate_results_{sample}_{region}/make_examples.tfrecord-{N}-of-{sharding}.gz.example_info.json',sharding=f'{config["shards"]:05}',N=[f'{i:05}' for i in range(config['shards'])],allow_missing=True),
    output:
        temp('{run}/deepvariant/intermediate_results_{sample}_{region}/call_variants_output-00000-of-00001.tfrecord.gz')
    params:
        examples = lambda wildcards,input: PurePath(input['examples'][0]).with_suffix('').with_suffix(f'.tfrecord@{config["shards"]}.gz'),
        model = get_checkpoint(config['model'])
    threads: config.get('resources',{}).get('call_variants',{}).get('threads',12)
    resources:
        mem_mb = config.get('resources',{}).get('call_variants',{}).get('mem_mb',1500),
        walltime = config.get('resources',{}).get('call_variants',{}).get('walltime','24h'),
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
        gvcf = expand('{run}/deepvariant/intermediate_results_{sample}_{region}/gvcf.tfrecord-{N}-of-{sharding}.gz',sharding=f'{config["shards"]:05}',N=[f'{i:05}' for i in range(config['shards'])],allow_missing=True)
    output:
        vcf = multiext('{run}/deepvariant/{sample}.{region}.vcf.gz','','.tbi'),
        gvcf = multiext('{run}/deepvariant/{sample}.{region}.g.vcf.gz','','.tbi')
    params:
        variants = lambda wildcards,input: PurePath(input.variants).with_name('call_variants_output@1.tfrecord.gz'),
        gvcf = lambda wildcards,input: PurePath(input.gvcf[0]).with_suffix('').with_suffix(f'.tfrecord@{config["shards"]}.gz'),
        handle_sex_chromosomes = f'--haploid_contigs "X,Y" --par_regions_bed "{config["haploid_sex"]}"' if 'haploid_sex' in config else ''
    threads: config.get('resources',{}).get('postprocess',{}).get('threads',2)
    resources:
        mem_mb = config.get('resources',{}).get('postprocess',{}).get('mem_mb',20000),
        walltime = config.get('resources',{}).get('postprocess',{}).get('walltime','4h'),
        scratch = "1G"
    priority: 100
    container: config.get('containers',{}).get('DV','docker://google/deepvariant:latest')
    shell:
        '''
        /opt/deepvariant/bin/postprocess_variants \
        --ref {input.reference[0]} \
        --infile {params.variants} \
        --outfile {output.vcf[0]} \
        --gvcf_outfile {output.gvcf[0]} \
        --nonvariant_site_tfrecord_path {params.gvcf} \
        --novcf_stats_report \
        {params.handle_sex_chromosomes} \
        --cpus {threads}
        '''

rule bcftools_scatter:
    input:
        gvcf = expand(rules.deepvariant_postprocess.output['gvcf'],region='all',allow_missing=True),
        regions = '/cluster/work/pausch/vcf_UCD/2023_07/regions.bed',
        fai = config['reference'] + ".fai"
    output:
        expand('{run}/deepvariant/{sample}.{region}.g.vcf.gz',region=regions,allow_missing=True),
        expand('{run}/deepvariant/{sample}.{region}.g.vcf.gz.tbi',region=regions,allow_missing=True)
    wildcard_constraints:
        region = "(?!all)"
    params:
        regions = regions,
        _dir = lambda wildcards, output: PurePath(output[0]).with_suffix('').with_suffix('').with_suffix('').with_suffix('')
    threads: 2
    resources:
        mem_mb = 2500,
        walltime = '1h'
    shell:
        '''
        bcftools +scatter {input.gvcf[0]} -o {params._dir} -Oz --threads {threads} -S {input.regions} -x unplaced --no-version
        for R in {params.regions}
        do 
          bcftools reheader -f {input.fai} {params._dir}/$R.vcf.gz > {params._dir}.$R.g.vcf.gz
          tabix -p vcf {params._dir}.$R.g.vcf.gz
        done
        '''

def get_GL_config(preset):
    match preset:
        case 'DeepVariantWGS' | 'DeepVariant_unfiltered' | 'DeepVariantWES_MED_DP':
            return preset
        case _:
            return config['GL_config'][preset]

rule GLnexus_merge:
    input:
        gvcfs = expand('{run}/deepvariant/{sample}.{region}.g.vcf.gz',sample=cohort_samples,allow_missing=True),
        tbi = expand('{run}/deepvariant/{sample}.{region}.g.vcf.gz.tbi',sample=cohort_samples,allow_missing=True)
    output:
        multiext('{run}/{region}.{preset}.vcf.gz','','.tbi')
    params:
        preset = lambda wildcards: get_GL_config(wildcards.preset),
        mem = lambda wildcards,threads,resources: threads*resources['mem_mb']/1024
    threads: config.get('resources',{}).get('merge',{}).get('threads',12),
    resources:
        mem_mb = config.get('resources',{}).get('merge',{}).get('mem_mb',6000),
        walltime = config.get('resources',{}).get('merge',{}).get('walltime','4h'),
        scratch = '50G'
    priority: 25
    shell:
        '''
        glnexus_cli \
        --dir $TMPDIR/GLnexus.DB \
        --config {params.preset} \
        --threads {threads} \
        --mem-gbytes {params.mem} \
        {input.gvcfs} |\
        bcftools view -@ {threads} -W=tbi {output[0]} 
        '''
