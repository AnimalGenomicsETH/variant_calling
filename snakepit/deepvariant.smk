
#NOTE: may need to be updated if deepvariant changes its internal parameters.
def make_custom_example_arguments(model):
    match model:
        case 'RNA' | 'WES':
            return '--channels \'\' --split_skip_reads'
        case 'PACBIO':
            return '--add_hp_channel --alt_aligned_pileup "diff_channels" --max_reads_per_partition "600" --min_mapping_quality "1" --parse_sam_aux_fields --partition_size "25000" --phase_reads --pileup_image_width "199" --norealign_reads --sort_by_haplotypes --track_ref_reads --vsc_min_fraction_indels "0.12"'
        case 'WGS':
            return '--channels "insert_size"'
        case _:
            return '--channels \'\' --split_skip_reads'

def get_checkpoint(model):
    match model:
        case 'PACBIO' | 'WGS':
            return f'/opt/models/{model.lower()}'
        case 'hybrid':
            return '/opt/models/hybrid_pacbio_illumina'
        case _:
            return config['model']

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

def get_sample_bam_path(wildcards):
    # query alignment_metadata for path?
    # and index
    return ''

#error can be caused by corrupted file, check gzip -t -v
rule deepvariant_make_examples:
    input:
        reference = multiext(config['reference'],'','.fai'),
        bam =  get_sample_bam_path
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
        mem_mb_per_cpu = config.get('resources',{}).get('make_examples',{}).get('mem_mb',6000),
        runtime = config.get('resources',{}).get('make_examples',{}).get('runtime','4h'),
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
        examples = expand('{run}/deepvariant/intermediate_results_{sample}_{region}/make_examples.tfrecord-{N}-of-{sharding}.gz',sharding=f'{config["shards"]:05}',N=[f'{i:05}' for i in range(config['shards'])],allow_missing=True),
        json = expand('{run}/deepvariant/intermediate_results_{sample}_{region}/make_examples.tfrecord-{N}-of-{sharding}.gz.example_info.json',sharding=f'{config["shards"]:05}',N=[f'{i:05}' for i in range(config['shards'])],allow_missing=True),
    output:
        temp('{run}/deepvariant/intermediate_results_{sample}_{region}/call_variants_output-00000-of-00001.tfrecord.gz')
    params:
        examples = lambda wildcards, input: PurePath(input['examples'][0]).with_suffix('').with_suffix(f'.tfrecord@{config["shards"]}.gz'),
        model = get_checkpoint(config['model'])
    threads: config.get('resources',{}).get('call_variants',{}).get('threads',12)
    resources:
        mem_mb_per_cpu = config.get('resources',{}).get('call_variants',{}).get('mem_mb',1500),
        runtime = config.get('resources',{}).get('call_variants',{}).get('runtime','24h'),
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
        gvcf = lambda wildcards,input: PurePath(input.gvcf[0]).with_suffix('').with_suffix(f'.tfrecord@{config["shards"]}.gz')
    threads: config.get('resources',{}).get('postprocess',{}).get('threads',2)
    resources:
        mem_mb_per_cpu = config.get('resources',{}).get('postprocess',{}).get('mem_mb',20000),
        runtime = config.get('resources',{}).get('postprocess',{}).get('runtime','4h'),
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
--cpus {threads}
        '''

#can we scatter off of a better BED file?
rule bcftools_scatter:
    input:
        gvcf = expand(rules.deepvariant_postprocess.output['gvcf'],region='all',allow_missing=True),
        regions = '/cluster/work/pausch/vcf_UCD/2023_07/regions.bed' # fix to be general
    output:
        expand('{run}/deepvariant/{sample}.{region}.g.vcf.gz',region=regions,allow_missing=True),
        expand('{run}/deepvariant/{sample}.{region}.g.vcf.gz.csi',region=regions,allow_missing=True)
    params:
        regions = regions,
        _dir = lambda wildcards, output: PurePath(output[0]).with_suffix('').with_suffix('').with_suffix('').with_suffix('')
    threads: 2
    resources:
        mem_mb_per_cpu = 2500,
        runtime = '1h'
    shell:
        '''
bcftools +scatter {input.gvcf[0]} -o {params._dir} -Oz --threads {threads} --write-index -S {input.regions} -x unplaced --no-version
for R in {params.regions}
do 
    mv {params._dir}/$R.vcf.gz {params._dir}.$R.g.vcf.gz
    mv {params._dir}/$R.vcf.gz.csi {params._dir}.$R.g.vcf.gz.csi
done
        '''

# add default to repo under config
def get_GL_config(preset):
    match preset:
        case 'DeepVariantWGS' | 'DeepVariant_unfiltered' | 'DeepVariantWES_MED_DP':
            return preset
        case _:
            return config['GL_config'][preset]

rule GLnexus_merge:
    input:
        #gvcfs = expand('{run}/deepvariant/{region}/{sample}.g.vcf.gz',sample=cohort_samples,allow_missing=True),
        #tbi = expand('{run}/deepvariant/{region}/{sample}.g.vcf.gz.csi',sample=cohort_samples,allow_missing=True)
        gvcfs = expand('{run}/deepvariant/{sample}.{region}.g.vcf.gz',sample=samples,allow_missing=True),
        tbi = expand('{run}/deepvariant/{sample}.{region}.g.vcf.gz.tbi',sample=samples,allow_missing=True)
    output:
        multiext('{run}/{region}.{preset}.vcf.gz','','.tbi')
    params:
        preset = lambda wildcards: get_GL_config(wildcards.preset),
        mem = lambda wildcards,threads,resources: threads*resources['mem_mb']/1024
    threads: config.get('resources',{}).get('merge',{}).get('threads',12),
    resources:
        mem_mb_per_cpu = config.get('resources',{}).get('merge',{}).get('mem_mb',6000),
        runtime = config.get('resources',{}).get('merge',{}).get('runtime','4h'),
    container: config.get('containers',{}).get('GLnexus','docker://ghcr.io/dnanexus-rnd/glnexus:v1.4.3')
    shell:
        '''
/usr/local/bin/glnexus_cli \
--dir $TMPDIR/GLnexus.DB \
--config {params.preset} \
--threads {threads} \
--mem-gbytes {params.mem} \
{input.gvcfs} |\
bcftools view - |\
bgzip -@ 8 -c > {output[0]} 

tabix -p vcf {output[0]}
        '''
