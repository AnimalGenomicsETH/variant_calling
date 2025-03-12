#NOTE: may need to be updated if deepvariant changes its internal parameters.
def make_custom_example_arguments(model,small=True): #TODO: allow switching on/off small model
    common_long_read_options = '--alt_aligned_pileup "diff_channels" --parse_sam_aux_fields --partition_size "25000" --phase_reads --norealign_reads --sort_by_haplotypes --track_ref_reads --trim_reads_for_pileup --vsc_min_fraction_indels "0.12" '
    small_model = '--small_model_snp_gq_threshold 25 --small_model_indel_gq_threshold 30 --small_model_vaf_context_window_size 51 ' #TODO: this does not work for ONTr10
    match model:
        case 'WGS':
            return small_model + '--channel_list BASE_CHANNELS,insert_size'
        case 'PACBIO':
            return common_long_read_options + small_model + '--max_reads_per_partition "600" --min_mapping_quality "1" --pileup_image_width "147"'
        case 'MASSEQ':
            return common_long_read_options + '--max_reads_per_partition 0 --min_mapping_quality 1 --pileup_image_width "199" --max_reads_for_dynamic_bases_per_region "1500"'
        case 'ONT_R104':
            return common_long_read_options + '--max_reads_per_partition "600" --min_mapping_quality "5" --pileup_image_width "99" --vsc_min_fraction_snps "0.08"'
        case 'RNA' | 'WES' | _:
            return '--channel_list BASE_CHANNELS --split_skip_reads'

#TODO: need to update models?
def get_checkpoint(model,small=True):
    path = f'/opt/{"small" if small else ""}models/'
    match model:
        case 'PACBIO' | 'MASSEQ' | 'WGS' | 'ONT_R104' :
            return path + model.lower()
        case 'hybrid':
            return path + 'hybrid_pacbio_illumina'
        case _:
            return config['model']

def get_regions(wildcards):
    if wildcards.region == 'all':
        return ''
    else:
        return '--regions "{' '.join(region_map[wildcards.region])}"'

## implictly assume indexed bam/cram files, snakemake won't complain but jobs may fail.
def get_sample_bam_path(wildcards):
    return alignment_metadata.filter(pl.col('sample ID')==wildcards.sample).get_column('bam path').to_list()[0]

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
        regions = get_regions
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
        gvcf = lambda wildcards,input: PurePath(input.gvcf[0]).with_suffix('').with_suffix(f'.tfrecord@{config["shards"]}.gz'),
        handle_sex_chromosomes = '' if 'PAR_regions' not in config else f'--haploid_contigs "X,Y" --par_regions_bed "{config['PAR_regions']}"'
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
{params.handle_sex_chromosomes} \
--cpus {threads}
        '''

rule bcftools_scatter:
    input:
        gvcf = expand(rules.deepvariant_postprocess.output['gvcf'],region='all',allow_missing=True),
        regions = config.get('regions',f"{config['reference']}.fai"),
        fai = config['reference'] + ".fai"
    output:
        expand('{run}/deepvariant/{sample}.{region}.g.vcf.gz',region=regions,allow_missing=True),
        expand('{run}/deepvariant/{sample}.{region}.g.vcf.gz.tbi',region=regions,allow_missing=True)
    wildcard_constraints:
        region = "(?!all)"
    params:
        regions = ' '.join(regions),
        region_cols = lambda wildcards, input: f"<(cut -f {2 if 'regions' in config else 1} {input.regions})",
        _dir = lambda wildcards, output: PurePath(output[0]).with_suffix('').with_suffix('').with_suffix('').with_suffix('')
    threads: 2
    resources:
        mem_mb_per_cpu = 2500,
        runtime = '1h'
    shell:
        '''
        bcftools +scatter {input.gvcf[0]} -o $TMPDIR -Oz --threads {threads} -S {params.region_cols} -x unplaced --no-version
        for R in {params.regions}
        do 
          bcftools reheader -f {input.fai} $TMPDIR/$R.vcf.gz > {params._dir}.$R.g.vcf.gz
          tabix -p vcf {params._dir}.$R.g.vcf.gz
        done
        '''

# see https://github.com/dnanexus-rnd/GLnexus/pull/310
def get_GL_config(preset):
    match preset:
        case 'DeepVariantWGS' | 'DeepVariant_unfiltered' | 'DeepVariantWES_MED_DP':
            return preset
        case _:
            return config['GL_config'][preset]

rule GLnexus_merge:
    input:
        gvcfs = expand('{run}/deepvariant/{sample}.{region}.g.vcf.gz',sample=samples,allow_missing=True),
        tbi = expand('{run}/deepvariant/{sample}.{region}.g.vcf.gz.tbi',sample=samples,allow_missing=True)
    output:
        multiext('{run}/{region}.Unrevised.vcf.gz','','.tbi')
    params:
        preset = lambda wildcards: get_GL_config('Unrevised'),#wildcards.preset),
        mem = lambda wildcards,threads,resources: threads*resources['mem_mb_per_cpu']/1024
    threads: config.get('resources',{}).get('merge',{}).get('threads',12),
    resources:
        mem_mb_per_cpu = config.get('resources',{}).get('merge',{}).get('mem_mb',6000),
        runtime = config.get('resources',{}).get('merge',{}).get('walltime','4h'),
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
