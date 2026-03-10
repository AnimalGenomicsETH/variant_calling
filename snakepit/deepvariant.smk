def get_resource(step, key, default):
    return config.get('resources', {}).get(step, {}).get(key, default)

def shard_expand(pattern, **kwargs):
    return expand(
        pattern,
        sharding=f'{config["shards"]:05}',
        N=[f'{i:05}' for i in range(config['shards'])],
        allow_missing=True,
        **kwargs
    )

def tfrecord_sharded_path(path):
    return str(Path(path).with_suffix('').with_suffix(f'.tfrecord@{config["shards"]}.gz'))

def get_checkpoint(model):
    path = f'/opt/models/'
    match model:
        case 'PACBIO' | 'MASSEQ' | 'WGS' | 'ONT_R104' | 'RNASEQ':
            return path + model.lower()
        case 'hybrid':
            return path + 'hybrid_pacbio_illumina'
        case _:
            return config['custom_model']

def get_regions(wildcards):
    if wildcards.region == 'all':
        return ''
    else:
        return f'--regions "{' '.join(regions[wildcards.region])}"'

## implictly assume indexed bam/cram files, snakemake won't complain but jobs may fail.
def get_sample_bam_path(wildcards):
    bam = alignment_metadata.filter(pl.col('sample ID')==wildcards.sample).get_column('bam path').to_list()[0]
    match Path(bam).suffix:
        case '.bam':
            if Path(f"{bam}.bai").exists():
                return bam, f"{bam}.bai"
            elif Path(f"{bam}.csi").exists():
                return bam, f"{bam}.csi"
        case '.cram':
            if Path(f"{bam}.crai").exists():
                return bam, f"{bam}.crai"
    raise Exception(f"BAM file ({bam}) is not indexed.")

def is_long_read_model():
    return config.get('model','WGS') in ['PACBIO','ONT_R104']

rule deepvariant_make_examples:
    input:
        reference = multiext(config['reference'],'','.fai'),
        bam =  get_sample_bam_path
    output:
        examples = temp('{run}/deepvariant/intermediate_results_{sample}_{region}/make_examples.tfrecord-{N}-of-{sharding}.gz'),
        small = temp('{run}/deepvariant/intermediate_results_{sample}_{region}/make_examples_call_variant_outputs.tfrecord-{N}-of-{sharding}.gz') if config.get('small_model',False) else [],
        phasing = temp('{run}/deepvariant/intermediate_results_{sample}_{region}/read-phasing_debug-{N}-of-{sharding}.tsv') if is_long_read_model() else [],
        gvcf = temp('{run}/deepvariant/intermediate_results_{sample}_{region}/gvcf.tfrecord-{N}-of-{sharding}.gz'),
        json = temp('{run}/deepvariant/intermediate_results_{sample}_{region}/make_examples.tfrecord-{N}-of-{sharding}.gz.example_info.json')
    params:
        examples = lambda wildcards, output: tfrecord_sharded_path(output['examples']),
        gvcf = lambda wildcards, output: tfrecord_sharded_path(output['gvcf']),
        phasing = lambda wildcards, output: f"--output_phase_info --output_local_read_phasing {Path(output['phasing']).parent / f'read-phasing_debug@{config["shards"]}.tsv'}" if is_long_read_model() else '',
        small = f'--{"" if config.get("small_model") else "no"}call_small_model_examples',
        model = get_checkpoint(config['model']),
        regions = get_regions
    threads: 1
    resources:
        mem_mb_per_cpu = get_resource('make_examples', 'mem_mb', 6000),
        runtime = get_resource('make_examples', 'runtime', '4h'),
    container: config.get('deepvariant_sif')
    shell:
        '''
/opt/deepvariant/bin/make_examples \
--mode calling \
--ref {input.reference[0]} \
--include_med_dp \
--checkpoint {params.model} \
--checkpoint_json {params.model}/model.example_info.json \
--reads {input.bam[0]} \
--examples {params.examples} \
--gvcf {params.gvcf} \
--sample_name {wildcards.sample} \
{params.regions} \
{params.phasing} \
{params.small} \
--task {wildcards.N}
        '''

rule deepvariant_call_variants:
    input:
        examples = shard_expand(rules.deepvariant_make_examples.output['examples']),
        json = shard_expand(rules.deepvariant_make_examples.output['json'])
    output:
        temp('{run}/deepvariant/intermediate_results_{sample}_{region}/call_variants_output-00000-of-00001.tfrecord.gz')
    params:
        examples = lambda wildcards, input: tfrecord_sharded_path(input['examples'][0]),
        model = get_checkpoint(config['model'])
    threads: get_resource('call_variants', 'threads', 4)
    resources:
        mem_mb_per_cpu = get_resource('call_variants', 'mem_mb', 1500),
        runtime = get_resource('call_variants', 'runtime', '4h'),
    container: config.get('deepvariant_sif')
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
        gvcf = shard_expand(rules.deepvariant_make_examples.output['gvcf']),
        small = shard_expand(rules.deepvariant_make_examples.output['small']),
        phasing = shard_expand(rules.deepvariant_make_examples.output['phasing']),
    output:
        vcf = multiext('{run}/deepvariant/{sample}.{region}.vcf.gz','','.tbi'),
        gvcf = multiext('{run}/deepvariant/{sample}.{region}.g.vcf.gz','','.tbi')
    params:
        variants = lambda wildcards,input: Path(input.variants).with_name('call_variants_output@1.tfrecord.gz'),
        gvcf = lambda wildcards, input: tfrecord_sharded_path(input.gvcf[0]),
        handle_sex_chromosomes = '' if 'PAR_regions' not in config else f'--haploid_contigs "X,Y" --par_regions_bed "{config["PAR_regions"]}"',
        small_model = lambda wildcards, input: f'--small_model_cvo_records {tfrecord_sharded_path(input.small[0])}' if config.get('small_model', False) else '',
        phasing = lambda wildcards, input: f"--phased_reads_input_path {Path(input['phasing']).parent / f'read-phasing_debug@{config["shards"]}.tsv'}" if is_long_read_model() else '',
        model = get_checkpoint(config['model']),
    threads: get_resource('postprocess', 'threads', 2)
    resources:
        mem_mb_per_cpu = get_resource('postprocess', 'mem_mb', 5000),
        runtime = get_resource('postprocess', 'runtime', '4h'),
    container: config.get('deepvariant_sif')
    shell:
        '''
/opt/deepvariant/bin/postprocess_variants \
--ref {input.reference[0]} \
--cpus {threads} \
--infile {params.variants} \
--outfile {output.vcf[0]} \
--checkpoint_json {params.model}/model.example_info.json \
--gvcf_outfile {output.gvcf[0]} \
--nonvariant_site_tfrecord_path {params.gvcf} \
--novcf_stats_report \
{params.small_model} \
{params.phasing} \
{params.handle_sex_chromosomes}
        '''

rule bcftools_scatter:
    input:
        gvcf = expand(rules.deepvariant_postprocess.output['gvcf'],region='all',allow_missing=True),
        regions = config.get('regions',f"{config['reference']}.fai"),
        fai = config['reference'] + ".fai"
    output:
        gvcf = temp(expand('{run}/deepvariant/{sample}.{region}.g.vcf.gz',region=regions,allow_missing=True)),
        tbi = temp(expand('{run}/deepvariant/{sample}.{region}.g.vcf.gz.tbi',region=regions,allow_missing=True))
    wildcard_constraints:
        region = "(?!all)"
    params:
        region_cols = lambda wildcards, input: f"<(cut -f {2 if 'regions' in config else 1} {input.regions})",
        _dir = lambda wildcards, output: Path(output['gvcf'][0]).parent
    threads: 2
    resources:
        mem_mb_per_cpu = 2500,
        runtime = '1h'
    shell:
        '''
bcftools +scatter {input.gvcf[0]} -o {params._dir} -W=tbi -Oz --threads {threads} -S {params.region_cols} --no-version --prefix {wildcards.sample}.
for R in {params._dir}/{wildcards.sample}.*
do
  mv "${{R}}" "${{R//.vcf.gz/.g.vcf.gz}}"
done
        '''

# see https://github.com/dnanexus-rnd/GLnexus/pull/310
def get_GL_config(preset):
    match preset:
        case 'Unrevised':
            return f"{workflow.basedir}/GLnexus_Unrevised.yml"
        case 'DeepVariantWGS' | 'DeepVariant_unfiltered' | 'DeepVariantWES_MED_DP':
            return preset
        case _:
            return config['GL_config'][preset]

rule GLnexus_merge:
    input:
        gvcfs = expand('{run}/deepvariant/{sample}.{region}.g.vcf.gz',sample=samples,allow_missing=True),
        tbi = expand('{run}/deepvariant/{sample}.{region}.g.vcf.gz.tbi',sample=samples,allow_missing=True)
    output:
        bcf = temp('{run}/{region}.Unrevised.bcf')
    params:
        preset = lambda wildcards: get_GL_config('Unrevised'),
        mem = lambda wildcards,threads,resources: threads*resources['mem_mb_per_cpu']/1024
    threads: get_resource('merge', 'threads', 4),
    resources:
        mem_mb_per_cpu = get_resource('merge', 'mem_mb', 6000),
        runtime = get_resource('merge', 'walltime', '4h')
    container: config.get('glnexus_sif')
    shell:
        '''
/usr/local/bin/glnexus_cli \
--dir $TMPDIR/GLnexus.DB \
--config {params.preset} \
--threads {threads} \
--mem-gbytes {params.mem} \
{input.gvcfs} > {output.bcf}
        '''

rule bcftools_view:
    input:
        bcf = rules.GLnexus_merge.output['bcf']
    output:
        vcf = multiext('{run}/{region}.Unrevised.vcf.gz','','.tbi')
    threads: 4
    resources:
        mem_mb_per_cpu = 2500,
        runtime = '1h'
    shell:
        '''
bcftools view \
--threads {threads} \
--write-index=tbi \
--output {output.vcf[0]} \
{input.bcf}
        '''
