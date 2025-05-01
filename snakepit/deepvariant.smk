#NOTE: may need to be updated if deepvariant changes its internal parameters.
#e.g., running Apptainer> run_deepvariant  --model_type ONT_R104 --ref GL_DV_raw.yml --reads README.md --output_vcf test --num_shards 10 --intermediate_results_dir inter --dry_run --disable_small_model=true

#use_multiallelic_model?
def make_custom_example_arguments(model):
    match model:
        case 'WGS':
            if config.get('small_model',False):
                return '--checkpoint "/opt/models/wgs" --call_small_model_examples --small_model_indel_gq_threshold "30" --small_model_snp_gq_threshold "25" --small_model_vaf_context_window_size "51" --trained_small_model_path "/opt/smallmodels/wgs"'
            else:
                return '--checkpoint "/opt/models/wgs"'
        case 'PACBIO':
            if config.get('small_model',False):
                return '--checkpoint "/opt/models/pacbio" --alt_aligned_pileup "diff_channels" --call_small_model_examples --max_reads_per_partition "600" --min_mapping_quality "1" --parse_sam_aux_fields --partition_size "25000" --phase_reads --pileup_image_width "147" --norealign_reads --small_model_indel_gq_threshold "30" --small_model_snp_gq_threshold "25" --small_model_vaf_context_window_size "51" --sort_by_haplotypes --track_ref_reads --trained_small_model_path "/opt/smallmodels/pacbio" --trim_reads_for_pileup --vsc_min_fraction_indels "0.12"'
            else:
                return '--checkpoint "/opt/models/pacbio" --alt_aligned_pileup "diff_channels" --max_reads_per_partition "600" --min_mapping_quality "1" --parse_sam_aux_fields --partition_size "25000" --phase_reads --pileup_image_width "147" --norealign_reads --sort_by_haplotypes --track_ref_reads --trim_reads_for_pileup --vsc_min_fraction_indels "0.12"'
        case 'MASSEQ': #no small model implemented
            return '--checkpoint "/opt/models/masseq" --alt_aligned_pileup "diff_channels" --max_reads_for_dynamic_bases_per_region "1500" --max_reads_per_partition "0" --min_mapping_quality "1" --parse_sam_aux_fields --partition_size "25000" --phase_reads --pileup_image_width "199" --norealign_reads --sort_by_haplotypes --track_ref_reads --trim_reads_for_pileup --vsc_min_fraction_indels "0.12"'
        case 'ONT_R104':
            if config.get('small_model',False):
                return '--checkpoint "/opt/models/ont_r104" --alt_aligned_pileup "diff_channels" --call_small_model_examples --max_reads_per_partition "600" --min_mapping_quality "5" --parse_sam_aux_fields --partition_size "25000" --phase_reads --pileup_image_width "99" --norea    lign_reads --small_model_indel_gq_threshold "25" --small_model_snp_gq_threshold "20" --small_model_vaf_context_window_size "51" --sort_by_haplotypes --track_ref_reads --trained_small_model_path "/opt/smallmodels/ont_r104" --trim_reads_for_pileup --vsc_min_fraction_indels "0.12"     --vsc_min_fraction_snps "0.08"'
            else:
                return '--checkpoint "/opt/models/ont_r104" --alt_aligned_pileup "diff_channels" --max_reads_per_partition "600" --min_mapping_quality "5" --parse_sam_aux_fields --partition_size "25000" --phase_reads --pileup_image_width "99" --norealign_reads --sort_by_haplotypes --track_ref_reads --trim_reads_for_pileup --vsc_min_fraction_indels "0.12" --vsc_min_fraction_snps "0.08"'
        case 'RNA' | 'WES' | _: #only valid for v1.4-v1.5
            return '--channels "" --split_skip_reads'

def get_checkpoint(model):
    path = f'/opt/models/'
    match model:
        case 'PACBIO' | 'MASSEQ' | 'WGS' | 'ONT_R104' :
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

rule deepvariant_make_examples:
    input:
        reference = multiext(config['reference'],'','.fai'),
        bam =  get_sample_bam_path
    output:
        examples = temp('{run}/deepvariant/intermediate_results_{sample}_{region}/make_examples.tfrecord-{N}-of-{sharding}.gz'),
        small = temp('{run}/deepvariant/intermediate_results_{sample}_{region}/make_examples_call_variant_outputs.tfrecord-{N}-of-{sharding}.gz') if config.get('small_model',False) else [],
        gvcf = temp('{run}/deepvariant/intermediate_results_{sample}_{region}/gvcf.tfrecord-{N}-of-{sharding}.gz'),
        json = temp('{run}/deepvariant/intermediate_results_{sample}_{region}/make_examples.tfrecord-{N}-of-{sharding}.gz.example_info.json')
    params:
        examples = lambda wildcards, output: Path(output['examples']).with_suffix('').with_suffix(f'.tfrecord@{config["shards"]}.gz'),
        gvcf = lambda wildcards, output: Path(output['gvcf']).with_suffix('').with_suffix(f'.tfrecord@{config["shards"]}.gz'),
        model_args = make_custom_example_arguments(config['model']),
        regions = get_regions
    threads: 1
    resources:
        mem_mb_per_cpu = config.get('resources',{}).get('make_examples',{}).get('mem_mb',6000),
        runtime = config.get('resources',{}).get('make_examples',{}).get('runtime','4h'),
    container: config.get('deepvariant_sif')
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
        examples = expand(rules.deepvariant_make_examples.output['examples'],sharding=f'{config["shards"]:05}',N=[f'{i:05}' for i in range(config['shards'])],allow_missing=True),
        json = expand(rules.deepvariant_make_examples.output['json'],sharding=f'{config["shards"]:05}',N=[f'{i:05}' for i in range(config['shards'])],allow_missing=True)
    output:
        temp('{run}/deepvariant/intermediate_results_{sample}_{region}/call_variants_output-00000-of-00001.tfrecord.gz')
    params:
        examples = lambda wildcards, input: Path(input['examples'][0]).with_suffix('').with_suffix(f'.tfrecord@{config["shards"]}.gz'),
        model = get_checkpoint(config['model'])
    threads: config.get('resources',{}).get('call_variants',{}).get('threads',4)
    resources:
        mem_mb_per_cpu = config.get('resources',{}).get('call_variants',{}).get('mem_mb',5000),
        runtime = config.get('resources',{}).get('call_variants',{}).get('runtime','4h'),
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
        gvcf = expand(rules.deepvariant_make_examples.output['gvcf'],sharding=f'{config["shards"]:05}',N=[f'{i:05}' for i in range(config['shards'])],allow_missing=True),
        small = expand(rules.deepvariant_make_examples.output['small'],sharding=f'{config["shards"]:05}',N=[f'{i:05}' for i in range(config['shards'])],allow_missing=True)
    output:
        vcf = multiext('{run}/deepvariant/{sample}.{region}.vcf.gz','','.tbi'),
        gvcf = multiext('{run}/deepvariant/{sample}.{region}.g.vcf.gz','','.tbi')
    params:
        variants = lambda wildcards,input: Path(input.variants).with_name('call_variants_output@1.tfrecord.gz'),
        gvcf = lambda wildcards,input: Path(input.gvcf[0]).with_suffix('').with_suffix(f'.tfrecord@{config["shards"]}.gz'),
        handle_sex_chromosomes = '' if 'PAR_regions' not in config else f'--haploid_contigs "X,Y" --par_regions_bed "{config['PAR_regions']}"',
        small_model = lambda wildcards, input: f'--small_model_cvo_records {Path(input.small[0]).with_suffix("").with_suffix(f".tfrecord@{config['shards']}.gz")}' if config.get('small_model',False) else '',
        #this is a workaround to allow older versions without parallel postprocessing to not complain
        cpus = lambda wildcards, threads: '--cpus {threads}' if config.get('resources',{}).get('postprocess',{}).get('threads',1) > 1 else ''
    threads: config.get('resources',{}).get('postprocess',{}).get('threads',2)
    resources:
        mem_mb_per_cpu = config.get('resources',{}).get('postprocess',{}).get('mem_mb',5000),
        runtime = config.get('resources',{}).get('postprocess',{}).get('runtime','4h'),
    container: config.get('deepvariant_sif')
    shell:
        '''
/opt/deepvariant/bin/postprocess_variants \
--ref {input.reference[0]} \
--infile {params.variants} \
--outfile {output.vcf[0]} \
--gvcf_outfile {output.gvcf[0]} \
--nonvariant_site_tfrecord_path {params.gvcf} \
--novcf_stats_report \
{params.small_model} \
{params.handle_sex_chromosomes} \
{params.cpus}
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
        _dir = lambda wildcards, output: Path(output[0]).with_suffix('').with_suffix('').with_suffix('').with_suffix('')
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
        bcf = '{run}/{region}.Unrevised.bcf'
    params:
        preset = lambda wildcards: get_GL_config('Unrevised'),
        mem = lambda wildcards,threads,resources: threads*resources['mem_mb_per_cpu']/1024
    threads: config.get('resources',{}).get('merge',{}).get('threads',4),
    resources:
        mem_mb_per_cpu = config.get('resources',{}).get('merge',{}).get('mem_mb',6000),
        runtime = config.get('resources',{}).get('merge',{}).get('walltime','4h')
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
--output {output[0]} \
{input.bcf}
        '''
