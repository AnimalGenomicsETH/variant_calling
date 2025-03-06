from pathlib import Path, PurePath
import polars as pl
from itertools import product

wildcard_constraints:
    region = r'\w+',
    run = r'\w+',
    N = r'\d',

if 'binding_paths' in config:
    for path in config['binding_paths']:
        workflow._singularity_args += f' -B {path}'

genome_regions = config.get('bed')
regions = list(map(str,range(1,30))) + ['X','Y','MT','unplaced']

regions = [line.strip().split()[0] for line in open('all_regions.bed')]

cohort_samples = config['samples'] if 'glob' not in config['samples'] else glob_wildcards(config["bam_path"] + config["samples"]["glob"]).sample

include: 'snakemake/deepvariant.smk'
include: 'snakemake/imputation.smk'
include: 'snakemake/medelian.smk'
include: 'snakemake/merfin.smk'

ruleorder: deepvariant_postprocess > bcftools_scatter
ruleorder: remove_tigs > GLnexus_merge_families

def get_files():
    if 'regions' not in config:
        config['regions'] == 'all'
    # handle per chromosome, imputed/unfiltered/etc
    # raw DV call
    # post-DV filter
    # raw impute
    # post-impute filter
    return []

rule all:
    input:
        expand('{name}/{region}.{preset}.vcf.gz',name=config.get('run_name','DeepVariant'),region=regions,preset=config['GL_config']),
        expand('{name}/{region}.{preset}.vcf.gz.tbi',name=config.get('run_name','DeepVariant'),region=regions,preset=config['GL_config']),
        expand('variants/contigs/{region}.beagle4.vcf.gz',region=regions),
        'merfin_isec.upset',
        'merfin_filtering.csv',
        expand('merfin/{sample}_merfin.intersections',sample=config['samples']),
        expand('merfin/{sample}.{vcf}.merfin.{fitted}.filter.vcf',sample=config['samples'],vcf=config['vcfs'],fitted=config['modes'])