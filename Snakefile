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

regions = list(map(str,range(1,30))) + ['X','Y','MT','unplaced']

#regions = [line.strip().split()[0] for line in open('all_regions.bed')]

samples = pl.read_csv(config['alignment_metadata']).get_column('sample')
#cohort_samples = config['samples'] if 'glob' not in config['samples'] else glob_wildcards(config["bam_path"] + config["samples"]["glob"]).sample

include: 'snakepit/deepvariant.smk'
include: 'snakepit/imputation.smk'
include: 'snakepit/mendelian.smk'
#include: 'snakepit/merfin.smk'

ruleorder: deepvariant_postprocess > bcftools_scatter
ruleorder: remove_tigs > GLnexus_merge_families

def get_files():
    targets = []
    for region in config.get('regions',['all']):
        targets.append(f"{config.get('run_name','DeepVariant')}/{region}.Unrevised.vcf.gz")
    
    # handle per chromosome, imputed/unfiltered/etc
    # raw DV call
    # post-DV filter
    # raw impute
    # post-impute filter
    return targets

rule all:
    input:
        get_files()
        #expand('{name}/{region}.{preset}.vcf.gz',name=config.get('run_name','DeepVariant'),region=regions,preset=config['GL_config']),
        #expand('{name}/{region}.{preset}.vcf.gz.tbi',name=config.get('run_name','DeepVariant'),region=regions,preset=config['GL_config']),
        #expand('variants/contigs/{region}.beagle4.vcf.gz',region=regions),
        #'merfin_isec.upset',
        #'merfin_filtering.csv',
        #expand('merfin/{sample}_merfin.intersections',sample=config['samples']),
        #expand('merfin/{sample}.{vcf}.merfin.{fitted}.filter.vcf',sample=config['samples'],vcf=config['vcfs'],fitted=config['modes'])