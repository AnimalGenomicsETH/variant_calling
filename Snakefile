from pathlib import Path, PurePath
import polars as pl
#from itertools import product

wildcard_constraints:
    region = r'[\w\.]+',
    sample = r'\w+',
    run = r'\w+',
    filtration = r'Unrevised_filtered|beagle4',
    N = r'\d+',

if 'binding_paths' in config:
    for path in config['binding_paths']:
        workflow._singularity_args += f' -B {path}'

from collections import defaultdict
def build_regions():
    region_map = defaultdict(list)

    for line in open(config.get('regions',f"{config['reference']}.fai")):
        parts = line.rstrip().split()
        if len(parts) == 2:
            region_map[parts[1]].append(parts[0])
        else: #if single entry or 4-column is assumed FAI format
           region_map[parts[0]].append(parts[0])

    return dict(region_map)

regions = build_regions()
alignment_metadata = pl.read_csv(config['alignment_metadata'])
samples = alignment_metadata.get_column('sample ID')

include: 'snakepit/deepvariant.smk'
ruleorder: deepvariant_postprocess > bcftools_scatter
include: 'snakepit/imputation.smk'
include: 'snakepit/mendelian.smk'
#include: 'snakepit/merfin.smk'

def get_files():
    targets = []

    postprocess_steps = {}
    for i,step in enumerate(config['variant_postprocess']):
        if i == 0:
            postprocess_steps[step] = 'Unrevised'
        else:
            postprocess_steps[step] = list(postprocess_steps.keys())[-1]

    for region in regions:
        targets.append(f"{config.get('run_name','DeepVariant')}/{region}.{config.get('imputed','Unrevised')}.vcf.gz")
    
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