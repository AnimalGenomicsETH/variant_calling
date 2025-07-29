from pathlib import Path
import polars as pl

wildcard_constraints:
    region = r'[\w\.]+',
    sample = r'[\w-]+',
    run = r'\w+',
    filtration = r'Unrevised_filtered|beagle4',
    N = r'\d+',

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
alignment_metadata = pl.read_csv(config['alignment_metadata'],comment_prefix='#')
samples = alignment_metadata.get_column('sample ID')

include: 'snakepit/deepvariant.smk'
if config.get('split_then_merge',False):
    ruleorder: bcftools_scatter > deepvariant_postprocess
else:
    ruleorder: deepvariant_postprocess > bcftools_scatter

include: 'snakepit/imputation.smk'
include: 'snakepit/mendelian.smk'
include: 'snakepit/SV_calling.smk'
#include: 'snakepit/merfin.smk'

def get_files():
    targets = []

    if 'SV' in config.get('targets',[]):
        for region in regions:
            targets.append(f'SVs/sawfish/cohort_{region}')

    if 'small' in config.get('targets',[]):
        postprocess_steps = {} #TODO: this is not good
        for i,step in enumerate(config.get('variant_postprocess',[])):
            if i == 0:
                postprocess_steps[step] = 'Unrevised'
            else:
                postprocess_steps[step] = list(postprocess_steps.keys())[-1]

        if 'regions' in config:
            for region in regions:
                targets.append(f"{config.get('run_name','DeepVariant')}/{region}.{config.get('imputed','Unrevised')}.vcf.gz")
        else:
            if config.get('split_then_merge',False):
                for region in regions:
                    targets.append(f"{config.get('run_name','DeepVariant')}/{region}.{config.get('imputed','Unrevised')}.vcf.gz")
            else:
                targets.append(f"{config.get('run_name','DeepVariant')}/all.{config.get('imputed','Unrevised')}.vcf.gz")

    return targets

rule all:
    input:
        get_files()
