from pathlib import PurePath
from itertools import product, combinations

wildcard_constraints:
    norm = r'normed|raw',
    collapse = r'none|all',
    cols = r'pos_only|bases',
    var = r'snp|indel'
    
#localrules: format_record_tsv, summarise_matches, variant_density, generate_karyotype

class Default(dict):
    def __missing__(self, key):
        return '{'+key+'}'

def get_dir(base='work',ext='', **kwargs):
    if base == 'concordance':
        base_dir = 'array_concordance'
    elif base == 'intersection':
        base_dir = 'intersection_{norm}_{collapse}_{caller1}_{caller2}'
    elif base == 'main':
        base_dir = ''
    else:
        raise Exception('Base not found')
    return str(Path(base_dir.format_map(Default(kwargs))) / ext.format_map(Default(kwargs)))

def capture_logic():
    targets = []

    for vcfs in combinations(config['vcfs'],2):
        for callset in (*vcfs,'shared'):
            targets.append(get_dir('intersection',f'{callset}.bases.snp.matches.count',norm='normed',collapse='all',caller1=vcfs[0],caller2=vcfs[1]))
    
    #for vcf in config['vcf']:
    #    targets.append('intersection_normed_all.DV1.DV2.pos_only.snp.summary.stats')
    
    for T in ('vcfs','imputed'):
        targets.append(get_dir('concordance',f'{T}.summary.txt'))
    return targets
    
    
rule all:
    input:
        #expand('intersection_{norm}_{collapse}.{cols}.{var}.summary.stats',norm=('normed','raw'),collapse=('none','all'),cols=('pos_only','bases'),var='snp'),
        #expand('intersection_{norm}_{collapse}.{cols}.{var}.summary.stats',norm=('normed','raw'),collapse=('none','all'),cols='pos_only',var='indel'),
        #expand('intersection_{norm}_{collapse}.{cols}.{var}.summary.stats',norm='normed',collapse='none',cols=('pos_only','bases'),var='snp'),
        #expand('intersection_{norm}_{collapse}.{cols}.{var}.summary.stats',norm='normed',collapse='none',cols='pos_only',var='indel'),
        #expand('{mode}.{var}.{norm}.{collapse}.ideogram.png',mode=('shared','missing','exclusive'),var='snp',norm='normed',collapse='none'),
        #expand('intersection_{norm}_{collapse}/{N}.qc',N=('DV','GATK','shared'),norm='normed',collapse='none'),
        #'array_concordance/summary.txt'
        capture_logic()

rule bcftools_norm:
    input:
        vcf = lambda wildcards: config[wildcards.imputed][wildcards.caller],
        ref = config['ref']
    output:
        get_dir('main','{caller}.{imputed}.normed.vcf.gz')
    threads: 4
    resources:
        mem_mb = 2000,
        walltime = '30'
    shell:
        '''
        bcftools norm --threads {threads} -m- -f {input.ref} -Oz -o {output} {input.vcf}
        tabix -p vcf {output}
        '''

rule bcftools_isec:
    input:
        caller1 = lambda wildcards: get_dir('main','{caller1}.normed.vcf.gz') if wildcards.norm == 'normed' else config['vcfs'][wildcards.caller1],
        caller2 = lambda wildcards: get_dir('main','{caller2}.normed.vcf.gz') if wildcards.norm == 'normed' else config['vcfs'][wildcards.caller2]
    output:
        (get_dir('intersection','{callset}.vcf',callset=n) for n in ('{caller1}','{caller2}','shared'))
    params:
        lambda wildcards, output: PurePath(output[0]).parent
    threads: 4
    resources:
        mem_mb = 2000,
        walltime = '45'
    shell:
        '''
        bcftools isec --threads {threads} -c {wildcards.collapse} -p {params} {input.caller1} {input.caller2}
        mv {params}/0000.vcf {params}/{wildcards.caller1}.vcf
        mv {params}/0001.vcf {params}/{wildcards.caller2}.vcf
        mv {params}/0002.vcf {params}/shared.vcf
        '''

rule average_quality:
    input:
        get_dir('intersection','{callset}.vcf')
    output:
        get_dir('intersection','{callset}.qv')
    shell:
        '''
        awk '!/#/ {{c+=$6;j+=1}} END {{print c/j}}' {input} > {output}
        '''
        
rule format_record_tsv:
    input:
        get_dir('intersection','{callset}.vcf')
    output:
        get_dir('intersection','{callset}.{cols}.{var}.tsv')
    params:
        target_cols = lambda wildcards: '$1"\t"$2"\t"$4"\t"$5' if wildcards.cols == 'bases' else '$1"\t"$2',
        v_type = lambda wildcards: '-m2 -M2 -v snps' if wildcards.var == 'snp' else '-m2 -M2 -v indels'
    resources:
        walltime = '30'
    shell:
        '''
        bcftools view {params.v_type} {input} | awk '$1!~/#/ {{print {params.target_cols}}}' > {output}
        '''

rule format_truth_tsv:
    input:
        config['truth']
    output:
        'truth.{cols}.tsv'
    params:
        target_cols = lambda wildcards: '$1"\t"$2"\t"$4"\t"$5' if wildcards.cols == 'bases' else '$1"\t"$2'
    shell:
        '''
        zcat {input} | awk '$1!~/#/ {{print {params.target_cols}}}' > {output}
        '''

rule count_matches:
    input:
        variants = get_dir('intersection','{callset}.{cols}.{var}.tsv'),
        truth = 'truth.{cols}.tsv' #lambda wildcards: config['truth'][wildcards.cols]
    output:
        get_dir('intersection','{callset}.{cols}.{var}.matches.count')
    threads: 8 #lambda wildcards, input: int(max(input.size_mb,1))
    resources:
        mem_mb = 5000,
        disk_scratch = 10,
        walltime = '30'
    shell:
        '''
        parallel --pipepart -a {input.variants} --jobs {threads} --block 1M LC_ALL=C fgrep -F -c -f {input.truth} | awk '{{c+=$1}}END{{print c}}' > {output}
        '''

rule truth_only_positions:
    input:
        variants = lambda wildcards: (get_dir('intersection','{callset}.{cols}.snp.tsv',callset=C) for C in ('shared',wildcards.caller1,wildcards.caller2)),
        truth = 'truth.{cols}.tsv' #lambda wildcards: config['truth'][wildcards.cols]
    output:
        missing = 'intersection_{norm}_{collapse}_{caller1}_{caller2}/missing.{cols}.snp.tsv',
        query = temp('intersection_{norm}_{collapse}_{caller1}_{caller2}/missing_temp.{cols}.var'),
        truth = temp('intersection_{norm}_{collapse}_{caller1}_{caller2}/truth_temp.{cols}.tsv')
    threads: 8
    resources:
        mem_mb = 3000,
        walltime = '60'
    shell:
        '''
        for i in {{1..29}}; do
          awk -v i=$i '$1 == i' {input.variants[0]} > {output.query}
          awk -v i=$i '$1 == i' {input.variants[1]} >> {output.query}
          awk -v i=$i '$1 == i' {input.variants[2]} >> {output.query}
          awk -v i=$i '$1 == i' {input.truth} > {output.truth}
          parallel --pipepart -a {output.truth} --jobs {threads} --block 100k LC_ALL=C fgrep -F -v -f {output.query} >> {output.missing}
        done
        '''

rule summarise_matches:
    input:
        tsvs = lambda wildcards: expand('intersection_{{norm}}_{{collapse}}_{{caller1}_{{caller2}}/{callset}.{{cols}}.{{var}}.tsv',callset=(wildcards.caller1,wildcards.caller2,'shared')),
        counts = lambda wildcards: expand('intersection_{{norm}}_{{collapse}}_{{caller1}_{{caller2}}/{callset}.{{cols}}.{{var}}.matches.count',callset=(wildcards.caller1,wildcards.caller2,'shared'))
    output:
        'intersection_{norm}_{collapse}.{caller1}.{caller2}.{cols}.{var}.summary.stats'
    shell:
        '''
        echo $(wc -l {input.tsvs[0]}) $(cat {input.counts[0]}) | awk '{{print "DV "$3/$1}}' >> {output}
        echo $(wc -l {input.tsvs[1]}) $(cat {input.counts[1]}) | awk '{{print "GATK "$3/$1}}' >> {output}
        echo $(wc -l {input.tsvs[2]}) $(cat {input.counts[2]}) | awk '{{print "shared "$3/$1}}' >> {output}
        '''

from collections import defaultdict

def get_density_counts(f_path,window):
    v_counts = defaultdict(int)
    with open(f_path,'r') as fin:
        for line in fin:
            chr, pos = line.rstrip().split()
            v_counts[(int(chr),int(pos)//window)] += 1
    return v_counts

rule variant_density:
    input:
        lambda wildcards: get_dir('intersection','{mode}.pos_only.{var}.tsv') if wildcards.mode != 'exclusive' else ('intersection_{norm}_{collapse}/DV.pos_only.{var}.tsv','intersection_{norm}_{collapse}/GATK.pos_only.{var}.tsv')
    output:
        get_dir('intersection','{mode}.{var}.ideogram.tsv')
    run:
        window = 100000
        with open(output[0],'w') as fout:
            if wildcards.mode != 'exclusive':
                fout.write('\t'.join(('Chr','Start','End','Value','Color')) + '\n')
                counts = get_density_counts(input[0],window)
                for (chr,pos),value in counts.items():
                    pos *= window
                    fout.write('\t'.join(map(str,(chr,pos,pos+window,value,'17A589' if wildcards.mode == 'shared' else 'A51733'))) + '\n')

            else:
                fout.write('\t'.join(('Chr','Start','End','Value_1','Color_1','Value_2','Color_2')) + '\n')
                counts_1 = get_density_counts(input[0],window)
                counts_2 = get_density_counts(input[1],window)
                for key in sorted(counts_1.keys() | counts_2.keys()):
                    chr,pos = key 
                    pos *= window
                    fout.write('\t'.join(map(str,(chr,pos,pos+window,counts_1[key],'2874A6',counts_2[key],'D4AC0D'))) + '\n')

rule generate_karyotype:
    input:
        config['ref'] + ".fai"
    output:
        'ref.karotype.tsv'
    shell:
        '''
        echo -e "Chr\tStart\tEnd" > {output}
        head -n 29 {input} | awk '{{print $1"\t1\t"$2}}' >> {output}
        '''

rule plot_density:
    input:
        karotype = 'ref.karotype.tsv', 
        features = get_dir('intersection','{mode}.{var}.ideogram.tsv'),
    output:
        '{mode}.{var}.{norm}.{collapse}.{caller1}_{caller2}.ideogram.png'
    params:
        lambda wildcards, output: PurePath(output[0]).with_suffix('')
    shell:
        'Rscript --vanilla {workflow.basedir}/../scripts/ideogram_plotter.R {input.karotype} {input.features} --output {params}'

rule get_snps:
    input:
        get_dir('main','{caller}.normed.vcf.gz')
    output:
        get_dir('main','{caller}.normed.snp.vcf')
    shell:
        'bcftools view -m2 -M2 -v snps -o {output} {input}'

#rule gtcomapre:
#vcfgtcompare.sh other cohort.29.imputed.vcf.gz ../imputed/29_beagle.vcf.gz | vcf2tsv >this_other_comparison.tsv

rule picard_concordance:
    input:
        query = get_dir('main','{caller}.{imputed}.normed.snp.vcf'),
        truth = lambda wildcards: config['chips'][wildcards.chip]['vcf']
    output:
        multiext(get_dir('concordance','{sample}.{chip}.{caller}.{imputed}'),'.conc','.metrics')
    params:
        ID = lambda wildcards: config['chips'][wildcards.chip][wildcards.sample],
        gc_out = lambda wildcards,output: PurePath(output[0]).parent / f'{wildcards.sample}_{wildcards.caller}'
    threads: 1
    resources:
        mem_mb = 25000,
        walltime = '60'
    envmodules:
        'gcc/8.2.0',
        'picard/2.25.7'
    shell:
        '''
        picard GenotypeConcordance CALL_VCF={input.query} TRUTH_VCF={input.truth} CALL_SAMPLE={wildcards.sample} TRUTH_SAMPLE={params.ID} O={params.gc_out}
        tail -n 4 {params.gc_out}.genotype_concordance_summary_metrics | awk '$1~/SNP/ {{print $2"\t{wildcards.caller}\t"$10"\t"$12"\t"2*$4*$5/($4+$5)"\t"2*$7*$8/($7+$8)"\t"2*$10*$11/($10+$11)"\t"$13"\t"$14}}' > {output[0]}
        awk '$1~/SNP/ {{print $3"\t{wildcards.chip}\t{wildcards.caller}\t"$4"\t"$5"\t"$6}}' {params.gc_out}.genotype_concordance_detail_metrics > {output[1]} 
        '''
        
def ID_concordances():
    targets = []
    for chip,keys in config['chips'].items():
        for sample in list(keys.keys())[1:]:
            targets.append(f'{sample}.{chip}')
    return targets

rule aggreate_concordance:
    input:
        #multiext(get_dir('concordance','{ID}.{caller}.{imputed}',ID=I,caller=C,imputed=wildcards.imputed) for (I,C) in product(ID_concordances(),config[wildcards.imputed]),'.conc','.metrics')
        concordance = lambda wildcards: (get_dir('concordance','{ID}.{caller}.{imputed}.conc',ID=I,caller=C,imputed=wildcards.imputed) for (I,C) in product(ID_concordances(),config[wildcards.imputed])),
        metrics = lambda wildcards: (get_dir('concordance','{ID}.{caller}.{imputed}.metrics',ID=I,caller=C,imputed=wildcards.imputed) for (I,C) in product(ID_concordances(),config[wildcards.imputed]))
    output:
        multiext(get_dir('concordance','{imputed}.summary'),'.txt','.df')
    shell:
        '''
        echo -e "sample\tcaller\tSensitivity\tSpecificity\tHET F1\tHOMV F1\tVAR F1\tGC\tnR-GC" > {output[0]}
        sort -k 13,13nr -k 11,11nr {input.concordance} >> {output[0]}
        echo -e "sample\tchip\tcaller\ttruth\tcall\tcount\timputed" > {output[1]}
        cat {input.metrics} | awk '{{print $0"\\t{wildcards.imputed}"}}' >> {output[1]}
        '''
