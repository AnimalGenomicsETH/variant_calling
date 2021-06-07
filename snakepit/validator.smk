from pathlib import PurePath
from itertools import product

wildcard_constraints:
    norm = r'normed|raw',
    collapse = r'none|all',
    cols = r'pos_only|bases',
    var = r'snp|indel'
    
localrules: format_record_tsv, summarise_matches, variant_density, generate_karyotype

def capture_logic():
    targets = []
    
    
rule all:
    input:
        #expand('intersection_{norm}_{collapse}.{cols}.{var}.summary.stats',norm=('normed','raw'),collapse=('none','all'),cols=('pos_only','bases'),var='snp'),
        #expand('intersection_{norm}_{collapse}.{cols}.{var}.summary.stats',norm=('normed','raw'),collapse=('none','all'),cols='pos_only',var='indel'),
        expand('intersection_{norm}_{collapse}.{cols}.{var}.summary.stats',norm='normed',collapse='none',cols=('pos_only','bases'),var='snp'),
        expand('intersection_{norm}_{collapse}.{cols}.{var}.summary.stats',norm='normed',collapse='none',cols='pos_only',var='indel'),
        expand('{mode}.{var}.{norm}.{collapse}.ideogram.png',mode=('shared','missing','exclusive'),var='snp',norm='normed',collapse='none'),
        expand('intersection_{norm}_{collapse}/{N}.qc',N=('DV','GATK','shared'),norm='normed',collapse='none'),
        'array_concordance/summary.txt'

rule bcftools_norm:
    input:
        vcf = lambda wildcards: config['vcf'][wildcards.caller],
        ref = config['ref']
    output:
        '{caller}.normed.vcf.gz'
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
        DV = lambda wildcards: 'DV.normed.vcf.gz' if wildcards.norm == 'normed' else config['vcf']['DV'],
        GATK = lambda wildcards: 'GATK.normed.vcf.gz' if wildcards.norm == 'normed' else config['vcf']['GATK']
    output:
        expand('intersection_{{norm}}_{{collapse}}/{N}.vcf',N=('DV','GATK','shared'))
    params:
        lambda wildcards, output: PurePath(output[0]).parent
    threads: 4
    resources:
        mem_mb = 2000,
        walltime = '45'
    shell:
        '''
        bcftools isec --threads {threads} -c {wildcards.collapse} -p {params} {input.DV} {input.GATK}
        mv {params}/0000.vcf {params}/DV.vcf
        mv {params}/0001.vcf {params}/GATK.vcf
        mv {params}/0002.vcf {params}/shared.vcf
        '''

rule average_quality:
    input:
        'intersection_{norm}_{collapse}/{N}.vcf'
    output:
        'intersection_{norm}_{collapse}/{N}.qc'
    shell:
        '''
        awk '!/#/ {{c+=$6;j+=1}} END {{print c/j}}' {input} > {output}
        '''
        
rule format_record_tsv:
    input:
        'intersection_{norm}_{collapse}/{callset}.vcf'
    output:
        'intersection_{norm}_{collapse}/{callset}.{cols}.{var}.tsv'
    params:
        target_cols = lambda wildcards: '$1"\t"$2"\t"$5' if wildcards.cols == 'bases' else '$1"\t"$2',
        v_type = lambda wildcards: '-m2 -M2 -v snps' if wildcards.var == 'snp' else '-m2 -M2 -v indels'
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
        variants = 'intersection_{norm}_{collapse}/{callset}.{cols}.{var}.tsv',
        truth = 'truth.{cols}.tsv' #lambda wildcards: config['truth'][wildcards.cols]
    output:
        'intersection_{norm}_{collapse}/{callset}.{cols}.{var}.matches.count'
    threads: 8 #lambda wildcards, input: int(max(input.size_mb,1))
    resources:
        mem_mb = 5000,
        disk_scratch = 10,
        walltime = '15'
    shell:
        '''
        parallel --pipepart -a {input.variants} --jobs {threads} --block 1M LC_ALL=C fgrep -F -c -f {input.truth} | awk '{{c+=$1}}END{{print c}}' > {output}
        '''

rule truth_only_positions:
    input:
        variants = expand('intersection_{{norm}}_{{collapse}}/{callset}.{{cols}}.snp.tsv',callset=('shared','DV','GATK')),
        truth = 'truth.{cols}.tsv' #lambda wildcards: config['truth'][wildcards.cols]
    output:
        missing = 'intersection_{norm}_{collapse}/missing.{cols}.snp.tsv',
        query = temp('intersection_{norm}_{collapse}/missing_temp.{cols}.var'),
        truth = temp('intersection_{norm}_{collapse}/truth_temp.{cols}.tsv')
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
        tsvs = expand('intersection_{{norm}}_{{collapse}}/{callset}.{{cols}}.{{var}}.tsv',callset=('DV','GATK','shared')),
        counts = expand('intersection_{{norm}}_{{collapse}}/{callset}.{{cols}}.{{var}}.matches.count',callset=('DV','GATK','shared'))
    output:
        'intersection_{norm}_{collapse}.{cols}.{var}.summary.stats'
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
        lambda wildcards: 'intersection_{norm}_{collapse}/{mode}.pos_only.{var}.tsv' if wildcards.mode != 'exclusive' else ('intersection_{norm}_{collapse}/DV.pos_only.{var}.tsv','intersection_{norm}_{collapse}/GATK.pos_only.{var}.tsv')
    output:
        'intersection_{norm}_{collapse}/{mode}.{var}.ideogram.tsv'
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
        features = 'intersection_{norm}_{collapse}/{mode}.{var}.ideogram.tsv',
    output:
        '{mode}.{var}.{norm}.{collapse}.ideogram.png'
    params:
        lambda wildcards, output: PurePath(output[0]).with_suffix('')
    shell:
        'Rscript --vanilla {workflow.basedir}/../scripts/ideogram_plotter.R {input.karotype} {input.features} --output {params}'

rule get_snps:
    input:
        '{caller}.normed.vcf.gz'
    output:
        '{caller}.normed.snp.vcf'
    shell:
        'bcftools view -m2 -M2 -v snps -o {output} {input}'

rule picard_concordance:
    input:
        query = lambda wildcards: f'{wildcards.caller}.normed.snp.vcf',
        truth = lambda wildcards: config['chips'][wildcards.chip]['vcf']
    output:
        concordance = 'array_concordance/{sample}.{chip}.{caller}.conc',
        metrics = 'array_concordance/{sample}.{chip}.{caller}.metrics'
    params:
        ID = lambda wildcards: config['chips'][wildcards.chip][wildcards.sample],
        gc_out = lambda wildcards,output: PurePath(output[0]).parent / f'{wildcards.sample}_{wildcards.caller}'
    threads: 1
    resources:
        mem_mb = 25000,
        walltime = '60'
    envmodules:
        'gcc/8.2.0',
        'picard/2.25.4'
    shell:
        '''
        picard GenotypeConcordance CALL_VCF={input.query} TRUTH_VCF={input.truth} CALL_SAMPLE={wildcards.sample} TRUTH_SAMPLE={params.ID} O={params.gc_out}
        tail -n 4 {params.gc_out}.genotype_concordance_summary_metrics | awk '$1~/SNP/ {{print $2"\t{wildcards.caller}\t"$10"\t"$12"\t"2*$4*$5/($4+$5)"\t"2*$7*$8/($7+$8)"\t"2*$10*$11/($10+$11)"\t"$13"\t"$14}}' > {output.concordance}
        awk '$1~/SNP/ {{print $3"\t{wildcards.chip}\t{wildcards.caller}\t"$4"\t"$5"\t"$6}}' {params.gc_out}.genotype_concordance_detail_metrics > {output.metrics} 
        '''

def ID_concordances():
    targets = []
    for chip,keys in config['chips'].items():
        for sample in list(keys.keys())[1:]:
            targets.append(f'{sample}.{chip}')
    return targets

rule aggreate_concordance:
    input:
        concordance = expand('array_concordance/{ID}.{caller}.conc',ID=ID_concordances(),caller=('DV','GATK')),
        metrics = expand('array_concordance/{ID}.{caller}.metrics',ID=ID_concordances(),caller=('DV','GATK'))
    output:
        summary = 'array_concordance/summary.txt',
        df = 'array_concordance/summary.df'
    shell:
        '''
        echo -e "sample\tcaller\tSensitivity\tSpecificity\tHET F1\tHOMV F1\tVAR F1\tGC\tnR-GC" > {output.summary}
        sort -k 13,13nr -k 11,11nr {input.concordance} >> {output.summary}
        echo -e "sample\tchip\tcaller\ttruth\tcall\tcount" > {output.df}
        cat {input.metrics} >> {output.df}
        '''

rule plot_concordance:
    input:
        'array_concordance/summary.df'
    output:
        'array_concordance/summary.png'
    run:
        import pandas as pd
        import seaborn as sns
        import matplotlib.pyplot as plt

        df = pd.read_csv(input,delimiter='\t')
        df_no_missing = df[(df['truth']!='MISSING')&(df['call']!='MISSING')]
        g = sns.catplot(data=df,x='call',y='count',col='truth',col_wrap=2,hue='caller',sharey=False)
        for ax in g.axes.flatten():
            ax.set_yscale('log')
        plt.savefig(output[0])
        #df2['conc']=df2['truth']==df2['call']
        #df2.groupby(['caller','conc']).sum()
