from pathlib import PurePath

rule all:
    input:
        str(PurePath(config['imputed']) / "cohort.autosomes.shapeit.vcf.gz")
        
rule beagle_imputation:
    input:
        str(PurePath(config['vcfs']) / "cohort.{chromosome}.vcf.gz")
    output:
        str(PurePath(config['imputed']) / "cohort.{chromosome}.beagle.vcf.gz")
    params:
        lambda wildcards, output: PurePath(output[0]).with_suffix('').with_suffix('')
    threads: 12
    resources:
        mem_mb = 3500,
        walltime = '4:00'
    shell:
        '''
        java -jar -Xss25m -Xmx40G {config[beagle]} gl={input} nthreads={threads} out={params}
        tmp_file=$(mktemp .XXXXXX)
        mv {output} $tmp_file
        bcftools reheader -f {config[reference]} -o {output} $tmp_file
        tabix -fp vcf {output}
        '''

#rule prepare_genetic_maps:
    output:
        str(PurePath(config['gmaps']) / "{chromosome}.ARS.gmap.gz")
    shell:
        '''
        for i in {{1..29}};
        do
          echo -e "pos\\tchr\\tcM" > ${{i}}.ARS.gmap; awk -v a=$i -F, '$1==a {{print sprintf("%.0f",$3*1000000)"\\t"$1"\\t"$5}}' {config[genetic_map]} | awk '$3=="NA"{$3=prev}{prev=$3}1' >> ${i}.ARS.gmap;
        done
        '''
#for i in {2..28}; do awk -i inplace '$3=="NA"{$3=prev}{prev=$3}1' ${i}.ARS.gmap; done


rule shapeit_imputation:
    input:
        vcf = str(PurePath(config['vcfs']) / "cohort.{chromosome}.vcf.gz"),
        gmap = str(PurePath(config['gmaps']) / "{chromosome}.ARS.gmap.gz")
    output:
        str(PurePath(config['imputed']) / "cohort.{chromosome}.shapeit.vcf.gz")
    threads: 4
    resources:
        mem_mb = 1000,
        walltime = '30'
    shell:
        '''
        /cluster/work/pausch/alex/software/shapeit4/bin/shapeit4.2 --input {input} --map {input.gmap} --region {wildcards.chromosome} --output {output} --thread {threads} --sequencing --mcmc-iterations 10b,1p,1b,1p,1b,1p,1b,1p,10m
        tmp_file=$(mktemp .XXXXXX)
        mv {output} $tmp_file
        bcftools reheader -f {config[reference]} -o {output} $tmp_file
        rm $tmp_file
        tabix -fp vcf {output}
        '''


rule bcftools_concat:
    input:
        (str(PurePath(config['imputed']) / f"cohort.{chromosome}.{{imputer}}.vcf.gz") for chromosome in range(1,30))
    output:
        multiext(str(PurePath(config['imputed']) / "cohort.autosomes.{imputer}.vcf.gz"),'','.tbi')
    threads: 4
    resources:
        mem_mb = 500,
        walltime = '10'
    shell:
        '''
        bcftools concat --threads {threads} -O z -o {output[0]} {input}
        tabix -p vcf {output[0]}
        '''
