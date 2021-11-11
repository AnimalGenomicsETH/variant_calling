from pathlib import PurePath

rule all:
    input:
        expand(str(PurePath(config['imputed']) / 'cohort.{chromosome}.imputed.vcf.gz'),chromosome=range(1,30))

rule beagle_imputation:
    input:
        str(PurePath(config['vcfs']) / "cohort.{chromosome}.vcf.gz")
    output:
        str(PurePath(config['imputed']) / "cohort.{chromosome}.imputed.vcf.gz")
    params:
        lambda wildcards, output: PurePath(output[0]).with_suffix('').with_suffix('')
    threads: 12
    resources:
        mem_mb = 2500,
        walltime = '2:00'
    shell:
        '''
        java -jar -Xss25m -Xmx40G {config[beagle]} gl={input} nthreads={threads} out={params}
        tmp_file=$(mktemp .XXXXXX)
        mv {output} $tmp_file
        bcftools reheader -f {config[reference]} -o {output} $tmp_file
        tabix -fp vcf {output}
        '''
