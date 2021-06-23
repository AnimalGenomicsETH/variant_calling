bcftools concat --threads 4 -Ou $(ls CHR*/ref_alt_beagle5.vcf.gz | sort -V) | bcftools norm --threads 4 -f /cluster/work/pausch/alex/REF_DATA/ARS-UCD1.2_Btau5.0.1Y.fa -c x -o $1.vcf.gz -Oz
tabix -p vcf $1.vcf.gz
