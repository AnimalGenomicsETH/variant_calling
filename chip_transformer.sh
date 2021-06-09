zgrep "#" CHR1/ref_alt_beagle5.vcf.gz | sed -e "9r /cluster/work/pausch/alex/ARS.vcf.contigs" > $1
for i in {1..29}; do
  zgrep -v "#" CHR${i}/ref_alt_beagle5.vcf.gz >> $1
done
bcftools norm --threads 4 -f /cluster/work/pausch/alex/REF_DATA/ARS-UCD1.2_Btau5.0.1Y.fa -c x -o $1.gz -Oz $1
tabix -p vcf $1.gz
