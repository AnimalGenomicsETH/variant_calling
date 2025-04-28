for B in *bam
do
  samtools bedcov -d 2 <(awk -v OFS='\t' '$1~/^[0-9XY]/ {print $1,0,$2}' /cluster/work/pausch/inputs/ref/BTA/UCD2.0/GCA_002263795.4_ARS-UCD2.0_genomic.fa.fai) $B | awk -v S=${B%.STAR.bam} '{print S,$1,$3,$4,$5}'
done
