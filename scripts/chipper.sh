module load plink
plink --bfile /nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/SUISAG/transfer_2021/SNP_2021/prepared/Suisag_${1} --recode vcf --threads 2 --memory 6000 --keep-allele-order  --out $1
awk '$1!="0"&&$1!="20"&&$1!="21"&&$1!="22"' $1.vcf | bgzip --threads 2 -c > $1.vcf.gz
tabix -fp vcf $1.vcf.gz
#bcftools view -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18 -O z FBFv2.vcf.gz | 
bcftools reheader -f /cluster/work/pausch/inputs/ref/SSC/11.1/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.fai $1.vcf.gz > $1.new.vcf.gz
tabix -fp vcf $1.new.vcf.gz

bcftools +fixref $1.new.vcf.gz -o $1.fixed.vcf.gz -- -f /cluster/work/pausch/inputs/ref/SSC/11.1/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa -m flip -d
tabix -p vcf $1.fixed.vcf.gz
#bcftools norm -f /cluster/work/pausch/inputs/ref/SSC/11.1/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa -c x $1.new.vcf.gz -o $1.normed.vcf.gz
#tabix -p vcf $1.normed.vcf.gz
