# BSW_analysis

## Notes on truth set validation

Merged autosomes from HD chip, other chips contributed very small amount of non-redundant sites

#### ID mapping
```
TRUTH sample -> bam sample
BSWCHEF120071057962 -> BSWCHEF120071057962
BSWCHEM110514060086 -> RM1894
BSWUSAM000000168052 -> RM1896
BSWUSAM000000186594 -> RM1899
```


#### Edits due to mismatch of reference and panel

```
sed -i '/2\t123539540/d' HD_chip.removed.vcf 
sed -i '/4\t36292062/d' HD_chip.removed.vcf 
sed -i '/9\t84677057/d' HD_chip.removed.vcf 
sed -i '/10\t11158957/d' HD_chip.removed.vcf 
sed -i '/11\t91524491/d' HD_chip.removed.vcf 
sed -i '/12\t63246637/d' HD_chip.removed.vcf 
sed -i '/14\t43010217/d' HD_chip.removed.vcf 
sed -i '/14\t70964024/d' HD_chip.removed.vcf 
sed -i '/15\t64749172/d' HD_chip.removed.vcf 
sed -i '/15\t69847415/d' HD_chip.removed.vcf 
sed -i '/19\t39399538/d' HD_chip.removed.vcf 
sed -i '/19\t55836350/d' HD_chip.removed.vcf 
sed -i '/20\t32338521/d' HD_chip.removed.vcf 
sed -i '/20\t43561728/d' HD_chip.removed.vcf 
sed -i '/20\t56749303/d' HD_chip.removed.vcf 
sed -i '/21\t8765426/d' HD_chip.removed.vcf 
sed -i '/26\t51144437/d' HD_chip.removed.vcf 
```


## Results

### Genotype concordance

Note that Specificity is calculated based on the quantity of false positives. I believe these count for calls from sequencing not found in the chip truth, i.e. we find many "FPs" because the sequencing approach generates more calls.

| Caller | Sensitivity | Specificity | F1       | GC       |
|--------|-------------|-------------|----------|----------|
| DV     | 0.989905    | 0.392515    | 0.11259  | 0.998023 |
| GATK   | 0.919109    | 0.408594    | 0.109807 | 0.998539 |
| % diff | 7.4%        | 4.0%        | 2.5%     | 0.05%    |

**NOTE** does this take into account the underlying prediction, or only the genotype? Appears to be latter.

### Call validation

Comparison to chip data as "truth"

| Call   | n Var    | avg Q | SNP Base | SNP Pos | Indel Pos |
|--------|----------|---------|-----------|---|---|
| Shared |    13.2m   | 56.7/7520 | 0.0407  | 0.0414  | 0.0014    |
| DV     |       839k   | 44.8 | 0.0341   | 0.0349  | 0.003     |
| GATK   |        441k  | 3629 | 0.00014  | 0.00092 | 0.0413737 |

DV and shared have similar profiles, while GATK is much lower.

Exclusive
![image](https://user-images.githubusercontent.com/29678761/119689647-50a7db00-be49-11eb-999d-f62daaec1876.png)
Missing
![image](https://user-images.githubusercontent.com/29678761/119689744-64ebd800-be49-11eb-953b-7e5f0a64b77a.png)
Shared
![image](https://user-images.githubusercontent.com/29678761/119689844-7b922f00-be49-11eb-883a-e79869652959.png)







