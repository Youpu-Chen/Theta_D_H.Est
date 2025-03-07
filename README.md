>**Theta_D_H.Est**
>
>https://github.com/panyuwen/Theta_D_H.Est 
>
>Calculator for statictics including Theta, D, H, and so on, based on the VCF input. 
>
>If you use the **Theta_D_H.Est** in your research work, please cite at least one of the following paper(s):
>
>- [Genomic diversity and post-admixture adaptation in the Uyghurs](https://academic.oup.com/nsr/article/9/3/nwab124/6368880) (*National Science Review*, 2021)
>- [Lineage-specific positive selection on *ACE2* contributes to the genetic susceptibility of COVID-19](https://academic.oup.com/nsr/advance-article/doi/10.1093/nsr/nwac118/6623880) (*National Science Review*, 2022)

---

### Notes

- for common users, the autosomes (and chrX) variants set can be applied directly, but for the haploid varinats dataset like chrY, the haploid-filling process is need which requires the sample file with two columns: 1) sampleID, 2) which haplotype to use
  - `hapFill.py` can help you to conduct this process
  - `ThetaEstDemo.py` is a numpy-based scripts which speed up the theta calculation. Only allowed for theta pi and theta k (theta watterson) calculation right now.

### Description: 

To calculate the following statistic within given regions:

- number of sequence
- number of genetic markers
- number of singleton
- Theta PI
- Theta K
- number of segregating site
- number of haplotype
- Haplotype diversity
- Fay & Wu's H
- normalized Fay & Wu's H
- Fu & Li's F(F*)
- Fu & Li's D(D*)
- Tajima's D (with p-values and BH-corrected P-values)

### Input:

1. phased VCF.gz file `required` 
2. target regions to be analyzed (file) `required` 
3. samples to be kept (file) `optional` 
4. haplotypes to be kept (file) `optional`
5. length of sliding windows (bp) and increment (bp) (string) `optional` 
6. whether the state of ancestral/derived allele is determined (string) `optional` 
7. name of output file (string) `optional` 

### Usage:
	$ ./Theta_D_H.Est -h 
	$ ./Theta_D_H.Est \
	    --gzvcf phased.vcf.gz \
	    --region region.txt \
	    --samples sample.info \
	    --haps haps.info \
	    --window_shift 50000@10000 \
	    --outgroup N \
	    --out output.txt

>**details about all these arguments** 
>
>**--gzvcf** (required): phased VCF.gz file, in GT format, like "1|0". Note: no duplicate physical positions, separate autosomes and X chromosome  
>**--region** (required): file including target regions to be analyzed. 4 columns: `region ID` `chrom ID` `start pos` `end pos`. no header line, tab or space delimited  
>**--samples** (optional): file including samples to be analyzed. 1 or 2 column: `sample ID` `gender, optional`. gender code, 1: male; 2: female. no header line. **default:** all samples in the VCF.gz file  
>**--haps** (optional): file including haplotypes to be analyzed. 2 columns: `sample ID` `haplotype index (1/2)`, 1: first haplotype; 2: second haplotype. no header line. **default:** all haplotypes in the VCF.gz file
>**--window_shift** (optional): `window_length@increment`, to partition the target region(s) into sliding windows of `window_length ` advanced by `increment`. Then calculate all of the statistic within the sliding windows. **default:** calculate all of the statistic within target region(s)  
>**--outgroup** (optional): `Y / N`, whether the ancestral/derived allele is determined, required for Fu&Li's and Fay&Wu's tests. **default:** N  
>**--out** (optional): output file name. output file would be gzipped. **default:** out.txt

### Tips & Notes:
1. to calculate statistic across the whole chromosome / genome, you can define the chromosome boundaries in the region file, 
	
	>ID1 21 1 48129895  
	>ID2 22 1 51304566  
	
	then specify the window size and increment (e.g., 100000@5000). This script will help you to calculate all of the statistic within sliding windows of 100KB in length shift by 5KB.   
	
2. heterozygotes are not allowed for male X chromosome. But this script can't help you to check the data format.   

3. all individuals would be considered as diploid, if no gender information is provided  

4. if both "--haps" and "--samples" are used, only samples provided by both arguments are remained for analysis 

5. if "--haps" is used, it would be probably  unnecessary to use "--samples", (for chrX, only if males are in homo format)

6. for chrX in males, only the first haplotypes would be used, no more second haplotypes  

7. if the ancestral/derived alleles are not determined, there won't output Fay & Wu's H, normalized Fay & Wu's H, Fu & Li's F, or Fu & Li's D (output Fu & Li's F\* and Fu & Li's D\* instead)

8. if the ancestral/derived alleles are determined, number of singletons will be estimated as the number of derived singletons  

9. chromosome IDs in the VCF file and region file should be coded in the same way (i.e., "1" is different from "chr1")

10. linear computational complexity. you are suggested to filter homozygous loci to speed up the program. it may takes <1h and ~2Gb memory for chromosome 1 of CHB (50000@20000, ~12.5K sliding windows, 103 individuals)

11. this program is compiled in centos7, older systems may not be supported. 

12. If you have problem using the compiled program, it can still be run in the following way: `python2 Theta_D_F_H.py2.py [--options]`; OR `python3 Theta_D_F_H.py3.py [--options]`. All the required packages are accessible in conda.

---
By: Yuwen Pan, 2021  
Contact: panyuwen.x@gmail.com

