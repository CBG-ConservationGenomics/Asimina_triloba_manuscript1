#https://speciationgenomics.github.io/filtering_vcfs/
#removed individuals with no coordinate or genotype data
bcftools view --force-samples --samples-file ^2remove all_called_snps_reheader4landscape.vcf > all_called_snps_reheader4landscape_rem.vcf
#Warn: exclude called for sample that does not exist in header: "VA_HighBT-1-1"... skipping
#Warn: exclude called for sample that does not exist in header: "VA_HighBT-1-2"... skipping
#Warn: exclude called for sample that does not exist in header: "VA_HighBT-1-3"... skipping
#Warn: exclude called for sample that does not exist in header: "VA_HighBT-1-4"... skipping
#Warn: exclude called for sample that does not exist in header: "VA_HighBT-1-5"... skipping
#Warn: exclude called for sample that does not exist in header: "PA_VenU-11-3"... skipping

#remove scaffolds and fungal contaminant contig (chr9)
grep -v Astri105S all_called_snps_reheader4landscape_rem.vcf > all_called_snps_reheader4landscape_rem2.vcf
grep -v Astri105C9 all_called_snps_reheader4landscape_rem2.vcf > all_called_snps_reheader4landscape_rem3.vcf

#filter for biallelic snps
vcftools --vcf all_called_snps_reheader4landscape_rem3.vcf --min-alleles 2 --max-alleles 2 --recode --out all_called_snps_reheader4landscape_rem4
#After filtering, kept 2072572 out of a possible 2167750 Sites

#filter on quality
vcftools --vcf all_called_snps_reheader4landscape_rem4.recode.vcf --recode --out all_called_snps_reheader4landscape_rem5 --minQ 30
#After filtering, kept 299 out of 299 Individuals
#After filtering, kept 2072572 out of a possible 2072572 Sites

#look at coverage per individual
vcftools --vcf all_called_snps_reheader4landscape_rem5.recode.vcf --depth --out individ_depth.txt

#filter on missing data
vcftools --vcf all_called_snps_reheader4landscape_rem5.recode.vcf --recode --out all_called_snps_reheader4landscape_rem6 --max-missing 0.9
#After filtering, kept 299 out of 299 Individuals #
#After filtering, kept 368602 out of a possible 2072572 Sites

#filter on min and max depth
#calc depth of coverage per individual
vcftools --vcf all_called_snps_reheader4landscape_rem6.recode.vcf --site-mean-depth --out site_mean_depth.txt
awk '{sum += $3; n++} END {if (n>0) print sum/n;}' site_mean_depth.txt.ldepth.mean
#=5.32


vcftools --vcf all_called_snps_reheader4landscape_rem6.recode.vcf --min-meanDP 5 --max-meanDP 11  --out all_called_snps_reheader4landscape7 --recode #max-meanDP=5.3*2
#After filtering, kept 299 out of 299 Individuals  
#After filtering, kept 212324 out of a possible 368602 Sites


#filter on minor allele frequency
#it is best practice to produce one dataset with a good MAF threshold and keep another without any MAF filtering at all.
#i didnt do this
