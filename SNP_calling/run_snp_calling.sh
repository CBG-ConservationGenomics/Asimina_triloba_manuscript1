#!/bin/sh                                                                                                                                                                                                 

###################################################################                                                                                                                                       
# pipeline for running snp-calling pipeline on re-sequencing data #                                                                                                                                       
#                                                                 #                                                                                                                                       
# usage:  snp_calling_April_13_2018.sh  $basedir $num_of_processors     #                                                                                                                                       
#                                                                 #                                                                                                                                       
###################################################################                                                                                                                                       
cd $1
CPU=$2

bowtie2-build ~/Projects/Strickler/Asimina_triloba/genome/Astri105_genome.fasta ~/Projects/Strickler/Asimina_triloba/genome/Astri105_genome

for file in `dir -d *_1.fq.gz ` ; do                                                                                                                      file2=`echo "$file" |sed 's/_1.fq.gz/_2.fq.gz/'`
    samfile=`echo "$file" | sed 's/_1.fq.gz/.sam/'`                                                                                                                                                       
    bowtie2 -p $CPU -x ~/Projects/Strickler/Asimina_triloba/genome/Astri105_genome -1 $file -2 $file2 -S $samfile                                                                  
done

#Do file compression and sorting
ls *.sam | parallel -j $CPU samtools view -bS -o {.}.bam {}                                                                                                                                              
ls *.bam | parallel -j $CPU samtools sort -o {.}.sort.bam {}

#Mark duplicate reads
ls *.sort.bam | parallel -j $CPU java -Djava.io.tmpdir=`pwd`~/tmp -jar /home/general/Programs/picard.jar MarkDuplicates INPUT={} OUTPUT={.}.md.bam METRICS_FILE={.}.metrics REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT

#Add read group info 
ls *.md.bam | parallel -j $CPU java -Djava.io.tmpdir=`pwd`~/tmp -jar /home/general/Programs/picard.jar AddOrReplaceReadGroups INPUT={} OUTPUT={}.rg.bam SORT_ORDER=coordinate RGID={} RGLB=1 RGPL=illumina RGPU=run RGSM={} RGCN=tom360 RGDS={}

#Index bam files
ls *.rg.bam | parallel -j $CPU samtools index {}

#Generate gvcf files
ls *.rg.bam |parallel -j $CPU java -jar /home/general/Programs/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar HaplotypeCaller -R ~/Projects/Strickler/Asimina_triloba/genome/Astri105_genome.fasta -I {} -O {.}.g.vcf -ERC GVCF

#CombineGVCFs and call SNPs
ls *.vcf > gvcf.list

java -Djava.io.tmpdir=`pwd`~/tmp -jar /home/general/Programs/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar CombineGVCFs -R ~/Projects/Strickler/Asimina_triloba/genome/Astri105_genome.fasta --variant gvcf.list -O combined_g.vcf

#Call SNPs
java -Djava.io.tmpdir=`pwd`~/tmp -jar /home/general/Programs/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar GenotypeGVCFs -R ~/Projects/Strickler/Asimina_triloba/genome/Astri105_genome.fasta -V combined_g.vcf -O output_snps.vcf
