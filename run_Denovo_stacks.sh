#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 16
#$ -v PATH
#$ -V
#$ -q main.q@@himem

date
hostname


#create a new folder for the quality check reports and files of the Raw data

mkdir -p -f   FASTQC_RD

#run fastqc to check the quality of reads

fastqc -o ./FASTQC_RD/  *.fastq
#clean the Raw data using cutadapt
mkdir -p -f   clean_data

#please note, that GBS_Rwadata is the folder where your raw data are saved 

for f in GBS_Rawdata\/*.fastq # for each sample

do
    n=${f%%.fastq} # strip part of file name
    cutadapt -a  AGATCGGAAGAGC -q 20,20  -m 45 -o /GIVE/YOUR/OUTPUT_PATH/clean_data/${n}_cut.fastq $1 ${n}.fastq 
                                                       
done

#create a new folder for the quality check reports and files

mkdir -p -f   FASTQC_CLND

#run fastqc to check the quality of reads

fastqc -o ./FASTQC_CLND/  *_cut.fastq

# run Multiqc with --interractive parameter and give the floder or the path where you saved your fastq reports
multiqc ./FASTQC_CLND/ --interractive
#automate stacks
#test from Khaoula
mkdir -p -f   output_stacks4
read -r -a samples < samples_ID.txt

for i in "${samples[@]}"
do
  echo "$i"
done
#ustacks
	
i=1
for sample in $samples
do
   ustacks -t fastq -f ./clean_data/${sample}.fastq -o ./output_stacks4/ -i $i -m 3 -M 2 -p 16 --force-diff-len 
   let "i+=1";
done

#cstacks
cstacks -n 6 -P ./output_stacks4/ -M ./popumap/popumap.txt -p 8 --k_len 15

#sstacks
sstacks -P ./output_stacks4/ -M ./popumap/popumap.txt -p 8


#tsv2bam
tsv2bam -P ./output_stacks4/ -M ./popumap/popumap.txt  -t 8

#gstacks
gstacks -P ./output_stack4/ -M ./popumap/popumap.txt  -t 8 

#populations
populations -P ./output_stacks4/ -M ./popumap/popumap.txt -r 0.65 --vcf --genepop --structure --fstats --hwe -t 8

#vcftools
vcftools --vcf ./output_stacks4/populations.snps.vcf  --max-missing 0.5 --maf 0.05  --recode --recode-INFO-all --out ./output_stacks1/filtered_popnovo_snps_005maf


#populations(convert)
populations -V ./output_stacks4/filtered_popnovo_snps_005maf.recode.vcf -O ./out_treemix/ -M ./popumap/popumap.txt --treemix --structure  --phylip --hwe  -t 8

#prepare data for plink analyses

grep -v "#" filtered_popnovo_snps_005maf.recode.vcf | sed -E 's/^/X_/g' > filtered_popnovo_snps_005maf.snps.vcf
grep "#" filtered_popnovo_snps_005maf.recode.vcf >  edit_head
cat edit_head filtered_popnovo_snps_005maf.snps.vcf >> filtered_popnovo_snps_005maf_plink.snps.vcf

#identify prune sites

./plink2 --vcf filtered_popnovo_snps_005maf_plink.snps.vcf --double-id --aec --set-missing-var-ids @:# --out popnovo_snps_005maf

#run  pc analyses with plink
#this will generate eigenvalue and eigenvectors files
./plink2 --vcf rfiltered_popnovo_snps_005maf_plink.snps.vcf --double-id --aec --set-missing-var-ids @:# --extract popnovo_snps_005maf.prune.in --make-bed --pca 10 --out plink_005maf
#for huge samples use the approx option to force generating pc analyses
.plink2 --vcf filtered_popnovo_snps_005maf_plink.snps.vcf --pca approx 10 --aec --set-missing-var-ids @:# --extract popnovo_snps_005maf.prune.in --make-bed  --out  plink_005maf 
