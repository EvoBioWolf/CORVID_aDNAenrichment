#!/bin/bash -l
#SBATCH -J backup
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --time=2:00:00

###############################################################################################
#                                             Part3                                           # 
# This script allows different combination of random SNP selection and population clustering. #
# Backup neutral snps for mybaits                                                             #
###############################################################################################                 
#                                            To run                                           # 
# sbatch 3.4_backup.sh cor1 cor2 cnx3 ori1 ori2 ori3 pec1 20000 60000 backup set6up           #
# Population clustering if desired: $1-$7=cor1-pec1                                           #
# Random snp selection: $8=intron=20000 / $9=intergene=60000                                  #
# $10=Name of this set  #organise your folders!                                               #          
# final LD prune + biallelic check                                                            #
###############################################################################################

conda activate renv
module load user_spack
module load plink2/1.9-beta6.10  
module load vcftools
module load bcftools/1.9 

dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/neutral"
ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/fst"
pan="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/snp_panel_summary/03_neutral"

cd $dat/neutral_perpop/summary
mkdir ./intron/10kbthin/${10}
mkdir ./intergene/10kbthin/${10}

#extract remaining snps not selected from the pool of all_reduced_10kbthin snps
cd $dat/neutral_perpop/summary/intron/10kbthin/${11}
vcftools --vcf ${11}_random.vcf --freq --out ${11}_random
awk 'NR>1 {print $1, $2}' ${11}_random.frq > ${11}_random.pos
grep -wiv -f ${11}_random.pos $dat/neutral_perpop/summary/intron/10kbthin/all_reduced_10kbthin_positions.txt > $dat/neutral_perpop/summary/intron/10kbthin/${10}/${10}_all_10kbthin.pos
cd $dat/neutral_perpop/summary/intron/10kbthin/${10}
vcftools --vcf $ref/christen_vcf/genome_crow.gatkHC.vqsrts99.sort5_incl.scaf5.repcont_SNP_ID.recode.vcf \
--positions ${10}_all_10kbthin.pos --keep $dat/poplist/all_subset_pop.txt --recode --out ${10}_all_10kbthin

#random selection 
plink2 --vcf ${10}_all_10kbthin.recode.vcf --allow-extra-chr -allow-no-sex --out ${10}_random --thin-count $8 --export vcf-4.2
#final 10kbprune
vcftools --thin 10000 --vcf ${10}_random.vcf --recode --out ${10}_random_10kbthin
vcftools --vcf ${10}_random_10kbthin.recode.vcf --freq --out ${10}_random_10kbthin
awk 'NR>1 {print $1, $2}' ${10}_random_10kbthin.frq > ${10}_random_10kbthin.pos

#extract sites from each population
for i in $(cat $dat/neutral_perpop/panelpop)
do
echo $i
vcftools --gzvcf $dat/neutral_perpop/genome_crow.gatkHC.vqsrts99.sort5_incl.scaf5.repcont_SNP_ID_$i.vcf.gz \
--positions ${10}_random_10kbthin.pos --keep $dat/poplist/$i.txt --recode --out $i
vcftools --vcf $i.recode.vcf --freq --out $i
done
vcftools --vcf $ref/christen_vcf/genome_crow.gatkHC.vqsrts99.sort5_incl.scaf5.repcont_SNP_ID.recode.vcf \
--positions ${10}_random_10kbthin.pos --keep $dat/poplist/all_subset_pop.txt \
--min-alleles 2 --max-alleles 2 --recode --out all_${10}_biallelic
vcftools --vcf all_${10}_biallelic.recode.vcf --freq --out all_${10}_biallelic
awk 'NR>1' all_${10}_biallelic.frq | awk -F "\t|:" 'BEGIN{print "Chrom", "Pos", "Ref", "Alt", "RefFreq", "AltFreq", "LocusName"} {print $1, $2, $5, $7, $6, $8, $9=$1":"$2}' > all_${10}_biallelic.frq.summary

## Intergene

cd $dat/neutral_perpop/summary/intergene/10kbthin/${11}    
vcftools --vcf ${11}_random.vcf --freq --out ${11}_random        
awk 'NR>1 {print $1, $2}' ${11}_random.frq > ${11}_random.pos
grep -wiv -f ${11}_random.pos $dat/neutral_perpop/summary/intergene/10kbthin/all_reduced_10kbthin_positions.txt > $dat/neutral_perpop/summary/intergene/10kbthin/${10}/${10}_all_10kbthin.pos
cd $dat/neutral_perpop/summary/intergene/10kbthin/${10}
vcftools --vcf $ref/christen_vcf/genome_crow.gatkHC.vqsrts99.sort5_incl.scaf5.repcont_SNP_ID.recode.vcf \
--positions ${10}_all_10kbthin.pos --keep $dat/poplist/all_subset_pop.txt --recode --out ${10}_all_10kbthin

plink2 --vcf ${10}_all_10kbthin.recode.vcf --allow-extra-chr -allow-no-sex --out ${10}_random --thin-count $9 --export vcf-4.2
vcftools --thin 10000 --vcf ${10}_random.vcf --recode --out ${10}_random_10kbthin
vcftools --vcf ${10}_random_10kbthin.recode.vcf --freq --out ${10}_random_10kbthin
awk 'NR>1 {print $1, $2}' ${10}_random_10kbthin.frq > ${10}_random_10kbthin.pos

for i in $(cat $dat/neutral_perpop/panelpop)
do
echo $i
vcftools --gzvcf $dat/neutral_perpop/genome_crow.gatkHC.vqsrts99.sort5_incl.scaf5.repcont_SNP_ID_$i.vcf.gz \
--positions ${10}_random_10kbthin.pos --keep $dat/poplist/$i.txt --recode --out $i
vcftools --vcf $i.recode.vcf --freq --out $i
done
vcftools --vcf $ref/christen_vcf/genome_crow.gatkHC.vqsrts99.sort5_incl.scaf5.repcont_SNP_ID.recode.vcf \
--positions ${10}_random_10kbthin.pos --keep $dat/poplist/all_subset_pop.txt \
--min-alleles 2 --max-alleles 2 --recode --out all_${10}_biallelic
vcftools --vcf all_${10}_biallelic.recode.vcf --freq --out all_${10}_biallelic
awk 'NR>1' all_${10}_biallelic.frq | awk -F "\t|:" 'BEGIN{print "Chrom", "Pos", "Ref", "Alt", "RefFreq", "AltFreq", "LocusName"} {print $1, $2, $5, $7, $6, $8, $9=$1":"$2}' > all_${10}_biallelic.frq.summary

## concat to generate backup neutral snps
cd $pan
awk 'FNR==1 && NR!=1{next;}{print}' $dat/neutral_perpop/summary/intergene/10kbthin/${10}/all_${10}_biallelic.frq.summary \
$dat/neutral_perpop/summary/intron/10kbthin/${10}/all_${10}_biallelic.frq.summary > ${10}_biallelic_reduced_positions.txt
awk 'NR>1 {print $1, $2}' ${10}_biallelic_reduced_positions.txt > ${10}_biallelic_pos.txt

Rscript ../outlier_genic_neutral.R ${10}
rm VennDiagram*.log
