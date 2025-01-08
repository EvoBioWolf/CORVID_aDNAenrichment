#!/bin/bash -l
#SBATCH -J ldpru
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --time=1-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/snp_panel_2/slurms/slurm-%j-%x.out

###############################################################################################
#                                             Part1                                           # 
# This is the first part of neutral SNP selection. Unlinked 4fold, intronic and intergenic    #
# sites were extracted from each population: only IRQ                                         #   
###############################################################################################

# module load user_spack
module load plink2
module load vcftools
module load bcftools

dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/snp_panel_2"
ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/fst"

#strict thinning 1 SNP per 10kb for intronic and intergenic sites 

#extract intronic & intergenic sites from the VCF of each pop
cd $dat/neutral_addpop
for i in $(cat panelpop)
do
echo $i
vcftools --gzvcf $i.vcf.gz \
--positions ./summary/intron/${i}_pos.txt --recode --out ./summary/intron/$i
vcftools --gzvcf $i.vcf.gz \
--positions ./summary/intergene/${i}_pos.txt --recode --out ./summary/intergene/$i
done

## Intron

cd $dat/neutral_addpop/summary/intron
for fname in *.vcf
do
base=${fname%.recode.vcf*}
vcftools --thin 10000 --vcf ${base}.recode.vcf --recode --out ./10kbthin/${base}_thin10kb
vcftools --vcf ./10kbthin/${base}_thin10kb.recode.vcf --freq --out ./10kbthin/${base}_thin10kb
awk 'NR>1 {print $1, $2}' ./10kbthin/${base}_thin10kb.frq > ./10kbthin/${base}_thin10kb.pos
done

## Intergene

cd $dat/neutral_addpop/summary/intergene
for fname in *.vcf
do
base=${fname%.recode.vcf*}
vcftools --thin 10000 --vcf ${base}.recode.vcf --recode --out ./10kbthin/${base}_thin10kb
vcftools --vcf ./10kbthin/${base}_thin10kb.recode.vcf --freq --out ./10kbthin/${base}_thin10kb
awk 'NR>1 {print $1, $2}' ./10kbthin/${base}_thin10kb.frq > ./10kbthin/${base}_thin10kb.pos
done

##4fold 

cd ${dat}/neutral_addpop
for i in $(cat panelpop)
do
echo $i
vcftools --gzvcf $i.vcf.gz \
--positions ./summary/4fold/${i}_pos.txt --recode --out ./summary/4fold/$i
done

cd $dat/neutral_addpop/summary/4fold
for i in $(cat ${dat}/neutral_addpop/panelpop)
do
base=${i%.recode.vcf*}
vcftools --thin 1000 --vcf ${base}.recode.vcf --recode --out ${base}_thin1kb
vcftools --vcf ${base}_thin1kb.recode.vcf --freq --out ${base}_thin1kb
awk 'NR>1 {print $1, $2}' ${base}_thin1kb.frq > ${base}_thin1kb.pos
done

awk 'NR==FNR{c[$1, $2]++;next} !($1, $2) in c' /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/snp_panel_summary/03_neutral/set6up_biallelic_pos.txt cnx6_thin1kb.pos > cnx6_thin1kb_nooverlap_pos.txt

