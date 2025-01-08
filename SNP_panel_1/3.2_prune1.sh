#!/bin/bash -l
#SBATCH -J ldpru
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=3:00:00

###############################################################################################
#                                             Part1                                           # 
# This is the first part of neutral SNP selection. Unlinked 4fold, intronic and intergenic    #
# sites were extracted from each population                                                   #   
###############################################################################################

module load user_spack
module load plink2/1.9-beta6.10  
module load vcftools
module load bcftools/1.9 

dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/neutral"
ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/fst"

#strict thinning 1 SNP per 10kb for 4fold, intronic and intergenic sites 
cd $dat/neutral_perpop/summary/4fold
sort *_pos.txt | uniq > all_4fold_pos.txt
vcftools --vcf ${ref}/genome_crow.gatkHC.vqsrts99.sort5_incl.scaf5.repcont_SNP_ID.recode.vcf \
--positions all_4fold_pos.txt --keep $dat/poplist/all_subset_pop.txt --recode --out all_4fold
vcftools --vcf ${base}.recode.vcf --freq --out ${base}
awk 'NR>1' all_4fold.frq | awk '$3==2' | awk '{print $1, $2}' > all_4fold_biallelic_pos.txt
vcftools --vcf ${ref}/genome_crow.gatkHC.vqsrts99.sort5_incl.scaf5.repcont_SNP_ID.recode.vcf \
--positions all_4fold_biallelic_pos.txt --keep ${dat}/poplist/all_subset_pop.txt --recode --out all_4fold_biallelic
vcftools --vcf all_4fold_biallelic.recode.vcf --thin 1000 --recode --out 4fold_biallelic_thin1000
vcftools --vcf 4fold_biallelic_thin1000.recode.vcf --freq --out 4fold_biallelic_thin1000
awk 'NR>1' 4fold_biallelic_thin1000.frq | awk -F "\t|:" 'BEGIN{print "Chrom", "Pos", "Ref", "Alt", "RefFreq", "AltFreq", "LocusName"} {print $1, $2, $5, $7, $6, $8, $9=$1":"$2}' > 4fold_biallelic_thin1000.frq.summary
awk 'NR>1 {print $1, $2}' 4fold_biallelic_thin1000.frq.summary > 4fold_biallelic_thin1000_pos.txt

## Intron

cd $dat/neutral_perpop/summary/intron
for fname in *.vcf
do
base=${fname%.recode.vcf*}
vcftools --thin 10000 --vcf ${base}.recode.vcf --recode --out ./10kbthin/${base}_thin10kb
vcftools --vcf ./10kbthin/${base}_thin10kb.recode.vcf --freq --out ./10kbthin/${base}_thin10kb
awk 'NR>1 {print $1, $2}' ./10kbthin/${base}_thin10kb.frq > ./10kbthin/${base}_thin10kb.pos
done

## Intergene

cd $dat/neutral_perpop/summary/intergene
for fname in *.vcf
do
base=${fname%.recode.vcf*}
vcftools --thin 10000 --vcf ${base}.recode.vcf --recode --out ./10kbthin/${base}_thin10kb
vcftools --vcf ./10kbthin/${base}_thin10kb.recode.vcf --freq --out ./10kbthin/${base}_thin10kb
awk 'NR>1 {print $1, $2}' ./10kbthin/${base}_thin10kb.frq > ./10kbthin/${base}_thin10kb.pos
done

