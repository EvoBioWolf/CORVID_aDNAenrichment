#!/bin/bash -l
#SBATCH -J ldpru2
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=3:00:00

###############################################################################################
#                                             Part2                                           # 
# This script allows different combination of random SNP selection and population clustering. #
# This will generate sets 6-8 and the final panel set6up. Sets 1-5 require different script   #
# because cnx was not downsampled and/ or random selection was carried out after second round # 
# of LD prune (different order).                                                              #
###############################################################################################                 
#                                            To run                                           # 
# sbatch 3.3_prune2.sh cor1 cor2 cnx3 ori1 ori2 ori3 pec1 17250 54900 set6up                  #
# Population clustering if desired: $1-$7=cor1-pec1                                           #
# Random snp selection: $8=intron=17250 / $9=intergene=54900                                  #
# ${10}=Name of this set  #organise your folders!                                             #          
# final LD prune + biallelic check                                                            #
# merge with outlier_genic regions                                                            #
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

#collapse all neutral snps of each pop into a single file: no pop clustering
cd $dat/neutral_perpop/summary/intron/10kbthin
sort $1_thin10kb.pos $2_thin10kb.pos $3_thin10kb.pos $4_thin10kb.pos \
$5_thin10kb.pos $6_thin10kb.pos $7_thin10kb.pos | uniq > all_reduced_10kbthin_positions.txt
wc -l all_reduced_10kbthin_positions.txt
vcftools --vcf $ref/christen_vcf/genome_crow.gatkHC.vqsrts99.sort5_incl.scaf5.repcont_SNP_ID.recode.vcf \
--positions all_reduced_10kbthin_positions.txt --keep $dat/poplist/all_subset_pop.txt --recode --out $dat/neutral_perpop/summary/intron/10kbthin/${10}/${10}_10kbthin

#random selection 
cd ${10}
plink2 --vcf ${10}_10kbthin.recode.vcf --allow-extra-chr -allow-no-sex --out ${10}_random --thin-count $8 --export vcf-4.2
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

cd $dat/neutral_perpop/summary/intergene/10kbthin
sort $1_thin10kb.pos $2_thin10kb.pos $3_thin10kb.pos \
$4_thin10kb.pos $5_thin10kb.pos $6_thin10kb.pos $7_thin10kb.pos | uniq > all_reduced_10kbthin_positions.txt
wc -l all_reduced_10kbthin_positions.txt
vcftools --vcf $ref/christen_vcf/genome_crow.gatkHC.vqsrts99.sort5_incl.scaf5.repcont_SNP_ID.recode.vcf \
--positions all_reduced_10kbthin_positions.txt --keep $dat/poplist/all_subset_pop.txt --recode --out $dat/neutral_perpop/summary/intron/10kbthin/${10}/${10}_10kbthin

cd ${10}
plink2 --vcf ${10}_10kbthin.recode.vcf --allow-extra-chr -allow-no-sex --out ${10}_random --thin-count $9 --export vcf-4.2
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

## concat to generate final neutral panel

cd $pan
awk 'FNR==1 && NR!=1{next;}{print}' $dat/neutral_perpop/summary/intergene/10kbthin/${10}/all_${10}_biallelic.frq.summary \
$dat/neutral_perpop/summary/intron/10kbthin/${10}/all_${10}_biallelic.frq.summary \
$dat/neutral_perpop/summary/4fold/4fold_biallelic_thin1000.frq.summary > ${10}_biallelic_reduced_positions.txt
awk 'NR>1 {print $1, $2}' ${10}_biallelic_reduced_positions.txt > ${10}_biallelic_pos.txt

## merge

Rscript ../outlier_genic_neutral.R ${10}
rm VennDiagram*.log
