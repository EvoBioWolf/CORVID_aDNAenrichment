#!/bin/bash -l
#SBATCH -J ldpru2
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --time=1-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/snp_panel_2/slurms/slurm-%j-%x.out

###############################################################################################
#                                             Part2                                           # 
# This script allows different combination of random SNP selection and population clustering. #
# This will generate sets 6-8 and the final panel set6up. Sets 1-5 require different script   #
# because cnx was not downsampled and/ or random selection was carried out after second round # 
# of LD prune (different order).                                                              #
###############################################################################################                 
#                                            To run                                           # 
# sbatch 1.2_prune2_addpop.sh cor3 cnx6 1614 4714 setadd03                                    #
# Has to be comparable to the previous design (1/7)                                           #
# Previous design (final, rd up): $8=intron=11300 / $9=intergene=33000                        #
# ${5}=Name of this set  #organise your folders!                                              #          
# final LD prune + biallelic check                                                            #
# merge with outlier_genic regions                                                            #
###############################################################################################

conda activate renv
module load plink2  
module load vcftools
module load bcftools

dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/snp_panel_2"
ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/05.1_recal/overlap/134inds_overlapped_filtered_norepeats.vcf.gz"
new="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"
hwe="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/05.1_recal/overlap/134inds_overlapped_filtered_norepeats_cnx3_hwedis.pos"

cd $dat/neutral_addpop/summary
mkdir ./intron/10kbthin/${5}
mkdir ./intergene/10kbthin/${5}

#collapse all neutral snps of each pop into a single file: no pop clustering
cd $dat/neutral_addpop/summary/intron/10kbthin
cp $dat/neutral_addpop/summary/intron/10kbthin/cnx6_thin10kb.recode.vcf $dat/neutral_addpop/summary/intron/10kbthin/${5}/${5}_10kbthin.recode.vcf

#random selection 
cd $dat/neutral_addpop/summary/intron/10kbthin/${5}
plink2 --vcf ${5}_10kbthin.recode.vcf --allow-extra-chr -allow-no-sex --out ${5}_random --thin-count ${3} --export vcf-4.2
#final 10kbprune
vcftools --thin 10000 --vcf ${5}_random.vcf --recode --out ${5}_random_10kbthin
vcftools --vcf ${5}_random_10kbthin.recode.vcf --freq --out ${5}_random_10kbthin
awk 'NR>1 {print $1, $2}' ${5}_random_10kbthin.frq > ${5}_random_10kbthin.pos

#extract sites from each population
vcftools --gzvcf ${ref} \
--positions ${5}_random_10kbthin.pos --keep $dat/poplist/all_subset_pop_cnx6.txt \
--min-alleles 2 --max-alleles 2 --recode --out all_${5}_biallelic
vcftools --vcf all_${5}_biallelic.recode.vcf --freq --out all_${5}_biallelic
awk 'NR>1' all_${5}_biallelic.frq | \
awk -F "\t|:" 'BEGIN{print "Chrom", "Pos", "Ref", "Alt", "RefFreq", "AltFreq", "LocusName"} {print $1, $2, $5, $7, $6, $8, $9=$1":"$2}' > all_${5}_biallelic.frq.summary

## Intergene
cd $dat/neutral_addpop/summary/intergene/10kbthin
cp $dat/neutral_addpop/summary/intergene/10kbthin/cnx6_thin10kb.recode.vcf $dat/neutral_addpop/summary/intergene/10kbthin/${5}/${5}_10kbthin.recode.vcf
cd $dat/neutral_addpop/summary/intergene/10kbthin/${5}
plink2 --vcf ${5}_10kbthin.recode.vcf --allow-extra-chr -allow-no-sex --out ${5}_random --thin-count ${4} --export vcf-4.2
#final 10kbprune
vcftools --thin 10000 --vcf ${5}_random.vcf --recode --out ${5}_random_10kbthin
vcftools --vcf ${5}_random_10kbthin.recode.vcf --freq --out ${5}_random_10kbthin
awk 'NR>1 {print $1, $2}' ${5}_random_10kbthin.frq > ${5}_random_10kbthin.pos

vcftools --gzvcf ${ref} \
--positions ${5}_random_10kbthin.pos --keep $dat/poplist/all_subset_pop_cnx6.txt \
--min-alleles 2 --max-alleles 2 --recode --out all_${5}_biallelic
vcftools --vcf all_${5}_biallelic.recode.vcf --freq --out all_${5}_biallelic
awk 'NR>1' all_${5}_biallelic.frq | awk -F "\t|:" 'BEGIN{print "Chrom", "Pos", "Ref", "Alt", "RefFreq", "AltFreq", "LocusName"} {print $1, $2, $5, $7, $6, $8, $9=$1":"$2}' > all_${5}_biallelic.frq.summary

## concat to generate final neutral panel
cd ${dat}/IRQ_neu
mkdir $5
cd $5
awk 'FNR==1 && NR!=1{next;}{print}' $dat/neutral_addpop/summary/intergene/10kbthin/${5}/all_${5}_biallelic.frq.summary $dat/neutral_addpop/summary/intron/10kbthin/${5}/all_${5}_biallelic.frq.summary > ${5}_biallelic_reduced_positions.txt
awk 'NR>1 {print $1, $2}' ${5}_biallelic_reduced_positions.txt > ${5}_biallelic_pos.txt

##adding 4fold sites

cat ${5}_biallelic_pos.txt ${dat}/neutral_addpop/summary/4fold/cnx6_thin1kb_nooverlap_pos.txt | sort -n > ${5}_inc4fold_biallelic_pos.txt

#check overlap with previous SNP panel 1
awk 'NR==FNR{c[$1, $2]++;next} !($1, $2) in c' /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/snp_panel_summary/05_final/snp_panel_104k_final.pos ${5}_inc4fold_biallelic_pos.txt > ${5}_biallelic_pos_nooverlap.txt

#check paralogs
awk 'NR==FNR{c[$1, $2]++;next} !($1, $2) in c' ${hwe} ${5}_biallelic_pos_nooverlap.txt > ${5}_biallelic_pos_nooverlap_hwe.txt

#plot pca for neutral snps including IRQ panel
cat /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/snp_panel_summary/05_final/neutral.pos ${5}_biallelic_pos_nooverlap_hwe.txt | sort -n > all_neutral_60K_$5.pos
vcftools --gzvcf ${ref} \
--positions all_neutral_60K_$5.pos --recode --out all_neutral_60K_$5
bcftools view -Oz -o all_neutral_60K_$5.vcf.gz all_neutral_60K_$5.recode.vcf
rm all_neutral_60K_$5.recode.vcf

#remove ori2
vcftools --gzvcf ${ref} \
--positions all_neutral_60K_$5.pos --keep ${dat}/poplist/poplist_131.txt --recode --out all_neutral_60K_$5_noori2
bcftools view -Oz -o all_neutral_60K_$5_noori2.vcf.gz all_neutral_60K_$5_noori2.recode.vcf
rm all_neutral_60K_$5_noori2.recode.vcf

Rscript ${new}/pca.R all_neutral_60K_$5 ${dat}/IRQ_neu/${5} ${dat}/IRQ_neu/${5} poplist_134
rm vcf.gds
Rscript ${new}/pca.R all_neutral_60K_$5_noori2 ${dat}/IRQ_neu/${5} ${dat}/IRQ_neu/${5} poplist_131
rm vcf.gds
