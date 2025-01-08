#!/bin/bash -l
#SBATCH -J singlehit
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --time=30:00

### sbatch 5.1_outgroup_singlehit.sh set6up ###
### sbatch 5.1_outgroup_singlehit.sh backup ###
### keeping only single hits mapped to outgroups New Caledonian crow and American crow ###

module load user_spack
module load vcftools

ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/fst"
ref2="/dss/dsslegfs01/pr53da/pr53da-dss-0018/assemblies/Corvus.cornix/genome/v2/v2.5"
scaff2chr="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/scaffold2chr"
neu="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/snp_panel_summary/03*"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/snp_panel_summary/04_combined"
final="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/snp_panel_summary/06_singlehits"

cd $dat
cd $1_snp_panel

#grep only single hits
grep -A 1 "# 1 hits found" bCorMon1_$1_mapped | grep -v "# 1 hits found\|--" > bCorMon1_mapped_singlehits
grep -A 1 "# 1 hits found" ASM69197v1_$1_mapped | grep -v "# 1 hits found\|--" > ASM69197v1_mapped_singlehits

#pull out additional snps removed due to multiple hits
awk -F ":|-" '{print $1, $3-60}' bCorMon1_mapped_singlehits > bCorMon1_mapped_singlehits_snp.pos
awk -F ":|-" '{print $1, $3-60}' ASM69197v1_mapped_singlehits > ASM69197v1_mapped_singlehits_snp.pos
sort *_mapped_singlehits_snp.pos | uniq > all_mapped_singlehits_snp.pos
wc -l all_mapped_singlehits_snp.pos
grep -wiv -f all_mapped_singlehits_snp.pos snp_panel_$1_trimmed_mapped_pos.txt > all_multihits_snps_removed.pos
wc -l all_multihits_snps_removed.pos
grep -wiv -f all_multihits_snps_removed.pos snp_panel_$1_trimmed_mapped_pos.txt > snp_panel_$1_trimmed_mapped_singlehits_pos.txt
wc -l snp_panel_$1_trimmed_mapped_singlehits_pos.txt #101399_this makes sure both reference must be mapped
cp snp_panel_$1_trimmed_mapped_singlehits_pos.txt ${final}_$1

#generate final positions
cd ${final}_$1
awk 'NR==FNR{c[$1, $2]++;next};c[$3, $4] > 0' snp_panel_$1_trimmed_mapped_singlehits_pos.txt ../02_genic/all_chrz_fst_dxy99.95_knief_genic_38K_biallelic_only.txt > snp_panel_outlier_genic.txt

#identify unqiue neutral positions - non-overlapping with outlier / genic
awk '{print $3, $4}' snp_panel_outlier_genic.txt > outlier_genic.pos
grep -wiv -f outlier_genic.pos snp_panel_$1_trimmed_mapped_singlehits_pos.txt > neutral.pos
cat <(awk '{print $1":"$2, "neutral", $1, $2}' neutral.pos) snp_panel_outlier_genic.txt > snp_panel_outlier_genic_neutral.txt

#insert chr info
awk '{print $3, $4, $2, $3, $1}' snp_panel_outlier_genic_neutral.txt > snp_panel_outlier_genic_neutral_chr.txt
cd $scaff2chr
for fname in *.scaffolds
do
base=${fname%.scaffolds*}
for i in $(cat $base.scaffolds)
do
echo $i $base
awk -v scaff=$i -v chr=$base '{gsub("^"scaff"$", chr, $4)} 1' ${final}_$1/snp_panel_outlier_genic_neutral_chr.txt > ${final}_$1/tmp.txt && mv -v ${final}_$1/tmp.txt ${final}_$1/snp_panel_outlier_genic_neutral$
done
done

cd $final_$1
#recall snps from vcf
vcftools --vcf $ref/christen_vcf/genome_crow.gatkHC.vqsrts99.sort5_incl.scaf5.repcont_SNP_ID.recode.vcf \
--positions snp_panel_$1_trimmed_mapped_pos.txt --recode --out snp_panel_$1_final

#outlier and genic only
vcftools --vcf $ref/christen_vcf/genome_crow.gatkHC.vqsrts99.sort5_incl.scaf5.repcont_SNP_ID.recode.vcf --positions outlier_genic.pos --recode --out snp_panel_outlier_genic_$1_final

#neutral only
vcftools --vcf $ref/christen_vcf/genome_crow.gatkHC.vqsrts99.sort5_incl.scaf5.repcont_SNP_ID.recode.vcf --positions neutral.pos --recode --out snp_panel_neutral_$1_final
