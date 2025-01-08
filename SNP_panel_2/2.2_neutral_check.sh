#!/bin/bash -l
#SBATCH -J neu_chk
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --time=1-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/snp_panel_2/slurms/slurm-%j-%x.out
#SBATCH --mem-per-cpu=4763mb

# sbatch 2.2_neutral_check.sh 9pop_overlapped_filtered_norepeats_AMcrow_biallele setadd03 outgroup_ascertained 9pop 8pop_noamcrow outgroup_ascertained4
# sbatch 2.2_neutral_check.sh 9pop_overlapped_filtered_norepeats_AMcrow_biallele setadd03 outgroup_NOTascertained 9pop 8pop_noamcrow outgroup_ascertained4

conda activate py2
module load vcftools/0.1.14-gcc8
module load bcftools
#picard.jar=2.25.7

dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/snp_panel_2"
neu="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/neutral"
new="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2"
ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/05.1_recal/overlap/134inds_overlapped_filtered_norepeats.vcf.gz"
hwe="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/05.1_recal/overlap/134inds_overlapped_filtered_norepeats_cnx3_hwedis.pos"

echo $(date)
STARTTIME=$(date +%s)

cd ${dat}/$6

#check overlap with previous design
awk 'NR==FNR{c[$1, $2]++;next} !($1, $2) in c' /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/snp_panel_summary/05_final/snp_panel_104k_final.pos neuall_$3.pos > neuall_$3_nooverlap.txt

#check overlap with selected IRQ SNPs
awk 'NR==FNR{c[$1, $2]++;next} !($1, $2) in c' ${dat}/IRQ_neu/${2}/${2}_biallelic_pos.txt neuall_${3}_nooverlap.txt > neuall_${3}_nooverlap_final.txt

#check transversion sites
vcftools --gzvcf ${1}.vcf.gz --positions neuall_${3}_nooverlap_final.txt --recode --out ${4}_neuall_${3}_nooverlap_final
vcftools --vcf ${4}_neuall_${3}_nooverlap_final.recode.vcf --keep ${dat}/poplist/${5}.pop --recode --out ${5}_neuall_${3}_nooverlap_final
vcftools --vcf ${4}_neuall_${3}_nooverlap_final.recode.vcf --TsTv-summary --out ${4}_neuall_${3}_nooverlap_final

# select transversion only
cat ${4}_neuall_${3}_nooverlap_final.recode.vcf | grep -v '#' | \
awk '(($4 == "A" && $5 == "T") || ($4 == "T" && $5 == "A") || ($4 == "C" && $5 == "G") || ($4 == "G" && $5 == "C") || ($4 == "C" && $5 == "A") || ($4 == "A" && $5 == "C") || ($4 == "G" && $5 == "T") || ($4 == "T" && $5 == "G") || $0 ~ /^#/ )' | \
awk '{print $1, $2}' >  neuall_${3}_nooverlap_final_TV.txt
vcftools --vcf ${4}_neuall_${3}_nooverlap_final.recode.vcf  --positions neuall_${3}_nooverlap_final_TV.txt --recode --out ${4}_neuall_${3}_nooverlap_final_TV
vcftools --vcf ${5}_neuall_${3}_nooverlap_final.recode.vcf  --positions neuall_${3}_nooverlap_final_TV.txt --recode --out ${5}_neuall_${3}_nooverlap_final_TV

#compare SFS of AMascertained sites and 3pops
vcftools --vcf ${4}_neuall_${3}_nooverlap_final.recode.vcf --freq --out ${4}_neuall_${3}_nooverlap_final
vcftools --vcf ${5}_neuall_${3}_nooverlap_final.recode.vcf --freq --out ${5}_neuall_${3}_nooverlap_final
vcftools --vcf ${4}_neuall_${3}_nooverlap_final.recode.vcf --counts --out ${4}_neuall_${3}_nooverlap_final
vcftools --vcf ${5}_neuall_${3}_nooverlap_final.recode.vcf --counts --out ${5}_neuall_${3}_nooverlap_final

awk 'NR>1' ${5}_neuall_${3}_nooverlap_final.frq.count |\
awk -F":|\t" 'BEGIN{print "Scaffold","Pos","Ref","Alt","RefCount","AltCount","TotalCount","RefFreq","AltFreq","LocusName"};{print $1,$2,$5,$7,$6,$8,$4,$6/$4,$8/$4,$1_$2}' > ${5}_neuall_${3}_nooverlap_final_AF_summary.txt 

#remove hwe disequilibrium sites 
vcftools --vcf ${5}_neuall_${3}_nooverlap_final.recode.vcf --exclude-positions ${hwe} --recode --out ${5}_neuall_${3}_nooverlap_final_hwe
vcftools --vcf ${5}_neuall_${3}_nooverlap_final_TV.recode.vcf --exclude-positions ${hwe} --recode --out ${5}_neuall_${3}_nooverlap_final_TV_hwe
vcftools --vcf ${5}_neuall_${3}_nooverlap_final_hwe.recode.vcf --freq --out ${5}_neuall_${3}_nooverlap_final_hwe
vcftools --vcf ${5}_neuall_${3}_nooverlap_final_hwe.recode.vcf --counts --out ${5}_neuall_${3}_nooverlap_final_hwe
vcftools --vcf ${5}_neuall_${3}_nooverlap_final_TV_hwe.recode.vcf --freq --out ${5}_neuall_${3}_nooverlap_final_TV_hwe
vcftools --vcf ${5}_neuall_${3}_nooverlap_final_TV_hwe.recode.vcf --counts --out ${5}_neuall_${3}_nooverlap_final_TV_hwe
awk 'NR>1' ${5}_neuall_${3}_nooverlap_final_TV_hwe.frq.count |\
awk -F":|\t" 'BEGIN{print "Scaffold","Pos","Ref","Alt","RefCount","AltCount","TotalCount","RefFreq","AltFreq","LocusName"};{print $1,$2,$5,$7,$6,$8,$4,$6/$4,$8/$4,$1_$2}' > ${5}_neuall_${3}_nooverlap_final_TV_hwe_AF_summary.txt
awk 'NR>1' ${5}_neuall_${3}_nooverlap_final_hwe.frq.count |\
awk -F":|\t" 'BEGIN{print "Scaffold","Pos","Ref","Alt","RefCount","AltCount","TotalCount","RefFreq","AltFreq","LocusName"};{print $1,$2,$5,$7,$6,$8,$4,$6/$4,$8/$4,$1_$2}' > ${5}_neuall_${3}_nooverlap_final_hwe_AF_summary.txt
awk 'NR>1 {print $1, $2}' ${5}_neuall_${3}_nooverlap_final_TV_hwe_AF_summary.txt > neuall_${3}_nooverlap_final_TV_hwe.txt
awk 'NR>1 {print $1, $2}' ${5}_neuall_${3}_nooverlap_final_hwe_AF_summary.txt > neuall_${3}_nooverlap_final_hwe.txt
vcftools --gzvcf ${ref} --positions neuall_${3}_nooverlap_final_hwe.txt --recode --out 134inds_neutral_${3}_${4}_hwe
bcftools view -Oz -o 134inds_neutral_${3}_${4}_hwe.vcf.gz 134inds_neutral_${3}_${4}_hwe.recode.vcf
vcftools --gzvcf ${ref} --positions neuall_${3}_nooverlap_final_TV_hwe.txt --recode --out 134inds_neutral_${3}_${4}_TV_hwe
bcftools view -Oz -o 134inds_neutral_${3}_${4}_TV_hwe.vcf.gz 134inds_neutral_${3}_${4}_TV_hwe.recode.vcf

#see pop structure in 134inds with these selected snps
vcftools --gzvcf ${ref} --positions neuall_${3}_nooverlap_final.txt --recode --out 134inds_neutral_${3}_${4}
bcftools view -Oz -o 134inds_neutral_${3}_${4}.vcf.gz 134inds_neutral_${3}_${4}.recode.vcf
vcftools --gzvcf 134inds_neutral_${3}_${4}.vcf.gz --counts --out 134inds_neutral_${3}_${4}
awk 'NR>1' 134inds_neutral_${3}_${4}.frq.count |\
awk -F":|\t" 'BEGIN{print "Scaffold","Pos","Ref","Alt","RefCount","AltCount","TotalCount","RefFreq","AltFreq","LocusName"};{print $1,$2,$5,$7,$6,$8,$4,$6/$4,$8/$4,$1_$2}' > 134inds_neutral_${3}_${4}_final_AF_summary.txt
vcftools --gzvcf ${ref} --positions neuall_${3}_nooverlap_final_TV.txt --recode --out 134inds_neutral_${3}_${4}_TV
bcftools view -Oz -o 134inds_neutral_${3}_${4}_TV.vcf.gz 134inds_neutral_${3}_${4}_TV.recode.vcf

conda activate renv
Rscript ${new}/pca.R 134inds_neutral_${3}_${4}_TV_hwe ${dat}/$6 ${dat}/$6 poplist_134
rm vcf.gds
Rscript ${new}/pca.R 134inds_neutral_${3}_${4}_hwe ${dat}/$6 ${dat}/$6 poplist_134
rm vcf.gds

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/60)) mins to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/3600)) hours to complete this task"


