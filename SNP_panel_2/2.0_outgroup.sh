#!/bin/bash -l
#SBATCH -J outgroup
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=6
#SBATCH --time=2-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/snp_panel_2/slurms/slurm-%j-%x.out
#SBATCH --mem-per-cpu=4763mb

sbatch 2.0_outgroup.sh 134inds_overlapped_filtered_norepeats_AMcrow_biallele 9pop outgroup_ascertained3

conda activate py2
module load vcftools/0.1.14-gcc8
module load bcftools
#picard.jar=2.25.7

ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/snp_panel_2"
neu="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/neutral"

echo $(date)
STARTTIME=$(date +%s)

cd ${dat}/$3
vcftools --keep ${dat}/poplist/${2}.txt --gzvcf /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/05.1_recal/overlap/${1}.vcf.gz --recode --out ${2}_overlapped_filtered_norepeats_AMcrow_biallele
bcftools view -Oz -o ${2}_overlapped_filtered_norepeats_AMcrow_biallele.vcf.gz ${2}_overlapped_filtered_norepeats_AMcrow_biallele.recode.vcf
bcftools index ${2}_overlapped_filtered_norepeats_AMcrow_biallele.vcf.gz
rm ${2}_overlapped_filtered_norepeats_AMcrow_biallele.recode.vcf
bcftools view --threads 6 -Oz -o ${2}_overlapped_filtered_norepeats_AMcrow_biallele_nomissing.vcf.gz --include "F_MISSING = 0" ${2}_overlapped_filtered_norepeats_AMcrow_biallele.vcf.gz

echo "count nomissing sites"
SitesG=$(zcat ${2}_overlapped_filtered_norepeats_AMcrow_biallele_nomissing.vcf.gz | grep -v '#' | wc -l)
echo $SitesG

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/60)) mins to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/3600)) hours to complete this task"


