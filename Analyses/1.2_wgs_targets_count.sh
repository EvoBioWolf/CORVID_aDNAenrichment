#!/bin/bash -l
#SBATCH -J wgstarcov
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --time=2-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/05_aDNA/slurms/slurm-%j-%x.out

# sbatch 1.2_wgs_targets_count.sh probes_104k_80bp
# sbatch 1.2_wgs_targets_count.sh probes_104k 
# sbatch 1.2_wgs_targets_count.sh probes_232015 
echo $(date)
STARTTIME=$(date +%s)

conda activate biotools
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/05_aDNA"
wgs2="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/05_aDNA/00_eager_wgs2/results/ext_mapped"
wgs3="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/05_aDNA/00_eager_wgs3/results/ext_mapped"

cd ${dat}
cd 00_baitscomparison
cd wgs_mapped_reads
#because I want to include duplicate also for target efficiency computation
#i will run on raw bam first

samples=("BRW001_WGS" "DVT014_WGS" "NCP002_WGS" "VKP001_WGS")
for sample in "${samples[@]}"
do
samtools index ${wgs2}/${sample}_2_sorted_noqc.bam 
samtools index ${wgs3}/${sample}_3_sorted_noqc.bam
bedtools multicov -bams ${wgs2}/${sample}_2_sorted_noqc.bam -bed ${dat}/${1}.bed > ${sample}_${1}_2.counts
bedtools multicov -bams ${wgs3}/${sample}_3_sorted_noqc.bam -bed ${dat}/${1}.bed > ${sample}_${1}_3.counts
done

for sample in "${samples[@]}"
do
awk '{sum += $NF} END {print sum}' ${sample}_${1}_2.counts | paste - <(awk '{sum += $NF} END {print sum}' ${sample}_${1}_3.counts) | awk '{print $1 + $2}' > ${sample}_${1}_total.countsdone

paste BRW001_WGS_${1}_total.counts DVT014_WGS_${1}_total.counts NCP002_WGS_${1}_total.counts VKP001_WGS_${1}_total.counts > combined_${1}_totals.counts


ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
