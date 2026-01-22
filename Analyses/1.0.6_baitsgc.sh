#!/bin/bash -l
#SBATCH -J gc
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=4
#SBATCH --time=2-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/05_aDNA/slurms/slurm-%j-%x.out

# sbatch 1.0.2_gc.sh 00_eager_mybaits 104K_80bp_panel probes_104k_80bp mybaits
# sbatch 1.0.2_gc.sh 00_eager_twist 104K_80bp_panel probes_104k_80bp twist

echo $(date)
STARTTIME=$(date +%s)

conda activate biotools
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/05_aDNA"
ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/05_aDNA/genome_HC_allpaths41687_v2.5.fasta"

cd $dat
cd ${1}/results/${2}/trimmed_bam

for i in *.trimmed.bam
do
base=${i%.trimmed.bam*}
# samtools view -@ 4 -bF 4 -L ${dat}/${3}.bed -b ${base}.trimmed.bam > ${base}.trimmed.104k80bp.bam
# samtools index -@ 4 ${base}.trimmed.104k80bp.bam
java -jar ${dat}/picard.jar CollectGcBiasMetrics \
    I=${base}.trimmed.104k80bp.bam \
    O=${dat}/00_baitscomparison/gc/${base}_gc_bias_metrics.txt \
    CHART=${dat}/00_baitscomparison/gc/${base}_gc_bias_chart.pdf \
    S=${dat}/00_baitscomparison/gc/${base}_summary_metrics.txt \
    R=${ref} \
    SCAN_WINDOW_SIZE=80 \
    VALIDATION_STRINGENCY=LENIENT 
done

#gc_at_dropout.py
# samtools depth -a -b ${dat}/probes_104k_SNPsite.bed $(echo *.trimmed.104k80bp.bam) > ${dat}/00_baitscomparison/coverage_104k_${4}.txt

# conda activate py3.8
# cd ${dat}/00_baitscomparison/gc
# gc_at_dropout.py

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"

