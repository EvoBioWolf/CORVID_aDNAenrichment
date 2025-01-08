#!/bin/bash -l
#SBATCH -J eagerTEbam
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=12
#SBATCH --time=2-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/05_aDNA/slurms/slurm-%j-%x.out

# for i in mybaits twist
# do
# sbatch 1.0.5_eagerTEbam.sh 00_eager_${i} /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/05_aDNA/genome_HC_allpaths41687_v2.5.fasta /dss/dsslegfs01/pr53da/pr53da-dss-0018/#projects/2020__ancientDNA/05_aDNA/acrow_eager_list_${i}_bam.tsv probes_232015_SNPsite
# sbatch 1.0.5_eagerTEbam.sh 00_eager_${i} /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/05_aDNA/genome_HC_allpaths41687_v2.5.fasta /dss/dsslegfs01/pr53da/pr53da-dss-0018/#projects/2020__ancientDNA/05_aDNA/acrow_eager_list_${i}_bam.tsv probes_104k_80bp
# sbatch 1.0.5_eagerTEbam.sh 00_eager_${i} /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/05_aDNA/genome_HC_allpaths41687_v2.5.fasta /dss/dsslegfs01/pr53da/pr53da-dss-0018/#projects/2020__ancientDNA/05_aDNA/acrow_eager_list_${i}_bam.tsv probes_104k_SNPsite
# sbatch 1.0.5_eagerTEbam.sh 00_eager_${i} /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/05_aDNA/genome_HC_allpaths41687_v2.5.fasta /dss/dsslegfs01/pr53da/pr53da-dss-0018/#projects/2020__ancientDNA/05_aDNA/acrow_eager_list_${i}_bam.tsv probes_232015
# sbatch 1.0.5_eagerTEbam.sh 00_eager_${i} /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/05_aDNA/genome_HC_allpaths41687_v2.5.fasta /dss/dsslegfs01/pr53da/pr53da-dss-0018/#projects/2020__ancientDNA/05_aDNA/acrow_eager_list_${i}_bam.tsv probes_104k
# done

echo $(date)
STARTTIME=$(date +%s)

#module load charliecloud
#module load purge
module load nextflow
conda activate biotools

dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/05_aDNA"

cd ${dat}
mkdir ${1}
cd $1

#with snpsite bed files and trimmed bam
unset DISPLAY
nextflow run nf-core/eager -profile conda -r 2.4.7 --input ${3} --fasta ${2} \
-c /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/05_aDNA/base_modified2.config \
--snpcapture_bed /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/05_aDNA/${4}.bed \
--mtnucratio_header chrM --run_mtnucratio \
--run_trim_bam --bamutils_clip_single_stranded_none_udg_left 5 --bamutils_clip_single_stranded_none_udg_right 5 

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
