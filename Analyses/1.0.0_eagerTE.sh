#!/bin/bash -l
#SBATCH -J eager_TE
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=8
#SBATCH --time=14-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/05_aDNA/slurms/slurm-%j-%x.out

# sbatch 1.0.0_eagerTE.sh 00_eager_twist /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta acrow_eager_list_twist
# sbatch 1.0.0_eagerTE.sh 00_eager_mybaits /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta acrow_eager_list_mybaits

echo $(date)
STARTTIME=$(date +%s)

module load charliecloud/0.30
module load nextflow
conda activate biotools

dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/05_aDNA"

cd ${dat}
mkdir ${1}
cd $1

nextflow run nf-core/eager -r 2.4.7 -profile conda --input /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/05_aDNA/${3}.tsv --fasta ${2} \
-c /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/05_aDNA/base_modified2.config \
--snpcapture_bed /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/05_aDNA/probes_232015.bed \
--clip_adapters_list /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/05_aDNA/adapterlist.txt 

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
