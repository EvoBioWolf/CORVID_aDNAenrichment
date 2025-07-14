#!/bin/bash -l
#SBATCH -J ANCangsd
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=12
#SBATCH --time=2-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/05_aDNA/slurms/slurm-%j-%x.out

# sbatch 1.1.1_angsd_baitscomparison.sh twist 104kpanel probes_104k_SNPsite.1based
# sbatch 1.1.1_angsd_baitscomparison.sh mybaits 104kpanel probes_104k_SNPsite.1based
# sbatch 1.1.1_angsd_baitscomparison.sh shotgun_noudg 104kpanel probes_104k_SNPsite.1based

echo $(date)
STARTTIME=$(date +%s)

dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/05_aDNA"
ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/05_aDNA/genome_HC_allpaths41687_v2.5.fasta"

conda activate biotools
cd ${dat}
cd 00_baitscomparison/angsd


#geno
# angsd -bam bam_${1}.filelist -GL 2 -doMaf 2 -doMajorMinor 1 -ref ${ref} -doGeno 4 -doPost 1 -doGlf 2 -doCounts 1 -doPlink 2 \
# -minMapQ 20 -minQ 20 -uniqueOnly 1 -remove_bads 1 \
# -nThreads 8 -rf ${dat}/${3} -out ${1}_${2}_geno_maxmis_q20

#geno dp3
angsd -bam bam_${1}_4.filelist -GL 2 -doMaf 2 -doMajorMinor 1 -ref ${ref} -doGeno 4 -doPost 1 -doGlf 2 -doCounts 1 -doPlink 2 \
-minMapQ 20 -minQ 20 -uniqueOnly 1 -remove_bads 1 -geno_minDepth 3 \
-nThreads 12 -rf ${dat}/${3} -out ${1}_${2}_geno_maxmis_q20_dp3

#get allele count for individual sample for just the 4 selected samples
angsd -out ${1}_${2}_geno_maxmis_q20_dp3 -doCounts 1 -dumpCounts 4 -bam bam_${1}_4.filelist -minMapQ 20 -minQ 20 -uniqueOnly 1 -remove_bads 1 \
-nThreads 8 -rf ${dat}/${3}


ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
