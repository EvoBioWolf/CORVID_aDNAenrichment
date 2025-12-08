#!/bin/bash -l
#SBATCH -J bedtools
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH -o PATH/01_probes/snp_panel_twist/slurms/slurm-%j-%x.out
#SBATCH --mem-per-cpu=4763mb

# sbatch 1.1_probes.sh snp_panel_combined_232k_0based

conda activate py3.8
module load bedtools2/2.27.1-gcc8 
#picard.jar=2.25.7

dat="PATH/01_probes/snp_panel_twist"
neu="PATH/01_probes/neutral"
new="PATH/04_fresh2"
ref="PATH/assemblies/Corvus.cornix/genome/v2/v2.5/genome_HC_allpaths41687_v2.5.fasta"

echo $(date)
STARTTIME=$(date +%s)

cd ${dat}

cat *_0based_final.txt | sort -n | uniq > ${1}.txt
awk -v OFS="\t" '{print $1, $2-40, $2}' $1.txt > $1_upstream.bed
awk -v OFS="\t" '{print $1, $2+1, $2+40}' $1.txt > $1_downstream.bed
cp $1.txt ${1}_final.txt 

bedtools getfasta -fi $ref -bed $1_upstream.bed -tab -fo $1_upstream.fa.bed
bedtools getfasta -fi $ref -bed $1_downstream.bed -tab -fo $1_downstream.fa.bed

python probes_generator.py ${dat}/ $1

for i in *_assignedallele.txt
do
fname=${i%.txt*}
awk 'BEGIN{print "scaffold", "0based_pos", "ref", "alt", "assigned"}1' ${fname}.txt > ${fname}_header.txt
done
 
ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/60)) mins to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/3600)) hours to complete this task"


