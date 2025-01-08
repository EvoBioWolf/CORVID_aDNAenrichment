#!/bin/bash -l
#SBATCH -J neutral
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=4
#SBATCH --time=1-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/snp_panel_2/slurms/slurm-%j-%x.out
#SBATCH --mem-per-cpu=4763mb

# sbatch 2.1_neutral_outgroup.sh 9pop_overlapped_filtered_norepeats_AMcrow_biallele 8pop_noamcrow outgroup_ascertained4

conda activate biotools
#conda activate py2
module load vcftools/0.1.14-gcc8
#picard.jar=2.25.7

ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/snp_panel_2"
neu="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/neutral"
frh="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/06_results/neutral"

echo $(date)
STARTTIME=$(date +%s)

cd ${dat}/${3}

vcftools --gzvcf ${1}.vcf.gz --keep ${dat}/poplist/${2}.pop --mac 1 --recode --stdout | gzip -c > $2.vcf.gz
zcat $2.vcf.gz | grep -v "INFO=<ID=\|FORMAT=<ID" | gzip -c > $2_edited.vcf.gz
rm $2.vcf.gz
mv $2_edited.vcf.gz $2.vcf.gz 

#for comparison: including invariant sites of AMcrow#
vcftools --gzvcf ${1}.vcf.gz --keep ${dat}/poplist/amcrow.pop --recode --stdout | gzip -c > $2.vcf.gz
zcat $2.vcf.gz | grep -v "INFO=<ID=\|FORMAT=<ID" | gzip -c > $2_edited.vcf.gz
rm $2.vcf.gz
mv $2_edited.vcf.gz $2.vcf.gz

#sort 
java -Xmx20g -Djava.io.tmpdir=$dat -jar /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/picard.jar SortVcf TMP_DIR=$dat \
      MAX_RECORDS_IN_RAM=10000 \
      I=$2.vcf.gz \
      O=$2_sorted.vcf
bgzip ${2}_sorted.vcf

conda activate py2
python ${neu}/scripts_rwilliamson/vcfSummarizer.py $2_sorted.vcf.gz ${neu}/genome_HC_allpaths41687_v2.5_annotated.sites -q 0 -d 0 -D 1000 -L 0 -N 0 > $2.summary

mkdir 0fold 4fold intergene intron exon

#identify synonymous sites on each pop vcf (christen)
awk 'NR<=12 || $8==3' $2.summary > ./0fold/${2}_0fold.txt
awk 'NR<=12 || $8==4' $2.summary > ./4fold/${2}_4fold.txt
awk 'NR<=12 || $8==0' $2.summary > ./intergene/${2}_intergene.txt
awk 'NR<=12 || $8==1' $2.summary > ./intron/${2}_intron.txt
awk 'NR<=12 || $8==2' $2.summary > ./exon/${2}_exon.txt

cd ${dat}/${3}/intron
awk 'NR>12 {print $1, $2}' ${2}_intron.txt > ${2}.pos
echo $(wc -l ${2}.pos)
cp ${dat}/outgroup_ascertained2/intron/amcrow.pos ./
cp ${dat}/outgroup_ascertained2/intron/amcrow_incinvariant.pos ./
cp ${frh}/intron/cnx6.pos ./
cp ${frh}/intron/cor1.pos ./
cp ${frh}/intron/ori3.pos ./

cd ${dat}/${3}/intergene
awk 'NR>12 {print $1, $2}' ${2}_intergene.txt > ${2}.pos
echo $(wc -l ${2}.pos)
cp ${dat}//outgroup_ascertained2/intergene/amcrow.pos ./
cp ${dat}//outgroup_ascertained2/intergene/amcrow_incinvariant.pos ./
cp ${frh}/intergene/cnx6.pos ./
cp ${frh}/intergene/cor1.pos ./
cp ${frh}/intergene/ori3.pos ./

cd ${dat}/${3}/4fold
awk 'NR>12 {print $1, $2}' ${2}_4fold.txt > ${2}.pos
echo $(wc -l ${2}.pos)
cp ${dat}//outgroup_ascertained2/4fold/amcrow.pos ./
cp ${dat}//outgroup_ascertained2/4fold/amcrow_incinvariant.pos ./
cp ${frh}/4fold/cnx6.pos ./
cp ${frh}/4fold/cor1.pos ./
cp ${frh}/4fold/ori3.pos ./

#once all done, do once to extract only sites heterozygous in AM and any of the combinedpop
cd ${dat}/${3}/intron
awk 'NR==FNR{c[$1, $2]++;next};c[$1, $2] > 0' ${2}.pos amcrow.pos > overlap_${2}_amcrow_pos.txt
awk 'NR==FNR{c[$1, $2]++;next};c[$1, $2] > 0' ${2}.pos amcrow_incinvariant.pos > overlap_${2}_amcrow_incinvariant_pos.txt
awk 'NR==FNR{c[$1, $2]++;next};c[$1, $2] > 0' ${2}.pos cor1.pos > overlap_${2}_cor1.tmp
awk 'NR==FNR{c[$1, $2]++;next};c[$1, $2] > 0' ${2}.pos cnx6.pos > overlap_${2}_cnx6.tmp
awk 'NR==FNR{c[$1, $2]++;next};c[$1, $2] > 0' ${2}.pos ori3.pos > overlap_${2}_ori3.tmp
cat overlap_${2}_cor1.tmp overlap_${2}_cnx6.tmp overlap_${2}_ori3.tmp | sort -n | uniq > overlap_${2}_pos.txt
echo $(wc -l overlap_${2}_pos.txt)

cd ${dat}/${3}/intergene
awk 'NR==FNR{c[$1, $2]++;next};c[$1, $2] > 0' ${2}.pos amcrow.pos > overlap_${2}_amcrow_pos.txt
awk 'NR==FNR{c[$1, $2]++;next};c[$1, $2] > 0' ${2}.pos amcrow_incinvariant.pos > overlap_${2}_amcrow_incinvariant_pos.txt
awk 'NR==FNR{c[$1, $2]++;next};c[$1, $2] > 0' ${2}.pos cor1.pos > overlap_${2}_cor1.tmp
awk 'NR==FNR{c[$1, $2]++;next};c[$1, $2] > 0' ${2}.pos cnx6.pos > overlap_${2}_cnx6.tmp
awk 'NR==FNR{c[$1, $2]++;next};c[$1, $2] > 0' ${2}.pos ori3.pos > overlap_${2}_ori3.tmp
cat overlap_${2}_cor1.tmp overlap_${2}_cnx6.tmp overlap_${2}_ori3.tmp | sort -n | uniq > overlap_${2}_pos.txt
echo $(wc -l overlap_${2}_pos.txt)

cd ${dat}/${3}/4fold   
awk 'NR==FNR{c[$1, $2]++;next};c[$1, $2] > 0' ${2}.pos amcrow.pos > overlap_${2}_amcrow_pos.txt
awk 'NR==FNR{c[$1, $2]++;next};c[$1, $2] > 0' ${2}.pos amcrow_incinvariant.pos > overlap_${2}_amcrow_incinvariant_pos.txt
awk 'NR==FNR{c[$1, $2]++;next};c[$1, $2] > 0' ${2}.pos cor1.pos > overlap_${2}_cor1.tmp
awk 'NR==FNR{c[$1, $2]++;next};c[$1, $2] > 0' ${2}.pos cnx6.pos > overlap_${2}_cnx6.tmp
awk 'NR==FNR{c[$1, $2]++;next};c[$1, $2] > 0' ${2}.pos ori3.pos > overlap_${2}_ori3.tmp
cat overlap_${2}_cor1.tmp overlap_${2}_cnx6.tmp overlap_${2}_ori3.tmp | sort -n | uniq > overlap_${2}_pos.txt
echo $(wc -l overlap_${2}_pos.txt)

cd ${dat}/${3}
cat intron/overlap_${2}_amcrow_pos.txt intergene/overlap_${2}_amcrow_pos.txt 4fold/overlap_${2}_amcrow_pos.txt | sort | uniq > neuall_outgroup_ascertained.pos
cat intron/overlap_${2}_amcrow_incinvariant_pos.txt intergene/overlap_${2}_amcrow_incinvariant_pos.txt 4fold/overlap_${2}_amcrow_incinvariant_pos.txt | sort | uniq > neuall_outgroup_NOTascertained.pos

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/60)) mins to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/3600)) hours to complete this task"


