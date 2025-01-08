#!/bin/bash -l
#SBATCH -J bupchk
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00

### sbatch 4.2_backup_check.sh backup ###
### map to outgroups New Caledonian crow and American crow ###
### only for neutral backup snps ###

conda activate renv
module load user_spack
module load plink2/1.9-beta6.10
module load vcftools
module load blast-plus/2.9.0

ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/fst"
ref2="/dss/dsslegfs01/pr53da/pr53da-dss-0018/assemblies/Corvus.cornix/genome/v2/v2.5"
scaff2chr="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/scaffold2chr"
neu="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/snp_panel_summary/03_neutral"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/snp_panel_summary/04_combined"
final="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/snp_panel_summary/05_final_backup"

cd $dat
mkdir $1_snp_panel
cd $1*

## remove escapee polymorphic sites

#recall snps from vcf
vcftools --vcf $ref/christen_vcf/genome_crow.gatkHC.vqsrts99.sort5_incl.scaf5.repcont_SNP_ID.recode.vcf \
--positions $neu/$1_biallelic_pos.txt --recode --out snp_panel_$1
#make sure no non-biallelic sites
vcftools --vcf snp_panel_$1.recode.vcf --freq --out snp_panel_$1
awk 'NR>1' snp_panel_$1.frq | awk '$3==3' 
awk 'NR>1' snp_panel_$1.frq | awk '$3==2' | awk '{print $1, $2}' > snp_panel_$1_trimmed_pos.txt  

## map

#prepare dummy 60bp probes for blasting
awk -v OFS='\t' '{print $1, $2-60, $2+60}' snp_panel_$1_trimmed_pos.txt | awk -v OFS='\t' '{if ($2<1) {$2=0;} print}' > snp_panel_$1_trimmed_pos.bed
bedtools getfasta -fi $ref2/genome_HC_allpaths41687_v2.5.fasta -bed snp_panel_$1_trimmed_pos.bed -fo snp_panel_$1_probesize.fa
bedtools getfasta -fi $ref2/genome_HC_allpaths41687_v2.5.fasta -bed snp_panel_$1_trimmed_pos.bed -fo snp_panel_$1_probesize.fa

#New Caledonian ref genome
blastn -query snp_panel_$1_probesize.fa -db ../bCorMon1 -outfmt "7 qacc sacc evalue qstart qend sstart send" -out bCorMon1_$1_mapped
blastn -query snp_panel_$1_probesize.fa -db ../ASM69197v1 -outfmt "7 qacc sacc evalue qstart qend sstart send" -out ASM69197v1_$1_mapped

for fname in *_mapped
do
base=${fname%_mapped*}
echo ${base}
grep -B 2 "# 0 hits found" ${base}_mapped | grep "# Query:" > ${base}_mapped_nohits 
sed -i -e 's/\# Query: //g' ${base}_mapped_nohits
awk -F ":|-" '{print $1, $3-60}' ${base}_mapped_nohits | awk '{print $1, $2, $1":"$2}' > ${base}_mapped_nohits_snp.pos 
wc -l ${base}_mapped_nohits_snp.pos 
awk 'NR==FNR{c[$1, $2]++;next};c[$3, $4] > 0' ${base}_mapped_nohits_snp.pos $neu/$1_outlier_genic_neutral.txt > ${base}_mapped_nohits_snp.pos.type
done
sort *_mapped_nohits_snp.pos | uniq > bCorMon1_ASM69197v1_$1_mapped_nohits_snp.pos
wc -l bCorMon1_ASM69197v1_$1_mapped_nohits_snp.pos

#generate final positions
awk '{print $1, $2, $3=$1":"$2}' snp_panel_$1_trimmed_pos.txt > snp_panel_$1_trimmed_pos.txt.tmp
awk '{print $3}' bCorMon1_ASM69197v1_$1_mapped_nohits_snp.pos > bCorMon1_ASM69197v1_$1_mapped_nohits_snp.pos.tmp
grep -wiv -f bCorMon1_ASM69197v1_$1_mapped_nohits_snp.pos.tmp snp_panel_$1_trimmed_pos.txt.tmp | awk '{print $1, $2}' > snp_panel_$1_trimmed_mapped_pos.txt
wc -l snp_panel_$1_trimmed_mapped_pos.txt 
rm *.tmp
cp snp_panel_$1_trimmed_mapped_pos.txt $final

## final check 

cd $final

awk 'NR==FNR{c[$1, $2]++;next};c[$3, $4] > 0' snp_panel_$1_trimmed_mapped_pos.txt ../02_genic/all_chrz_fst_dxy99.95_knief_genic_38K_biallelic_only.txt > snp_panel_outlier_genic.txt
#a few neutral snps also overlap with outier/ genic snps

#identify unqiue neutral positions - non-overlapping with outlier / genic
awk '{print $3, $4}' snp_panel_outlier_genic.txt > outlier_genic.pos
grep -wiv -f outlier_genic.pos snp_panel_$1_trimmed_mapped_pos.txt > neutral.pos
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
awk -v scaff=$i -v chr=$base '{gsub("^"scaff"$", chr, $4)} 1' ${final}/snp_panel_outlier_genic_neutral_chr.txt > ${final}/tmp.txt && mv -v ${final}/tmp.txt ${final}/snp_panel_outlier_genic_neutral_chr.txt
done
done

cd $final
#recall snps from vcf
vcftools --vcf $ref/christen_vcf/genome_crow.gatkHC.vqsrts99.sort5_incl.scaf5.repcont_SNP_ID.recode.vcf \
--positions snp_panel_$1_trimmed_mapped_pos.txt --recode --out snp_panel_$1_final

#create overall AF summary and heterozygosity
vcftools --vcf snp_panel_$1_final.recode.vcf --freq --out snp_panel_$1_final
vcftools --vcf snp_panel_$1_final.recode.vcf --het --out snp_panel_$1_final
awk 'NR>1' snp_panel_$1_final.frq | awk -F "\t|:" 'BEGIN{print "Chrom", "Pos", "Ref", "Alt", "RefFreq", "AltFreq", "LocusName"} {print $1, $2, $5, $7, $6, $8, $9=$1":"$2}' > snp_panel_$1_final.frq.summary
awk '{print $1, $2, $3, $4, $7}' OFS='\t' snp_panel_$1_final.frq.summary > snp_panel_$1_final_biocat.txt

#outlier and genic only
vcftools --vcf $ref/christen_vcf/genome_crow.gatkHC.vqsrts99.sort5_incl.scaf5.repcont_SNP_ID.recode.vcf --positions outlier_genic.pos --recode --out snp_panel_outlier_genic_$1_final
vcftools --vcf snp_panel_outlier_genic_$1_final.recode.vcf --freq --out snp_panel_outlier_genic_$1_final
awk 'NR>1' snp_panel_outlier_genic_$1_final.frq | awk -F "\t|:" 'BEGIN{print "Chrom", "Pos", "Ref", "Alt", "LocusName"} {print $1, $2, $5, $7, $1":"$2}' OFS='\t' > snp_panel_outlier_genic_$1_final_biocat.txt

#neutral only
vcftools --vcf $ref/christen_vcf/genome_crow.gatkHC.vqsrts99.sort5_incl.scaf5.repcont_SNP_ID.recode.vcf --positions neutral.pos --recode --out snp_panel_neutral_$1_final
vcftools --vcf snp_panel_neutral_$1_final.recode.vcf --freq --out snp_panel_neutral_$1_final
vcftools --vcf snp_panel_neutral_$1_final.recode.vcf --het --out snp_panel_neutral_$1_final
awk 'NR>1' snp_panel_neutral_$1_final.frq | awk -F "\t|:" 'BEGIN{print "Chrom", "Pos", "Ref", "Alt", "LocusName"} {print $1, $2, $5, $7, $1":"$2}' OFS='\t' > snp_panel_neutral_$1_final_biocat.txt

#check transversion and transition
vcftools --vcf snp_panel_$1_final.recode.vcf --TsTv-summary --out snp_panel_$1_final
vcftools --vcf snp_panel_neutral_$1_final.recode.vcf --TsTv-summary --out snp_panel_neutral_$1_final

#chrplot
awk '{print $1, $2, $2+1, $3, $4, $5}' snp_panel_outlier_genic_neutral_chr.txt | awk -F "\t" 'BEGIN{print "Chrom", "Start", "End", "Group", "ChrNo", "LocusName"}1' > snp_panel_outlier_genic_neutral_chrplot.txt

#fst
vcftools --vcf snp_panel_$1_final.recode.vcf --weir-fst-pop $ref/G01.txt --weir-fst-pop $ref/G02.txt --out G01_G02
vcftools --vcf snp_panel_$1_final.recode.vcf --weir-fst-pop $ref/G06.txt --weir-fst-pop $ref/G07.txt --out G06_G07
vcftools --vcf snp_panel_$1_final.recode.vcf --weir-fst-pop $ref/G08.txt --weir-fst-pop $ref/G09.txt --out G08_G09

for fname in *.weir.fst
do
base=${fname%.weir.fst}
nl ${base}.weir.fst | awk '{print $1, $2, $3, $4, $2}' > ${base}.weir.fst.chr
done

cd ${scaff2chr}
for fname in *.scaffolds
do
base=${fname%.scaffolds*}
for i in $(cat $base.scaffolds)
do 
echo $i $base
awk -v scaff=$i -v chr=$base '{gsub("^"scaff"$", chr, $2)}1' ${final}/G01_G02.weir.fst.chr > ${final}/tmp.txt && mv -v ${final}/tmp.txt ${final}/G01_G02.weir.fst.chr
done
done

for fname in *.scaffolds
do
base=${fname%.scaffolds*}
for i in $(cat $base.scaffolds)
do 
echo $i $base
awk -v scaff=$i -v chr=$base '{gsub("^"scaff"$", chr, $2)}1' ${final}/G06_G07.weir.fst.chr > ${final}/tmp2.txt && mv -v ${final}/tmp2.txt ${final}/G06_G07.weir.fst.chr
done
done

for fname in *.scaffolds
do
base=${fname%.scaffolds*}
for i in $(cat $base.scaffolds)
do 
echo $i $base
awk -v scaff=$i -v chr=$base '{gsub("^"scaff"$", chr, $2)}1' ${final}/G08_G09.weir.fst.chr > ${final}/tmp3.txt && mv -v ${final}/tmp3.txt ${final}/G08_G09.weir.fst.chr
done
done

## plot

cd $final
Rscript ../final_plot.R $1 G01_G02 G06_G07 G08_G09 05_final_backup
rm *.gds*
mkdir plots
mv *.pdf *.jpg ./plots
