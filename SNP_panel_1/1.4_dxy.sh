### obtain persite dxy of G01vsG02 (only the main hybrid zone considered here) using Simon Martin's script ###

conda activate renv
module load bcftools/1.9

ref="dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/fst"
dat="refdss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/dxy"
scaff="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/scaffold2chr"

cd $dat

#preparing readable input for Simon Martin's script
#script will not work unless it's biallelic!
bcftools query -f '%CHROM\t%POS[\t%TGT]\n' $ref/christen_vcf/genome_crow.gatkHC.vqsrts99.sort5_incl.scaf5.repcont_SNP_ID.recode.vcf \
> genome_crow.gatkHC.vqsrts99.sort5_incl.scaf5.repcont_SNP_ID_christen.recode.geno
cat header.txt genome_crow.gatkHC.vqsrts99.sort5_incl.scaf5.repcont_SNP_ID_christen.recode.geno > genome_crow.gatkHC.vqsrts99.sort5_incl.scaf5.repcont_SNP_ID_christen.recode.final.geno 
awk '{gsub(/\.\/\./, "N\/N")}' genome_crow.gatkHC.vqsrts99.sort5_incl.scaf5.repcont_SNP_ID_christen.recode.final.geno > genome_crow.gatkHC.vqsrts99.sort5_incl.scaf5.repcont_SNP_ID_christen.recode.final.final.geno
#keeping only diploid sites
for i in $(cat ${scaff}/chrZ.scaffolds)
do
grep -v $i genome_crow.gatkHC.vqsrts99.sort5_incl.scaf5.repcont_SNP_ID_christen.recode.final.final.geno > temp mv temp genome_crow.gatkHC.vqsrts99.sort5_incl.scaf5.repcont_SNP_ID_christen.recode.final.final.geno
done 
for i in $(cat ${scaff}/allchrUN.txt)
do
grep -v $i genome_crow.gatkHC.vqsrts99.sort5_incl.scaf5.repcont_SNP_ID_christen.recode.final.final.geno > temp2 && mv temp2 genome_crow.gatkHC.vqsrts99.sort5_incl.scaf5.repcont_SNP_ID_christen.recode.final.final.geno
done
#remove polymorphic sites
awk 'length($NF)==3 && length($(NF-117))==3' genome_crow.gatkHC.vqsrts99.sort5_incl.scaf5.repcont_SNP_ID_christen.recode.final.final.geno > genome_crow.gatkHC.vqsrts99.sort5_incl.scaf5.repcont_SNP_ID_christen.recode.final.final.final.geno

rm genome_crow.gatkHC.vqsrts99.sort5_incl.scaf5.repcont_SNP_ID_christen.recode.final.geno genome_crow.gatkHC.vqsrts99.sort5_incl.scaf5.repcont_SNP_ID_christen.recode.final.final.geno 

#Simon Martin's script 
python genomics_general/popgenWindows.py --popsFile pop115.txt -g genome_crow.gatkHC.vqsrts99.sort5_incl.scaf5.repcont_SNP_ID_christen.recode.final.final.final.geno \
-o dxypersite_G01.G02_christen.csv \
-f phased -w 1 -p G01 -p G02 -T 12

awk -F ',' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $1}' dxypersite_G01.G02_christen.csv > dxypersite_G01.G02_christen.csv.chr

cd $scaff
for fname in *.scaffolds
do
base=${fname%.scaffolds*}
for i in $(cat $base.scaffolds)
do
awk -v scaff=$i -v chr=$base '{gsub("^"scaff"$", chr, $1)} 1' ${dat}/dxypersite_G01.G02_christen.csv.chr > ${dat}/tmp1.txt && mv -v ${dat}/tmp1.txt ${dat}/dxypersite_G01.G02_christen.csv.chr
done
done
