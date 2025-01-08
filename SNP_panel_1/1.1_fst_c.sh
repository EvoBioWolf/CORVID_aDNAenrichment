### Obtain FSTs of the 3 hybrid zones using Christen's VCF ###

module load vcftools

#directories
scr="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/fst"
refc="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/fst/christen_vcf"
scaff="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/scaffold2chr"

cd ${refc}
#hybrid zone: cornix_P+S vs corone_D
vcftools --vcf genome_crow.gatkHC.vqsrts99.sort5_incl.scaf5.repcont_SNP_ID.recode.vcf --weir-fst-pop G01.txt --weir-fst-pop G02.txt --out vcftools_weir_fst/hybridzone/G01_G02.fst
#hybriz zone: Russia_cornix vs Russia_orientalis
vcftools --vcf genome_crow.gatkHC.vqsrts99.sort5_incl.scaf5.repcont_SNP_ID.recode.vcf --weir-fst-pop G06.txt --weir-fst-pop G07.txt --out vcftools_weir_fst/hybridzone/G06_G07.fst
#hybrid zone: orientalis vs pectoralis
vcftools --vcf genome_crow.gatkHC.vqsrts99.sort5_incl.scaf5.repcont_SNP_ID.recode.vcf --weir-fst-pop G08.txt --weir-fst-pop G09.txt --out vcftools_weir_fst/hybridzone/G08_G09.fst

cd ./vcftools_weir_fst/hybridzone
for fname in *.fst.weir.fst
do
base=${fname%.fst.weir.fst.*}
nl ${base}.fst.weir.fst. | awk '{print $1, $2, $3, $4, $2, LocusName=$2":"$3}' > ${base}.fst.chr
done

cd $scaff
for fname in *.scaffolds
do
base=${fname%.scaffolds*}
for i in $(cat $base.scaffolds)
do 
echo $i $base
awk -v scaff=$i -v chr=$base '{gsub("^"scaff"$", chr, $2)} 1' ${refc}/vcftools_weir_fst/hybridzone/G01_G02.fst.chr > ${refc}/tmp.txt && mv -v ${refc}/tmp.txt ${refc}/vcftools_weir_fst/hybridzone/G01_G02.fst.chr
done
done

for fname in *.scaffolds
do
base=${fname%.scaffolds*}
for i in $(cat $base.scaffolds)
do 
echo $i $base
awk -v scaff=$i -v chr=$base '{gsub("^"scaff"$", chr, $2)} 1' ${refc}/vcftools_weir_fst/hybridzone/G06_G07.fst.chr > ${refc}/tmp2.txt && mv -v ${refc}/tmp2.txt ${refc}/vcftools_weir_fst/hybridzone/G06_G07.fst.chr
done
done

for fname in *.scaffolds
do
base=${fname%.scaffolds*}
for i in $(cat $base.scaffolds)
do 
echo $i $base
awk -v scaff=$i -v chr=$base '{gsub("^"scaff"$", chr, $2)} 1' ${refc}/vcftools_weir_fst/hybridzone/G08_G09.fst.chr > ${refc}/tmp3.txt && mv -v ${refc}/tmp3.txt ${refc}/vcftools_weir_fst/hybridzone/G08_G09.fst.chr
done
done

Rscript $scr/cmplot.R /christen_vcf/vcftools_weir_fst/hybridzone christen G01_G02 G06_G07 G08_G09
