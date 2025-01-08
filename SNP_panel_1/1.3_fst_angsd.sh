### FSTs from ANGSD by Vijay for the 3 hybrid zones ###

#directories
scr="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/fst"
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/reference/2016_Vijay_NatComm_archive/ANGSD"
refa="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/fst/angsd_fst"
scaff="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/scaffold2chr"

cd ${dat}
for fname in *.fst
do
base=${fname%.fst*}
nl ${base}.fst | awk 'BEGIN{print "snp", "CHROM", "POS", "FST", "NA", "CHROM.1", "LocusName"} {print $1, $2, $3, $5, $6, $2, $2":"$3}' > ${refa}/${base}.fst.chr
done

cd $scaff
for fname in *.scaffolds
do
base=${fname%.scaffolds*}
for i in $(cat $base.scaffolds)
do 
echo $i $base
awk -v scaff=$i -v chr=$base '{gsub("^"scaff"$", chr, $2)} 1' ${refa}/G01_G02.fst.chr > ${refa}/tmp.txt && mv -v ${refa}/tmp.txt ${refa}/G01_G02.fst.chr
done
done

for fname in *.scaffolds
do
base=${fname%.scaffolds*}
for i in $(cat $base.scaffolds)
do 
echo $i $base
awk -v scaff=$i -v chr=$base '{gsub("^"scaff"$", chr, $2)} 1' ${refa}/G06_G07.fst.chr > ${refa}/tmp2.txt && mv -v ${refa}/tmp2.txt ${refa}/G06_G07.fst.chr
done
done

for fname in *.scaffolds
do
base=${fname%.scaffolds*}
for i in $(cat $base.scaffolds)
do 
echo $i $base
awk -v scaff=$i -v chr=$base '{gsub("^"scaff"$", chr, $2)} 1' ${refa}/G08_G09.fst.chr > ${refa}/tmp3.txt && mv -v ${refa}/tmp3.txt ${refa}/G08_G09.fst.chr
done
done

Rscript $scr/cmplot.R /angsd_fst angsd G01_G02 G06_G07 G08_G09
