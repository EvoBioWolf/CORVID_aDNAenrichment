#!/bin/bash -l
#SBATCH -J syn2
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_highmem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --time=1-00:00:00

### extract 4fold, intronic and intergenic sites per population ###
# the reference and annotations (GCF_000738735.1) were downloaded directly from NCBI instead of the one stored # in the assembly folder
# they should correspond to v2.5.

module load vcftools
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/neutral"
ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/assemblies/Corvus.cornix/genome/v2/v2.5"
refc="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/fst/christen_vcf"
scaff="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/scaffold2chr"

cd $dat
python ./scripts_rwilliamson/NewAnnotateRef.py GCF_000738735.1_Hooded_Crow_genome_genomic.fna GCF_000738735.1_Hooded_Crow_genome_genomic.gff > genome_HC_allpaths41687_v2.5_annotated.sites

#convert scaffold names of GCF_000738735.1.gff to the same as our reference
awk 'NR==FNR{A[$2]=$1;next}A[$1]{sub($1,A[$1]);print}' ncbi_scafftransform.txt GCF_000738735.1_Hooded_Crow_genome_genomic.gff > GCF_000738735.1_Hooded_Crow_genome_genomic.transformed.gff 

#Generate population VCF
cd ./neutral_perpop
for i in $(cat panel)
do
vcftools --vcf ${refc}/genome_crow.gatkHC.vqsrts99.sort5_incl.scaf5.repcont_SNP_ID.recode.vcf \
--keep ../poplist/$i.txt --recode --mac 1 -c | gzip -c > $i.vcf.gz
done

#Generate a vcf summary using vcfSummarizer.py by Robert Williamson (https://github.com/fabbyrob/science/blob/master/pileup_analyzers/vcfSummarizer.py)
cd ../neutral_perpop
for fname in *.vcf.gz
do
base=${fname%.vcf.gz*}
python ../scripts_rwilliamson/vcfSummarizer.py ./${base}.vcf.gz \
$dat/genome_HC_allpaths41687_v2.5_annotated.sites -q 0 -d 0 -D 1000 -L 0 -N 0 > ${base}.summary
done

#identify synonymous sites on each pop vcf (christen)
for fname in *.summary
do
base=${fname%.summary*}
awk 'NR<=12 || $8==3' ${base}.summary >> ./summary/${base}_0fold.txt
awk 'NR<=12 || $8==4' ${base}.summary >> ./summary/${base}_4fold.txt
awk 'NR<=12 || $8==0' ${base}.summary >> ./summary/${base}_intergene.txt
awk 'NR<=12 || $8==1' ${base}.summary >> ./summary/${base}_intron.txt
awk 'NR<=12 || $8==2' ${base}.summary >> ./summary/${base}_exon.txt
done

#print SAF of each pop and remove chrZ
cd ./summary
for fname in *.txt
do
base=${fname%.txt*}
awk 'NR>12' ${base}.txt | \
awk '{if ($5 + 0 !=0) $10 = $5 / $7; else $10=0} {if ($6 + 0 !=0) $11 = $6 / $7; else $11=0}1' | \
awk 'BEGIN{print "Scaffold", "Pos", "Ref", "Alt", "RefCount", "AltCount", "TotalCount", "RefFreq", "AltFreq", "LocusName"} {print $1, $2, $3, $4, $5, $6, $7, $10, $11, $12=$1":"$2}' | \
grep -v -f ${scaff}/chrZ.scaffolds > ${base}_AF_nochrz.txt
done

mv *0fold*nochrz.txt ./0fold
mv *4fold*nochrz.txt ./4fold
mv *intergene*nochrz.txt ./intergene
mv *intron*nochrz.txt ./intron
mv *exon*nochrz.txt ./exon
mv *.txt ./original

#print SNP positions of intronic & intergenic region of each pop
cd $dat/neutral_perpop/summary/intron
for fname in *_nochrz.txt
do
base=${fname%_intron_AF_nochrz.txt*} 
awk 'NR>1 {print $1, $2}' ${base}_intron_AF_nochrz.txt > ${base}_pos.txt
done
cd $dat/neutral_perpop/summary/intergene  
for fname in *_nochrz.txt
do
base=${fname%_intergene_AF_nochrz.txt*} 
awk 'NR>1 {print $1, $2}' ${base}_intergene_AF_nochrz.txt > ${base}_pos.txt
done

#extract intronic & intergenic sites from the VCF of each pop
cd $dat/neutral_perpop
for i in $(cat panelpop)
do
echo $i
vcftools --gzvcf genome_crow.gatkHC.vqsrts99.sort5_incl.scaf5.repcont_SNP_ID_$i.vcf.gz \
--positions ./summary/intron/$i_pos.txt --recode --out ./summary/intron/$i
vcftools --gzvcf genome_crow.gatkHC.vqsrts99.sort5_incl.scaf5.repcont_SNP_ID_$i.vcf.gz \
--positions ./summary/intergene/$i_pos.txt --recode --out ./summary/intergene/$i
done
