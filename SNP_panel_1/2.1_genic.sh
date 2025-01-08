### extract variants between G01&G02 in selected genes involving in melanogenesis pathway ###
### merge with outlier region then filter out polymorphic sites ###
# the reference and annotations (GCF_000738735.1) were downloaded directly from NCBI instead of the one stored in the assembly 
# they should correspond to v2.5, but the scaffold names are different     
                                             

module load bedtools2/2.27.1
module load vcftools/0.1.14

dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/genic"
ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/fst"
des="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/scaffold2chr"
pan="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/snp_panel_summary"

cd ${dat}
for i in $(cat genelist_group1.txt)
do
grep -w $i /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/neutral/GCF_000738735.1_Hooded_Crow_genome_genomic.transformed.gff > ${i}_group1gene
done
cat *_group1gene > group1gene.GFF
mv *_group1gene ./group1genes

vcftools --vcf ${ref}/christen_vcf/genome_crow.gatkHC.vqsrts99.sort5_incl.scaf5.repcont_SNP_ID.recode.vcf \ 
--keep G01G02.txt --recode --mac 1 \
--out genome_crow.gatkHC.vqsrts99.sort5_incl.scaf5.repcont_SNP_ID_G01G02

bedtools intersect -a genome_crow.gatkHC.vqsrts99.sort5_incl.scaf5.repcont_SNP_ID_G01G02.recode.vcf -b group1gene.GFF > genic_variants_group1genes_G01G02.vcf
sort -u -k 3,3 genic_variants_group1genes_G01G02.vcf > genic_variants_group1genes_unique_G01G02.vcf

cd ${des}
for fname in *.scaffolds
do
base=${fname%.scaffolds*}
for i in $(cat $base.scaffolds)
do
echo $i $base
awk -v scaff=$i -v chr=$base '{gsub("^"scaff"$", chr, $2)} 1' ${dat}/genic_variants_group1genes_LocusName_G01G02.txt > tmp1.txt && mv -v tmp1.txt ${dat}/genic_variants_group1genes_LocusName_G01G02.txt
done
done

cp ${dat}/genic_variants_group1genes_LocusName_G01G02.txt ${pan}/02_genic/group1genes_G01G02_14k.txt

## merge outlier and genic regions
cd $pan
Rscript outlier_genic.R

## bialleic check
cd $pan/02*
awk 'NR>1 {print $3, $4}' all_chrz_fst_dxy99.95_knief_genic_46k.txt > outlier_genic_pos.txt

vcftools --vcf $ref/christen_vcf/genome_crow.gatkHC.vqsrts99.sort5_incl.scaf5.repcont_SNP_ID.recode.vcf \
--positions outlier_genic_pos.txt --recode --out outlier_genic
vcftools --vcf outlier_genic.recode.vcf --min-alleles 2 --max-alleles 2 --recode --out outlier_genic_biallelic
vcftools --vcf outlier_genic_biallelic.recode.vcf --freq --out outlier_genic_biallelic

awk 'NR>1 {print $1, $2}' outlier_genic_biallelic.frq > outlier_genic_biallelic_pos.txt
awk 'NR==FNR{c[$1$2]++;next};c[$3$4] > 0' outlier_genic_biallelic_pos.txt all_chrz_fst_dxy99.95_knief_genic_46k.txt > all_chrz_fst_dxy99.95_knief_genic_38K_biallelic_only.txt
wc -l all_chrz_fst_dxy99.95_knief_genic_38K_biallelic_only.txt #38885

#genic_biallele
awk 'NR==FNR{c[$3$4]++;next};c[$4$3] > 0' all_chrz_fst_dxy99.95_knief_genic_38K_biallelic_only.txt group1genes_G01G02_14k.txt | \
awk 'BEGIN{print "LocusName", "CHROM", "Pos", "Scaffold"}'1 > group1genes_G01G02_14k_biallelic_only.txt
wc -l group1genes_G01G02_14k_biallelic_only.txt #14197 (dropped from 14280)

#outlier_biallele
cd $pan/01*
awk 'NR==FNR{c[$3$4]++;next};c[$3$4] > 0' $dat/02_genic/all_chrz_fst_dxy99.95_knief_genic_38K_biallelic_only.txt all_chrz_fst_dxy99.95_knief_strict_31k.txt | \
awk 'BEGIN{print "LocusName", "OverlapType", "Scaffold", "Pos"}1' > all_chrz_fst_dxy99.95_knief_strict_25k_biallelic.txt
wc -l all_chrz_fst_dxy99.95_knief_strict_25k_biallelic.txt #24985 (reduced from 31615)

