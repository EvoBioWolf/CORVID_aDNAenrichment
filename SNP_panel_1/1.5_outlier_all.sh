### compile outlier fst and dxy from different datasets and hybrid zones into a single final output ###
### sbatch 1.5_outlier_all.sh ###

scr="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/fst"

#compile fsts from three different datasets into a single file for each hybrid zone
Rscript $scr/overlap_per_hybridzone.R G01_G02 G06_G07 G08_G09 chr_autosomes chr_chrZ

#replace dummy LocusName with corresponding Fst value of each dataset
cd /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/01_probes/fst/fst_summary
awk 'NR==FNR{A[$6]=$4;next}FNR==1{print $0; next}{if($3=A[$1])sub($3,A[$1],$3);else $3="NA";print}' fst_christen_G01_G02_chr_autosomes_only.txt G01_G02_chr_autosomes_persite_fst99.9.txt > tmp.txt && mv -v tmp.txt G01_G02_chr_autosomes_persite_fst99.9.txt
awk 'NR==FNR{A[$6]=$4;next}FNR==1{print $0; next}{if($4=A[$1])sub($4,A[$1],$4);else $4="NA";print}' fst_verena_G01_G02_chr_autosomes_only.txt G01_G02_chr_autosomes_persite_fst99.9.txt > tmp.txt && mv -v tmp.txt G01_G02_chr_autosomes_persite_fst99.9.txt
awk 'NR==FNR{A[$8]=$5;next}FNR==1{print $0; next}{if($5=A[$1])sub($5,A[$1],$5);else $5="NA";print}' fst_angsd_G01_G02_chr_autosomes_only.txt G01_G02_chr_autosomes_persite_fst99.9.txt > tmp.txt && mv -v tmp.txt G01_G02_chr_autosomes_persite_fst99.9.txt

awk 'NR==FNR{A[$6]=$4;next}FNR==1{print $0; next}{if($3=A[$1])sub($3,A[$1],$3);else $3="NA";print}' fst_christen_G06_G07_chr_autosomes_only.txt G06_G07_chr_autosomes_persite_fst99.9.txt > tmp.txt && mv -v tmp.txt G06_G07_chr_autosomes_persite_fst99.9.txt
awk 'NR==FNR{A[$6]=$4;next}FNR==1{print $0; next}{if($4=A[$1])sub($4,A[$1],$4);else $4="NA";print}' fst_verena_G06_G07_chr_autosomes_only.txt G06_G07_chr_autosomes_persite_fst99.9.txt > tmp.txt && mv -v tmp.txt G06_G07_chr_autosomes_persite_fst99.9.txt
awk 'NR==FNR{A[$8]=$5;next}FNR==1{print $0; next}{if($5=A[$1])sub($5,A[$1],$5);else $5="NA";print}' fst_angsd_G06_G07_chr_autosomes_only.txt G06_G07_chr_autosomes_persite_fst99.9.txt > tmp.txt && mv -v tmp.txt G06_G07_chr_autosomes_persite_fst99.9.txt

awk 'NR==FNR{A[$6]=$4;next}FNR==1{print $0; next}{if($3=A[$1])sub($3,A[$1],$3);else $3="NA";print}' fst_christen_G08_G09_chr_autosomes_only.txt G08_G09_chr_autosomes_persite_fst99.9.txt > tmp.txt && mv -v tmp.txt G08_G09_chr_autosomes_persite_fst99.9.txt
awk 'NR==FNR{A[$6]=$4;next}FNR==1{print $0; next}{if($4=A[$1])sub($4,A[$1],$4);else $4="NA";print}' fst_verena_G08_G09_chr_autosomes_only.txt G08_G09_chr_autosomes_persite_fst99.9.txt > tmp.txt && mv -v tmp.txt G08_G09_chr_autosomes_persite_fst99.9.txt
awk 'NR==FNR{A[$8]=$5;next}FNR==1{print $0; next}{if($5=A[$1])sub($5,A[$1],$5);else $5="NA";print}' fst_angsd_G08_G09_chr_autosomes_only.txt G08_G09_chr_autosomes_persite_fst99.9.txt > tmp.txt && mv -v tmp.txt G08_G09_chr_autosomes_persite_fst99.9.txt

for fname in *_chr_chrz_persite_fst99.9.txt
do
base=${fname%_chr_chrz_persite_fst99.9.txt*} 
awk '{print $1}' ${base}_chr_chrz_persite_fst99.9.txt > ${base}_chrz.pos
done
sort *_chrz.pos | uniq > allhybridzones_chrz_fst99.9.txt
wc -l allhybridzones_chrz_fst99.9.txt

#compile all outliers into a single file and save in snp_panel_summary
Rscript $scr/overlap_all_hybridzones.R G01_G02 G06_G07 G08_G09
