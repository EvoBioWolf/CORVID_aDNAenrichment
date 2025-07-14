# Pipeline for bait comparison

All scripts in ./Analyses

# nfcore/Eager2

I ran [Eager](https://nf-co.re/eager/2.4.7) version 2.4.7 on BioHPC with conda profile. The Eager pipeline works with nextflow and is designed for ancient DNA analysis.
I used the pipeline for both whole genome resequencing (wgs) and target-enriched data, but providing a SNP bed file for the latter in addition to the reference genome v2.5 with chrW.

I broke the Eager pipeline into two parts. First raw reads are processed up to the part just before mapping. I mapped the trimmed reads to the reference genome with my own script using BWA-aln to multithread more efficiently. I then resume the Eager pipeline with bam input (sorted but duplicates have not been removed).

---

## Part 1: raw reads processing

Run `1.0.0_eagerTE.sh` 

```
#!/bin/bash -l
#SBATCH -J eager_TE
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=8
#SBATCH --time=14-00:00:00
#SBATCH -o /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/05_aDNA/slurms/slurm-%j-%x.out

# sbatch 1.0.0_eagerTE.sh 00_eager_twist /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta acrow_eager_list_twist
# sbatch 1.0.0_eagerTE.sh 00_eager_mybaits /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta acrow_eager_list_mybaits

echo $(date)
STARTTIME=$(date +%s)

module load charliecloud/0.30
module load nextflow
conda activate biotools

dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/05_aDNA"

cd ${dat}
mkdir ${1}
cd $1

nextflow run nf-core/eager -r 2.4.7 -profile conda --input /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/05_aDNA/${3}.tsv --fasta ${2} \
-c /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/05_aDNA/base_modified2.config \
--snpcapture_bed /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/05_aDNA/probes_232015.bed \
--clip_adapters_list /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/05_aDNA/adapterlist.txt

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
```

`adapterlist.txt` is provided to specify 7 pairs of adapter combinations to be trimmed off from R1 and R2. The read 1 adapter used to make our [ssDNA libraries](https://www.nature.com/articles/s41596-020-0338-0) is 5bp shorter than the default Illumina adapter, so we need to provide this unique adapter pair to AdapterRemoval (shown in the 1st and 2nd rows, while the 3rd and 4th rows show the default pair of adapters). We also included the 6th and 7th pairs of seqeunces as Fastqc flagged over-representation of these unkown sequences in multiple samples after trimming. I have blasted these sequences and they seem to be bacteria or cloning vector.

```
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC  GGAAGAGCGTCGTGTAGGGAAAGAGTGT
AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG  GGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNNAGATCTCGGTGGTCGCCGTATCATT
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC  AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA
AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG  AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNNAGATCTCGGTGGTCGCCGTATCATT
AGATCGGAAGA GGAAGAGCGTCG
ATTCAGCTCCGGTTCCCAACGATCAAGGCGAGTTACATGAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC  TCTTCCGATCTGGAAGAGCGTCGTGTAGGGAAAGAGTGT
ATTCAGCTCCGGTTCCCAACGATCAAGGCGAGTTACATGA    TCATGTAACTCGCCTTGATCGTTGGGAACCGGAGCTG
```

AdapterRemoval will also merge read1 and read2 together to form longer reads. The merged reads are found for each sequenced sample in the folder `./results/adapterremoval/output`. 

---

## Part 2: BWA-aln mapping

Run `1.0.1_mapTE.sh` Seeding is disabled to [optimized mapped reads](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3468387/), but more computation time is required. This is the main reason why I stopped the Eager pipeline to map each sample as a separate job - increasing the number of CPUs on nextflow config file did not improve mapping efficiency

`1.0.1_mapTE.sh`

```
# cd adapterremoval/output
# for i in $(ls *.fq.gz | awk -v FS="_" '{print $1"_"$2}'); do sbatch /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/05_aDNA/1.0.1_mapTE.sh $i 00_eager_twist adapterremoval/output; done
# for i in $(ls *.fq.gz | awk -v FS="_" '{print $2"_"$3}'); do sbatch /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/05_aDNA/1.0.1_mapTE.sh $i 00_eager_mybaits adapterremoval/output; done

echo $(date)
STARTTIME=$(date +%s)

conda activate biotools
dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/05_aDNA"
ref="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/05_aDNA/genome_HC_allpaths41687_v2.5.fasta"

cd ${dat}
cd ${2}/results
mkdir ext_mapped
cd ${3}

echo ${1}

for fname in *${1}*.fq.gz
do
bwa aln -l 1024 -n 0.04 -o 2 -t 8 ${ref} ${fname} > ${dat}/${2}/results/ext_mapped/${1}.sai
done
#default seeding for -l is 32bp   
#-n changed from 0.02 to 0.04 (stricter)

cd ${dat}/${2}/results/ext_mapped
for fname in ${1}.sai
do
bwa samse ${ref} ${1}.sai ${dat}/${2}/results/${3}/*${1}*.fq.gz > ${1}.sam
done

#no quality threshold
samtools view -bh -@ 8 -bS ${1}.sam > ${1}_noqc.bam
samtools sort -@ 8 ${1}_noqc.bam -o ${1}_sorted_noqc.bam
rm ${1}_noqc.bam
```

---

## Part 3: BAM processing

Run `1.0.5_eagerTEbam.sh`

Continuing the Eager pipeline, but with sorted BAM files instead of fastq files this time. Final mutliqc output consolidated to `baitscomparison23.txt`

`1.0.5_eagerTEbam.sh` 
I ran this for 5 sets of probe combination for myBaits and Twist spaarately: (i)probes_232015_SNPsite, (ii)probes_232015, (iii)probes_104k_80bp, (iv)probes_104k, (v)probes_104k_SNPsite

```
# for i in mybaits twist
# do
# sbatch 1.0.5_eagerTEbam.sh 00_eager_${i} /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/05_aDNA/genome_HC_allpaths41687_v2.5.fasta /dss/dsslegfs01/pr53da/pr53da-dss-0018/#projects/2020__ancientDNA/05_aDNA/acrow_eager_list_${i}_bam.tsv probes_232015_SNPsite
# sbatch 1.0.5_eagerTEbam.sh 00_eager_${i} /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/05_aDNA/genome_HC_allpaths41687_v2.5.fasta /dss/dsslegfs01/pr53da/pr53da-dss-0018/#projects/2020__ancientDNA/05_aDNA/acrow_eager_list_${i}_bam.tsv probes_104k_80bp
# sbatch 1.0.5_eagerTEbam.sh 00_eager_${i} /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/05_aDNA/genome_HC_allpaths41687_v2.5.fasta /dss/dsslegfs01/pr53da/pr53da-dss-0018/#projects/2020__ancientDNA/05_aDNA/acrow_eager_list_${i}_bam.tsv probes_104k_SNPsite
# sbatch 1.0.5_eagerTEbam.sh 00_eager_${i} /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/05_aDNA/genome_HC_allpaths41687_v2.5.fasta /dss/dsslegfs01/pr53da/pr53da-dss-0018/#projects/2020__ancientDNA/05_aDNA/acrow_eager_list_${i}_bam.tsv probes_232015
# sbatch 1.0.5_eagerTEbam.sh 00_eager_${i} /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/05_aDNA/genome_HC_allpaths41687_v2.5.fasta /dss/dsslegfs01/pr53da/pr53da-dss-0018/#projects/2020__ancientDNA/05_aDNA/acrow_eager_list_${i}_bam.tsv probes_104k
# done

echo $(date)
STARTTIME=$(date +%s)

module load nextflow
conda activate biotools

dat="/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/05_aDNA"

cd ${dat}
mkdir ${1}
cd $1

#with snpsite bed files and trimmed bam
unset DISPLAY
nextflow run nf-core/eager -profile conda -r 2.4.7 --input ${3} --fasta ${2} \
-c /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/05_aDNA/base_modified2.config \
--snpcapture_bed /dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/05_aDNA/${4}.bed \
--mtnucratio_header chrM --run_mtnucratio \
--run_trim_bam --bamutils_clip_single_stranded_none_udg_left 5 --bamutils_clip_single_stranded_none_udg_right 5 
```
---

## Variant calling

Run `1.1.1_angsd_baitscomparison.sh`

```
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
angsd -bam bam_${1}_4.filelist -GL 2 -doMaf 2 -doMajorMinor 1 -ref ${ref} -doGeno 4 -doPost 1 -doGlf 2 -doCounts 1 -doPlink 2 \
-minMapQ 20 -minQ 20 -uniqueOnly 1 -remove_bads 1 -geno_minDepth 3 \
-nThreads 12 -rf ${dat}/${3} -out ${1}_${2}_geno_maxmis_q20_dp3

#parameters controlling missingness allowed for each sample and site are removed

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
```

## Comparison
comparison done on R. See `baitscomparison23.R`