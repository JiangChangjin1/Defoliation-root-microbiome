#!/bin/sh                                                                       

##########################################################         
#     This BASH script was used to cluster OTUs using 
#     UEASRCH, VSEARCH and Qiime1.
#    
# Changjin Jiang
# 2023.4.9
##########################################################  
#BSUB -J def
#BSUB -n 3
#BSUB -R span[hosts=1]
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -q normal


THREADS=3
REF=reference_data
#PERL=$(which perl)
USEARCH=./tools/usearch11
VSEARCH=./tools/vsearch-2.7.1/bin/vsearch
SEQ_DIR=RawData
OUT_DIR=OutputFolder
FASTQC=./tools/FastQC/FastQC/fastqc
RAW_QUALITY=OriginalQuality
FILTED_QUALITY=QualityAfterTreatment
taxonomy=taxonomy_2023
date


# Process samples  

sample_id=$(ls $SEQ_DIR/*R1.fq |  cut -d  "_"  -f1,2,3,4,5 |cut -d  "/"  -f2 | uniq)
 
for Sample in $sample_id; do

    echo
    echo ====================================
    echo Processing sample $Sample
    echo ====================================
	echo raw reads quality
	chmod  777 $SEQ_DIR/*fq 
	$FASTQC -o $RAW_QUALITY -t 3 $SEQ_DIR/${Sample}_R1.fq $SEQ_DIR/${Sample}_R2.fq
	
    echo Merge paired reads
    $USEARCH  -fastq_mergepairs  $SEQ_DIR/${Sample}_R1.fq \
	-reverse  $SEQ_DIR/${Sample}_R2.fq -fastqout  $OUT_DIR/${Sample}_merged.fq 
	
	echo Removing adaptor
	cutadapt -g AACMGGATTAGATACCCKG \
	-O 19 -e 0.1 -u -29 \
	-o $OUT_DIR/${Sample}_trimmed.fq  $OUT_DIR/${Sample}_merged.fq \
	--minimum-length 150 \
	--discard-untrimmed 
	
	echo
    echo Quality filtering
    $USEARCH -fastq_filter $OUT_DIR/${Sample}_trimmed.fq \
	-fastq_maxee_rate  0.01 \
	-fastaout $OUT_DIR/${Sample}_trimmed_filted.fa \
	-fastqout $FILTED_QUALITY/${Sample}_trimmed_filted.fq \
	-relabel ${Sample}.
	
	echo filted reads quality
    $FASTQC -o $FILTED_QUALITY -t 3 $FILTED_QUALITY/${Sample}_trimmed_filted.fq 
done

echo
echo merge quality result
multiqc  $RAW_QUALITY -o $RAW_QUALITY 
multiqc  $FILTED_QUALITY -o $FILTED_QUALITY  

### open $RAW_QUALITY/multiqc_report.html or $FILTED_QUALITY/multiqc_report.html 
### to check quality


echo
echo ====================================
echo Processing all samples together
echo ====================================

echo
echo Merge all samples

cat $OUT_DIR/*_trimmed_filted.fa > $OUT_DIR/all_filted.fa

#fastp=/public/home/cjjiang/wangfei2022N/all_modify_data/fastp 

echo
echo Dereplicate across samples

$VSEARCH --threads $THREADS \
    --derep_fulllength $OUT_DIR/all_filted.fa \
    --minuniquesize 8 \
    --sizeout \
    --output $OUT_DIR/all_filted_derep.fa

# 2. two piplines were used to obtain OTU table
### unoise3 pipline 
echo
echo Unoise the sequences by unoise3 
# Please see the USEARCH webpage for more information. http://www.drive5.com/usearch/manual/cmd_unoise3.html

$USEARCH -unoise3 $OUT_DIR/all_filted_derep.fa \
    -zotus $OUT_DIR/derep_zotus.fa

echo
echo Reference chimera detection

$VSEARCH  --uchime_ref $OUT_DIR/derep_zotus.fa \
  --db $REF/rdp_gold.fa \
  --nonchimeras $OUT_DIR/nonchimeras_zotus.fa
  
echo
echo generate ZOTU table

$VSEARCH  --usearch_global  $OUT_DIR/all_filted.fa \
  --db  $OUT_DIR/nonchimeras_zotus.fa \
  --id 0.97 \
  --threads 3 \
  --otutabout $OUT_DIR/zotus_tab.txt 
  
biom convert -i $OUT_DIR/zotus_tab.txt \
  -o  $OUT_DIR/zotus_tab.biom \
  --to-hdf5 \
  --table-type="OTU table" 


# 3. qiime1 annotation for unoise3 pipline

###for unoise3 pipline

echo
echo source activate qiime1

source activate qiime1

echo
echo assign taxonomy

assign_taxonomy.py -i  $OUT_DIR/nonchimeras_zotus.fa   \
  -r $REF/silva123_97_otus_16S.fasta \
  -t $REF/silva123_97_raw_taxonomy.txt \
  -o $taxonomy

echo
echo biom processing 

biom add-metadata -i $OUT_DIR/zotus_tab.biom \
    --observation-metadata-fp $taxonomy/nonchimeras_zotus_tax_assignments.txt \
    -o $OUT_DIR/zotus_tab_tax.biom \
    --sc-separated taxonomy --observation-header OTUID,taxonomy  
	
biom convert -i  $OUT_DIR/zotus_tab_tax.biom  \
    -o $OUT_DIR/zotus_tab_tax.tsv \
	--to-tsv \
	--header-key taxonomy 

biom convert -i $OUT_DIR/zotus_tab_tax.tsv \
    -o $OUT_DIR/zotus_tab_tax_final.biom \
	--to-json --table-type="OTU table" \
	--process-obs-metadata taxonomy  


