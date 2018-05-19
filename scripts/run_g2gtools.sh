#!/bin/bash
source activate g2gtools
g2gtools logo 

#
# Note: Put all vcf.gz files at VCF_DIR
#
ARRAY=($(ls ${VCF_DIR:-.}/*.vcf.gz))
EXPANDED=()
for E in "${ARRAY[@]}"; do
    EXPANDED+=("-i ${E}")
done

echo "[g2gtools::vcf2vci] Compiling variant calls..."
g2gtools vcf2vci -o ${SAMPLE_ID}.vci -s ${SAMPLE_ID} --diploid -p ${NUM_CORES} ${EXPANDED[@]}

echo "[g2gtools::patch] Incorporating snps..."
g2gtools patch -i ${REF_SEQ} -c ${SAMPLE_ID}.vci.gz -o ${SAMPLE_ID}.patched.fa -p ${NUM_CORES}

echo "[g2gtools::transform] Incorporating indels..."
g2gtools transform -i ${SAMPLE_ID}.patched.fa -c ${SAMPLE_ID}.vci.gz -o ${SAMPLE_ID}.fa -p ${NUM_CORES}

echo "[g2gtools::convert] Lifting over gtf file..."
g2gtools convert -i ${REF_GTF} -c ${SAMPLE_ID}.vci.gz -o ${SAMPLE_ID}.gtf

echo "[g2gtools::gtf2db] Converting gtf to DB..."
g2gtools gtf2db -i ${SAMPLE_ID}.gtf -o ${SAMPLE_ID}.gtf.db

echo "[g2gtools::extract] Extracting genes..."
g2gtools extract -i ${SAMPLE_ID}.fa -db ${SAMPLE_ID}.gtf.db --genes > ${SAMPLE_ID}.genes.fa

echo "[g2gtools::extract] Extracting transcripts..."
g2gtools extract -i ${SAMPLE_ID}.fa -db ${SAMPLE_ID}.gtf.db --transcripts > ${SAMPLE_ID}.transcripts.fa

echo "[g2gtools::extract] Extracting exons..."
g2gtools extract -i ${SAMPLE_ID}.fa -db ${SAMPLE_ID}.gtf.db --exons > ${SAMPLE_ID}.exons.fa

echo "[g2gtools] Done." 
