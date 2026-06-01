#!/bin/bash
#SBATCH --job-name=germline_calling
#SBATCH --output=/home/data/results/Cognito/kuznetsovads/test/logs/germline_%j.log
#SBATCH --error=/home/data/results/Cognito/kuznetsovads/test/logs/germline_%j.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8


#Поиск зародышевых мутаций для ERR034529 и HCC1143_normal
#ERR034529: BED файл Twist Bioscience 
#HCC1143_normal: BED файл Broad ICE (exome_calling_regions.v1.bed)

set -eo pipefail

GATK_SIF="/home/kuznetsovads/containers/gatk_4.6.2.0.sif"
GATK="singularity exec ${GATK_SIF} gatk"

REF="/home/kuznetsovads/d1/ref_bwa_mem2/hg38.analysisSet.fa"
DBSNP="/home/kuznetsovads/databases/gatk_bundle/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
MILLS="/home/kuznetsovads/databases/gatk_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
G1000="/home/kuznetsovads/databases/gatk_bundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz"

#BED файлы
BED_ERR="/home/kuznetsovads/d1/exome/hg38_Twist_Bioscience_for_Illumina_Exome_2_5_Mito.bed"
BED_HCC="/home/kuznetsovads/test_task/exome_bed/exome_calling_regions.v1.bed"

#BAM
ALIGN_DIR="/home/kuznetsovads/test_task/alignment"
INSTALL_DIR="/home/kuznetsovads/test_task/align_install"

ERR_BAM="${ALIGN_DIR}/aligned_ERR034529.dedup.bam"
# hcc1143_N_clean.bam уже прошёл BQSR — используем напрямую
N_BAM="${INSTALL_DIR}/hcc1143_N_clean.bam"
N_SM="HCC1143_normal"   # SM из @RG заголовка

#dыходные директории
BQSR_DIR="/home/kuznetsovads/test_task/bqsr"
GERM_DIR="/home/kuznetsovads/test_task/germline"
LOG_DIR="/home/data/results/Cognito/kuznetsovads/test/logs"

mkdir -p "${BQSR_DIR}" "${GERM_DIR}" "${LOG_DIR}"

#BQSR только для ERR034529


echo "[$(date '+%H:%M:%S')] === BQSR ERR034529 START ==="

ERR_RECAL="${BQSR_DIR}/ERR034529.recal.table"
ERR_BQSR="${BQSR_DIR}/ERR034529.bqsr.bam"

if [[ ! -f "${ERR_BQSR}" ]]; then
    echo "  BaseRecalibrator: ERR034529..."
    ${GATK} --java-options "-Xmx40g -Djava.io.tmpdir=/tmp" \
        BaseRecalibrator \
        -R "${REF}" \
        -I "${ERR_BAM}" \
        --known-sites "${DBSNP}" \
        --known-sites "${MILLS}" \
        --known-sites "${G1000}" \
        -L "${BED_ERR}" \
        -O "${ERR_RECAL}" \
        2>> "${LOG_DIR}/ERR034529_bqsr.log"

    echo "  ApplyBQSR: ERR034529..."
    ${GATK} --java-options "-Xmx40g -Djava.io.tmpdir=/tmp" \
        ApplyBQSR \
        -R "${REF}" \
        -I "${ERR_BAM}" \
        --bqsr-recal-file "${ERR_RECAL}" \
        -O "${ERR_BQSR}" \
        2>> "${LOG_DIR}/ERR034529_applybqsr.log"
else
    echo "  [SKIP] Already exists: ${ERR_BQSR}"
fi

echo "[$(date '+%H:%M:%S')] === BQSR ERR034529 DONE ==="

#HaplotypeCaller - каждый образец отдельно

echo "[$(date '+%H:%M:%S')] === HAPLOTYPECALLER START ==="

#ERR034529
ERR_GVCF="${GERM_DIR}/ERR034529_twist.g.vcf.gz"
if [[ ! -f "${ERR_GVCF}" ]]; then
    echo "  HaplotypeCaller: ERR034529 (Twist BED)..."
    ${GATK} --java-options "-Xmx40g -Djava.io.tmpdir=/tmp -XX:ParallelGCThreads=4" \
        HaplotypeCaller \
        -R "${REF}" \
        -I "${ERR_BQSR}" \
        -ERC GVCF \
        --dbsnp "${DBSNP}" \
        -L "${BED_ERR}" \
        --interval-padding 100 \
        --native-pair-hmm-threads 8 \
        -O "${ERR_GVCF}" \
        2>> "${LOG_DIR}/ERR034529_twist_haplotypecaller.log"
    echo "  Done: ERR034529 -> ${ERR_GVCF}"
else
    echo "  [SKIP] Already exists: ${ERR_GVCF}"
fi

# HCC1143_normal
N_GVCF="${GERM_DIR}/HCC1143_normal.g.vcf.gz"
if [[ ! -f "${N_GVCF}" ]]; then
    echo "  HaplotypeCaller: HCC1143_normal..."
    ${GATK} --java-options "-Xmx40g -Djava.io.tmpdir=/tmp -XX:ParallelGCThreads=4" \
        HaplotypeCaller \
        -R "${REF}" \
        -I "${N_BAM}" \
        -ERC GVCF \
        --dbsnp "${DBSNP}" \
        -L "${BED_HCC}" \
        --interval-padding 100 \
        -O "${N_GVCF}" \
        2>> "${LOG_DIR}/HCC1143_normal_haplotypecaller.log"
    echo "  Done: HCC1143_normal -> ${N_GVCF}"
else
    echo "  [SKIP] Already exists: ${N_GVCF}"
fi

echo "[$(date '+%H:%M:%S')] === HAPLOTYPECALLER DONE ==="

#GenotypeGVCFs - финальный VCF для каждого образца

echo "[$(date '+%H:%M:%S')] === GENOTYPE GVCFS START ==="

for SAMPLE in "ERR034529_twist" "HCC1143_normal"; do

    if [[ "${SAMPLE}" == "ERR034529_twist" ]]; then
        GVCF="${ERR_GVCF}"
        BED="${BED_ERR}"
    else
        GVCF="${N_GVCF}"
        BED="${BED_HCC}"
    fi

    RAW_VCF="${GERM_DIR}/${SAMPLE}_raw.vcf.gz"

    if [[ -f "${RAW_VCF}" ]]; then
        echo "  [SKIP] Already exists: ${RAW_VCF}"
        continue
    fi

    echo "  GenotypeGVCFs: ${SAMPLE}..."
    ${GATK} --java-options "-Xmx40g -Djava.io.tmpdir=/tmp" \
        GenotypeGVCFs \
        -R "${REF}" \
        -V "${GVCF}" \
        --dbsnp "${DBSNP}" \
        -L "${BED}" \
        -O "${RAW_VCF}" \
        2>> "${LOG_DIR}/${SAMPLE}_genotype.log"

    echo "  Done: ${SAMPLE} -> ${RAW_VCF}"
done

echo "[$(date '+%H:%M:%S')] === GENOTYPE GVCFS DONE ==="

#Hard Filtering

echo "[$(date '+%H:%M:%S')] === FILTERING START ==="

for SAMPLE in "ERR034529_twist" "HCC1143_normal"; do

    RAW_VCF="${GERM_DIR}/${SAMPLE}_raw.vcf.gz"
    SNP_VCF="${GERM_DIR}/${SAMPLE}_snp.vcf.gz"
    SNP_FILT="${GERM_DIR}/${SAMPLE}_snp.filtered.vcf.gz"
    INDEL_VCF="${GERM_DIR}/${SAMPLE}_indel.vcf.gz"
    INDEL_FILT="${GERM_DIR}/${SAMPLE}_indel.filtered.vcf.gz"
    FINAL_VCF="${GERM_DIR}/${SAMPLE}_germline_final.vcf.gz"

    if [[ -f "${FINAL_VCF}" ]]; then
        echo "  [SKIP] Already exists: ${FINAL_VCF}"
        continue
    fi

    echo "  Filtering SNP: ${SAMPLE}..."
    ${GATK} SelectVariants -R "${REF}" -V "${RAW_VCF}" \
        --select-type-to-include SNP -O "${SNP_VCF}" \
        2>> "${LOG_DIR}/${SAMPLE}_filter.log"

    ${GATK} VariantFiltration -R "${REF}" -V "${SNP_VCF}" \
        --filter-expression "QD < 2.0"             --filter-name "QD2"              \
        --filter-expression "QUAL < 30.0"           --filter-name "QUAL30"           \
        --filter-expression "FS > 60.0"             --filter-name "FS60"             \
        --filter-expression "SOR > 3.0"             --filter-name "SOR3"             \
        --filter-expression "MQ < 40.0"             --filter-name "MQ40"             \
        --filter-expression "MQRankSum < -12.5"     --filter-name "MQRankSum-12.5"   \
        --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
        -O "${SNP_FILT}" 2>> "${LOG_DIR}/${SAMPLE}_filter.log"

    echo "  Filtering Indel: ${SAMPLE}..."
    ${GATK} SelectVariants -R "${REF}" -V "${RAW_VCF}" \
        --select-type-to-include INDEL \
        --select-type-to-include MIXED \
        -O "${INDEL_VCF}" 2>> "${LOG_DIR}/${SAMPLE}_filter.log"

    ${GATK} VariantFiltration -R "${REF}" -V "${INDEL_VCF}" \
        --filter-expression "QD < 2.0"              --filter-name "QD2"               \
        --filter-expression "QUAL < 30.0"            --filter-name "QUAL30"            \
        --filter-expression "FS > 200.0"             --filter-name "FS200"             \
        --filter-expression "SOR > 10.0"             --filter-name "SOR10"             \
        --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
        -O "${INDEL_FILT}" 2>> "${LOG_DIR}/${SAMPLE}_filter.log"

    echo "  MergeVcfs: ${SAMPLE}..."
    ${GATK} MergeVcfs \
        -I "${SNP_FILT}" \
        -I "${INDEL_FILT}" \
        -O "${FINAL_VCF}" 2>> "${LOG_DIR}/${SAMPLE}_filter.log"

    echo "  Done: ${SAMPLE} -> ${FINAL_VCF}"
done

echo "[$(date '+%H:%M:%S')] end og filtering"
echo "[$(date '+%H:%M:%S')] the end"


