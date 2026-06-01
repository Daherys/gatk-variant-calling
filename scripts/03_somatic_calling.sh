#!/bin/bash
#SBATCH --job-name=somatic_calling
#SBATCH --output=/home/data/results/Cognito/kuznetsovads/test/logs/somatic_%j.log
#SBATCH --error=/home/data/results/Cognito/kuznetsovads/test/logs/somatic_%j.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8


#gоиск соматических мутаций: HCC1143_tumor vs HCC1143_normal

set -eo pipefail


#пути


GATK_SIF="/home/kuznetsovads/containers/gatk_4.6.2.0.sif"
GATK="singularity exec ${GATK_SIF} gatk"

REF="/home/kuznetsovads/d1/ref_bwa_mem2/hg38.analysisSet.fa"
DBSNP="/home/kuznetsovads/databases/gatk_bundle/Homo_sapiens_assembly38.dbsnp138.vcf.gz"

#BED файл Broad ICE exome kit
BED="/home/kuznetsovads/test_task/exome_bed/exome_calling_regions.v1.bed"

#входные BAM (уже предобработаны)
INSTALL_DIR="/home/kuznetsovads/test_task/align_install"
T_BAM="${INSTALL_DIR}/hcc1143_T_clean.bam"
N_BAM="${INSTALL_DIR}/hcc1143_N_clean.bam"

#имена сэмплов из @RG заголовков
T_SM="HCC1143_tumor"
N_SM="HCC1143_normal"

#ресурсы для Mutect2

EXAC="/home/kuznetsovads/databases/gatk_bundle/small_exac_common_3.hg38.vcf.gz"

#выходные директории
SOM_DIR="/home/kuznetsovads/test_task/somatic"
LOG_DIR="/home/data/results/Cognito/kuznetsovads/test/logs"

mkdir -p "${SOM_DIR}" "${LOG_DIR}"

#Mutect2 — tumor vs normal

echo "[$(date '+%H:%M:%S')] === MUTECT2 START ==="

RAW_VCF="${SOM_DIR}/somatic_raw.vcf.gz"
F1R2="${SOM_DIR}/f1r2.tar.gz"

if [[ ! -f "${RAW_VCF}" ]]; then
    echo "  Mutect2: ${T_SM} vs ${N_SM}..."
    ${GATK} --java-options "-Xmx40g -Djava.io.tmpdir=/tmp" \
        Mutect2 \
        -R "${REF}" \
        -I "${T_BAM}" \
        -I "${N_BAM}" \
        --tumor-sample "${T_SM}" \
        --normal-sample "${N_SM}" \
        --germline-resource "${DBSNP}" \
        -L "${BED}" \
        --interval-padding 100 \
        --f1r2-tar-gz "${F1R2}" \
        -O "${RAW_VCF}" \
        2>> "${LOG_DIR}/mutect2.log"
    echo "  Done: ${RAW_VCF}"
else
    echo "  [SKIP] Already exists: ${RAW_VCF}"
fi

echo "[$(date '+%H:%M:%S')] === MUTECT2 DONE ==="


#Оценка артефактов ориентации ридов (F1R2)
# Учитывает FFPE-подобные артефакты C->T/G->A

echo "[$(date '+%H:%M:%S')] === ORIENTATION MODEL START ==="

OB_MODEL="${SOM_DIR}/read_orientation_model.tar.gz"

if [[ ! -f "${OB_MODEL}" ]]; then
    echo "  LearnReadOrientationModel..."
    ${GATK} --java-options "-Xmx16g -Djava.io.tmpdir=/tmp" \
        LearnReadOrientationModel \
        -I "${F1R2}" \
        -O "${OB_MODEL}" \
        2>> "${LOG_DIR}/orientation_model.log"
    echo "  Done: ${OB_MODEL}"
else
    echo "  [SKIP] Already exists: ${OB_MODEL}"
fi

echo "[$(date '+%H:%M:%S')] === ORIENTATION MODEL DONE ==="


#Оценка контаминации

echo "[$(date '+%H:%M:%S')] === CONTAMINATION ESTIMATION START ==="

T_PILEUP="${SOM_DIR}/tumor_pileups.table"
N_PILEUP="${SOM_DIR}/normal_pileups.table"
CONTAM="${SOM_DIR}/contamination.table"
SEGMENTS="${SOM_DIR}/tumor_segments.table"

#Проверяем наличие ExAC файла
if [[ ! -f "${EXAC}" ]]; then
    echo "  [WARN] small_exac_common not found: ${EXAC}"
    echo "  Downloading small_exac_common..."
    wget -q https://storage.googleapis.com/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz \
        -O "${EXAC}"
    wget -q https://storage.googleapis.com/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz.tbi \
        -O "${EXAC}.tbi"
fi

if [[ ! -f "${T_PILEUP}" ]]; then
    echo "  GetPileupSummaries: tumor..."
    ${GATK} --java-options "-Xmx16g -Djava.io.tmpdir=/tmp" \
        GetPileupSummaries \
        -I "${T_BAM}" \
        -V "${EXAC}" \
        -L "${EXAC}" \
        -O "${T_PILEUP}" \
        2>> "${LOG_DIR}/pileup_tumor.log"
fi

if [[ ! -f "${N_PILEUP}" ]]; then
    echo "  GetPileupSummaries: normal..."
    ${GATK} --java-options "-Xmx16g -Djava.io.tmpdir=/tmp" \
        GetPileupSummaries \
        -I "${N_BAM}" \
        -V "${EXAC}" \
        -L "${EXAC}" \
        -O "${N_PILEUP}" \
        2>> "${LOG_DIR}/pileup_normal.log"
fi

if [[ ! -f "${CONTAM}" ]]; then
    echo "  CalculateContamination..."
    ${GATK} --java-options "-Xmx16g -Djava.io.tmpdir=/tmp" \
        CalculateContamination \
        -I "${T_PILEUP}" \
        --matched-normal "${N_PILEUP}" \
        --tumor-segmentation "${SEGMENTS}" \
        -O "${CONTAM}" \
        2>> "${LOG_DIR}/contamination.log"
fi

echo "[$(date '+%H:%M:%S')] === CONTAMINATION ESTIMATION DONE ==="


#FilterMutectCalls — финальная фильтрация

echo "[$(date '+%H:%M:%S')] === FILTER MUTECT CALLS START ==="

FILT_VCF="${SOM_DIR}/somatic_filtered.vcf.gz"

if [[ ! -f "${FILT_VCF}" ]]; then
    echo "  FilterMutectCalls..."
    ${GATK} --java-options "-Xmx16g -Djava.io.tmpdir=/tmp" \
        FilterMutectCalls \
        -R "${REF}" \
        -V "${RAW_VCF}" \
        --contamination-table "${CONTAM}" \
        --tumor-segmentation "${SEGMENTS}" \
        --ob-priors "${OB_MODEL}" \
        -O "${FILT_VCF}" \
        2>> "${LOG_DIR}/filter_mutect.log"
    echo "  Done: ${FILT_VCF}"
else
    echo "  [SKIP] Already exists: ${FILT_VCF}"
fi

# Оставляем только варианты с PASS
PASS_VCF="${SOM_DIR}/somatic_pass.vcf.gz"
if [[ ! -f "${PASS_VCF}" ]]; then
    echo "  SelectVariants PASS only..."
    ${GATK} SelectVariants \
        -R "${REF}" \
        -V "${FILT_VCF}" \
        --exclude-filtered \
        -O "${PASS_VCF}" \
        2>> "${LOG_DIR}/select_pass.log"
fi

echo "[$(date '+%H:%M:%S')] === FILTER MUTECT CALLS DONE ==="

# Статистика финального VCF
echo ""
echo "  Статистика PASS вариантов:"
singularity exec ${GATK_SIF} gatk CountVariants -V "${PASS_VCF}" 2>/dev/null

echo ""
echo "========================================="
echo "  ФИНАЛЬНЫЕ РЕЗУЛЬТАТЫ (соматические):"
echo "  Все варианты: ${FILT_VCF}"
echo "  Только PASS:  ${PASS_VCF}"
echo "========================================="
echo "[$(date '+%H:%M:%S')] === ALL DONE ==="


