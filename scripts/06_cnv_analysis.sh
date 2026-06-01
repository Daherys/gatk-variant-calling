#!/bin/bash
#SBATCH --job-name=gatk_cnv
#SBATCH --output=/home/data/results/Cognito/kuznetsovads/test/logs/cnv_%j.log
#SBATCH --error=/home/data/results/Cognito/kuznetsovads/test/logs/cnv_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8


# Соматический CNV анализ HCC1143 tumor vs normal
# Пайплайн: CollectReadCounts -> PoN -> DenoiseReadCounts ->
#           CollectAllelicCounts -> ModelSegments ->
#           CallCopyRatioSegments -> PlotModeledSegments

set -eo pipefail

source /home/kuznetsovads/miniconda3/etc/profile.d/conda.sh


GATK_SIF="/home/kuznetsovads/containers/gatk_4.6.2.0.sif"
GATK="singularity exec ${GATK_SIF} gatk"

REF="/home/kuznetsovads/d1/ref_bwa_mem2/hg38.analysisSet.fa"
DICT="/home/kuznetsovads/d1/ref_bwa_mem2/hg38.analysisSet.dict"

T_BAM="/home/kuznetsovads/test_task/align_install/hcc1143_T_clean.bam"
N_BAM="/home/kuznetsovads/test_task/align_install/hcc1143_N_clean.bam"

BED="/home/kuznetsovads/test_task/exome_bed/exome_calling_regions.v1.bed"
INTERVALS="/home/kuznetsovads/test_task/exome_bed/exome_calling_regions.interval_list"

GNOMAD="/home/kuznetsovads/databases/gatk_bundle/small_exac_common_3.hg38.vcf.gz"

CNV_DIR="/home/kuznetsovads/test_task/cnv"
LOG_DIR="/home/data/results/Cognito/kuznetsovads/test/logs"

mkdir -p "${CNV_DIR}" "${LOG_DIR}"

#промежуточные файлы
T_COUNTS="${CNV_DIR}/tumor.counts.hdf5"
N_COUNTS="${CNV_DIR}/normal.counts.hdf5"
PON="${CNV_DIR}/pon.hdf5"
DENOISED_CR="${CNV_DIR}/tumor.denoisedCR.tsv"
STANDARDIZED_CR="${CNV_DIR}/tumor.standardizedCR.tsv"
T_ALLELIC="${CNV_DIR}/tumor.allelicCounts.tsv"
N_ALLELIC="${CNV_DIR}/normal.allelicCounts.tsv"
SEGMENTS_PREFIX="${CNV_DIR}/tumor"
CALLED_SEGS="${CNV_DIR}/tumor.called.seg"
PLOT_DIR="${CNV_DIR}/plots"

mkdir -p "${PLOT_DIR}"

#BED -> interval_list

echo "[$(date '+%H:%M:%S')] === ШАГ 0: BedToIntervalList ==="

if [[ ! -f "${INTERVALS}" ]]; then
    ${GATK} BedToIntervalList \
        -I "${BED}" \
        -O "${INTERVALS}" \
        -SD "${DICT}" \
        2>> "${LOG_DIR}/cnv_step0.log"
    echo "  Done: ${INTERVALS}"
else
    echo "  [SKIP] Already exists: ${INTERVALS}"
fi


#CollectReadCounts

echo "[$(date '+%H:%M:%S')] === ШАГ 1: CollectReadCounts ==="

if [[ ! -f "${T_COUNTS}" ]]; then
    echo "  tumor..."
    ${GATK} --java-options "-Xmx40g -Djava.io.tmpdir=/tmp" \
        CollectReadCounts \
        -I "${T_BAM}" \
        -L "${INTERVALS}" \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O "${T_COUNTS}" \
        2>> "${LOG_DIR}/cnv_step1_tumor.log"
    echo "  Done: ${T_COUNTS}"
else
    echo "  [SKIP] Already exists: ${T_COUNTS}"
fi

if [[ ! -f "${N_COUNTS}" ]]; then
    echo "  normal..."
    ${GATK} --java-options "-Xmx40g -Djava.io.tmpdir=/tmp" \
        CollectReadCounts \
        -I "${N_BAM}" \
        -L "${INTERVALS}" \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O "${N_COUNTS}" \
        2>> "${LOG_DIR}/cnv_step1_normal.log"
    echo "  Done: ${N_COUNTS}"
else
    echo "  [SKIP] Already exists: ${N_COUNTS}"
fi

echo "[$(date '+%H:%M:%S')] === ШАГ 1 DONE ==="


#CreateReadCountPanelOfNormals

echo "[$(date '+%H:%M:%S')] === ШАГ 2: CreateReadCountPanelOfNormals ==="

if [[ ! -f "${PON}" ]]; then
    ${GATK} --java-options "-Xmx40g -Djava.io.tmpdir=/tmp" \
        CreateReadCountPanelOfNormals \
        -I "${N_COUNTS}" \
        --minimum-interval-median-percentile 5.0 \
        -O "${PON}" \
        2>> "${LOG_DIR}/cnv_step2.log"
    echo "  Done: ${PON}"
else
    echo "  [SKIP] Already exists: ${PON}"
fi

echo "[$(date '+%H:%M:%S')] === ШАГ 2 DONE ==="

#DenoiseReadCounts
echo "[$(date '+%H:%M:%S')] === ШАГ 3: DenoiseReadCounts ==="

if [[ ! -f "${DENOISED_CR}" ]]; then
    ${GATK} --java-options "-Xmx40g -Djava.io.tmpdir=/tmp" \
        DenoiseReadCounts \
        -I "${T_COUNTS}" \
        --count-panel-of-normals "${PON}" \
        --standardized-copy-ratios "${STANDARDIZED_CR}" \
        --denoised-copy-ratios "${DENOISED_CR}" \
        2>> "${LOG_DIR}/cnv_step3.log"
    echo "  Done: ${DENOISED_CR}"
else
    echo "  [SKIP] Already exists: ${DENOISED_CR}"
fi

echo "[$(date '+%H:%M:%S')] === ШАГ 3 DONE ==="

# CollectAllelicCounts

echo "[$(date '+%H:%M:%S')] === ШАГ 4: CollectAllelicCounts ==="

if [[ ! -f "${T_ALLELIC}" ]]; then
    echo "  tumor..."
    ${GATK} --java-options "-Xmx28g -Djava.io.tmpdir=/tmp" \
        CollectAllelicCounts \
        -I "${T_BAM}" \
        -L "${GNOMAD}" \
        -R "${REF}" \
        -O "${T_ALLELIC}" \
        2>> "${LOG_DIR}/cnv_step4_tumor.log"
    echo "  Done: ${T_ALLELIC}"
else
    echo "  [SKIP] Already exists: ${T_ALLELIC}"
fi

if [[ ! -f "${N_ALLELIC}" ]]; then
    echo "  normal..."
    ${GATK} --java-options "-Xmx28g -Djava.io.tmpdir=/tmp" \
        CollectAllelicCounts \
        -I "${N_BAM}" \
        -L "${GNOMAD}" \
        -R "${REF}" \
        -O "${N_ALLELIC}" \
        2>> "${LOG_DIR}/cnv_step4_normal.log"
    echo "  Done: ${N_ALLELIC}"
else
    echo "  [SKIP] Already exists: ${N_ALLELIC}"
fi

echo "[$(date '+%H:%M:%S')] === ШАГ 4 DONE ==="

ModelSegments

echo "[$(date '+%H:%M:%S')] === ШАГ 5: ModelSegments ==="

if [[ ! -f "${SEGMENTS_PREFIX}.modelFinal.seg" ]]; then
    ${GATK} --java-options "-Xmx28g -Djava.io.tmpdir=/tmp" \
        ModelSegments \
        --denoised-copy-ratios "${DENOISED_CR}" \
        --allelic-counts "${T_ALLELIC}" \
        --normal-allelic-counts "${N_ALLELIC}" \
        --output "${CNV_DIR}" \
        --output-prefix tumor \
        2>> "${LOG_DIR}/cnv_step5.log"
    echo "  Done: ${SEGMENTS_PREFIX}.modelFinal.seg"
else
    echo "  [SKIP] Already exists"
fi

echo "[$(date '+%H:%M:%S')] === ШАГ 5 DONE ==="


#CallCopyRatioSegments

echo "[$(date '+%H:%M:%S')] === ШАГ 6: CallCopyRatioSegments ==="

if [[ ! -f "${CALLED_SEGS}" ]]; then
    ${GATK} --java-options "-Xmx16g -Djava.io.tmpdir=/tmp" \
        CallCopyRatioSegments \
        -I "${SEGMENTS_PREFIX}.cr.seg" \
        -O "${CALLED_SEGS}" \
        2>> "${LOG_DIR}/cnv_step6.log"
    echo "  Done: ${CALLED_SEGS}"
else
    echo "  [SKIP] Already exists: ${CALLED_SEGS}"
fi

echo "[$(date '+%H:%M:%S')] === ШАГ 6 DONE ==="


# PlotModeledSegments

echo "[$(date '+%H:%M:%S')] === ШАГ 7: PlotModeledSegments ==="

if [[ ! -f "${PLOT_DIR}/tumor.modeled.png" ]]; then
    ${GATK} --java-options "-Xmx16g -Djava.io.tmpdir=/tmp" \
        PlotModeledSegments \
        --denoised-copy-ratios "${DENOISED_CR}" \
        --allelic-counts "${SEGMENTS_PREFIX}.hets.tsv" \
        --segments "${SEGMENTS_PREFIX}.modelFinal.seg" \
        --sequence-dictionary "${DICT}" \
        --minimum-contig-length 46709983 \
        --output "${PLOT_DIR}" \
        --output-prefix tumor \
        2>> "${LOG_DIR}/cnv_step7.log"
    echo "  Done: plots в ${PLOT_DIR}"
else
    echo "  [SKIP] Already exists"
fi

echo "[$(date '+%H:%M:%S')] === ШАГ 7 DONE ==="

#статистика

echo ""
echo "=== Статистика CNV ==="
TOTAL=$(grep -v "^@\|^CONTIG" "${CALLED_SEGS}" | wc -l)
AMP=$(grep -v "^@\|^CONTIG" "${CALLED_SEGS}" | awk '$6=="+"' | wc -l)
DEL=$(grep -v "^@\|^CONTIG" "${CALLED_SEGS}" | awk '$6=="-"' | wc -l)
NEU=$(grep -v "^@\|^CONTIG" "${CALLED_SEGS}" | awk '$6=="0"' | wc -l)

echo "Всего сегментов: ${TOTAL}"
echo "Амплификации (+): ${AMP}"
echo "Делеции (-): ${DEL}"
echo "Нейтральные (0): ${NEU}"

echo ""
echo "Топ-10 амплификаций (по log2CR, >10 точек):"
grep -v "^@\|^CONTIG" "${CALLED_SEGS}" | awk '$6=="+" && $4>10' | \
    sort -k5 -rn | head -10 | \
    awk '{printf "%s:%s-%s\tlog2CR=%.3f\tpoints=%s\n", $1,$2,$3,$5,$4}'

echo ""
echo "Топ-10 делеций (по log2CR, >10 точек):"
grep -v "^@\|^CONTIG" "${CALLED_SEGS}" | awk '$6=="-" && $4>10' | \
    sort -k5 -n | head -10 | \
    awk '{printf "%s:%s-%s\tlog2CR=%.3f\tpoints=%s\n", $1,$2,$3,$5,$4}'

echo ""
echo "========================================="
echo "  РЕЗУЛЬТАТЫ CNV:"
echo "  Сегменты:  ${CALLED_SEGS}"
echo "  Графики:   ${PLOT_DIR}/"
echo "========================================="
echo "[$(date '+%H:%M:%S')] === ALL DONE ==="


