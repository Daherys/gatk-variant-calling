#!/bin/bash
#SBATCH --job-name=hc_nobed_ERR034529
#SBATCH --output=/home/data/results/Cognito/kuznetsovads/test/logs/hc_nobed_%j.log
#SBATCH --error=/home/data/results/Cognito/kuznetsovads/test/logs/hc_nobed_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8


set -eo pipefail

GATK_SIF="/home/kuznetsovads/containers/gatk_4.6.2.0.sif"
REF="/home/kuznetsovads/d1/ref_bwa_mem2/hg38.analysisSet.fa"
DBSNP="/home/kuznetsovads/databases/gatk_bundle/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
IN_BAM="/home/kuznetsovads/test_task/bqsr/ERR034529.bqsr.bam"
GERM_DIR="/home/kuznetsovads/test_task/germline"
LOG_DIR="/home/data/results/Cognito/kuznetsovads/test/logs"

mkdir -p "${GERM_DIR}" "${LOG_DIR}"

#HaplotypeCaller без BED файла


GVCF="${GERM_DIR}/ERR034529_nobed.g.vcf.gz"
RAW_VCF="${GERM_DIR}/ERR034529_nobed_raw.vcf.gz"
FINAL_VCF="${GERM_DIR}/ERR034529_nobed_final.vcf.gz"

echo "[$(date '+%H:%M:%S')] === HAPLOTYPECALLER (no BED) START ==="

if [[ ! -f "${GVCF}" ]]; then
    singularity exec ${GATK_SIF} gatk \
        --java-options "-Xmx40g -Djava.io.tmpdir=/tmp" \
        HaplotypeCaller \
        -R "${REF}" \
        -I "${IN_BAM}" \
        -ERC GVCF \
        --dbsnp "${DBSNP}" \
        --native-pair-hmm-threads 8 \
        -O "${GVCF}" \
        2>> "${LOG_DIR}/ERR034529_nobed_hc.log"
    echo "  Done: ${GVCF}"
else
    echo "  [SKIP] Already exists: ${GVCF}"
fi

echo "[$(date '+%H:%M:%S')] === HAPLOTYPECALLER DONE ==="

echo "[$(date '+%H:%M:%S')] === GENOTYPEGVCFS START ==="

if [[ ! -f "${RAW_VCF}" ]]; then
    singularity exec ${GATK_SIF} gatk \
        --java-options "-Xmx40g -Djava.io.tmpdir=/tmp" \
        GenotypeGVCFs \
        -R "${REF}" \
        -V "${GVCF}" \
        --dbsnp "${DBSNP}" \
        -O "${RAW_VCF}" \
        2>> "${LOG_DIR}/ERR034529_nobed_genotype.log"
    echo "  Done: ${RAW_VCF}"
else
    echo "  [SKIP] Already exists: ${RAW_VCF}"
fi

echo "[$(date '+%H:%M:%S')] === GENOTYPEGVCFS DONE ==="


#Hard Filtering
echo "[$(date '+%H:%M:%S')] === FILTERING START ==="

SNP_VCF="${GERM_DIR}/ERR034529_nobed_snp.vcf.gz"
SNP_FILT="${GERM_DIR}/ERR034529_nobed_snp.filtered.vcf.gz"
INDEL_VCF="${GERM_DIR}/ERR034529_nobed_indel.vcf.gz"
INDEL_FILT="${GERM_DIR}/ERR034529_nobed_indel.filtered.vcf.gz"

singularity exec ${GATK_SIF} gatk SelectVariants \
    -R "${REF}" -V "${RAW_VCF}" \
    --select-type-to-include SNP -O "${SNP_VCF}" \
    2>> "${LOG_DIR}/ERR034529_nobed_filter.log"

singularity exec ${GATK_SIF} gatk VariantFiltration \
    -R "${REF}" -V "${SNP_VCF}" \
    --filter-expression "QD < 2.0"             --filter-name "QD2"              \
    --filter-expression "QUAL < 30.0"           --filter-name "QUAL30"           \
    --filter-expression "FS > 60.0"             --filter-name "FS60"             \
    --filter-expression "SOR > 3.0"             --filter-name "SOR3"             \
    --filter-expression "MQ < 40.0"             --filter-name "MQ40"             \
    --filter-expression "MQRankSum < -12.5"     --filter-name "MQRankSum-12.5"   \
    --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O "${SNP_FILT}" 2>> "${LOG_DIR}/ERR034529_nobed_filter.log"


singularity exec ${GATK_SIF} gatk SelectVariants \
    -R "${REF}" -V "${RAW_VCF}" \
    --select-type-to-include INDEL \
    --select-type-to-include MIXED \
    -O "${INDEL_VCF}" 2>> "${LOG_DIR}/ERR034529_nobed_filter.log"

singularity exec ${GATK_SIF} gatk VariantFiltration \
    -R "${REF}" -V "${INDEL_VCF}" \
    --filter-expression "QD < 2.0"              --filter-name "QD2"               \
    --filter-expression "QUAL < 30.0"            --filter-name "QUAL30"            \
    --filter-expression "FS > 200.0"             --filter-name "FS200"             \
    --filter-expression "SOR > 10.0"             --filter-name "SOR10"             \
    --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
    -O "${INDEL_FILT}" 2>> "${LOG_DIR}/ERR034529_nobed_filter.log"

singularity exec ${GATK_SIF} gatk MergeVcfs \
    -I "${SNP_FILT}" \
    -I "${INDEL_FILT}" \
    -O "${FINAL_VCF}" 2>> "${LOG_DIR}/ERR034529_nobed_filter.log"
echo "[$(date '+%H:%M:%S')] === FILTERING DONE ==="

#cтатистика для сравнения с BED результатом
# ============================================================
echo ""
echo "=== СТАТИСТИКА (no BED) ==="
singularity exec ${GATK_SIF} gatk CountVariants -V "${FINAL_VCF}" 2>/dev/null

conda activate bcftools 2>/dev/null || true
bcftools stats "${FINAL_VCF}" | grep -E "^SN|^TSTV" | grep -v "^#"
echo ""
echo "========================================="
echo "  РЕЗУЛЬТАТ: ${FINAL_VCF}"
echo "  Сравните Ti/Tv с предыдущим результатом:"
echo "  С BED:    32,629 вариантов, Ti/Tv=2.57"
echo "  Без BED:  см. выше"
echo "  WGS -> Ti/Tv ~2.0-2.1, ~200-400K вариантов"
echo "  WES -> Ti/Tv ~2.8-3.0, ~80-100K вариантов"
echo "========================================="
echo "[$(date '+%H:%M:%S')] === ALL DONE ==="



