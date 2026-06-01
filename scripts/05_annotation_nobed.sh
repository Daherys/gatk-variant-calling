#!/bin/bash
#SBATCH --job-name=vep_ERR_wgs
#SBATCH --output=/home/data/results/Cognito/kuznetsovads/test/logs/vep_ERR_wgs_%j.log
#SBATCH --error=/home/data/results/Cognito/kuznetsovads/test/logs/vep_ERR_wgs_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8


set -eo pipefail

source /home/kuznetsovads/miniconda3/etc/profile.d/conda.sh
conda activate bcftools


VEP_SIF="/home/kuznetsovads/containers/vep_115.2.sif"
VEP_DATA="/home/kuznetsovads/vep_data"
PLUGIN_DIR="/home/kuznetsovads/miniconda3/pkgs/ensembl-vep-115.2-pl5321h2a3209d_1/share/ensembl-vep-115.2-1"
REVEL="/home/kuznetsovads/databases/revel/revel_grch38.tsv.gz"
CLINVAR="/home/kuznetsovads/databases/clinvar/clinvar_20260208.vcf.gz"

IN_VCF="/home/kuznetsovads/test_task/germline/ERR034529_nobed_final.vcf.gz"
ANNO_DIR="/home/kuznetsovads/test_task/annotation"
LOG_DIR="/home/data/results/Cognito/kuznetsovads/test/logs"

NOCHR_VCF="${ANNO_DIR}/ERR034529_wgs_nochr.vcf.gz"
VEP_VCF="${ANNO_DIR}/ERR034529_wgs_nochr.vep.vcf.gz"
CLINVAR_VCF="${ANNO_DIR}/ERR034529_wgs_nochr.vep.clinvar.vcf.gz"
FINAL_VCF="${ANNO_DIR}/ERR034529_wgs_annotated.vcf.gz"

mkdir -p "${ANNO_DIR}" "${LOG_DIR}"

#для конвертации хромосом
CHR_MAP="/tmp/chr_rename.txt"
CHR_MAP_BACK="/tmp/chr_rename_back.txt"

for i in {1..22} X Y M; do echo "chr${i} ${i}"; done > "${CHR_MAP}"
for i in {1..22} X Y M; do echo "${i} chr${i}"; done > "${CHR_MAP_BACK}"


echo "[$(date '+%H:%M:%S')] [1/4] Конвертация chr -> без chr..."

bcftools annotate \
    --rename-chrs "${CHR_MAP}" \
    "${IN_VCF}" \
    -Oz -o "${NOCHR_VCF}"
bcftools index --threads 4 "${NOCHR_VCF}"


VEP аннотация

echo "[$(date '+%H:%M:%S')] [2/4] VEP аннотация (324K вариантов, ~2-4 часа)..."

singularity exec \
    --bind "${VEP_DATA}:${VEP_DATA}" \
    --bind "${PLUGIN_DIR}:${PLUGIN_DIR}" \
    --bind "/home/kuznetsovads/databases:/home/kuznetsovads/databases" \
    --bind "${ANNO_DIR}:${ANNO_DIR}" \
    "${VEP_SIF}" vep \
    --dir_cache "${VEP_DATA}" \
    --dir_plugins "${PLUGIN_DIR}" \
    --fasta "${VEP_DATA}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" \
    --cache --offline \
    --species homo_sapiens --assembly GRCh38 \
    --check_ref --check_existing \
    --symbol --canonical --biotype --mane \
    --af_gnomad --max_af \
    --hgvs \
    --sift b --polyphen b --protein \
    --plugin REVEL,file="${REVEL}" \
    --plugin NMD \
    --input_file "${NOCHR_VCF}" \
    --output_file "${VEP_VCF}" \
    --vcf --compress_output bgzip --force_overwrite \
    --fork 8 \
    --stats_file "${ANNO_DIR}/ERR034529_wgs_vep_stats.html" \
    2>> "${LOG_DIR}/ERR034529_wgs_vep.log"

bcftools index --threads 4 "${VEP_VCF}"
echo "[$(date '+%H:%M:%S')] VEP done"

#ClinVar

echo "[$(date '+%H:%M:%S')] [3/4] ClinVar аннотация..."

bcftools annotate \
    -a "${CLINVAR}" \
    -c INFO/CLNSIG,INFO/CLNREVSTAT,INFO/CLNDN,INFO/CLNVC \
    "${VEP_VCF}" \
    -Oz -o "${CLINVAR_VCF}"
bcftools index --threads 4 "${CLINVAR_VCF}"

#1 -> chr1
echo "[$(date '+%H:%M:%S')] [4/4] Конвертация обратно -> chr..."

bcftools annotate \
    --rename-chrs "${CHR_MAP_BACK}" \
    "${CLINVAR_VCF}" \
    -Oz -o "${FINAL_VCF}"
bcftools index --threads 4 "${FINAL_VCF}"

# Удаляем промежуточные файлы
rm -f "${NOCHR_VCF}" "${NOCHR_VCF}.tbi" \
      "${VEP_VCF}" "${VEP_VCF}.tbi" \
      "${CLINVAR_VCF}" "${CLINVAR_VCF}.tbi"

#статистика

echo ""
echo "=== Статистика ERR034529 WGS ==="
bcftools stats "${FINAL_VCF}" | grep "^SN" | grep -v "no-ALTs"

echo ""
echo "=== ClinVar значимости ==="
bcftools view "${FINAL_VCF}" | grep -v "^#" | \
    grep -oP 'CLNSIG=[^;]+' | sort | uniq -c | sort -rn | head -10

echo ""
echo "========================================="
echo "  Финальный файл: ${FINAL_VCF}"
echo "  Сравнение с BED результатом:"
echo "  С BED:    32,629 вариантов, Ti/Tv=2.57"
echo "  Без BED: 324,148 вариантов, Ti/Tv=2.09 -> WGS подтверждён"
echo "========================================="
echo "[$(date '+%H:%M:%S')] === ALL DONE ==="


