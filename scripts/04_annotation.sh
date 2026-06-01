#!/bin/bash
#SBATCH --job-name=vep_annotation
#SBATCH --output=/home/data/results/Cognito/kuznetsovads/test/logs/vep_%j.log
#SBATCH --error=/home/data/results/Cognito/kuznetsovads/test/logs/vep_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8


set -eo pipefail

source /home/kuznetsovads/miniconda3/etc/profile.d/conda.sh
conda activate bcftools


VEP_SIF="$HOME/containers/vep_115.2.sif"
VEP_DATA="$HOME/vep_data"
PLUGIN_DIR="/home/kuznetsovads/miniconda3/pkgs/ensembl-vep-115.2-pl5321h2a3209d_1/share/ensembl-vep-115.2-1"
REVEL="/home/kuznetsovads/databases/revel/revel_grch38.tsv.gz"
CLINVAR="/home/kuznetsovads/databases/clinvar/clinvar_20260208.vcf.gz"

GERM_DIR="/home/kuznetsovads/test_task/germline"
SOM_DIR="/home/kuznetsovads/test_task/somatic"
ANNO_DIR="/home/kuznetsovads/test_task/annotation"
LOG_DIR="/home/data/results/Cognito/kuznetsovads/test/logs"

mkdir -p "${ANNO_DIR}" "${LOG_DIR}"

# для конвертации хромосом
CHR_MAP="/tmp/chr_rename.txt"
CHR_MAP_BACK="/tmp/chr_rename_back.txt"

for i in {1..22} X Y M; do
    echo "chr${i} ${i}"
done > "${CHR_MAP}"

for i in {1..22} X Y M; do
    echo "${i} chr${i}"
done > "${CHR_MAP_BACK}"

#функция для аннотации

annotate_vcf() {
    local SAMPLE=$1
    local IN_VCF=$2
    local TYPE=$3

    echo "[$(date '+%H:%M:%S')] === Аннотация: ${SAMPLE} (${TYPE}) ==="

    local NOCHR_VCF="${ANNO_DIR}/${SAMPLE}_nochr.vcf.gz"
    local VEP_VCF="${ANNO_DIR}/${SAMPLE}_nochr.vep.vcf.gz"
    local CLINVAR_VCF="${ANNO_DIR}/${SAMPLE}_nochr.vep.clinvar.vcf.gz"
    local FINAL_VCF="${ANNO_DIR}/${SAMPLE}_annotated.vcf.gz"

    # chr1 -> 1 
    echo "  [1/4] Конвертация chr -> без chr..."
    bcftools annotate \
        --rename-chrs "${CHR_MAP}" \
        "${IN_VCF}" \
        -Oz -o "${NOCHR_VCF}"
    bcftools index --threads 4 "${NOCHR_VCF}"

    #VEP-
    echo "  [2/4] VEP аннотация (с REVEL, NMD, gnomAD, SIFT, PolyPhen)..."
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
        --stats_file "${ANNO_DIR}/${SAMPLE}_vep_stats.html" \
        2>> "${LOG_DIR}/${SAMPLE}_vep.log"

    bcftools index --threads 4 "${VEP_VCF}"
    echo "  VEP done"

    #ClinVar 
    echo " ClinVar аннотация..."
    bcftools annotate \
        -a "${CLINVAR}" \
        -c INFO/CLNSIG,INFO/CLNREVSTAT,INFO/CLNDN,INFO/CLNVC \
        "${VEP_VCF}" \
        -Oz -o "${CLINVAR_VCF}"
    bcftools index --threads 4 "${CLINVAR_VCF}"

    #1 -> chr1 ---
    echo "  [4/4] Конвертация обратно -> chr..."
    bcftools annotate \
        --rename-chrs "${CHR_MAP_BACK}" \
        "${CLINVAR_VCF}" \
        -Oz -o "${FINAL_VCF}"
    bcftools index --threads 4 "${FINAL_VCF}"

    # delete промежуточные файлы
    rm -f "${NOCHR_VCF}" "${NOCHR_VCF}.tbi" \
          "${VEP_VCF}" "${VEP_VCF}.tbi" \
          "${CLINVAR_VCF}" "${CLINVAR_VCF}.tbi"

    #cтатистика 
    echo ""
    echo "  === Статистика: ${SAMPLE} ==="
    bcftools stats "${FINAL_VCF}" | grep "^SN" | grep -v "no-ALTs"

    echo ""
    echo "  Топ ClinVar значимостей:"
    bcftools view "${FINAL_VCF}" | grep -v "^#" | \
        grep -oP 'CLNSIG=[^;]+' | sort | uniq -c | sort -rn | head -5

    echo "  Финальный файл: ${FINAL_VCF}"
    echo "[$(date '+%H:%M:%S')] === ${SAMPLE} DONE ==="
    echo ""
}


# ЗАРОДЫШЕВЫЕ:

echo "[$(date '+%H:%M:%S')] ===== GERMLINE ANNOTATION START ====="

annotate_vcf "ERR034529_germline" \
    "${GERM_DIR}/ERR034529_germline_final.vcf.gz" "germline"

annotate_vcf "HCC1143_normal_germline" \
    "${GERM_DIR}/HCC1143_normal_germline_final.vcf.gz" "germline"

echo "[$(date '+%H:%M:%S')] ===== GERMLINE ANNOTATION DONE ====="


# СОМАТИЧЕСКИЕ:

echo "[$(date '+%H:%M:%S')] ===== SOMATIC ANNOTATION START ====="

annotate_vcf "HCC1143_somatic" \
    "${SOM_DIR}/somatic_pass.vcf.gz" "somatic"

echo "[$(date '+%H:%M:%S')] ===== SOMATIC ANNOTATION DONE ====="

echo ""
echo "========================================="
echo "  ФИНАЛЬНЫЕ АННОТИРОВАННЫЕ ФАЙЛЫ:"
echo "  ${ANNO_DIR}/ERR034529_germline_annotated.vcf.gz"
echo "  ${ANNO_DIR}/HCC1143_normal_germline_annotated.vcf.gz"
echo "  ${ANNO_DIR}/HCC1143_somatic_annotated.vcf.gz"
echo ""
echo "  Базы данных:"
echo "  - VEP 115 cache (gnomAD, SIFT, PolyPhen, MANE, canonical)"
echo "  - REVEL v1.3 (патогенность missense вариантов)"
echo "  - ClinVar $(basename ${CLINVAR} .vcf.gz) (клиническая значимость)"
echo "  - NMD плагин (nonsense-mediated decay)"
echo "========================================="
echo "[$(date '+%H:%M:%S')] === ALL DONE ==="


