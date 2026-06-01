#!/bin/bash
#SBATCH --job-name=revel_prep
#SBATCH --output=/home/data/results/Cognito/kuznetsovads/test/logs/revel_prep_%j.log
#SBATCH --error=/home/data/results/Cognito/kuznetsovads/test/logs/revel_prep_%j.err
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

set -eo pipefail

source /home/kuznetsovads/miniconda3/etc/profile.d/conda.sh
conda activate bcftools

cd /home/kuznetsovads/databases/revel

echo "[$(date '+%H:%M:%S')] Конвертация CSV -> TSV..."
cat revel_with_transcript_ids \
    | tr "," "\t" \
    | awk 'NR==1{print "#chr\thg19_pos\tgrch38_pos\tref\talt\taaref\taaalt\tREVEL\tEnsembl_transcriptid"; next} {print}' \
    > revel_grch38.tsv

echo "[$(date '+%H:%M:%S')] Сжатие bgzip..."
bgzip revel_grch38.tsv

echo "[$(date '+%H:%M:%S')] Индексирование tabix..."
tabix -s 1 -b 3 -e 3 revel_grch38.tsv.gz

echo "[$(date '+%H:%M:%S')] Done!"
ls -lh revel_grch38.tsv.gz


