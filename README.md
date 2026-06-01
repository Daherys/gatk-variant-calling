# GATK Variant Calling Pipeline

Биоинформатический пайплайн для поиска зародышевых и соматических мутаций и анализа CNV на основе GATK Best Practices.

## Образцы
- ERR034529 (NA12892) — WES, здоровый индивид, 1000 Genomes Project
- HCC1143_tumor — WES, рак молочной железы (TNBC), Broad Institute
- HCC1143_normal — WES, парная норма, Broad Institute

## Скрипты
- `01_germline_calling.sh` — зародышевые мутации (HaplotypeCaller, Twist Bioscience BED)
- `02_germline_calling_nobed.sh` — зародышевые мутации без BED
- `03_somatic_calling.sh` — соматические мутации (Mutect2)
- `04_annotation.sh` — аннотация VEP + ClinVar
- `05_annotation_nobed.sh` — аннотация зародышевых без BED
- `06_cnv_analysis.sh` — CNV анализ (GATK somatic CNV pipeline)
- `07_prep_annotation.sh` — подготовка данных к аннотации

## Инструменты
- GATK 4.6.2.0 (Singularity)
- BWA-MEM2 2.2.1
- fastp 0.23.2
- VEP 115 (Singularity)
- AnnotSV 3.4

## Запуск
```bash
sbatch scripts/01_germline_calling.sh
sbatch scripts/03_somatic_calling.sh
sbatch scripts/06_cnv_analysis.sh
sbatch scripts/04_annotation.sh
```
