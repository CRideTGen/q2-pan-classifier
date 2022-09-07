#!/bin/bash
#SBATCH --job-name=rescript_get
#SBATCH --output=slurm_out/slurm_%x_%A_%a.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5000
#SBATCH --time=00:10:00
#SBATCH --array=0-9

module load anaconda3
conda activate qiime2-2022.2
output_dir="taxonomy_split/tax_out_0${SLURM_ARRAY_TASK_ID}"

srun qiime rescript get-ncbi-data --m-accession-ids-file accession_split/entero_ids0${SLURM_ARRAY_TASK_ID} --output-dir ${output_dir};


