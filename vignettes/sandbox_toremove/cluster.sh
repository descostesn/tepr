#!/bin/sh
#SBATCH --array=1
#SBATCH --nodes=1
#SBATCH --mem=60gb
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --job-name=bgmean5cpu
#SBATCH --output=slurm_%x_%A_%a.out

module purge
eval "$(conda shell.bash hook)"
conda activate R-4.4.1

srun Rscript /g/romebioinfo/Projects/tepr/vignettes/sandbox_toremove/testpreprocessing.R

echo "Done"
