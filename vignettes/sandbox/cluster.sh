#!/bin/sh
#SBATCH --array=1
#SBATCH --nodes=1
#SBATCH --mem=300gb
#SBATCH --time=20:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --job-name=pre15cpu
#SBATCH --output=slurm_%x_%A_%a.out

module purge
eval "$(conda shell.bash hook)"
conda activate R-4.4.1

srun Rscript /g/romebioinfo/Projects/tepr/vignettes/sandbox/testpreprocessing.R

echo "Done"
