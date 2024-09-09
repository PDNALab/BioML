#!/bin/bash
#SBATCH --job-name=af2
#SBATCH --output=af2.out
#SBATCH --error=af2.err
##SBATCH --mail-type=ALL
##SBATCH --mail-user=g.jokent@gmail.com
#SBATCH --time=100:00:00
#SBATCH --ntasks=1
#SBATCH --distribution=cyclic:cyclic
#SBATCH --mem-per-cpu=40000
#SBATCH --qos=alberto.perezant
#SBATCH --cpus-per-gpu=1
#SBATCH --gpus-per-task=1
#SBATCH --partition=gpu
#SBATCH --constraint=a100

##########################################################

export PATH="/blue/alberto.perezant/liweichang/dev/miniconda3/envs/colabfold_casp16/bin/:$PATH"
export LD_LIBRARY_PATH="/blue/alberto.perezant/liweichang/dev/miniconda3/envs/colabfold_casp16/bin/:$LD_LIBRARY_PATH"
export PYTHONPATH=/blue/alberto.perezant/liweichang/dev/miniconda3/envs/colabfold_casp16/lib/python3.11/site-packages/:$PYTHONPATH

module load ufrc cuda/12 gcc

directories=($(ls -vd *motif_*))

# Loop through the directories and change to each directory
for dir in "${directories[@]}"; do
  dir="${dir%/}"
  cd "$dir"
  mkdir output_nomsa
  colabfold_batch seq_custom.a3m output_nomsa  --model-type 'alphafold2_ptm'
  cd .. 
done

