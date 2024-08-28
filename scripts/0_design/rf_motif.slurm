#!/bin/bash
#SBATCH --job-name=new_rf_motif
#SBATCH --output=out.out
#SBATCH --error=err.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=g.jokent@gmail.com
#SBATCH --time=100:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=8000
#SBATCH --qos=alberto.perezant
#SBATCH --partition=gpu
#SBATCH --gpus=a100:1

##############################################

module purge
module load cuda/11.1.0
export PATH="/blue/alberto.perezant/liweichang/dev/miniconda3/envs/RFdiffusion/bin/:$PATH"
export LD_LIBRARY_PATH="/blue/alberto.perezant/liweichang/dev/miniconda3/envs/RFdiffusion/bin/:$LD_LIBRARY_PATH"

rf_diff="/orange/alberto.perezant/jokent.gaza/rfdiffusion/RFdiffusion"

python $rf_diff/scripts/run_inference.py \
    inference.output_prefix=motifs/motif \
    inference.input_pdb=../input_pdbs/NEW_motif.pdb \
    'contigmap.contigs=[A1-68/0 0-50/B69-77/0-50 ]' \
    contigmap.length=80-80 inference.num_designs=3000 \
    inference.ckpt_override_path=$rf_diff/models/Complex_beta_ckpt.pt
