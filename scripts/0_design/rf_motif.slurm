#!/bin/bash
#SBATCH --job-name=rf_motif
#SBATCH --output=out.out
#SBATCH --error=err.err
##SBATCH --mail-type=ALL
##SBATCH --mail-user=g.jokent@gmail.com
#SBATCH --time=300:00:00
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
    inference.input_pdb=NEW_motif.pdb \
    'contigmap.contigs=[A1-165/0 B166-330/0 0-50/C386-391/0-50 ]' \
    'ppi.hotspot_res=[A122,A123,A124,A125,A126,A127,A128,A129,A130,A131,A132,A133,A134,A135,A136,A137]' \
    contigmap.length=80-80 inference.num_designs=3000 \
    inference.ckpt_override_path=$rf_diff/models/Complex_beta_ckpt.pt
