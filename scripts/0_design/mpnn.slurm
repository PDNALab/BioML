#!/bin/bash
#SBATCH --job-name=MPNN
#SBATCH --output=MPNN_s.out
#SBATCH --error=MPNN_s.err
##SBATCH --mail-type=ALL
##SBATCH --mail-user=g.jokent@gmail.com
#SBATCH --time=200:00:00
#SBATCH --ntasks=1
#SBATCH --distribution=cyclic:cyclic
#SBATCH --mem-per-cpu=50000
#SBATCH --qos=alberto.perezant
#SBATCH --cpus-per-gpu=1
#SBATCH --gpus-per-task=1
#SBATCH --partition=gpu
#SBATCH --constraint=a100

##########################################################

module purge

module load cuda/11.1.0
export PATH="/orange/alberto.perezant/jokent.gaza/dl_binder_design/include/envs/dl_binder_design/bin/:$PATH"
export LD_LIBRARY_PATH="/orange/alberto.perezant/jokent.gaza/dl_binder_design/include/envs/dl_binder_design/bin/:$LD_LIBRARY_PATH"

dl_binder="/orange/alberto.perezant/jokent.gaza/dl_binder_design"

# retain sequence of motif
python $dl_binder/helper_scripts/addFIXEDlabels.py                      \
    --pdbdir ./motifs --trbdir ./motifs --verbose

# convert to silent files
silentfrompdbs ./motifs/*.pdb > ./motifs/backbone_designs.silent

# proteinmpnn
python $dl_binder/mpnn_fr/dl_interface_design.py                        \
    -silent ./motifs/backbone_designs.silent                            \
    -outsilent sequence_designs.silent -debug

# run af2 and compare
python $dl_binder/af2_initial_guess/predict.py                          \
    -silent sequence_designs.silent                                     \
    -outpdbdir ./af2                                                    \
    -outsilent ./af2/designed_af.silent                                 \
    -scorefilename ./af2/designed_score.sc

