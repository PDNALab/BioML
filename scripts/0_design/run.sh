#!/bin/bash

motif=$1
model=$2
scripts_dir=$(pwd)

# Prepare Files
mkdir $scripts_dir\/../../submissions/$motif; cd $scripts_dir\/../../submissions/$motif
seq_dir=$(pwd)
echo "Working directory: $seq_dir"
mkdir motifs; mkdir af2;
cp $scripts_dir/*.slurm .

# Replace RFDiff parameters based on input_pdb/*motif.pdb file
input_pdbs="$scripts_dir/../input_pdbs/$motif"_"motif.pdb"
sed -i "s|NEW_motif.pdb|$input_pdbs|g" rf_motif.slurm

A_ter=$(grep 'TER' $input_pdbs | awk 'NR == 1{print $5}')
B_ter=$(grep 'TER' $input_pdbs | awk 'NR == 2{print $5}')

sed -i "s|165|$A_ter|g" rf_motif.slurm
sed -i "s|166|$((A_ter+1))|g" rf_motif.slurm
sed -i "s|330|$B_ter|g" rf_motif.slurm

C_nter=$(awk 'substr($0, 22, 1) == "C"' $input_pdbs | awk '{print $6}' | sort | uniq | head -n 1)
C_cter=$(awk 'substr($0, 22, 1) == "C"' $input_pdbs | awk '{print $6}' | sort | uniq | tail -n 1)

sed -i "s|386|$C_nter|g" rf_motif.slurm
sed -i "s|391|$C_cter|g" rf_motif.slurm

model_name="Complex_$model"_"ckpt"
sed -i "s|Complex_beta_ckpt|$model_name|g" rf_motif.slurm

if [ "$motif" == "OF" ]; then
  hotspot1=106; hotspot2=127
elif [ "$motif" == "OB" ]; then
  hotspot1=121; hotspot2=136
else
  # RT numbering scheme  
  hotspot1=122; hotspot2=137
fi

hotspot_resides=""

for ((i=hotspot1; i<=hotspot2; i++)); do
  if [ $i -eq $((hotspot2)) ]; then
    hotspot_residues+="A$i"  
  else
    hotspot_residues+="A$i,"  
  fi
done

sed -i "s|hotspot-residues|$hotspot_residues|g" rf_motif.slurm

# Run
#job_id=$(sbatch rf_motif.slurm | awk '{print $NF}')
#sbatch --dependency=afterok:$job_id mpnn.slurm

cd $scripts_dir
