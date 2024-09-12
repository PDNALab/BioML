#!/bin/bash

motif=$1
scripts_dir=$(pwd)

# Prepare Files
cd $scripts_dir\/../../submissions/$motif
seq_dir=$(pwd)
echo "Working directory: $seq_dir"
mkdir sanity; mkdir sanity/monomer; mkdir sanity/complex;

# Get list of successful designs (pAE < 10)
awk '
$5 > 85 {
    match($9, /motif_[0-9]+_/)
    if (RSTART > 0) {
        print $5 "\t" substr($9, RSTART, RLENGTH)
    }
}' af2/designed_score.sc | \
sort -n | \
awk '{print $2}' > sanity/motifs_label.dat

while IFS= read -r line; do
    grep "$line" sequence_designs.silent | grep "ANNOTATED" | \
    awk '
    {
        # Remove everything inside brackets, including the brackets themselves
        gsub(/\[[^\]]*\]/, "")
        # Get binder sequence
        print substr($2, 1, 80)
    }'
done < sanity/motifs_label.dat > sanity/seqs.dat

# Monomer

cd sanity/monomer
cp $scripts_dir/nomsa.slurm .; cp $scripts_dir/write_monomer.py .
python write_monomer.py
sbatch nomsa.slurm
rm af2.err af2.out nomsa.slurm write_monomer.py
cd ../../

# Complex
cd sanity/complex               
cp $scripts_dir/nomsa.slurm .; cp $scripts_dir/write_complex.py .; cp $scripts_dir/CD20.msa .
python write_complex.py
sbatch nomsa.slurm
rm af2.err af2.out nomsa.slurm write_complex.py CD20.msa
cd ../../

cd $scripts_dir
