import mdtraj as md
import os
import numpy as np
import glob
import re

pep1_rmsd = []

monomer_motifs = [dir for dir in os.listdir("monomer") if '-motif' in dir]
monomer_motifs.sort(key=lambda dir: int(dir.split('-')[0]))

with open('rmsd.dat', 'w') as f:
    f.write("RMSD of monomer vs full CAR\n")

# Iterate over each motif
for motif in monomer_motifs:
    complex_dir = [dir for dir in os.listdir("full") if dir == motif][0]

    complex_pep=glob.glob(f"full/{complex_dir}/output_nomsa/seq_custom_unrelaxed_rank_001_*.pdb")
    monomer_pep=glob.glob(f"monomer/{motif}/output_nomsa/seq_custom_unrelaxed_rank_001_*.pdb")

    traj1 = md.load(complex_pep[0])
    chainC = traj1.topology.select('resid 21 to 100')

    traj2 = md.load(monomer_pep[0])
    chainA = traj2.topology.select('chainid 0')

    rmsd = md.rmsd(traj1, traj2, 0, atom_indices=chainC, ref_atom_indices=chainA)
    pep1_rmsd.append(rmsd[0])
    
    with open('rmsd.dat', 'a') as f:
        f.write("{0}\n".format(rmsd[0]))


pep1_rog = []

monomer_motifs = [dir for dir in os.listdir("monomer") if '-motif' in dir]
monomer_motifs.sort(key=lambda dir: int(dir.split('-')[0]))

with open('rog.dat', 'w') as f:
    f.write("Radius of Gyration of monomer\n")

# Iterate over each motif
for motif in monomer_motifs:
    complex_dir = [dir for dir in os.listdir("complex") if dir == motif][0]

    monomer_pep=glob.glob(f"monomer/{motif}/output_nomsa/seq_custom_unrelaxed_rank_001_*.pdb")

    traj2 = md.load(monomer_pep[0])
    chainA = traj2.topology.select('chainid 0')

    rmsd = md.compute_rg(traj2)
    pep1_rog.append(rmsd[0])
    
    with open('rog.dat', 'a') as f:
        f.write("{0}\n".format(rmsd[0]))

with open('passed.dat', 'w') as f:
    f.write('List of successful proteins\n')

for i, (rog, rmsd) in enumerate(zip(pep1_rog, pep1_rmsd)):
    if rog < 1.4 and rmsd < 0.2:
        with open('passed.dat', 'a') as f:
            f.write(monomer_motifs[i] + '\n')
