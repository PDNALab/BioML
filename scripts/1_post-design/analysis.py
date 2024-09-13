import mdtraj as md
import os
import numpy as np
import glob
import re

pep1_rmsd = []

monomer_motifs = [dir for dir in os.listdir("monomer") if '-motif' in dir]
monomer_motifs.sort(key=lambda dir: int(dir.split('-')[0]))

with open('rmsd.dat', 'w') as f:
    f.write("RMSD of monomer vs in-complex\n")

# Iterate over each motif
for motif in monomer_motifs:
    complex_dir = [dir for dir in os.listdir("complex") if dir == motif][0]

    complex_pep=glob.glob(f"complex/{complex_dir}/output_nomsa/seq_custom_unrelaxed_rank_001_*.pdb")
    monomer_pep=glob.glob(f"monomer/{motif}/output_nomsa/seq_custom_unrelaxed_rank_001_*.pdb")

    traj1 = md.load(complex_pep[0])
    chainC = traj1.topology.select('chainid 2')

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


pep1_iptm = []

with open('ipTM.dat', 'w') as f:
    f.write("ipTM values of complex\n")

# Iterate over each motif
for i in range(1, len(monomer_motifs)+1):
    filepath = f"complex/{i}-motif_*/output_nomsa/log.txt"
    
    for filename in glob.glob(filepath):
        with open(filename, 'r') as file:
            for line in file:
                if "rank_001_" in line:
                    fields = line.split()
                    if len(fields) >= 6:
                        match = re.search(r'ipTM=([^ ]+)', fields[5])
                        if match:
                            pep1_iptm.append(float(match.group(1)))
                            with open('ipTM.dat', 'a') as f:
                                f.write("{0}\n".format(match.group(1)))

with open('passed.dat', 'w') as f:
    f.write('List of successful proteins\n')

for i, (rog, rmsd, iptm) in enumerate(zip(pep1_rog, pep1_rmsd, pep1_iptm)):
    if rog < 1.4 and rmsd < 0.2 and iptm > 7.0:
        with open('passed.dat', 'a') as f:
            f.write(monomer_motifs[i] + '\n')
