import os
import re
import itertools
from itertools import combinations

def list_all_pairs(input_list):
    all_pairs = list(combinations(input_list, 2))
    return all_pairs

seq_dict = {}
name_all = []
label_all = []

CD20='MRESKTLGAVQIMNGLFHIALGGLLMIPAGIYAPICVTVWYPLWGGIMYIISGSLLAATEKNSRKCLVKGKMIMNSLSLFAAISGMILSIMDILNIKISHFLKMESLNFIRAHTPYINIYNCEPANPSEKNSPSTQYCYSIQSLFLGILSVMLIFAFFQELVIAG'

with open('../seqs.dat', 'r') as f:
  for line in f:
    names = line.strip("\n")
    name_all.append(names)

with open('../motifs_label.dat', 'r') as f:
  for line in f:
    parts = line.split('_')
    names = '_'.join(parts[0:2]).strip("\n")
    label_all.append(names)

for i, (j, k) in enumerate(zip(name_all, label_all), 1):
  af2_file = open("seq_custom.a3m", "w")
  full_seq="{0}".format(j)
  af2_file.write("#80\t1\n")
  af2_file.write(">101\n{0}".format(full_seq))

  os.system("mkdir {0}-{1}".format(i, k))
  af2_file.close()
  os.system("mv seq_custom.a3m {0}-{1}".format(i, k))
