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
  os.system("mkdir {0}-{1}".format(i,k))

  head = "#165,80 2,1\n>101\t102\n{0}{1}\n".format(CD20,str(j))
  tail = ">102\n---------------------------------------------------------------------------------------------------------------------------------------------------------------------{0}".format(str(j))

  with open("CD20.msa", "r") as f:
      msa = f.read()

  with open("{0}-{1}/seq_custom.a3m".format(i,k), "w") as f:
      f.write(head+msa+tail)

  f.close()
