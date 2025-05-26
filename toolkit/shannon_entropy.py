import os
import math
import argparse
import warnings
from collections import Counter
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib import style
import numpy as np
from pyfamsa import Aligner, Sequence
import time
warnings.filterwarnings('ignore')
parser = argparse.ArgumentParser(description='Please provide following arguments.') 
parser.add_argument("-i", "--input", type=str, required=True, help="Input: FASTA format file of protein sequence/sequences")
args = parser.parse_args()

pfile = args.input
if not os.path.exists('Entropy_output'):
    os.mkdir('Entropy_output')

name = pfile.replace(".fasta", "")

try:
    with open(pfile, 'r') as fh:
        data = fh.readlines()
except IOError:
    print("Unable to open the file. Try again.")
    exit()
err = open("./Entropy_output/error_output.txt","w")
if(">" not in data[0]):
    err.write("Invalid input\nMissing '>' from your fasta file\n")
    exit()
if(">" not in data[0]):
    err.write("Invalid input\nMissing '>' from your fasta file\n")
    exit()
acc=[]
for i in data:
    if(">" in i):
        a = i.replace(">","").replace("\n","").replace("|","_").replace(" ","_")
        if len(a)>30:
            acc.append(a[:30])
        else:
            acc.append(a)
    

pos=[]
for i in range(len(data)):
    if(">" in data[i]):
        pos.append(i)
pos.append(len(data))

seq=[]

for i in range(len(pos)-1):   
    for j in range(pos[i]+1,pos[i+1]):
        seq.append(data[j])
    seq.append("\t")

protein_list=[]
gr=[]
for i in seq:
    if i != "\t":
        gr.append(i)
    else:
        protein_list.append("".join(gr).replace(' ', '').replace('\n', ''))
        gr=[]
if gr:
    protein_list.append("".join(gr).replace(' ', '').replace('\n', ''))
fasta_seq = []
for n in range(len(acc)):
    fasta_seq.append(Sequence(acc[n].encode(), protein_list[n].encode()))
aligner = Aligner()
msa = aligner.align(fasta_seq)

whm = open(f"./Entropy_output/{name}_msa.fasta","w")

for seqq in msa:
    whm.write(f">{seqq.id.decode()}\n{seqq.sequence.decode()}\n")
whm.close()


try:
    with open(f"./Entropy_output/{name}_msa.fasta", 'r') as fh2:
        data2 = fh2.readlines()
except IOError:
    print("Unable to open the file. Try again.")
    exit()
acc2=[]
for i in data2:
    if(">" in i):
        a = i.replace(">","").replace("\n","")
        acc2.append(a)
    

pos2=[]
for i in range(len(data2)):
    if(">" in data2[i]):
        pos2.append(i)
pos2.append(len(data2))

seq2=[]

for i in range(len(pos2)-1):   
    for j in range(pos2[i]+1,pos2[i+1]):
        seq2.append(data2[j])
    seq2.append("\t")

protein_list2=[]
gr2=[]
for i in seq2:
    if i != "\t":
        gr2.append(i)
    else:
        protein_list2.append("".join(gr2).replace(' ', '').replace('\n', ''))
        gr2=[]
if gr2:
    protein_list2.append("".join(gr2).replace(' ', '').replace('\n', ''))

matrix_file = f"./Entropy_output/{name}_matrix.tsv"
with open(matrix_file, "w") as wh1:
    for seq in protein_list2:
        wh1.write("\t".join(seq) + "\n")


with open(matrix_file, 'r') as fh1:
    rows = [line.strip().split('\t') for line in fh1 if line.strip()]


columns = list(map(list, zip(*rows)))

aa_single_here = "ACDEFGHIKLMNPQRSTVWY-"
entro = []
consensus = []
calc = 0
x_ax = []
output_file = "./Entropy_output/Shannon_entropy.tsv"
with open(output_file, "w") as wh:
    for column in columns:
        calc = calc + 1
        x_ax.append(calc)
        counter = Counter(column)
        most_common_element = counter.most_common(1)[0][0]  
        total = len(column)

        entropy = sum(
            (-count/total) * math.log(count/total) for aa, count in counter.items() if count > 0
        )
        
        wh.write(f"{most_common_element}\t{entropy}\n")
        consensus.append(most_common_element)
        entro.append(entropy)

consensus_seq = ''.join(map(str, consensus))
wh2 = open("./Entropy_output/Consensus_Sequence.fasta","w")
wh2.write(f">Consensus_Sequence\n{consensus_seq}")

ax=plt.axes()
ax.axhline(0,color='black')

ax.plot(x_ax,entro,color="blue")

plt.title("Graphical Representation of Entropy Value of Protein Sequences")
plt.xlabel(f"Sequence Position")
plt.ylabel("Entropy Value")
plt.savefig(f"./Entropy_output/Shannon_entropy.svg",bbox_inches="tight")