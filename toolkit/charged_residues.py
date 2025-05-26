import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib import style
import numpy as np
import os
import argparse  
import warnings
warnings.filterwarnings('ignore')
parser = argparse.ArgumentParser(description='Please provide following arguments.') 

parser.add_argument("-i", "--input", type=str, required=True, help="Input: FASTA format file of protein sequence/sequences")
args = parser.parse_args()

pfile=args.input

try:
    with open(pfile, 'r') as fh:
        data = fh.readlines()
except IOError:
    print("Unable to open the file. Try again.")
    exit()
if(">" not in data[0]):
        print("Invalid input\nMissing '>' from your fasta file")
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


for ac in range(len(acc)):
    plen=len(protein_list[ac])
    e=protein_list[ac].count("E")
    d=protein_list[ac].count("D")
    r=protein_list[ac].count("R")
    k=protein_list[ac].count("K")
    data={f"Total Positive Charged Residues Present is:[{r+k}]":r+k,f"Total Negative Charged Residues Present is:[{e+d}]":e+d}
    key=list(data.keys())
    value=list(data.values())

    explode = [0.009,0.009]
    plt.pie(value,colors="green,red",labels=key,explode=explode,radius = 1.2)
    plt.text(-2.5,-1, f"Total amount of charged residues of\n{acc[ac]} is {r+k+e+d}", fontsize=10 , ha='center')
    if not os.path.exists("neg_pos_output"):
        os.mkdir("neg_pos_output")
    plt.savefig(f"./neg_pos_output/{acc[ac]}.svg",bbox_inches="tight")
    plt.clf()

