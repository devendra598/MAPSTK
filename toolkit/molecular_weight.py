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

if not os.path.exists('Molecular_Weight_output'):
        os.mkdir('Molecular_Weight_output')
mw = {
      'A': 71.08, 'R': 156.20, 'N': 114.11, 'D': 115.09, 'C': 103.14,
      'E': 129.12, 'Q': 128.41, 'G': 57.06, 'H': 137.15, 'I': 113.17,
      'L': 113.17, 'K': 128.1741, 'M': 131.21, 'F': 147.18, 'P': 97.12,
      'S': 87.08, 'T': 101.11, 'W': 186.21, 'Y': 163.18, 'V': 99.14
    }

molwt=[]

if(len(acc)<=2000):
    with open("./Molecular_Weight_output/Molecular_weight.tsv","w") as wh:
        wh.write("Accession_ID\tMol_Wt\n")
        for ac in range(len(acc)):
            plen=len(protein_list[ac])
            add=0
            for b in range(0,plen):
                base=protein_list[ac][b]
                try:
                    add=add+mw[f"{base}"]
                except(KeyError):
                    add=add

            mwprot=add+18.01524
            molwt.append(mwprot)
            wh.write(f"{acc[ac]}\t{mwprot/1000:.2f}\n")

   

    tick_spicing=1
    if len(acc)<20:
        fig = plt.figure(figsize=(5,5))
    else:
        fig = plt.figure(figsize=(len(acc)/5,5))
    # fig.patch.set_facecolor('xkcd:dark')
    ax=plt.axes()
    ax.axhline(0,color='black')
    
    colors=np.random.rand(len(acc),3)
    
    ax.bar(acc,molwt,color=colors) 
    plt.xticks(acc,rotation=90)
    # plt.yticks(color='white')
    ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_spicing))

    plt.title("Graphical Representation of Molecular Weights of Protein Sequences")
    plt.xlabel(f"FASTA Idnetifier")
    plt.ylabel("Molecular Weights")
    plt.savefig(f"./Molecular_Weight_output/molecular_weight.svg",bbox_inches="tight")
   
else:
    with open("./Molecular_Weight_output/Molecular_weight.tsv","w") as wh:
        wh.write("FASTA Idnetifier\tAA_length\tMol_Wt\n")
        for ac in range(len(acc)):
            plen=len(protein_list[ac])
            add=0
            for b in range(0,plen):
                base=protein_list[ac][b]
                add=add+mw[f"{base}"]

            mwprot=add+18.01524
            wh.write(f"{acc[ac]}\t{plen}\t{mwprot}\n")

