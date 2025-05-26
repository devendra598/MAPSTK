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

seq = []

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


if not os.path.exists('Aliphatic_Index_output'):
        os.mkdir('Aliphatic_Index_output')

ALI=[]


if(len(acc)<=2000):
    with open("./Aliphatic_Index_output/aliphatic_index.tsv","w") as wh:
        wh.write("Accession_ID\tAliphatic_Index\n")
        for j in range(len(acc)):
            plen=len(protein_list[j])
            a = protein_list[j].count("A")
            i = protein_list[j].count("I")
            l = protein_list[j].count("L")
            v = protein_list[j].count("V")

            MolPA=a/plen*100
            MolPI=i/plen*100
            MolPL=l/plen*100
            MolPV=v/plen*100
            rvV=2.9
            rvIL=3.9
            AI=MolPA + (rvV*MolPV) + (rvIL*(MolPI+MolPL))
            ALI.append(AI)
            wh.write(f"{acc[j]}\t{AI}\n")
            
    tick_spicing=1
    if len(acc)<20:
        fig = plt.figure(figsize=(5,5))
    else:
        fig = plt.figure(figsize=(len(acc)/5,5))
    
    ax=plt.axes()
    ax.axhline(0,color='black')
    
    colors=np.random.rand(len(acc),3)
    
    ax.bar(acc,ALI,color=colors) 
    plt.xticks(acc,rotation=90,color='black')
    plt.yticks(color='black')
    
    ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_spicing))

    plt.title("Graphical Representation of Aliphatic Index of Protein Sequences", color='black')
    plt.xlabel(f"FASTA Idnetifier")
    plt.ylabel("Aliphatic Index Value")
    plt.savefig(f"./Aliphatic_Index_output/aliphatic_index.svg",bbox_inches="tight")
   
else:
    with open("./Aliphatic_Index_output/aliphatic_index.tsv","w") as wh:
        wh.write("FASTA Idnetifier\tAliphatic_Index\n")
        for j in range(len(acc)):
            plen=len(protein_list[j])
            a = protein_list[j].count("A")
            i = protein_list[j].count("I")
            l = protein_list[j].count("L")
            v = protein_list[j].count("V")

            MolPA=a/plen*100
            MolPI=i/plen*100
            MolPL=l/plen*100
            MolPV=v/plen*100
            rvV=2.9
            rvIL=3.9
            AI=MolPA + (rvV*MolPV) + (rvIL*(MolPI+MolPL))
            wh.write(f"{acc[j]}\t{AI}\n")
    
