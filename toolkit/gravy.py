import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib import style
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
    print("No such file in the directory.")
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

hyd = {
    'A': 1.8,  # Alanine
    'R': -4.5, # Arginine
    'N': -3.5, # Asparagine
    'D': -3.5, # Aspartic acid
    'C': 2.5,  # Cysteine
    'Q': -3.5, # Glutamine
    'E': -3.5, # Glutamic acid
    'G': -0.4, # Glycine
    'H': -3.2, # Histidine
    'I': 4.5,  # Isoleucine
    'L': 3.8,  # Leucine
    'K': -3.9, # Lysine
    'M': 1.9,  # Methionine
    'F': 2.8,  # Phenylalanine
    'P': -1.6, # Proline
    'S': -0.8, # Serine
    'T': -0.7, # Threonine
    'W': -0.9, # Tryptophan
    'Y': -1.3, # Tyrosine
    'V': 4.2   # Valine
}

yax=[]
yax1=[]
if not os.path.exists('Gravy_output'):
        os.mkdir('Gravy_output')

if(len(acc)<=2000):
    with open("./Gravy_output/gravy.tsv","w") as wh:
        wh.write("Accession_ID\tgravy_value\n")
        for ac in range(len(acc)):
            plen=len(protein_list[ac])
            add=0
            
            for b in range(0,plen):
                base=protein_list[ac][b]
                try:
                    add=add+hyd[f"{base}"]
                except(KeyError):
                    add=add
            gravy=add/plen
            wh.write(f"{acc[ac]}\t{gravy}\n")
            if(gravy>0):
                yax.append(gravy)
                yax1.append(0)
            else:
                yax.append(0)
                yax1.append(gravy)

    
    if len(acc)<20:
        plt.figure(figsize=(5,5))
    else:
        plt.figure(figsize=(5,len(acc)/5))
    tick_spicing=0.1
    
    ax=plt.axes()
    ax.axvline(0,color='blue')
    
    ax.barh(acc,yax,color="green")
    ax.barh(acc,yax1,color="red")

    plt.xticks(yax,rotation=90)
    
    ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_spicing))

    plt.title("Graphical Representation of GRAVY Value of Protein Sequences")
    plt.xlabel("Gravy Value")
    plt.ylabel("FASTA Idnetifier")
    plt.savefig(f"./Gravy_output/gravy.svg",bbox_inches="tight")
    

else:
    with open("./Gravy_output/gravy.tsv","w") as wh:
        wh.write("FASTA Idnetifier\tgravy_value\n")
        for ac in range(len(acc)):
            plen=len(protein_list[ac])
            add=0
            
            for b in range(0,plen):
                base=protein_list[ac][b]
                try:
                    add=add+hyd[f"{base}"]
                except(KeyError):
                    add=add
            gravy=add/plen
            wh.write(f"{acc[ac]}\t{gravy}\n")

    
