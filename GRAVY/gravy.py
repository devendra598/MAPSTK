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
        try:
            a = i.split("|")
            acc.append(a[1][:15])
        except(IndexError):
            a = i.split("\n")
            acc.append(a[0][:15].replace(">",""))
    

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
    "A": 0.62, "R": -2.5, "N": -0.78, "D": -0.9, "C": 0.29,
    "E": -0.74, "Q": -0.85, "G": 0.48, "H": -0.40, "I": 1.4,
    "L": 1.1, "K": -1.5, "M": 0.64, "F": 1.2, "P": 0.12,
    "S": -0.18, "T": -0.05, "W": 0.81, "Y": 0.26, "V": 1.1
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
    plt.ylabel("Accession ID")
    plt.savefig(f"./Gravy_output/gravy.pdf",format="pdf",bbox_inches="tight")
    

else:
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

print("The outputs are saved in Gravy_output folder")
    
