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
# print(len(pos))
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


for ac in range(len(acc)):
    plen=len(protein_list[ac])
    xax=[]
    hypho=[]
    hyphi=[]
    aa=[]
    hydro=[]
    for b in range(0,plen):
        xax.append(f"{protein_list[ac][b]}{b+1}")
        base=protein_list[ac][b]
        aa.append(base)
        hydro.append(hyd[f"{base}"])
        if(hyd[f"{base}"]>0):
            hypho.append(hyd[f"{base}"])
            hyphi.append(0)
        else:
            hyphi.append(hyd[f"{base}"])
            hypho.append(0)
    if len(protein_list[ac])<=2000:
        tick_spicing=1
        
        
        plt.figure(figsize=(len(aa)/5,8))

        ax=plt.axes()
        ax.axhline(0,color='black')
        ax.plot(xax,hypho,color='blue')
        ax.fill_between(xax, hypho, color='lightblue', alpha=0.5)
        ax.plot(xax,hyphi,color='red')
        ax.fill_between(xax, hyphi, color='lightcoral', alpha=0.5)
        # ax.bar(ind,h,color='skyblue')
        ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_spicing))
        plt.margins(0)
        plt.title(f"Hydropathy Graph of Protein {acc[ac]}",fontsize=50)
        plt.xlabel("Positions of Amino Acids",fontsize=40)
        plt.ylabel("Hydropathy Scale",fontsize=40)
        plt.xticks(xax,rotation=90)
        if not os.path.exists("hydropathy_output"):
            os.mkdir("hydropathy_output")
        plt.savefig(f"./hydropathy_output/{acc[ac]}.svg",bbox_inches="tight")
        plt.clf()
    with open(f"./hydropathy_output/{acc[ac]}_hydropathy.tsv","w") as wh:
        wh.write("Amino_Acids_w_positions\tHydropathy_values\n")
        for i in range(len(aa)):
            wh.write(f"{xax[i]}\t{hydro[i]}\n")



