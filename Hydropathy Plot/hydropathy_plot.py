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
    "A": 0.62, "R": -2.5, "N": -0.78, "D": -0.9, "C": 0.29,
    "E": -0.74, "Q": -0.85, "G": 0.48, "H": -0.40, "I": 1.4,
    "L": 1.1, "K": -1.5, "M": 0.64, "F": 1.2, "P": 0.12,
    "S": -0.18, "T": -0.05, "W": 0.81, "Y": 0.26, "V": 1.1
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
        plt.title("Hydropathy Graph of Protein",fontsize=50)
        plt.xlabel("Positions of Amino Acids",fontsize=40)
        plt.ylabel("Hydropathy Scale",fontsize=40)
        plt.xticks(xax,rotation=90)
        if not os.path.exists("hydropathy_output"):
            os.mkdir("hydropathy_output")
        plt.savefig(f"./hydropathy_output/{acc[ac]}.pdf",format="pdf",bbox_inches="tight")
        plt.clf()
    with open(f"./hydropathy_output/{acc[ac]}_hydropathy.tsv","w") as wh:
        wh.write("Amino_Acids_w_positions\tHydropathy_values\n")
        for i in range(len(aa)):
            wh.write(f"{xax[i]}\t{hydro[i]}\n")


print(f"Output tsv files and figures are saved in hydropathy_output folder as their accession_id\nexample: {acc[0]}.pdf\t{acc[0]}.tsv")

