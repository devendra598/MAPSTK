import os
import math
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

mw = {
      'A': 71.0788, 'R': 156.1875, 'N': 114.1038, 'D': 115.0886, 'C': 103.1388,
      'E': 129.115, 'Q': 128.1307, 'G': 57.0519, 'H': 137.1411, 'I': 113.1594,
      'L': 113.1594, 'K': 128.1741, 'M': 131.1926, 'F': 147.1766, 'P': 97.1167,
      'S': 87.0782, 'T': 101.1051, 'W': 186.2132, 'Y': 163.1760, 'V': 99.1326
    }

if not os.path.exists('Entropy_output'):
        os.mkdir('Entropy_output')

entro=[]
if(len(acc)<=2000):
    with open("./Entropy_output/Shannon_Entropy.tsv","w") as wh:
        
        wh.write("Accession_ID\tEntropy\n")
        for ac in range(len(acc)):
            plen=len(protein_list[ac])
            aa_single_here	= ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"]

            a=c=d=e=f=g=h=i=k=l=m=n=p=q=r=s=t=v=w=y=u=0
            for aa in protein_list[ac]:
                if(aa == 'A'):
                    a=a+1
                elif(aa == 'R'):
                    r=r+1
                elif(aa == 'N'):
                    n=n+1
                elif(aa == 'D'):
                    d=d+1
                elif(aa == 'C'):
                    c=c+1
                elif(aa == 'Q'):
                    q=q+1
                elif(aa == 'E'):
                    e=e+1
                elif(aa == 'G'):
                    g=g+1
                elif(aa == 'H'):
                    h=h+1
                elif(aa == 'I'):
                    i=i+1
                elif(aa == 'L'):
                    l=l+1
                elif(aa == 'K'):
                    k=k+1
                elif(aa == 'M'):
                    m=m+1
                elif(aa == 'F'):
                    f=f+1
                elif(aa == 'S'):
                    s=s+1
                elif(aa == 'T'):
                    t=t+1
                elif(aa == 'W'):
                    w=w+1
                elif(aa == 'Y'):
                    y=y+1
                elif(aa == 'V'):
                    v=v+1
                elif(aa == 'P'):
                    p=p+1
                else:
                    u=u+1


            data={f'A':a,f'C':c,f'D':d,f'E':e,f'F':f,f'G':g,f'H':h,f'I':i,f'K':k,f'L':l,f'M':m,f'N':n,f'P':p,f'Q':q,f'R':r,f'S':s,f'T':t,f'V':v,f'W':w,f'Y':y}

            entropy=0
            
            for aa in aa_single_here:
                if data[aa] != 0:
                    fraction=data[aa]/plen
                    ent	= -fraction* (math.log(fraction) / math.log(2.0))
                    entropy = entropy + ent
            wh.write(f"{acc[ac]}\t{entropy}\n")
            entro.append(entropy)
        

        tick_spicing=1
        if len(acc)<20:
            fig = plt.figure(figsize=(5,5))
        else:
            fig = plt.figure(figsize=(len(acc)/5,5))

        ax=plt.axes()
        ax.axhline(0,color='black')

        
        colors=np.random.rand(len(acc),3)
        
        ax.bar(acc,entro,color=colors)
        plt.xticks(acc,rotation=90)

        ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_spicing))
        
        plt.title("Graphical Representation of Entropy Value of Protein Sequences")
        plt.xlabel(f"Accession ID")
        plt.ylabel("Entropy Value")
        plt.savefig(f"./Entropy_output/Shannon_entropy.pdf",format="pdf",bbox_inches="tight")

else:
    with open("./Entropy_output/Shannon_Entropy.tsv","w") as wh:
        
        wh.write("Accession_ID\tEntropy\n")
        for ac in range(len(acc)):
            plen=len(protein_list[ac])
            aa_single_here	= ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"]

            a=c=d=e=f=g=h=i=k=l=m=n=p=q=r=s=t=v=w=y=u=0
            for aa in protein_list[ac]:
                if(aa == 'A'):
                    a=a+1
                elif(aa == 'R'):
                    r=r+1
                elif(aa == 'N'):
                    n=n+1
                elif(aa == 'D'):
                    d=d+1
                elif(aa == 'C'):
                    c=c+1
                elif(aa == 'Q'):
                    q=q+1
                elif(aa == 'E'):
                    e=e+1
                elif(aa == 'G'):
                    g=g+1
                elif(aa == 'H'):
                    h=h+1
                elif(aa == 'I'):
                    i=i+1
                elif(aa == 'L'):
                    l=l+1
                elif(aa == 'K'):
                    k=k+1
                elif(aa == 'M'):
                    m=m+1
                elif(aa == 'F'):
                    f=f+1
                elif(aa == 'S'):
                    s=s+1
                elif(aa == 'T'):
                    t=t+1
                elif(aa == 'W'):
                    w=w+1
                elif(aa == 'Y'):
                    y=y+1
                elif(aa == 'V'):
                    v=v+1
                elif(aa == 'P'):
                    p=p+1
                else:
                    u=u+1


            data={f'A':a,f'C':c,f'D':d,f'E':e,f'F':f,f'G':g,f'H':h,f'I':i,f'K':k,f'L':l,f'M':m,f'N':n,f'P':p,f'Q':q,f'R':r,f'S':s,f'T':t,f'V':v,f'W':w,f'Y':y}

            entropy=0
            
            for aa in aa_single_here:
                if data[aa] != 0:
                    fraction=data[aa]/plen
                    ent	= -fraction* (math.log(fraction) / math.log(2.0))
                    entropy = entropy + ent
            wh.write(f"{acc[ac]}\t{entropy}\n")

print("The entropy output is saved in Entropy_output folder")               
            
            
