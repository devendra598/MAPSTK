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

if not os.path.exists('Theoritical_pI_output'):
        os.mkdir('Theoritical_pI_output')
thpi=[]

if(len(acc)<=2000):
    with open("./Theoritical_pI_output/theoritical_pi.tsv","w") as wh:
        wh.write("Accession_ID\tTheoritical_pI\n")
        for ac in range(len(acc)):
            plen=len(protein_list[ac])

            d=e=c=y=h=k=r=u=0
            for aa in protein_list[ac]:
                if(aa == 'D'):
                    d=d+1
                elif(aa == 'E'):
                    e=e+1
                elif(aa == 'C'):
                    c=c+1
                elif(aa == 'Y'):
                    y=y+1
                elif(aa == 'H'):
                    h=h+1
                elif(aa == 'K'):
                    k=k+1
                elif(aa == 'R'):
                    r=r+1
                else:
                    u=u+1

            ph=0
            end=0
            eqd=eqe=eqc=eqy=eqh=eqk=eqr=eqcooh=eqnh=0.0
            while(end==0):
                eqd=-d/(1+(10**(3.65-ph)))
                eqe=-e/(1+(10**(4.25-ph)))
                eqc=-c/(1+(10**(8.18-ph)))
                eqy=-y/(1+(10**(10.7-ph)))
                eqh=h/(1+(10**(ph-6.0)))
                eqk=k/(1+(10**(ph-10.53)))
                eqr=r/(1+(10**(ph-12.48)))
                eqcooh=-1/(1+(10**(2.34-ph)))
                eqnh=1/(1+(10**(ph-9.69)))
                tpi=eqd+eqe+eqc+eqy+eqh+eqk+eqr+eqcooh+eqnh
                if(tpi>0):
                    ph=ph+0.01
                else:
                    end+=1

            wh.write(f"{acc[ac]}\t{ph}\n")
            thpi.append(ph)

    if len(acc)<20:
        plt.figure(figsize=(5,5))
    else:
        plt.figure(figsize=(len(acc)/5,5))
    tick_spicing=1
    
    ax=plt.axes()
    ax.axhline(0,color='black')
    ax.axhline(7,color='seagreen')
    
    acid=[]
    base=[]
    neutral=[]
    for i in range(len(thpi)):
        if(thpi[i]>7):
            base.append(thpi[i])
            acid.append(0)
            neutral.append(0)
        elif(thpi[i]==7):
            base.append(0)
            acid.append(0)
            neutral.append(thpi[i])
        else:
            base.append(0)
            acid.append(thpi[i])
            neutral.append(0)
    
    cl = ['blue','red','green']
    labels = ['Basic','Acidic','Neutral']
    ax.bar(acc,base,color="blue")
    ax.bar(acc,acid,color="red")
    ax.bar(acc,neutral,color="green")
    plt.xticks(acc,rotation=90)
    patches = [plt.Rectangle((0,0),1,1,fc=color, edgecolor='none') for color in cl]
    plt.legend(patches, labels, loc=(1,1.1))

    ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_spicing))

    plt.title("Graphical Representation of Theoretical Isoelectric Point of Protein Sequences")
    plt.xlabel(f"FASTA Idnetifier")
    plt.ylabel("Isoelectric Point")
    plt.savefig(f"./Theoritical_pI_output/theoritical_pI.svg",bbox_inches="tight")
   

else:
    with open("./Theoritical_pI_output/theoritical_pi.tsv","w") as wh:
        wh.write("FASTA Idnetifier\tTheoritical_pI\n")
        for ac in range(len(acc)):
            plen=len(protein_list[ac])

            d=e=c=y=h=k=r=u=0
            for aa in protein_list[ac]:
                if(aa == 'D'):
                    d=d+1
                elif(aa == 'E'):
                    e=e+1
                elif(aa == 'C'):
                    c=c+1
                elif(aa == 'Y'):
                    y=y+1
                elif(aa == 'H'):
                    h=h+1
                elif(aa == 'K'):
                    k=k+1
                elif(aa == 'R'):
                    r=r+1
                else:
                    u=u+1

            ph=0
            end=0
            eqd=eqe=eqc=eqy=eqh=eqk=eqr=eqcooh=eqnh=0.0
            while(end==0):
                eqd=-d/(1+(10**(3.65-ph)))
                eqe=-e/(1+(10**(4.25-ph)))
                eqc=-c/(1+(10**(8.18-ph)))
                eqy=-y/(1+(10**(10.7-ph)))
                eqh=h/(1+(10**(ph-6.0)))
                eqk=k/(1+(10**(ph-10.53)))
                eqr=r/(1+(10**(ph-12.48)))
                eqcooh=-1/(1+(10**(2.34-ph)))
                eqnh=1/(1+(10**(ph-9.69)))
                tpi=eqd+eqe+eqc+eqy+eqh+eqk+eqr+eqcooh+eqnh
                if(tpi>0):
                    ph=ph+0.01
                else:
                    end+=1

            wh.write(f"{acc[ac]}\t{ph}\n")
            

