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


if not os.path.exists('Instability_Index_output'):
        os.mkdir('Instability_Index_output')

yax=[]
unyax=[]
normal=[]

if(len(acc)<=2000):
    with open("./Instability_Index_output/Instability_index.tsv","w") as wh:
        wh.write("Accession_ID\tInstability_Index\tStability\n")
        for ac in range(len(acc)):
            plen=len(protein_list[ac])
            out=[]
            for i in range(0,plen-1):
                try:
                    if(protein_list[ac][i]=="W"):
                        value={ 'W': 1.0, 'C': 1.0, 'M': 24.68, 'H': 24.68, 'Y': 1.0, 'F': 1.0, 'Q': 1.0, 'N': 13.34, 'I': 1.0, 'R': 1.0, 'D': 1.0, 'P': 1.0, 'T': -14.03, 'K': 1.0, 'E': 1.0, 'V': -7.49, 'S': 1.0, 'G': -9.37, 'A': -14.03, 'L': 13.34 }
                        out.append(value[f"{protein_list[ac][i+1]}"])
                    elif(protein_list[ac][i]=="C"):
                        value={ 'W': 24.68, 'C': 1.0, 'M': 33.6, 'H': 33.6, 'Y': 1.0, 'F': 1.0, 'Q': -6.54, 'N': 1.0, 'I': 1.0, 'R': 1.0, 'D': 20.26, 'P': 20.26, 'T': 33.6, 'K': 1.0, 'E': 1.0, 'V': -6.54, 'S': 1.0, 'G': 1.0, 'A': 1.0, 'L': 20.26 }
                        out.append(value[f"{protein_list[ac][i+1]}"])
                    elif(protein_list[ac][i]=="M"):
                        value={'W':1.0,'C':1.0,'M':-1.88,'H':58.28,'Y':24.68,'F':1.0,'Q':-6.54,'N':1.0,'I':1.0,'R':-6.54,'D':1.0,'P':44.94,'T':-1.88,'K':1.0,'E':1.0,'V':1.0,'S':44.94,'G':1.0,'A':13.34,'L':1.0}
                        out.append(value[f"{protein_list[ac][i+1]}"])
                    elif(protein_list[ac][i]=="H"):
                        value={'W':-1.88,'C':1.0,'M':1.0,'H':1.0,'Y':44.94,'F':-9.37,'Q':1.0,'N':24.68,'I':44.94,'R':1.0,'D':1.0,'P':-1.88,'T':-6.54,'K':24.68,'E':1.0,'V':1.0,'S':1.0,'G':-9.37,'A':1.0,'L':1.0}
                        out.append(value[f"{protein_list[ac][i+1]}"])
                    elif(protein_list[ac][i]=="Y"):
                        value={'W':-9.37,'C':1.0,'M':44.94,'H':13.34,'Y':13.34,'F':1.0,'Q':1.0,'N':1.0,'I':1.0,'R':-15.91,'D':24.68,'P':13.34,'T':-7.49,'K':1.0,'E':-6.54,'V':1.0,'S':1.0,'G':-7.49,'A':24.68,'L':1.0}
                        out.append(value[f"{protein_list[ac][i+1]}"])
                    elif(protein_list[ac][i]=="F"):
                        value={'W':1.0,'C':1.0,'M':1.0,'H':1.0,'Y':33.6,'F':1.0,'Q':1.0,'N':1.0,'I':1.0,'R':1.0,'D':13.34,'P':20.26,'T':1.0,'K':-14.03,'E':1.0,'V':1.0,'S':1.0,'G':1.0,'A':1.0,'L':1.0}
                        out.append(value[f"{protein_list[ac][i+1]}"])
                    elif(protein_list[ac][i]=="Q"):
                        value={'W':1.0,'C':-6.54,'M':1.0,'H':1.0,'Y':-6.54,'F':-6.54,'Q':20.26,'N':1.0,'I':1.0,'R':1.0,'D':20.26,'P':20.26,'T':1.0,'K':1.0,'E':20.26,'V':-6.54,'S':44.94,'G':1.0,'A':1.0,'L':1.0}
                        out.append(value[f"{protein_list[ac][i+1]}"])
                    elif(protein_list[ac][i]=="N"):
                        value={'W':-9.37,'C':-1.88,'M':1.0,'H':1.0,'Y':1.0,'F':-14.03,'Q':-6.54,'N':1.0,'I':44.94,'R':1.0,'D':1.0,'P':-1.88,'T':-7.49,'K':24.68,'E':1.0,'V':1.0,'S':1.0,'G':-14.03,'A':1.0,'L':1.0}
                        out.append(value[f"{protein_list[ac][i+1]}"])
                    elif(protein_list[ac][i]=="I"):
                        value={'W':1.0,'C':1.0,'M':1.0,'H':13.34,'Y':1.0,'F':1.0,'Q':1.0,'N':1.0,'I':1.0,'R':1.0,'D':1.0,'P':-1.88,'T':1.0,'K':-7.49,'E':44.94,'V':-7.49,'S':1.0,'G':1.0,'A':1.0,'L':20.26}
                        out.append(value[f"{protein_list[ac][i+1]}"])
                    elif(protein_list[ac][i]=="R"):
                        value={'W':58.28,'C':1.0,'M':1.0,'H':20.26,'Y':-6.54,'F':1.0,'Q':20.26,'N':13.34,'I':1.0,'R':58.28,'D':1.0,'P':20.26,'T':1.0,'K':1.0,'E':1.0,'V':1.0,'S':44.94,'G':-7.49,'A':1.0,'L':1.0}
                        out.append(value[f"{protein_list[ac][i+1]}"])
                    elif(protein_list[ac][i]=="D"):
                        value={'W':1.0,'C':1.0,'M':1.0,'H':1.0,'Y':1.0,'F':-6.54,'Q':1.0,'N':1.0,'I':1.0,'R':-6.54,'D':1.0,'P':1.0,'T':-14.03,'K':-7.49,'E':1.0,'V':1.0,'S':20.26,'G':1.0,'A':1.0,'L':1.0}
                        out.append(value[f"{protein_list[ac][i+1]}"])
                    elif(protein_list[ac][i]=="P"):
                        value={'W':-1.88,'C':-6.54,'M':-6.54,'H':1.0,'Y':1.0,'F':20.26,'Q':20.26,'N':1.0,'I':1.0,'R':-6.54,'D':-6.54,'P':20.26,'T':1.0,'K':1.0,'E':18.38,'V':20.26,'S':20.26,'G':1.0,'A':20.26,'L':1.0}
                        out.append(value[f"{protein_list[ac][i+1]}"])
                    elif(protein_list[ac][i]=="T"):
                        value={'W':-14.03,'C':1.0,'M':1.0,'H':1.0,'Y':1.0,'F':13.34,'Q':-6.54,'N':-14.03,'I':1.0,'R':1.0,'D':1.0,'P':1.0,'T':1.0,'K':1.0,'E':20.26,'V':1.0,'S':1.0,'G':-7.49,'A':1.0,'L':1.0}
                        out.append(value[f"{protein_list[ac][i+1]}"])
                    elif(protein_list[ac][i]=="K"):
                        value={'W':1.0,'C':1.0,'M':33.6,'H':1.0,'Y':1.0,'F':1.0,'Q':24.68,'N':1.0,'I':-7.49,'R':33.6,'D':1.0,'P':-6.54,'T':1.0,'K':1.0,'E':1.0,'V':-7.49,'S':1.0,'G':-7.49,'A':1.0,'L':-7.49}
                        out.append(value[f"{protein_list[ac][i+1]}"])
                    elif(protein_list[ac][i]=="E"):
                        value={'W':-14.03,'C':44.94,'M':1.0,'H':-6.54,'Y':1.0,'F':1.0,'Q':20.26,'N':1.0,'I':20.26,'R':1.0,'D':20.26,'P':20.26,'T':1.0,'K':1.0,'E':33.6,'V':1.0,'S':20.26,'G':1.0,'A':1.0,'L':1.0}
                        out.append(value[f"{protein_list[ac][i+1]}"])
                    elif(protein_list[ac][i]=="V"):
                        value={'W':1.0,'C':1.0,'M':1.0,'H':1.0,'Y':-6.54,'F':1.0,'Q':1.0,'N':1.0,'I':1.0,'R':1.0,'D':-14.03,'P':20.26,'T':-7.49,'K':-1.88,'E':1.0,'V':1.0,'S':1.0,'G':-7.49,'A':1.0,'L':1.0}
                        out.append(value[f"{protein_list[ac][i+1]}"])
                    elif(protein_list[ac][i]=="S"):
                        value={'W':1.0,'C':33.6,'M':1.0,'H':1.0,'Y':1.0,'F':1.0,'Q':20.26,'N':1.0,'I':1.0,'R':20.26,'D':1.0,'P':44.94,'T':1.0,'K':1.0,'E':20.26,'V':1.0,'S':20.26,'G':1.0,'A':1.0,'L':1.0}
                        out.append(value[f"{protein_list[ac][i+1]}"])
                    elif(protein_list[ac][i]=="G"):
                        value={'W':13.34,'C':1.0,'M':1.0,'H':1.0,'Y':-7.49,'F':1.0,'Q':1.0,'N':-7.49,'I':-7.49,'R':1.0,'D':1.0,'P':1.0,'T':-7.49,'K':-7.49,'E':-6.54,'V':1.0,'S':1.0,'G':13.34,'A':-7.49,'L':1.0}
                        out.append(value[f"{protein_list[ac][i+1]}"])
                    elif(protein_list[ac][i]=="A"):
                        value={'W':1.0,'C':44.94,'M':1.0,'H':-7.49,'Y':1.0,'F':1.0,'Q':1.0,'N':1.0,'I':1.0, 'R':1.0,'D':-7.49,'P':20.26,'T':1.0,'K':1.0,'E':1.0,'V':1.0,'S':1.0,'G':1.0,'A':1.0,'L':1.0}
                        out.append(value[f"{protein_list[ac][i+1]}"])
                    elif(protein_list[ac][i]=="L"):
                        value={'W':24.68,'C':1.0,'M':1.0,'H':1.0,'Y':1.0,'F':1.0,'Q':33.6,'N':1.0,'I':1.0,'R':20.26,'D':1.0,'P':20.26,'T':1.0,'K':-7.49,'E':1.0,'V':1.0,'S':1.0,'G':1.0,'A':1.0,'L':1.0}
                        out.append(value[f"{protein_list[ac][i+1]}"])
                except(KeyError):
                    out.append(0)
                    
            add=sum(out)

            II=(10/plen)*add
            
            
            if(II<40):
                yax.append(II)
                unyax.append(0)
                normal.append(0)
                wh.write(f"{acc[ac]}\t{II}\tStable\n")
            else:
                wh.write(f"{acc[ac]}\t{II}\tUnstable\n")
                normal.append(40)
                yax.append(40)
                unyax.append(II)
    
    
    if len(acc)<20:
        plt.figure(figsize=(5,5))
    else:
        plt.figure(figsize=(len(acc)/5,5))
    tick_spicing=1
    
    ax=plt.axes()
    ax.axhline(0,color='black')
    ax.axhline(40,color='indianred')
    
    colors=np.random.rand(len(acc),3)
    
    cl = ['orangered','chocolate','teal']
    labels = ['Amount of Instability','Unstable Protein','Stable Protein']
    ax.bar(acc,unyax,color="orangered") 
    ax.bar(acc,yax,color="teal")
    ax.bar(acc,normal,color="chocolate")
    

    patches = [plt.Rectangle((0,0),1,1,fc=color, edgecolor='none') for color in cl]
    plt.legend(patches, labels, loc=(1,1.1))
    
    
    plt.xticks(acc,rotation=90)
    
    ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_spicing))

    plt.title("Graphical Representation of Stability of Protein Sequences")
    plt.xlabel(f"FASTA Idnetifier")
    plt.ylabel("Instability Index Value")
    
    plt.savefig(f"./Instability_Index_output/instability_index.svg",bbox_inches="tight")
    
else:
    with open("./Instability_Index_output/Instability_index.tsv","w") as wh:
        wh.write("FASTA Idnetifier\tInstability_Index\tStability\n")
        for ac in range(len(acc)):
            plen=len(protein_list[ac])
            out=[]
            for i in range(0,plen-1):
                try:
                    if(protein_list[ac][i]=="W"):
                        value={ 'W': 1.0, 'C': 1.0, 'M': 24.68, 'H': 24.68, 'Y': 1.0, 'F': 1.0, 'Q': 1.0, 'N': 13.34, 'I': 1.0, 'R': 1.0, 'D': 1.0, 'P': 1.0, 'T': -14.03, 'K': 1.0, 'E': 1.0, 'V': -7.49, 'S': 1.0, 'G': -9.37, 'A': -14.03, 'L': 13.34 }
                        out.append(value[f"{protein_list[ac][i+1]}"])
                    elif(protein_list[ac][i]=="C"):
                        value={ 'W': 24.68, 'C': 1.0, 'M': 33.6, 'H': 33.6, 'Y': 1.0, 'F': 1.0, 'Q': -6.54, 'N': 1.0, 'I': 1.0, 'R': 1.0, 'D': 20.26, 'P': 20.26, 'T': 33.6, 'K': 1.0, 'E': 1.0, 'V': -6.54, 'S': 1.0, 'G': 1.0, 'A': 1.0, 'L': 20.26 }
                        out.append(value[f"{protein_list[ac][i+1]}"])
                    elif(protein_list[ac][i]=="M"):
                        value={'W':1.0,'C':1.0,'M':-1.88,'H':58.28,'Y':24.68,'F':1.0,'Q':-6.54,'N':1.0,'I':1.0,'R':-6.54,'D':1.0,'P':44.94,'T':-1.88,'K':1.0,'E':1.0,'V':1.0,'S':44.94,'G':1.0,'A':13.34,'L':1.0}
                        out.append(value[f"{protein_list[ac][i+1]}"])
                    elif(protein_list[ac][i]=="H"):
                        value={'W':-1.88,'C':1.0,'M':1.0,'H':1.0,'Y':44.94,'F':-9.37,'Q':1.0,'N':24.68,'I':44.94,'R':1.0,'D':1.0,'P':-1.88,'T':-6.54,'K':24.68,'E':1.0,'V':1.0,'S':1.0,'G':-9.37,'A':1.0,'L':1.0}
                        out.append(value[f"{protein_list[ac][i+1]}"])
                    elif(protein_list[ac][i]=="Y"):
                        value={'W':-9.37,'C':1.0,'M':44.94,'H':13.34,'Y':13.34,'F':1.0,'Q':1.0,'N':1.0,'I':1.0,'R':-15.91,'D':24.68,'P':13.34,'T':-7.49,'K':1.0,'E':-6.54,'V':1.0,'S':1.0,'G':-7.49,'A':24.68,'L':1.0}
                        out.append(value[f"{protein_list[ac][i+1]}"])
                    elif(protein_list[ac][i]=="F"):
                        value={'W':1.0,'C':1.0,'M':1.0,'H':1.0,'Y':33.6,'F':1.0,'Q':1.0,'N':1.0,'I':1.0,'R':1.0,'D':13.34,'P':20.26,'T':1.0,'K':-14.03,'E':1.0,'V':1.0,'S':1.0,'G':1.0,'A':1.0,'L':1.0}
                        out.append(value[f"{protein_list[ac][i+1]}"])
                    elif(protein_list[ac][i]=="Q"):
                        value={'W':1.0,'C':-6.54,'M':1.0,'H':1.0,'Y':-6.54,'F':-6.54,'Q':20.26,'N':1.0,'I':1.0,'R':1.0,'D':20.26,'P':20.26,'T':1.0,'K':1.0,'E':20.26,'V':-6.54,'S':44.94,'G':1.0,'A':1.0,'L':1.0}
                        out.append(value[f"{protein_list[ac][i+1]}"])
                    elif(protein_list[ac][i]=="N"):
                        value={'W':-9.37,'C':-1.88,'M':1.0,'H':1.0,'Y':1.0,'F':-14.03,'Q':-6.54,'N':1.0,'I':44.94,'R':1.0,'D':1.0,'P':-1.88,'T':-7.49,'K':24.68,'E':1.0,'V':1.0,'S':1.0,'G':-14.03,'A':1.0,'L':1.0}
                        out.append(value[f"{protein_list[ac][i+1]}"])
                    elif(protein_list[ac][i]=="I"):
                        value={'W':1.0,'C':1.0,'M':1.0,'H':13.34,'Y':1.0,'F':1.0,'Q':1.0,'N':1.0,'I':1.0,'R':1.0,'D':1.0,'P':-1.88,'T':1.0,'K':-7.49,'E':44.94,'V':-7.49,'S':1.0,'G':1.0,'A':1.0,'L':20.26}
                        out.append(value[f"{protein_list[ac][i+1]}"])
                    elif(protein_list[ac][i]=="R"):
                        value={'W':58.28,'C':1.0,'M':1.0,'H':20.26,'Y':-6.54,'F':1.0,'Q':20.26,'N':13.34,'I':1.0,'R':58.28,'D':1.0,'P':20.26,'T':1.0,'K':1.0,'E':1.0,'V':1.0,'S':44.94,'G':-7.49,'A':1.0,'L':1.0}
                        out.append(value[f"{protein_list[ac][i+1]}"])
                    elif(protein_list[ac][i]=="D"):
                        value={'W':1.0,'C':1.0,'M':1.0,'H':1.0,'Y':1.0,'F':-6.54,'Q':1.0,'N':1.0,'I':1.0,'R':-6.54,'D':1.0,'P':1.0,'T':-14.03,'K':-7.49,'E':1.0,'V':1.0,'S':20.26,'G':1.0,'A':1.0,'L':1.0}
                        out.append(value[f"{protein_list[ac][i+1]}"])
                    elif(protein_list[ac][i]=="P"):
                        value={'W':-1.88,'C':-6.54,'M':-6.54,'H':1.0,'Y':1.0,'F':20.26,'Q':20.26,'N':1.0,'I':1.0,'R':-6.54,'D':-6.54,'P':20.26,'T':1.0,'K':1.0,'E':18.38,'V':20.26,'S':20.26,'G':1.0,'A':20.26,'L':1.0}
                        out.append(value[f"{protein_list[ac][i+1]}"])
                    elif(protein_list[ac][i]=="T"):
                        value={'W':-14.03,'C':1.0,'M':1.0,'H':1.0,'Y':1.0,'F':13.34,'Q':-6.54,'N':-14.03,'I':1.0,'R':1.0,'D':1.0,'P':1.0,'T':1.0,'K':1.0,'E':20.26,'V':1.0,'S':1.0,'G':-7.49,'A':1.0,'L':1.0}
                        out.append(value[f"{protein_list[ac][i+1]}"])
                    elif(protein_list[ac][i]=="K"):
                        value={'W':1.0,'C':1.0,'M':33.6,'H':1.0,'Y':1.0,'F':1.0,'Q':24.68,'N':1.0,'I':-7.49,'R':33.6,'D':1.0,'P':-6.54,'T':1.0,'K':1.0,'E':1.0,'V':-7.49,'S':1.0,'G':-7.49,'A':1.0,'L':-7.49}
                        out.append(value[f"{protein_list[ac][i+1]}"])
                    elif(protein_list[ac][i]=="E"):
                        value={'W':-14.03,'C':44.94,'M':1.0,'H':-6.54,'Y':1.0,'F':1.0,'Q':20.26,'N':1.0,'I':20.26,'R':1.0,'D':20.26,'P':20.26,'T':1.0,'K':1.0,'E':33.6,'V':1.0,'S':20.26,'G':1.0,'A':1.0,'L':1.0}
                        out.append(value[f"{protein_list[ac][i+1]}"])
                    elif(protein_list[ac][i]=="V"):
                        value={'W':1.0,'C':1.0,'M':1.0,'H':1.0,'Y':-6.54,'F':1.0,'Q':1.0,'N':1.0,'I':1.0,'R':1.0,'D':-14.03,'P':20.26,'T':-7.49,'K':-1.88,'E':1.0,'V':1.0,'S':1.0,'G':-7.49,'A':1.0,'L':1.0}
                        out.append(value[f"{protein_list[ac][i+1]}"])
                    elif(protein_list[ac][i]=="S"):
                        value={'W':1.0,'C':33.6,'M':1.0,'H':1.0,'Y':1.0,'F':1.0,'Q':20.26,'N':1.0,'I':1.0,'R':20.26,'D':1.0,'P':44.94,'T':1.0,'K':1.0,'E':20.26,'V':1.0,'S':20.26,'G':1.0,'A':1.0,'L':1.0}
                        out.append(value[f"{protein_list[ac][i+1]}"])
                    elif(protein_list[ac][i]=="G"):
                        value={'W':13.34,'C':1.0,'M':1.0,'H':1.0,'Y':-7.49,'F':1.0,'Q':1.0,'N':-7.49,'I':-7.49,'R':1.0,'D':1.0,'P':1.0,'T':-7.49,'K':-7.49,'E':-6.54,'V':1.0,'S':1.0,'G':13.34,'A':-7.49,'L':1.0}
                        out.append(value[f"{protein_list[ac][i+1]}"])
                    elif(protein_list[ac][i]=="A"):
                        value={'W':1.0,'C':44.94,'M':1.0,'H':-7.49,'Y':1.0,'F':1.0,'Q':1.0,'N':1.0,'I':1.0, 'R':1.0,'D':-7.49,'P':20.26,'T':1.0,'K':1.0,'E':1.0,'V':1.0,'S':1.0,'G':1.0,'A':1.0,'L':1.0}
                        out.append(value[f"{protein_list[ac][i+1]}"])
                    elif(protein_list[ac][i]=="L"):
                        value={'W':24.68,'C':1.0,'M':1.0,'H':1.0,'Y':1.0,'F':1.0,'Q':33.6,'N':1.0,'I':1.0,'R':20.26,'D':1.0,'P':20.26,'T':1.0,'K':-7.49,'E':1.0,'V':1.0,'S':1.0,'G':1.0,'A':1.0,'L':1.0}
                        out.append(value[f"{protein_list[ac][i+1]}"])
                except(KeyError):
                    out.append(0)
            add=sum(out)

            II=(10/plen)*add
            
            if(II<40):
                wh.write(f"{acc[ac]}\t{II}\tStable\n")
            else:
                wh.write(f"{acc[ac]}\t{II}\tUnstable\n")

