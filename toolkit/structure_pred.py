import os
import argparse  
import warnings
import requests
from RamachanDraw import phi_psi, plot
import matplotlib.pyplot as plt
from modeller import *
import subprocess
import time
import gzip
import shutil
warnings.filterwarnings('ignore')
parser = argparse.ArgumentParser(description='Please provide following arguments.') 


parser.add_argument("-i", "--input", type=str, required=True, help="Input: FASTA format file of protein sequence/sequences")
parser.add_argument("-pi","--identity", type=float, default=90.0, help="Identity: Minimum percent of identity for finding templates in structure prediction Default = 90.0")
parser.add_argument("-e","--evalue", type=float, default=0.01, help="E-value: E-value for BLAST in structure prediction Default = 0.01")
parser.add_argument("-os","--ossystem", type=str, default="l", help="OSsystem: Enter the operating system you are using (l = Linux, w = Windows) Default = l")

args = parser.parse_args()

pfile=args.input
pi = args.identity
eval = args.evalue  
ossys = args.ossystem 

if not os.path.exists(f'./structure_predict_out'):
    os.mkdir(f'./structure_predict_out')

if ossys == "l":
    os.system("cp ./toolkit/aln_model.py ./structure_predict_out") 
    os.system("cp ./toolkit/build_model.py ./structure_predict_out")
elif ossys == "w":
    os.system("copy toolkit\\aln_model.py structure_predict_out") 
    os.system("copy toolkit\\build_model.py structure_predict_out")


try:
    with open(pfile, 'r') as fh:
        data = fh.readlines()
except IOError:
    print("Unable to open the file. Try again.")
    exit()
if(">" not in data[0]):
        print("Invalid input\nMissing '>' from your fasta file")
        exit()
accs=[]
for i in data:
    if(">" in i):
        a = i.replace(">","").replace("\n","").replace("|","_").replace(" ","_")
        if len(a)>30:
            accs.append(a[:30])
        else:
            accs.append(a)

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

i = 0

for i in range(len(accs)):

    with open(f"./structure_predict_out/{accs[i]}.fasta","w") as wh:
        wh.write(f">{accs[i]}\n{protein_list[i]}")
    wrt = open(f"./structure_predict_out/{accs[i]}.ali","w")
    wrt.write(f">P1;{accs[i]}\nsequence:{accs[i]}:::::::0.00: 0.00\n{protein_list[i]}*")
    wrt.close()
    if ossys == "l" or ossys == "m":
        subprocess.run(f"diamond blastp --query ./structure_predict_out/{accs[i]}.fasta --db ./datasets/daimond_database.dmnd -o ./structure_predict_out/{accs[i]}_blast.tsv --evalue {eval} --quiet",shell=True)
    elif ossys == "w":
        cmd = rf"datasets\diamond.exe blastp --query structure_predict_out\{accs[i]}.fasta --db datasets\daimond_database.dmnd -o structure_predict_out\{accs[i]}_blast.tsv --evalue {eval} --quiet"
        subprocess.run(cmd, shell=True)


for i in range(len(accs)):
    
    print(os.getcwd())
    print(accs[i])
    try:
        with open(f'./structure_predict_out/{accs[i]}_blast.tsv', 'r') as fh:
            data = fh.readlines()
    except IOError:
        print("Unable to open the tsv file. Try again.")
        continue
    ac = []
    ident = []
    for m in range(len(data)):
        a = data[m].split("\n")
        b = a[0].split("\t")
        ident.append(float(b[2]))
        try:
            c = b[1].split('|')
            ac.append(c[1])
        except(IndexError):
            c = b[1].split('_')
            ac.append(c[1])

    acc = []
    try:
        if max(ident) > pi:
            for cut in range(len(ident)):

                if ident[cut] > pi:
                    acc.append(ac[cut])
        else:
            err = open(f"./structure_predict_out/{accs[i]}_empty.txt","w")
            err.write(f"No templates avialable for {accs[i]} when percent of identity is {pi}")
            err.close()
            continue
    except(ValueError):
        err = open(f"./structure_predict_out/{accs[i]}_empty.txt","w")
        err.write(f"No templates avialable for {accs[i]} when percent of identity is {pi}")
        err.close()
        continue
    # alph=[]
    uni =[]
    templates = []
    for id in range(len(acc)):
        alph = f"https://alphafold.ebi.ac.uk/files/AF-{acc[id]}-F1-model_v4.pdb"
        a=requests.get(alph)
        f2=a.text
        if "<Error>" not in f2:
            f=open(f"./structure_predict_out/{accs[i]}_{acc[id]}_template.pdb","w")
            f.write(f2)
            f.close()
            templates.append(f"{accs[i]}_{acc[id]}_template")
        
    if len(templates) == 0:
        err = open(f"./structure_predict_out/{accs[i]}_empty.txt","w")
        err.write(f"No templates avialable for {accs[i]} when percent of identity is {pi}")
        err.close()
        continue
        
  
    jn = " ".join(templates)
    os.chdir("./structure_predict_out")
    os.system(f"python aln_model.py -ali {accs[i]} -tmp {jn}")
    os.system(f"python build_model.py -ali {accs[i]}-multiple_templates -tmp {jn}")

    os.chdir("../")
