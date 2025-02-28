import os
import argparse  
import warnings
import requests
from RamachanDraw import phi_psi, plot
import matplotlib.pyplot as plt
import subprocess
import time
import gzip
import shutil
warnings.filterwarnings('ignore')
parser = argparse.ArgumentParser(description='Please provide following arguments.') 


parser.add_argument("-i", "--input", type=str, required=True, help="Input: FASTA format file of protein sequence/sequences")
parser.add_argument("-p","--identity", type=float, default=90.0, help="Identity: Minimum percent of identity Default = 90.0")
parser.add_argument("-e","--evalue", type=float, default=0.01, help="E-value: E-value for BLAST Default = 0.01")
# parser.add_argument("-o","--output", type=str, default='./structure_predict_out', help="Output: Output directory name Default = ./structure_predict_out")
args = parser.parse_args()

pfile=args.input

pi = args.identity
eval = args.evalue  
# outf = args.output  
outf = "./structure_predict_out"

whh1 = open("total_template.csv","a")
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
        try:
            a = i.split("|")
            accs.append(a[1])
        except(IndexError):
            a = i.split("\n")
            accs.append(a[0][1:])

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
if not os.path.exists(f'./structure_predict_out'):
            os.mkdir(f'./structure_predict_out')
i = 0
for i in range(len(accs)):

    with open(f"./structure_predict_out/{accs[i]}.fasta","w") as wh:
        
        wh.write(f">{accs[i]}\n{protein_list[i]}")
    subprocess.run(f"./diamond.exe blastp --query ./structure_predict_out/{accs[i]}.fasta --db ./uniprot_rev.dmnd -o ./structure_predict_out/{accs[i]}_blast.tsv --evalue {eval} --quiet",shell=True)

def ramch(accession):

    ax = plot(f"./structure_predict_out/{accession}_predicted.pdb")
   
    
    plt.savefig(f"./structure_predict_out/{accession}_predicted_ramachandran.png")

    pdb_file = f"./structure_predict_out/{accession}_predicted.pdb"

    phi_psi_dict = phi_psi(pdb_file)
    aa = phi_psi_dict.keys()
    aa=list(aa)
    ph_ps = phi_psi_dict.values()
    ph_ps=list(ph_ps)

    with open(f"./structure_predict_out/{accession}_predicted_phi_psi.tsv","w") as wh:
        wh.write("Amino_acids\tphi_values\tpsi_values\n")
        for l in range(len(aa)):
            wh.write(f"{aa[l]}\t{ph_ps[l][0]}\t{ph_ps[l][1]}\n")

for i in range(len(accs)):
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
        c = b[1].split('|')
        ac.append(c[1])

    acc = []
    if max(ident) > pi:
        for cut in range(len(ident)):

            if ident[cut] == max(ident):
                acc.append(ac[cut])
    else:
        err = open(f"{accs[i]}_empty.txt","w")
        err.write(f"No templates avialable for {accs[i]} when percent of identity is {pi}")
        continue
    alph=[]
    uni =[]
    
    alph.append(f"https://alphafold.ebi.ac.uk/files/AF-{acc[0]}-F1-model_v4.pdb")
    a=requests.get(alph[0])
    f2=a.text
    if "<Error>" not in f2:
        f=open(f"./structure_predict_out/{accs[i]}_{acc[0]}_template.pdb","w")
        f.write(f2)
        f.close()
    if max(ident) == 100:
        try:
            src = f"./structure_predict_out/{accs[i]}_{acc[0]}_template.pdb"
            dst = f"./structure_predict_out/{accs[i]}_predicted.pdb"
            shutil.copyfile(src = src , dst = dst)
            ramch(accs[i])
            continue
        except(FileNotFoundError):
            continue

            
        
    token = "d44dd29a1fd66dc103c3eb3c414cb2e9dfa879fb"

    with open(f"./structure_predict_out/{accs[i]}.fasta","r") as fa:
        fas = fa.readlines()

    with open(f"./structure_predict_out/{accs[i]}_{acc[0]}_template.pdb", "r") as f:
        template_coordinates = f.read()
 
    response = requests.post(
        "https://swissmodel.expasy.org/user_template",
        headers={ "Authorization": f"Token {token}" },
        json={
            "target_sequences": fas[1],
            "template_coordinates": template_coordinates,
            "project_title":f"{accs[i]}_{acc[0]}"
            })
    project_id = response.json()["project_id"]


    while True:
        
        time.sleep(10)

         
        response = requests.get(
            f"https://swissmodel.expasy.org/project/{ project_id }/models/summary/", 
            headers={ "Authorization": f"Token {token}" })

        
        status = response.json()["status"]

        print('Job status is now', status)

        if status in ["COMPLETED", "FAILED"]:
            break
    response_object = response.json()

    if response_object['status'] == 'COMPLETED':
        for model in response_object['models']:
            pdb_file = model['coordinates_url'] 

            
            b = requests.get(pdb_file)
            if b.status_code == 200:
                
                with open(f"./structure_predict_out/{accs[i]}_predicted.pdb.gz", "wb") as f1:
                    f1.write(b.content)

                
                with gzip.open(f"./structure_predict_out/{accs[i]}_predicted.pdb.gz", 'rb') as gz_file:
                    with open(f"./structure_predict_out/{accs[i]}_predicted.pdb", 'wb') as pdb_file:
                        shutil.copyfileobj(gz_file, pdb_file)
            else:
                print(f"Failed to download the PDB file. Status code: {b.status_code}")
        os.remove(f"./structure_predict_out/{accs[i]}_predicted.pdb.gz")
        ramch(accs[i])
    if response_object['status'] == 'FAILED':
        err = open(f"{accs[i]}_empty.txt","w")
        err.write(f"No templates avialable for {accs[i]}")
    
