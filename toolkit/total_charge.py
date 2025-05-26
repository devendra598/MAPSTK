import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib import style
import numpy as np
import seaborn as sns
import pandas as pd
import os
import argparse  
import warnings
warnings.filterwarnings('ignore')
parser = argparse.ArgumentParser(description='Please provide following arguments.') 

parser.add_argument("-i", "--input", type=str, required=True, help="Input: FASTA format file of protein sequence/sequences")
parser.add_argument("-ph", "--ph", type=float, default= 7.0, required=True, help="pH: Provide pH at which you want to analyse the total charge")
args = parser.parse_args()

pfile=args.input
ph = args.ph
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

if not os.path.exists('Total_charge'):
        os.mkdir('Total_charge')

def calculate_charge(protein_seq, ph):
    # pKa values for amino acids and termini
    pKa_values = {
        'D': 3.65,  # Aspartic acid
        'E': 4.25,  # Glutamic acid
        'C': 8.18,  # Cysteine
        'Y': 10.07, # Tyrosine
        'H': 6.0,   # Histidine
        'K': 10.53, # Lysine
        'R': 12.48, # Arginine
        'COOH': 2.34, # C-terminal
        'NH': 9.69   # N-terminal
    }

    # Count amino acids
    aa_counts = {aa: protein_seq.count(aa) for aa in pKa_values if aa not in ['COOH', 'NH']}

    # Calculate charge contributions
    charge = 0
    for aa, count in aa_counts.items():
        if aa in ['D', 'E', 'C', 'Y']:
            charge += -count / (1 + 10**(pKa_values[aa] - ph))
        elif aa in ['H', 'K', 'R']:
            charge += count / (1 + 10**(ph - pKa_values[aa]))

    # Add termini charges
    charge += -1 / (1 + 10**(pKa_values['COOH'] - ph))  # C-terminal
    charge += 1 / (1 + 10**(ph - pKa_values['NH']))     # N-terminal

    return charge
CH = []
if(len(acc)<=2000):
    with open("./Total_charge/total_charge.tsv","w") as wh:
        wh.write("FASTA Idnetifier\tTotal_charge\n")
        for seq in range(len(protein_list)): 
            charge = calculate_charge(protein_list[seq], float(ph))
            CH.append(charge)
            wh.write(f"{acc[seq]}\t{charge}\n")

combined = list(zip(acc, CH))
sorted_combined = sorted(combined, key=lambda x: x[1])

# Unzip the sorted pairs back into separate lists
sorted_labels, sorted_values = zip(*sorted_combined)

# Convert tuples back to lists (if needed)
sorted_labels = list(sorted_labels)
sorted_values = list(sorted_values)

# print("Sorted Labels:", sorted_labels)
# print("Sorted Values:", sorted_values)

data = pd.DataFrame({'FASTA Idnetifier': sorted_labels, 'Total Charge': sorted_values})
data = data.set_index('FASTA Idnetifier')

# Plot the heatmap
if len(acc)<20:
    plt.figure(figsize=(23,22))
else:
    plt.figure(figsize=(len(acc)/2,60))
sns.heatmap(
    data.T, 
    annot=False,  # Disable annotations
    fmt='.3f', 
    cmap='RdYlGn', 
    cbar=True, 
    linewidths=1, 
    linecolor='black'
)
# Increase x and y labels fontsize
plt.xlabel('FASTA Idnetifier', fontsize=18)


plt.savefig("./Total_charge/total_charge.svg")

#plt.show()
