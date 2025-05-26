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
    
halflifeM = {
        "A": "4.4 hour", "R": "1 hour", "N": "1.4 hour", "D": "1.1 hour", "C": "1.1 hour",
        "Q": "0.8 hour", "E": "1 hour", "G": "30 hour", "H": "3.5 hour", "I": "20 hour",
        "L": "5.5 hour", "K": "1.3 hour", "M": "30 hour", "F": "1.1 hour", "P": ">20 hour",
        "S": "1.9 hour", "T": "7.2 hour", "W": "2.8 hour", "Y": "2.8 hour", "V": "100 hour"
      }
halflifeY = {
        "A": ">20 hour", "R": "2 min", "N": "3 min", "D": "3 min", "C": ">20 hour",
        "Q": "10 min", "E": "30 min", "G": ">20 hour", "H": "10 min", "I": "30 min",
        "L": "3 min", "K": "3 min", "M": ">20 hour", "F": "3 min", "P": ">20 hour",
        "S": ">20 hour", "T": ">20 hour", "W": "3 min", "Y": "10 min", "V": ">20 hour"
      }
      
halflifeE = {
        "A": ">10 hour", "R": "2 min", "N": ">10 hour", "D": ">10 hour", "C": ">10 hour",
        "Q": ">10 hour", "E": ">10 hour", "G": ">10 hour", "H": ">10 hour", "I": ">10 hour",
        "L": "2 min", "K": "2 min", "M": ">10 hour", "F": "2 min", "P": "...", "S": ">10 hour",
        "T": ">10 hour", "W": "2 min", "Y": "2 min", "V": ">10 hour"
      }
if not os.path.exists('Half_Life_output'):
        os.mkdir('Half_Life_output')


with open("./Half_Life_output/half_life.tsv","w") as wh:
    wh.write("FASTA Idnetifier\tHalf_life_Mammal\tHalf_life_Yeast\tHalf_life_E.coli\n")
    for i in range(len(acc)):
        wh.write(f"{acc[i]}\t{halflifeM[protein_list[i][0]]}\t{halflifeY[protein_list[i][0]]}\t{halflifeE[protein_list[i][0]]}\n")


