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

mw = {
      'A': 71.0788, 'R': 156.1875, 'N': 114.1038, 'D': 115.0886, 'C': 103.1388,
      'E': 129.115, 'Q': 128.1307, 'G': 57.0519, 'H': 137.1411, 'I': 113.1594,
      'L': 113.1594, 'K': 128.1741, 'M': 131.1926, 'F': 147.1766, 'P': 97.1167,
      'S': 87.0782, 'T': 101.1051, 'W': 186.2132, 'Y': 163.1760, 'V': 99.1326
    }

if not os.path.exists('ext_coeff_output'):
        os.mkdir('ext_coeff_output')
  
with open("./ext_coeff_output/Extinction_coefficient.tsv","w") as wh:
    wh.write("ECprotCC : Exinction Co-efficient of protein while considering Cystine residues\nABprotCC : Absorbance of protein while considering Cystine residues at 280nm\nConcProtCC : Concentration of protein while considering Cystine residues where path length L=1cm\n")
    wh.write("ECprotC : Exinction Co-efficient of protein while considering Cysteine residues\nABprotCC : Absorbance of protein while considering Cysteine residues at 280nm\nConcProtC : Concentration of protein while considering Cysteine residues where path length L=1cm\n")
    wh.write("FASTA Identifier\tECprotCC(M^-1 cm^-1)\tABprotCC(mg/ml)\tConcProtCC(mol/ml)\tECprotC(M^-1 cm^-1)\tABprotC(mg/ml)\tConcProtC(mol/ml)\n")
    for ac in range(len(acc)):
        plen=len(protein_list[ac])
        c = protein_list[ac].count("C")
        w = protein_list[ac].count("W")
        y = protein_list[ac].count("Y")
        u = protein_list[ac].count("U")

        add=0
        for b in range(0,plen):
            base=protein_list[ac][b]
            add=add+mw[f"{base}"]
        mwprot=add+18.01524

        ecy=1490
        ecw=5500
        ecc=125

        if(c%2 ==0):
            ecprot=y*ecy+w*ecw+(c/2)*ecc
        elif(c%2 !=0):
            ecprot=y*ecy+w*ecw+((c-1)/2)*ecc
        elif(c<=1):
            ecprot=y*ecy+w*ecw

        if(ecprot==0):
            abprot=0
            ecprot1=0
            abprot1=0
            cprot=0
            cprot1=0
        else:
            abprot=ecprot/mwprot
            cprot=(abprot/ecprot)
            ecprot1=y*ecy+w*ecw
            abprot1=ecprot1/mwprot
            cprot1=(abprot1/ecprot1)
        wh.write(f"{acc[ac]}\t{ecprot}\t{round(abprot,3)}\t{cprot}\t{ecprot1}\t{round(abprot1,3)}\t{cprot1}\n")
            
         
            
            
