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


for ac in range(len(acc)):
    plen=len(protein_list[ac])
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


    data={f'A [{a}]':a,f'C [{c}]':c,f'D [{d}]':d,f'E [{e}]':e,f'F [{f}]':f,f'G [{g}]':g,f'H [{h}]':h,f'I [{i}]':i,f'K [{k}]':k,f'L [{l}]':l,f'M [{m}]':m,f'N [{n}]':n,f'P [{p}]':p,f'Q [{q}]':q,f'R [{r}]':r,f'S [{s}]':s,f'T [{t}]':t,f'V [{v}]':v,f'W [{w}]':w,f'Y [{y}]':y}
    key=list(data.keys())
    value=list(data.values())
    skolors= [0,2,1,3,5,4,6,8,7,9,11,10,12,14,13,15,17,16,18,20,19]


    explode = [0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005]
    cmap=plt.get_cmap('tab20')
    colors=cmap(skolors)

    plt.pie(value,colors=colors,labels=key,explode=explode,radius = 1.2)
    plt.pie(value,colors='w',radius = 0.7)
    plt.text(0, 0, f"{acc[ac]} [{plen}]", fontsize=12, ha='center')
    if not os.path.exists("aa_composition_output"):
        os.mkdir("aa_composition_output")
    plt.savefig(f"./aa_composition_output/{acc[ac]}.svg",bbox_inches="tight")
    plt.clf()


