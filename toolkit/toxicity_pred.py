import argparse  
import warnings
import os
import numpy as np
import pandas as pd
import joblib
import matplotlib.pyplot as plt
from matplotlib import ticker

warnings.filterwarnings('ignore')
parser = argparse.ArgumentParser(description='Please provide following arguments.')
parser.add_argument("-i", "--input", type=str, required=True, help="Input: FASTA format file of protein sequence/sequences")
parser.add_argument("-t","--threshold", type=float, default=0.38, help="Threshold: The threshold value for toxicity prediction between 0 to 1 by default 0.38")
parser.add_argument("-s","--segments",type=int, default=50, help="segments: Segment length for toxicity prediction by default 50")


args = parser.parse_args()

pfile= args.input        

Threshold= float(args.threshold)

Parts=int(args.segments)

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

if not os.path.exists('./toxicity_output/'):
    os.mkdir('./toxicity_output/')

standard = list("ACDEFGHIKLMNPQRSTVWY")


for ac in range(len(accs)):
    print(accs[ac])
    sequence =''
    seqid=[]
    data1 = protein_list[ac]

    for each in data1:
        sequence=sequence+each.replace('\n','')
    

    start=[]
    end=[]
 
    lngth = []
    for i in range(len(sequence)):
        lngth.append(i+1)
    seq = [sequence[i:i + Parts] for i in range(len(sequence)-(Parts-1))]
    flen = [lngth[i:i + Parts] for i in range(len(lngth)-(Parts-1))]
    for i in flen:
        start.append(i[0])
        end.append(i[-1])
    for i in range (1,len(seq)+1):
        seqid.append("Seg_"+str(i))



    dtfm1 = pd.DataFrame(seq, columns=["Seq"])

    d1 = []
    for j in dtfm1['Seq']:
        c = []
        for i in standard:
            count = 0
            for k in j:
                temp1 = k
                if temp1 == i:
                    count += 1
                comp = (count/len(j))*100
            c.append(comp)
        d1.append(c)
    dtfm2 = pd.DataFrame(d1)
    head = []
    for i in standard:
        head.append('AA_'+i)
    dtfm2.columns = head
    dtfm2.to_csv('./toxicity_output/aa_comp.csv', index=None, header=False)



    q=1
    zz = dtfm1.Seq
    d2 = []
    for i in range(0,len(zz)):
        c = []
        for j in standard:
            for k in standard:
                count = 0
                temp = j+k
                for m in range(0,len(zz[i])-q):
                    b = zz[i][m:m+q+1:q]
                    b.upper()
                    if b == temp:
                        count += 1
                    comp = (count/(len(zz[i])-(q)))*100
                c.append(comp)
        d2.append(c)
    dtfm3 = pd.DataFrame(d2)
    head = []
    for i in standard:
        for j in standard:
            head.append("DP"+str(q)+"_"+i+j)
    dtfm3.columns = head
    dtfm3.to_csv('./toxicity_output/dipep_comp.csv', index=None, header=False)

    dtfm = pd.DataFrame()
    a=[]
    file = './toxicity_output/aa_comp.csv'
    file1 = './datasets/toxinpred3.0_model/toxinpred3.0_model.pkl'
    file2 = './toxicity_output/dipep_comp.csv'
    imp_model = joblib.load(file1)

    test1 = np.loadtxt(file, delimiter=',')
    test2 = np.loadtxt(file2, delimiter=',')
    test3 = np.concatenate([test1,test2], axis=1)
    X_test = test3
    y_score=imp_model.predict_proba(X_test)

    y_s=y_score.tolist()
    dtfm = pd.DataFrame(y_s)
    dtfm_final = dtfm.iloc[:,-1]

    os.remove('./toxicity_output/aa_comp.csv')
    os.remove('./toxicity_output/dipep_comp.csv')
    list = dtfm_final.tolist()
    
    mlscore=[]
    toxic=[]
    acc=[]
    normal=[]
    pred = []
    name=accs[ac]


    for i in range(0, len(seq)):
        mlscore.append(float(list[i]))
        if float(mlscore[i]) >= Threshold:
            normal.append(Threshold)
            toxic.append(float(mlscore[i]))
            pred.append("Toxic")
            acc.append(f"Seg_{i+1}")
        else:
            normal.append(0)
            toxic.append(0)
            pred.append("Non-toxic")
            acc.append(f"Seg_{i+1}")
    with open(f'./toxicity_output/{name}.tsv', 'w') as wh:
        wh.write(f"Segments\tSequences\tPositions\tML_score\tToxicity_ml_score\tPrediction\n")
        for b in range(0, len(seq)):
            wh.write(f"{acc[b]}\t{seq[b]}\t{start[b]}-{end[b]}\t{mlscore[b]}\t{toxic[b]}\t{pred[b]}\n")
   
   
    if len(seq) < 2000:
        if len(seq)<20:
            plt.figure(figsize=(5,5))
        else:
            plt.figure(figsize=(len(seq)/5,5))
        tick_spicing=1
        ax=plt.axes()
        ax.axhline(0,color='black')
        ax.axhline(Threshold,color='blue')

        cl = ['red','purple','green']
        labels = ['Amount of Toxicity','Toxic Protein','Non-Toxic Protein']
        ax.bar(acc,mlscore,color='green')
        ax.bar(acc,toxic,color='red') 
        ax.bar(acc,normal,color='purple')
        patches = [plt.Rectangle((0,0),1,1,fc=color, edgecolor='none') for color in cl]
        plt.legend(patches, labels, loc=(1,1.1))

        ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_spicing))
        plt.xticks(acc,rotation=90)
        plt.title(f"Finding Possible Toxic and Non-toxic Region of Peptides of {name}", fontsize=40)
        plt.xlabel(f"Segments of Amino Acids", fontsize=30)
        plt.ylabel("ML Score", fontsize=30)
        plt.savefig(f"./toxicity_output/{name}.svg",bbox_inches="tight")
    

 
 