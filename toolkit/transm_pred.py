import requests
import time
import re
import sys
import os
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib import style
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

if not os.path.exists('./transm_output/'):
    os.mkdir('./transm_output/')


for ac in range(len(accs)):
	sequence=protein_list[ac]

	alph = f'http://cctop.ttk.hu/direct/submit?id={accs[ac]}&seq={sequence}'
	a=requests.get(alph)
	f2=a.text
	print("Accession Id of input sequence is:",accs[ac])
	print("job id:",f2)

	while 1:
		alph2=f'http://cctop.ttk.hu/direct/poll?hash={f2}'
		d=requests.get(alph2)
		f3=d.text
		print("status:",f3)
		alph3 = f'https://cctop.ttk.hu/job/results/{f2}/xml'
		e=requests.get(alph3)
		f4=e.text
		if('<CCTOPItems>\n' in f4):
			outf='cctop.txt'
			with open(outf,'w') as wh:
				wh.write(f4)
			break
		# time.sleep(10)
		if 'Invalid' not in f3:
			time.sleep(10)
		else:
			break


	seq_len = len(sequence)	


	try:
		with open(outf, 'r') as fh:
			data = fh.readlines()
	except NameError:
		quit()

	os.remove(outf)
	try:
		for i in range(len(data)):
			if 'name="TMHMM">\n' in data[i]:
					x=i
			if '</CCTOPItem>\n'  in data[i]:
					y=i


		m=[]
		o=[]
		ii=[]
		for j in range(x+1,y-2):
			if 'loc="M"/>\n' in data[j]:
				m.append(data[j])
			elif 'loc="O"/>\n' in data[j]:
				o.append(data[j])
			elif 'loc="I"/>\n' in data[j]:
				ii.append(data[j])
	except(NameError):
		ww = open(f"./transm_output/{accs[ac]}_non-membrane.txt","w")
		ww.write("This sequence have no transmembrane segemnt")
		ww.close()
		continue
	
	if len(m) == 0:
		continue
		
	inn=[]
	out=[]
	mm=[]

	for k in range(len(m)):
		int_values = [int(match.group(1)) for match in re.finditer(r'(\d+)', m[k])]
		mm.append(int_values)
	for k in range(len(o)):
		int_values = [int(match.group(1)) for match in re.finditer(r'(\d+)', o[k])]
		out.append(int_values)
	for k in range(len(ii)):
		int_values = [int(match.group(1)) for match in re.finditer(r'(\d+)', ii[k])]
		inn.append(int_values)


	result1 = [element for sublist in mm for element in range(sublist[0]-1, sublist[1])]
	result2 = [element for sublist in out for element in range(sublist[0]-1, sublist[1])]
	result3 = [element for sublist in inn for element in range(sublist[0]-1, sublist[1])]

	xax=[]
	yax=[]
	xm=[]
	xi=[]
	xo=[]
	for a in range(0,seq_len):
		flag=0
		xax.append(f"{sequence[a]}{a+1}")
		
		for am1 in result1:
			if(a==am1):
				flag=1
		for am2 in result2:
			if(a==am2):
				flag=2
		if(flag==1):
			xm.append(int(a))
			xo.append(0)
			xi.append(0)
		elif(flag==2):
			xo.append(int(a))
			xm.append(0)
			xi.append(0)
		else:
			xi.append(int(a))
			xm.append(0)
			xo.append(0)
		
	name = accs[ac]
	with open(f'./transm_output/{name}.tsv', 'w') as wh:
		wh.write(f"Type\tPositions\n")
		for b in range(0,len(mm)):
			wh.write(f"Transmebrane_{b+1}\t{mm[b]}\n")
		for b in range(0,len(out)):
			wh.write(f"Out_side_{b+1}\t{out[b]}\n")
		for b in range(0,len(inn)):
			wh.write(f"In_side_{b+1}\t{inn[b]}\n")
	# plt.style.use('dark_background')
	plt.figure(figsize=(seq_len/5,8))
	ax=plt.axes()
	tick_spicing=1
	ax.axhline(0,color='black')

	cl = ['orange','teal','skyblue']
	labels = ['Transmembrane_Segments','Inside','Outside']

	
	ax.bar(xax,xo,color='teal')
	ax.bar(xax,xi,color='skyblue')
	ax.bar(xax,xm,color='orange')

	patches = [plt.Rectangle((0,0),1,1,fc=color, edgecolor='none') for color in cl]
	plt.legend(patches, labels, fontsize='30', loc='best')

	ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_spicing))
	plt.xticks(xax,rotation=90)
	plt.margins(0)
	plt.title(f"Finding Possible Transmembrane Region of Peptides of {name}", fontsize=40)
	plt.xlabel("Positions of Amino Acids", fontsize=30)
	plt.ylabel("Positions of Regions", fontsize=30)
	plt.savefig(f"./transm_output/{name}.svg",bbox_inches="tight")




