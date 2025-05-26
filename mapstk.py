import subprocess
import argparse  
import warnings
import os
warnings.filterwarnings('ignore')
parser = argparse.ArgumentParser(description='Please provide following arguments.')
parser.add_argument("-f", "--input", type=str, required=True, help="Input: FASTA format file of protein sequence/sequences")
parser.add_argument("-op","--ossystem", type=str, default="l", help="ossystem: Enter the operating system you are using (l = Linux, w = Windows and m = MacOS) Default = l")
parser.add_argument("-t","--tools", nargs='+', default= ["all_tools"], help='Tools: Choose the tools you want to run: aa_composition, aliphatic_index, charged_residues, exitinction_coefficient, gravy, half_life, hydropathy_plot, instability_index, isoelectric_point, molecular_weight, shannon_entropy, structure_pred, total_charge, toxicity_pred, transm_pred and by default: all_tools')
parser.add_argument("-pi","--identity", type=float, default=90.0, help="Identity: Minimum percent of identity for finding templates in structure prediction Default = 90.0")
parser.add_argument("-e","--evalue", type=float, default=0.01, help="E-value: E-value for BLAST in structure prediction Default = 0.01")
parser.add_argument("-db", "--database", type=str, default= "uniprot_rev.fasta.gz", help="Database: Provide a protein sequence database in '.fasta' or '.fasta.gz' format for BLASTP executation Default = uniprot_sport.fasta.gz")
parser.add_argument("-ph", "--ph", type=float, default= 7.0, help="pH: Provide pH at which you want to analyse the total charge by default 7.0")
parser.add_argument("-th","--threshold", type=float, default=0.38, help="Threshold: The threshold value for toxicity prediction between 0 to 1 by default 0.38")
parser.add_argument("-s","--segments",type=int, default=50, help="Parts: Segment length for toxicity prediction Default = 50")
args = parser.parse_args()

pfile=args.input
tools = args.tools
ident = args.identity
evalue = args.evalue
db = args.database
ph = args.ph
threshold = args.threshold
segment = args.segments
ossys = args.ossystem
if ossys == "w":
    os.system(f"copy {pfile} toolkit")
elif ossys == "l":
    os.system(f"cp ./{pfile} ./toolkit")
elif ossys == "m":
    os.system(f"cp ./{pfile} ./toolkit")
else:
    print("Your operating system is not Windows or Linux or MacOS")
    exit()

if "aa_composition" in tools or "all_tools" in tools:
    print("***aa_composition tool is running***")
    subprocess.run(f"python ./toolkit/aa_composition.py -i {pfile}",shell=True)

if "aliphatic_index" in tools or "all_tools" in tools:
    print("***aliphatic_index tool is running***")
    subprocess.run(f"python ./toolkit/aliphatic_index.py -i {pfile}",shell=True)

if "charged_residues" in tools or "all_tools" in tools:
    print("***charged_residues tool is running***")
    subprocess.run(f"python ./toolkit/charged_residues.py -i {pfile}",shell=True)

if "exitinction_coefficient" in tools or "all_tools" in tools:
    print("***exitinction_coefficient tool is running***")
    subprocess.run(f"python ./toolkit/exitinction_coefficient.py -i {pfile}",shell=True)

if "gravy" in tools or "all_tools" in tools:
    print("***gravy tool is running***")
    subprocess.run(f"python ./toolkit/gravy.py -i {pfile}",shell=True)

if "half_life" in tools or "all_tools" in tools:
    print("***half_life tool is running***")
    subprocess.run(f"python ./toolkit/half_life.py -i {pfile}",shell=True)

if "hydropathy_plot" in tools or "all_tools" in tools:
    print("***hydropathy_plot tool is running***")
    subprocess.run(f"python ./toolkit/hydropathy_plot.py -i {pfile}",shell=True)

if "instability_index" in tools or "all_tools" in tools:
    print("***instability_index tool is running***")
    subprocess.run(f"python ./toolkit/instability_index.py -i {pfile}",shell=True)

if "isoelectric_point" in tools or "all_tools" in tools:
    print("***isoelectric_point tool is running***")
    subprocess.run(f"python ./toolkit/isoelectric_point.py -i {pfile}",shell=True)

if "molecular_weight" in tools or "all_tools" in tools:
    print("***molecular_weight tool is running***")
    subprocess.run(f"python ./toolkit/molecular_weight.py -i {pfile}",shell=True)

if "shannon_entropy" in tools or "all_tools" in tools:
    print("***shannon_entropy tool is running***")
    subprocess.run(f"python ./toolkit/shannon_entropy.py -i {pfile}",shell=True)

if "structure_pred" in tools or "all_tools" in tools:
    print("***structure_pred tool is running***")
    if ossys == "l" or ossys == "m":
        os.system(f"diamond makedb --in datasets/{db} -d datasets/daimond_database")
    elif ossys == "w":
        os.system(f"datasets/diamond.exe makedb --in datasets/{db} -d datasets/daimond_database")
    subprocess.run(f"python ./toolkit/structure_pred.py -i {pfile} -pi {ident} -e {evalue} -os {ossys}",shell=True)

if "total_charge" in tools or "all_tools" in tools:
    print("***total_charge tool is running***")
    subprocess.run(f"python ./toolkit/total_charge.py -i {pfile} -ph {ph}",shell=True)

if "toxicity_pred" in tools or "all_tools" in tools:
    print("***toxicity_pred tool is running***")
    subprocess.run(f"python ./toolkit/toxicity_pred.py -i {pfile} -t {threshold} -s {segment}",shell=True)

if "transm_pred" in tools or "all_tools" in tools:
    print("***transm_pred tool is running***")
    subprocess.run(f"python ./toolkit/transm_pred.py -i {pfile}",shell=True)
os.remove(f"toolkit/{pfile}")