MAPSTK: Multi-parameter Analysis Protein Sequence Toolkit
A Toolkit for Multiple Protein Sequence Analysis

Introduction

MAPSTK have 14 different tools (Amino Acid Composition, Total Charged Residues, In-vivo Half-life, Molecular Weight, Aliphatic Index, Instability Index, Iso-electric Point, Extinction Coefficient, Hydropathy Plot, GRAVY, Entropy, Transmembrane Prediction, Toxicity and Structure Prediction.) which can extract the values of more than 20 properties of proteins from its sequence. It is obvious that each of the 14 different tools are made for a specific purpose. These tools will generate data both in tabular TSV and graphical format as well. Matplotlib module was used for generating the graph. 

Standalone

The standalone version of MAPSTK is developed in Python3 and requires the following libraries for successful execution:

matplotlib: pip install matplotlib
requests: pip install requests
urllib3: pip install urllib3==1.26.6
pandas: pip install pandas
joblib: pip install joblib
numpy: pip install numpy
scikit-learn: pip install scikit-learn==1.0.2
ramachandraw: pip install ramachandraw==0.2.3

Important Note

Due to the large size of the model file in Toxicity tool, it has been compressed. Before using the code or model, it is essential to unzip the file. The compressed file must be extracted to its original form to ensure the proper functioning of the code.

For Linux users, in the Structure Prediction Tool of MAPSTK, extract the "diamond.tar.xz" file to obtain the executable "diamond" tool. Additionally, remove ".exe" from "./diamond.exe" in the script at line 79 to ensure compatibility. 

USAGE of ToolKit:

Amino Acid Composition:
usage: python aa_composition.py [-h] -i INPUT
Please provide following arguments.
options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
Input: FASTA format file of protein sequence/sequences

Charged Residues:
usage: charged_residues.py [-h] -i INPUT
Please provide following arguments.
options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
Input: FASTA format file of protein sequence/sequences

Aliphatic Index:
usage: aliphatic_index.py [-h] -i INPUT
Please provide following arguments.
options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
Input: FASTA format file of protein sequence/sequences

Hydropathy Plot:
usage: hydropathy_plot.py [-h] -i INPUT
Please provide following arguments.
options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
Input: FASTA format file of protein sequence/sequences

GRAVY:
usage: gravy.py [-h] -i INPUT
Please provide following arguments.
options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
Input: FASTA format file of protein sequence/sequences

Transmembrane Prediction:
usage: transm_pred.py [-h] -i INPUT
Please provide following arguments.
options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
Input: FASTA format file of protein sequence/sequences

Instability Index:
usage: instability_index.py [-h] -i INPUT
Please provide following arguments.
options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
Input: FASTA format file of protein sequence/sequences

Isoelectric Point:
usage: isoelectric_point.py [-h] -i INPUT
Please provide following arguments.
options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
Input: FASTA format file of protein sequence/sequences

Molecular Weight:
usage: molecular_weight.py [-h] -i INPUT
Please provide following arguments.
options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
Input: FASTA format file of protein sequence/sequences

Shannon Entropy:
usage: shannon_entropy.py [-h] -i INPUT
Please provide following arguments.
options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
Input: FASTA format file of protein sequence/sequences

In-vivo Half-life:
usage: half_life.py [-h] -i INPUT
Please provide following arguments.
options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
Input: FASTA format file of protein sequence/sequences

Extinction Co-efficient:
usage: exitinction_coefficient.py [-h] -i INPUT
Please provide following arguments.
options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
Input: FASTA format file of protein sequence/sequences

Toxicity:
usage: toxicity.py [-h] -i INPUT [-t THRESHOLD] [-p PARTS]
Please provide following arguments.
options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input: FASTA format file of protein sequence/sequences
  -t THRESHOLD, --threshold THRESHOLD
                        Threshold: Value between 0 to 1 by default 0.38
  -p PARTS, --parts PARTS
                        Parts: Segment length Default = 50

Structure Pred:
usage: structure_pred.py [-h] -i INPUT [-p IDENTITY] [-e EVALUE]
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input: FASTA format file of protein sequence/sequences
  -p IDENTITY, --identity IDENTITY
                        Identity: Minimum percent of identity Default = 90.0
-e EVALUE, --evalue EVALUE
                        E-value: E-value for BLAST Default = 0.01

Input File: Users can upload input files in the FASTA format. Both multi-sequence and single-sequence FASTA files can be submitted by the user. 

Output Folder: Every tool in the MAPSTK automatically generates an output folder and stores the result files within that folder.
