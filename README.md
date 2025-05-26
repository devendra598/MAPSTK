# MAPSTK: Multi-parameter Analysis Protein Sequence Toolkit
## A Toolkit for Multiple Protein Sequence Analysis

### Introduction

MAPSTK provides 15 different tools for protein sequence analysis, which includes:

- Amino Acid Composition
- Total Charged Residues
- In-vivo Half-life
- Molecular Weight
- Aliphatic Index
- Instability Index
- Iso-electric Point
- Total Charge
- Extinction Coefficient
- Hydropathy Plot
- GRAVY
- Entropy
- Transmembrane Prediction
- Toxicity
- Structure Prediction

Each tool serves a specific purpose, generating results in both tabular TSV format and graphical format.

---

## Standalone Installation

The standalone version of MAPSTK is developed in Python3 

### Download Python3

You can download and install Python3 from the official Python website:  
🔗 [Python Download](https://www.python.org/downloads)

### Download Conda

You can download and install Conda from official website:
🔗 [Conda Download](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)

### Creating a Conda Environment with Python 3.7.0

Some tools may encounter errors with recent versions of Python. To ensure compatibility, it is recommended to use **Python 3.7.0**.
You can create an environment with Python 3.7.0 using the following command:
```
conda create -n mapstk_env python==3.7.0
```
Activate the Conda environment using the following command:
```
conda activate mapstk_env
```

### Download and Install Diamond with Conda

For Linux and MacOS environment:
```
conda install bioconda::diamond
```

For Windows Command Line:
Download `diamond.exe` from the official website:
🔗 [Diamond Download](https://github.com/bbuchfink/diamond/releases/download/v2.1.11/diamond-windows.zip)
Extract `diamond.exe` from the zip-folder and place it in the "datasets" folder for proper blastp execution. 

### Download and Install MODELLER with Conda

To add salilab channels to conda:
```
conda config --add channels salilab
```
Install modeller through conda:

In the command provided below, replace "XXXX" with MODELLER license key.
You can obtain MODELLER license key from this website:
🔗 [Get MODELLER license key](https://salilab.org/modeller/registration.html)

For Linux and MacOS environment:
```
KEY_MODELLER=XXXX conda install modeller
```

For Windows Command Line:
Set MODELLER license key
```
set KEY_MODELLER=XXXX
```
Install Modeller:
```
conda install modeller
```

## Install dependensies 

MAPSTK depends on the following packages:
```
pip install matplotlib requests==2.31.0 pandas joblib numpy scikit-learn==1.0.2 ramachandraw==0.2.3 pyfamsa seaborn urllib3==1.26.6
```

## Important Note

- The **Toxicity tool** contains a large model file that is present in compressed format in **datasets** folder. Before using it, **unzip** the file to ensure proper execution.
- The Structure Prediction and Transmembrane tools rely on **API**, requiring users to stay connected to the **internet** while running them.
- The **Structure Prediction** tool requires a protein sequence database to operate. You can use the default uniprot reviewed database **uniprot_sport.fasta.gz**, the download link for this database is here 🔗 [Uniprot_sport Database](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz) or you can also provide any other protein sequence database from Uniprot in **.fasta** or **.fasta.gz** format. Ensure that the database file is placed in the **datasets** folder for proper execution. 
You may access your desired UniProt database from this site:
🔗 [Uniprot Database](https://www.uniprot.org/help/downloads)

---

## Usage

Each tool in MAPSTK requires a **FASTA format file** as input. Below are the usage instructions for each tool:

### **To see the help box**
```
python mapstk.py -h
```

### **To run all the tools at once**
```
python mapstk.py -f [FASTA file name] -op [Operating system name] -ph [pH Value] -th [THRESHOLD] -s [SEGMENTS] -pi [IDENTITY] -e [EVALUE] -db [DATABASE]
```
- `-f INPUT` Input: FASTA format file of protein sequence/sequences.
- `-op ossystem` (default = l) Enter the operating system you are using (l = Linux and w = Windows)
- `-ph PH` (default: 7.0) specifies the pH value (0 to 14) applicable only for Total Charge Tool.
- `-th THRESHOLD` (default: 0.38) specifies the threshold value (0 to 1) applicable only for Toxicity Tool.
- `-s SEGMENTS` (default: 50) defines segment length applicable only for Toxicity Tool.
- `-pi IDENTITY` (default: 90.0) sets the minimum percent identity applicable only for Structure Prediction Tool.
- `-e EVALUE` (default: 0.01) defines the e-value for BLAST applicable only for Structure Prediction Tool.
- `-db DATABASE` (default: uniprot_rev.fasta.gz) Provide a protein sequence database in '.fasta' or '.fasta.gz' format for BLASTP executation

### **To run multiple tools at once**
```
python mapstk.py -f [FASTA file name] -op [Operating system name] -t [1st Tool Name] [2nd Tool Name] ... [nth Tool Name]
```

### **Amino Acid Composition**
```
python mapstk.py -f [FASTA file name] -op [Operating system name] -t aa_composition
```

### **Charged Residues**
```
python mapstk.py -f [FASTA file name] -op [Operating system name] -t charged_residues
```

### **Aliphatic Index**
```
python mapstk.py -f [FASTA file name] -op [Operating system name] -t aliphatic_index
```

### **Hydropathy Plot**
```
python mapstk.py -f [FASTA file name] -op [Operating system name] -t hydropathy_plot
```

### **GRAVY**
```
python mapstk.py -f [FASTA file name] -op [Operating system name] -t gravy
```

### **Transmembrane Prediction**
```
python mapstk.py -f [FASTA file name] -op [Operating system name] -t transm_pred
```

### **Instability Index**
```
python mapstk.py -f [FASTA file name] -op [Operating system name] -t instability_index
```

### **Isoelectric Point**
```
python mapstk.py -f [FASTA file name] -op [Operating system name] -t isoelectric_point
```

### **Total Charge**
```
python mapstk.py -f [FASTA file name] -op [Operating system name] -t total_charge
```

### **Molecular Weight**
```
python mapstk.py -f [FASTA file name] -op [Operating system name] -t molecular_weight
```

### **Shannon Entropy**
```
python mapstk.py -f [FASTA file name] -op [Operating system name] -t shannon_entropy -ph [pH Value]
```
- `-ph PH` (default: 7.0) specifies the pH value (0 to 14).

### **In-vivo Half-life**
```
python mapstk.py -f [FASTA file name] -op [Operating system name] -t half_life
```

### **Extinction Coefficient**
```
python mapstk.py -f [FASTA file name] -op [Operating system name] -t exitinction_coefficient
```

### **Toxicity**
```
python mapstk.py -f [FASTA file name] -op [Operating system name] -t toxicity_pred -th [THRESHOLD] -s [SEGMENTS]
```
- `-th THRESHOLD` (default: 0.38) specifies the threshold value (0 to 1).
- `-s SEGMENTS` (default: 50) defines segment length.

### **Structure Prediction**
```
python mapstk.py -f [FASTA file name] -op [Operating system name] -t structure_pred -pi [IDENTITY] -e [EVALUE]
```
- `-pi IDENTITY` (default: 90.0) sets the minimum percent identity.
- `-e EVALUE` (default: 0.01) defines the e-value for BLAST.
- `-db DATABASE` (default: uniprot_rev.fasta.gz) Provide a protein sequence database in '.fasta' or '.fasta.gz' format for BLASTP executation

---

## Input & Output

- **Input File:** Users can upload input files in the FASTA format. Both multi-sequence and single-sequence FASTA files can be submitted by the user. 
- **Output Folder:** Every tool in the MAPSTK automatically generates an output folder and stores the result files within that folder.

---

For any further queries, please contact:
devbioinfo@gmail.com,
rajrohanrnr@gmail.com
