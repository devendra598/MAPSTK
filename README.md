# MAPSTK: Multi-parameter Analysis Protein Sequence Toolkit
## A Toolkit for Multiple Protein Sequence Analysis

### Introduction

MAPSTK provides 14 different tools for protein sequence analysis, which includes:

- Amino Acid Composition
- Total Charged Residues
- In-vivo Half-life
- Molecular Weight
- Aliphatic Index
- Instability Index
- Iso-electric Point
- Extinction Coefficient
- Hydropathy Plot
- GRAVY
- Entropy
- Transmembrane Prediction
- Toxicity
- Structure Prediction

These tools can extract over 20 properties from protein sequences. Each tool serves a specific purpose, generating results in both tabular TSV format and graphical format using the `matplotlib` module.

---

## Standalone Installation

The standalone version of MAPSTK is developed in Python3 

## Python 3.7.0 Compatibility Notice

Some tools may encounter errors with recent versions of Python. To ensure compatibility, it is recommended to use **Python 3.7.0**.

## Download Python 3.7.0

You can download Python 3.7.0 from the official Python website:  
🔗 [Python 3.7.0 Download](https://www.python.org/downloads/release/python-370/)

**OR**

## Creating a Conda Environment with Python 3.7.0

If you are using Conda, you can create an environment with Python 3.7.0 using the following command:

```
conda create -n myenv python==3.7.0
```
Activate the Conda environment using the following command:
```
conda activate myenv
```

## Install dependensies 

MAPSTK depends on the following packages:
```
pip install matplotlib requests==2.31.0 pandas joblib numpy scikit-learn==1.0.2 ramachandraw==0.2.3 urllib3==1.26.6
```

### Important Note

- The **Toxicity tool** contains a large model file that is compressed. Before using it, **unzip** the file to ensure proper execution.
- The **Structure Prediction** tool requires the **uniprot_rev.dmnd** database to operate. You can find the download link for this database in the **uniprot_rev.txt** file. After downloading, ensure that the **uniprot_rev.dmnd** file is placed in the same directory as the **structure_pred.py** script for proper execution.
- For **Linux users**, in the Structure Prediction Tool, extract the `diamond.tar.xz` file to obtain the `diamond` executable. Additionally, remove `.exe` from `diamond.exe` in the script at **line 79** for compatibility.
- The Structure Prediction and Transmembrane tools rely on **API**, requiring users to stay connected to the **internet** while running them.


---

## Usage

Each tool in MAPSTK requires a **FASTA format file** as input. Below are the usage instructions for each tool:

### **Amino Acid Composition**
```
python aa_composition.py -i gpr_final.fasta
```

### **Charged Residues**
```
python charged_residues.py -i gpr_final.fasta
```

### **Aliphatic Index**
```
python aliphatic_index.py -i gpr_final.fasta
```

### **Hydropathy Plot**
```
python hydropathy_plot.py -i gpr_final.fasta
```

### **GRAVY**
```
python gravy.py -i gpr_final.fasta
```

### **Transmembrane Prediction**
```
python transm_pred.py -i gpr_final.fasta
```

### **Instability Index**
```
python instability_index.py -i gpr_final.fasta
```

### **Isoelectric Point**
```
python isoelectric_point.py -i gpr_final.fasta
```

### **Molecular Weight**
```
python molecular_weight.py -i gpr_final.fasta
```

### **Shannon Entropy**
```
python shannon_entropy.py -i gpr_final.fasta
```

### **In-vivo Half-life**
```
python half_life.py -i gpr_final.fasta
```

### **Extinction Coefficient**
```
python extinction_coefficient.py -i gpr_final.fasta
```

### **Toxicity**
```
python toxicity.py -i gpr_final.fasta [-t THRESHOLD] [-p PARTS]
```
- `-t THRESHOLD` (default: 0.38) specifies the threshold value (0 to 1).
- `-p PARTS` (default: 50) defines segment length.

### **Structure Prediction**
```
python structure_pred.py -i gpr_final.fasta [-p IDENTITY] [-e EVALUE]
```
- `-p IDENTITY` (default: 90.0) sets the minimum percent identity.
- `-e EVALUE` (default: 0.01) defines the e-value for BLAST.

---

## Input & Output

- **Input File:** Users can upload input files in the FASTA format. Both multi-sequence and single-sequence FASTA files can be submitted by the user. 
- **Output Folder:** Every tool in the MAPSTK automatically generates an output folder and stores the result files within that folder.

---

For any further queries, please contact:
devbioinfo@gmail.com,
rajrohanrnr@gmail.com
