# MAPSTK: Multi-parameter Analysis Protein Sequence Toolkit
## A Toolkit for Multiple Protein Sequence Analysis

### Introduction

MAPSTK provides 14 different tools for protein sequence analysis, including:

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

The standalone version of MAPSTK is developed in Python 3 and requires the following dependencies:

```
pip install matplotlib requests urllib3==1.26.6 pandas joblib numpy scikit-learn==1.0.2 ramachandraw==0.2.3
```

### Important Note

- The **Toxicity tool** contains a large model file that is compressed. Before using it, **unzip** the file to ensure proper execution.
- For **Linux users**, in the Structure Prediction Tool, extract the `diamond.tar.xz` file to obtain the `diamond` executable. Additionally, remove `.exe` from `./diamond.exe` in the script at **line 79** for compatibility.

---

## Usage

Each tool in MAPSTK requires a **FASTA format file** as input. Below are the usage instructions for each tool:

### **Amino Acid Composition**
```
python aa_composition.py -i input.fasta
```

### **Charged Residues**
```
python charged_residues.py -i input.fasta
```

### **Aliphatic Index**
```
python aliphatic_index.py -i input.fasta
```

### **Hydropathy Plot**
```
python hydropathy_plot.py -i input.fasta
```

### **GRAVY**
```
python gravy.py -i input.fasta
```

### **Transmembrane Prediction**
```
python transm_pred.py -i input.fasta
```

### **Instability Index**
```
python instability_index.py -i input.fasta
```

### **Isoelectric Point**
```
python isoelectric_point.py -i input.fasta
```

### **Molecular Weight**
```
python molecular_weight.py -i input.fasta
```

### **Shannon Entropy**
```
python shannon_entropy.py -i input.fasta
```

### **In-vivo Half-life**
```
python half_life.py -i input.fasta
```

### **Extinction Coefficient**
```
python extinction_coefficient.py -i input.fasta
```

### **Toxicity**
```
python toxicity.py -i input.fasta [-t THRESHOLD] [-p PARTS]
```
- `-t THRESHOLD` (default: 0.38) specifies the threshold value (0 to 1).
- `-p PARTS` (default: 50) defines segment length.

### **Structure Prediction**
```
python structure_pred.py -i input.fasta [-p IDENTITY] [-e EVALUE]
```
- `-p IDENTITY` (default: 90.0) sets the minimum percent identity.
- `-e EVALUE` (default: 0.01) defines the e-value for BLAST.

---

## Input & Output

- **Input File:** Users must provide **FASTA** format files (single or multi-sequence).
- **Output Folder:** Each tool automatically generates an output folder where result files are stored.

---

For any issues, please contact the us.
