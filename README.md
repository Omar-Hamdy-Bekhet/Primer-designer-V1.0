# Primer-designer-V1.0
An end-to-end Python in-silico pipeline for codon optimization, automated primer design, QC checks, and peptide physicochemical profiling using RDKit.

# 🧬 Bioinformatics & Molecular Design Toolkit
> **An end-to-end Python pipeline integrating Bioinformatics, Primer Engineering, Cheminformatics, and 3D Visualization.**

---

## 📌 Overview
This project is an experimental and educational toolkit designed to explore how different computational techniques can be integrated into a unified molecular biology workflow.

The pipeline starts from an optimized DNA sequence and executes a full multi-stage computational analysis:

$$ \text{DNA} \longrightarrow \text{Primer Design} \longrightarrow \text{Tm/GC QC} \longrightarrow \text{Smart Optimization} \longrightarrow \text{RNA} \longrightarrow \text{Protein} \longrightarrow \text{RDKit Profiling} \longrightarrow \text{3D Rendering} $$

---

## ⚠️ Scientific Disclaimer
> **Note:** This project is intended for educational, exploratory, and computational purposes. Its results should be treated as in-silico predictions rather than experimental laboratory validations.

---

## ✨ Key Features

### 🧬 1. DNA Sequence Analysis
* Calculates optimized gene sequence length.
* Adds flanking sequences and restriction enzyme recognition sites ($NdeI$ & $XhoI$).
* Generates a complete construct sequence and reverse-complement strands.

| Enzyme | Recognition Site | Position |
| :--- | :--- | :--- |
| **NdeI** | `CATATG` | 5' End |
| **XhoI** | `CTCGAG` | 3' End |

---

### 🧪 2. Automated Primer Design & QC
Generates primers flanked with restriction sites and calculates exact Nearest-Neighbor melting temperatures ($T_m$) using **Biopython**:
* **Forward Primer:** `Flank + NdeI + Gene Binding Region`
* **Reverse Primer:** `Flank + XhoI + Reverse Complement Region`
* **QC Assessment:** Evaluates GC content against ideal range ($40\% - 60\%$) with status warnings (**Excellent** / **Accepted** / **Warning**).

---

### 🤖 3. Smart Primer Optimization
* Iterative search algorithm spanning **18–30 bp**.
* Dynamically optimizes binding lengths based on $T_m \ge 54^\circ\text{C}$ and optimal GC thresholds.

---

### 📊 4. Data & Chart Visualization
* Graphical analysis of Forward/Reverse $T_m$ and GC% using **Matplotlib** and **Seaborn**.
* Comparative analysis against manual external input (e.g., Benchling).

---

### 🧫 5. Central Dogma Translation Pipeline
* **DNA $\rightarrow$ mRNA:** Automated transcription ($T \rightarrow U$).
* **mRNA $\rightarrow$ Protein:** Custom codon table translation into single-letter amino acid sequences.

---

### ⚗️ 6. Cheminformatics with RDKit
Connects sequence biology to chemical properties by constructing continuous peptide bonds (`Chem.MolFromSequence`) and calculating physicochemical descriptors:
* Molecular Weight (g/mol)
* LogP (Lipophilicity)
* H-Bond Donors & Acceptors
* Rotatable Bonds & TPSA ($\text{\AA}^2$)
* **Gasteiger Partial Charges** (Max Positive / Min Negative)

---

### 🧊 7. Interactive 3D Visualization
Dual 3D rendering powered by **py3Dmol**:
* **Enzyme Level:** Cartoon representation with VDW surfaces for **enhanced PETase** (PDB: `6EQM`).
* **Molecule Level:** Interactive stick-and-sphere 3D conformers for RDKit molecules.

---

## 🛠️ Tech Stack & Dependencies

```bash
pip install py3Dmol rdkit biopython matplotlib seaborn pandas numpy
LibraryPurposeBiopythonSequence manipulation & Thermodynamic $T_m$ calculationRDKitCheminformatics, peptide building, and partial chargespy3DmolInteractive 3D structural renderingMatplotlib / SeabornData plotting & Quality Control charts🔬 Scientific Pipeline Architecture                    ┌─────────────────────────┐
                    │  Optimized DNA Sequence │
                    └────────────┬────────────┘
                                 │
                                 ▼
                    ┌─────────────────────────┐
                    │   Primer Design & QC    │
                    └────────────┬────────────┘
                                 │
                     ┌───────────┴───────────┐
                     ▼                       ▼
            ┌─────────────────┐     ┌─────────────────┐
            │   Tm (°C) QC    │     │      GC %       │
            └────────┬────────┘     └────────┬────────┘
                     │                       │
                     └───────────┬───────────┘
                                 │
                                 ▼
                    ┌─────────────────────────┐
                    │   Smart Optimization    │
                    └────────────┬────────────┘
                                 │
                                 ▼
                    ┌─────────────────────────┐
                    │    RNA & Translation    │
                    └────────────┬────────────┘
                                 │
                                 ▼
                    ┌─────────────────────────┐
                    │  RDKit Cheminformatics  │
                    └────────────┬────────────┘
                                 │
                                 ▼
                    ┌─────────────────────────┐
                    │ 3D Structure & Graphics │
                    └─────────────────────────┘
⚠️ Known Limitations & Future RoadmapTranscription Model: Assumes coding strand sequence input.Translation Frame: Currently translates from base 0; future releases will include explicit CDS / ORF boundary extraction.Primer Thermodynamics: Primary calculations use binding regions; full-length overhang thermodynamics planned for V2.Drug-likeness Context: Lipinski rules are provided as RDKit demonstrations rather than validated metrics for full macromolecules.👨‍💻 Author & StatusStatus: Experimental / Active Educational DevelopmentDeveloper: Omar Hamdy BekhetBiopharmaceutical & In-Silico Computational Biology Researcher
py3Dmol

Their respective documentation and scientific communities made this type of computational workflow possible.

Built as an experimental exploration of computational biology, molecular design, and scientific programming. 🧬🐍⚗️
