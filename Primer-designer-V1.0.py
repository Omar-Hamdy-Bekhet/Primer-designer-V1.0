import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import py3Dmol
from Bio.SeqUtils import MeltingTemp as mt
from time import sleep
from os import wait
import sys
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors, Draw, Lipinski


def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement[base] for base in reversed(seq))
#_______________________________________________________________________________
def calculate_exact_Tm(seq, primer_conc=500, na_conc=50):
    return mt.Tm_NN(seq, dnac1=primer_conc, Na=na_conc)
#_______________________________________________________________________________
def calculate_gc(seq):
    g_count = seq.count('G')
    c_count = seq.count('C')
    gc_content = ((g_count + c_count) / len(seq)) * 100
    return gc_content
#_______________________________________________________________________________

print(" ")
print("--------------------------------")
print(" ")
#-------------------------------------- #simple Calculates & restriction enzymes.
optm_seq = 'ATGCAAACTAAACCATTTCGCCGTACTATGCGCAAGCTGGCGACCACCCTGATGTTTACCACCCTGGCGGCAAGCGCGGCGACCGCGTGTGATAGCAACCAGACCGCCCCGGTGGGCATGGCGTATTTTGCGAACGGCCATACCATGTGCCACAGCATTGCGTGGAGCATTACCAGCCTGGGCACCGCGAACGATGTGCTGACCGCGAGCCTGAACAGCGCGGTTAACCTGCCACAAAGCACCCAGGATGCGAGCATTGGCATGAGCATGGGCGGCGGCGCGGCGCTGGGCACCCTGCTGATTAGCCCG'
flank = 'GCTA'

print("Optmized Gene sequence length : ", len(optm_seq))

ndeI_site = "CATATG"   # 5' NdeI site
xhoI_site = "CTCGAG"   # 3' XhoI site
full_Optm_seq = ndeI_site + optm_seq + xhoI_site
print("Fully Optmized Gene sequence : ", full_Optm_seq)
print("Fully Optmized Gene sequence length : ", len(full_Optm_seq))

#--------------------------------------
print(" ")
print("--------------------------------")
print(" ")
#-------------------------------------- # primer Calcs.

fwd_pri = flank + ndeI_site + optm_seq[:18]
rev_pri = flank + xhoI_site + reverse_complement(optm_seq[-18:])

print("Forward Primer (5' -> 3'):", fwd_pri)
print("Forward Primer Length:", len(fwd_pri))
print("Reverse Primer (5' -> 3'):", rev_pri)
print("Reverse Primer Length:", len(rev_pri))

#--------------------------------------
print(" ")
print("--------------------------------")
print(" ")
#-------------------------------------- # Annealing Temperature Calcs + QC & QA For The Primer.

fwd_bind = optm_seq[:18]
rev_bind = reverse_complement(optm_seq[-18:])
fwd_tm = calculate_exact_Tm(fwd_bind)
rev_tm = calculate_exact_Tm(rev_bind)

ta_temp = min(fwd_tm, rev_tm) - 5

print(f"Forward Binding Tm: {fwd_tm:.2f} °C")
print(f"Reverse Binding Tm: {rev_tm:.2f} °C")
print(f"Recommended PCR Annealing Temp (Ta): {ta_temp:.2f} °C")


fwd_gc = calculate_gc(fwd_pri)
rev_gc = calculate_gc(rev_pri)

print("Forward Primer GC Content:", round(fwd_gc, 2), "%")
print("Reverse Primer GC Content:", round(rev_gc, 2), "%")

print(" ")
print("Test 1 Conditions :")
if 40 <= fwd_gc <= 60:
    print("\033[92mStatus: Forward Primer is in Excellent condition, Ready for synthesis\033[0m")
else:
    print("\033[91mStatus: Warning! Forward Primer content is outside the ideal 40-60% range.\033[0m")
if 40 <= rev_gc <= 60:
    print("\033[92mStatus: Reverse Primer is in Excellent condition, Ready for synthesis\033[0m")
else:
    print("\033[91mStatus: Warning! Reverse Primer content is outside the ideal 40-60% range.\033[0m")

# Variables to hold current primer values for plotting, initialized with initial calculations
current_fwd_tm = fwd_tm
current_rev_tm = rev_tm
current_fwd_gc = fwd_gc
current_rev_gc = rev_gc

#--------------------------------------
print(" ")
print("--------------------------------")
print(" ")
#-------------------------------------- # Smart Optmization Mode (Optinal).

min_length = 18
max_length = 30
ELN = input("Do you Want to start the smart Optmization mode ? | y/n (Heavliy Recommended to be ON): ")

#FWDOPM
if ELN == 'y' or ELN == 'Y' or ELN == 'غ':
    print("excuteing Smart Optmization Proccess ...")
    sleep(5)
    print("12%")
    sleep(0.3)
    print("48%")
    sleep(0.3)
    print("73%")
    sleep(0.3)
    print("96%")
    sleep(1.5)
    print("100%")
    for length in range(min_length, max_length + 1):
        test_bind = optm_seq[:length]
        test_tm = calculate_exact_Tm(test_bind)
        test_gc = calculate_gc(flank + ndeI_site + test_bind)
        if test_tm >= 54 and test_gc >= 40:
            opt_fwd_bind = test_bind
            opt_fwd_pri = flank + ndeI_site + opt_fwd_bind
            print(" ")
            print(f"Optimal Forward Primer Length found: {len(opt_fwd_pri)} bp")
            print(f"New Forward Binding Tm: {test_tm:.2f} °C")
            print(f"New Forward GC Content: {round(test_gc, 2)} %")
            current_fwd_tm = test_tm
            current_fwd_gc = test_gc
            break
    #REVOPM
    for length in range(min_length, max_length + 1):
        test_bind2 = reverse_complement(optm_seq[-length:])
        test_tm2 = calculate_exact_Tm(test_bind2)
        test_gc2 = calculate_gc(flank + xhoI_site + test_bind2)
        if test_tm2 >= 54 and test_gc2 <= 61.5:
            opt_rev_bind = test_bind2
            opt_rev_pri = flank + xhoI_site + opt_rev_bind
            print(f"Optimal Reverse Primer Length found: {len(opt_rev_pri)} bp")
            print(f"New Reverse Binding Tm: {test_tm2:.2f} °C")
            print(f"New Reverse GC Content: {round(test_gc2, 2)} %")
            current_rev_tm = test_tm2
            current_rev_gc = test_gc2
            break

    print(" ")
    print("Test 2 Conditions :")
    if 40 <= current_fwd_gc <= 60:
        print("\033[92mStatus: Forward Primer is in Excellent condition, Ready for synthesis\033[0m")
    elif 35 <= current_fwd_gc <= 65:
      print("\033[32mStatus: Forward Primer is in Accepted condition, Good for synthesis\033[0m")
    else:
        print("\033[91mStatus: Warning! Forward Primer content is outside the ideal 40-60% range.\033[0m")
    if 40 <= current_rev_gc <= 60:
        print("\033[92mStatus: Reverse Primer is in Excellent condition, Ready for synthesis\033[0m")
    elif 35 <= current_rev_gc <= 65:
      print("\033[32mStatus: Reverse Primer is in Accepted condition, Good for synthesis\033[0m")
    else:
        print("\033[91mStatus: Warning! Reverse Primer content is outside the ideal 40-60% range.\033[0m")

    ta_temp2 = min(current_fwd_tm, current_rev_tm) - 5
    print(f"New Recommended PCR Annealing Temp (Ta): {ta_temp2:.2f} °C")

else:
  print("Smart Optmization Proccess is Skipped")

#--------------------------------------
print(" ")
print("--------------------------------")
print(" ")
#-------------------------------------- # Manual Input From Benchling.

ELNBE = input("Do you want to manually input Benchling values for plotting ? (Enter n if you dont have outer values) | y/n : ")
if ELNBE == 'y' or ELNBE == 'Y' or ELNBE == 'غ':
    fwd_tm_plot = float(input("Enter Forward Primer Tm from Benchling (°C): "))
    rev_tm_plot = float(input("Enter Reverse Primer Tm from Benchling (°C): "))
    fwd_gc_plot = float(input("Enter Forward Primer GC% from Benchling: "))
    rev_gc_plot = float(input("Enter Reverse Primer GC% from Benchling: "))
else:
    fwd_tm_plot, rev_tm_plot = current_fwd_tm, current_rev_tm
    fwd_gc_plot, rev_gc_plot = current_fwd_gc, current_rev_gc
    print("Manual Input is Skipped")

#--------------------------------------
print(" ")
print("--------------------------------")
print(" ")
#-------------------------------------- # Estimated Protien 3D Structure view.

ELN3D = input("Do you want to see the Estimated Protien 3D Structure ? | y/n : ")
if ELN3D == 'y' or ELN3D == 'Y' or ELN3D == 'غ':
    print("Visualizng a 3D Structure View ...")
    viewer = py3Dmol.view(query='pdb:6EQM') # >>>>>>>>>>>>>>>>>>>>>>>>>>> CHANGE PDB HERE ! <<<<<<<<<<<<<<<<<<<<<<<<< (Current enzyme is : enhanced PETase)
    viewer.setStyle({'cartoon': {'color': 'spectrum'}})
    viewer.addSurface(py3Dmol.VDW, {'opacity': 0.7, 'color': 'white'})
    viewer.zoomTo()
    viewer.show()
else:
    print("3D Structure view is Skipped")

#--------------------------------------
print(" ")
print("--------------------------------")
print(" ")
#-------------------------------------- # Chart Visualzation.

ELNCHA = input("Do you want to start a chart Visualzation ? | y/n : ")
if ELNCHA == 'y' or ELNCHA == 'Y' or ELNCHA == 'غ':
    print("Visualizng a Chart Visualzation ...")
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    axes[0].bar(['Forward Primer', 'Reverse Primer'], [fwd_tm_plot, rev_tm_plot], color=['#3498db', '#e74c3c'])
    axes[0].axhline(54, color='gray', linestyle='--', label='Min Ideal Tm (54°C)')
    axes[0].set_ylabel('Temperature (°C)')
    axes[0].set_title('Primer Binding Melting Temperatures (Tm)')
    axes[0].set_ylim(0, 80)
    axes[0].legend()

    axes[1].bar(['Forward Primer', 'Reverse Primer'], [fwd_gc_plot, rev_gc_plot], color=['#2ecc71', '#9b59b6'])
    axes[1].axhline(40, color='orange', linestyle='--', label='Min GC (40%)')
    axes[1].axhline(60, color='red', linestyle='--', label='Max GC (60%)')
    axes[1].set_ylabel('Percentage (%)')
    axes[1].set_title('Primer GC Content (%)')
    axes[1].set_ylim(0, 100)
    axes[1].legend()

    plt.tight_layout()
    plt.show()

else:
  print("Chart Visualzation is Skipped")

#--------------------------------------
print(" ")
print("--------------------------------")
print(" ")
#-------------------------------------- # DNA to mRNA.

full_Optm_mRNA = full_Optm_seq.replace('T', 'U')

ELNRNA = input("Do you want to Translate DNA sequence to RNA ? | y/n : ")

if ELNRNA == 'y' or ELNRNA == 'Y' or ELNRNA == 'غ':
    print("Fully Optmized mRNA sequence : ", full_Optm_mRNA)
else:
    print("Translation to RNA is Skipped")

#--------------------------------------
print(" ")
print("--------------------------------")
print(" ")
#-------------------------------------- # mRNA To Protein.

codon_table = {'AUA':'I', 'AUC':'I', 'AUU':'I', 'AUG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACU':'T',
               'AAC':'N', 'AAU':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGU':'S', 'AGA':'R', 'AGG':'R',
               'CUA':'L', 'CUC':'L', 'CUG':'L', 'CUU':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCU':'P',
               'CAC':'H', 'CAU':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGU':'R',
               'GUA':'V', 'GUC':'V', 'GUG':'V', 'GUU':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCU':'A',
               'GAC':'D', 'GAU':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGU':'G',
               'UCA':'S', 'UCC':'S', 'UCG':'S', 'UCU':'S', 'UUC':'F', 'UUU':'F', 'UUA':'L', 'UUG':'L',
               'UAC':'Y', 'UAU':'Y', 'UGC':'C', 'UGU':'C', 'UGG':'W', 'UGA':'!!_*', 'UAA':'!!_*', 'UAG':'!!_*'}
amino_array = []

for i in range(0, len(full_Optm_mRNA), 3):
    codon = full_Optm_mRNA[i:i+3]
    amino_acid = codon_table.get(codon, 'X')
    amino_array.append(amino_acid)
    amino_seq = ''.join(amino_array)

ELNPRO = input("Do you want to Translate mRNA sequence to Amino Acids ? | y/n : ")
if ELNPRO == 'y' or ELNPRO == 'Y' or ELNPRO == 'غ':
    print("Fully Optmized Amino Acids sequence : ", amino_seq)
else:
    print("Translation to Amino Acids is Skipped")


#--------------------------------------
print(" ")
print("--------------------------------")
print(" ")
#-------------------------------------- # integreating Biology with Chemistry Using RdKit library.

amino_acid_smiles = {
    'A': 'CC(C(=O)O)N',                     # Alanine (Ala)
    'R': 'C(CC(C(=O)O)N)CN=C(N)N',          # Arginine (Arg)
    'N': 'C(C(C(=O)O)N)C(=O)N',             # Asparagine (Asn)
    'D': 'C(C(C(=O)O)N)C(=O)O',             # Aspartate (Asp)
    'C': 'C(C(C(=O)O)N)S',                  # Cysteine (Cys)
    'E': 'CCC(C(=O)O)NC(=O)O',              # Glutamate (Glu)
    'Q': 'CCC(C(=O)O)NC(=O)N',              # Glutamine (Gln)
    'G': 'NCC(=O)O',                        # Glycine (Gly)
    'H': 'C1=C(NC=N1)CC(C(=O)O)N',          # Histidine (His)
    'I': 'CCC(C)C(C(=O)O)N',                # Isoleucine (Ile)
    'L': 'CC(C)CC(C(=O)O)N',                # Leucine (Leu)
    'K': 'NCCCC[C@@H](C(=O)O)N',            # Lysine (Lys)
    'M': 'CSCCC(C(=O)O)N',                  # Methionine (Met)
    'F': 'c1ccc(cc1)CC(C(=O)O)N',           # Phenylalanine (Phe)
    'P': 'C1CC(NC1)C(=O)O',                 # Proline (Pro)
    'S': 'C(C(C(=O)O)N)O',                  # Serine (Ser)
    'T': 'CC(C(C(=O)O)N)O',                 # Threonine (Thr)
    'W': 'c1ccc2c(c1)c(c[nH]2)CC(C(=O)O)N', # Tryptophan (Trp)
    'Y': 'c1cc(ccc1CC(C(=O)O)N)O',          # Tyrosine (Tyr)
    'V': 'CC(C)C(C(=O)O)N'                  # Valine (Val)
}
smiles_seq=[]

for amino_acid in amino_seq:
    if amino_acid in amino_acid_smiles:
        smiles_seq.append(amino_acid_smiles[amino_acid])
    else:
        smiles_seq.append('X')

mol_connected = Chem.MolFromSequence(amino_seq)

if mol_connected is not None:
    mol = Chem.AddHs(mol_connected)

mw = Descriptors.MolWt(mol)
logp = Descriptors.MolLogP(mol)
HYD = Lipinski.NumHDonors(mol)
HYA = Lipinski.NumHAcceptors(mol)
BOR = Descriptors.NumRotatableBonds(mol)
TPSA = Descriptors.TPSA(mol)
AllChem.ComputeGasteigerCharges(mol)

charges = [float(atom.GetProp('_GasteigerCharge')) for atom in mol.GetAtoms() if atom.HasProp('_GasteigerCharge')]
max_pos_charge = max(charges)
min_neg_charge = min(charges)

print(f"Molecular Weight : {mw:.1f} g/mol")
print(f"LogP             : {logp:.1f}")
print(f"H-Bond Donors    : {HYD}")
print(f"H-Bond Acceptors : {HYA}")
print(f"Rotatable Bonds  : {BOR}")
print(f"TPSA             : {TPSA:.1f} Å²")
print(f"Max Positive Charge (+): {max_pos_charge:.1f}")
print(f"Max Negative Charge (-): {min_neg_charge:.1f}")


LIR = (mw <= 500) and (logp <= 5) and (HYD <= 5) and (HYA <= 10)
print(f"Passes Lipinski Rule ? : {LIR}")

Draw.MolToImage(mol)

mb = Chem.MolToMolBlock(mol)

print("\033[93mWarning! the next step needs a Med~HI PC recources, This step is heavily Unrecommended for weak GPUs\033[0m")
ELNFUS = input("Do you want to start a full Molucular Visualzation ? | y/n : ")
if ELNFUS == 'y' or ELNFUS == 'Y' or ELNFUS == 'غ':
    print("Visualizng a Full Molucular Visualzation ...")
    view = py3Dmol.view(width=400, height=400)
    view.addModel(mb, 'mol')
    view.setStyle({'sphere': {'radius': 0.25, 'colorscheme': 'Jmol'}, 'stick': {'radius': 0.15}})

    #ELNELC = input("Do you want to show electronegativty? | y/n : " )
    #if ELNELC == 'y' or ELNELC == 'Y' or ELNELC == 'غ':
    #    AllChem.MMFFOptimizeMolecule(mol)
    #    for i, atom in enumerate(mol.GetAtoms()):
    #        if atom.HasProp('_GasteigerCharge'):
    #            chg = float(atom.GetProp('_GasteigerCharge'))
    #            pos = mol.GetConformer().GetAtomPosition(i)
    #            view.addLabel(f"{chg:+.2f}", {'position': {'x': pos.x, 'y': pos.y, 'z': pos.z}, 'backgroundColor': 'black', 'fontColor': 'white', 'fontSize': 10, 'backgroundOpacity': 0.7})

    view.zoomTo()
    view.show()
else:
    print("Full Molucular Visualzation is Skipped")

#--------------------------------------
print(" ")
print("--------------------------------")
print(" ")
#-------------------------------------- #
