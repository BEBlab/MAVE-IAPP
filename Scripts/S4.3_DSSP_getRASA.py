#This script is needed to the relative ASA (accessible surface area) of each residue resolved in a PDB structure. 

#Command to run it from the terminal: python3 S4.3_DSSP_getRASA.py 
#PDB files where the relative ASA needs to be extracted must be in the same folder together with this script. 
#This script will generate separate .csv files for each PDB structure in the folder and one .csv file "merged_ss_sasa.csv" where all the relative ASAs can be found.
#"merged_ss_sasa.csv" is the input file needed to run the script S4.4_Correlation of nucleation scores and relative ASA.R.

import os
import pandas as pd
import numpy as np
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.DSSP import DSSP
import warnings

warnings.filterwarnings('ignore')

# Define a dictionary for amino acids
d_aa = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H',
        'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
        'TYR': 'Y', 'VAL': 'V'}

# Initialize the PDB parser
parser = PDBParser()
structures_names = [i for i in os.listdir() if i.endswith('pdb')]

def change_bfactor(pdb, new_bfactors, out_name, bfactor):
    try:
        # Parse the pdb file
        structure = parser.get_structure(pdb, pdb)
        for chain in structure[0].get_chains():
            for residue in chain:
                for atom in residue:
                    atom.set_bfactor(0)
                for i in new_bfactors.index:
                    if d_aa[residue.resname] + str(residue.id[1]) == new_bfactors.loc[i, 'AAPos']:
                        for atom in residue:
                            atom.set_bfactor(new_bfactors.loc[i, bfactor])
        io = PDBIO()
        io.set_structure(structure)
        io.save(out_name + ".pdb")
        return structure
    except Exception as e:
        print(f"Error processing {pdb}: {e}")

def aa_ss_sasa_pdb(pdb):
    try:
        structure = parser.get_structure(pdb, pdb)
        dssp = DSSP(structure[0], pdb)
        AA, Pos, AAPos, sec_str, sasa = [], [], [], [], []
        
        for residue in structure[0]['A']:
            AA.append(d_aa[residue.resname])
            Pos.append(residue.id[1])
            AAPos.append(d_aa[residue.resname] + str(residue.id[1]))
        
        for pos in Pos:
            sasa_in = []
            for chain in structure[0].get_chains():
                try:
                    sasa_in.append(dssp[(chain.id, (' ', pos, ' '))][3])
                except KeyError:
                    pass
            sec_str.append(dssp[('A', (' ', pos, ' '))][2])
            sasa.append(np.mean(sasa_in) if sasa_in else np.nan)
        
        df = pd.DataFrame({'AAPos': AAPos, f'SecStr_{pdb.split(".pdb")[0]}': sec_str, f'SASA_{pdb.split(".pdb")[0]}': sasa})
        df.to_csv(f'ss_sasa_{pdb.split(".pdb")[0]}.csv', index=False)
        return df
    except Exception as e:
        print(f"Error processing {pdb}: {e}")

# Calculate ss and sasa for each structure
for pdb in structures_names:
    aa_ss_sasa_pdb(pdb)

# Load the first CSV file
df = pd.read_csv(f'ss_sasa_{structures_names[0].split(".pdb")[0]}.csv')

# Merge the remaining CSV files
for pdb in structures_names[1:]:
    df_other = pd.read_csv(f'ss_sasa_{pdb.split(".pdb")[0]}.csv')
    df = df.merge(df_other, how='outer', on='AAPos')

df.to_csv('merged_ss_sasa.csv', index=False)
print(df)