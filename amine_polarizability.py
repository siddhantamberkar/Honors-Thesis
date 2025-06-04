# Example: Python3 ./polarizability.py "CCCN" 
    ##   [INFO] Matched: CCCN is a primary aliphatic amine.
    ##   [RESULT] Isotropic polarizability: 40.68418 Å³
# Note: If the inputed SMILE string does not correspond to an amine, the calculation will be canceled. 

import sys
import os
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem

# SMARTS patterns for amines
patterns = {
    "primary aliphatic amine": "[NX3;H2][CX4]",
    "primary aromatic amine": "[NX3;H2][c]",
    "secondary aliphatic amine": "[NX3;H1]([CX4])[CX4]",
    "secondary aromatic amine": "[NX3;H1]([CX4])[c]",
    "tertiary aliphatic amine": "[NX3;H0]([CX4])([CX4])[CX4]",
    "tertiary aromatic amine": "[NX3;H0]([CX4])([CX4])[c]",
    "quaternary amine": "[NX4+]"
}

def identify_amine(smiles):
    mol = Chem.MolFromSmiles(smiles)
    for name, smarts in patterns.items():
        pattern = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(pattern):
            print(f"[INFO] Matched: {smiles} is a {name}.")
            return True
    print(f"[INFO] No amine match found for {smiles}.")
    return False

def smiles_to_xyz(smiles, name="mol"):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.UFFOptimizeMolecule(mol)

    conf = mol.GetConformer()
    xyz = ""
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        xyz += f"{atom.GetSymbol()} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}\n"
    return xyz

def write_orca_input(xyz_block, filename="mol.inp"):
    with open(filename, "w") as f:
        f.write("! HF def2-SVP TightSCF\n")
        f.write("%elprop Polar 1 end\n")
        f.write("* xyz 0 1\n")
        f.write(xyz_block)
        f.write("*\n")

def run_orca(input_file="mol.inp", output_file="mol.out", orca_path="orca"):
    try:
        with open(output_file, "w") as out:
            subprocess.run([orca_path, input_file], stdout=out, stderr=subprocess.STDOUT, check=True)
    except FileNotFoundError:
        print("ERROR: ORCA executable not found. Please set the ORCA_PATH environment variable or add ORCA to your system PATH.")
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        print("ORCA calculation failed:")
        print(e.stderr)
        sys.exit(1)

def validate_orca_run(output_file="mol.out"):
    with open(output_file) as f:
        for line in f:
            if "ORCA TERMINATED NORMALLY" in line:
                return True
    return False

def extract_isotropic_polarizability(orca_output_file):
    with open(orca_output_file, 'r') as f:
        for line in f:
            if "Isotropic polarizability" in line:
                return line.strip().split()[-1]
    return None

def main():
    if len(sys.argv) != 2:
        print("Usage: ./polarizability.py 'SMILES_STRING'")
        sys.exit(1)

    smiles = sys.argv[1]

    if not identify_amine(smiles):
        sys.exit(0)  # Exit cleanly if not an amine

    xyz = smiles_to_xyz(smiles)
    write_orca_input(xyz)

    orca_path = os.getenv("ORCA_PATH", "orca")
    run_orca(orca_path=orca_path)

    if not validate_orca_run():
        print("[ERROR] ORCA did not terminate normally.")
        sys.exit(1)

    polar = extract_isotropic_polarizability("mol.out")
    if polar:
        print(f"[RESULT] Isotropic polarizability: {polar} Å³")
    else:
        print("[ERROR] Polarizability not found.")

if __name__ == "__main__":
    main()
