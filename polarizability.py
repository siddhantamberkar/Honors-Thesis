# Example: Python3 ./polarizability.py "CCCN" 
    ##    Isotropic polarizability: 40.68418 Å³

import os
import sys
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem

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

def write_orca_input(xyz, filename="mol.inp"):
    with open(filename, "w") as f:
        f.write("! HF def2-SVP TightSCF\n")
        f.write("%elprop Polar 1 end\n")
        f.write(f"* xyz 0 1\n{xyz}\n*\n")

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
                parts = line.strip().split()
                return parts[-1]
    return None

def main():
    if len(sys.argv) != 2:
        print("Usage: python3 polarizability.py <SMILES>")
        sys.exit(1)

    smiles = sys.argv[1]
    print(f"[INFO] SMILES: {smiles}")

    xyz = smiles_to_xyz(smiles)
    print(f"[INFO] XYZ block:\n{xyz}")

    write_orca_input(xyz)

  # Ensure ORCA is in your PATH or modify this
    orca_path = os.getenv("ORCA_PATH", "orca")
    print(f"[INFO] Running ORCA from: {orca_path}")

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