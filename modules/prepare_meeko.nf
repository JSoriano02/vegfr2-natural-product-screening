
process PREPARE_MEEKO_3D {
    // Allocate 2 CPU threads for computationally intensive 3D geometry calculations
    cpus 2

    // Define conda environment with RDKit, Meeko, and Python dependencies
    conda "conda-forge::rdkit conda-forge::meeko conda-forge::python=3.10"

    input:
    // Batch file containing SMILES of safe molecules from ADMET filtering
    path safe_smi

    output:
    // All generated PDBQT files from this processing batch
    path "*.pdbqt", emit: ligands_pdbqt

    script:
    """
    python3 <<EOF
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation

# Stage 1: Load safe molecules from current batch
# Read SMILES file and create molecule objects in RDKit
suppl = Chem.SmilesMolSupplier('${safe_smi}', delimiter='\\t', smilesColumn=0, nameColumn=1, titleLine=False, sanitize=True)

# Stage 2: Initialize Meeko for PDBQT conversion with partial charges
preparator = MoleculePreparation()

# Extract batch identifier from filename to ensure molecule tracking
# Prevents mixing molecules from different processing batches
batch_id = str('${safe_smi}').replace('.smi', '').replace('clinically_safe_', '')

# Counter for sequential numbering of molecules in this batch
molecule_counter = 0

# Stage 3: Process each molecule through 3D structure generation and optimization
for mol in suppl:
    if mol is not None:
        try:
            # Step 1: Add hydrogen atoms (critical for accurate hydrogen bonding interactions)
            mol = Chem.AddHs(mol)

            # Step 2: Generate 3D coordinates using ETKDGv3 algorithm
            # This algorithm finds a 3D conformation that minimizes steric clashes
            result_code = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())

            # Step 3: Check if 3D embedding was successful (result_code == 0)
            if result_code == 0:
                # Step 4: Minimize molecular energy with MMFF94 force field
                # Adjusts atomic positions to reduce steric collisions and strain
                AllChem.MMFFOptimizeMolecule(mol)

                # Step 5: Convert optimized structure to PDBQT format with Gasteiger charges
                preparator.prepare(mol)
                pdbqt_string = preparator.write_pdbqt_string()

                # Step 6: Generate filename using molecule name or batch identifier
                # Format: ligand_[batch_id]_[counter].pdbqt
                mol_name = mol.GetProp("_Name") if mol.HasProp("_Name") and mol.GetProp("_Name") else f"ligand_{batch_id}_{molecule_counter}"

                # Step 7: Write PDBQT file for docking
                with open(f"{mol_name}.pdbqt", "w") as f:
                    f.write(pdbqt_string)

                molecule_counter += 1
        except Exception as e:
            # Silently skip molecules that cannot form valid 3D geometries
            # This includes molecules with impossible steric arrangements
            pass
EOF
    """
}