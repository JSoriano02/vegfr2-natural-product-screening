process PREPARE_MEEKO_3D {
    // Allocate 2 CPU threads for computationally intensive 3D geometry calculations
    cpus 2

    // Define conda environment with RDKit, Meeko, Python, and OpenBabel
    conda "conda-forge::rdkit conda-forge::meeko conda-forge::openbabel conda-forge::python=3.10"

    input:
    // Batch file containing SMILES of safe molecules
    path safe_smi

    output:
    // All generated PDBQT files from this processing batch
    path "*.pdbqt", emit: ligands_pdbqt, optional: true

    script:
    """
    obabel ${safe_smi} -O protonated.smi -p 7.4 2>/dev/null || cp ${safe_smi} protonated.smi

    python3 <<EOF
import sys
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
from meeko import MoleculePreparation

# Suppress RDKit C++ warnings to keep the Nextflow console output clean
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

# Stage 2: Initialize Meeko for PDBQT conversion (REQUIREMENT 2: Gasteiger charges)
preparator = MoleculePreparation()

# Extract batch identifier to ensure molecule tracking
batch_id = str('${safe_smi}').replace('.smi', '').replace('clinically_safe_', '')

total_processed = 0
successful_conversions = 0

with open('protonated.smi', 'r') as f:
    for line in f:
        parts = line.strip().split()
        if not parts:
            continue
            
        smiles_string = parts[0]
        # Assign LOTUS ID if available, else use sequential generic ID
        raw_name = parts[1] if len(parts) > 1 else f"ligand_{batch_id}_{total_processed}"
        
        mol = Chem.MolFromSmiles(smiles_string)
        total_processed += 1
        
        if mol is not None:
            try:
                # REQUIREMENT 1: Add explicit hydrogen atoms
                mol = Chem.AddHs(mol)

                # Generate 3D coordinates using ETKDGv3 algorithm (Original User Code)
                result_code = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())

                if result_code == 0:
                    # REQUIREMENT 3: Minimize molecular energy with MMFF94 force field
                    AllChem.MMFFOptimizeMolecule(mol)

                    # Convert optimized structure to PDBQT format
                    preparator.prepare(mol)
                    pdbqt_string = preparator.write_pdbqt_string()

                    # Sanitize filename to prevent Linux filesystem errors
                    safe_filename = "".join(c if c.isalnum() or c in "-_" else "_" for c in raw_name) + ".pdbqt"

                    with open(safe_filename, "w") as out:
                        out.write(pdbqt_string)
                        
                    successful_conversions += 1
            except Exception as e:
                # Silently skip molecules that cannot form valid 3D geometries
                pass

if successful_conversions == 0:
    print(f"CRITICAL ERROR: Batch {batch_id} processed {total_processed} molecules but 0 were successfully converted to 3D. Check formatting.", file=sys.stderr)
    sys.exit(1)

EOF
    """
}