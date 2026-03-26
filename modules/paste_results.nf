
process CONSOLIDATE_RESULTS {
    // Publish final results to designated output directory
    publishDir "results/top_final", mode: 'copy'

    // Include RDKit for SMILES extraction from 3D structures
    conda "conda-forge::python=3.10 conda-forge::rdkit"

    input:
    // Multiple CSV files from docking batches containing scoring results
    path csv_files
    // Multiple SDF files containing successful docking poses
    path sdf_files

    output:
    // Master CSV with top candidates and all metadata
    path "top_candidates_master.csv", emit: master_csv
    // Directory containing SDF files of top ranked poses
    path "top_poses/*.sdf", emit: top_sdf

    script:
    """
    mkdir -p top_poses

    python3 -c "
import glob
import csv
import os
import shutil
from rdkit import Chem

// Stage 1: Aggregate all docking results from batch CSV files
all_rows = []
header = []

// Read all CSV files and combine into single dataset
for f in glob.glob('*.csv'):
    with open(f, 'r') as file:
        reader = csv.reader(file)
        current_header = next(reader)
        if not header:
            header = current_header
        for row in reader:
            if len(row) >= 4:
                all_rows.append(row)

// Stage 2: Apply ultra-strict filtering for top-tier candidates
// Filter molecules based on binding affinity and CNN confidence
filtered_rows = []
for r in all_rows:
    try:
        vina = float(r[1])
        cnn = float(r[2])
        // Criteria: very strong affinity (vina <= -8.0 kcal/mol) AND high AI confidence (cnn >= 0.6)
        // These are strict filters for high-quality candidates only
        if vina <= -8.0 and cnn >= 0.6:
            filtered_rows.append(r)
    except ValueError:
        pass

// Stage 3: Sort candidates by binding affinity (most negative = strongest binding)
filtered_rows.sort(key=lambda x: float(x[1]))

// Stage 4: Select top 50 candidates from all successful predictions
top_molecules = filtered_rows[:50]

// Stage 5: Prepare output by adding SMILES column
// Add header for SMILES structures extracted from 3D poses
header.append('SMILES_Estructura')

// Stage 6: Rescue 3D structures and extract SMILES for each top candidate
for row in top_molecules:
    basename = row[0]
    sdf_name = f'{basename}_docked.sdf'
    smiles_string = 'ERROR_NO_ENCONTRADO'  // Keeps original error message format

    if os.path.exists(sdf_name):
        // Copy 3D structure file to VIP results directory
        shutil.copy(sdf_name, f'top_poses/{sdf_name}')

        // Extract SMILES representation from 3D SDF structure using RDKit
        try:
            // Load 3D molecule structure from SDF file
            suppl = Chem.SDMolSupplier(sdf_name)
            mol = suppl[0]
            if mol is not None:
                // Convert 3D structure to canonical SMILES notation
                smiles_string = Chem.MolToSmiles(mol)
        except:
            // Skip SMILES extraction if RDKit fails
            pass

    // Append extracted SMILES to molecule's data row
    row.append(smiles_string)

// Stage 7: Write comprehensive master CSV with all information
with open('top_candidates_master.csv', 'w', newline='') as out:
    writer = csv.writer(out)
    writer.writerow(header)
    writer.writerows(top_molecules)
"
    """
}