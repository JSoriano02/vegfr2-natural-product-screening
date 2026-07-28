process CONSOLIDATE_RESULTS {
    // Publish final results to designated output directory
    publishDir "results/top_final", mode: 'copy'

    // Include RDKit and Pandas for SMILES extraction and ADMET data handling
    conda "conda-forge::python=3.10 conda-forge::rdkit conda-forge::pandas"

    input:
    // Multiple CSV files from docking batches containing scoring results
    path csv_files
    // Multiple SDF files containing successful docking poses
    path sdf_files
    // Multiple CSV files containing ADMET prediction data for safe molecules
    path admet_data_files

    output:
    // Master CSV with top candidates and all metadata
    path "top_candidates_master.csv", emit: master_csv
    // Directory containing SDF files of top ranked poses
    path "top_poses/*.sdf",optional: true , emit: top_sdf

    script:
    """
    mkdir -p top_poses

    python3 -c "
import glob
import csv
import os
import shutil
import pandas as pd
from rdkit import Chem

def safe_get(row, col_name):
    if col_name in row and pd.notna(row[col_name]):
        try:
            return round(float(row[col_name]), 4)
        except:
            return row[col_name]
    return 'N/A'

# Stage 1: Load ADMET data into a dictionary for quick access
admet_dict = {}
for f in glob.glob('admet_data_*.csv'):
    df_admet = pd.read_csv(f)
    for _, row in df_admet.iterrows():
        mol_id = str(row['name'])
        
        # Intentamos buscar el nombre antiguo y el nuevo por si acaso
        carc_val = safe_get(row, 'Carcinogens_Lagunin')
        if carc_val == 'N/A':
            carc_val = safe_get(row, 'Carcinogens')

        # Store critical clinical safety properties safely
        admet_dict[mol_id] = {
            'BBB_Martins': safe_get(row, 'BBB_Martins'),
            'AMES_Toxicity': safe_get(row, 'AMES'),
            'DILI_Risk': safe_get(row, 'DILI'),
            'Carcinogens': carc_val
        }

# Stage 2: Aggregate all docking results from batch CSV files
all_rows = []
header = []

# Read all CSV files and combine into single dataset
for f in glob.glob('*.csv'):
    # Skip ADMET data files to prevent mixing headers with docking results
    if 'admet_data' in f:
        continue
        
    with open(f, 'r') as file:
        reader = csv.reader(file)
        try:
            current_header = next(reader)
        except StopIteration:
            continue # Salta archivos CSV vacíos si los hay
            
        if not header:
            header = current_header
        for row in reader:
            if len(row) >= 4:
                all_rows.append(row)

# Stage 3: Apply ultra-strict filtering for top-tier candidates
filtered_rows = []
for r in all_rows:
    try:
        vina = float(r[1])
        cnn = float(r[2])
        # Criteria: very strong affinity (vina <= -8.0 kcal/mol) AND high AI confidence (cnn >= 0.6)
        if vina <= -8.0 and cnn >= 0.6:
            filtered_rows.append(r)
    except ValueError:
        pass

# Stage 4: Sort candidates by binding affinity (most negative = strongest binding)
filtered_rows.sort(key=lambda x: float(x[1]))

# Stage 5: Select top 50 candidates from all successful predictions
top_molecules = filtered_rows[:50]

# Stage 6: Prepare output headers
header.extend(['BBB_Martins', 'AMES_Toxicity', 'DILI_Risk', 'Carcinogens', 'Extracted_SMILES'])

# Stage 7: Rescue 3D structures, extract SMILES, and merge ADMET data for each top candidate
for row in top_molecules:
    basename = row[0]
    sdf_name = f'{basename}_docked.sdf'
    smiles_string = 'ERROR_NOT_FOUND'

    if os.path.exists(sdf_name):
        # Copy 3D structure file to VIP results directory
        shutil.copy(sdf_name, f'top_poses/{sdf_name}')

        # Extract SMILES representation from 3D SDF structure using RDKit
        try:
            suppl = Chem.SDMolSupplier(sdf_name)
            mol = suppl[0]
            if mol is not None:
                smiles_string = Chem.MolToSmiles(mol)
        except:
            pass

    # Fetch ADMET properties from dictionary using the molecule ID
    admet_info = admet_dict.get(basename, {
        'BBB_Martins': 'N/A', 
        'AMES_Toxicity': 'N/A', 
        'DILI_Risk': 'N/A', 
        'Carcinogens': 'N/A'
    })

    # Append ADMET properties and extracted SMILES to molecule's data row
    row.extend([
        admet_info['BBB_Martins'],
        admet_info['AMES_Toxicity'],
        admet_info['DILI_Risk'],
        admet_info['Carcinogens'],
        smiles_string
    ])

# Stage 8: Write comprehensive master CSV with all information
with open('top_candidates_master.csv', 'w', newline='') as out:
    writer = csv.writer(out)
    writer.writerow(header)
    writer.writerows(top_molecules)
"
    """
}