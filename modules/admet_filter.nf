process FILTER_ADMET {
    // Use GPU for intensive ADMET property prediction
    label 'gpu_intensive'
    // Define conda environment with Python and data analysis tools
    conda "conda-forge::python=3.10 conda-forge::pip conda-forge::pandas"

    input:
    // Batch file containing RDKit-filtered viable SMILES molecules
    path viable_smi

    output:
    // Batch file containing only clinically safe molecules after ADMET filtering and their corresponding ADMET prediction data
    path "clinically_safe_${viable_smi}", emit: safe_candidates
    path "admet_data_*.csv", emit: admet_data 

    script:
    """
    # 1. Force install an older setuptools (<70) because hyperopt requires pkg_resources
    pip install "setuptools<70.0.0" --quiet
    
    # 2. Install ADMET-AI
    pip install admet-ai --quiet

    # Stage 1: Prepare input with proper CSV header format
    # ADMET-AI expects a header row named "smiles"
    python3 <<EOF
import pandas as pd
df_in = pd.read_csv('${viable_smi}', sep='\\t', header=None, names=['smiles', 'name'])
df_in.to_csv('input_with_header.csv', index=False)
EOF

    # Stage 2: Run ADMET-AI prediction on all molecules
    # Generates predictions for multiple toxicity and bioavailability properties
    admet_predict \\
        --data_path input_with_header.csv \\
        --save_path predictions.csv

    # Stage 3: Filter molecules based on OPHTHALMIC clinical safety thresholds
    python3 <<EOF
import pandas as pd

# Load predictions from ADMET-AI
df = pd.read_csv('predictions.csv')

# Dynamic column naming (handles slight variations in ADMET-AI output names)
col_eye = 'Eye_Irritation' if 'Eye_Irritation' in df.columns else ('Eye_toxicity' if 'Eye_toxicity' in df.columns else 'Eye_Corrosion')
col_neuro = 'Neurotoxicity'
col_immuno = 'Immunotoxicity'
col_ames = 'AMES'
col_hepa = 'Hepatotoxicity' if 'Hepatotoxicity' in df.columns else 'DILI'
col_bbb = 'BBB_Martins'

# Apply strict filters for INTRAVITREAL INJECTION clinical safety:
# - BBB_Martins: Used as proxy for Blood-Retinal Barrier (must NOT leak to brain/blood)
# - Eye/Neuro/Immunotoxicity: Critical local dealbreakers for the eye
# - AMES/Hepatotoxicity: Systemic safety parameters
mask = (
    (df.get(col_eye, 0) < 0.5) & 
    (df.get(col_neuro, 0) < 0.5) & 
    (df.get(col_immuno, 0) < 0.5) &
    (df.get(col_ames, 0) < 0.5) & 
    (df.get(col_hepa, 0) < 0.5) &
    (df.get(col_bbb, 0) < 0.5)
)

# Apply filters to get safe molecules
df_safe = df[mask]

# Write output file containing only SMILES of safe molecules for the next Nextflow stage
df_safe[['smiles', 'name']].to_csv('clinically_safe_${viable_smi}', sep='\\t', index=False, header=False)

# Define the comprehensive pharmacokinetic parameters we want to save for the final report
columns_to_keep = [
    'smiles', 'name', 'Log(D)', col_bbb, 'CYP2D6_Inhibitor', 
    'CYP3A4_Inhibitor', 'Clearance', col_eye, col_hepa, 
    col_neuro, 'Nephrotoxicity', col_immuno, col_ames
]

# Keep only columns that actually exist to avoid KeyError
existing_cols = [c for c in columns_to_keep if c in df_safe.columns]

# Save the rich ADMET prediction data for these safe molecules
clean_name = "${viable_smi}".replace('.smi', '')
df_safe[existing_cols].to_csv(f'admet_data_{clean_name}.csv', index=False)
EOF
    """
}