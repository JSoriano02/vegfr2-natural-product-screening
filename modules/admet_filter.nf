process FILTER_ADMET {
    // Use GPU for intensive ADMET property prediction
    label 'gpu_intensive'
    // Define conda environment with Python and data analysis tools
    conda "conda-forge::python=3.10 conda-forge::pandas conda-forge::pip"

    input:
    // Batch file containing RDKit-filtered viable SMILES molecules
    path viable_smi

    output:
    // Batch file containing only clinically safe molecules after ADMET filtering
    path "clinically_safe_${viable_smi}", emit: safe_candidates

    script:
    """
    // Install ADMET-AI package for toxicity and bioavailability prediction
    pip install admet-ai --quiet

    // Stage 1: Prepare input with proper CSV header format
    // ADMET-AI expects a header row named "smiles"
    echo "smiles" > input_with_header.csv
    cat ${viable_smi} >> input_with_header.csv

    // Stage 2: Run ADMET-AI prediction on all molecules
    // Generates predictions for multiple toxicity and bioavailability properties
    admet_predict \\
        --data_path input_with_header.csv \\
        --save_path predictions.csv

    // Stage 3: Filter molecules based on safety thresholds
    // Keep only molecules predicted to be safe across all critical ADMET properties
    python3 <<EOF
import pandas as pd

// Load predictions from ADMET-AI
df = pd.read_csv('predictions.csv')

// Apply strict filters for clinical safety:
// - BBB_Martins: Blood-brain barrier permeability (< 0.5 means limited brain exposure)
// - AMES: Mutagenicity test (< 0.5 means non-mutagenic)
// - DILI: Drug-induced liver injury risk (< 0.5 means low risk)
// - Carcinogens_Lagunin: Carcinogenicity prediction (< 0.5 means non-carcinogenic)
// Note: ADMET-AI returns probability of BEING toxic (class 1), so we filter for < 0.5
mask = (df['BBB_Martins'] < 0.5) & (df['AMES'] < 0.5) & (df['DILI'] < 0.5) & (df['Carcinogens_Lagunin'] < 0.5)

// Apply filters to get safe molecules
df_safe = df[mask]

// Write output file containing only SMILES of safe molecules
df_safe.iloc[:, 0].to_csv('clinically_safe_${viable_smi}', index=False, header=False)
print("Filtering complete. Remaining safe molecules: " + str(len(df_safe)))
EOF
    """
}