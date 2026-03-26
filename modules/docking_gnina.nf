
process DOCKING_GNINA {
    // Use GPU for intensive docking computations
    label 'gpu_intensive'
    // Allocate 4 CPU threads for parallel operations
    cpus 4

    // Continue processing other tasks even if individual dockings fail
    errorStrategy 'ignore'
    // Define conda environment with Python for results parsing
    conda "conda-forge::python=3.10"

    input:
    // Multiple ligand files in PDBQT format to be docked
    path ligands
    // Receptor protein structure in PDBQT format
    path receptor
    // Control ligand for defining docking box coordinates
    path control_ligand

    output:
    // CSV files containing docking scores for each batch processed
    path "results_batch_*.csv", emit: reporte_csv
    // SDF files with successful docking poses (optional, may not exist)
    path "successful_poses/*.sdf", emit: poses_3d, optional: true

    script:
    """
    // Stage 1: Download and prepare GNINA docking engine
    // GNINA combines Vina docking with CNN-based pose scoring
    wget -qO gnina https://github.com/gnina/gnina/releases/download/v1.0.3/gnina
    chmod +x gnina

    // Stage 2: Initialize output directories and CSV file
    mkdir -p successful_poses
    // Create results file with timestamp and process ID to avoid conflicts
    CSV_FILE="results_batch_\${HOSTNAME}_\$\$.csv"
    // Write CSV header with docking and scoring metrics
    echo "SMILES_ID,Vina_Affinity_kcal_mol,CNN_Pose_Score,CNN_Affinity_pKd" > "\$CSV_FILE"

    // Stage 3: Dock each ligand against the receptor
    for lig in ${ligands}; do
        // Extract molecule identifier from filename
        basename=\$(basename "\$lig" .pdbqt)

        // Run GNINA docking with CNN-based scoring
        // --exhaustiveness 8: balance between speed and accuracy
        // --cnn_scoring rescore: use CNN to rescore poses from Vina
        ./gnina -r ${receptor} \\
              -l "\$lig" \\
              --autobox_ligand ${control_ligand} \\
              --exhaustiveness 8 \\
              --cnn_scoring rescore \\
              --out "\${basename}_docked.sdf" > "\${basename}_log.txt" 2>&1 || true

        // Stage 4: Parse docking results and extract scoring metrics
        // Python script for processing GNINA output logs and results
        python3 -c "
import os
log_file = '\${basename}_log.txt'
sdf_file = '\${basename}_docked.sdf'
csv_file = '\$CSV_FILE'
basename = '\${basename}'

try:
    // Read log file to extract docking results
    if os.path.exists(log_file):
        with open(log_file, 'r') as f:
            lines = f.readlines()

        // Find start of scoring table (after header lines with column names)
        start_idx = 0
        for i, line in enumerate(lines):
            // Look for table header row containing mode and affinity columns
            if 'mode |' in line and 'affinity' in line:
                start_idx = i + 3  // Skip header separator lines
                break

        // Extract top pose (first scoring line after header)
        if start_idx > 0 and start_idx < len(lines):
            top_pose = lines[start_idx].split()

            // Extract metrics from the parsed line
            // Index 1: Vina binding affinity (kcal/mol)
            // Index 2: CNN pose scoring confidence
            // Index 3: CNN binding affinity prediction (pKd)
            vina_affinity = float(top_pose[1])
            cnn_score = float(top_pose[2])
            cnn_affinity = float(top_pose[3])

            // Write results to CSV file
            with open(csv_file, 'a') as f_csv:
                f_csv.write(f'{basename},{vina_affinity},{cnn_score},{cnn_affinity}\\n')

            // Stage 5: Keep successful poses that pass quality filters
            // Criteria: Vina affinity <= -7.0 kcal/mol AND CNN confidence >= 0.5
            if vina_affinity <= -7.0 and cnn_score >= 0.5:
                if os.path.exists(sdf_file):
                    // Move successful pose to dedicated directory
                    os.rename(sdf_file, f'successful_poses/{basename}_docked.sdf')
            else:
                // Discard unsuccessful poses to save disk space
                if os.path.exists(sdf_file):
                    os.remove(sdf_file)
except Exception as e:
    // Silently skip any parsing errors
    pass
"
        // Cleanup temporary files
        rm -f "\${basename}_log.txt" "\${basename}_docked.sdf"
    done
    """
}