
process DOCKING_GNINA {
    label 'gpu_intensive'
    cpus 4
    
    publishDir "results/docking_final", mode: 'copy'
    
    conda "conda-forge::gnina conda-forge::python=3.10"

    input:
    tuple path(ligandos), path(receptor), path(ligando_control)

    output:
    path "resultados_lote_*.csv", emit: reporte_csv
    path "poses_exitosas/*.sdf", emit: poses_3d, optional: true

    script:
    """
    mkdir -p poses_exitosas
    CSV_FILE="resultados_lote_\${HOSTNAME}_\$\$.csv"
    echo "SMILES_ID,Vina_Affinity_kcal_mol,CNN_Pose_Score,CNN_Affinity_pKd" > "\$CSV_FILE"

    # Iteramos sobre cada ligando del lote
    for lig in ${ligandos}; do
        basename=\$(basename "\$lig" .pdbqt)

        gnina -r ${receptor} \\
              -l "\$lig" \\
              --autobox_ligand ${ligando_control} \\
              --exhaustiveness 8 \\
              --cnn_scoring rescore \\
              --out "\${basename}_docked.sdf" > "\${basename}_log.txt"

        python3 -c "
import os
log_file = '\${basename}_log.txt'
sdf_file = '\${basename}_docked.sdf'
csv_file = '\$CSV_FILE'
basename = '\${basename}'

try:
    with open(log_file, 'r') as f:
        lines = f.readlines()
        
    start_idx = 0
    for i, line in enumerate(lines):
        if 'mode |   affinity |  dist from best mode |  CNN |   CNN' in line:
            start_idx = i + 3
            break
            
    if start_idx > 0 and start_idx < len(lines):
        top_pose = lines[start_idx].split()
        
        vina_affinity = float(top_pose[1])
        cnn_score = float(top_pose[4])
        cnn_affinity = float(top_pose[5])
        
        if vina_affinity <= -7.0 and cnn_score >= 0.5:
            with open(csv_file, 'a') as f_csv:
                f_csv.write(f'{basename},{vina_affinity},{cnn_score},{cnn_affinity}\\n')
            os.rename(sdf_file, f'poses_exitosas/{basename}_docked.sdf')
        else:
            if os.path.exists(sdf_file):
                os.remove(sdf_file)
except Exception as e:
    pass
"
        rm -f "\${basename}_log.txt" "\${basename}_docked.sdf"
    done
    """
}