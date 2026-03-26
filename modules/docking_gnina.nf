
process DOCKING_GNINA {
    label 'gpu_intensive'
    cpus 4
    
    errorStrategy 'ignore'
    conda "conda-forge::python=3.10"

    input:
    path ligandos
    path receptor
    path ligando_control

    output:
    path "resultados_lote_*.csv", emit: reporte_csv
    path "poses_exitosas/*.sdf", emit: poses_3d, optional: true

    script:
    """
    wget -qO gnina https://github.com/gnina/gnina/releases/download/v1.0.3/gnina
    chmod +x gnina

    mkdir -p poses_exitosas
    CSV_FILE="resultados_lote_\${HOSTNAME}_\$\$.csv"
    echo "SMILES_ID,Vina_Affinity_kcal_mol,CNN_Pose_Score,CNN_Affinity_pKd" > "\$CSV_FILE"

    for lig in ${ligandos}; do
        basename=\$(basename "\$lig" .pdbqt)

        ./gnina -r ${receptor} \\
              -l "\$lig" \\
              --autobox_ligand ${ligando_control} \\
              --exhaustiveness 8 \\
              --cnn_scoring rescore \\
              --out "\${basename}_docked.sdf" > "\${basename}_log.txt" 2>&1 || true

        # EL SCRIPT DE PYTHON CORREGIDO
        python3 -c "
import os
log_file = '\${basename}_log.txt'
sdf_file = '\${basename}_docked.sdf'
csv_file = '\$CSV_FILE'
basename = '\${basename}'

try:
    if os.path.exists(log_file):
        with open(log_file, 'r') as f:
            lines = f.readlines()
            
        start_idx = 0
        for i, line in enumerate(lines):
            # Ahora busca la cabecera real sin importar las columnas intermedias
            if 'mode |' in line and 'affinity' in line:
                start_idx = i + 3
                break
                
        if start_idx > 0 and start_idx < len(lines):
            top_pose = lines[start_idx].split()
            
            # Índices corregidos basándonos en tu captura de pantalla
            vina_affinity = float(top_pose[1])
            cnn_score = float(top_pose[2])
            cnn_affinity = float(top_pose[3])
            
            with open(csv_file, 'a') as f_csv:
                f_csv.write(f'{basename},{vina_affinity},{cnn_score},{cnn_affinity}\\n')
            
            if vina_affinity <= -7.0 and cnn_score >= 0.5:
                if os.path.exists(sdf_file):
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