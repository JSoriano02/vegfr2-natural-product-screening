
process CONSOLIDAR_RESULTADOS {
    publishDir "results/top_final", mode: 'copy'
    
    // AÑADIDO: Le damos RDKit a este proceso para poder extraer los SMILES
    conda "conda-forge::python=3.10 conda-forge::rdkit"

    input:
    path csv_files
    path sdf_files

    output:
    path "top_candidates_master.csv", emit: master_csv
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

all_rows = []
header = []

# 1. Leer todos los CSV de los lotes
for f in glob.glob('*.csv'):
    with open(f, 'r') as file:
        reader = csv.reader(file)
        current_header = next(reader)
        if not header:
            header = current_header
        for row in reader:
            if len(row) >= 4:
                all_rows.append(row)

# 2. El Filtro Ultra-Estricto
filtered_rows = []
for r in all_rows:
    try:
        vina = float(r[1])
        cnn = float(r[2])
        # Filtro exigente: Afinidad muy fuerte y confianza alta de la IA
        if vina <= -8.0 and cnn >= 0.6:
            filtered_rows.append(r)
    except ValueError:
        pass

# 3. Ordenar alfabéticamente por energía (el más negativo es el primero)
filtered_rows.sort(key=lambda x: float(x[1]))

# 4. Nos quedamos SOLO con los 50 mejores absolutos
top_molecules = filtered_rows[:50]

# 5. AÑADIDO: Añadir la columna SMILES a la cabecera
header.append('SMILES_Estructura')

# 6. Misión de Rescate 3D y Extracción de SMILES
for row in top_molecules:
    basename = row[0]
    sdf_name = f'{basename}_docked.sdf'
    smiles_string = 'ERROR_NO_ENCONTRADO'
    
    if os.path.exists(sdf_name):
        # A) Copiamos el archivo 3D a la carpeta VIP
        shutil.copy(sdf_name, f'top_poses/{sdf_name}')
        
        # B) Usamos RDKit para leer el 3D y extraer el SMILES
        try:
            suppl = Chem.SDMolSupplier(sdf_name)
            mol = suppl[0]
            if mol is not None:
                smiles_string = Chem.MolToSmiles(mol)
        except:
            pass
            
    # Añadimos el SMILES recuperado al final de la fila de esa molécula
    row.append(smiles_string)

# 7. Guardar el archivo maestro con toda la información
with open('top_candidates_master.csv', 'w', newline='') as out:
    writer = csv.writer(out)
    writer.writerow(header)
    writer.writerows(top_molecules)
"
    """
}