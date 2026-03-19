
process FILTRADO_ADMET {
    label 'gpu_intensive'
    conda "conda-forge::python=3.10 conda-forge::pandas conda-forge::pip"

    input:
    path viables_smi

    output:
    path "seguros_clinicos_${viables_smi}", emit: candidatos_seguros

    script:
    """
    pip install admet-ai --quiet

    # 1. Preparar entrada
    echo "smiles" > input_con_cabecera.csv
    cat ${viables_smi} >> input_con_cabecera.csv

    # 2. Predicción de IA
    admet_predict \\
        --data_path input_con_cabecera.csv \\
        --save_path predicciones.csv

    # 3. Filtrado dinámico (IMPORTANTE: Pegado al margen izquierdo)
    python3 <<EOF
import pandas as pd
df = pd.read_csv('predicciones.csv')

# Vemos qué columnas ha generado realmente la IA (saldrá en tu log)
print("Columnas detectadas:", df.columns.tolist())

# Definimos los filtros buscando el nombre que contenga la palabra clave
# Esto evita el KeyError si el nombre cambia ligeramente
col_bbb = [c for c in df.columns if 'BBB' in c][0]
col_ames = [c for c in df.columns if 'AMES' in c][0]

# Para hepatotoxicidad y carcinogenicidad, buscamos coincidencias comunes
col_hepato = [c for c in df.columns if 'H-HT' in c or 'Hepatotoxicity' in c][0]
col_carcino = [c for c in df.columns if 'Carcinogenicity' in c][0]

mask = (df[col_bbb] == 0) & (df[col_ames] == 0) & (df[col_hepato] == 0) & (df[col_carcino] == 0)
df_seguro = df[mask]

df_seguro.iloc[:, 0].to_csv('seguros_clinicos_${viables_smi}', index=False, header=False)
print(f"Filtrado completado. Quedan {len(df_seguro)} moleculas seguras.")
EOF
    """
}