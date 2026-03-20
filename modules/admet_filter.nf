

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

    # 3. Filtrado con columnas exactas y umbral de probabilidad (Pegado al margen izquierdo)
    python3 <<EOF
import pandas as pd
df = pd.read_csv('predicciones.csv')

# Nombres exactos extraídos de tu terminal
# Usamos < 0.5 porque ADMET-AI devuelve la probabilidad de que sea tóxico (clase 1)
mask = (df['BBB_Martins'] < 0.5) & (df['AMES'] < 0.5) & (df['DILI'] < 0.5) & (df['Carcinogens_Lagunin'] < 0.5)

df_seguro = df[mask]

# Guardamos el resultado final
df_seguro.iloc[:, 0].to_csv('seguros_clinicos_${viables_smi}', index=False, header=False)
print("Filtrado completado. Quedan " + str(len(df_seguro)) + " moleculas seguras.")
EOF
    """
}