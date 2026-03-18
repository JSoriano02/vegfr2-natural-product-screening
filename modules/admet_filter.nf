
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

    admet_predict \\
        --data_path ${viables_smi} \\
        --save_path predicciones.csv \\
        --smiles_column 0

    # 2. Filtrado de Seguridad Clínica
    python3 <<EOF
        import pandas as pd
        df = pd.read_csv('predicciones.csv')

        # Filtros: No cruza cerebro, no mutagénico, no hepatotóxico, no carcinogénico
        mask = (df['BBB_Martins'] == 0) & (df['AMES'] == 0) & (df['H-HT'] == 0) & (df['Rat_Carcinogenicity'] == 0)
        df_seguro = df[mask]

        # Guardamos el resultado final (SMILES está en la primera columna)
        df_seguro.iloc[:, 0].to_csv('seguros_clinicos_${viables_smi}', index=False, header=False)
        EOF
    """
}