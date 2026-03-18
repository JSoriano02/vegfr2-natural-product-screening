process FRACCIONAR_LOTUS {
    conda "conda-forge::coreutils"

    input:
    path base_datos

    output:
    path "lote_*"

    script:
    """
    # Desglose de cadenas de texto SMILES en bloques asíncronos en disco
    split -l ${params.chunk_size} ${base_datos} lote_
    """
}

process FILTRADO_RDKIT {
    conda "conda-forge::rdkit conda-forge::python=3.10"

    input:
    path lote_smi

    output:
    path "${lote_smi}_viable.smi"

    script:
    """
    #!/usr/bin/env python3
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski

    out_file = open('${lote_smi}_viable.smi', 'w')
    
    suppl = Chem.MultithreadedSmilesMolSupplier('${lote_smi}', delimiter='\\t', smilesColumn=0, titleLine=False)

    for mol in suppl:
        if mol is not None:
            if (Descriptors.MolWt(mol) <= 500 and 
                Descriptors.MolLogP(mol) <= 5 and 
                Lipinski.NumHDonors(mol) <= 5 and 
                Lipinski.NumHAcceptors(mol) <= 10 and 
                Lipinski.NumRotatableBonds(mol) < 10):
                
                out_file.write(Chem.MolToSmiles(mol) + '\\n')
                
    out_file.close()
    """
}