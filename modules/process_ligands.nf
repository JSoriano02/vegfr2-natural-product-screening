process CHUNK_LOTUS {
    // Lightweight conda environment for Unix utilities
    conda "conda-forge::coreutils"

    input:
    // Complete SMILES database to be split into batches
    path database

    output:
    // Multiple batch files containing subsets of ligands
    path "batch_*"

    script:
    """
    # Split SMILES strings into manageable batch chunks stored asynchronously on disk
    # This parallelizes filteringfor faster processing of large datasets
    split -l ${params.chunk_size} ${database} batch_
    """
}

process FILTER_RDKIT {
    // Conda environment with RDKit and Python for molecular filtering
    conda "conda-forge::rdkit conda-forge::python=3.10"

    input:
    // Single batch file containing SMILES strings to filter
    path batch_smi

    output:
    // Batch file containing only viable molecules that pass Lipinski's rules
    path "${batch_smi}_viable.smi"

    script:
    """
    #!/usr/bin/env python3
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski

    # Open output file for writing viable molecules
    out_file = open('${batch_smi}_viable.smi', 'w')

    # Load all molecules from input batch
    suppl = Chem.MultithreadedSmilesMolSupplier('${batch_smi}', delimiter='\\t', smilesColumn=0, nameColumn=1, titleLine=False)

    # Iterate through molecules and apply Lipinski's Rule of Five filters
    for mol in suppl:
        if mol is not None:
            # Verify molecular weight, lipophilicity, hydrogen donors/acceptors, and rotatable bonds
            if (Descriptors.MolWt(mol) <= 500 and
                Descriptors.MolLogP(mol) <= 5 and
                Lipinski.NumHDonors(mol) <= 5 and
                Lipinski.NumHAcceptors(mol) <= 10 and
                Lipinski.NumRotatableBonds(mol) < 10):

                # Write viable molecule to SMILES format
                out_file.write(f"{Chem.MolToSmiles(mol)}\\t{mol_id}\\n")

    out_file.close()
    """
}