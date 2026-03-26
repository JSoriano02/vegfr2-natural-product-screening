process OBTAIN_DATA_RAW {
    // Define conda environment with required tools for data acquisition
    conda "conda-forge::wget conda-forge::coreutils"

    // Save data in the raw_data folder, which is not tracked by Git and is ignored in .gitignore
    publishDir "${params.raw_dir}", mode: 'copy'

    input:
    // URL pointing to LOTUS Natural Products Database SMILES file
    val lotus_url
    // Protein Data Bank (PDB) identifier code for the receptor structure
    val pdb_id

    output:
    // Full LOTUS ligand database in SMILES format
    path "lotus_full.smi", emit: lotus_smi
    // Receptor structure file from Protein Data Bank
    path "${pdb_id}.pdb", emit: receptor_pdb

    script:
    """
    // Download the complete LOTUS Natural Products Database containing all available ligands
    wget -O lotus_full.smi ${lotus_url}

    // Download the receptor 3D structure from the Protein Data Bank (PDB)
    wget -O ${pdb_id}.pdb "https://files.rcsb.org/download/${pdb_id}.pdb"
    """
}