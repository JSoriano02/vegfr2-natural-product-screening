process OBTAIN_DATA_RAW {
    conda "conda-forge::wget conda-forge::coreutils"
    
    // Save data in the raw_data folder, which is not tracked by Git and is ignored in .gitignore
    publishDir "${params.raw_dir}", mode: 'copy'

    input:
    val lotus_url
    val pdb_id

    output:
    path "lotus_full.smi", emit: lotus_smi
    path "${pdb_id}.pdb", emit: receptor_pdb

    script:
    """
    # Obtention of the full dataset from LOTUS Natural Products Database
    wget -O lotus_full.smi ${lotus_url}
    
    # Obtention of the receptor structure from the Protein Data Bank (PDB)
    wget -O ${pdb_id}.pdb "https://files.rcsb.org/download/${pdb_id}.pdb"
    """
}