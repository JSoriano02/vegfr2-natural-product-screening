process PREPARE_RECEPTOR {
    // Use single CPU thread for sequential processing
    cpus 1

    // Use Bioconda channel for the official and stable MGLTools package
    conda "bioconda::mgltools"

    input:
    // Receptor structure file from Protein Data Bank
    path receptor_pdb

    output:
    // Receptor in PDBQT format (AutoDock compatible atom types and Kollman charges)
    path "ready_receptor.pdbqt", emit: receptor_pdbqt
    // Co-crystalized ligand extracted from receptor structure
    path "control_ligand.pdb", emit: control_ligand

    script:
    """
    # Stage 1: Extract and prepare co-crystalized ligand
    awk '/^HETATM/ && !/HOH/' ${receptor_pdb} > control_ligand.pdb

    # Stage 2: Extract protein structure
    awk '/^ATOM/' ${receptor_pdb} > receptor_clean.pdb

    # Stage 3: Prepare AutoDock receptor using the official MGLTools script
    prepare_receptor4.py -r receptor_clean.pdb -o ready_receptor.pdbqt -A hydrogens -U nphs_lps_waters
    """
}