
process PREPARE_RECEPTOR {
    // Use single CPU thread for sequential processing
    cpus 1

    // Use OpenBabel - robust conda package for molecular format conversion
    conda "conda-forge::openbabel"

    input:
    // Receptor structure file from Protein Data Bank
    path receptor_pdb

    output:
    // Receptor in PDBQT format (AutoDock compatible with charges)
    path "ready_receptor.pdbqt", emit: receptor_pdbqt
    // Co-crystalized ligand extracted from receptor structure
    path "control_ligand.pdb", emit: control_ligand

    script:
    """
    // Stage 1: Extract and prepare co-crystalized ligand
    // Extract all HETATM (non-water ligand) records, excluding HOH (water molecules)
    // This control ligand is used to define the docking binding box
    awk '/^HETATM/ && !/HOH/' ${receptor_pdb} > control_ligand.pdb

    // Stage 2: Extract protein structure
    // Keep only ATOM records (protein backbone and sidechains)
    // Exclude all other records (HETATM, CONECT, etc.)
    awk '/^ATOM/' ${receptor_pdb} > receptor_clean.pdb

    // Stage 3: Prepare receptor for docking with OpenBabel
    // -p 7.4    : Add hydrogens calculated for physiological pH 7.4
    // -xr        : AutoDock format flag - calculates Gasteiger partial charges, rigid PDBQT format
    obabel receptor_clean.pdb -O ready_receptor.pdbqt -p 7.4 -xr
    """
}