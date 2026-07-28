#!/usr/bin/env nextflow

// Main workflow orchestration for VGFR2 molecular docking pipeline
// Processes include: data acquisition, ligand filtering, receptor preparation, and molecular docking

include { OBTAIN_DATA_RAW }              from './modules/obtain_data.nf'
include { CHUNK_LOTUS; FILTER_RDKIT }   from './modules/process_ligands.nf'
include { FILTER_ADMET }                from './modules/admet_filter.nf'
include { PREPARE_MEEKO_3D }            from './modules/prepare_meeko.nf'
include { PREPARE_RECEPTOR }            from './modules/prepare_receptor.nf'
include { DOCKING_GNINA }               from './modules/docking_gnina.nf'
include { CONSOLIDATE_RESULTS }         from './modules/paste_results.nf'

workflow {
    // Stage 1: Data acquisition from LOTUS database and PDB structure database
    OBTAIN_DATA_RAW(params.lotus_url, params.pdb_id)
    lotus_smi_ch    = OBTAIN_DATA_RAW.out.lotus_smi
    receptor_pdb_ch = OBTAIN_DATA_RAW.out.receptor_pdb

    // Stage 2: Ligand database chunking and filtering
    // Chunk large LOTUS dataset into manageable batches
    raw_batches_ch = CHUNK_LOTUS(lotus_smi_ch)

    // Apply RDKit filtering to remove molecules violating Lipinski's rules
    viable_molecules_ch = FILTER_RDKIT(raw_batches_ch.flatten())

    // Filter viable molecules by ADMET properties for safety
    admet_results = FILTER_ADMET(viable_molecules_ch)
    safe_candidates_ch = admet_results.safe_candidates
    admet_data_ch      = admet_results.admet_data

    // Stage 3: 3D structure generation and optimization
    meeko_results = PREPARE_MEEKO_3D(safe_candidates_ch)
    ligands_3d_ch = meeko_results.ligands_pdbqt 

    // Stage 4: Receptor preparation
    receptor_preparation = PREPARE_RECEPTOR(receptor_pdb_ch)
    
    ready_receptor_ch = receptor_preparation.receptor_pdbqt
    control_ligand_ch = receptor_preparation.control_ligand

    // Stage 5: Molecular docking with GNINA
    final_results_ch = DOCKING_GNINA(
        ligands_3d_ch,
        ready_receptor_ch, 
        control_ligand_ch  
    )

    // Stage 6: Results consolidation
    // Collect, filter top candidates, and extract SMILES representations
    CONSOLIDATE_RESULTS(
        final_results_ch.reporte_csv.collect(),
        final_results_ch.poses_3d.collect(),
        admet_data_ch.collect()
   )
}