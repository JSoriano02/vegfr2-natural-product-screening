#!/usr/bin/env nextflow

include { OBTAIN_DATA_RAW }                     from './modules/obtain_data.nf'
include { FRACCIONAR_LOTUS; FILTRADO_RDKIT }    from './modules/process_ligands.nf'
include { FILTRADO_ADMET }                      from './modules/admet_filter.nf'
include { PREPARACION_MEEKO_3D }                from './modules/prepare_meeko.nf'
include { PREPARAR_RECEPTOR }                   from './modules/prepare_receptor.nf'

workflow {
    // 1. Ingesta Inicial
    OBTAIN_DATA_RAW(params.lotus_url, params.pdb_id)
    
    lotus_smi_ch    = OBTAIN_DATA_RAW.out.lotus_smi
    receptor_pdb_ch = OBTAIN_DATA_RAW.out.receptor_pdb

    // 2. Fraccionamiento
    lotes_crudos_ch = FRACCIONAR_LOTUS(lotus_smi_ch)

    // 3. Filtrado RDKit
    moleculas_viables_ch = FILTRADO_RDKIT(lotes_crudos_ch.flatten())

    // 4. Filtrado ADMET-AI
    candidatos_seguros_ch = FILTRADO_ADMET(moleculas_viables_ch)

    // 5. Preparación de Ligandos 3D
    ligandos_3d_ch = PREPARACION_MEEKO_3D(candidatos_seguros_ch)

    // 6. Preparación del Receptor y Rescate del Control
    preparacion_receptor = PREPARAR_RECEPTOR(receptor_pdb_ch)
    
    // Separamos las dos salidas en canales distintos para usarlos en el Docking
    receptor_listo_ch  = preparacion_receptor.receptor_pdbqt
    ligando_control_ch = preparacion_receptor.ligando_control
}