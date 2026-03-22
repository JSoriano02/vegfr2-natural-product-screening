#!/usr/bin/env nextflow

include { OBTAIN_DATA_RAW }                     from './modules/obtain_data.nf'
include { FRACCIONAR_LOTUS; FILTRADO_RDKIT }    from './modules/process_ligands.nf'
include { FILTRADO_ADMET }                      from './modules/admet_filter.nf'
include { PREPARACION_MEEKO_3D }                from './modules/prepare_meeko.nf'
include { PREPARAR_RECEPTOR }                   from './modules/prepare_receptor.nf'
include { DOCKING_GNINA }                       from './modules/docking_gnina.nf'

workflow {
    OBTAIN_DATA_RAW(params.lotus_url, params.pdb_id)
    lotus_smi_ch    = OBTAIN_DATA_RAW.out.lotus_smi
    receptor_pdb_ch = OBTAIN_DATA_RAW.out.receptor_pdb

    lotes_crudos_ch = FRACCIONAR_LOTUS(lotus_smi_ch)
    moleculas_viables_ch = FILTRADO_RDKIT(lotes_crudos_ch.flatten())
    candidatos_seguros_ch = FILTRADO_ADMET(moleculas_viables_ch)
    ligandos_3d_ch = PREPARACION_MEEKO_3D(candidatos_seguros_ch)

    // 6. Preparación del Receptor y Rescate del Control
    preparacion_receptor = PREPARAR_RECEPTOR(receptor_pdb_ch)
    
    receptor_listo_ch  = preparacion_receptor.receptor_pdbqt
    ligando_control_ch = preparacion_receptor.ligando_control

    // 7. EL DOCKING FINAL
    resultados_finales_ch = DOCKING_GNINA(
        ligandos_3d_ch, 
        receptor_listo_ch, 
        ligando_control_ch
    )
}