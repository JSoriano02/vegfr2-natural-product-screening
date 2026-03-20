#!/usr/bin/env nextflow

include { OBTAIN_DATA_RAW }                     from './modules/obtain_data.nf'
include { FRACCIONAR_LOTUS; FILTRADO_RDKIT }    from './modules/process_ligands.nf'
include { FILTRADO_ADMET }                      from './modules/admet_filter.nf'
include { PREPARACION_MEEKO_3D }                from './modules/prepare_meeko.nf'
workflow {
    // 1. Ingesta Inicial
    OBTAIN_DATA_RAW(params.lotus_url, params.pdb_id)
    
    // 2. Extracción de canales
    lotus_smi_ch    = OBTAIN_DATA_RAW.out.lotus_smi
    receptor_pdb_ch = OBTAIN_DATA_RAW.out.receptor_pdb

    // 3. Fraccionamiento del archivo gigante
    lotes_crudos_ch = FRACCIONAR_LOTUS(lotus_smi_ch)

    // 4. Filtrado Farmacológico inicial (Lipinski)
    moleculas_viables_ch = FILTRADO_RDKIT(lotes_crudos_ch.flatten())

    // 5. Filtrado de Seguridad Clínica con IA (ADMET-AI)
    candidatos_seguros_ch = FILTRADO_ADMET(moleculas_viables_ch)

    // 6. Conversión a 3D y Preparación para Docking (RDKit + Meeko)
    ligandos_3d_ch = PREPARACION_MEEKO_3D(candidatos_seguros_ch)
}