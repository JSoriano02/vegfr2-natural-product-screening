#!/usr/bin/env nextflow

include { OBTAIN_DATA_RAW }                   from './modules/obtain_data.nf'
include { FRACCIONAR_LOTUS; FILTRADO_RDKIT }    from './modules/process_ligands.nf'
include { FILTRADO_ADMET }                      from './modules/admet_filter.nf'

workflow {
    // 1. Ingesta Inicial. Ejecutamos el proceso directamente
    OBTAIN_DATA_RAW(params.lotus_url, params.pdb_id)
    
    // 2. Extraemos los canales usando la sintaxis estricta ".out"
    lotus_smi_ch    = OBTAIN_DATA_RAW.out.lotus_smi
    receptor_pdb_ch = OBTAIN_DATA_RAW.out.receptor_pdb

    // 3. Fraccionamiento Asíncrono
    lotes_crudos_ch = FRACCIONAR_LOTUS(lotus_smi_ch)

    // 4. Filtrado Farmacológico en paralelo
    moleculas_viables_ch = FILTRADO_RDKIT(lotes_crudos_ch.flatten())

    // 5. Filtrado de Seguridad Clínica con IA (ADMET-AI)
    candidatos_seguros_ch = FILTRADO_ADMET(moleculas_viables_ch)

}