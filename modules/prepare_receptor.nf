
process PREPARAR_RECEPTOR {
    cpus 1
    
    // Usamos OpenBabel, que es moderno y nunca falla en Conda
    conda "conda-forge::openbabel"

    input:
    path receptor_pdb

    output:
    path "receptor_listo.pdbqt", emit: receptor_pdbqt
    path "ligando_control.pdb", emit: ligando_control

    script:
    """
    # 1. RESCATE DEL CONTROL Y LIMPIEZA BASICA (Usando 'awk' que es a prueba de fallos)
    # Guardamos los ligandos ignorando el agua
    awk '/^HETATM/ && !/HOH/' ${receptor_pdb} > ligando_control.pdb
    
    # Guardamos SOLO la proteína (líneas ATOM), descartando todo lo demás
    awk '/^ATOM/' ${receptor_pdb} > receptor_clean.pdb

    # 2. PREPARACIÓN FINAL CON OPENBABEL
    # -p 7.4 : Añade hidrógenos calculados para un pH fisiológico de 7.4
    # -xr    : Flag especial de AutoDock (calcula cargas Gasteiger y formatea como PDBQT rígido)
    obabel receptor_clean.pdb -O receptor_listo.pdbqt -p 7.4 -xr
    """
}