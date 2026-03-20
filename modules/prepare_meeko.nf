
process PREPARACION_MEEKO_3D {
    // Le damos 2 hilos de CPU porque calcular geometría 3D requiere matemáticas intensivas
    cpus 2
    
    // Entorno Conda con RDKit y Meeko
    conda "conda-forge::rdkit conda-forge::meeko conda-forge::python=3.10"

    input:
    path seguros_smi

    output:
    // Emitimos todos los archivos .pdbqt que se generen en este lote
    path "*.pdbqt", emit: ligandos_pdbqt

    script:
    """
    python3 <<EOF
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation

# 1. Leer las moléculas seguras del lote actual
suppl = Chem.SmilesMolSupplier('${seguros_smi}', titleLine=False, sanitize=True)

# 2. Inicializar el motor de Meeko
preparator = MoleculePreparation()

# Extraemos un identificador del nombre del archivo para no mezclar moléculas de distintos lotes
lote_id = str('${seguros_smi}').replace('.smi', '').replace('seguros_clinicos_', '')

contador = 0
for mol in suppl:
    if mol is not None:
        try:
            # Añadir hidrógenos (vital para simular puentes de hidrógeno)
            mol = Chem.AddHs(mol)
            
            # Plegado 3D (Algoritmo ETKDGv3)
            res = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
            
            # Si res == 0, el algoritmo encontró una geometría 3D estable
            if res == 0:
                # Minimización de energía (Acomodar los átomos para que no choquen)
                AllChem.MMFFOptimizeMolecule(mol)
                
                # Conversión a formato PDBQT con cargas parciales
                preparator.prepare(mol)
                pdbqt_string = preparator.write_pdbqt_string()
                
                # Nombramos el archivo: ej. ligando_lote_ab_0.pdbqt
                mol_name = mol.GetProp("_Name") if mol.HasProp("_Name") and mol.GetProp("_Name") else f"ligando_{lote_id}_{contador}"
                
                with open(f"{mol_name}.pdbqt", "w") as f:
                    f.write(pdbqt_string)
                
                contador += 1
        except Exception as e:
            # Si la molécula es físicamente imposible de doblar, la descartamos en silencio
            pass
EOF
    """
}