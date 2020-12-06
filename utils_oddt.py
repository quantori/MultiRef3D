import oddt
from oddt.shape import electroshape
from rdkit.Chem import AllChem as Chem

import numpy as np

max_oddt_bounds =  [10] * 15
min_oddt_bounds = [-10] * 15

def fp_norm( fps, max_bounds = max_oddt_bounds, min_bounds = min_oddt_bounds ) :
    return [(f - min_bounds[i]) / (max_bounds[i] - min_bounds[i] + 1e-8) for i, f in enumerate(fps)]

def tanimoto_float(fp1, fp2, threshold = 0.5) :
    t = 1.0 - np.sqrt( np.mean( np.abs( np.array(fp1) - np.array(fp2) ) ) )
    if t < threshold:
        t = 0
    return t

def fp_oddt( cpd_struct=None, numConfs = 100, numThreads = 0 ):
    '''cpd_struct must be either inchi (preferred) or SMILES (preferably isomeric)'''
    
    if cpd_struct is None:
        raise ValueError("Must supply either smiles or inchi as cpd_struct!")
     
    try:
        if cpd_struct.startswith("InChI="):
            mol_rdkit = Chem.MolFromInchi( cpd_struct ) 
        else:
            mol_rdkit = Chem.MolFromSmiles( cpd_struct )
                    
        mol_rdkit = Chem.AddHs(mol_rdkit)
        Chem.EmbedMultipleConfs(mol_rdkit, numConfs = numConfs, randomSeed = 1, numThreads = numThreads)
        Chem.MMFFOptimizeMoleculeConfs(mol_rdkit, maxIters=10000, numThreads = numThreads)
    except:
        print("Unable to produce any converged conformers for cpd_struct =", cpd_struct)
        return []
     
    oddt_fps = []
    for conf in mol_rdkit.GetConformers():
        conf = Chem.Conformer(conf)
        conf.SetId(0)
        mol_oddt = mol_rdkit
        mol_oddt.RemoveAllConformers()
        mol_oddt.AddConformer(conf)
     
        conf_oddt = Chem.MolToMolBlock(mol_oddt)
        conf_toolkit_rdkit = oddt.toolkit.readstring('sdf', conf_oddt)
        oddt_fps.append( fp_norm(electroshape(conf_toolkit_rdkit) ) )
    
    return oddt_fps