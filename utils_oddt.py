import oddt
from oddt.shape import electroshape
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Descriptors # Descriptors.MolWt( Chem.MolFromInchi( cpd_Olaparib.inchi ) )
from rdkit.Chem.MolStandardize import rdMolStandardize

import numpy as np

#%%
FPS_DIM = 15
max_oddt_bounds =  [10] * FPS_DIM
min_oddt_bounds = [-10] * FPS_DIM

def fp_norm( fps, max_bounds = max_oddt_bounds, min_bounds = min_oddt_bounds ) :
    return [(f - min_bounds[i]) / (max_bounds[i] - min_bounds[i] + 1e-8) for i, f in enumerate(fps)]

def tanimoto_float(fp1, fp2, threshold = 0.5) :
    t = 1.0 - np.sqrt( np.mean( np.abs( np.array(fp1) - np.array(fp2) ) ) )
    if t < threshold:
        t = 0
    return t

def fp_oddt( cpd_struct=None, normalize = True, numConfs = 100, numThreads = 0, largest_fragment = True, mol_weight_thresh = 1000 ):
    '''cpd_struct must be either inchi (preferred) or SMILES (preferably isomeric)'''
    
    if cpd_struct is None:
        raise ValueError("Must supply either smiles or inchi as cpd_struct!")
     
    try:
        if cpd_struct.startswith("InChI="):
            mol_rdkit = Chem.MolFromInchi( cpd_struct ) 
        else:
            mol_rdkit = Chem.MolFromSmiles( cpd_struct )
            
        if largest_fragment:
            largest_Fragment = rdMolStandardize.LargestFragmentChooser()
            mol_rdkit = largest_Fragment.choose(mol_rdkit)
            
        mol_weight = Descriptors.MolWt(mol_rdkit)
        if mol_weight > mol_weight_thresh:
            print("   MW exceeds threshold of", mol_weight_thresh)
            return []
                    
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
        conf_toolkit_rdkit.calccharges()
        fp_float = electroshape(conf_toolkit_rdkit)
        # print('fp_float =', fp_float)
        
        if normalize: # BE CAREFUL (min/max bounds are tricky! use usr_similarity() from oddt whenever possible!)
            oddt_fps.append( fp_norm( fp_float ) )
        else:
            oddt_fps.append( fp_float )
    
    return oddt_fps