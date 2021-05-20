import os
import numpy as np
import pandas as pd
from collections import defaultdict

import pubchempy as pcp
from oddt.shape import usr_similarity
from rdkit.Chem import AllChem as Chem
import oddt
from oddt.shape import electroshape
from utils_oddt import fp_oddt

#%%
def load_fps_REFS( fname_dict_fps_REFS ):
    '''NOTE: current reference list is limited by names/PubChem IDs below (for demo purposes).
    '''
    dict_fps_REFS = dict()
    
    if os.path.exists( fname_dict_fps_REFS ):
        print("Reading from existing file -->", fname_dict_fps_REFS)
        dict_fps_REFS = pd.read_pickle( fname_dict_fps_REFS )

    else:
        print("Fetching and fingerprinting compounds from PubChem...")
        cpd_ref = pcp.Compound.from_cid( 23725625 )
        dict_fps_REFS["Olaparib"] = fp_oddt( cpd_struct = cpd_ref.inchi, normalize = False )
        
        cpd_ref = pcp.Compound.from_cid( 110635 )
        dict_fps_REFS["Tadalafil"] = fp_oddt( cpd_struct = cpd_ref.inchi, normalize = False )
                       
        cpd_ref = pcp.Compound.from_cid( 8223 )
        dict_fps_REFS["Ergotamine"] = fp_oddt( cpd_struct = cpd_ref.inchi, normalize = False )
               
        cpd_ref = pcp.Compound.from_cid( 121304016 )
        dict_fps_REFS["Remdesivir"] = fp_oddt( cpd_struct = cpd_ref.inchi, normalize = False )
         
        pd.to_pickle(dict_fps_REFS, fname_dict_fps_REFS)
        
    print("dict_fps_REFS contains", len(dict_fps_REFS), "compounds")
    return dict_fps_REFS


def load_fps_universe( fname_fp_universe ) :
    '''Read in fp_file to a dict( id -> list of conf fps ) 
    '''
    cpd_fps_universe = defaultdict( list )
    
    with open(fname_fp_universe, "r") as f_in:
        for idx_line, line in enumerate(f_in):
            str_id_conf, str_fp = line.split("\t")
            if "nan" in str_fp:
                continue
            
            fp = np.array( eval(str_fp) )
            id_cpd, str_idx_conf = str_id_conf.split("-")
            str_idx_conf.replace("conf","")
            
            cpd_fps_universe[id_cpd].append( fp )
            
            # if idx_line > 99:
            #     break
        
    cpd_fps_universe = dict( cpd_fps_universe )
    print("Loaded", len(cpd_fps_universe), "FP'd compounds from", fname_fp_universe)
    
    return cpd_fps_universe


def score_qfps_vs_rfps( q_fps, r_fps, tanimoto_thresh = 0.5 ):
    TOTAL_SCORE = 0
    
    for idx_q, q_fp in enumerate( q_fps ):
        if idx_q % 10 == 0:
            # print("Evaluating q_conformer", idx_q+1, "of 100...")
            pass
            
        local_MAX = 0
        for idx_r, r_fp in enumerate( r_fps ):
            # t_score = tanimoto_float( q_fp, r_fp )
            t_score = usr_similarity( q_fp, r_fp )
            if t_score > local_MAX:
                local_MAX = t_score
            
        if local_MAX < tanimoto_thresh:
            local_MAX = 0
            
        TOTAL_SCORE += local_MAX
        
    return TOTAL_SCORE

def screen_fps_lib( cpd_fps_universe = None, 
                   list_REFS_NAMES = None,
                   dict_fps_REFS = None,
                   TANIMOTO_THRESH = 0.5 
                   ) :
    if cpd_fps_universe is None or list_REFS_NAMES is None:
        raise ValueError("MUST SUPPLY both pd_fps_universe, list_REFS and dict_fps_REFS!")
        
    if dict_fps_REFS is None:
        dict_fps_REFS = pd.read_pickle( "dict_fps_REFS.pkl" )
    list_REFS = [dict_fps_REFS[ list_REFS_NAMES[i] ] for i in range( len(list_REFS_NAMES) )]
    
    screen_dict = defaultdict( list ) # dict()
    cpd_count = -1
    for cpd_id, cpd_fps in cpd_fps_universe.items():
        cpd_count += 1
        if cpd_count % 100 == 0:
            print("Screening cpd#", cpd_count)
        
        for fps_ref in list_REFS:
            t_score = score_qfps_vs_rfps( cpd_fps, fps_ref, tanimoto_thresh = TANIMOTO_THRESH )
            # print("t_score =", t_score)
            
            screen_dict[ cpd_id ].append( t_score )
            
        # break
    screen_dict = dict( screen_dict )
    
    screen_results = dict()
    for cpd_id, t_scores in screen_dict.items():
        screen_results[ cpd_id ] = sum( t_scores )
    
    screen_sorted = sorted(screen_results.items(), key=lambda item: item[1], reverse=True)
    # print( screen_sorted[:5] )
    print("Screened", len(screen_dict), "compounds")
    
    # Create df with sorted (Top) scores/compounds
    TopN = len(screen_dict) # 20
    dict_TopHits = dict()
    dict_TopHits[ "cpd_id" ] = [""]*TopN
    dict_TopHits[ "cpd_name" ] = [""]*TopN
    dict_TopHits[ "TotalScore" ] = [0.0]*TopN
    dict_TopHits = {**dict_TopHits, **{ list_REFS_NAMES[i] : [0.0]*TopN for i in range( len(list_REFS_NAMES) ) }}
    
    df_TopHits = pd.DataFrame.from_dict( dict_TopHits )
    
    for idx_top in range(TopN):
        cpd_id = screen_sorted[idx_top][0]
        cpd_score = screen_sorted[idx_top][1]
        
        cpd_name = cpd_id
        
        # if 'molecule_pref_name' in dict_universe[ cpd_id ]:
        #     cpd_name = dict_universe[ cpd_id ][ 'molecule_pref_name' ].capitalize()
        # else:
        #     cpd_name = cpd_id
        
        # print("Cpd:", cpd_id, "|", cpd_name, "| Score:", round(cpd_score,2) ) # 
        # print("Cpd:", cpd_id, "|", cpd_name, "| Score:", str([round(x,4) for x in screen_dict[str(id_test)]]) ) # round(cpd_score,2))
        df_TopHits.at[ idx_top, "cpd_id" ] = cpd_id
        df_TopHits.at[ idx_top, "cpd_name" ] = cpd_name
        df_TopHits.at[ idx_top, "TotalScore" ] = cpd_score
        
        sim_scores = [ round(x,4) for x in screen_dict[ cpd_id] ]
        for i in range(len(sim_scores)):
            df_TopHits.at[ idx_top, list_REFS_NAMES[i] ] = sim_scores[i]
    
    return df_TopHits, screen_results, screen_dict

def gen_fps_universe( input_SMILES, filter_for_tox=True ):
    import QDD_MolGen    
    novel_SMILES = QDD_MolGen.generate_novel_SMILES( input_SMILES, filter_for_tox=filter_for_tox )
    print("   Generated", len(novel_SMILES), "novel compounds.")
    
    cpd_fps_universe = dict()
    print("   Fingerprinting generated universe...")
    for smi in novel_SMILES:
        fps = fp_oddt( cpd_struct=smi, normalize = False )
        cpd_fps_universe[ smi ] = fps
        
    return cpd_fps_universe

#%% EOF