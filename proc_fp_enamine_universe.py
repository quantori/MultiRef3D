#%%
'''This script fingerprints the universe for screening (either "anitviral-like" or "Diverse Discovery")
'''
#%% imports:
import pandas as pd
from pathlib import Path
import sys
sys.stdout.flush()

from utils_oddt import fp_oddt

#%% Load structures for the universe:
fname_cpd_universe = "Enamine_Discovery_Diversity_Set_10K_Library.tsv" # downloaded from Enamine
# fname_cpd_universe = "Enamine_Antivirial_Library.tsv" # downloaded from Enamine
fname_fps = "fp__" + fname_cpd_universe

with Path( fname_fps ) as myfile:
    myfile.touch(exist_ok=True)

df_universe = pd.read_csv( fname_cpd_universe, sep = "\t" )

df_universe.rename(columns={'idnumber': 'id_cpd'}, inplace=True)
dict_universe = df_universe.set_index('id_cpd').T.to_dict() # 'records')

#%% Now read cpd names that are already FP'd:
cpd_processed = set()
try:
    df_fpd = pd.read_csv( fname_fps, sep = "\t", header = None )
    for conf_name in list( set( df_fpd[0] ) ):
        cpd_processed.add( conf_name.split("-")[0] )
except:
    print("unable to read_csv() file -->", fname_fps)
    print("presuming empty...")

#%% fp them / save to file:
idx_cpd = -1
for id_cpd, cpd_rec in dict_universe.items():
    idx_cpd += 1
    print("Processing cpd#", idx_cpd+1, " | cid =", id_cpd)
    
    if id_cpd in cpd_processed:
        print("   Skipping (already FPd) -->", id_cpd)
        continue
    
    cpd_struct = cpd_rec[ "smiles" ]
    fp_original = fp_oddt( cpd_struct=cpd_struct, normalize = False )

    with open(fname_fps, "a") as f_out:
        if len(fp_original) == 0:
            f_out.write( str(id_cpd) + "-conf0" + "\t" + str( [0]*15 ) + "\n" )
            f_out.flush()
            continue
        for id_fp, fp in enumerate(fp_original):
            id_conf = str(id_cpd) + "-conf" + str(id_fp)
            str_fp =  str([f for f in fp])
            # print( id_conf, ":", str_fp )
            f_out.write( id_conf + "\t" + str_fp + "\n" )
            f_out.flush()
    
print("\Done!")