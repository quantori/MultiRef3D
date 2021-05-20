#%% imports:
import pandas as pd
from utils_MultiRef3D import load_fps_universe, screen_fps_lib

#%% Load chembl cpds with known MoA
fname_fps = "fp__Enamine_Discovery_Diversity_Set_10K_Library.tsv" # "fp__Enamine_Antiviral_Library.tsv" # "fp__Enamine_Discovery_Diversity_Set_10K_Library.tsv"
# fname_fps = "fp__Enamine_Antiviral_Library.tsv" 
# fname_fps = "fp__chembl_moa.tsv"
cpd_fps_universe = load_fps_universe( fname_fps )

fname_cpd = fname_fps.replace("fp__", "") # "ETL0_dump_moa_via_sqlite.txt"
df_universe = pd.read_csv( fname_cpd, sep = "\t" )
if "idnumber" in df_universe.columns: # Enamine db dump:
    dict_universe = df_universe.rename(columns={'idnumber': 'id_cpd'}).set_index('id_cpd').T.to_dict()
elif "molecule_chembl_id" in df_universe.columns: # ChEMBL db dump:
    dict_universe = df_universe.rename(columns={'molecule_chembl_id': 'id_cpd'}).set_index('id_cpd').T.to_dict()
else:
    raise ValueError("UNKNOWN compounds universe data format!")
print("Loaded fingerprinted universe of", len(dict_universe), "compounds.")

#%% screen cpd universe wrt ref cpds:
list_REFS_NAMES = ["Olaparib", "Tadalafil", "Ergotamine", "Remdesivir"]

# screen fingerprinted lib against chosen refs
df_TopHits, screen_results, screen_dict = screen_fps_lib( 
    cpd_fps_universe = cpd_fps_universe, 
    list_REFS_NAMES = list_REFS_NAMES
)
print("\ndf_TopHits =\n", df_TopHits[:10])

pd.to_pickle(df_TopHits, fname_fps.replace("fp__", "TopHits-" + str(len(list_REFS_NAMES)) + "refs--").replace(".tsv", ".pkl"))

#%% EOF