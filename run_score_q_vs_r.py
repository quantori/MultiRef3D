import argparse
from utils_oddt import fp_oddt, tanimoto_float

'''This is how to get structural info on reference compounds:'''
'''COMMENT OUT IF USING YOUR OWN QUERY AND REFERENCE COMPOUNDS!'''

DEMO_GET_REFS_FROM_PUBCHEM = False
if DEMO_GET_REFS_FROM_PUBCHEM:
    import pubchempy as pcp
    print("Getting structural info for reference compounds...")
    cpd_cloroquine  = pcp.Compound.from_cid( 2719 ) # CCN(CC)CCCC(C)NC1=C2C=CC(=CC2=NC=C1)Cl
    cpd_remdesivir  = pcp.Compound.from_cid( 121304016 ) # CCC(CC)COC(=O)C(C)NP(=O)(OCC1C(C(C(O1)(C#N)C2=CC=C3N2N=CN=C3N)O)O)OC4=CC=CC=C4
    cpd_favipiravir = pcp.Compound.from_cid( 492405 ) # C1=C(N=C(C(=O)N1)C(=O)N)F
    cpd_jq1         = pcp.Compound.from_cid( 46907787 ) # CC1=C(SC2=C1C(=NC(C3=NN=C(N32)C)CC(=O)OC(C)(C)C)C4=CC=C(C=C4)Cl)C
    cpd_apicidin    = pcp.Compound.from_cid( 6918328 ) # CCC(C)C1C(=O)N2CCCCC2C(=O)NC(C(=O)NC(C(=O)N1)CC3=CN(C4=CC=CC=C43)OC)CCCCCC(=O)CC
    cpd_haloperidol = pcp.Compound.from_cid( 3559 ) # C1CN(CCC1(C2=CC=C(C=C2)Cl)O)CCCC(=O)C3=CC=C(C=C3)F
    
    '''Example of how to extract SMILES or InChI from fetched compound object:'''
    smiles_jq1 = cpd_jq1.isomeric_smiles # isomeric SMILES
    inchi_jq1  = cpd_jq1.inchi # InChI

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    
    '''If using InChi as an input from command prompt, you might need to enclose it in double quotes (depending on the OS)''' 
    parser.add_argument('-q', '--query_structure', type=str, help='SMILES or InChI for the Query compound', required = True)
    parser.add_argument('-r', '--ref_structure', type=str, help='SMILES or InChI for the Reference compound', required = True)
    
    args = parser.parse_args()
    q_STRUCT = args.query_structure # O=C(O)c1cccnc1Oc1cccnc1
    r_STRUCT = args.ref_structure   # CCN(CC)CCCC(C)NC1=C2C=CC(=CC2=NC=C1)Cl
    
    TOTAL_SCORE = 0
    
    print("Fingerprinting query structure...")
    q_fps = fp_oddt( q_STRUCT )
    print("Fingerprinting reference structure...")
    r_fps = fp_oddt( r_STRUCT )
    
    for idx_q, q_fp in enumerate( q_fps ):
        if idx_q % 10 == 0:
            print("Evaluating q_conformer", idx_q+1, "of 100...")
            
        local_MAX = 0
        for idx_r, r_fp in enumerate( r_fps ):
            t_score = tanimoto_float( q_fp, r_fp )
            if t_score > local_MAX:
                local_MAX = t_score
            
        TOTAL_SCORE += local_MAX
            
    print("\nTOTAL_SCORE =", round(TOTAL_SCORE,2) )

