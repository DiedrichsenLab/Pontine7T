import pandas as pd
from pathlib import Path
import numpy as np
import torch as pt
import nibabel as nb
import SUITPy as suit
import matplotlib.pyplot as plt
import seaborn as sb
from copy import copy,deepcopy
import Functional_Fusion.atlas_map as am
import Functional_Fusion.dataset as ds
import Functional_Fusion.matrix as matrix
import HierarchBayesParcel.full_model as fm
import HierarchBayesParcel.spatial as sp
import HierarchBayesParcel.arrangements as ar
import HierarchBayesParcel.emissions as em
import HierarchBayesParcel.evaluation as ev
import DCBC.utilities as ut
from scripts import indiv_parc
from DCBC.dcbc import compute_DCBC 

wk_dir = '/Volumes/diedrichsen_data$/data/Cerebellum/Pontine7T/atlases/thalamus'
data_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion_new/MDTB'

def evaluate_dcbc(wta_indiv_data,wta_indiv_group,U_group,atlas='MNISymThalamus1',max_dist=40):
    
    atlas_obj, _ = am.get_atlas(atlas)
    coords = atlas_obj.world.T

    wta_group = np.argmax(U_group, axis=0) + 1

    dist = ut.compute_dist(coords,1.5, backend='numpy')

    T = pd.DataFrame()

    subj_info = pd.read_csv(f'{data_dir}/participants.tsv',sep='\t')
    
    for i, s in enumerate(subj_info['participant_id']):
        
        # Clean up re
        data,info,_ = ds.get_dataset(base_dir,'MDTB',session='sc2',atlas='MNI',type='CondAll')
        
        dcbc_group = compute_DCBC(maxDist=max_dist, binWidth=1.5, parcellation=wta_group,
                                func= data, dist=dist, weighting=True, backend='numpy')

        D1 = {}
        D1['type'] = ['group']
        D1['runs'] = [0]
        D1['dcbc'] = [dcbc_group['DCBC']]
        D1['subject'] = [i + 1]
        T = pd.concat([T, pd.DataFrame(D1)])

        for r in range(len(wta_indiv_data)):
            dcbc_indiv_data = compute_DCBC(maxDist=max_dist, binWidth=1.5, parcellation=wta_indiv_data[r],
                                            func= data, dist=dist, weighting=True, backend='numpy')
            dcbc_indiv_group = compute_DCBC(maxDist=max_dist, binWidth=1.5, parcellation=wta_indiv_group[r],
                                            func= data, dist=dist, weighting=True, backend='numpy')
            
            D1 = {}
            D1['type'] = ['data']
            D1['runs'] = [r + 1]
            D1['dcbc'] = [dcbc_indiv_data['DCBC']]
            D1['subject'] = [i + 1]
            T = pd.concat([T, pd.DataFrame(D1)])
            
            D1 = {}
            D1['type'] = ['data and group']
            D1['runs'] = [r + 1]
            D1['dcbc'] = [dcbc_indiv_group['DCBC']]
            D1['subject'] = [i + 1]
            T = pd.concat([T, pd.DataFrame(D1)])

    return T 

if __name__ == "__main__":


    data_full = np.load(f'{wk_dir}/indiv_parcellations/MDTB-ses2/MDTB_ses2_data.npy')
    data_section = data_full[:3]
    U_indiv_data = np.load(f'{wk_dir}/indiv_parcellations/MDTB-ses1_flat_prior/all_indiv_parcels.npy')
    U_indiv_group = np.load(f'{wk_dir}/indiv_parcellations/MDTB-ses1/all_indiv_parcels.npy')
    U_group = np.load(f'{wk_dir}/Prob_thalamus.npy')

    D = evaluate_dcbc(U_indiv_data,U_indiv_group,U_group,atlas='MNISymThalamus1',max_dist=60)


    








