import pandas as pd
import os 
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
import HierarchBayesParcel.util as ut_hierarch
import HierarchBayesParcel.evaluation as ev
import DCBC.utilities as ut
from scripts import indiv_parc
from DCBC.dcbc import compute_DCBC 

wk_dir = '/Volumes/diedrichsen_data$/data/Cerebellum/Pontine7T/atlases/thalamus'
data_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion_new/MDTB'

def get_group_atlas(wk_dir = wk_dir, atlas_name = 'MNISymThalamus1', roi = 'thalamus'):
    atlas, _ = am.get_atlas(atlas_name)
    atlas_group = np.load(f"{wk_dir}/Prob_{roi}.npy")
    ar_model = ar.build_arrangement_model(atlas_group, prior_type='prob', atlas=atlas)

    return ar_model, atlas 

def get_data(data_dir = data_dir, atlas_name = 'MNISymThalamus1', session = 'localizerfm'):
    subj_info = pd.read_csv(f'{data_dir}/participants.tsv',sep='\t')
    data = []

    for i, s in enumerate(subj_info['participant_id']):

        file_name = f'{data_dir}/derivatives/ffextract/{s}/{s}_space-{atlas_name}_ses-{session}_CondRun.dscalar.nii'

        #MDTB - {s}_space-{atlas_name.name}_ses-s1_CondRun.dscalar.nii (or ses-s2)
        #Language - {s}_space-{atlas_name.name}_ses-localizerfm_CondRun.dscalar.nii
        #MDTB-high-res - {s}_space-{atlas_name.name}_ses-s1_CondRun.dscalar.nii
        
        if not os.path.exists(file_name):
            print(f"File {file_name} does not exist. Skipping subject {s}.")
            continue

        datafile = nb.load(file_name)
        data.append(datafile.get_fdata())

    data_np = np.stack(data)

    # data is: numsubj x numcond (repetitions x conditions) x numvoxel tensor

    data = pt.tensor(data_np, dtype=pt.get_default_dtype())


    return data

def get_info_emission(data_dir = data_dir, sample_subj = 'sub_01', session = 'localizerfm', dataset = 'Language'):

    info = pd.read_csv(f'{data_dir}/derivatives/ffextract/{sample_subj}/{sample_subj}_ses-{session}_CondRun.tsv',sep='\t')

    #MDTB - sub-02_ses-s1_CondRun.tsv (cond_name)
    #Language - sub-01_ses-localizerfm_CondRun.tsv (task_name)
    #MDTB-high-res - sub-01_ses-s1_CondRun.tsv (task_name)

    if dataset in ['Language', 'MDTB-high-res']:
        cond_v = info['task_name'].values
    if dataset in ['MDTB']:
        cond_v = info['cond_name'].values 

    part_v = info['run'].values

    #data = ds.remove_baseline(data,part_v) #for language, pontine, do this separately 

    return cond_v, part_v

def build_indiv_parc_runwise(ar_model, atlas, data, cond_v, part_v, wk_dir = wk_dir):
    K = ar_model.K
    
    X= ut_hierarch.indicator(cond_v)
    
    em_model = em.MixVMF(K=K,P=atlas.P, X=X,part_vec=part_v)
    
    em_model.V = pt.load(f"{wk_dir}/V_cerebcortex_MDTB(ses1).pt")
    
    em_model.set_param_list(['kappa'])
    
    M = fm.FullMultiModel(ar_model, [em_model])
    
    M.initialize([data])
    
    M, _, _, _ = M.fit_em(iter=200, tol=0.01,
                          fit_arrangement=False,fit_emission=True,first_evidence=False)
    
    Uhat_indiv_data = []
    Uhat_indiv_group = []

    runs = np.unique(part_v)

    for r in runs:

        ind = part_v <=r
        M.emissions[0].X = pt.tensor(ut_hierarch.indicator(cond_v[ind]), dtype=pt.get_default_dtype())
        M.emissions[0].part_vec = pt.tensor(part_v[ind], dtype=pt.get_default_dtype())
        M.emissions[0].initialize(data[:,ind,:])

        emloglik = M.emissions[0].Estep()
        U_indiv_data_run = pt.softmax(emloglik, dim=1)
        Uhat_indiv_data.append(U_indiv_data_run)
        
        U_indiv_group_run, _ = M.arrange.Estep(emloglik)
        Uhat_indiv_group.append(U_indiv_group_run)

    return Uhat_indiv_data, Uhat_indiv_group 

def evaluate_dcbc(U_indiv_data,U_indiv_group,U_group,atlas='MNISymThalamus1',max_dist=40):
    
    atlas_obj, _ = am.get_atlas(atlas)
    coords = atlas_obj.world.T

    wta_group = np.argmax(U_group, axis=0) + 1
    wta_indiv_data = [pt.argmax(r, dim=1) + 1 for r in U_indiv_data]
    wta_indiv_group = [pt.argmax(r, dim=1) + 1 for r in U_indiv_group]

    dist = ut.compute_dist(coords,1.5, backend='numpy')

    T = pd.DataFrame()

    subj_info = pd.read_csv(f'{data_dir}/participants.tsv',sep='\t')
    
    for i, s in enumerate(subj_info['participant_id']):
        
        data = get_data(data_dir = data_dir, atlas_name = atlas, session = 's2')[i].T.numpy()
        
        dcbc_group = compute_DCBC(maxDist=max_dist, binWidth=1.5, parcellation=wta_group,
                                func= data, dist=dist, weighting=True, backend='numpy')

        D1 = {}
        D1['type'] = ['group']
        D1['runs'] = [0]
        D1['dcbc'] = [dcbc_group['DCBC']]
        D1['subject'] = [i + 1]
        T = pd.concat([T, pd.DataFrame(D1)])

        for r in range(len(wta_indiv_data)):
            dcbc_indiv_data = compute_DCBC(maxDist=max_dist, binWidth=1.5, parcellation=wta_indiv_data[r][i],
                                            func= data, dist=dist, weighting=True, backend='numpy')
            dcbc_indiv_group = compute_DCBC(maxDist=max_dist, binWidth=1.5, parcellation=wta_indiv_group[r][i],
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

    ar_model, atlas = get_group_atlas(wk_dir = wk_dir, atlas_name = 'MNISymThalamus1', roi = 'thalamus')

    data = get_data(data_dir = data_dir, atlas_name = 'MNISymThalamus1', session = 's1')

    cond_v, part_v = get_info_emission(data_dir = data_dir, sample_subj = 'sub-02', session = 's1', dataset = 'MDTB')

    Uhat_indiv_data, Uhat_indiv_group = build_indiv_parc_runwise(ar_model, atlas, data, cond_v, part_v, wk_dir = wk_dir)

    D = evaluate_dcbc(Uhat_indiv_data,Uhat_indiv_group,
                      np.load(f'{wk_dir}/Prob_thalamus.npy'),atlas='MNISymThalamus1',max_dist=40)

    print("Check")

    #U_group = np.load(f'{wk_dir}/Prob_thalamus.npy')

    #D = evaluate_dcbc(U_indiv_data,U_indiv_group,U_group,atlas='MNISymThalamus1',max_dist=60)


    








