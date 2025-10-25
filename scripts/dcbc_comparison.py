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
<<<<<<< HEAD
base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion_new'
=======
ff_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion_new'
>>>>>>> 0ec1f61c8d66c1ada1feda1513c8f38f491d906f

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

def build_indiv_parc_runwise_globalk(ar_model, atlas, data, cond_v, part_v, wk_dir = wk_dir):
    #global kappa 
    
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

def build_indiv_parc_runwise_runk(ar_model, atlas, data, cond_v, part_v, wk_dir = wk_dir):

    #runwise kappas

    K = ar_model.K

    V_fixed = pt.load(f"{wk_dir}/V_cerebcortex_MDTB(ses1).pt")

    Uhat_indiv_data_cumulative = []
    Uhat_indiv_group_cumulative = []

    runs = np.unique(part_v)

    for r in runs:
        print(f"Processing run data up to run {r}...")

        ind = part_v <=r

        cumulative_data = data[:,ind,:]    
        cumulative_cond_v = cond_v[ind]
        cumulative_part_v = part_v[ind]   

        X_cumulative = ut_hierarch.indicator(cumulative_cond_v)
        em_model_cumulative = em.MixVMF(K=K,P=atlas.P, X=X_cumulative,part_vec=cumulative_part_v)

        em_model_cumulative.V = V_fixed
        em_model_cumulative.set_param_list(['kappa'])

        M_cumulative = fm.FullMultiModel(ar_model, [em_model_cumulative])

        M_cumulative.initialize([cumulative_data])

        M_cumulative, _, _, _ = M_cumulative.fit_em(iter=200, tol=0.01,
                              fit_arrangement=False,fit_emission=True,first_evidence=False)


        emloglik = M_cumulative.emissions[0].Estep()
        U_indiv_data_run = pt.softmax(emloglik, dim=1)
        Uhat_indiv_data_cumulative.append(U_indiv_data_run)
        
        U_indiv_group_run, _ = M_cumulative.arrange.Estep(emloglik)
        Uhat_indiv_group_cumulative.append(U_indiv_group_run)

    return Uhat_indiv_data_cumulative, Uhat_indiv_group_cumulative 

def evaluate_dcbc(U_indiv_data,U_indiv_group,U_group,atlas='MNISymThalamus1',max_dist=40):
    
    atlas_obj, _ = am.get_atlas(atlas)
    coords = atlas_obj.world.T

    wta_group = np.argmax(U_group, axis=0) + 1
    wta_indiv_data = [pt.argmax(r, dim=1) + 1 for r in U_indiv_data]
    wta_indiv_group = [pt.argmax(r, dim=1) + 1 for r in U_indiv_group]

    dist = ut.compute_dist(coords, backend='numpy')

    T = pd.DataFrame()

    subj_info = pd.read_csv(f'{data_dir}/participants.tsv',sep='\t')
    
    for i, s in enumerate(subj_info['participant_id']):

        data, info, ds_obj = ds.get_dataset(base_dir,'MDTB',atlas="MNISymThalamus1",sess='ses-s2', subj=None, 
                                type='CondRun')
        
        datai = data[i].T
        
        #data = get_data(data_dir = data_dir, atlas_name = atlas, session = 's2')[i].T.numpy()
        
        dcbc_group = compute_DCBC(maxDist=max_dist, binWidth=1.5, parcellation=wta_group,
                                func= datai, dist=dist, weighting=True, backend='numpy')

        D1 = {}
        D1['type'] = ['group']
        D1['runs'] = [0]
        D1['dcbc'] = [dcbc_group['DCBC']]
        D1['subject'] = [i + 1]
        T = pd.concat([T, pd.DataFrame(D1)])

        for r in range(len(wta_indiv_data)):

            print(f"Processing run data up to run {r}...")


            dcbc_indiv_data = compute_DCBC(maxDist=max_dist, binWidth=1.5, parcellation=wta_indiv_data[r][i],
                                            func= datai, dist=dist, weighting=True, backend='numpy')
            dcbc_indiv_group = compute_DCBC(maxDist=max_dist, binWidth=1.5, parcellation=wta_indiv_group[r][i],
                                            func= datai, dist=dist, weighting=True, backend='numpy')
            
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

<<<<<<< HEAD
    data, info, ds_obj = ds.get_dataset(base_dir,'MDTB',atlas="MNISymThalamus1",sess='ses-s1', subj=None, 
                                type='CondRun')

    #cond_v, part_v = get_info_emission(data_dir = data_dir, sample_subj = 'sub-02', session = 's1', dataset = 'MDTB')

    Uhat_indiv_data, Uhat_indiv_group = build_indiv_parc_runwise_globalk(ar_model, atlas, data, info.cond_num, info.run, wk_dir = wk_dir)
=======
    data,info,ds = ds.get_dataset(ff_dir,'MDTB', atlas = 'MNISymThalamus1', sess = 'ses-s1',type='CondRun')


    Uhat_indiv_data, Uhat_indiv_group = build_indiv_parc_runwise_globalk(ar_model, atlas, data, 
                                                                            cond_v = info.cond_num, 
                                                                            part_v = info.run, 
                                                                            wk_dir = wk_dir)
>>>>>>> 0ec1f61c8d66c1ada1feda1513c8f38f491d906f

    D = evaluate_dcbc(Uhat_indiv_data,Uhat_indiv_group,
                      np.load(f'{wk_dir}/Prob_thalamus.npy'),atlas='MNISymThalamus1',max_dist=40)

    print("Check")

    #U_group = np.load(f'{wk_dir}/Prob_thalamus.npy')

    #D = evaluate_dcbc(U_indiv_data,U_indiv_group,U_group,atlas='MNISymThalamus1',max_dist=60)


    








