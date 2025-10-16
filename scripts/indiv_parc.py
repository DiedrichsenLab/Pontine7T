import numpy as np
import torch as pt
import nibabel as nb
import nitools as nt
import pandas as pd
import glob 
import matplotlib.pyplot as plt
import Functional_Fusion.atlas_map as am
import Functional_Fusion.dataset as ds
import HierarchBayesParcel.arrangements as ar
import HierarchBayesParcel.emissions as em
import HierarchBayesParcel.full_model as fm
import HierarchBayesParcel.util as ut
import SUITPy as suit
import Functional_Fusion.plot as plot
import os 

wk_dir = '/Volumes/diedrichsen_data$/data/Cerebellum/Pontine7T/atlases/thalamus'
data_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion_new/Pontine7T'

def get_group_atlas(wk_dir = wk_dir, atlas_name = 'MNISymThalamus1', roi = 'thalamus'):
    atlas, _ = am.get_atlas(atlas_name)
    atlas_group = np.load(f"{wk_dir}/Prob_{roi}.npy")
    ar_model = ar.build_arrangement_model(atlas_group, prior_type='prob', atlas=atlas)

    return ar_model, atlas 

def get_data(data_dir = data_dir, atlas_name = 'MNISymThalamus1', session = 'localizerfm'):
    subj_info = pd.read_csv(f'{data_dir}/participants.tsv',sep='\t')
    data = []

    for i, s in enumerate(subj_info['participant_id']):

        file_name = f'{data_dir}/derivatives/ffextract/{s}/{s}_space-{atlas_name.name}_ses-{session}_CondRun.dscalar.nii'

        #MDTB - {s}_space-{atlas_name.name}_ses-s1_CondRun.dscalar.nii (or ses-s2)
        #Language - {s}_space-{atlas_name.name}_ses-localizerfm_CondRun.dscalar.nii
        #MDTB-high-res - {s}_space-{atlas_name.name}_ses-s1_CondRun.dscalar.nii
        
        if not os.path.exists(file_name):
            print(f"File {file_name} does not exist. Skipping subject {s}.")
            continue

        datafile = nb.load(file_name)
        data.append(datafile.get_fdata())

    data = np.stack(data)

    # data is: numsubj x numcond (repetitions x conditions) x numvoxel tensor

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

def build_indiv_parc(ar_model, atlas, data, cond_v, part_v, wk_dir = wk_dir):
    K = ar_model.K
    
    X= ut.indicator(cond_v)
    
    em_model = em.MixVMF(K=K,P=atlas.P, X=X,part_vec=part_v)
    
    em_model.V = pt.load(f"{wk_dir}/V_cerebcortex_Pontine.pt")
    
    em_model.set_param_list(['kappa'])
    
    M = fm.FullMultiModel(ar_model, [em_model])
    
    M.initialize([data])
    
    M, _, _, _ = M.fit_em(iter=200, tol=0.01,
                          fit_arrangement=False,fit_emission=True,first_evidence=False)

    emloglik = M.emissions[0].Estep()

    U_indiv_data = pt.softmax(emloglik, dim=1) 

    U_indiv_group, _ = M.arrange.Estep(emloglik)

    return U_indiv_data, U_indiv_group 

if __name__ == '__main__':

    ar_model, atlas = get_group_atlas()
    
    data = get_data()
    
    cond_v, part_v = get_info_emission()
    
    U_indiv = build_indiv_parc(ar_model, atlas, data, cond_v, part_v)
    
