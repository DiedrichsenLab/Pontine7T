import numpy as np 
import nibabel as nb
import Functional_Fusion.dataset as ds
import numpy as np
import nibabel as nb
import Functional_Fusion.atlas_map as am
import torch as pt

wk_dir = '/Users/incehusain/fs_projects'
base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion_new' 

def get_Vs(subj =['sub-01'], dataset = 'Language', session = 'ses-localizerfm', map_file = f"{wk_dir}/pseg_Language7T/sub-01_4d_thalamus_prob_map.nii.gz", map_type='indiv'):

    #weighted average of data times the probability maps to get V matrices of size Tasks x Parcels 

    for sub in subj:
        data, info, ds_obj = ds.get_dataset(base_dir, dataset ,atlas='MNISymThalamus1', 
                                            sess=session, 
                                            subj=[sub], 
                                            type='CondHalf')
        
        atlas, _ = am.get_atlas('MNISymThalamus1')
        
        map_img = nb.load(map_file)
        map_data = atlas.read_data(map_img) 
        
        Vs = data[0] @ map_data.T
        
        parcel_weight_sums_subj = map_data.sum(axis=1)
        parcel_weight_sums_subj[parcel_weight_sums_subj == 0] = 1
        
        Vs_avg = Vs/ parcel_weight_sums_subj  
        
        #pt.save(Vs_avg,f"{wk_dir}/V_matrices/V_{map_type}_{subj}_unnorm.pt")
        
        col_norms = np.linalg.norm(Vs_avg, axis=0)
        col_norms[col_norms == 0] = 1
        Vs_subj_normalized = Vs_avg / col_norms
        
        pt.save(Vs_subj_normalized,f"{wk_dir}/V_matrices_{dataset}/{session}_V_{map_type}_{sub}_norm.pt")

def get_Vs_from_data(data, map_file):
    atlas, _ = am.get_atlas('MNISymThalamus1')
    map_img = nb.load(map_file)
    map_data = atlas.read_data(map_img)
    
    Vs_list = []
    
    for subj_idx in range(data.shape[0]):
        subj_data = data[subj_idx]
        Vs = subj_data @ map_data.T
        
        parcel_weight_sums = map_data.sum(axis=1)
        parcel_weight_sums[parcel_weight_sums == 0] = 1
        Vs_avg = Vs / parcel_weight_sums
        
        col_norms = np.linalg.norm(Vs_avg, axis=0)
        col_norms[col_norms == 0] = 1
        Vs_subj_normalized = Vs_avg / col_norms
        
        Vs_list.append(Vs_subj_normalized)
        
    return np.stack(Vs_list, axis=0)