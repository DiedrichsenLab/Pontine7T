import numpy as np 
import ants 
import nibabel as nb
import os
import Functional_Fusion.dataset as ds
import numpy as np
import nibabel as nb
import Functional_Fusion.atlas_map as am
import torch as pt
import HierarchBayesParcel.util as ut
import pandas as pd 
import HierarchBayesParcel.emissions as em
import HierarchBayesParcel.arrangements as ar
import HierarchBayesParcel.full_model as fm
import matplotlib.pyplot as plt
import pandas as pd 
import recenter_data as rd 

wk_dir = '/Users/incehusain/fs_projects'
base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion_new' 

def get_subj_prob_maps(subj =['sub-01'], dataset='Language7T'):
    
    #Get individual subject anatomical thalamus probability maps in MNI space; assumes freesurfer thalamus segmentations already done
    
    template = f"{wk_dir}/tpl-MNI152NLin2009cSym_res-1_T1w.nii"
    template_img = ants.image_read(template)

    lut_file = f"{wk_dir}/compressionLookupTable.txt"

    lut_label_to_name = {}
    with open(lut_file, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                parts = line.strip().split()
                label = int(parts[0])
                name = parts[1]
                lut_label_to_name[label] = name
                
    thalamus_labels = sorted(list(lut_label_to_name.keys()))
    
    K = len(thalamus_labels)    
    
    for sub in subj:
        segmentation_path = f"{wk_dir}/dseg_segmentations_{dataset}/{sub}_thalamus_dseg.nii.gz"
        seg_img = ants.image_read(segmentation_path)
        T1_nii = f"{wk_dir}/freesurfer_{dataset}/{sub}/mri/T1.mgz"  
        T1 = ants.image_read(T1_nii)
        
        output_prefix = f"{wk_dir}/{sub}_to_MNI_"
        registration = ants.registration(fixed=template_img, moving=T1, type_of_transform='SyN', outprefix=output_prefix)

        warped_prob_maps_list = []
        for i, label in enumerate(thalamus_labels):
            print(f"Processing label {label} ({i+1}/{K})")
            mask_3d = (seg_img == label).astype('float32')
            warped_mask_3d = ants.apply_transforms(fixed=template_img, moving=mask_3d, transformlist=registration['fwdtransforms'], interpolator='linear')
            warped_prob_maps_list.append(warped_mask_3d)
        warped_prob_4d_img = ants.merge_channels(warped_prob_maps_list)
        
        save_data_5d = f"{wk_dir}/{sub}_5d_thalamus_prob_map.nii.gz" #fsleyes adds an extra dimension because it is assuming time series data, so we need to convert back to 4d
        ants.image_write(warped_prob_4d_img, save_data_5d)
        img_5d = nb.load(save_data_5d)
        data_5d = img_5d.get_fdata()    
        data_4d = np.squeeze(data_5d)
        
        img_4d = nb.Nifti1Image(data_4d, img_5d.affine)
        data_4d_name = f"{wk_dir}/{sub}_4d_thalamus_prob_map.nii.gz"        

        nb.save(img_4d, data_4d_name)


def get_group_prob_map(subject_file ='sub-01_4d_thalamus_prob_map.nii.gz'):

    #averages subject probability maps to get group map

    N= len(subject_file)
    first_img = nb.load(f"{wk_dir}/{subject_file[0]}")
    sum_of_data = first_img.get_fdata().astype('float64')

    for i in range(1, N):
        img = nb.load(f"{wk_dir}/{subject_file[i]}")
        sum_of_data += img.get_fdata()

    mean_data = sum_of_data / N

    mean_img = nb.Nifti1Image(mean_data.astype('float32'), first_img.affine, first_img.header)

    nb.save(mean_img, f"{wk_dir}/group_mean_thalamus_prob_map.nii.gz")


def get_Vs(subj =['sub-01'], dataset = 'Language', session = 'ses-localizerfm', map_file = f"{wk_dir}/pseg_Language7T/sub-01_4d_thalamus_prob_map.nii.gz", map_type='indiv'):

    #weighted average of data times the probability maps to get V matrices of size Tasks x Parcels 

    for sub in subj:
        data, info, ds_obj = ds.get_dataset(base_dir, dataset ,atlas='MNISymThalamus1', 
                                            sess=session, 
                                            subj=[sub], 
                                            type='CondAll')
        
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

def get_Vs_from_data(data, info, map_file, map_type='indiv'):
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

def build_emissions(K,P,atlas='MNISymThalamus1'):

    data, info,ds_obj = ds.get_dataset(base_dir,'MDTB',atlas=atlas,sess='ses-s1',  subj=None, type='CondAll')
    #cond_v = info['cond_num_uni']
    #part_v = info['run']
    #X= ut.indicator(cond_v)
    em_model_mdtb1 = em.MixVMF(K=K,P=P)
    em_model_mdtb1.initialize(data)

    data2, info2,ds_obj2 = ds.get_dataset(base_dir,'MDTB',atlas=atlas,sess='ses-s2',subj=None, type='CondAll')
    #cond_v2 = info2['cond_num_uni']
    #part_v2 = info2['run']
    #X2= ut.indicator(cond_v2)
    em_model_mdtb2 = em.MixVMF(K=K,P=P)
    em_model_mdtb2.initialize(data2)

    data3, info3, ds_obj3 = ds.get_dataset(base_dir,'Language',atlas=atlas, sess='ses-localizerfm', subj=None, type='CondRun')
    #cond_v3 = info3['task']
    #part_v3 = info3['run']
    #X3= ut.indicator(cond_v3)
    em_model_language = em.MixVMF(K=K,P=P)
    em_model_language.initialize(data3)

    #data4, info4, ds_obj4 = ds.get_dataset(base_dir,'Pontine',atlas=atlas, sess='ses-s1', subj=None, type='CondRun')    
    #cond_v4 = info4['task']
    #part_v4 = info4['run']
    #X4= ut.indicator(cond_v4)
    #em_model_pontine = em.MixVMF(K=K,P=P, X=X4,part_vec=part_v4)
    #em_model_pontine(data4)

    return em_model_mdtb1, em_model_mdtb2, em_model_language

def estimate_new_atlas(): 
     
    # Get atlas for dentate 

     atlas, _ = am.get_atlas('MNISymThalamus1')    

     ar_model = ar.ArrangeIndependent(K=32, P=atlas.P, spatial_specific=True, remove_redundancy=False)

     emissions = build_emissions(ar_model.K,atlas.P,atlas='MNISymThalamus1')

     emissions[0].V = pt.load(f"{wk_dir}/V_cerebcortex_MDTB(ses1).pt")
     emissions[1].V = pt.load(f"{wk_dir}/V_cerebcortex_MDTB(ses2).pt")
     emissions[2].V = pt.load(f"{wk_dir}/V_cerebcortex_Language.pt")
     #emissions[3].V = pt.load(f"{wk_dir}/V_cerebcortex_Pontine.pt")

     for em_model in emissions:
            em_model.set_param_list(['kappa'])
    
     M= fm.FullMultiModel(ar_model, [emissions[0], emissions[1], emissions[2]])

     M.initialize()

     M, ll, _, _ = M.fit_em(iter=200, tol=0.01,
        fit_arrangement=True,fit_emission=True,first_evidence=True)
     
     Prob = M.arrange.marginal_prob().numpy()
     np.save(f"{wk_dir}/Prob_thalamus.npy",Prob)

     return M

def calc_cosine_similarity(subj = ['sub-01', 'sub-02'], type='group', dataset='Language', sess = 'ses-localizerfm'):

    #computes average cosine similarities between subjects for each same parcel 

    avg_similarity_per_parcel = []
    similarity_distributions_per_parcel = []
    
    for k in range(58):
        parcel_vectors_subj = []

        for sub in subj:
            V_indiv = pt.load(f"{wk_dir}/V_matrices_{dataset}/{sess}_V_{type}_{sub}_norm.pt")
            parcel_vector = V_indiv[:, k]
            parcel_vectors_subj.append(parcel_vector)
            
            X = np.stack(parcel_vectors_subj, axis=0)

        cosine_similarity_matrix = X@X.T

        triangle_indices = np.triu_indices_from(cosine_similarity_matrix, k=1)
        cosine_similarity_matrix_values = cosine_similarity_matrix[triangle_indices]

        similarity_distributions_per_parcel.append(cosine_similarity_matrix_values)

        mean_similarity = np.mean(cosine_similarity_matrix_values) 

        avg_similarity_per_parcel.append(mean_similarity) 
    
    return avg_similarity_per_parcel, similarity_distributions_per_parcel   

def calc_cosine_similarity_within(subj = ['sub-01', 'sub-02'], type='group', dataset='Language', sess = 'ses-localizerfm'):

#calculates cosine similarities within each subject between all parcel pairs
    
    avg_similarity_per_subject = []
    subj_similarity_matrix = {}
    
    for sub in subj:
        V = pt.load(f"{wk_dir}/V_matrices_{dataset}/{sess}_V_{type}_{sub}_norm.pt")
        X = V.T
        
        cosine_similarity_matrix = X@X.T

        triangle_indices = np.triu_indices_from(cosine_similarity_matrix, k=1)
        cosine_similarity_matrix_values = cosine_similarity_matrix[triangle_indices]

        mean_similarity = np.mean(cosine_similarity_matrix_values) 

        avg_similarity_per_subject.append(mean_similarity) 
        subj_similarity_matrix[sub] = cosine_similarity_matrix

    return avg_similarity_per_subject, subj_similarity_matrix

def avg_similarity_matrices(subj_similarity_matrices):

    all_matrices = []
    
    for sub in subj_similarity_matrices.keys():
        all_matrices.append(subj_similarity_matrices[sub])
    
    mean_similarity_matrix = np.mean(np.array(all_matrices), axis=0)
    
    return mean_similarity_matrix

if __name__ == '__main__':
    
    #1: Get individual subject probability maps in MNI space
    #get_subj_prob_maps(subj=subjects, dataset='Language7T')
    
    #2: Get group average probability map
    #subject_files = [f"{sub}_4d_thalamus_prob_map.nii.gz" for sub in subjects]
    #get_group_prob_map(subject_file=subject_files)
    
    # Step 3: Get V matrices for individual and group maps
    #for sub in subjects:

    dataset= 'WMFS'
    sessions = ['ses-01', 'ses-02']

    base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion_new' 

    #sub_list = ['sub-01','sub-02','sub-03', 'sub-04', 'sub-05','sub-06']
    
    sub_list = ['sub-01','sub-02','sub-03', 'sub-04', 'sub-05','sub-06','sub-07','sub-08','sub-09', 'sub-10','sub-11','sub-12','sub-13','sub-14','sub-15','sub-16']

    #sub_list = ['sub-02','sub-03','sub-04', 'sub-06','sub-08','sub-09', 'sub-10', 'sub-12','sub-14','sub-15','sub-17', 'sub-18',
     #         'sub-19','sub-20', 'sub-21', 'sub-22','sub-24', 'sub-25','sub-26','sub-27', 'sub-28', 'sub-29','sub-30', 'sub-31']
        
    data_ses1, info_ses1, ds_obj_ses1 = ds.get_dataset(base_dir, dataset ,atlas='MNISymThalamus1', 
                                            sess=sessions[0], 
                                            subj=sub_list, 
                                            type='CondAll')
    
    data_ses2, info_ses2, ds_obj_ses2 = ds.get_dataset(base_dir, dataset ,atlas='MNISymThalamus1', 
                                            sess=sessions[1], 
                                            subj=sub_list, 
                                            type='CondAll')
    
    info_ses1_path = f'{base_dir}/{dataset}/derivatives/ffextract/sub-02/sub-02_{sessions[0]}_CondAll.tsv'
    info_ses1 = pd.read_csv(info_ses1_path, sep='\t')

    info_ses2_path = f'{base_dir}/{dataset}/derivatives/ffextract/sub-02/sub-02_{sessions[1]}_CondAll.tsv'
    info_ses2 = pd.read_csv(info_ses2_path, sep='\t')   

    #recentered_ses1 = rd.recenter_data(data_ses1, info_ses1, center_full_code='rest_task', keep_center=True)

    #recentered_ses2 = rd.recenter_data(data_ses2, info_ses2, center_full_code='rest_task', keep_center=True)

    #merged_data, merged_info = rd.merge_sessions([recentered_ses1[0], recentered_ses2[0]], [recentered_ses1[1], recentered_ses2[1]])

    merged_data, merged_info = rd.merge_sessions([data_ses1, data_ses2], [info_ses1, info_ses2])   
    
    baseline_all = np.ones(merged_data.shape[1])
    mean_centered = ds.remove_baseline(merged_data, baseline_all)

    map_file=f"{wk_dir}/group_mean_thalamus_prob_map.nii.gz"
    Vs_merged = get_Vs_from_data(mean_centered, merged_info, map_file, map_type='merged')
    
    for subj_idx, sub in enumerate(sub_list):
        pt.save(Vs_merged[subj_idx], f"{wk_dir}/V_matrices_{dataset}/V_concat_{sub}_norm.pt")


    #-------------------------------------#
    
    dataset = 'IBC'
    sess = 'ses-tom'

    sub_list = ['sub-01', 'sub-04', 'sub-05','sub-06', 'sub-07', 'sub-08', 'sub-09', 'sub-11', 'sub-12', 'sub-13', 'sub-14', 'sub-15']

    for sub in sub_list:        
        cosine = calc_cosine_similarity(subj=sub_list, type='group', dataset='IBC', sess='ses-tom')
        pt.save(cosine[0],f"{wk_dir}/cosine_similarities_{dataset}/{sess}_cosine_btwn_subj_avg_parcels_group.pt")
        pt.save(cosine[1],f"{wk_dir}/cosine_similarities_{dataset}/{sess}_cosine_btwn_subj_parcels_group.pt")

        subj_matrices = calc_cosine_similarity_within(subj=sub_list, type='group', dataset='IBC', sess='ses-tom')
        pt.save(subj_matrices[0],f"{wk_dir}/cosine_similarities_{dataset}/{sess}_cosine_within_subj_groupV.pt")

        mean_similarity_matrix = avg_similarity_matrices(subj_matrices[1])
        pt.save(mean_similarity_matrix,f"{wk_dir}/cosine_similarities_{dataset}/{sess}_cosine_avg_all_parcels.pt")

        #get_Vs(subj=sub_list, dataset = 'Social', session = 'ses-social', map_file = f"{wk_dir}/group_mean_thalamus_prob_map.nii.gz", map_type='group')
        #get_Vs(subj=sub_list, dataset = 'Language', session = 'ses-localizerfm', map_file=f"{wk_dir}/group_mean_thalamus_prob_map.nii.gz", map_type='group')
        #emissions = build_emissions(K=58,P=24291,atlas='MNISymThalamus1')

    
    #MDTB: sub_list = ['sub-02','sub-03','sub-04', 'sub-06','sub-08','sub-09', 'sub-10', 'sub-12','sub-14','sub-15','sub-17', 'sub-18',
             # 'sub-19','sub-20', 'sub-21', 'sub-22','sub-24', 'sub-25','sub-26','sub-27', 'sub-28', 'sub-29','sub-30', 'sub-31']
        
    #Language:  sub_list=['sub-01', 'sub-02', 'sub-04', 'sub-06', 'sub-07', 'sub-08', 'sub-09',
               #'sub-10', 'sub-12', 'sub-13', 'sub-14', 'sub-15', 'sub-16', 'sub-17', 'sub-18', 'sub-19']
        
    #Social:  sub_list= ['sub-04', 'sub-05',  'sub-07', 'sub-08', 'sub-09', 'sub-10',
           #    'sub-11', 'sub-12', 'sub-13', 'sub-14',  
            #   'sub-16', 'sub-17', 'sub-18', 'sub-19', 'sub-20', 'sub-21', 'sub-22', 'sub-24',
             #  'sub-25', 'sub-26', 'sub-27']
    
    #WMFS: sub_list = ['sub-01','sub-02','sub-03', 'sub-04', 'sub-05','sub-06','sub-07','sub-08','sub-09', 'sub-10','sub-11','sub-12','sub-13','sub-14','sub-15','sub-16']

    #Nishimoto: sub_list = ['sub-01','sub-02','sub-03', 'sub-04', 'sub-05','sub-06']

    #IBC: sub_list = ['sub-01','sub-02', 'sub-04', 'sub-05','sub-06', 'sub-07', 'sub-08', 'sub-09', 'sub-11', 'sub-12', 'sub-13', 'sub-14', 'sub-15']
    ## IBS sessions: 'ses-archi', 'ses-clips4',  'ses-enumeration','ses-hcp1', 'ses-hcp2','ses-lyon1', 'ses-lyon2','ses-mathlang', 'ses-mtt1', 'ses-mtt2',  'ses-preference', 'ses-rsvplanguage', 'ses-spatialnavigation', 'ses-tom']




