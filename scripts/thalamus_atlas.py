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


def get_Vs(subj =['sub-01'], map_file = f"{wk_dir}/pseg_Language7T/sub-01_4d_thalamus_prob_map.nii.gz", map_type='indiv'):

    #weighted average of data times the probability maps to get V matrices of size Tasks x Parcels 

    for sub in subj:
        data, info, ds_obj = ds.get_dataset(base_dir,'Language',atlas='MNISymThalamus1', 
                                            sess='ses-localizerfm', 
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
        
        pt.save(Vs_subj_normalized,f"{wk_dir}/V_matrices/V_{map_type}_{sub}_norm.pt")


def calc_cosine_similarity(subj = ['sub-01', 'sub-02'], type='group'):

    avg_similarity_per_parcel = []
    similarity_distributions_per_parcel = []
    
    for k in range(58):
        parcel_vectors_subj = []

        for sub in subj:
            V_indiv = pt.load(f"{wk_dir}/V_matrices/V_{type}_{sub}_norm.pt")
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

def calc_cosine_similarity_within(subj = ['sub-01', 'sub-02'], type='group'):

    avg_similarity_per_subject = []
    
    for sub in subj:
        V = pt.load(f"{wk_dir}/V_matrices/V_{type}_{sub}_norm.pt")
        X = V.T
        
        cosine_similarity_matrix = X@X.T

        triangle_indices = np.triu_indices_from(cosine_similarity_matrix, k=1)
        cosine_similarity_matrix_values = cosine_similarity_matrix[triangle_indices]

        mean_similarity = np.mean(cosine_similarity_matrix_values) 

        avg_similarity_per_subject.append(mean_similarity) 

    return avg_similarity_per_subject


if __name__ == '__main__':

    #subjects = ['sub-01','sub-02','sub-03','sub-04','sub-05',
                #'sub-06','sub-07','sub-08','sub-09','sub-10']
    
    #1: Get individual subject probability maps in MNI space
    #get_subj_prob_maps(subj=subjects, dataset='Language7T')
    
    #2: Get group average probability map
    #subject_files = [f"{sub}_4d_thalamus_prob_map.nii.gz" for sub in subjects]
    #get_group_prob_map(subject_file=subject_files)
    
    # Step 3: Get V matrices for individual and group maps
    #for sub in subjects:
                
    sub_list=['sub-01','sub-02','sub-03','sub-04',
              'sub-06','sub-07','sub-08','sub-09']

    for sub in sub_list:
        cosine = calc_cosine_similarity_within(subj=sub_list, type='subj')
        get_Vs(subj=sub_list, map_file = f"{wk_dir}/pseg_Language7T/{sub}_4d_thalamus_prob_map.nii.gz", map_type='indiv')
        get_Vs(subj='group', map_file=f"{wk_dir}/group_mean_thalamus_prob_map.nii.gz", map_type='group')

