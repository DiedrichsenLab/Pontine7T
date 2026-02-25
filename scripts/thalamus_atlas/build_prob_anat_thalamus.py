import numpy as np 
import ants 
import nibabel as nb
import numpy as np
import nibabel as nb

"This script builds anatomical thalamus probability maps for each subject and then averages them to create a group map. It applies a threshold to the group probability map to create a binary mask of the thalamus."
"It assumes that individual subject thalamus segmentations have already been done using freesurfer, and that the resulting segmentation files are in MNI space."

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

def thresh_mask(thalamus_mask = '/Users/incehusain/fs_projects/thalamus_masks/thalamus_mask.nii', threshold=0.004):
    img = nb.load(thalamus_mask)
    data = img.get_fdata()
    
    binary_mask = (data >= threshold).astype(np.uint8)
    
    binary_mask_img = nb.Nifti1Image(binary_mask, img.affine, img.header)
    
    nb.save(binary_mask_img, '/Users/incehusain/fs_projects/thalamus_mask_thresholded_025.nii')
