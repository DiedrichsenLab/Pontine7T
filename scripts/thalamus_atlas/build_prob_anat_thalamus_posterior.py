import numpy as np 
import ants 
import nibabel as nb
import numpy as np
import nibabel as nb
import os 
import subprocess
import glob


"This script builds anatomical thalamus probability maps for each subject and then averages them to create a group map. It applies a threshold to the group probability map to create a binary mask of the thalamus."
"It assumes that individual subject thalamus segmentations have already been done using freesurfer, and that the resulting segmentation files are in MNI space."

wk_dir = '/Users/incehusain/fs_projects/'

def convert_mgz_to_nii(input_dir, output_dir=None):
    """
    Convert all .mgz files in a directory to .nii.gz
    Before running, in terminal, write: 
    export FREESURFER_HOME=/Applications/freesurfer/8.2.0
    source $FREESURFER_HOME/SetUpFreeSurfer.sh 
    """
    
    os.makedirs(output_dir, exist_ok=True)

    mgz_files = [f for f in os.listdir(input_dir) if f.startswith("posterior") and f.endswith(".mgz")]

    print(f"Found {len(mgz_files)} .mgz files")

    for i, fname in enumerate(sorted(mgz_files)):
        in_path = os.path.join(input_dir, fname)
        
        out_name = fname.replace(".mgz", ".nii.gz")
        out_path = os.path.join(output_dir, out_name)

        print(f"[{i+1}/{len(mgz_files)}] Converting {fname}")

        cmd = ["/Applications/freesurfer/8.2.0/bin/mri_convert", in_path, out_path]
        subprocess.run(cmd, check=True)

    print("Done.")

def get_subj_prob_maps(subj =['sub-01'], dataset='Social'):

    wk_dir = '/Users/incehusain/fs_projects/'
    
    #Get individual subject anatomical thalamus probability maps in MNI space; assumes freesurfer thalamus segmentations already done
    
    template = f"{wk_dir}/tpl-MNI152NLin2009cSym_res-1_T1w.nii"
    template_img = ants.image_read(template)
    
    for sub in subj:
        print(f"\nProcessing subject: {sub}")
        
        T1_nii = f"{wk_dir}/{dataset}/{sub}/mri/T1.mgz"  
        T1 = ants.image_read(T1_nii)
        
        output_prefix = f"{wk_dir}/{sub}_to_MNI_"
        registration = ants.registration(fixed=template_img, moving=T1, type_of_transform='SyN', outprefix=output_prefix)

        pseg_dir = f"{wk_dir}/{dataset}_pseg/{sub}/"
        pseg_files = sorted([
            os.path.join(pseg_dir, f)
            for f in os.listdir(pseg_dir)
            if f.endswith("v13.T1.nii.gz") 
        ])

        warped_prob_maps_list = []

        for i, pseg_file in enumerate(pseg_files):
            print(f"Warping {i+1}/{len(pseg_files)}: {os.path.basename(pseg_file)}")
            prob_img = ants.image_read(pseg_file)
            warped_prob = ants.apply_transforms(fixed=template_img, moving=prob_img, transformlist=registration['fwdtransforms'], interpolator='linear')
            warped_prob_maps_list.append(warped_prob)
        warped_prob_4d = ants.merge_channels(warped_prob_maps_list)
        
        tmp_5d = f"{wk_dir}/{dataset}_pseg/{sub}/{sub}_5d_tmp_prob_map.nii.gz" #fsleyes adds an extra dimension because it is assuming time series data, so we need to convert back to 4d
        ants.image_write(warped_prob_4d, tmp_5d)
        img_5d = nb.load(tmp_5d)
        data_5d = img_5d.get_fdata()    
        data_4d = np.squeeze(data_5d)
        
        data_4d_name = f"{wk_dir}/{dataset}_pseg/{sub}/{sub}_4d_thalamus_prob_map.nii.gz" 
        img_4d = nb.Nifti1Image(data_4d, img_5d.affine)

        nb.save(img_4d, data_4d_name)


def get_group_prob_map(wk_dir):

    pattern = os.path.join(wk_dir, "*_pseg", "sub-*", "sub-*_4d_thalamus_prob_map.nii.gz")
    subject_files = glob.glob(pattern)

    #averages subject probability maps to get group map

    N= len(subject_files)
    first_img = nb.load(f"{subject_files[0]}")
    sum_of_data = first_img.get_fdata().astype('float64')

    for i in range(1, N):
        img = nb.load(f"{subject_files[i]}")
        sum_of_data += img.get_fdata()

    mean_data = sum_of_data / N

    mean_img = nb.Nifti1Image(mean_data.astype('float32'), first_img.affine, first_img.header)

    nb.save(mean_img, f"{wk_dir}/group_mean_thalamus_prob_map.nii.gz")

def get_thalamus_mask(wk_dir, threshold=0.25):

    in_path = os.path.join(wk_dir, "67_group_mean_thalamus_prob_map.nii.gz")
    out_path = os.path.join(wk_dir, f"group_thalamus_mask_thr{threshold}.nii.gz")

    img = nb.load(in_path)
    data = img.get_fdata()  # shape: (X, Y, Z, N)

    # collapse across nuclei, determine max value for voxel across nuclei 
    max_data = data.max(axis=3)

    # threshold
    binary_data = (max_data > threshold).astype(np.uint8)

    mask_img = nb.Nifti1Image(binary_data, img.affine, img.header)
    nb.save(mask_img, out_path)

    print(f"Saved 3D thalamus mask to: {out_path}")


import nibabel as nb
import numpy as np
import os

def symmetrize_mask(wk_dir, threshold=25):

    in_path = os.path.join(wk_dir, f"group_thalamus_mask_thr0{threshold}.nii.gz")
    out_path = os.path.join(wk_dir, f"group_thalamus_mask_thr0{threshold}_sym.nii.gz")

    img = nb.load(in_path)
    data = img.get_fdata()

    # flip left-right (assumes axis 0 is L-R)
    flipped = np.flip(data, axis=0)

    # combine original + flipped
    sym_data = np.maximum(data, flipped)

    # ensure binary
    sym_data = (sym_data > 0).astype(np.uint8)

    sym_img = nb.Nifti1Image(sym_data, img.affine, img.header)
    nb.save(sym_img, out_path)

    print(f"Saved symmetric mask to: {out_path}")

def symmetrize_copy_right_to_left(wk_dir):

    in_path = os.path.join(wk_dir, "corr_group_thalamus_mask_thr02_sym_copy.nii.gz")
    out_path = os.path.join(wk_dir, "group_thalamus_mask_thr02_sym_RtoL.nii.gz")

    img = nb.load(in_path)
    data = img.get_fdata()

    # Check orientation
    axcodes = nb.aff2axcodes(img.affine)
    print("Axis orientation:", axcodes)

    # Find left-right axis
    lr_axis = None
    for i, ax in enumerate(axcodes):
        if ax in ['L', 'R']:
            lr_axis = i
            break

    if lr_axis is None:
        raise ValueError("Could not find left-right axis.")

    # Flip across left-right axis
    flipped = np.flip(data, axis=lr_axis)

    # Create output copy
    sym_data = data.copy()

    # Determine midpoint
    mid = data.shape[lr_axis] // 2

    # Replace LEFT half with mirrored RIGHT half
    slicer = [slice(None)] * 3
    slicer[lr_axis] = slice(0, mid)

    sym_data[tuple(slicer)] = flipped[tuple(slicer)]

    # Ensure binary
    sym_data = (sym_data > 0).astype(np.uint8)

    sym_img = nb.Nifti1Image(sym_data, img.affine, img.header)
    nb.save(sym_img, out_path)

    print(f"Saved right→left symmetric mask to: {out_path}")

if __name__ == '__main__':

    right_to_left = symmetrize_copy_right_to_left(wk_dir)

    #sym_thalamus_mask = symmetrize_mask(wk_dir, threshold=2)

    #thalamus_mask = get_thalamus_mask(wk_dir, threshold=0.18)

    #group_map = get_group_prob_map(wk_dir)
    
    #data_dir = f"{wk_dir}/Language"
    
    #subjects = sorted(
    #[d for d in os.listdir(data_dir)
    #if d.startswith("sub-") and os.path.isdir(os.path.join(data_dir, d))])

    #subjects = ['sub-29']
    #dataset = 'MDTB'

    #for sub in subjects:
     #   convert_mgz_to_nii(input_dir=f"{wk_dir}/{dataset}/{sub}/mri/", output_dir=f"{wk_dir}/{dataset}_pseg/{sub}/")
      #  get_subj_prob_maps(subj=[sub], dataset=dataset)
    
    #convert_mgz_to_nii(input_dir=f"{wk_dir}/Social/sub-03/mri/", output_dir=f"{wk_dir}/Social_pseg/sub-03/")
    #get_subj_prob_maps(subj=['sub-03'],dataset='Social')


