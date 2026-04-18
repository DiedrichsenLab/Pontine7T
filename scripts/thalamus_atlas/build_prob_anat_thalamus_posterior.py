import numpy as np 
import ants 
import nibabel as nb
import numpy as np
import nibabel as nb
import os 
import subprocess


"This script builds anatomical thalamus probability maps for each subject and then averages them to create a group map. It applies a threshold to the group probability map to create a binary mask of the thalamus."
"It assumes that individual subject thalamus segmentations have already been done using freesurfer, and that the resulting segmentation files are in MNI space."

wk_dir = '/Users/incehusain/fs_projects/'

def convert_mgz_to_nii(input_dir, output_dir=None):
    """
    Convert all .mgz files in a directory to .nii.gz
    Before running, in terminal, write: 
    export FREESURFER_HOME=/Applications/freesurfer_7.4.1/7.4.1
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

        cmd = ["/Applications/freesurfer_7.4.1/7.4.1/bin/mri_convert", in_path, out_path]
        subprocess.run(cmd, check=True)

    print("Done.")

def get_subj_prob_maps(subj =['sub-01'], dataset='Social'):
    
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

if __name__ == '__main__':
    
    #data_dir = f"{wk_dir}/Social"
    
    #subjects = sorted(
    #[d for d in os.listdir(data_dir)
    #if d.startswith("sub-") and os.path.isdir(os.path.join(data_dir, d))])

    subjects = ['sub-10','sub-11','sub-12','sub-13','sub-14','sub-15','sub-16','sub-17','sub-18','sub-19','sub-20',
                'sub-21','sub-22','sub-23','sub-24','sub-25','sub-26','sub-27']

    for sub in subjects:
        convert_mgz_to_nii(input_dir=f"{wk_dir}/Social/{sub}/mri/", output_dir=f"{wk_dir}/Social_pseg/{sub}/")
        get_subj_prob_maps(subj=[sub], dataset='Social')
    
    #convert_mgz_to_nii(input_dir=f"{wk_dir}/Social/sub-03/mri/", output_dir=f"{wk_dir}/Social_pseg/sub-03/")
    #get_subj_prob_maps(subj=['sub-03'],dataset='Social')


