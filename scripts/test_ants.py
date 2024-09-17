import nibabel as nb
import numpy as np
import ants 

base_dir = '/Volumes/diedrichsen_data$/data/Cerebellum/Pontine7T' 
wkdir = ''





if __name__ == '__main__':
    # Load data
    subj = 'S17'
    filename = f'{base_dir}/imaging_data/{subj}/{subj}_mean_bold.nii'
    img = nb.load(filename)
    data = img.get_fdata()
    print(data.shape)

    # Convert to ants image
    ants_img = ants.image_read(filename)
    print(ants_img)
    ants.plot(ants_img)
    # Save to nifti
    pass 