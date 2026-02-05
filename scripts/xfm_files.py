import SUITPy.normalization as suit_norm
import ants

wk_dir = '/Users/incehusain/fs_projects'

def adjust_masks(source_file, mask_file, sub):
    
    source_img = ants.image_read(source_file)
    mask_img = ants.image_read(mask_file)
        
    mask_new = ants.resample_image_to_target(mask_img, source_img, interp_type='nearestNeighbor')

    mask_new.shape == source_img.shape
    mask_new.spacing == source_img.spacing
    mask_new.origin == source_img.origin
    mask_new.direction == source_img.direction
        
    out_file = f"{wk_dir}/xfm_files/Social/{sub}_desc-cereb_mask_adjusted.nii"
    mask_new.to_filename(out_file)
    
    return out_file

def make_xfm_files(source_file, mask_file, sub):
        
    suit_norm.normalize(source_file, mask_file, space='MNISym')
        
    template_file = f"{wk_dir}/tpl-MNI152NLin2009cSym_res-1_T1w.nii"
    displacement_file = f"{wk_dir}/xfm_files/Social/{sub}_T1w_to-MNISym_mode-image_xfm.nii.gz"
    deformation_file = f"{wk_dir}/xfm_files/Social/{sub}_space-MNI152NLin2009cSym_xfm.nii.gz"
        
    suit_norm.deformation_from_displacement(template_file, displacement_file, deformation_file)

if __name__ == '__main__':
    
    sub_list= ['sub-06', 'sub-07']

    for sub in sub_list:
        source_file = f"{wk_dir}/xfm_files/Social/{sub}_T1w.nii"
        mask_file = f"{wk_dir}/xfm_files/Social/{sub}_desc-cereb_mask.nii"

        mask_adjusted = adjust_masks(source_file, mask_file, sub)
        make_xfm = make_xfm_files(source_file, mask_adjusted, sub)



