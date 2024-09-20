import nibabel as nb
import numpy as np
import ants 
import pandas as pd
import matplotlib.pyplot as plt
import glob
import os

base_dir = '/Volumes/diedrichsen_data$/data/Cerebellum/Pontine7T' 

pinfo = pd.read_csv(f'{base_dir}/participants.tsv', sep='\t') 


def align_subjs(src_sub,trg_sub):
    """uses Syn registration to morph the src subject into the trg subject"""
    src = f'{base_dir}/imaging_data/{src_sub}/{src_sub}_mean_bold.nii'
    trg = f'{base_dir}/imaging_data/{trg_sub}/{trg_sub}_mean_bold.nii'
    # Prefix for the resulting outfile
    prefix = f'{base_dir}/bold_normalization/transforms/xfm_{src_sub}_{trg_sub}_'

    src_img = ants.image_read(src)
    trg_img = ants.image_read(trg)
    mytx = ants.registration(fixed=trg_img, moving=src_img, type_of_transform='SyN',outprefix=prefix,write_composite_transform=False)
    wsrc_img = mytx['warpedmovout']
    return src_img, trg_img, wsrc_img

def avrg_images(pattern,outname):
    image_list = glob.glob(pattern)
    data = []
    for img in image_list:
        img = nb.load(img)
        data.append(img.get_fdata())
    data = np.array(data)
    avrg_data = np.mean(data,axis=0)
    avrg_img = nb.Nifti1Image(avrg_data, img.affine)
    nb.save(avrg_img,outname)
    return avrg_img

def show_plot():
    ants.plot(tmp_img,overlay= anat_img*mask_img,axis=1,slices=[40])

def test_plt():
    A = np.random.normal(0,1,(30,30))
    plt.imshow(A,cmap='gray')
    plt.show()
    pass 

def normalize_bold_all(name,template,mask=None,tot='SyN'): 
    subj = pinfo['participant_id']
    subj = subj[pinfo.good==1]
    def_dir = f'{base_dir}/bold_normalization/deformed/{name}'
    xfm_dir = f'{base_dir}/bold_normalization/transforms/{name}'
    if not os.path.isdir(def_dir):
        os.mkdir(def_dir)
    if not os.path.isdir(xfm_dir):
        os.mkdir(xfm_dir)
    trg = f'{base_dir}/bold_normalization/templates/{template}'
    trg_img = ants.image_read(trg)
    for src_sub in subj:
        src = f'{base_dir}/bold_normalization/individual/{src_sub}_mean_sbref.nii'
        if not os.path.isfile(src):
            src = f'{base_dir}/bold_normalization/individual/{src_sub}_mean_bold.nii'
            print(f'Using bold image for {src_sub}')
        else:
            print(f'Using sbref image for {src_sub}')
        prefix = f'{xfm_dir}/xfm_{src_sub}_'
        src_img = ants.image_read(src)
        mytx = ants.registration(fixed=trg_img , moving=src_img,type_of_transform=tot,outprefix=prefix,write_composite_transform=False)
        fname_w = f'{def_dir}/w{src_sub}.nii'
        wsrc_img = mytx['warpedmovout']
        # ants.image_write(src_img,fname_s)
        ants.image_write(wsrc_img,fname_w)
    pattern  = f'{def_dir}/w*.nii'
    outname = f'{def_dir}/wavrg.nii'
    avrg_images(pattern,outname) # Load data

def normalize_suit_all():
    template  = f'{base_dir}/bold_normalization/templates/tpl-SUIT_T1w.nii'
    trg_img = ants.image_read(template)
    subj = pinfo['participant_id']
    subj = subj[pinfo.good==1]
    for src_sub in ["S17"]:
        print(f'Processing {src_sub}')
        src = f'{base_dir}/suit/anatomicals/{src_sub}/c_{src_sub}_T1w.nii'
        msk = f'{base_dir}/suit/anatomicals/{src_sub}/c_{src_sub}_T1w_pcereb_corr.nii'
        prefix = f'{base_dir}/bold_normalization/transforms/xfm_{src_sub}_SUIT_'
        src_img = ants.image_read(src)
        mask_img = ants.image_read(msk)
        mytx = ants.registration(fixed=trg_img , moving=src_img*mask_img,type_of_transform='SyN',outprefix=prefix,write_composite_transform=False)
        fname_w = f'{base_dir}/bold_normalization/deformed/wsuit_{src_sub}.nii'
        wsrc_img = mytx['warpedmovout']
        # ants.image_write(src_img,fname_s)
        ants.image_write(wsrc_img,fname_w)

def make_initial_template(src_sub='S03'): 
    template  = f'{base_dir}/bold_normalization/templates/tpl-MNI152NLin2009cSym_res-1_T2w.nii.gz'
    src = f'{base_dir}/anatomicals/{src_sub}/{src_sub}_T1w.nii'

    prefix = f'{base_dir}/bold_normalization/transforms/MNI2009c_T1w/xfm_{src_sub}_'
    trg_img = ants.image_read(template)
    src_img = ants.image_read(src) 
    mytx = ants.registration(fixed=trg_img , moving=src_img,type_of_transform='SyN',outprefix=prefix,write_composite_transform=True)
    fname_w = f'{base_dir}/bold_normalization/deformed/MNI2009c_T1w/w{src_sub}_T1w.nii'
    wsrc_img = mytx['warpedmovout']
    ants.image_write(wsrc_img,fname_w)
    
    # Read the transform and apply it to the bold image
    tx=ants.read_transform(mytx['fwdtransforms'])
    bold_img=ants.image_read(f'{base_dir}/bold_normalization/individual/{src_sub}_mean_sbref.nii')
    ntemp_img  =tx.apply_to_image(bold_img,reference=trg_img)
    ntemp_name = f'{base_dir}/bold_normalization/templates/tpl-{src_sub}_bold.nii'
    ants.image_write(ntemp_img,ntemp_name)

    pass

def symmetrize_image(src,outname):
    img = nb.load(src)
    data = img.get_fdata()
    data = (data + np.flip(data,axis=0))/2
    img = nb.Nifti1Image(data, img.affine)
    nb.save(img,outname)
    return img

def symmetrize_template(src_sub='S03'):
    inf = f'{base_dir}/bold_normalization/templates/tpl-SyN_bold.nii'
    outf = f'{base_dir}/bold_normalization/templates/tpl-SyNSym_bold.nii'
    symmetrize_image(inf,outf)

if __name__ == '__main__':
    # src_sub='S16'
    # make_initial_template()
    symmetrize_template()
    normalize_bold_all('S03Sym_CC','tpl-S03Sym_bold.nii',tot='SyNCC')
    normalize_bold_all('S03Sym_Agg','tpl-S03Sym_bold.nii',tot='SyNAggro')
    # Save to nifti
    pass 