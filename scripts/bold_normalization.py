import nibabel as nb
import numpy as np
import ants 
import pandas as pd
import matplotlib.pyplot as plt
import glob
import os
import nitools as nt

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

def normalize_bold_all(name,template,mask=None,tot='SyN',kwargs={},subj=None): 
    if subj is None: 
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
        mytx = ants.registration(fixed=trg_img , 
            moving=src_img,
            type_of_transform=tot,
            outprefix=prefix,
            write_composite_transform=True,
            **kwargs)
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

def make_initial_template_step1(src_sub='S03'):
    """ Apply Bias field correction the the T1 image"""
    src = f'{base_dir}/anatomicals/{src_sub}/{src_sub}_T1w.nii'
    src_img = ants.image_read(src) 
    fname_w = f'{base_dir}/bold_normalization/templateS03_building/{src_sub}_T1w_n3.nii'
    src_img4 = ants.n3_bias_field_correction(src_img)
    ants.image_write(src_img4,fname_w)

def make_initial_template_step2(src_sub='S03'): 
    """ Get the Transform from the bias corrected s03_T1 to the MNI2009c template
    This does not work well with mattes loss function, but much better with cross-correlation"""
    template  = f'{base_dir}/bold_normalization/templates/tpl-MNI152NLin2009cSym_res-1_T1w.nii'
    src = f'{base_dir}/bold_normalization/templateS03_building/{src_sub}_T1w_n3.nii'

    prefix = f'{base_dir}/bold_normalization/templateS03_building/T1to2009Sym_n3cc_{src_sub}_'
    trg_img = ants.image_read(template)
    src_img = ants.image_read(src) 
    mytx = ants.registration(fixed=trg_img , moving=src_img,type_of_transform='SyNCC',outprefix=prefix,write_composite_transform=True)
    fname_w = f'{base_dir}/bold_normalization/templateS03_building/w{src_sub}_T1w_n3cc.nii'
    wsrc_img = mytx['warpedmovout']
    ants.image_write(wsrc_img,fname_w)
    
def make_initial_template_step2a(src_sub='S03'):
    """ Attempt for better alignment of S03_mean_sbref to S03_T1
    This does not work well in the current attempts, so decide to skip this step"""
    src = f'{base_dir}/bold_normalization/individual/{src_sub}_mean_sbref.nii'
    template = f'{base_dir}/anatomicals/{src_sub}/{src_sub}_T1w.nii'
    prefix = f'{base_dir}/bold_normalization/templateS03_building/boldToT1_{src_sub}_'
    trg_img = ants.image_read(template)
    src_img = ants.image_read(src) 
    mytx = ants.registration(fixed=trg_img , moving=src_img,type_of_transform='SyNBold',outprefix=prefix,write_composite_transform=True)
    fname_w = f'{base_dir}/bold_normalization/templateS03_building/{src_sub}_mean_bold_corrected.nii'
    wsrc_img = mytx['warpedmovout']
    ants.image_write(wsrc_img,fname_w)

def make_initial_template_step3(src_sub='S03'):
    """ Apply T1to2009cSym transform to the bold image - and symmetrize it"""
    xfm_name = f'{base_dir}/bold_normalization/templateS03_building/T1to2009Sym_n3cc_{src_sub}_Composite.h5'
    # Read the transform and apply it to the bold image
    tx=ants.read_transform(xfm_name)
    bold_img=ants.image_read(f'{base_dir}/bold_normalization/individual/{src_sub}_mean_sbref.nii')

    template  = f'{base_dir}/bold_normalization/templates/tpl-MNI152NLin2009cSym_res-1_T1w.nii'
    src = f'{base_dir}/anatomicals/{src_sub}/{src_sub}_T1w.nii'

    prefix = f'{base_dir}/bold_normalization/templateS03_building/T1to2009Sym_{src_sub}_'
    trg_img = ants.image_read(template)
    trg_crp = ants.crop_indices( trg_img, (0,22,0), (193,160,164) )
    ntemp_img  =tx.apply_to_image(bold_img,reference=trg_crp)
    ntemp_name = f'{base_dir}/bold_normalization/templateS03_building/tpl-{src_sub}_bold.nii'
    ants.image_write(ntemp_img,ntemp_name)
    outf = f'{base_dir}/bold_normalization/templates/tpl-{src_sub}Sym_bold.nii'
    symmetrize_image(ntemp_name,outf)
    pass

def symmetrize_image(src,outname):
    img = nb.load(src)
    data = img.get_fdata()
    data = (data + np.flip(data,axis=0))/2
    img = nb.Nifti1Image(data, img.affine)
    nb.save(img,outname)
    return img

def symmetrize_template(src_sub='S03'):
    inf = f'{base_dir}/bold_normalization/deformed/S03Sym_CC/wavrg.nii'
    outf = f'{base_dir}/bold_normalization/templates/tpl-SyNSym_bold.nii'
    symmetrize_image(inf,outf)

def normalize_to_MNI2009c():
    tmp_dir = f'{base_dir}/bold_normalization/templates'
    src = f'{tmp_dir}/tpl-SyNSym_bold.nii'
    trg = f'{tmp_dir}/tpl-MNI152NLin2009cSym_res-1_T2w.nii.gz'
    mask_trg=f'{tmp_dir}/tpl-MNI152NLin2009cSymC_desc-cereb_mask.nii'
    prefix = f'{tmp_dir}/xfm_SyNSym_to_2009c/xfm_SynSym_to_2009cT2'
    trg_img = ants.image_read(trg)
    src_img = ants.image_read(src) 
    msk_img = ants.image_read(mask_trg)
    mytx = ants.registration(fixed=trg_img , moving=src_img, mask = msk_img,type_of_transform='SyN',outprefix=prefix,write_composite_transform=True)
    fname_w = f'{tmp_dir}/xfm_SyNSym_to_2009c/w_SynSym_to_2009cT2m.nii'
    wsrc_img = mytx['warpedmovout']
    ants.image_write(wsrc_img,fname_w)

def ants_transform_to_deformation_field(xfm,src_img,ref_img):
    """Converts the ANTs transform to a deformation field"""
    ref_ants = ants.image_read(ref_img)
    src_img = nb.load(src_img)

    # Make a image in src space with the voxel coordinates 
    Xv,Yv,Zv = np.meshgrid(range(src_img.shape[0]),
                           range(src_img.shape[1]),
                           range(src_img.shape[2]),indexing='ij')
    X,Y,Z = nt.affine_transform(Xv,Yv,Zv,src_img.affine)
    data = np.stack([X,Y,Z],axis=3)
    coord_img  = nb.Nifti1Image(data,src_img.affine)
    coord_ants = ants.from_nibabel(coord_img)

    # Apply the transform to this coordinate image
    y_ants = ants.apply_transforms(ref_ants, coord_ants, xfm, 
                          interpolator='linear', 
                          imagetype=3, 
                          defaultvalue=np.nan)
    # Insert a dummy 4th dimension
    y_img = y_ants.to_nibabel()
    data = y_img.get_fdata()
    data = np.expand_dims(data,axis=3)
    y_img  = nb.Nifti1Image(data,y_img.affine)
    return y_img


def make_deformation_fields(name,template,subj=None):
    """ Generates SPM-compatible deformation fields from the ANTs transforms"""
    if subj is None: 
        subj = pinfo['participant_id']
        subj = subj[pinfo.good==1]
    xfm_dir = f'{base_dir}/bold_normalization/transforms/{name}'
    ref_img = f'{base_dir}/bold_normalization/templates/{template}'

    for src_sub in subj:
        src = f'{base_dir}/bold_normalization/individual/{src_sub}_mean_sbref.nii'
        if not os.path.isfile(src):
            src = f'{base_dir}/bold_normalization/individual/{src_sub}_mean_bold.nii'

        xfm = f'{xfm_dir}/xfm_{src_sub}_composite.h5'
        def_field = ants_transform_to_deformation_field([xfm],src,ref_img)
        fname=f'{xfm_dir}/y_{src_sub}.nii'
        nb.save(def_field,fname)
        pass

def test_bold_images(subj=None):
    if subj is None: 
        subj = pinfo['participant_id']
        subj = subj[pinfo.good==1]

    for src_sub in subj:
        src = f'{base_dir}/bold_normalization/individual/{src_sub}_mean_sbref.nii'
        if not os.path.isfile(src):
            src = f'{base_dir}/bold_normalization/individual/{src_sub}_mean_bold.nii'

        nb_img = nb.load(src)

        SA,sa=nb_img.header.get_sform(coded=True)
        SQ,sq=nb_img.header.get_qform(coded=True)

        print(f'{src_sub} sform: {sa} qform: {sq}')
        if (np.isclose(SA,SQ,1e-4,1e-4)).all(): 
            print('sform and qform are the same')
        else:
            print('sform and qform are different')
            nb_img.header.set_qform(nb_img.header.get_sform(),code='aligned')
            nb.save(nb_img,src)

if __name__ == '__main__':
    #  make_initial_template_step3()
    # src_sub='S16'
    # make_initial_template()
    # test_bold_images()
    
    # kw={'verbose':True}
    # normalize_bold_all('MNI2009c_T1bold','tpl-MNI152NLin2009cSym_res-1_T1w.nii',tot='SyNBold',kwargs=kw)
    # normalize_bold_all('S03Sym_Syn','tpl-S03Sym_bold.nii',tot='SyN')
    # normalize_bold_all('S03Sym_CC','tpl-S03Sym_bold.nii',tot='SyNCC')
    # normalize_bold_all('SynSym_CC','tpl-SyNSym_bold.nii',tot='SyNCC',subj=['S08'])
    # make_deformation_fields('SynSym_CC','tpl-SynSym_bold.nii',subj=['S08'])
    # make_deformation_fields('S03Sym_Syn','tpl-S03Sym_bold.nii')
    #  Save to nifti
    normalize_to_MNI2009c()
    pass 