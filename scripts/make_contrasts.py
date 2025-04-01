import nibabel as nb
import numpy as np
import pandas 
import seaborn as sns
import matplotlib.pyplot as plt
import Functional_Fusion.dataset as ds
import Functional_Fusion.atlas_map as am 
import pandas as pd
from scripts import decomposing_variances

base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion' 
atlas_dir = base_dir + '/Atlases/tpl-MNI152NLin2009cSymC'

def make_pontine_contrasts(atlas='MNISymDentate1'):
    data, info, ds_obj = ds.get_dataset(base_dir,'Pontine',atlas=atlas,
                                        type='CondRun', sess='ses-s1', 
                                        subj=None)
    cond_v = info['task']
    part_v = info['run']

    data = ds.remove_baseline(data,part_v)

    flat_data = decomposing_variances.flat2ndarray(data, cond_v, part_v)

    flat_data = np.nan_to_num(flat_data)

    cond_avg = np.mean(flat_data, axis=1)

    num_subj = cond_avg.shape[0]
    num_cond = cond_avg.shape[1]

    contrast_per_subj = np.zeros((num_subj,num_cond, cond_avg.shape[-1]))

    #make contrast vectors

    T = pd.read_csv(f"{base_dir}/Pontine/derivatives/sub-01/data/sub-01_ses-s1_CondRun.tsv", sep='\t')

    contrast_names_all = T['taskName'].tolist()

    contrast_names = list(set(contrast_names_all))

    contrast_names = list(dict.fromkeys(contrast_names_all))

    contrast = []

    for i in range(num_cond):  # Loop over all conditions (tasks)
        c = [1 if j == i else -1 / (num_cond - 1) for j in range(num_cond)]
        contrast.append(c)

    contrast = np.array(contrast)

    for i,c in enumerate(contrast):
        for s in range(num_subj):
            contrast_per_subj[s,i,:] = np.dot(c, cond_avg[s, :, :])

    CON = np.mean(contrast_per_subj,axis=0)

    STD = np.std(contrast_per_subj,axis=0)
        
    t = CON/(STD/np.sqrt(num_subj))
        
    # Get one example cifti-file for the header 
    
    ref_img=nb.load(f"{base_dir}/Pontine/derivatives/sub-01/data/sub-01_space-MNISymDentate1_ses-s1_CondRun.dscalar.nii")
    
    bm = ref_img.header.get_axis(1)

    row_axis = nb.cifti2.ScalarAxis(contrast_names)
    header = nb.Cifti2Header.from_axes((row_axis,bm))

    con_img = nb.Cifti2Image(dataobj=CON[:, :], header=header)
    t_img = nb.Cifti2Image(dataobj=t[:, :], header=header)

    con_filename = f'/Volumes/diedrichsen_data$/data/Cerebellum/Pontine7T/atlases/dentate/contrasts/condavg_contrast.dscalar.nii'
    t_filename = f'/Volumes/diedrichsen_data$/data/Cerebellum/Pontine7T/atlases/dentate/contrasts/condavg_Tstat.dscalar.nii'

        # Save the contrast and T-statistic images
    nb.save(con_img, con_filename)
    nb.save(t_img, t_filename)

def make_language_contrasts(atlas='MNISymDentate1'):
    data, info, ds_obj = ds.get_dataset(base_dir,'Language',atlas=atlas,
                                        type='CondRun', sess='ses-localizer_cond_fm', 
                                        subj=None)
    cond_v = info['task']
    part_v = info['run']

    data = ds.remove_baseline(data,part_v)

    flat_data = decomposing_variances.flat2ndarray(data, cond_v, part_v)

    flat_data = np.nan_to_num(flat_data)

    cond_avg = np.mean(flat_data, axis=1)

    num_subj = cond_avg.shape[0]
    num_cond = cond_avg.shape[1]

    contrast_per_subj = np.zeros((num_subj,num_cond, cond_avg.shape[-1]))

    #make contrast vectors

    T = pd.read_csv(f"{base_dir}/Language/derivatives/sub-01/data/sub-01_ses-localizer_cond_fm_CondRun.tsv", sep='\t')

    contrast_names_all = T['taskName'].tolist()

    contrast_names = list(set(contrast_names_all))

    contrast_names = list(dict.fromkeys(contrast_names_all))

    contrast = []

    for i in range(num_cond):  # Loop over all conditions (tasks)
        c = [1 if j == i else -1 / (num_cond - 1) for j in range(num_cond)]
        contrast.append(c)

    contrast = np.array(contrast)

    for i,c in enumerate(contrast):
        for s in range(num_subj):
            contrast_per_subj[s,i,:] = np.dot(c, cond_avg[s, :, :])

    CON = np.mean(contrast_per_subj,axis=0)

    STD = np.std(contrast_per_subj,axis=0)
        
    t = CON/(STD/np.sqrt(num_subj))
        
    # Get one example cifti-file for the header 
    
    ref_img=nb.load(f"{base_dir}/Language/derivatives/sub-01/data/sub-01_space-MNISymDentate1_ses-localizer_cond_fm_CondRun.dscalar.nii")
    
    bm = ref_img.header.get_axis(1)

    row_axis = nb.cifti2.ScalarAxis(contrast_names)
    header = nb.Cifti2Header.from_axes((row_axis,bm))

    con_img = nb.Cifti2Image(dataobj=CON[:, :], header=header)
    t_img = nb.Cifti2Image(dataobj=t[:, :], header=header)

    con_filename = f'/Volumes/diedrichsen_data$/data/Cerebellum/Pontine7T/atlases/dentate/contrasts/condavg_contrast.dscalar.nii'
    t_filename = f'/Volumes/diedrichsen_data$/data/Cerebellum/Pontine7T/atlases/dentate/contrasts/condavg_Tstat.dscalar.nii'

        # Save the contrast and T-statistic images
    nb.save(con_img, con_filename)
    nb.save(t_img, t_filename)


def make_mdtb_contrasts(atlas='MNISymDentate1'):
    data, info, ds_obj = ds.get_dataset(base_dir,'MDTB',atlas=atlas,
                                        type='CondRun', sess='ses-s1', 
                                        subj=None)
    cond_v = info['cond_num_uni']
    part_v = info['run']

    flat_data = decomposing_variances.flat2ndarray(data, cond_v, part_v)

    flat_data = np.nan_to_num(flat_data)

    cond_avg = np.mean(flat_data, axis=1)

    num_subj = cond_avg.shape[0]
    num_cond = cond_avg.shape[1]

    contrast_per_subj = np.zeros((num_subj,num_cond, cond_avg.shape[-1]))

    #make contrast vectors

    T = pd.read_csv(f"{base_dir}/MDTB/derivatives/sub-02/data/sub-02_ses-s1_CondRun.tsv", sep='\t')

    contrast_names_all = T['task_name'].tolist()

    contrast_names = list(set(contrast_names_all))

    contrast_names = list(dict.fromkeys(contrast_names_all))

    contrast = []

    for i in range(num_cond):  # Loop over all conditions (tasks)
        c = [1 if j == i else -1 / (num_cond - 1) for j in range(num_cond)]
        contrast.append(c)

    contrast = np.array(contrast)

    for i,c in enumerate(contrast):
        for s in range(num_subj):
            contrast_per_subj[s,i,:] = np.dot(c, cond_avg[s, :, :])

    CON = np.mean(contrast_per_subj,axis=0)

    STD = np.std(contrast_per_subj,axis=0)
        
    t = CON/(STD/np.sqrt(num_subj))
        
    # Get one example cifti-file for the header 
    
    ref_img=nb.load(f"{base_dir}/MDTB/derivatives/sub-02/data/sub-02_space-MNISymDentate1_ses-s1_CondRun.dscalar.nii")
    
    bm = ref_img.header.get_axis(1)

    row_axis = nb.cifti2.ScalarAxis(contrast_names)
    header = nb.Cifti2Header.from_axes((row_axis,bm))

    con_img = nb.Cifti2Image(dataobj=CON[:, :], header=header)
    t_img = nb.Cifti2Image(dataobj=t[:, :], header=header)

    con_filename = f'/Volumes/diedrichsen_data$/data/Cerebellum/Pontine7T/atlases/dentate/contrasts/condavg_contrast.dscalar.nii'
    t_filename = f'/Volumes/diedrichsen_data$/data/Cerebellum/Pontine7T/atlases/dentate/contrasts/condavg_Tstat.dscalar.nii'

        # Save the contrast and T-statistic images
    nb.save(con_img, con_filename)
    nb.save(t_img, t_filename)



    #beyond this point is garbage 

def make_contrast_vectors_handedness():

    contrast = np.array([[-1,1],[1,-1]])

    return contrast 

def group_analysis_handedness(contrast):
    
    Y = get_structure_data(structure='cereb_gray', data_dir='/Volumes/diedrichsen_data$/data/Cerebellum/Pontine7T/RegionOfInterest_BOLDMNI/data/group_smoothed') 
    cond_vec = np.tile(np.arange(1,11),16)
    part_vec = np.repeat(np.arange(1,17), 10)
    
    Y_array = flat2ndarray(Y, cond_vec, part_vec)

    left_fingseq = Y_array[:, 4:12, 2:3, :]  #this selects tasks with finger sequence (axis 2) for runs 5 through 12 (4:12)

    left_fingseq_avg =  np.mean(left_fingseq, axis=1)

    left_rest = Y_array[:, 4:12, 8:9, :]  #this selects rest for left hand 

    left_rest_avg = np.mean(left_rest, axis=1) #shape is subj x cond x voxels 

    right_fingseq = Y_array[:, np.r_[0:4, 12:16], 2:3, :] 

    right_fingseq_avg = np.mean(right_fingseq, axis=1)

    right_rest = Y_array[:, np.r_[0:4, 12:16], 8:9, :]

    right_rest_avg = np.mean(right_rest, axis=1)

    combined_array = np.stack([left_fingseq_avg, left_rest_avg, right_fingseq_avg, right_rest_avg], axis=1 )
    combined_array = np.squeeze(combined_array)

    num_subj = combined_array.shape[0]

    contrast_per_subj = np.zeros((num_subj,1,combined_array.shape[-1]))

    for s in range(num_subj):
        contrast_per_subj[s,:,:] = np.dot(contrast, combined_array[s, :, :])

    CON = np.mean(contrast_per_subj,axis=0)

    STD = np.std(contrast_per_subj,axis=0)
        
    t = CON/(STD/np.sqrt(num_subj))

    ref_img = nb.load("/Volumes/diedrichsen_data$/data/Cerebellum/Pontine7T/RegionOfInterest_BOLDMNI/data/group_smoothed/ref_cereb.dscalar.nii")
    bm = ref_img.header.get_axis(1)

    row_axis = nb.cifti2.ScalarAxis(["Handedness contrast"])
    header = nb.Cifti2Header.from_axes((row_axis,bm))

    con_img = nb.Cifti2Image(dataobj=CON[:, :], header=header)
    t_img = nb.Cifti2Image(dataobj=t[:, :], header=header)

    con_filename = f'/Volumes/diedrichsen_data$/data/Cerebellum/Pontine7T/RegionOfInterest_BOLDMNI/data/group_avg/cond_handedness_contrast.dscalar.nii'
    t_filename = f'/Volumes/diedrichsen_data$/data/Cerebellum/Pontine7T/RegionOfInterest_BOLDMNI/data/group_avg/cond_handedness_Tstat.dscalar.nii'

        # Save the contrast and T-statistic images
    nb.save(con_img, con_filename)
    nb.save(t_img, t_filename)
        


def make_contrast_vectors(): 

    T = pandas.read_csv('/Volumes/diedrichsen_data$/data/Cerebellum/Pontine7T/GLM_firstlevel_2/S17/SPM_info.tsv', sep='\t')

    contrast_names_all = T['taskName'].tolist()

    contrast_names = list(set(contrast_names_all))

    contrast_names = list(dict.fromkeys(contrast_names_all))

    num_cond = 10

    contrast = []
   
   #a 160 x 160 matrix where every row is a 
    
    inst_vector = [int(1) if j % 10 ==0
                   else -1/(num_cond-1) for j in range(0,num_cond)]
    contrast.append(inst_vector)
    
    for i in range(1,num_cond):
        c = [1-1/(num_cond-1) if j == i and j % 10 != 0
             else int(0) if j == 0
             else -1/(num_cond-1) if j % 10 !=0 
             else 0 for j in range(0,num_cond)]
        
        contrast.append(c)

    contrast = np.array(contrast)
    
    return contrast, contrast_names

def group_analysis(contrast,contrast_names):

    # Does group analysis 
    
    Y = get_structure_data(structure='dentate', data_dir='/Volumes/diedrichsen_data$/data/Cerebellum/Pontine7T/RegionOfInterest_BOLDMNI/data/group_smoothed') 

    cond_vec = np.tile(np.arange(1,11),16)

    part_vec = np.repeat(np.arange(1,17), 10)
    
    Y_array = flat2ndarray(Y, cond_vec, part_vec)

    Y_array_cond_avg = np.mean(Y_array, axis=1)

    num_subj = Y_array_cond_avg.shape[0]
    num_cond = Y_array_cond_avg.shape[1]

    contrast_per_subj = np.zeros((num_subj,num_cond, Y_array_cond_avg.shape[-1]))


    for i,c in enumerate(contrast):
        for s in range(num_subj):
            contrast_per_subj[s,i,:] = np.dot(c, Y_array_cond_avg[s, :, :])

    CON = np.mean(contrast_per_subj,axis=0)

    STD = np.std(contrast_per_subj,axis=0)
        
    t = CON/(STD/np.sqrt(num_subj))
        
    # Get one exmample cifti-file for the header 
    
    ref_img = nb.load("/Volumes/diedrichsen_data$/data/Cerebellum/Pontine7T/RegionOfInterest_BOLDMNI/data/group_smoothed/ref_dentate.dscalar.nii")
    bm = ref_img.header.get_axis(1)

    row_axis = nb.cifti2.ScalarAxis(contrast_names)
    header = nb.Cifti2Header.from_axes((row_axis,bm))

    con_img = nb.Cifti2Image(dataobj=CON[:, :], header=header)
    t_img = nb.Cifti2Image(dataobj=t[:, :], header=header)

    con_filename = f'/Volumes/diedrichsen_data$/data/Cerebellum/Pontine7T/RegionOfInterest_BOLDMNI/data/group_avg/condavg_contrast.dscalar.nii'
    t_filename = f'/Volumes/diedrichsen_data$/data/Cerebellum/Pontine7T/RegionOfInterest_BOLDMNI/data/group_avg/condavg_Tstat.dscalar.nii'

        # Save the contrast and T-statistic images
    nb.save(con_img, con_filename)
    nb.save(t_img, t_filename)

if __name__ == '__main__':

    pontine = make_mdtb_contrasts()

    print("YO")