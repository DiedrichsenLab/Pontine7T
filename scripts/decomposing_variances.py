import nibabel as nb
import numpy as np
import pandas 
import seaborn as sns
import matplotlib.pyplot as plt
from Functional_Fusion.dataset import decompose_pattern_into_group_indiv_noise

#data_dir = '/Volumes/diedrichsen_data$/data/Cerebellum/Pontine7T/RegionOfInterest_SUIT/data/group'
 
 #appending files: 10 conditions, 16 runs; converting to tensor subj x cond x voxels

#
def get_structure_data(structure='pontine', data_dir='/Volumes/diedrichsen_data$/data/Cerebellum/Pontine7T/RegionOfInterest_BOLD/data/group'):
    T = pandas.read_csv('/Volumes/diedrichsen_data$/data/Cerebellum/Pontine7T/participants.tsv', sep='\t')
    A = []
    for i, good_value in zip(T.participant_id, T.good):
        if good_value==1:
            file_path = f'{data_dir}/beta_glm2_{structure}_{i}.dscalar.nii'
            cifti = nb.load(file_path)
            A.append(cifti.get_fdata())
        to_tensor = np.array(A)
    return to_tensor

def flat2ndarray(flat_data, cond_vec, part_vec):
    """
    convert flat data (n_subjects x n_trials x n_voxels) into a 4d ndarray (n_subjects x n_conditions x n_partitions x n_voxels)

    Args:
        flat_data:
        part_vec:
        cond_vec:

    Returns:
        data

    """

    [n_subjects, n_trials, n_voxels] = flat_data.shape

    unique_conditions = np.unique(cond_vec)
    unique_partitions = np.unique(part_vec)

    n_conditions = unique_conditions.size
    n_partitions = unique_partitions.size
   
    data = np.zeros((n_subjects, n_partitions, n_conditions, n_voxels))

    for p,partI in enumerate(unique_partitions):
        for c,condI in enumerate(unique_conditions):
            trial_inds = np.where(np.logical_and(cond_vec == condI, part_vec == partI))
            data[:, p, c, :] = flat_data[:, trial_inds, :].squeeze()
            
    return data


def make_contrast_vectors(): 

    T = pandas.read_csv('/Volumes/diedrichsen_data$/data/Cerebellum/Pontine7T/pontine_taskConds_GLM_reordered.tsv', sep='\t')

    contrast_names_all = T['taskNames'].tolist()

    contrast_names = list(set(contrast_names_all))

    contrast_names = list(dict.fromkeys(contrast_names))

    num_conditions = 10

    contrast = []
    
    for i in range(num_conditions):
        c = [1 if j == i else -1/num_conditions for j in range(num_conditions)]
        contrast.append(c)

    contrast = np.array(contrast)
    
    return contrast, contrast_names



def group_analysis(contrast,contrast_names):

    # Does group analysis 
    
    Y = get_structure_data(structure='dentate', data_dir='/Volumes/diedrichsen_data$/data/Cerebellum/Pontine7T/RegionOfInterest_BOLD/data/group_smoothed') 

    cond_vec = np.tile(np.arange(1,11),16)

    part_vec = np.repeat(np.arange(1,17), 10)
    
    Y_flat = flat2ndarray(Y, cond_vec, part_vec)

    Y_flat_cond_avg = np.mean(Y_flat, axis=1)

    num_subj = Y_flat_cond_avg.shape[0]
    num_cond = Y_flat_cond_avg.shape[1]

    contrast_per_subj = np.zeros((num_subj,num_cond, Y_flat_cond_avg.shape[-1]))


    for i,c in enumerate(contrast):
        for s in range(num_subj):
            contrast_per_subj[s,i,:] = np.dot(c, Y_flat_cond_avg[s, :, :])

    CON = np.mean(contrast_per_subj,axis=0)

    STD = np.std(contrast_per_subj,axis=0)
        
    t = CON/(STD*np.sqrt(num_subj))
        
    # Get one exmample cifti-file for the header 
    
    ref_img = nb.load("/Volumes/diedrichsen_data$/data/Cerebellum/Pontine7T/RegionOfInterest_BOLD/group_avg/ref.dscalar.nii")
    bm = ref_img.header.get_axis(1)

    row_axis = nb.cifti2.ScalarAxis(contrast_names)
    header = nb.Cifti2Header.from_axes((row_axis,bm))

    con_img = nb.Cifti2Image(dataobj=CON[:, :], header=header)
    t_img = nb.Cifti2Image(dataobj=t[:, :], header=header)

    con_filename = f'/Volumes/diedrichsen_data$/data/Cerebellum/Pontine7T/RegionOfInterest_BOLD/group_avg/condavg_contrast.dscalar.nii'
    t_filename = f'/Volumes/diedrichsen_data$/data/Cerebellum/Pontine7T/RegionOfInterest_BOLD/group_avg/condavg_Tstat.dscalar.nii'

        # Save the contrast and T-statistic images
    nb.save(con_img, con_filename)
    nb.save(t_img, t_filename)


if __name__=='__main__':
    flat_data = get_structure_data()

    cond_vec = numpy.tile(numpy.arange(1,11),16)

    part_vec = numpy.repeat(numpy.arange(1,17), 10)

    tensor_4d = flat2ndarray(flat_data, cond_vec, part_vec)

    tensor_no_nans = numpy.nan_to_num(tensor_4d)

    tensor_std_cond = numpy.std(tensor_no_nans, axis=1, keepdims=1)


    variances= decompose_pattern_into_group_indiv_noise(tensor_no_nans, criterion='subject-wise')

    #locating voxels with missing data:

    missing_data = []

    for i in range(0, tensor_4d.shape[3]):
        if numpy.any(numpy.isnan(tensor_4d[:, :, :, i])):
            missing = i
            missing_data.append(missing)


            
