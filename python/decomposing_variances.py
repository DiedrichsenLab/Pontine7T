import nibabel
import numpy 
import pandas 
from Functional_Fusion.dataset import decompose_pattern_into_group_indiv_noise

data_dir = '/Volumes/diedrichsen_data$/data/Cerebellum/Pontine7T/RegionOfInterest/data/group'
 
 #appending files: 10 conditions, 16 runs; converting to tensor subj x cond x voxels

#
def get_structure_data(structure='cerebellum_gray'):
    T = pandas.read_csv('/Volumes/diedrichsen_data$/data/Cerebellum/Pontine7T/participants_no_s99.tsv', sep='\t')
    A = []
    for i in T.participant_id:
        file_path = f'{data_dir}/beta_glm2_{structure}_{i}.dscalar.nii'
        cifti = nibabel.load(file_path)
        A.append(cifti.get_fdata())
    to_tensor = numpy.array(A)
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

    unique_conditions = numpy.unique(cond_vec)
    unique_partitions = numpy.unique(part_vec)

    n_conditions = unique_conditions.size
    n_partitions = unique_partitions.size
   
    data = numpy.zeros((n_subjects, n_partitions, n_conditions, n_voxels))

    for p,partI in enumerate(unique_partitions):
        for c,condI in enumerate(unique_conditions):
            trial_inds = numpy.where(numpy.logical_and(cond_vec == condI, part_vec == partI))
            data[:, p, c, :] = flat_data[:, trial_inds, :].squeeze()
            
    return data

if __name__=='__main__':
    flat_data = get_structure_data()

    cond_vec = numpy.tile(numpy.arange(1,11),16)

    part_vec = numpy.repeat(numpy.arange(1,17), 10)

    tensor_4d = flat2ndarray(flat_data, cond_vec, part_vec)

    tensor_no_nans = numpy.nan_to_num(tensor_4d)

    tensor_avg_cond = tensor_no_nans.mean(axis=1, keepdims=1)

    tensor_subtract = tensor_no_nans - tensor_avg_cond

    variances = decompose_pattern_into_group_indiv_noise(tensor_subtract, criterion='global')
    print("global variances:", variances)

    #locating voxels with missing data

    missing_data = []

    for i in range(0, tensor_4d.shape[3]):
        if numpy.any(numpy.isnan(tensor_4d[:, :, :, i])):
            missing = i
            missing_data.append(missing)

   # print(missing_data)

    print("done")