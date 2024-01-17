import nibabel
import numpy 
import pandas 
import sys
import site 
import itertools 

path_to_func_fusion = '/Users/incehusain/Documents/GitHub/Functional_Fusion'

site.addsitedir(path_to_func_fusion)

from Functional_Fusion.dataset import decompose_pattern_into_group_indiv_noise

#appending 2 files: 160 conditions, 2385 voxels each; and converting to tensor of form S x N x P

def get_data(structure='dentate'):
    A = []
    for i in range(3,5):
        file_path = f'/Users/incehusain/decomposing_variances/beta_glm2_{structure}_S0{i}.dscalar.nii'
        cifti = nibabel.load(file_path)
        A.append(cifti.get_fdata())
    to_tensor = numpy.array(A)
    print("tensor:", A)
    print("shape of tensor:", to_tensor.shape)
    print("number of dimensions:", numpy.ndim(A))
    print('number of voxels:', len(A[0][0]))
    print('number of conditions:', len(A[0]))
    print('subject 2, condition 2, voxel 3:', A[1][1][2])
    print('subject 1, condition 1, voxel 1:', A[0][0][0])
    return to_tensor

get_data()

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
   
    data = numpy.zeros((n_subjects, n_conditions, n_partitions, n_voxels))

    for c,condI in enumerate(unique_conditions):
        for p,partI in enumerate(unique_partitions):
            trial_inds = numpy.where(numpy.logical_and(cond_vec == condI, part_vec == partI))
            data[:, c, p, :] = numpy.nanmean(flat_data[:, trial_inds, :], axis=1).squeeze()
            
    return data

flat_data = get_data()

cond_vec = numpy.tile(numpy.arange(1,11),16)

part_vec = numpy.repeat(numpy.arange(1,17), 10)

result = flat2ndarray(flat_data, cond_vec, part_vec)
print("result subject x condition x repetition x voxel:", result)

print("last datapoint; subject 2, condition 10, repetition 16, voxel 2384:", result[1][9][15][2384])
print("first datapoint; subject 1, condition 1, repetition 1, voxel 1:", result[0][0][0][0])

variances = decompose_pattern_into_group_indiv_noise(result, criterion='voxel_wise')
print("variances:", variances)
print("done")

#there's an issue with the variance values: negative vg..
#per condition and global won't compute (give nans as results)
#the result data has nans in it: numpy.sum(numpy.isnan(result)) gives 1760/ 763200