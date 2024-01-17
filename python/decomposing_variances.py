import nibabel
import numpy 
import pandas 
import sys
import site 

path_to_func_fusion = '/Users/incehusain/Documents/GitHub/Functional_Fusion'

site.addsitedir(path_to_func_fusion)

from Functional_Fusion.dataset import decompose_pattern_into_group_indiv_noise

#appending 2 files: 160 conditions, 2385 voxels each; and converting to tensor of form S x N x P

def get_data(structure='dentate'):
    A = []
    B = A
    for i in range(3,5):
        file_path = f'/Users/incehusain/decomposing_variances/beta_glm2_{structure}_S0{i}.dscalar.nii'
        cifti = nibabel.load(file_path)
        A.append(cifti.get_fdata())
        to_tensor = numpy.array(B)
    print("tensor:", B)
    print("shape of tensor:", to_tensor.shape)
    print("number of dimensions:", numpy.ndim(B))
    print('number of voxels:', len(B[0][0]))
    print('number of conditions:', len(B[0]))
    print('subject 2, condition 1, voxel 3:', B[1][0][2])
    print('subject 1, condition 1, voxel 1:', B[0][0][0])

get_data()



