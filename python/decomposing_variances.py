import nibabel
import numpy 
import pandas 
import sys
import site 

path_to_func_fusion = '/Users/incehusain/Documents/GitHub/Functional_Fusion'

site.addsitedir(path_to_func_fusion)

from Functional_Fusion.dataset import decompose_pattern_into_group_indiv_noise

#loading data

cifti = nibabel.load('/Users/incehusain/decomposing_variances/beta_glm2_cerebellum_gray_S03.dscalar.nii')


print(cifti.shape)

data = cifti.get_fdata()

unique_data = numpy.unique(data)



