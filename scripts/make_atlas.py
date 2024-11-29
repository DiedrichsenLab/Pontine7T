import nibabel as nb
import torch as pt
import numpy as np
import Functional_Fusion.atlas_map as am
import Functional_Fusion.dataset as ds
import HierarchBayesParcel.arrangements as arr
import HierarchBayesParcel.emissions as em
import HierarchBayesParcel.arrangements as ar
import HierarchBayesParcel.emissions as em
import HierarchBayesParcel.full_model as fm
import HierarchBayesParcel.util as ut
import matplotlib.pyplot as plt
import SUITPy as suit 
import nitools as nt 
import Functional_Fusion.plot as plot
from matplotlib.colors import ListedColormap
from scripts import decomposing_variances as ac

base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion' 
atlas_dir = base_dir + '/Atlases/tpl-MNI152NLin2009cSymC'


def build_emission_mdtb(K,P,atlas='MNISymCereb2'):
    data, info,ds_obj = ds.get_dataset(base_dir,'MDTB',atlas=atlas,type='CondRun',sess='all',subj=None)
    cond_v = info['cond_num_uni']
    part_v = info['run']

    # Make a design matrix
    X= ut.indicator(cond_v)
    # Build an emission model
    em_model = em.MixVMF(K=K,P=P, X=X,part_vec=part_v)
    em_model.initialize(data)
    return em_model

def build_emission_language(K,P,atlas='MNISymCereb2'):
    data, info, ds_obj = ds.get_dataset(base_dir,'Language',atlas=atlas,type='CondRun', sess='ses-localizer_cond_fm', subj=None)
    cond_v = info['task']
    part_v = info['run']

    # Make a design matrix
    X= ut.indicator(cond_v)

    data = ds.remove_baseline(data,part_v)

    # Build an emission model
    em_model = em.MixVMF(K=K,P=P, X=X,part_vec=part_v)
    em_model.initialize(data)
    return em_model


def build_emission_pontine(K,P,atlas='MNISymCereb2'):
    data_dir = '/Volumes/diedrichsen_data$/data/Cerebellum/Pontine7T/RegionOfInterest_BOLDMNI/data/group'
    if atlas=='MNISymCereb2': 
        pontine_flat_data = ac.get_structure_data(structure='cereb_gray', data_dir=data_dir )
    elif atlas == 'MNISymDentate1':
        pontine_flat_data = ac.get_structure_data(structure='cereb_gray',  data_dir=data_dir)
 
    cond_v = np.tile(np.arange(1,11),16)

    part_v = np.repeat(np.arange(1,17), 10)

    tensor_pontine = ac.flat2ndarray(pontine_flat_data, cond_v, part_v)

    tensor_no_nans = np.nan_to_num(tensor_pontine)

    tensor_avg_cond = tensor_no_nans.mean(axis=3, keepdims=1) #this is the mean activity pattern 

    data_p = tensor_no_nans - tensor_avg_cond

    data = data_p.reshape(16,160,18290)

    X = ut.indicator(cond_v)

    em_model = em.MixVMF(K=K,P=P,X=X,part_vec=part_v)
    em_model.initialize(data)
    return em_model

def estimate_emission_models():
    """ Estimate emission models for dataset
    """
    # Get the atlas for the cerebellar cortex
    atlas, _ = am.get_atlas('MNISymCereb2')

    # Load the probabilistic atlas into that space 
    atlas_fname = 'atl-NettekovenSym32_space-MNI152NLin2009cSymC_probseg.nii'
    U = atlas.read_data(atlas_dir + '/' + atlas_fname).T
    
    # Build the arrangement model 
    ar_model = ar.build_arrangement_model(U, prior_type='prob', atlas=atlas)
        
    # Load the data using functional fusion 

    em_model1 = build_emission_mdtb(ar_model.K,atlas.P,atlas='MNISymCereb2')
    em_model2 = build_emission_language(ar_model.K,atlas.P,atlas='MNISymCereb2')
    em_model3  = build_emission_pontine(ar_model.K,atlas.P,atlas='MNISymCereb2')
 
    # Build the full model: The emission models are passed as a list, as usually we have multiple data sets
    M = fm.FullMultiModel(ar_model, [em_model1, em_model2, em_model3])

    # Attach the data to the model - this is done for speed
    # The data is passed as a list with on element per data set
    M.initialize() # Initialize the default that each em has different subjects

    # Now we can run the EM algorithm
    M, _, _, _ = M.fit_em(iter=200, tol=0.01,
        fit_arrangement=False,fit_emission=True,first_evidence=False)
    
    Vs = [i.V for i in M.emissions]

    return Vs    

def estimate_new_atlas(cereb_Vs): 
     
    # Get atlas for dentate 

     atlas, _ = am.get_atlas('MNISymDentate1')    

     ar_model = ar.ArrangeIndependent(K=32, P=atlas.P, spatial_specific=True, remove_redundancy=False)
         
     em_model1 = build_emission_mdtb(ar_model.K,atlas.P,atlas='MNISymCereb2')
     em_model2 = build_emission_language(ar_model.K,atlas.P,atlas='MNISymCereb2')
     em_model3  = build_emission_pontine(ar_model.K,atlas.P,atlas='MNISymCereb2')

     em_model1.V = cereb_Vs[0]
     em_model2.V = cereb_Vs[1]
     em_model3.V = cereb_Vs[2]

     em_model1.set_param_list(['kappa'])
     em_model2.set_param_list(['kappa'])
     em_model3.set_param_list(['kappa'])

     M= fm.FullMultiModel(ar_model, [em_model1, em_model2, em_model3])

     M.initialize()

     M, _, _, _ = M.fit_em(iter=200, tol=0.01,
        fit_arrangement=True,fit_emission=True,first_evidence=True)
     
     return M


if __name__ == '__main__':
    all_Vs = estimate_emission_models()

    dentate_group_map = estimate_new_atlas(all_Vs)

    # Load colormap and labels

    lid,cmap,names = nt.read_lut('/Volumes/diedrichsen_data$/data/FunctionalFusion/Atlases/tpl-MNI152NLin2009cSymC/atl-NettekovenSym32.lut')

    cmap_string = ListedColormap(cmap)

    #plot group map 
    data = dentate_group_map.arrange.marginal_prob().numpy()
    #wta = np.argmax(data,axis=0)
    wta = np.argmax(data, axis=0) 
    wta += 1

    wta_int32 = wta.astype(np.int32)
    
    dentate_parcellation = plot.plot_dentate(wta_int32,cscale=[0,32],cmap=cmap_string)

    
    
    Vs = [em.V for em in M.emissions]
    #plt.imshow(Vs[0])

   # plt.yticks(info.cond_name[info.run==1].values)
   # plt.y()


   # arrangement = estimate_new_atlas(ems)
    
    # arrangement.marginal_prob() will give the probability of each parcel for each voxel in the new atlas. (KxP)