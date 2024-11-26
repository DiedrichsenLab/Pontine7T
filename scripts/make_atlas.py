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


base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion' 
atlas_dir = base_dir + '/Atlases/tpl-MNI152NLin2009cSymC'

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
        
    # Load the data using functional fusion (for pontine, follow step 4)
    data, info,ds_obj = ds.get_dataset(base_dir,'MDTB',atlas='MNISymCereb2',type='CondRun',sess='all',subj=[0,1,2])
    cond_v = info['cond_num_uni']
    part_v = info['run']

    # Step 5: Fit new emission model to the data

    # K is the number of parcels
    K = ar_model.K
    # Make a design matrix
    X= ut.indicator(cond_v)
    # Build an emission model
    em_model1 = em.MixVMF(K=K,P=atlas.P, X=X,part_vec=part_v)

    # Build the full model: The emission models are passed as a list, as usually we have multiple data sets
    M = fm.FullMultiModel(ar_model, [em_model1])

    # Attach the data to the model - this is done for speed
    # The data is passed as a list with on element per data set

    M.initialize([data])

    # Now we can run the EM algorithm
    M, _, _, _ = M.fit_em(iter=200, tol=0.01,
        fit_arrangement=False,fit_emission=True,first_evidence=False)
    
    Vs = [i.V for i in M.emissions]

    return Vs[0]    

def estimate_new_atlas(cereb_Vs): 
    # Get atlas for dentate 

     dentate_atlas, _ = am.get_atlas('MNISymDentate1')    
    
    # Get all three datasets + info for the dentate

     data, info,ds_obj = ds.get_dataset(base_dir,'MDTB',atlas='MNISymDentate1',type='CondRun',sess='all',subj=[0,1,2])
     cond_v = info['cond_num_uni']
     part_v = info['run']
     
     ar_model = ar.ArrangeIndependent(K=32, P=dentate_atlas.P, spatial_specific=True, remove_redundancy=False)

    
    # Make a design matrix
     
     X= ut.indicator(cond_v)
     
     em_model1 = em.MixVMF(K=32, P=dentate_atlas.P, X=X, part_vec=part_v)

     for k in range(len(cereb_Vs)):
         em_model1.V[k] = cereb_Vs[k]

     em_model1.set_param_list(['kappa'])

     M= fm.FullMultiModel(ar_model, [em_model1])

     M.initialize([data])

     M, _, _, _ = M.fit_em(iter=200, tol=0.01,
        fit_arrangement=True,fit_emission=True,first_evidence=True)
     
     return ar_model
     
     #corresponding with the dataset

    # Set the parameter list of the emission models to kappa only to avoid fitting of V. 
     
    #for em in em_models:
     #   em.set_param_list(['kappa'])
    #M = fm.FullMultiModel(ar_model, em_models)


    #M.initialize([data])

    #M, _, _, _ = M.fit_em(iter=200, tol=0.01,
     #   fit_arrangement=True,fit_emission=True,first_evidence=True)

    #return ar_model

if __name__ == '__main__':
    Vs = estimate_emission_models()

    dentate_group_map = estimate_new_atlas(Vs)

    
    Vs = [em.V for em in M.emissions]
    #plt.imshow(Vs[0])

   # plt.yticks(info.cond_name[info.run==1].values)
   # plt.y()


   # arrangement = estimate_new_atlas(ems)
    
    # arrangement.marginal_prob() will give the probability of each parcel for each voxel in the new atlas. (KxP)