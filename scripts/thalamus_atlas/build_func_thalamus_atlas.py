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
from scripts import decomposing_variances as dv

base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion_new' 
wk_dir = '/Users/incehusain/fs_projects'

def build_emissions(atlas='MNISymThalamus1'):
    
    K=atlas.K
    P=atlas.P

    data, info, ds_obj = ds.get_dataset(base_dir,'Language',atlas=atlas, sess='ses-localizerfm', subj=None, type='CondRun')
    cond_vec = info['task']
    part_vec = info['run']
    X= ut.indicator(cond_vec)
    em_model_language = em.MixVMF(K=K,P=P,X3=X,part_v3=part_vec)
    em_model_language.initialize(data)

    return em_model_language

def estimate_emission_models():
    """ Estimate emission models for dataset
    """
    atlas, _ = am.get_atlas('MNISymThalamus1')

    atlas_fname = 'group_mean_thalamus_prob_map.nii.gz'

    U = atlas.read_data(wk_dir + '/' + atlas_fname)
    
    ar_model = ar.build_arrangement_model(U, prior_type='prob', atlas=atlas)
            
    emissions = build_emissions(ar_model.K,atlas.P,atlas='MNISymThalamus1')
 
    M = fm.FullMultiModel(ar_model, [emissions])

    # Attach the data to the model - this is done for speed
    # The data is passed as a list with on element per data set

    M.initialize() 

    M, _, _, _ = M.fit_em(iter=200, tol=0.01,
        fit_arrangement=False,fit_emission=True,first_evidence=False)

    pt.save(M.emissions[1].V,f"{wk_dir}/V_Iglesias_Language.pt")

def estimate_new_atlas(): 
     
     atlas, _ = am.get_atlas('MNISymThalamus1')    

     ar_model = ar.ArrangeIndependent(K=32, P=atlas.P, spatial_specific=True, remove_redundancy=False)

     emissions = build_emissions(ar_model.K,atlas.P,atlas='MNISymThalamus1')

     emissions.V = pt.load(f"{wk_dir}/V_Iglesias_Language.pt")

     for em_model in emissions:
            em_model.set_param_list(['kappa'])
    
     M= fm.FullMultiModel(ar_model, [emissions])

     M.initialize()

     M, ll, _, _ = M.fit_em(iter=200, tol=0.01,
        fit_arrangement=True,fit_emission=True,first_evidence=True)
     
     Prob = M.arrange.marginal_prob().numpy()
     np.save(f"{wk_dir}/Prob_thalamus.npy",Prob)

     return M

if __name__ == '__main__':

    language = build_emissions()
     
    estimate_emission_models()
    thalamus = estimate_new_atlas()
     

     
