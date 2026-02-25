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
import plot_thalamus 

base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion_new' 
wk_dir = '/Users/incehusain/fs_projects'

def build_emissions(K,atlas='MNISymThalamus1'):

    atlas_read, _ = am.get_atlas('MNISymThalamus1')

    P = atlas_read.P

    data, info, ds_obj = ds.get_dataset(base_dir,'Social',atlas=atlas, sess='ses-social', type='CondRun')
    cond_vec = info['task_name']
    part_vec = info['run']
    X= ut.indicator(cond_vec)
    em_model = em.MixVMF(K=K,P=P,X=X,part_vec=part_vec)
    em_model.initialize(data)

    return em_model

def estimate_emission_models():
    """ Estimate emission models for dataset
    """
    atlas, _ = am.get_atlas('MNISymThalamus1')

    atlas_fname = 'group_mean_thalamus_prob_map.nii.gz'

    U = atlas.read_data(wk_dir + '/' + atlas_fname)

    non_empty = U.sum(axis=1) > 0
    U = U[non_empty]
    
    ar_model = ar.build_arrangement_model(U, prior_type='prob', atlas=atlas)
            
    emissions = build_emissions(ar_model.K,'MNISymThalamus1')
 
    M = fm.FullMultiModel(ar_model, [emissions])

    # Attach the data to the model - this is done for speed
    # The data is passed as a list with on element per data set

    M.initialize() 

    M, _, _, _ = M.fit_em(iter=200, tol=0.01,
        fit_arrangement=False,fit_emission=True,first_evidence=False)

    pt.save(M.emissions[0].V,f"{wk_dir}/V_Iglesias_Social.pt")

def estimate_new_atlas(): 
     
     atlas, _ = am.get_atlas('MNISymThalamus1')    

     ar_model = ar.ArrangeIndependent(K=47, P=atlas.P, spatial_specific=True, remove_redundancy=False)

     emissions = build_emissions(ar_model.K,'MNISymThalamus1')

     emissions.V = pt.load(f"{wk_dir}/V_Iglesias_Social.pt")

     #for em_model in emissions:
     emissions.set_param_list(['kappa'])
    
     M= fm.FullMultiModel(ar_model, [emissions])

     M.initialize()

     M, ll, _, _ = M.fit_em(iter=200, tol=0.01,
        fit_arrangement=True,fit_emission=True,first_evidence=True)
     
     Prob = M.arrange.marginal_prob().numpy()
     np.save(f"{wk_dir}/Prob_thalamus_Social.npy",Prob)

     return M

if __name__ == '__main__':    

    pmap = np.load(f"{wk_dir}/Prob_thalamus_Social.npy")
    pmap2 = np.load(f"{wk_dir}/Prob_Language_thalamus.npy")

    # Load colormap and labels #is there a .lut file I can use for the Iglesias atlas? If I did: what would this mean?
    lid,cmap_lut,names = nt.read_lut(f'{wk_dir}/thalamus_atlas.lut')

    wta = np.argmax(pmap, axis=0) 
    wta += 1
    wta_int32 = wta.astype(np.int32)
    
    fig, axes = plot_thalamus.plot_thalamus_wta(
    wta_data=wta_int32,
    bg_img=None,       # or your structural Nifti
    lut_cmap=cmap_lut,
    z_coords=[-7,-5,-3,1,8,12,15])
    
    plt.show()
     

     
