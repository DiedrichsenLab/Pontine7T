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

def build_emissions(K=58,atlas='MNISymThalamus1'):

    #use K=58 

    atlas_read, _ = am.get_atlas('MNISymThalamus1')

    P = atlas_read.P

    data, info, ds_obj = ds.get_dataset(base_dir,'Nishimoto',atlas=atlas, sess='ses-01', type='CondRun')
    cond_vec = info['task_name']
    part_vec = info['run']
    X= ut.indicator(cond_vec)
    em_model = em.MixVMF(K=K,P=P,X=X,part_vec=part_vec)
    em_model.initialize(data)

    data2, info2, ds_obj2 = ds.get_dataset(base_dir,'Nishimoto',atlas=atlas, sess='ses-02', type='CondRun')
    cond_vec2 = info2['task_name']
    part_vec2 = info2['run']
    X2= ut.indicator(cond_vec2)
    em_model2 = em.MixVMF(K=K,P=P,X=X2,part_vec=part_vec2)
    em_model2.initialize(data=data2)

    data3, info3, ds_obj3 = ds.get_dataset(base_dir,'MDTB',atlas=atlas, sess='ses-s1', type='CondRun')
    cond_vec3 = info3['cond_name']
    part_vec3 = info3['run']
    X3= ut.indicator(cond_vec3)
    em_model3 = em.MixVMF(K=K,P=P,X=X3,part_vec=part_vec3)
    em_model3.initialize(data=data3)

    data4, info4, ds_obj4 = ds.get_dataset(base_dir,'MDTB',atlas=atlas, sess='ses-s2', type='CondRun')
    cond_vec4 = info4['cond_name']
    part_vec4 = info4['run']
    X4= ut.indicator(cond_vec4)
    em_model4 = em.MixVMF(K=K,P=P,X=X4,part_vec=part_vec4)
    em_model4.initialize(data=data4)

    data5, info5, ds_obj5 = ds.get_dataset(base_dir,'Social',atlas=atlas, sess='ses-social', type='CondRun')
    cond_vec5 = info5['task_name']
    part_vec5 = info5['run']
    X5= ut.indicator(cond_vec5)
    em_model5 = em.MixVMF(K=K,P=P,X=X5,part_vec=part_vec5)
    em_model5.initialize(data=data5)

    data6, info6, ds_obj6 = ds.get_dataset(base_dir,'Language',atlas=atlas, sess='ses-localizerfm', type='CondRun')
    cond_vec6 = info6['task_name']
    part_vec6 = info6['run']
    X6= ut.indicator(cond_vec6)
    em_model6 = em.MixVMF(K=K,P=P,X=X6,part_vec=part_vec6)
    em_model6.initialize(data=data6)

    return em_model, em_model2, em_model3, em_model4, em_model5, em_model6

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
 
    M = fm.FullMultiModel(ar_model, [emissions[2]])

    # Attach the data to the model - this is done for speed
    # The data is passed as a list with on element per data set

    M.initialize() 

    M, _, _, _ = M.fit_em(iter=200, tol=0.01,
        fit_arrangement=False,fit_emission=True,first_evidence=False)

    pt.save(M.emissions[0].V,f"{wk_dir}/V_Iglesias_MDTB_ses1.pt")

def estimate_new_atlas(): 
     
     atlas, _ = am.get_atlas('MNISymThalamus1')    

     ar_model = ar.ArrangeIndependent(K=47, P=atlas.P, spatial_specific=True, remove_redundancy=False)

     emissions = build_emissions(ar_model.K,'MNISymThalamus1')

     emissions[0].V = pt.from_numpy(np.load(f"{wk_dir}/V_matrices_Nishimoto/V_group_avg.pt"))
     #emissions[1].V = pt.load(f"{wk_dir}/V_Iglesias_Nishimoto_ses2.pt")
     #emissions[2].V = pt.load(f"{wk_dir}/V_Iglesias_MDTB_ses1.pt")
     #emissions[3].V = pt.load(f"{wk_dir}/V_Iglesias_MDTB_ses2.pt")
     emissions[4].V = pt.from_numpy(np.load(f"{wk_dir}/V_matrices_Social/V_group_avg.pt"))
     #emissions[5].V = pt.load(f"{wk_dir}/V_matrices_Language/V_group_avg.pt")
     emissions[5].V = pt.from_numpy(np.load(f"{wk_dir}/V_matrices_Language/V_group_avg.pt")
)

     #for em_model in emissions:
     emissions[0].set_param_list(['kappa'])
     #emissions[1].set_param_list(['kappa'])
     #emissions[2].set_param_list(['kappa'])
     #emissions[3].set_param_list(['kappa'])
     emissions[4].set_param_list(['kappa'])
     emissions[5].set_param_list(['kappa'])
    
     M= fm.FullMultiModel(ar_model, [emissions[0], emissions[4], emissions[5]])

     M.initialize()

     M, ll, _, _ = M.fit_em(iter=200, tol=0.01,
        fit_arrangement=True,fit_emission=True,first_evidence=True)
     
     Prob = M.arrange.marginal_prob().numpy()
     np.save(f"{wk_dir}/Prob_thalamus_first-stepV_nMDTB.npy",Prob)

     return M


def estimate_new_atlas_de_novo(): 
          
     atlas, _ = am.get_atlas('MNISymThalamus1')    

     ar_model = ar.ArrangeIndependent(K=47, P=atlas.P, spatial_specific=True, remove_redundancy=False)

     emissions = build_emissions(ar_model.K,'MNISymThalamus1')
    
     M= fm.FullMultiModel(ar_model, [emissions[0], emissions[1], emissions[2], emissions[3], emissions[4], emissions[5]])

     M.initialize()

     M, ll, _, _ = M.fit_em(iter=200, tol=0.01,
        fit_arrangement=True,fit_emission=True,first_evidence=True)
     
     Prob = M.arrange.marginal_prob().numpy()
     np.save(f"{wk_dir}/Prob_5_datasets_de_novo.npy",Prob)

     return M



if __name__ == '__main__':    
    
    estimate_new_atlas()

     # Load probability

    pmap = np.load(f"{wk_dir}/Prob_thalamus_MDTBses2.npy")
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
     

     
