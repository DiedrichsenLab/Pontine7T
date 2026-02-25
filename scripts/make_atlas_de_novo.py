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
atlas_dir = base_dir + '/Atlases/tpl-MNI152NLin2009cSymC'
wk_dir = '/Volumes/diedrichsen_data$/data/Cerebellum/Pontine7T/atlases/thalamus'

def build_emission_mdtb(K,P,atlas='MNISymThalamus1'):
    data, info,ds_obj = ds.get_dataset(base_dir,'MDTB',atlas=atlas,sess='ses-s1',  subj=None, type='CondRun')
    cond_v = info['cond_num_uni']
    part_v = info['run']

    # Make a design matrix
    X= ut.indicator(cond_v)

    # Build an emission model
    em_model = em.MixVMF(K=K,P=P, X=X,part_vec=part_v)
    em_model.initialize(data)
    return em_model

def build_emission_mdtb_ses2(K,P,atlas='MNISymThalamus1'):
    data, info,ds_obj = ds.get_dataset(base_dir,'MDTB',atlas=atlas,sess='ses-s2',subj=None, type='CondRun')

    cond_v = info['cond_num_uni']
    part_v = info['run']

    # Make a design matrix
    X= ut.indicator(cond_v)
    # Build an emission model
    em_model = em.MixVMF(K=K,P=P, X=X,part_vec=part_v)
    em_model.initialize(data)
    return em_model

def build_emission_language(K,P,atlas='MNISymThalamus1'):

    #data is in form: (subjects, run x tasks, voxels)
    data, info, ds_obj = ds.get_dataset(base_dir,'Language',atlas=atlas, sess='ses-localizer_cond_fm', subj=None, type='CondRun')
    cond_v = info['task']
    part_v = info['run']

    # Make a design matrix
    X= ut.indicator(cond_v)

    data = ds.remove_baseline(data,part_v)

    # Build an emission model
    em_model = em.MixVMF(K=K,P=P, X=X,part_vec=part_v)
    em_model.initialize(data)
    return em_model


def build_emission_pontine(K,P,atlas='MNISymThalamus1'):
    data, info, ds_obj = ds.get_dataset(base_dir,'Pontine',atlas=atlas, sess='ses-s1', subj=None, type='CondRun')
    cond_v = info['task']
    part_v = info['run']

    # Make a design matrix
    X= ut.indicator(cond_v)

    data = ds.remove_baseline(data,part_v)

    # Build an emission model
    em_model = em.MixVMF(K=K,P=P, X=X,part_vec=part_v)
    em_model.initialize(data)
    return em_model

#add more build_emission functions for remaining datasets: HCP, Social, Somatotopic, Working Memory, IBC, Demand, Nishimito 

def estimate_emission_models():
    """ Estimate emission models for dataset
    """
    # Get the atlas for the thalamus
    atlas, _ = am.get_atlas('MNISymThalamus1')

    # Load the probabilistic atlas into that space - this will be replaced with the Iglesias thalamus atlas if it is of form "KxP" 

    atlas_fname = 'atl-NettekovenSym32_space-MNI152NLin2009cSymC_probseg.nii'

    U = atlas.read_data(atlas_dir + '/' + atlas_fname)
    
    # Build the arrangement model 
    ar_model = ar.build_arrangement_model(U, prior_type='prob', atlas=atlas)
        
    # Load the data using functional fusion for all 10 datasets

    em_model = build_emission_mdtb(ar_model.K,atlas.P,atlas='MNISymThalamus1')

    em_model1 = build_emission_mdtb_ses2(ar_model.K,atlas.P,atlas='MNISymThalamus1')

    em_model2 = build_emission_language(ar_model.K,atlas.P,atlas='MNISymThalamus1')

    em_model3  = build_emission_pontine(ar_model.K,atlas.P,atlas='MNISymThalamus1')
 
    # Build the full model: inlcude all 10 datasets

    M = fm.FullMultiModel(ar_model, [em_model, em_model1, em_model2, em_model3])

    # Attach the data to the model - this is done for speed
    # The data is passed as a list with on element per data set

    M.initialize() # Initialize the default that each em has different subjects

    # Now we can run the EM algorithm
    M, _, _, _ = M.fit_em(iter=200, tol=0.01,
        fit_arrangement=False,fit_emission=True,first_evidence=False)
    
    #Prob = M.arrange.marginal_prob().numpy()
    
    #np.save(f"{wk_dir}/Prob_cereb_grey(mdtb_ses1).npy",Prob)

    pt.save(M.emissions[0].V,f"{wk_dir}/V_Iglesias_MDTB(ses2).pt")
    pt.save(M.emissions[1].V,f"{wk_dir}/V_Iglesias_Language.pt")
    pt.save(M.emissions[0].V,f"{wk_dir}/V_Iglesias_Pontine(set1).pt")

def estimate_new_atlas(): 

     atlas, _ = am.get_atlas('MNISymThalamus1')    

     ar_model = ar.ArrangeIndependent(K=23, P=atlas.P, spatial_specific=True, remove_redundancy=False)
    
    #add all 10 datasets
     
     em_model = build_emission_mdtb_ses2(ar_model.K,atlas.P,atlas='MNISymThalamus1')
     em_model1 = build_emission_mdtb(ar_model.K,atlas.P,atlas='MNISymThalamus1')
     em_model2 = build_emission_language(ar_model.K,atlas.P,atlas='MNISymThalamus1')
     em_model3  = build_emission_pontine(ar_model.K,atlas.P,atlas='MNISymThalamus1')

     em_model.V = pt.load(f"{wk_dir}/V_Iglesias_MDTB(ses2).pt")
     em_model1.V = pt.load(f"{wk_dir}/V_Iglesias_MDTB(ses1).pt")
     em_model2.V = pt.load(f"{wk_dir}/V_Iglesias_Language.pt")
     em_model3.V = pt.load(f"{wk_dir}/V_Iglesias_Pontine.pt")

     em_model.set_param_list(['kappa'])
     em_model1.set_param_list(['kappa'])
     em_model2.set_param_list(['kappa'])
     em_model3.set_param_list(['kappa'])

     M= fm.FullMultiModel(ar_model, [em_model, em_model1, em_model2, em_model3])

     M.initialize()

     M, ll, _, _ = M.fit_em(iter=200, tol=0.01,
        fit_arrangement=True,fit_emission=True,first_evidence=True)
     
     Prob = M.arrange.marginal_prob().numpy()
     np.save(f"{wk_dir}/Prob_thalamus_Iglesias.npy",Prob)

     return M


if __name__ == '__main__':

    thalamus = estimate_new_atlas()

    # Load probability 
    
    pmap = np.load(f"{wk_dir}/Prob_thalamus_Iglesias.npy")

    # Load colormap and labels #is there a .lut file I can use for the Iglesias atlas? If I did: what would this mean?
    lid,cmap,names = nt.read_lut('/Volumes/diedrichsen_data$/data/FunctionalFusion/Atlases/tpl-MNI152NLin2009cSymC/atl-NettekovenSym32.lut')

    wta = np.argmax(pmap, axis=0) 
    wta += 1
    wta_int32 = wta.astype(np.int32)
    
    parcellation = plot.plot_thalamus(wta_int32,cscale=[0,47],cmap=cmap[0:48])

    #pass 
    
    #Vs = [em.V for em in M.emissions]
    #plt.imshow(Vs[0])

    #plt.yticks(info.cond_name[info.run==1].values)
    #plt.y()


    #arrangement = estimate_new_atlas(ems)
    
    # arrangement.marginal_prob() will give the probability of each parcel for each voxel in the new atlas. (KxP)