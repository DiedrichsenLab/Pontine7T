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

base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion' 
atlas_dir = base_dir + '/Atlases/tpl-MNI152NLin2009cSymC'
wk_dir = '/Volumes/diedrichsen_data$/data/Cerebellum/Pontine7T/dentate_atlas'

def build_emission_mdtb(K,P,atlas='MNISymCereb2'):
    data, info,ds_obj = ds.get_dataset(base_dir,'MDTB',atlas=atlas,type='CondRun',sess='ses-s1',subj=None)
    cond_v = info['cond_num_uni']
    part_v = info['run']

    # Make a design matrix
    X= ut.indicator(cond_v)
    # Build an emission model
    em_model = em.MixVMF(K=K,P=P, X=X,part_vec=part_v)
    em_model.initialize(data)
    return em_model

def build_emission_language(K,P,atlas='MNISymCereb2'):

    #data is in form: (subjects, run x tasks, voxels)
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


def build_emission_pontine(K,P,atlas='MNISymCereb1'):
    data_dir = '/Volumes/diedrichsen_data$/data/Cerebellum/Pontine7T/RegionOfInterest_BOLDMNI/data/group'
    if atlas=='MNISymCereb1': 
        at,_ = am.get_atlas('MNISymCereb1')
        data_pontine = dv.get_structure_data(structure='cereb_gray', data_dir=data_dir )
        
    elif atlas == 'MNISymDentate1':
        at,_ = am.get_atlas('MNISymDentate1')
        data_pontine = dv.get_structure_data(structure='dentate',  data_dir=data_dir)

    elif atlas == 'MNISymOlive11':
        at,_ = am.get_atlas('MNISymOlive1')
        data_pontine = dv.get_structure_data(structure='olive',  data_dir=data_dir)

    elif atlas == 'MNISymRedNucleus1':
        at,_ = am.get_atlas('MNISymRedNucleus1')
        data_pontine = dv.get_structure_data(structure='rednucleus',  data_dir=data_dir)

    elif atlas == 'MNISymPontine1':
        at,_ = am.get_atlas('MNISymPontine1')
        data_pontine = dv.get_structure_data(structure='pontine',  data_dir=data_dir)

    elif atlas == 'MNISymThalamus1':
        at,_ = am.get_atlas('MNISymThalamus1')
        data_pontine = dv.get_structure_data(structure='thalamus',  data_dir=data_dir)
 
    cifti_objects = data_pontine[4]
    data = []

    for cifti in cifti_objects:
        subject_data = at.read_data(cifti).T
        data.append(subject_data)

    data = np.stack(data)

    cond_v = np.tile(np.arange(1,11),16)

    part_v = np.repeat(np.arange(1,17), 10)

    data = ds.remove_baseline(data,part_v)

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

    #I had to remove the "T" for it to work 

    U = atlas.read_data(atlas_dir + '/' + atlas_fname)
    
    # Build the arrangement model 
    ar_model = ar.build_arrangement_model(U, prior_type='prob', atlas=atlas)
        
    # Load the data using functional fusion 

    em_model1 = build_emission_mdtb(ar_model.K,atlas.P,atlas='MNISymCereb2')

    em_model2 = build_emission_language(ar_model.K,atlas.P,atlas='MNISymCereb2')

   # em_model3  = build_emission_pontine(ar_model.K,atlas.P,atlas='MNISymCereb1')
 
    # Build the full model: The emission models are passed as a list, as usually we have multiple data sets
    M = fm.FullMultiModel(ar_model, [em_model1])

    # Attach the data to the model - this is done for speed
    # The data is passed as a list with on element per data set
    M.initialize() # Initialize the default that each em has different subjects

    # Now we can run the EM algorithm
    M, _, _, _ = M.fit_em(iter=200, tol=0.01,
        fit_arrangement=False,fit_emission=True,first_evidence=False)
    
    Prob = M.arrange.marginal_prob().numpy()
    
    np.save(f"{wk_dir}/Prob_cereb_grey_lang+mdtb(ses1).npy",Prob)

    pt.save(M.emissions[0].V,f"{wk_dir}/V_cerebcortex_MDTB(ses1).pt")
    #pt.save(M.emissions[1].V,f"{wk_dir}/V_cerebcortex_res2_Language.pt")
   # pt.save(M.emissions[2].V,f"{wk_dir}/V_cerebcortex_Pontine.pt")

def estimate_new_atlas(): 
     
    # Get atlas for dentate 

     atlas, _ = am.get_atlas('MNISymDentate1')    

     ar_model = ar.ArrangeIndependent(K=32, P=atlas.P, spatial_specific=True, remove_redundancy=False)
         
     em_model1 = build_emission_mdtb(ar_model.K,atlas.P,atlas='MNISymDentate1')
     #em_model2 = build_emission_language(ar_model.K,atlas.P,atlas='MNISymDentate1')
     #em_model3  = build_emission_pontine(ar_model.K,atlas.P,atlas='MNISymDentate1')

     em_model1.V = pt.load(f"{wk_dir}/V_cerebcortex_MDTB(ses1).pt")
     #em_model2.V = pt.load(f"{wk_dir}/V_cerebcortex_Language.pt")
     #em_model3.V = pt.load(f"{wk_dir}/V_cerebcortex_Pontine.pt")

     em_model1.set_param_list(['kappa'])
     #em_model2.set_param_list(['kappa'])
     #em_model3.set_param_list(['kappa'])

     M= fm.FullMultiModel(ar_model, [em_model1])

     M.initialize()

     M, ll, _, _ = M.fit_em(iter=200, tol=0.01,
        fit_arrangement=True,fit_emission=True,first_evidence=True)
     
     Prob = M.arrange.marginal_prob().numpy()
     np.save(f"{wk_dir}/Prob_dentate_mdtb(ses1).npy",Prob)

     return M


if __name__ == '__main__':

    #lang_data = build_emission_pontine(32,18290,'MNISymCereb1')

    #all_Vs = estimate_emission_models()

    #dentateM = estimate_new_atlas()

    # Load probability 
    pmap = np.load(f"{wk_dir}/Prob_dentate_language.npy")
    # Load colormap and labels
    lid,cmap,names = nt.read_lut('/Volumes/diedrichsen_data$/data/FunctionalFusion/Atlases/tpl-MNI152NLin2009cSymC/atl-NettekovenSym32.lut')


    #wta = np.argmax(data,axis=0)
    wta = np.argmax(pmap, axis=0) 
    wta += 1

    wta_int32 = wta.astype(np.int32)
    
    dentate_parcellation = plot.plot_dentate(wta_int32,cscale=[0,32],cmap=cmap)

    #pass 
    
    Vs = [em.V for em in M.emissions]
    #plt.imshow(Vs[0])

    #plt.yticks(info.cond_name[info.run==1].values)
    #plt.y()


    #arrangement = estimate_new_atlas(ems)
    
    # arrangement.marginal_prob() will give the probability of each parcel for each voxel in the new atlas. (KxP)