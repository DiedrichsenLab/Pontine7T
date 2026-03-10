import numpy as np 
from DCBC.dcbc import compute_DCBC
import DCBC.utilities as ut
import Functional_Fusion.dataset as ds
import Functional_Fusion.atlas_map as am


wk_dir = '/Users/incehusain/fs_projects'
base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion_new'

def evaluate_atlases_with_dcbc(atlases, datasets, coords,
                               maxDist=40, binWidth=1.5,
                               backend='numpy'):

    # compute voxel distance matrix once
    dist = ut.compute_dist(coords, backend=backend)

    signs =  np.sign(coords[:,0])  #extracts sign of x-coord voxels (right is positive, left is negative)
    
    cross_hemisph_voxels = signs[:, None] != signs[None, :]  #True where signs are opposite 
    
    dist[cross_hemisph_voxels]=100 #changing cross-correlated voxel distances to 100 

    results = []

    for i in range(len(atlases)):

        atlas = atlases[i]       # probabilistic atlas
        func = datasets[i]       # left-out dataset

        # convert probabilistic atlas → hard parcellation
        parcellation = np.argmax(atlas, axis=0) + 1 

        D = compute_DCBC(
            maxDist=maxDist,
            binWidth=binWidth,
            parcellation=parcellation,
            func=func,
            dist=dist,
            weighting=True,
            backend=backend
        )

        results.append(D["DCBC"])

    return results

if __name__ == "__main__":

    atlas_obj, _ = am.get_atlas("MNISymThalamus1")
    coords = atlas_obj.world.T

    data_mdtb1, info_mdtb1, ds_obj_mdtb1 = ds.get_dataset(base_dir,'MDTB',atlas="MNISymThalamus1",sess='ses-s1', subj=None, 
                                type='CondRun')
    data_mdtb1_avg = np.mean(data_mdtb1, axis=0)
    dataset_MDTB_ses1 = data_mdtb1_avg.T

    data_mdtb2, info_mdtb2, ds_obj_mdtb2 = ds.get_dataset(base_dir,'MDTB',atlas="MNISymThalamus1",sess='ses-s2', subj=None, 
                                type='CondRun')
    data_mdtb2_avg = np.mean(data_mdtb2, axis=0)
    dataset_MDTB_ses2 = data_mdtb2_avg.T

    data_Social, info_Social, ds_obj_Social = ds.get_dataset(base_dir,'Social',atlas="MNISymThalamus1",sess='ses-social', subj=None,
                                type='CondRun')
    data_Social_avg = np.mean(data_Social, axis=0)
    dataset_Social = data_Social_avg.T

    data_Language, info_Language, ds_obj_Language = ds.get_dataset(base_dir,'Language',atlas="MNISymThalamus1",sess='ses-localizerfm', subj=None,
                                type='CondRun')
    data_Language_avg = np.mean(data_Language, axis=0)
    dataset_Language = data_Language_avg.T

    data_Nishimoto, info_Nishimoto, ds_obj_Nishimoto = ds.get_dataset(base_dir,'Nishimoto',atlas="MNISymThalamus1",sess='ses-01', subj=None,
                                type='CondRun')
    data_Nishimoto_avg = np.mean(data_Nishimoto, axis=0)
    dataset_Nishimoto_ses1 = data_Nishimoto_avg.T

    
    atlases = [np.load(f"{wk_dir}/Prob_thalamus_nMDTB-ses1.npy"),
               np.load(f"{wk_dir}/Prob_thalamus_nMDTB-ses2.npy"),
               np.load(f"{wk_dir}/Prob_thalamus_nSocial.npy"),
               np.load(f"{wk_dir}/Prob_thalamus_nLanguage.npy"),
               np.load(f"{wk_dir}/Prob_thalamus_nNishimoto-ses1.npy")]

    datasets = [dataset_MDTB_ses1, dataset_MDTB_ses2, dataset_Social, dataset_Language, dataset_Nishimoto_ses1]
    
    dcbc_results = evaluate_atlases_with_dcbc(atlases, datasets, coords)

    print("done")