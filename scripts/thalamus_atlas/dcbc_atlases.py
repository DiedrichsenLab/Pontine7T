import numpy as np 
from DCBC.dcbc import compute_DCBC
import DCBC.utilities as ut
import Functional_Fusion.dataset as ds
import Functional_Fusion.atlas_map as am
import pandas as pd
from scripts import recenter_data as rd

wk_dir = '/Users/incehusain/fs_projects'
base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion_new'

def evaluate_atlases_subject_dcbc(atlases, atlas_names,
                                  datasets, dataset_names,
                                  coords,
                                  maxDist=40, binWidth=1.5,
                                  backend='numpy'):

    # Compute voxel distance matrix once
    dist = ut.compute_dist(coords, backend=backend)

    signs = np.sign(coords[:,0])
    cross_hemisph_voxels = signs[:, None] != signs[None, :]
    dist[cross_hemisph_voxels] = 100

    rows = []

    for a, atlas in enumerate(atlases):

        atlas_name = atlas_names[a]

        # winner-take-all parcellation
        parcellation = np.argmax(atlas, axis=0) + 1

        for d, data in enumerate(datasets):

            dataset_name = dataset_names[d]

            n_subj = data.shape[0]

            for s in range(n_subj):

                # subject functional data
                func = data[s].T   
                
                D = compute_DCBC(
                    maxDist=maxDist,
                    binWidth=binWidth,
                    parcellation=parcellation,
                    func=func,
                    dist=dist,
                    weighting=True,
                    backend=backend
                )

                rows.append({
                    "Test dataset": dataset_name,
                    "Test subject": s + 1,
                    "Atlas": atlas_name,
                    "DCBC": D["DCBC"]
                })

    df = pd.DataFrame(rows)

    return df

if __name__ == "__main__":

    #atlas, _ = am.get_atlas('MNISymThalamus1')
    #atlas_fname = 'group_mean_thalamus_prob_map.nii.gz'
    #U = atlas.read_data(wk_dir + '/' + atlas_fname)
    #non_empty = U.sum(axis=1) > 0
    #U = U[non_empty]

    atlas_obj, _ = am.get_atlas("MNISymThalamus1")
    coords = atlas_obj.world.T

    data_ses1, info_ses1, ds_obj_ses1 = ds.get_dataset(base_dir, "Nishimoto" ,atlas='MNISymThalamus1', 
                                            sess='ses-01', 
                                            subj=None, 
                                            type='CondRun')
    
    data_ses2, info_ses2, ds_obj_ses2 = ds.get_dataset(base_dir, "Nishimoto" ,atlas='MNISymThalamus1', 
                                            sess="ses-02", 
                                            subj=None, 
                                            type='CondRun')
    
    info_ses1_path = f'{base_dir}/Nishimoto/derivatives/ffextract/sub-02/sub-02_ses-01_CondRun.tsv'
    info_ses1 = pd.read_csv(info_ses1_path, sep='\t')

    info_ses2_path = f'{base_dir}/Nishimoto/derivatives/ffextract/sub-02/sub-02_ses-02_CondRun.tsv'
    info_ses2 = pd.read_csv(info_ses2_path, sep='\t')   

    recentered_ses1 = rd.recenter_data(data_ses1, info_ses1, center_full_code='rest_close', keep_center=True)

    recentered_ses2 = rd.recenter_data(data_ses2, info_ses2, center_full_code='rest_close', keep_center=True)

    merged_data, merged_info = rd.merge_sessions([recentered_ses1[0], recentered_ses2[0]], [recentered_ses1[1], recentered_ses2[1]])


    #data_Social, _, _ = ds.get_dataset(base_dir,'Social',
     #                                  atlas="MNISymThalamus1",
      #                                 sess='ses-social',
       #                                subj=None,
        #                               type='CondRun')

    #data_Language, _, _ = ds.get_dataset(base_dir,'Language',
     #                                    atlas="MNISymThalamus1",
      #                                   sess='ses-localizerfm',
       #                                  subj=None,
        #                                 type='CondRun')

    datasets = [
        merged_data
    ]

    dataset_names = [
        "Nishimoto"
    ]

    atlases = [
        np.load(f"{wk_dir}/Prob_thalamus_nNishimoto.npy"),
    ]

    #atlases = [U]

    atlas_names = [
        "Atlas_no_Nishimoto"
    ]

    df = evaluate_atlases_subject_dcbc(
        atlases,
        atlas_names,
        datasets,
        dataset_names,
        coords
    )

    #df.to_csv(f"{wk_dir}/subject_level_dcbc2.tsv",
     #     sep="\t",
      #    mode="a",
       #   header=False,
        #  index=False)

    df.to_csv(f"{wk_dir}/subject_level_dcbc_Nishimoto.tsv", sep="\t", index=False)

    print("DCBC evaluation finished.")