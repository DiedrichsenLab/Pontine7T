import numpy as np 
from DCBC.dcbc import compute_DCBC
import DCBC.utilities as ut
import Functional_Fusion.dataset as ds
import Functional_Fusion.atlas_map as am
import pandas as pd

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

    atlas, _ = am.get_atlas('MNISymThalamus1')
    atlas_fname = 'group_mean_thalamus_prob_map.nii.gz'
    U = atlas.read_data(wk_dir + '/' + atlas_fname)
    non_empty = U.sum(axis=1) > 0
    U = U[non_empty]

    atlas_obj, _ = am.get_atlas("MNISymThalamus1")
    coords = atlas_obj.world.T

    data_Social, _, _ = ds.get_dataset(base_dir,'Social',
                                       atlas="MNISymThalamus1",
                                       sess='ses-social',
                                       subj=None,
                                       type='CondRun')

    data_Language, _, _ = ds.get_dataset(base_dir,'Language',
                                         atlas="MNISymThalamus1",
                                         sess='ses-localizerfm',
                                         subj=None,
                                         type='CondRun')

    datasets = [
        data_Social,
        data_Language,
    ]

    dataset_names = [
        "Social",
        "Language",
    ]

    #atlases = [
     #   np.load(f"{wk_dir}/Prob_thalamus_nSocial.npy"),
      #  np.load(f"{wk_dir}/Prob_thalamus_nLanguage.npy"),
    #]

    atlases = [U]

    atlas_names = [
        "Anatomical"
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

    df.to_csv(f"{wk_dir}/subject_level_dcbc_anat.tsv", sep="\t", index=False)

    print("DCBC evaluation finished.")