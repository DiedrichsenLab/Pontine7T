import ProbabilisticParcellation.util as ut
import Functional_Fusion.dataset as ds
import matplotlib.pyplot as plt
import seaborn as sb
import pandas as pd
from scipy import stats
import glob
import numpy as np
import os
import time
import torch as pt
import Functional_Fusion.atlas_map as am



def calc_test_dcbc(
    parcels, testdata, dist, max_dist=110, bin_width=5, trim_nan=False, verbose=True
):
    """DCBC: evaluate the resultant parcellation using DCBC
    Args:
        parcels (np.ndarray): the input parcellation:
            either group parcellation (1-dimensional: P)
            individual parcellation (num_subj x P )
        testdata (np.ndarray): the functional test dataset,
                                shape (num_sub, N, P)
        dist (<AtlasVolumetric>): the class object of atlas
        max_dist (int): the maximum distance between voxel pairs
        bin_width (int): the bin width of distance
        trim_nan (boolean): if true, the nan voxel label will be
                            removed from DCBC calculation. Otherwise,
                            we treat nan voxels are in the same parcel
                            which is label 0 by default.
    Returns:
        dcbc_values (np.ndarray): the DCBC values of subjects
    """

    #
    # if trim_nan:  # mask the nan voxel pairs distance to nan
    #     dist[np.where(np.isnan(parcels))[0], :] = np.nan
    #     dist[:, np.where(np.isnan(parcels))[0]] = np.nan

    dcbc_values = []
    for sub in range(testdata.shape[0]):
        if verbose:
            print(f"Subject {sub}", end=":")
        tic = time.perf_counter()
        if parcels.ndim == 1:
            D = ut.compute_DCBC(
                maxDist=max_dist,
                binWidth=bin_width,
                parcellation=parcels,
                dist=dist,
                func=testdata[sub].T,
            )
        else:
            D = ut.compute_DCBC(
                maxDist=max_dist,
                binWidth=bin_width,
                parcellation=parcels[sub],
                dist=dist,
                func=testdata[sub].T,
            )
        dcbc_values.append(D["DCBC"])
        toc = time.perf_counter()
        if verbose:
            print(f"{toc-tic:0.4f}s")
    return pt.stack(dcbc_values)


if __name__ == "__main__":

    atlas, _ = am.get_atlas('MNISymThalamus1')

    base_dir = "/Volumes/diedrichsen_data$/data/Cerebellum/Pontine7T/atlases/thalamus/indiv_parcellations/MDTB-ses1"

    indiv_parcellation = np.load(f"{base_dir}/all_indiv_parcels.npy")
    indiv_parcellation_torch = pt.tensor(indiv_parcellation, dtype=pt.int32)

    dataset = np.load(f"{base_dir}/MDTB_ses1_data.npy")
    dataset_torch = pt.tensor(dataset, dtype=pt.float32)

    dcbc_values = calc_test_dcbc(parcels = indiv_parcellation_torch, 
                                 testdata = dataset_torch, dist = atlas, 
                                 max_dist=110, bin_width=5, trim_nan=False, verbose=True)

    dcbc_values_npy = dcbc_values.cpu().numpy() if pt.is_tensor(dcbc_values) else dcbc_values

    np.save(f"{base_dir}/dcbc.npy", dcbc_values.cpu().numpy())
    
    """DCBC: evaluate the resultant parcellation using DCBC
    Args:
        parcels (np.ndarray): the input parcellation:
            either group parcellation (1-dimensional: P)
            individual parcellation (num_subj x P )
        dist (<AtlasVolumetric>): the class object of atlas
        testdata (np.ndarray): the functional test dataset,
                                shape (num_sub, N, P)
        trim_nan (boolean): if true, the nan voxel label will be
                            removed from DCBC calculation. Otherwise,
                            we treat nan voxels are in the same parcel
                            which is label 0 by default.
    Returns:
        dcbc_values (np.ndarray): the DCBC values of subjects
    """


