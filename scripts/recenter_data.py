import numpy as np
import Functional_Fusion.dataset as ds
import pandas as pd


def recenter_data(data, info, center_full_code='rest_task', keep_center=False):
    """
    Recenter  data by subtracting the center condition from each partition.
    Works with CondAll, CondHalf, and CondRun data types.
    
    Parameters:
        data (np.ndarray):  data of shape (n_subjects, n_conditions x repititions, n_voxels)
        info (pd.DataFrame): condition information with columns like 'task_code', 'cond_code', 'half', 'run'
        center_full_code (str): full task condition code to center around as 'taskcode_condcode' (default: 'rest_task')
        keep_center (bool): if True, keep center condition as zeros; if False, remove it
    
    Returns:
        processed_data (np.ndarray): recentered fMRI data
        updated_info (pd.DataFrame): updated info
    """
    n_subjects, _, _ = data.shape
    
    # Split on first underscore only
    center_task, center_cond = center_full_code.split('_', 1)
    
    # Determine partition column
    if 'half' in info.columns:
        partition_col = 'half'
        partitions = info['half'].unique()
    elif 'run' in info.columns:
        partition_col = 'run'
        partitions = info['run'].unique()
    else:
        partition_col = None
        partitions = [None]
    
    # Find center condition indices (match both task_code and cond_code)
    center_mask = (info['task_code'] == center_task) & (info['cond_code'] == center_cond)
    
    if not center_mask.any():
        raise ValueError(f"Center condition with task_code='{center_task}' and cond_code='{center_cond}' not found in data")
    
    # Process each subject
    processed_data_list = []
    for subject_idx in range(n_subjects):
        subject_data = data[subject_idx]
        subject_recentered = []
        
        for partition in partitions:
            # Get partition mask
            if partition_col is not None:
                partition_mask = (info[partition_col] == partition)
            else:
                partition_mask = np.ones(len(info), dtype=bool)
            
            # Find center condition for this partition
            center_idx_partition = np.where(partition_mask & center_mask)[0]
            
            if len(center_idx_partition) == 0:
                raise ValueError(f"No center condition '{center_full_code}' found for partition {partition}")
            
            # Get center condition data
            center_data = subject_data[center_idx_partition[0]]
            
            if keep_center:
                # Keep all conditions including center (as zeros)
                partition_indices = np.where(partition_mask)[0]
                recentered = subject_data[partition_indices] - center_data
            else:
                # Remove center condition
                task_mask = partition_mask & ~center_mask
                partition_indices = np.where(task_mask)[0]
                recentered = subject_data[partition_indices] - center_data
            
            subject_recentered.append(recentered)
        
        processed_data_list.append(np.vstack(subject_recentered))
    
    processed_data = np.stack(processed_data_list, axis=0)
    
    # Update info
    if keep_center:
        updated_info = info.copy()
    else:
        updated_info = info[~center_mask].reset_index(drop=True)
    
    return processed_data, updated_info


def merge_sessions(data_list, info_list):
    
    """
    Merge multiple sessions (already recentered).
    
    Parameters
    ----------
    data_list : list of np.ndarray
    Each array has shape (n_subjects, n_conditions, n_voxels)
    info_list : list of pd.DataFrame
    Matching info DataFrames for each session
    Must contain 'task_code' and 'cond_code'
    
    Returns
    -------
    merged_data : np.ndarray
    Shape (n_subjects, merged_conditions, n_voxels)
    merged_info : pd.DataFrame
    Info corresponding to merged_data
    """

    #Build condition codes for matching ---
    
    full_codes_list = []

    for info in info_list:
        full_codes = info['task_code'] + '_' + info['cond_code']
        full_codes_list.append(full_codes.values)

    # Collect all unique condition codes 
    
    all_codes = sorted(set(np.concatenate(full_codes_list)))

    merged_data_dict = {}
    merged_info_rows = []

    for code in all_codes:
        condition_data = []
        info_row = None

        for data, info, codes in zip(data_list, info_list, full_codes_list):
            matches = np.where(codes == code)[0]

            if len(matches) > 0:
                idx = matches[0]
                condition_data.append(data[:, idx, :])
                info_row = info.iloc[idx]

        stacked = np.stack(condition_data, axis=0)
        averaged = stacked.mean(axis=0)
        
        merged_data_dict[code] = averaged
        merged_info_rows.append(info_row)    

    merged_data = np.stack([merged_data_dict[c] for c in all_codes],axis=1)

    merged_info = pd.DataFrame(merged_info_rows).reset_index(drop=True)

    return merged_data, merged_info 

if __name__ == '__main__':

    #center_full_code = 'rest_close', 'rest_task'

    dataset= 'WMFS'
    sessions = ['ses-01', 'ses-02']

    base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion_new' 

    #Nishimoto: sub_list = ['sub-01','sub-02','sub-03', 'sub-04', 'sub-05','sub-06']
    
    sub_list = ['sub-01','sub-02','sub-03', 'sub-04', 'sub-05','sub-06','sub-07','sub-08','sub-09', 'sub-10','sub-11','sub-12','sub-13','sub-14','sub-15','sub-16']

    #MDTB: sub_list = ['sub-02','sub-03','sub-04', 'sub-06','sub-08','sub-09', 'sub-10', 'sub-12','sub-14','sub-15','sub-17', 'sub-18',
     #         'sub-19','sub-20', 'sub-21', 'sub-22','sub-24', 'sub-25','sub-26','sub-27', 'sub-28', 'sub-29','sub-30', 'sub-31']
        
    data_ses1, info_ses1, ds_obj_ses1 = ds.get_dataset(base_dir, dataset ,atlas='MNISymThalamus1', 
                                            sess=sessions[0], 
                                            subj=sub_list, 
                                            type='CondAll')
    
    data_ses2, info_ses2, ds_obj_ses2 = ds.get_dataset(base_dir, dataset ,atlas='MNISymThalamus1', 
                                            sess=sessions[1], 
                                            subj=sub_list, 
                                            type='CondAll')
    
    info_ses1_path = f'{base_dir}/{dataset}/derivatives/ffextract/sub-02/sub-02_{sessions[0]}_CondAll.tsv'
    info_ses1 = pd.read_csv(info_ses1_path, sep='\t')

    info_ses2_path = f'{base_dir}/{dataset}/derivatives/ffextract/sub-02/sub-02_{sessions[1]}_CondAll.tsv'
    info_ses2 = pd.read_csv(info_ses2_path, sep='\t')   

    #recentered_ses1 = recenter_data(data_ses1, info_ses1, center_full_code='rest_close', keep_center=True)

    #recentered_ses2 = recenter_data(data_ses2, info_ses2, center_full_code='rest_close', keep_center=True)

    #merged_data, merged_info = merge_sessions([recentered_ses1[0], recentered_ses2[0]], [recentered_ses1[1], recentered_ses2[1]])

    merged_data, merged_info = merge_sessions([data_ses1, data_ses2], [info_ses1, info_ses2])   
    
    baseline_all = np.ones(merged_data.shape[1])

    mean_centered = ds.remove_baseline(merged_data, baseline_all)

    print("DONE!")


