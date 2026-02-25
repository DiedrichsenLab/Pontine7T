import numpy as np 
import torch as pt

wk_dir = '/Users/incehusain/fs_projects'
base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion_new' 

def calc_cosine_similarity_btwn(subj = ['sub-01', 'sub-02'], type='group', dataset='Language', sess = 'ses-localizerfm'):

    #computes average cosine similarities between subjects for each same parcel 

    avg_similarity_per_parcel = []
    similarity_distributions_per_parcel = []
    
    for k in range(58):
        parcel_vectors_subj = []

        for sub in subj:
            V_indiv = pt.load(f"{wk_dir}/V_matrices_{dataset}/{sess}_V_{type}_{sub}_norm.pt")
            parcel_vector = V_indiv[:, k]
            parcel_vectors_subj.append(parcel_vector)
            
            X = np.stack(parcel_vectors_subj, axis=0)

        cosine_similarity_matrix = X@X.T

        triangle_indices = np.triu_indices_from(cosine_similarity_matrix, k=1)
        cosine_similarity_matrix_values = cosine_similarity_matrix[triangle_indices]

        similarity_distributions_per_parcel.append(cosine_similarity_matrix_values)

        mean_similarity = np.mean(cosine_similarity_matrix_values) 

        avg_similarity_per_parcel.append(mean_similarity) 
    
    return avg_similarity_per_parcel, similarity_distributions_per_parcel   

def calc_cosine_similarity_within(subj = ['sub-01', 'sub-02'], type='group', dataset='Language', sess = 'ses-localizerfm'):

#calculates cosine similarities within each subject between all parcel pairs
    
    avg_similarity_per_subject = []
    subj_similarity_matrix = {}
    
    for sub in subj:
        V = pt.load(f"{wk_dir}/V_matrices_{dataset}/V_{sess}_{sub}_norm.pt")
        X = V.T
        
        cosine_similarity_matrix = X@X.T

        triangle_indices = np.triu_indices_from(cosine_similarity_matrix, k=1)
        cosine_similarity_matrix_values = cosine_similarity_matrix[triangle_indices]

        mean_similarity = np.mean(cosine_similarity_matrix_values) 

        avg_similarity_per_subject.append(mean_similarity) 
        subj_similarity_matrix[sub] = cosine_similarity_matrix

    return avg_similarity_per_subject, subj_similarity_matrix


def calc_parcelwise_cosine_reliability(subj=['sub-01', 'sub-02'],
                                       dataset='Language',
                                       sess='ses-localizerfm'
                                       ):
    """
    Computes parcel-wise cross-run cosine reliability per subject.

    For each subject:
        - Split trials into first half and second half
        - Compute parcel x parcel cosine similarity matrix for each half
        - For each parcel, compute mean similarity to all other parcels (excluding self)
        - Compute cosine similarity between the two mean-similarity vectors (parcel-wise reliability)

    Returns:
        parcelwise_reliability: dict of subject -> array of length n_parcels
            reliability score per parcel
        S_matrices: dict of subject -> (S1, S2)
            full cosine similarity matrices per half
    """

    parcelwise_reliability = {}
    S_matrices = {}

    for sub in subj:
        # Load V matrix: (trials, parcels)
        V = pt.load(f"{wk_dir}/V_matrices_{dataset}/ses2_V_{sess}_{sub}_norm.pt")
        n_trials, n_parcels = V.shape

        half = n_trials // 2
        V1 = V[:half, :]
        V2 = V[half:, :]

        # Cosine similarity matrices: parcel x parcel
        S1 = V1.T @ V1  # parcels x parcels
        S2 = V2.T @ V2

        # Compute mean similarity for each parcel (excluding self)
        #mean_sim_S1 = (S1.sum(dim=1) - 1) / (n_parcels - 1)
        #mean_sim_S2 = (S2.sum(dim=1) - 1) / (n_parcels - 1)

        # Compute cross-run cosine similarity per parcel (reliability)
        # Each parcel has a scalar: how similar its profile is across runs
        # Note: mean similarity vectors are already 1D tensors of length n_parcels

        # But for full profile reliability, you might want parcel-to-parcel similarity vector
        # Here, we use the row of S1 and S2 (excluding diagonal) for each parcel

        parcel_reliability = pt.zeros(n_parcels)

        for i in range(n_parcels):
            # similarity profile for parcel i (exclude self)
            profile1 = np.concatenate([S1[i, :i], S1[i, i+1:]])
            profile2 = np.concatenate([S2[i, :i], S2[i, i+1:]])

            profile1 = profile1 / np.linalg.norm(profile1)
            profile2 = profile2 / np.linalg.norm(profile2)

            parcel_reliability[i] = np.dot(profile1, profile2)

        parcelwise_reliability[sub] = parcel_reliability
        S_matrices[sub] = (S1, S2)

    return parcelwise_reliability, S_matrices

def avg_similarity_matrices(subj_similarity_matrices):

    all_matrices = []
    
    for sub in subj_similarity_matrices.keys():
        all_matrices.append(subj_similarity_matrices[sub])
    
    mean_similarity_matrix = np.mean(np.array(all_matrices), axis=0)
    
    return mean_similarity_matrix