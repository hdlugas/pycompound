
# import libraries
import numpy as np
import pandas as pd
import os
import sys

# add src directory to python's module search path and import all functions from processing and similarity_measures
sys.path.append(os.path.abspath(os.path.join(os.getcwd(), '../src')))
from processing import *
from similarity_measures import *

# set seed for reproducibility
np.random.seed(1)

# display entire pandas dataframe
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)



###################################### LC-MS/MS Example ######################################
print('\n\n\nLC-MS/MS example:\n')

# import query and reference libraries
df_query = pd.read_csv(os.path.join(os.getcwd(), '../data') + '/lcms_query_library.csv')
df_reference = pd.read_csv(os.path.join(os.getcwd(), '../data') + '/lcms_reference_library.csv')

# get unique spectrum IDs
unique_query_ids = df_query.iloc[:,0].unique()
unique_reference_ids = df_reference.iloc[:,0].unique()

# compute the similarity score between each query spectrum/spectra and all reference spectra and record the predicted compound for each query along with its corresponding similarity score
preds_cosine = []
preds_shannon = []
preds_tsallis = []
scores_cosine = []
scores_shannon = []
scores_tsallis = []
for query_idx in range(0,len(unique_query_ids)): 
    q_idxs_tmp = np.where(df_query.iloc[:,0] == unique_query_ids[query_idx])[0]
    q_spec_tmp = np.asarray(pd.concat([df_query.iloc[q_idxs_tmp,1], df_query.iloc[q_idxs_tmp,2]], axis=1).reset_index(drop=True))
    q_spec = q_spec_tmp

    scores_cosine_tmp = []
    scores_shannon_tmp = []
    scores_tsallis_tmp = []
    for ref_idx in range(0,len(unique_reference_ids)):
        if ref_idx % 100 == 0:
            print(f'Query spectrum #{query_idx} has had its similarity with {ref_idx} reference library spectra computed')
        r_idxs_tmp = np.where(df_reference.iloc[:,0] == unique_reference_ids[ref_idx])[0]
        r_spec = np.asarray(pd.concat([df_reference.iloc[r_idxs_tmp,1], df_reference.iloc[r_idxs_tmp,2]], axis=1).reset_index(drop=True))

        # perform weight factor transformation
        q_spec[:,1] = wf_transform(spec_mzs=q_spec[:,0], spec_ints=q_spec[:,1], wf_mz=1.1, wf_int=0.9)
        r_spec[:,1] = wf_transform(spec_mzs=r_spec[:,0], spec_ints=r_spec[:,1], wf_mz=1.1, wf_int=0.9)

        # perform low-entropy transformation 
        q_spec[:,1] = LE_transform(intensity=q_spec[:,1], thresh=3, normalization_method='standard')
        r_spec[:,1] = LE_transform(intensity=r_spec[:,1], thresh=3, normalization_method='standard')

        # remove noise (i.e. remove low-intensity ion fragments with intensity below max(intensity)*noise_threshold
        q_spec = remove_noise(spec=q_spec, nr=0.01)
        r_spec = remove_noise(spec=r_spec, nr=0.01)

        # centroid spectra
        q_spec = centroid_spectrum(spec=q_spec, window_size=0.05) 
        r_spec = centroid_spectrum(spec=r_spec, window_size=0.05) 

        # match m/z values so spectra intensities can be represented by vectors of the same length
        m_spec = match_peaks_in_spectra(spec_a=q_spec, spec_b=r_spec, window_size=0.05)
        q_spec = m_spec[:,0:2]
        r_spec = m_spec[:,[0,2]]

        q_ints = q_spec[:,1]
        r_ints = r_spec[:,1]

        # normalize intensities so they sum to 1 to represent a probability distribution (only necessary for entropy-based similarity measures)
        q_ints = normalize(intensities=q_ints, method='standard')
        r_ints = normalize(intensities=r_ints, method='standard')

        # compute similarity using three different similarity measures
        sim_cosine = S_cos(ints_a=q_ints, ints_b=r_ints)
        sim_shannon = S_shannon(ints_a=q_ints, ints_b=r_ints)
        sim_tsallis = S_tsallis(ints_a=q_ints, ints_b=r_ints, q=1.1)
        scores_cosine_tmp.append(sim_cosine)
        scores_shannon_tmp.append(sim_shannon)
        scores_tsallis_tmp.append(sim_tsallis)

    # the predicted compound is the compound from the reference library with the largest similarity score
    preds_cosine.append(unique_reference_ids[np.where(scores_cosine_tmp == max(scores_cosine_tmp))[0][0]])
    preds_shannon.append(unique_reference_ids[np.where(scores_shannon_tmp == max(scores_shannon_tmp))[0][0]])
    preds_tsallis.append(unique_reference_ids[np.where(scores_tsallis_tmp == max(scores_tsallis_tmp))[0][0]])
    scores_cosine.append(max(scores_cosine_tmp))
    scores_shannon.append(max(scores_shannon_tmp))
    scores_tsallis.append(max(scores_tsallis_tmp))

df_lcms = pd.DataFrame({'QUERY_ID':unique_query_ids, 'PREDICTED_COMPOUND_COSINE':preds_cosine, 'SIMILARITY_SCORE_COSINE':scores_cosine, 'PREDICTED_COMPOUND_SHANNON':preds_shannon, 'SIMILARITY_SCORE_SHANNON':scores_shannon, 'PREDICTED_COMPOUND_TSALLIS':preds_tsallis, 'SIMILARITY_SCORE_TSALLIS':scores_tsallis})

print('\n Example LC-MS/MS results:')
print(df_lcms)






###################################### GC-MS Example ######################################
print('\n\n\nGC-MS example:\n')

# import query and reference libraries
df_query = pd.read_csv(os.path.join(os.getcwd(), '../data') + '/gcms_query_library.csv')
df_reference = pd.read_csv(os.path.join(os.getcwd(), '../data') + '/gcms_reference_library.csv')

# get unique spectrum IDs
unique_query_ids = df_query.iloc[:,0].unique()
unique_reference_ids = df_reference.iloc[:,0].unique()

# get the range of m/z values
min_mz = np.min([np.min(df_query.iloc[:,1]), np.min(df_reference.iloc[:,1])])
max_mz = np.max([np.max(df_query.iloc[:,1]), np.max(df_reference.iloc[:,1])])
mzs = np.linspace(min_mz,max_mz,(max_mz-min_mz+1))

# compute the similarity score between each query spectrum/spectra and all reference spectra and record the predicted compound for each query along with its corresponding similarity score
preds_cosine = []
preds_shannon = []
preds_tsallis = []
scores_cosine = []
scores_shannon = []
scores_tsallis = []
for query_idx in range(0,len(unique_query_ids)):
    q_idxs_tmp = np.where(df_query.iloc[:,0] == unique_query_ids[query_idx])[0]
    q_spec_tmp = np.asarray(pd.concat([df_query.iloc[q_idxs_tmp,1], df_query.iloc[q_idxs_tmp,2]], axis=1).reset_index(drop=True))
    q_spec_tmp = convert_spec(q_spec_tmp,mzs)

    scores_cosine_tmp = []
    scores_shannon_tmp = []
    scores_tsallis_tmp = []
    for ref_idx in range(0,len(unique_reference_ids)):
        q_spec = q_spec_tmp
        if ref_idx % 1000 == 0:
            print(f'Query spectrum #{query_idx} has had its similarity with {ref_idx} reference library spectra computed')
        r_idxs_tmp = np.where(df_reference.iloc[:,0] == unique_reference_ids[ref_idx])[0]
        r_spec_tmp = np.asarray(pd.concat([df_reference.iloc[r_idxs_tmp,1], df_reference.iloc[r_idxs_tmp,2]], axis=1).reset_index(drop=True))
        r_spec = convert_spec(r_spec_tmp,mzs)

        # perform weight factor transformation
        q_spec[:,1] = wf_transform(spec_mzs=q_spec[:,0], spec_ints=q_spec[:,1], wf_mz=1.1, wf_int=0.9)
        r_spec[:,1] = wf_transform(spec_mzs=r_spec[:,0], spec_ints=r_spec[:,1], wf_mz=1.1, wf_int=0.9)

        # perform low-entropy transformation 
        q_spec[:,1] = LE_transform(intensity=q_spec[:,1], thresh=3, normalization_method='standard')
        r_spec[:,1] = LE_transform(intensity=r_spec[:,1], thresh=3, normalization_method='standard')

        # remove noise (i.e. remove low-intensity ion fragments with intensity below max(intensity)*noise_threshold
        q_spec = remove_noise(spec=q_spec, nr=0.01)
        r_spec = remove_noise(spec=r_spec, nr=0.01)

        q_ints = q_spec[:,1]
        r_ints = r_spec[:,1]

        # normalize intensities so they sum to 1 to represent a probability distribution (only necessary for entropy-based similarity measures)
        q_ints = normalize(intensities=q_ints, method='standard')
        r_ints = normalize(intensities=r_ints, method='standard')

        # compute similarity using three different similarity measures
        sim_cosine = S_cos(ints_a=q_ints, ints_b=r_ints)
        sim_shannon = S_shannon(ints_a=q_ints, ints_b=r_ints)
        sim_tsallis = S_tsallis(ints_a=q_ints, ints_b=r_ints, q=1.1)
        scores_cosine_tmp.append(sim_cosine)
        scores_shannon_tmp.append(sim_shannon)
        scores_tsallis_tmp.append(sim_tsallis)

    # the predicted compound is the compound from the reference library with the largest similarity score
    preds_cosine.append(unique_reference_ids[np.where(scores_cosine_tmp == max(scores_cosine_tmp))[0][0]])
    preds_shannon.append(unique_reference_ids[np.where(scores_shannon_tmp == max(scores_shannon_tmp))[0][0]])
    preds_tsallis.append(unique_reference_ids[np.where(scores_tsallis_tmp == max(scores_tsallis_tmp))[0][0]])
    scores_cosine.append(max(scores_cosine_tmp))
    scores_shannon.append(max(scores_shannon_tmp))
    scores_tsallis.append(max(scores_tsallis_tmp))

df_gcms = pd.DataFrame({'QUERY_ID':unique_query_ids, 'PREDICTED_COMPOUND_COSINE':preds_cosine, 'SIMILARITY_SCORE_COSINE':scores_cosine, 'PREDICTED_COMPOUND_SHANNON':preds_shannon, 'SIMILARITY_SCORE_SHANNON':scores_shannon, 'PREDICTED_COMPOUND_TSALLIS':preds_tsallis, 'SIMILARITY_SCORE_TSALLIS':scores_tsallis})

print('\n Example GC-MS/MS results:')
print(df_gcms)




