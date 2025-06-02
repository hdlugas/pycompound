
# this script performs spectral library matching to identify unknown query compound(s) from GC-MS data

# load libraries
from processing import *
from similarity_measures import *
import pandas as pd
import argparse
from pathlib import Path
import sys


# create argparse object so command-line input can be imported
parser = argparse.ArgumentParser()

# import optional command-line arguments
parser.add_argument('--query_data', type=str, metavar='\b', help='CSV file of query mass spectrum/spectra to be identified. Each row should correspond to a mass spectrum, the left-most column should contain an identifier, and each of the other columns should correspond to a single mass/charge ratio. Mandatory argument.')
parser.add_argument('--reference_data', type=str, metavar='\b', help='CSV file of the reference mass spectra. Each row should correspond to a mass spectrum, the left-most column should contain in identifier (i.e. the CAS registry number or the compound name), and the remaining column should correspond to a single mass/charge ratio. Mandatory argument.')
parser.add_argument('--likely_reference_IDs', type=str, metavar='\b', help='CSV file with one column containing the IDs of a subset of all compounds in the reference_data to be used in spectral library matching. Each ID in this file must be an ID in the reference library. Default: none (i.e. default is to use entire reference library)')
parser.add_argument('--similarity_measure', type=str, default='cosine', metavar='\b', help='Similarity measure: options are \'cosine\', \'shannon\', \'renyi\', and \'tsallis\'. Default: cosine.')
parser.add_argument('--chromatography_platform', type=str, metavar='\b', help='Chromatography platform: options are \'HRMS\' and \'NRMS\'. Mandatory argument.')
parser.add_argument('--spectrum_preprocessing_order', type=str, metavar='\b', help='The LC-MS/MS spectrum preprocessing transformations and the order in which they are to be applied. Note that these transformations are applied prior to computing similarity scores. Format must be a string with 2-6 characters chosen from C, F, M, N, L, W representing centroiding, filtering based on mass/charge and intensity values, matching, noise removal, low-entropy trannsformation, and weight-factor-transformation, respectively. For example, if \'WCM\' is passed, then each spectrum will undergo a weight factor transformation, then centroiding, and then matching. Note that if an argument is passed, then \'M\' must be contained in the argument, since matching is a required preprocessing step in spectral library matching of LC-MS/MS data. Furthermore, \'C\' must be performed before matching since centroiding can change the number of ion fragments in a given spectrum. Default: FCNMWL for HRMS, FNLW for NRMS')
parser.add_argument('--high_quality_reference_library', type=str, default='False', metavar='\b', help='True/False flag indicating whether the reference library is considered to be of high quality. If True, then the spectrum preprocessing transformations of filtering and noise removal are performed only on the query spectrum/spectra. If False, all spectrum preprocessing transformations specified will be applied to both the query and reference spectra. Default: False')
parser.add_argument('--mz_min', type=int, default=0, metavar='\b', help='Remove all peaks with mass/charge less than mz_min in each spectrum. Default: 0')
parser.add_argument('--mz_max', type=int, default=999999999999, metavar='\b', help='Remove all peaks with mass/charge greater than mz_max in each spectrum. Default: 999999999999')
parser.add_argument('--int_min', type=float, default=0, metavar='\b', help='Remove all peaks with intensity less than int_min in each spectrum. Default: 0')
parser.add_argument('--int_max', type=float, default=999999999999, metavar='\b', help='Remove all peaks with intensity greater than int_max in each spectrum. Default: 999999999999')
parser.add_argument('--window_size_centroiding', type=float, default=0.5, metavar='\b', help='Window size parameter used in centroiding a given spectrum. Only for HRMS. Default: 0.5')
parser.add_argument('--window_size_matching', type=float, default=0.5, metavar='\b', help='Window size parameter used in matching a query spectrum and a reference library spectrum. Only for HRMS. Default: 0.5')
parser.add_argument('--noise_threshold', type=float, default=0, metavar='\b', help='Ion fragments (i.e. points in a given mass spectrum) with intensity less than max(intensities)*noise_threshold are removed. Default: 0')
parser.add_argument('--wf_mz', type=float, default=0, metavar='\b', help='Mass/charge weight factor parameter. Default: 0.')
parser.add_argument('--wf_intensity', type=float, default=1, metavar='\b', help='Intensity weight factor parameter. Default: 1.')
parser.add_argument('--LET_threshold', type=float, default=0, metavar='\b', help='Low-entropy transformation threshold parameter. Spectra with Shannon entropy less than LET_threshold are transformed according to intensitiesNew=intensitiesOriginal^{(1+S)/(1+LET_threshold)}. Default: 0.')
parser.add_argument('--entropy_dimension', type=float, default=1.1, metavar='\b', help='Entropy dimension parameter. Must have positive value other than 1. When the entropy dimension is 1, then Renyi and Tsallis entropy are equivalent to Shannon entropy. Therefore, this parameter only applies to the renyi and tsallis similarity measures. This parameter will be ignored if similarity measure cosine or shannon is chosen. Default: 1.1.')
parser.add_argument('--normalization_method', type=str, default='standard', metavar='\b', help='Method used to normalize the intensities of each spectrum so that the intensities sum to 1. Since the objects entropy quantifies the uncertainy of must be probability distributions, the intensities of a given spectrum must sum to 1 prior to computing the entropy of the given spectrum intensities. Options: \'standard\' and \'softmax\'. Default: standard.')
parser.add_argument('--n_top_matches_to_save', type=int, default=1, metavar='\b', help='The number of top matches to report. For example, if n_top_matches_to_save=5, then for each query spectrum, the five reference spectra with the largest similarity with the given query spectrum will be reported. Default: 1.')
parser.add_argument('--print_id_results', type=str, default='False', metavar='\b', help='Flag that prints identification results if True. Default: False')
parser.add_argument('--output_identification', type=str, default=f'{Path.cwd()}/output_identification.csv', metavar='\b', help='Output CSV file containing the most-similar reference spectra for each query spectrum along with the corresponding similarity scores. Default is to save identification output in current working directory (i.e. same directory this script is contained in) with filename \'output_identification.csv\'.')
parser.add_argument('--output_similarity_scores', type=str, default=f'{Path.cwd()}/output_all_similarity_scores.csv', metavar='\b', help='Output CSV file containing similarity scores between all query spectrum/spectra and all reference spectra. Each row corresponds to a query spectrum, the left-most column contains the query spectrum/spectra identifier, and the remaining column contain the similarity scores with respect to all reference library spectra. If no argument passed, then this CSV file is written to the current working directory with filename \'output_all_similarity_scores\'.csv.')


# parse the user-input arguments
args = parser.parse_args()


print('\nPerforming spectral library matching')


# import the query library
if args.query_data is not None:
    df_query = pd.read_csv(args.query_data)
else:
    print('\nError: No argument passed to query_data. To view usage, run \"python plot_spectra.py -h\".')
    sys.exit()


# load the reference library
if args.reference_data is not None:
    df_reference = pd.read_csv(args.reference_data)
else:
    print('\nError: No argument passed to reference_data. To view usage, run \"python plot_spectra.py -h\".')
    sys.exit()


# load the IDs of a subset of the reference library for use in spectral library matching
if args.likely_reference_IDs is not None:
    if args.likely_reference_IDs != 'None' and args.likely_reference_IDs != 'none':
        likely_reference_IDs = pd.read_csv(args.likely_reference_IDs, header=None)
        df_reference = df_reference.loc[df_reference.iloc[:,0].isin(likely_reference_IDs.iloc[:,0].tolist())]


# throw error if similarity measure is not one of the available similarity measures
similarity_measure = args.similarity_measure
if similarity_measure not in ['cosine','shannon','renyi','tsallis']:
    print('\nError: similarity_measure must be either \'cosine\', \'shannon\', \'tsallis\'')
    sys.exit()


# specify the chromatography platform
if args.chromatography_platform is not None:
    chromatography_platform = args.chromatography_platform
else:
    print('\nError: No argument passed to chromatography_platform. To view usage, run \"python plot_spectra.py -h\".')
    sys.exit()

if chromatography_platform not in ['HRMS','NRMS']:
    print('\nError: chromatography_platform must be either \'HRMS\' or \'NRMS\'')
    sys.exit()


# get the spectrum preprocessing order
if chromatography_platform == 'HRMS':
    preprocessing_error_message1 = 'Error: \'M\' must be a character in spectrum_preprocessing_order.'
    preprocessing_error_message2 = 'Error: \'C\' must come before \'M\' in spectrum_preprocessing_order.'
    preprocessing_error_message3 = 'Error: spectrum_preprocessing_order must contain only \'C\', \'F\', \'M\', \'N\', \'L\', \'W\' for HRMS'

    if args.spectrum_preprocessing_order is not None:
        spectrum_preprocessing_order = list(args.spectrum_preprocessing_order)
    else:
        spectrum_preprocessing_order = ['F', 'C', 'N', 'M', 'W', 'L']

    if 'M' not in spectrum_preprocessing_order:
        print(f'\n{preprocessing_error_message1}\n')
        sys.exit()

    if 'C' in spectrum_preprocessing_order:
        if spectrum_preprocessing_order.index('C') > spectrum_preprocessing_order.index('M'):
            print(f'\n{preprocessing_error_message2}\n')
            sys.exit()

    if set(spectrum_preprocessing_order) - {'F','C','N','M','W','L'}:
        print(f'\n{preprocessing_error_message3}')
        sys.exit()

elif chromatography_platform == 'NRMS':
    preprocessing_error_message3 = 'Error: spectrum_preprocessing_order must contain only \'F\', \'N\', \'L\', \'W\' for NRMS'
    if args.spectrum_preprocessing_order is not None:
        spectrum_preprocessing_order = list(args.spectrum_preprocessing_order)
    else:
        spectrum_preprocessing_order = ['F', 'N', 'L', 'W']

    if set(spectrum_preprocessing_order) - {'F','N','W','L'}:
        print(f'\n{preprocessing_error_message3}')
        sys.exit()


# load the flag indicating whether the reference library is considered to be of high quality
high_quality_reference_library = args.high_quality_reference_library


# load the filtering parameters
mz_min = int(args.mz_min)
if mz_min < 0:
    print('\nError: mz_min should be a non-negative integer')
    sys.exit()

mz_max = int(args.mz_max)
if mz_max <= 0:
    print('\nError: mz_max should be a positive integer')
    sys.exit()

int_min = float(args.int_min)
if int_min < 0:
    print('\nError: int_min should be a non-negative float')
    sys.exit()

int_max = float(args.int_max)
if int_max <= 0:
    print('\nError: int_max should be a positive float')
    sys.exit()


# load the centroiding window size parameter
if args.window_size_centroiding is not None:
    try:
        window_size_centroiding = float(args.window_size_centroiding)
    except ValueError:
        print('\nError: window_size_centroiding should be a positive float')
        sys.exit()
else:
    window_size_centroiding = 0.5


# load the matching window size parameter
if args.window_size_matching is not None:
    try:
        window_size_matching = float(args.window_size_matching)
    except ValueError:
        print('\nError: window_size_matching should be a positive float')
        sys.exit()
else:
    window_size_matching = 0.5


# load the noise removal parameter
if args.noise_threshold is not None:
    try:
        noise_threshold = float(args.noise_threshold)
    except ValueError:
        print('\nError: noise_threshold should be a positive float')
        sys.exit()
else:
    noise_threshold = 0


# load the weight factor parameters
if args.wf_mz is not None:
    try:
        wf_mz = float(args.wf_mz)
    except ValueError:
        print('\nError: wf_mz should be a float')
        sys.exit()
else:
    wf_mz = 0

if args.wf_intensity is not None:
    try:
        wf_intensity = float(args.wf_intensity)
    except ValueError:
        print('\nError: wf_intensity should be a float')
        sys.exit()
else:
    wf_intensity = 1


# load the low-entropy transformation threshold
LET_threshold = float(args.LET_threshold)


# load the entropy dimension parameter (if applicable)
if args.similarity_measure == 'renyi' or args.similarity_measure == 'tsallis':
    q = float(args.entropy_dimension)
    if q <= 0:
        print('\nError: entropy_dimension should be a positive float')
        sys.exit()


# load the normalization method
normalization_method = args.normalization_method
if normalization_method not in ['softmax','standard']:
    print('\nError: normalization_method must be either \'softmax\' or \'standard\'')
    sys.exit()


# get the number of most-similar reference library spectra to report for each query spectrum
n_top_matches_to_save = int(args.n_top_matches_to_save)
if n_top_matches_to_save <= 0:
    print('\nError: n_top_matches_to_save should be a positive integer')
    sys.exit()


# get the flag to determine whether or not to print identification results
print_id_results = str(args.print_id_results)
if print_id_results not in ['True','False']:
    print('\nError: print_id_results must be either \'True\' or \'False\'')
    sys.exit()


# consider the cases of HRMS and NRMS separately
if chromatography_platform == 'HRMS':

    # get unique query/reference library IDs; each query/reference ID corresponds to exactly one query/reference mass spectrum
    unique_query_ids = df_query.iloc[:,0].unique()
    unique_reference_ids = df_reference.iloc[:,0].unique()

    # compute the similarity score between each query library spectrum/spectra and all reference library spectra
    all_similarity_scores =  []
    for query_idx in range(0,len(unique_query_ids)):
        q_idxs_tmp = np.where(df_query.iloc[:,0] == unique_query_ids[query_idx])[0]
        q_spec_tmp = np.asarray(pd.concat([df_query.iloc[q_idxs_tmp,1], df_query.iloc[q_idxs_tmp,2]], axis=1).reset_index(drop=True))

        # compute the similarity score between the given query spectrum and all spectra in the reference library
        similarity_scores = []
        for ref_idx in range(0,len(unique_reference_ids)):
            #if ref_idx % 100 == 0:
            #    print(f'Query spectrum #{query_idx} has had its similarity with {ref_idx} reference library spectra computed')
            q_spec = q_spec_tmp
            r_idxs_tmp = np.where(df_reference.iloc[:,0] == unique_reference_ids[ref_idx])[0]
            r_spec = np.asarray(pd.concat([df_reference.iloc[r_idxs_tmp,1], df_reference.iloc[r_idxs_tmp,2]], axis=1).reset_index(drop=True))

            # apply spectrum preprocessing transformation in the order specified by user
            is_matched = False
            for transformation in spectrum_preprocessing_order:
                if np.isinf(q_spec[:,1]).sum() > 0:
                    q_spec[:,1] = np.zeros(q_spec.shape[0])
                if np.isinf(r_spec[:,1]).sum() > 0:
                    r_spec[:,1] = np.zeros(r_spec.shape[0])
                if transformation == 'C' and q_spec.shape[0] > 1 and r_spec.shape[1] > 1: # centroiding
                    q_spec = centroid_spectrum(q_spec, window_size=window_size_centroiding) 
                    r_spec = centroid_spectrum(r_spec, window_size=window_size_centroiding) 
                if transformation == 'M' and q_spec.shape[0] > 1 and r_spec.shape[1] > 1: # matching
                    m_spec = match_peaks_in_spectra(spec_a=q_spec, spec_b=r_spec, window_size=window_size_matching)
                    q_spec = m_spec[:,0:2]
                    r_spec = m_spec[:,[0,2]]
                    is_matched = True
                if transformation == 'W' and q_spec.shape[0] > 1 and r_spec.shape[1] > 1: # weight factor transformation
                    q_spec[:,1] = wf_transform(q_spec[:,0], q_spec[:,1], wf_mz, wf_intensity)
                    r_spec[:,1] = wf_transform(r_spec[:,0], r_spec[:,1], wf_mz, wf_intensity)
                if transformation == 'L' and q_spec.shape[0] > 1 and r_spec.shape[1] > 1: # low-entropy tranformation
                    q_spec[:,1] = LE_transform(q_spec[:,1], LET_threshold, normalization_method=normalization_method)
                    r_spec[:,1] = LE_transform(r_spec[:,1], LET_threshold, normalization_method=normalization_method)
                if transformation == 'N' and q_spec.shape[0] > 1 and r_spec.shape[1] > 1: # noise removal
                    q_spec = remove_noise(q_spec, nr = noise_threshold)
                    if high_quality_reference_library == False:
                        r_spec = remove_noise(r_spec, nr = noise_threshold)
                if transformation == 'F' and q_spec.shape[0] > 1 and r_spec.shape[1] > 1: # filter with respect to mz and/or intensity
                    q_spec = filter_spec_lcms(q_spec, mz_min = mz_min, mz_max = mz_max, int_min = int_min, int_max = int_max, is_matched = is_matched)
                    if high_quality_reference_library == False:
                        r_spec = filter_spec_lcms(r_spec, mz_min = mz_min, mz_max = mz_max, int_min = int_min, int_max = int_max, is_matched = is_matched)

            # query and reference spectrum intensities
            q_ints = q_spec[:,1]
            r_ints = r_spec[:,1]

            if np.sum(q_ints) != 0 and np.sum(r_ints) != 0 and q_spec.shape[0] > 1 and r_spec.shape[1] > 1:
                if similarity_measure == 'cosine':
                    similarity_score = S_cos(q_ints, r_ints)
                else:
                    q_ints = normalize(q_ints, method = normalization_method)
                    r_ints = normalize(r_ints, method = normalization_method)

                    if similarity_measure == 'shannon':
                        similarity_score = S_shannon(q_ints, r_ints)
                    elif similarity_measure == 'renyi':
                        similarity_score = S_renyi(q_ints, r_ints, q)
                    elif similarity_measure == 'tsallis':
                        similarity_score = S_tsallis(q_ints, r_ints, q)
            else:
                similarity_score = 0

            similarity_scores.append(similarity_score)
        all_similarity_scores.append(similarity_scores)


elif chromatography_platform == 'NRMS':

    # get unique query/reference library IDs; each query/reference ID corresponds to exactly one query/reference mass spectrum
    unique_query_ids = df_query.iloc[:,0].unique()
    unique_reference_ids = df_reference.iloc[:,0].unique()

    # get the range of m/z values
    min_mz = np.min([np.min(df_query.iloc[:,1]), np.min(df_reference.iloc[:,1])])
    max_mz = np.max([np.max(df_query.iloc[:,1]), np.max(df_reference.iloc[:,1])])
    mzs = np.linspace(min_mz,max_mz,(max_mz-min_mz+1))

    # compute the similarity score between each query library spectrum/spectra and all reference library spectra
    # for each query spectrum, compute its similarity with all reference spectra 
    all_similarity_scores =  []
    for query_idx in range(0,len(unique_query_ids)):
        q_idxs_tmp = np.where(df_query.iloc[:,0] == unique_query_ids[query_idx])[0]
        q_spec_tmp = np.asarray(pd.concat([df_query.iloc[q_idxs_tmp,1], df_query.iloc[q_idxs_tmp,2]], axis=1).reset_index(drop=True))
        q_spec_tmp = convert_spec(q_spec_tmp,mzs)

        similarity_scores = []
        for ref_idx in range(0,len(unique_reference_ids)):
            q_spec = q_spec_tmp
            if ref_idx % 1000 == 0:
                print(f'Query spectrum #{query_idx} has had its similarity with {ref_idx} reference library spectra computed')
            r_idxs_tmp = np.where(df_reference.iloc[:,0] == unique_reference_ids[ref_idx])[0]
            r_spec_tmp = np.asarray(pd.concat([df_reference.iloc[r_idxs_tmp,1], df_reference.iloc[r_idxs_tmp,2]], axis=1).reset_index(drop=True))
            r_spec = convert_spec(r_spec_tmp,mzs)


            # apply spectrum preprocessing transformation in the order specified by user
            for transformation in spectrum_preprocessing_order:
                if np.isinf(q_spec[:,1]).sum() > 0:
                    q_spec[:,1] = np.zeros(q_spec.shape[0])
                if np.isinf(r_spec[:,1]).sum() > 0:
                    r_spec[:,1] = np.zeros(r_spec.shape[0])
                if transformation == 'W': # weight factor transformation
                    q_spec[:,1] = wf_transform(q_spec[:,0], q_spec[:,1], wf_mz, wf_intensity)
                    r_spec[:,1] = wf_transform(r_spec[:,0], r_spec[:,1], wf_mz, wf_intensity)
                if transformation == 'L': # low-entropy transformation
                    q_spec[:,1] = LE_transform(q_spec[:,1], LET_threshold, normalization_method=normalization_method)
                    r_spec[:,1] = LE_transform(r_spec[:,1], LET_threshold, normalization_method=normalization_method)
                if transformation == 'N': # noise removal
                    q_spec = remove_noise(q_spec, nr = noise_threshold)
                    if high_quality_reference_library == False:
                        r_spec = remove_noise(r_spec, nr = noise_threshold)
                if transformation == 'F': # filter with respect to mz and/or intensity
                    q_spec = filter_spec_gcms(q_spec, mz_min = mz_min, mz_max = mz_max, int_min = int_min, int_max = int_max)
                    if high_quality_reference_library == False:
                        r_spec = filter_spec_gcms(r_spec, mz_min = mz_min, mz_max = mz_max, int_min = int_min, int_max = int_max)

            # query and reference spectrum intensities
            q_ints = q_spec[:,1]
            r_ints = r_spec[:,1]

            # if there are no non-zero intensities in the query or reference spectrum, their similarity is 0
            if np.sum(q_ints) != 0 and np.sum(r_ints) != 0:
                if similarity_measure == 'cosine':
                    similarity_score = S_cos(q_ints, r_ints)
                else:
                    # normalize intensities of each spectrum so they sum to 1 so that they represent a probability distribution
                    q_ints = normalize(q_ints, method = normalization_method)
                    r_ints = normalize(r_ints, method = normalization_method)

                    if similarity_measure == 'shannon':
                        similarity_score = S_shannon(q_ints, r_ints)
                    elif similarity_measure == 'renyi':
                        similarity_score = S_renyi(q_ints, r_ints, q)
                    elif similarity_measure == 'tsallis':
                        similarity_score = S_tsallis(q_ints, r_ints, q)
            else:
                similarity_score = 0

            similarity_scores.append(similarity_score)
        all_similarity_scores.append(similarity_scores)




# create pandas dataframe containing all similarity scores computed with one row for each query spectrum and one column for each reference spectrum
df_scores = pd.DataFrame(all_similarity_scores, columns = unique_reference_ids)
df_scores.index = unique_query_ids
df_scores.index.names = ['Query Spectrum ID']


# get predicted identity/identities of each query spectrum and the corresponding maximum similarity score
preds = []
scores = []
for i in range(0, df_scores.shape[0]):
    df_scores_tmp = df_scores
    preds_tmp = []
    scores_tmp = []
    for j in range(0, n_top_matches_to_save):
        top_ref_specs_tmp = df_scores_tmp.iloc[i,np.where(df_scores_tmp.iloc[i,:] == np.max(df_scores_tmp.iloc[i,:]))[0]]
        cols_to_keep = np.where(df_scores_tmp.iloc[i,:] != np.max(df_scores_tmp.iloc[i,:]))[0]
        df_scores_tmp = df_scores_tmp.iloc[:,cols_to_keep]

        #preds_tmp.append(';'.join(top_ref_specs_tmp.index.to_list()))
        preds_tmp.append(';'.join(map(str,top_ref_specs_tmp.index.to_list())))
        if len(top_ref_specs_tmp.values) == 0:
            scores_tmp.append(0)
        else:
            scores_tmp.append(top_ref_specs_tmp.values[0])
    preds.append(preds_tmp)
    scores.append(scores_tmp)

preds = np.array(preds)
scores = np.array(scores)
out = np.c_[preds,scores]

# get column names for a pandas dataframe with the n_top_matches_to_save top-matches for each query spectrum
cnames_preds = []
cnames_scores = []
for i in range(0,n_top_matches_to_save):
    cnames_preds.append(f'RANK.{i+1}.PRED')
    cnames_scores.append(f'RANK.{i+1}.SIMILARITY.SCORE')

# get pandas dataframe with identifcation results with each row corresponding to a query spectrum, n_top_matches_to_save columns for the top predictions, and n_top_matches_to_save columns for the similarity scores corresponding to the predictions
df_top_ref_specs = pd.DataFrame(out, columns = [*cnames_preds, *cnames_scores])
df_top_ref_specs.index = unique_query_ids
df_top_ref_specs.index.names = ['Query Spectrum ID']

# print the identification results if the user desires
if print_id_results == 'True':
    print(df_top_ref_specs.to_string())


# write spectral library matching results to disk
df_top_ref_specs.to_csv(args.output_identification)


# write all similarity scores to disk
df_scores.columns = ['Reference Spectrum ID: ' + col for col in  list(map(str,df_scores.columns.tolist()))]
df_scores.to_csv(args.output_similarity_scores)

