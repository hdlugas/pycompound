
# this script's function runs spectral library matching to identify unknown query compound(s)

from pycompound_fy7392.build_library import build_library_from_raw_data
from .processing import *
from .similarity_measures import *
import pandas as pd
from pathlib import Path
import sys


def run_spec_lib_matching_on_HRMS_data(query_data=None, reference_data=None, likely_reference_IDs=None, similarity_measure='cosine', spectrum_preprocessing_order='FCNMWL', high_quality_reference_library=False, mz_min=0, mz_max=9999999, int_min=0, int_max=9999999, window_size_centroiding=0.5, window_size_matching=0.5, noise_threshold=0.0, wf_mz=0.0, wf_intensity=1.0, LET_threshold=0.0, entropy_dimension=1.1, normalization_method='standard', n_top_matches_to_save=1, print_id_results=False, output_identification=None, output_similarity_scores=None):
    '''
    runs spectral library matching on high-resolution mass spectrometry data

    --query_data: CSV file of query mass spectrum/spectra to be identified. Each row should correspond to a mass spectrum, the left-most column should contain an identifier, and each of the other columns should correspond to a single mass/charge ratio. Mandatory argument.
    --reference_data: CSV file of the reference mass spectra. Each row should correspond to a mass spectrum, the left-most column should contain in identifier (i.e. the CAS registry number or the compound name), and the remaining column should correspond to a single mass/charge ratio. Mandatory argument.
    --likely_reference_IDs: CSV file with one column containing the IDs of a subset of all compounds in the reference_data to be used in spectral library matching. Each ID in this file must be an ID in the reference library. Default: None (i.e. default is to use entire reference library)
    --similarity_measure: \'cosine\', \'shannon\', \'renyi\', and \'tsallis\'. Default: cosine.
    --spectrum_preprocessing_order: The spectrum preprocessing transformations and the order in which they are to be applied. Note that these transformations are applied prior to computing similarity scores. Format must be a string with 2-6 characters chosen from C, F, M, N, L, W representing centroiding, filtering based on mass/charge and intensity values, matching, noise removal, low-entropy trannsformation, and weight-factor-transformation, respectively. For example, if \'WCM\' is passed, then each spectrum will undergo a weight factor transformation, then centroiding, and then matching. Note that if an argument is passed, then \'M\' must be contained in the argument, since matching is a required preprocessing step in spectral library matching of HRMS data. Furthermore, \'C\' must be performed before matching since centroiding can change the number of ion fragments in a given spectrum. Default: FCNMWL')
    --high_quality_reference_library: True/False flag indicating whether the reference library is considered to be of high quality. If True, then the spectrum preprocessing transformations of filtering and noise removal are performed only on the query spectrum/spectra. If False, all spectrum preprocessing transformations specified will be applied to both the query and reference spectra. Default: False')
    --mz_min: Remove all peaks with mass/charge value less than mz_min in each spectrum. Default: 0
    --mz_max: Remove all peaks with mass/charge value greater than mz_max in each spectrum. Default: 9999999
    --int_min: Remove all peaks with intensity value less than int_min in each spectrum. Default: 0
    --int_max: Remove all peaks with intensity value greater than int_max in each spectrum. Default: 9999999
    --window_size_centroiding: Window size parameter used in centroiding a given spectrum. Default: 0.5
    --window_size_matching: Window size parameter used in matching a query spectrum and a reference library spectrum. Default: 0.5
    --noise_threshold: Ion fragments (i.e. points in a given mass spectrum) with intensity less than max(intensities)*noise_threshold are removed. Default: 0.0
    --wf_mz: Mass/charge weight factor parameter. Default: 0.0
    --wf_intensity: Intensity weight factor parameter. Default: 0.0
    --LET_threshold: Low-entropy transformation threshold parameter. Spectra with Shannon entropy less than LET_threshold are transformed according to intensitiesNew=intensitiesOriginal^{(1+S)/(1+LET_threshold)}. Default: 0.0
    --entropy_dimension: Entropy dimension parameter. Must have positive value other than 1. When the entropy dimension is 1, then Renyi and Tsallis entropy are equivalent to Shannon entropy. Therefore, this parameter only applies to the renyi and tsallis similarity measures. This parameter will be ignored if similarity measure cosine or shannon is chosen. Default: 1.1
    --normalization_method: Method used to normalize the intensities of each spectrum so that the intensities sum to 1. Since the objects entropy quantifies the uncertainy of must be probability distributions, the intensities of a given spectrum must sum to 1 prior to computing the entropy of the given spectrum intensities. Options: \'standard\' and \'softmax\'. Default: standard.
    --n_top_matches_to_save: The number of top matches to report. For example, if n_top_matches_to_save=5, then for each query spectrum, the five reference spectra with the largest similarity with the given query spectrum will be reported. Default: 1
    --print_id_results: Flag that prints identification results if True. Default: False
    --output_identification: Output CSV file containing the most-similar reference spectra for each query spectrum along with the corresponding similarity scores. Default is to save identification output in current working directory with filename \'output_identification.csv\'.
    --output_similarity_scores: Output CSV file containing similarity scores between all query spectrum/spectra and all reference spectra. Each row corresponds to a query spectrum, the left-most column contains the query spectrum/spectra identifier, and the remaining column contain the similarity scores with respect to all reference library spectra. If no argument passed, then this CSV file is written to the current working directory with filename \'output_all_similarity_scores\'.csv.')
    '''

    # load query and reference libraries
    if query_data is None:
        print('\nError: No argument passed to the mandatory query_data. Please pass the path to the CSV file of the query data.')
        sys.exit()
    else:
        extension = query_data.rsplit('.',1)
        extension = extension[(len(extension)-1)]
        if extension == 'mgf' or extension == 'MGF' or extension == 'mzML' or extension == 'mzml' or extension == 'MZML' or extension == 'cdf' or extension == 'CDF':
            output_path_tmp = query_data[:-3] + 'csv'
            build_library_from_raw_data(input_path=query_data, output_path=output_path_tmp, is_reference=False)
            df_query = pd.read_csv(output_path_tmp)
        if extension == 'csv' or extension == 'CSV':
            df_query = pd.read_csv(query_data)
        unique_query_ids = df_query.iloc[:,0].unique()

    if reference_data is None:
        print('\nError: No argument passed to the mandatory reference_data. Please pass the path to the CSV file of the reference data.')
        sys.exit()
    else:
        extension = reference_data.rsplit('.',1)
        extension = extension[(len(extension)-1)]
        if extension == 'mgf' or extension == 'MGF' or extension == 'mzML' or extension == 'mzml' or extension == 'MZML' or extension == 'cdf' or extension == 'CDF':
            output_path_tmp = reference_data[:-3] + 'csv'
            build_library_from_raw_data(input_path=reference_data, output_path=output_path_tmp, is_reference=True)
            df_reference = pd.read_csv(output_path_tmp)
        if extension == 'csv' or extension == 'CSV':
            df_reference = pd.read_csv(reference_data)

        if likely_reference_IDs is not None:
            likely_reference_IDs = pd.read_csv(likely_reference_IDs, header=None)
            df_reference = df_reference.loc[df_reference.iloc[:,0].isin(likely_reference_IDs.iloc[:,0].tolist())]
        unique_reference_ids = df_reference.iloc[:,0].unique()


    ##### process input parameters and ensure they are in a valid format #####
    if spectrum_preprocessing_order is not None:
        spectrum_preprocessing_order = list(spectrum_preprocessing_order)
    else:
        spectrum_preprocessing_order = ['F', 'C', 'N', 'M', 'W', 'L']
    if 'M' not in spectrum_preprocessing_order:
        print(f'Error: \'M\' must be a character in spectrum_preprocessing_order.')
        sys.exit()
    if 'C' in spectrum_preprocessing_order:
        if spectrum_preprocessing_order.index('C') > spectrum_preprocessing_order.index('M'):
            print(f'Error: \'C\' must come before \'M\' in spectrum_preprocessing_order.')
            sys.exit()
    if set(spectrum_preprocessing_order) - {'F','C','N','M','W','L'}:
        print(f'Error: spectrum_preprocessing_order must contain only \'C\', \'F\', \'M\', \'N\', \'L\', \'W\'.')
        sys.exit()

    if similarity_measure not in ['cosine','shannon','renyi','tsallis']:
        print('\nError: similarity_measure must be either \'cosine\', \'shannon\', \'renyi\', or \'tsallis\'')
        sys.exit()

    if isinstance(int_min,int) is True:
        int_min = float(int_min)
    if isinstance(int_max,int) is True:
        int_max = float(int_max)
    if isinstance(mz_min,int) is False or isinstance(mz_max,int) is False or isinstance(int_min,float) is False or isinstance(int_max,float) is False:
        print('Error: mz_min must be a non-negative integer, mz_max must be a positive integer, int_min must be a non-negative float, and int_max must be a positive float')
        sys.exit()
    if mz_min < 0:
        print('\nError: mz_min should be a non-negative integer')
        sys.exit()
    if mz_max <= 0:
        print('\nError: mz_max should be a positive integer')
        sys.exit()
    if int_min < 0:
        print('\nError: int_min should be a non-negative float')
        sys.exit()
    if int_max <= 0:
        print('\nError: int_max should be a positive float')
        sys.exit()

    if isinstance(window_size_centroiding,float) is False or window_size_centroiding <= 0.0:
        print('Error: window_size_centroiding must be a positive float.')
        sys.exit()
    if isinstance(window_size_matching,float) is False or window_size_matching<= 0.0:
        print('Error: window_size_matching must be a positive float.')
        sys.exit()

    if isinstance(noise_threshold,int) is True:
        noise_threshold = float(noise_threshold)
    if isinstance(noise_threshold,float) is False or noise_threshold < 0:
        print('Error: noise_threshold must be a positive float.')
        sys.exit()

    if isinstance(wf_intensity,int) is True:
        wf_intensity = float(wf_intensity)
    if isinstance(wf_mz,int) is True:
        wf_mz = float(wf_mz)
    if isinstance(wf_intensity,float) is False or isinstance(wf_mz,float) is False:
        print('Error: wf_mz and wf_intensity must be integers or floats')
        sys.exit()

    if entropy_dimension <= 0:
        print('\nError: entropy_dimension should be a positive float')
        sys.exit()
    else:
        q = entropy_dimension

    if normalization_method not in ['softmax','standard']:
        print('\nError: normalization_method must be either \'softmax\' or \'standard\'')
        sys.exit()

    if n_top_matches_to_save <= 0 or isinstance(n_top_matches_to_save,int)==False:
        print('\nError: n_top_matches_to_save should be a positive integer')
        sys.exit()

    if isinstance(print_id_results,bool)==False:
        print('\nError: print_id_results must be either True or False')
        sys.exit()
    
    if output_identification is None:
        output_identification = f'{Path.cwd()}/output_identification.csv'
        print(f'Warning: writing identification output to {output_identification}')

    if output_similarity_scores is None:
        output_similarity_scores = f'{Path.cwd()}/output_all_similarity_scores.csv'
        print(f'Warning: writing similarity scores to {output_similarity_scores}')


    ####################################### begin spectral library matching #######################################
    # compute the similarity score between each query library spectrum/spectra and all reference library spectra
    all_similarity_scores =  []
    print(len(unique_query_ids))
    for query_idx in range(0,len(unique_query_ids)):
        print(f'query spectrum #{query_idx} is being identified')
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
    if print_id_results == True:
        print(df_top_ref_specs.to_string())

    # write spectral library matching results to disk
    df_top_ref_specs.to_csv(output_identification)

    # write all similarity scores to disk
    df_scores.columns = ['Reference Spectrum ID: ' + col for col in  list(map(str,df_scores.columns.tolist()))]
    df_scores.to_csv(output_similarity_scores)





def run_spec_lib_matching_on_NRMS_data(query_data=None, reference_data=None, likely_reference_IDs=None, spectrum_preprocessing_order='FNLW', similarity_measure='cosine', high_quality_reference_library=False, mz_min=0, mz_max=9999999, int_min=0, int_max=9999999, noise_threshold=0.0, wf_mz=0.0, wf_intensity=1.0, LET_threshold=0.0, entropy_dimension=1.1, normalization_method='standard', n_top_matches_to_save=1, print_id_results=False, output_identification=None, output_similarity_scores=None):
    '''
    runs spectral library matching on nominal-resolution mass spectrometry data

    --query_data: CSV file of query mass spectrum/spectra to be identified. Each row should correspond to a mass spectrum, the left-most column should contain an identifier, and each of the other columns should correspond to a single mass/charge ratio. Mandatory argument.
    --reference_data: CSV file of the reference mass spectra. Each row should correspond to a mass spectrum, the left-most column should contain in identifier (i.e. the CAS registry number or the compound name), and the remaining column should correspond to a single mass/charge ratio. Mandatory argument.
    --likely_reference_IDs: CSV file with one column containing the IDs of a subset of all compounds in the reference_data to be used in spectral library matching. Each ID in this file must be an ID in the reference library. Default: None (i.e. default is to use entire reference library)
    --similarity_measure: \'cosine\', \'shannon\', \'renyi\', and \'tsallis\'. Default: cosine.
    --spectrum_preprocessing_order: The spectrum preprocessing transformations and the order in which they are to be applied. Note that these transformations are applied prior to computing similarity scores. Format must be a string with 2-4 characters chosen from F, N, L, W representing filtering based on mass/charge and intensity values, noise removal, low-entropy trannsformation, and weight-factor-transformation, respectively. For example, if \'WN\' is passed, then each spectrum will undergo a weight factor transformation and then noise removal. Default: FNLW')
    --high_quality_reference_library: True/False flag indicating whether the reference library is considered to be of high quality. If True, then the spectrum preprocessing transformations of filtering and noise removal are performed only on the query spectrum/spectra. If False, all spectrum preprocessing transformations specified will be applied to both the query and reference spectra. Default: False')
    --mz_min: Remove all peaks with mass/charge value less than mz_min in each spectrum. Default: 0
    --mz_max: Remove all peaks with mass/charge value greater than mz_max in each spectrum. Default: 9999999
    --int_min: Remove all peaks with intensity value less than int_min in each spectrum. Default: 0
    --int_max: Remove all peaks with intensity value greater than int_max in each spectrum. Default: 9999999
    --noise_threshold: Ion fragments (i.e. points in a given mass spectrum) with intensity less than max(intensities)*noise_threshold are removed. Default: 0.0
    --wf_mz: Mass/charge weight factor parameter. Default: 0.0
    --wf_intensity: Intensity weight factor parameter. Default: 0.0
    --LET_threshold: Low-entropy transformation threshold parameter. Spectra with Shannon entropy less than LET_threshold are transformed according to intensitiesNew=intensitiesOriginal^{(1+S)/(1+LET_threshold)}. Default: 0.0
    --entropy_dimension: Entropy dimension parameter. Must have positive value other than 1. When the entropy dimension is 1, then Renyi and Tsallis entropy are equivalent to Shannon entropy. Therefore, this parameter only applies to the renyi and tsallis similarity measures. This parameter will be ignored if similarity measure cosine or shannon is chosen. Default: 1.1
    --normalization_method: Method used to normalize the intensities of each spectrum so that the intensities sum to 1. Since the objects entropy quantifies the uncertainy of must be probability distributions, the intensities of a given spectrum must sum to 1 prior to computing the entropy of the given spectrum intensities. Options: \'standard\' and \'softmax\'. Default: standard.
    --n_top_matches_to_save: The number of top matches to report. For example, if n_top_matches_to_save=5, then for each query spectrum, the five reference spectra with the largest similarity with the given query spectrum will be reported. Default: 1
    --print_id_results: Flag that prints identification results if True. Default: False
    --output_identification: Output CSV file containing the most-similar reference spectra for each query spectrum along with the corresponding similarity scores. Default is to save identification output in current working directory with filename \'output_identification.csv\'.
    --output_similarity_scores: Output CSV file containing similarity scores between all query spectrum/spectra and all reference spectra. Each row corresponds to a query spectrum, the left-most column contains the query spectrum/spectra identifier, and the remaining column contain the similarity scores with respect to all reference library spectra. If no argument passed, then this CSV file is written to the current working directory with filename \'output_all_similarity_scores\'.csv.')
    '''

    # load query and reference libraries
    if query_data is None:
        print('\nError: No argument passed to the mandatory query_data. Please pass the path to the CSV file of the query data.')
        sys.exit()
    else:
        extension = query_data.rsplit('.',1)
        extension = extension[(len(extension)-1)]
        if extension == 'mgf' or extension == 'MGF' or extension == 'mzML' or extension == 'mzml' or extension == 'MZML' or extension == 'cdf' or extension == 'CDF':
            output_path_tmp = query_data[:-3] + 'csv'
            build_library_from_raw_data(input_path=query_data, output_path=output_path_tmp, is_reference=False)
            df_query = pd.read_csv(output_path_tmp)
        if extension == 'csv' or extension == 'CSV':
            df_query = pd.read_csv(query_data)
        unique_query_ids = df_query.iloc[:,0].unique()

    if reference_data is None:
        print('\nError: No argument passed to the mandatory reference_data. Please pass the path to the CSV file of the reference data.')
        sys.exit()
    else:
        extension = reference_data.rsplit('.',1)
        extension = extension[(len(extension)-1)]
        if extension == 'mgf' or extension == 'MGF' or extension == 'mzML' or extension == 'mzml' or extension == 'MZML' or extension == 'cdf' or extension == 'CDF':
            output_path_tmp = reference_data[:-3] + 'csv'
            build_library_from_raw_data(input_path=reference_data, output_path=output_path_tmp, is_reference=True)
            df_reference = pd.read_csv(output_path_tmp)
        if extension == 'csv' or extension == 'CSV':
            df_reference = pd.read_csv(reference_data)
            if likely_reference_IDs is not None:
                likely_reference_IDs = pd.read_csv(likely_reference_IDs, header=None)
                df_reference = df_reference.loc[df_reference.iloc[:,0].isin(likely_reference_IDs.iloc[:,0].tolist())]
            unique_reference_ids = df_reference.iloc[:,0].unique()


    ##### process input parameters and ensure they are in a valid format #####
    if spectrum_preprocessing_order is not None:
        spectrum_preprocessing_order = list(spectrum_preprocessing_order)
    else:
        spectrum_preprocessing_order = ['F','N','W','L']
    if set(spectrum_preprocessing_order) - {'F','N','W','L'}:
        print(f'Error: spectrum_preprocessing_order must contain only \'F\', \'N\', \'W\', \'L\'.')
        sys.exit()

    if similarity_measure not in ['cosine','shannon','renyi','tsallis']:
        print('\nError: similarity_measure must be either \'cosine\', \'shannon\', \'renyi\', or \'tsallis\'')
        sys.exit()

    if isinstance(int_min,int) is True:
        int_min = float(int_min)
    if isinstance(int_max,int) is True:
        int_max = float(int_max)
    if isinstance(mz_min,int) is False or isinstance(mz_max,int) is False or isinstance(int_min,float) is False or isinstance(int_max,float) is False:
        print('Error: mz_min must be a non-negative integer, mz_max must be a positive integer, int_min must be a non-negative float, and int_max must be a positive float')
        sys.exit()
    if mz_min < 0:
        print('\nError: mz_min should be a non-negative integer')
        sys.exit()
    if mz_max <= 0:
        print('\nError: mz_max should be a positive integer')
        sys.exit()
    if int_min < 0:
        print('\nError: int_min should be a non-negative float')
        sys.exit()
    if int_max <= 0:
        print('\nError: int_max should be a positive float')
        sys.exit()

    if isinstance(noise_threshold,int) is True:
        noise_threshold = float(noise_threshold)
    if isinstance(noise_threshold,float) is False or noise_threshold < 0:
        print('Error: noise_threshold must be a positive float.')
        sys.exit()

    if isinstance(wf_intensity,int) is True:
        wf_intensity = float(wf_intensity)
    if isinstance(wf_mz,int) is True:
        wf_mz = float(wf_mz)
    if isinstance(wf_intensity,float) is False or isinstance(wf_mz,float) is False:
        print('Error: wf_mz and wf_intensity must be integers or floats')
        sys.exit()

    if entropy_dimension <= 0:
        print('\nError: entropy_dimension should be a positive float')
        sys.exit()
    else:
        q = entropy_dimension

    if normalization_method not in ['softmax','standard']:
        print('\nError: normalization_method must be either \'softmax\' or \'standard\'')
        sys.exit()

    if n_top_matches_to_save <= 0 or isinstance(n_top_matches_to_save,int)==False:
        print('\nError: n_top_matches_to_save should be a positive integer')
        sys.exit()

    if isinstance(print_id_results,bool)==False:
        print('\nError: print_id_results must be either True or False')
        sys.exit()
    
    if output_identification is None:
        output_identification = f'{Path.cwd()}/output_identification.csv'
        print(f'Warning: writing identification output to {output_identification}')

    if output_similarity_scores is None:
        output_similarity_scores = f'{Path.cwd()}/output_all_similarity_scores.csv'
        print(f'Warning: writing similarity scores to {output_similarity_scores}')



    ####################################### begin spectral library matching #######################################
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
    if print_id_results == True:
        print(df_top_ref_specs.to_string())

    # write spectral library matching results to disk
    df_top_ref_specs.to_csv(output_identification)

    # write all similarity scores to disk
    df_scores.columns = ['Reference Spectrum ID: ' + col for col in  list(map(str,df_scores.columns.tolist()))]
    df_scores.to_csv(output_similarity_scores)


