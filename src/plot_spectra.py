
# this script plots a given query spectrum against a given reference spectrum before and after spectrum preprocessing transformations

# load libraries
from processing import *
from similarity_measures import *
import pandas as pd
import argparse
from pathlib import Path
import sys
import matplotlib.pyplot as plt


# create argparse object so command-line input can be imported
parser = argparse.ArgumentParser()

# import optional command-line arguments
parser.add_argument('--query_data', type=str, metavar='\b', help='CSV file of query mass spectrum/spectra to be identified. Each row should correspond to a mass spectrum, the left-most column should contain an identifier, and each of the other columns should correspond to a single mass/charge ratio. Mandatory argument.')
parser.add_argument('--reference_data', type=str, metavar='\b', help='CSV file of the reference mass spectra. Each row should correspond to a mass spectrum, the left-most column should contain in identifier (i.e. the CAS registry number or the compound name), and the remaining column should correspond to a single mass/charge ratio. Mandatory argument.')
parser.add_argument('--spectrum_ID1', type=str, metavar='\b', help='The identifier of one (probably query) spectrum to be plotted. Default: first query spectrum in query_data.')
parser.add_argument('--spectrum_ID2', type=str, metavar='\b', help='The identifier of one (probably reference) spectrum to be plotted. Default: first reference spectrum in reference_data.')
parser.add_argument('--similarity_measure', type=str, default='cosine', metavar='\b', help='Similarity measure: options are \'cosine\', \'shannon\', \'renyi\', and \'tsallis\'. Default: cosine.')
parser.add_argument('--chromatography_platform', type=str, metavar='\b', help='Chromatography platform: options are \'HRMS\' and \'NRMS\'. Mandatory argument.')
parser.add_argument('--spectrum_preprocessing_order', type=str, metavar='\b', help='The spectrum preprocessing transformations and the order in which they are to be applied. Note that these transformations are applied prior to computing similarity scores. Format must be a string with 2-6 characters chosen from C, F, M, N, L, W representing centroiding, filtering based on mass/charge and intensity values, matching, noise removal, low-entropy transformation, and weight-factor-transformation, respectively. For example, if \'WCM\' is passed, then each spectrum will undergo a weight factor transformation, then centroiding, and then matching. Note that if an argument is passed, then \'M\' must be contained in the argument, since matching is a required preprocessing step in spectral library matching of LC-MS/MS data. Furthermore, \'C\' must be performed before matching since centroiding can change the number of ion fragments in a given spectrum. Default: FCNMWL for HRMS, FNLW for NRMS')
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
parser.add_argument('--y_axis_transformation', type=str, default='normalized', metavar='\b', help='Transformation to apply to y-axis (i.e. intensity axis) of plots. Options: \'normalized\', \'none\', \'log10\', and \'sqrt\'. Default: normalized.')
parser.add_argument('--save_plots', type=str, metavar='\b', help='Output PDF file containing the plots of the spectra before and after preprocessing transformations. If no argument is passed, then the plots will be saved to the PDF ./spectrum1_{spectrum_ID1}_spectrum2_{spectrum_ID2}_plot.pdf in the current working directory.')

# parse the user-input arguments
args = parser.parse_args()


# import the query library
if args.query_data is not None:
    df_query = pd.read_csv(args.query_data)
else:
    print('Error: No argument passed to query_data. To view usage, run \"python plot_spectra.py -h\".')
    sys.exit()


# load the reference library
if args.reference_data is not None:
    df_reference = pd.read_csv(args.reference_data)
else:
    print('Error: No argument passed to reference_data. To view usage, run \"python plot_spectra.py -h\".')
    sys.exit()


# import the identifier of the query spectrum to be plotted
if args.spectrum_ID1 is not None:
    spectrum_ID1 = str(args.spectrum_ID1)
else:
    spectrum_ID1 = str(df_query.iloc[0,0])
    print('No argument passed to spectrum_ID1; using the first spectrum in query_data.')


# import the identifier of the reference spectrum to be plotted
if args.spectrum_ID2 is not None:
    spectrum_ID2 = str(args.spectrum_ID2)
else:
    spectrum_ID2 = str(df_reference.iloc[0,0])
    print('No argument passed to spectrum_ID2; using the first spectrum in reference_data.')


# throw error if similarity measure is not one of the available similarity measures
similarity_measure = args.similarity_measure
if similarity_measure not in ['cosine','shannon','renyi','tsallis']:
    print('\nError: similarity_measure must be either \'cosine\', \'shannon\', \'tsallis\'')
    sys.exit()


# specify the chromatography platform
if args.chromatography_platform is not None:
    chromatography_platform = args.chromatography_platform
else:
    print('Error: No argument passed to chromatography_platform. To view usage, run \"python plot_spectra.py -h\".')
    sys.exit()

if chromatography_platform not in ['HRMS','NRMS']:
    print('\nError: chromatography_platform must be either \'HRMS\' or \'NRMS\'')
    sys.exit()


# ensure that the transformation to be applied to the y-axis (i.e. intensity axis), if any, if a valid option
y_axis_transformation = args.y_axis_transformation
if y_axis_transformation not in ['normalized','none','log10','sqrt']:
    print('Error: y_axis_transformation must be either \'normalized\', \'none\', \'log10\', or \'sqrt\'. To view usage, run \"python plot_spectra -h\".')
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


# load the path which the output PDF file should be written to
if args.save_plots is not None:
    path_output = args.save_plots
else:
    path_output = f'{Path.cwd()}/spectrum1_{spectrum_ID1}_spectrum2_{spectrum_ID2}_plot.pdf'


# create the figure
fig, axes = plt.subplots(nrows=2, ncols=1)

# note the workflow is slightly different for HRMS and NRMS data
if chromatography_platform == 'HRMS':

    # get unique query/reference library IDs; each query/reference ID corresponds to exactly one query/reference mass spectrum
    unique_query_ids = df_query.iloc[:,0].unique()
    unique_reference_ids = df_reference.iloc[:,0].unique()
    unique_query_ids = list(map(str,unique_query_ids))
    unique_reference_ids = list(map(str,unique_reference_ids))
    common_IDs = np.intersect1d(unique_query_ids,unique_reference_ids)
    if len(common_IDs) > 0:
        print(f'Error: the query and reference library have overlapping IDs: {common_IDs}')
        sys.exit()


    df_query = df_query.astype(object)
    df_reference = df_reference.astype(object)

    df_query.iloc[:,0] = df_query.iloc[:,0].astype(str)
    df_reference.iloc[:,0] = df_reference.iloc[:,0].astype(str)

    if spectrum_ID1 in unique_query_ids and spectrum_ID2 in unique_query_ids:
        query_idx = unique_query_ids.index(spectrum_ID1)
        reference_idx = unique_query_ids.index(spectrum_ID2)
        q_idxs_tmp = np.where(df_query.iloc[:,0] == unique_query_ids[query_idx])[0]
        r_idxs_tmp = np.where(df_query.iloc[:,0] == unique_query_ids[reference_idx])[0]
        q_spec = np.asarray(pd.concat([df_query.iloc[q_idxs_tmp,1], df_query.iloc[q_idxs_tmp,2]], axis=1).reset_index(drop=True))
        r_spec = np.asarray(pd.concat([df_query.iloc[r_idxs_tmp,1], df_query.iloc[r_idxs_tmp,2]], axis=1).reset_index(drop=True))
    elif spectrum_ID1 in unique_reference_ids and spectrum_ID2 in unique_reference_ids:
        query_idx = unique_reference_ids.index(spectrum_ID1)
        reference_idx = unique_reference_ids.index(spectrum_ID2)
        q_idxs_tmp = np.where(df_reference.iloc[:,0] == unique_reference_ids[query_idx])[0]
        r_idxs_tmp = np.where(df_reference.iloc[:,0] == unique_reference_ids[reference_idx])[0]
        q_spec = np.asarray(pd.concat([df_reference.iloc[q_idxs_tmp,1], df_reference.iloc[q_idxs_tmp,2]], axis=1).reset_index(drop=True))
        r_spec = np.asarray(pd.concat([df_reference.iloc[r_idxs_tmp,1], df_reference.iloc[r_idxs_tmp,2]], axis=1).reset_index(drop=True))
    else:
        if spectrum_ID1 in unique_reference_ids and spectrum_ID2 in unique_query_ids:
            spec_tmp = spectrum_ID1
            spectrum_ID1 = spectrum_ID2
            spectrum_ID2 = spec_tmp
        query_idx = unique_query_ids.index(spectrum_ID1)
        reference_idx = unique_reference_ids.index(spectrum_ID2)
        q_idxs_tmp = np.where(df_query.iloc[:,0] == unique_query_ids[query_idx])[0]
        r_idxs_tmp = np.where(df_reference.iloc[:,0] == unique_reference_ids[reference_idx])[0]
        q_spec = np.asarray(pd.concat([df_query.iloc[q_idxs_tmp,1], df_query.iloc[q_idxs_tmp,2]], axis=1).reset_index(drop=True))
        r_spec = np.asarray(pd.concat([df_reference.iloc[r_idxs_tmp,1], df_reference.iloc[r_idxs_tmp,2]], axis=1).reset_index(drop=True))


    q_spec_pre_trans = q_spec.copy()
    r_spec_pre_trans = r_spec.copy()
    q_spec_pre_trans[:,1] = q_spec_pre_trans[:,1].astype(float)
    r_spec_pre_trans[:,1] = r_spec_pre_trans[:,1].astype(float)

    # apply transformation to y-axis if relevant
    if y_axis_transformation == 'normalized':
        q_spec_pre_trans[:,1] = q_spec_pre_trans[:,1] / np.max(q_spec_pre_trans[:,1])
        r_spec_pre_trans[:,1] = r_spec_pre_trans[:,1] / np.max(r_spec_pre_trans[:,1])
    elif y_axis_transformation == 'log10':
        q_spec_pre_trans[:,1] = np.log10(np.array(q_spec_pre_trans[:,1]+1,dtype=float))
        r_spec_pre_trans[:,1] = np.log10(np.array(r_spec_pre_trans[:,1]+1,dtype=float))
    elif y_axis_transformation == 'sqrt':
        q_spec_pre_trans[:,1] = np.sqrt(np.array(q_spec_pre_trans[:,1],dtype=float))
        r_spec_pre_trans[:,1] = np.sqrt(np.array(r_spec_pre_trans[:,1],dtype=float))

    # plot the untransformed spectra
    plt.subplot(2,1,1)
    plt.vlines(x=q_spec_pre_trans[:,0], ymin=[0]*q_spec_pre_trans.shape[0], ymax=q_spec_pre_trans[:,1], linewidth=3, color='blue', label=f'Spectrum ID 1: {spectrum_ID1}')
    plt.vlines(x=r_spec_pre_trans[:,0], ymin=[0]*r_spec_pre_trans.shape[0], ymax=-r_spec_pre_trans[:,1], linewidth=3, color='red', label=f'Spectrum ID 2: {spectrum_ID2}')
    plt.xlabel('m/z',fontsize=8)
    plt.ylabel('Intensity', fontsize=8)
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=8)
    plt.title('Untransformed Spectra', fontsize=12)


    # perform the spectrum preprocessing transformations in the order specified
    is_matched = False
    for transformation in spectrum_preprocessing_order:
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
        if transformation == 'L' and q_spec.shape[0] > 1 and r_spec.shape[1] > 1: # low-entropy transformation
            q_spec[:,1] = LE_transform(q_spec[:,1], LET_threshold, normalization_method=normalization_method)
            r_spec[:,1] = LE_transform(r_spec[:,1], LET_threshold, normalization_method=normalization_method)
        if transformation == 'N' and q_spec.shape[0] > 1 and r_spec.shape[1] > 1: # noise removal
            q_spec = remove_noise(q_spec, nr = noise_threshold)
            r_spec = remove_noise(r_spec, nr = noise_threshold)
        if transformation == 'F' and q_spec.shape[0] > 1 and r_spec.shape[1] > 1: # filtering
            q_spec = filter_spec_lcms(q_spec, mz_min = mz_min, mz_max = mz_max, int_min = int_min, int_max = int_max, is_matched = is_matched)
            r_spec = filter_spec_lcms(r_spec, mz_min = mz_min, mz_max = mz_max, int_min = int_min, int_max = int_max, is_matched = is_matched)

    # intensities of query and reference library
    q_ints = q_spec[:,1]
    r_ints = r_spec[:,1]

    # if there is at least one non-zero intensity ion fragment in either spectra, compute their similarity
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

    # plot the transformed spectra
    plt.subplot(2,1,2)

    # display warning message if either spectra are empty or have no non-zero intensity ion fragments
    if q_spec.shape[0] > 1:
        if np.max(q_spec[:,1]) == 0 or np.max(r_spec[:,1]) == 0:
            plt.text(0.5, 0.5, 'The query and/or reference spectrum has no non-zero intensities after transformations.\n Change transformation parameters.', ha='center', va='center', fontsize=7, color='black')
            plt.xticks([])
            plt.yticks([])
        else:
            # apply transformation to y-axis if relevant
            if y_axis_transformation == 'normalized':
                q_spec[:,1] = q_spec[:,1] / np.max(q_spec[:,1])
                r_spec[:,1] = r_spec[:,1] / np.max(r_spec[:,1])
            elif y_axis_transformation == 'log10':
                q_spec[:,1] = np.log10(q_spec[:,1]+1)
                r_spec[:,1] = np.log10(r_spec[:,1]+1)
            elif y_axis_transformation == 'sqrt':
                q_spec[:,1] = np.sqrt(q_spec[:,1])
                r_spec[:,1] = np.sqrt(r_spec[:,1])
            plt.vlines(x=q_spec[:,0], ymin=[0]*q_spec.shape[0], ymax=q_spec[:,1], linewidth=3, color='blue')
            plt.vlines(x=r_spec[:,0], ymin=[0]*r_spec.shape[0], ymax=-r_spec[:,1], linewidth=3, color='red')
            plt.xlabel('m/z', fontsize=8)
            plt.ylabel('Intensity', fontsize=8)
            plt.xticks(fontsize=8)
            plt.yticks(fontsize=8)
            plt.title(f'Transformed Spectra\n Similarity Score: {round(similarity_score,4)}', fontsize=12)
    else:
        plt.text(0.5, 0.5, 'All points in the spectra were removed during preprocessing. \nChange the spectrum_preprocesing_order and/or change other spectrum-preprocessing parameters.', ha='center', va='center', fontsize=7, color='black')
        plt.xticks([])
        plt.yticks([])


elif chromatography_platform == 'NRMS':

    # get m/z values
    min_mz = np.min([np.min(df_query.iloc[:,1]), np.min(df_reference.iloc[:,1])])
    max_mz = np.max([np.max(df_query.iloc[:,1]), np.max(df_reference.iloc[:,1])])
    mzs = np.linspace(min_mz,max_mz,(max_mz-min_mz+1))

    # get unique query/reference library IDs; each query/reference ID corresponds to exactly one query/reference mass spectrum
    unique_query_ids = df_query.iloc[:,0].unique().tolist()
    unique_reference_ids = df_reference.iloc[:,0].unique().tolist()
    unique_query_ids = [str(ID) for ID in unique_query_ids]
    unique_reference_ids = [str(ID) for ID in unique_reference_ids]
    common_IDs = np.intersect1d([str(ID) for ID in unique_query_ids], [str(ID) for ID in unique_reference_ids])
    if len(common_IDs) > 0:
        print(f'Error: the query and reference library have overlapping IDs: {common_IDs}')
        sys.exit()

    if spectrum_ID1 in unique_query_ids and spectrum_ID2 in unique_query_ids:
        q_idxs_tmp = np.where(df_query.iloc[:,0].astype(str) == spectrum_ID1)[0]
        r_idxs_tmp = np.where(df_query.iloc[:,0].astype(str) == spectrum_ID2)[0]
        q_spec = np.asarray(pd.concat([df_query.iloc[q_idxs_tmp,1], df_query.iloc[q_idxs_tmp,2]], axis=1).reset_index(drop=True))
        r_spec = np.asarray(pd.concat([df_query.iloc[r_idxs_tmp,1], df_query.iloc[r_idxs_tmp,2]], axis=1).reset_index(drop=True))
    elif spectrum_ID1 in unique_reference_ids and spectrum_ID2 in unique_reference_ids:
        q_idxs_tmp = np.where(df_reference.iloc[:,0].astype(str) == spectrum_ID1)[0]
        r_idxs_tmp = np.where(df_reference.iloc[:,0].astype(str) == spectrum_ID2)[0]
        q_spec = np.asarray(pd.concat([df_reference.iloc[q_idxs_tmp,1], df_reference.iloc[q_idxs_tmp,2]], axis=1).reset_index(drop=True))
        r_spec = np.asarray(pd.concat([df_reference.iloc[r_idxs_tmp,1], df_reference.iloc[r_idxs_tmp,2]], axis=1).reset_index(drop=True))
    else:
        if spectrum_ID1 in unique_reference_ids and spectrum_ID2 in unique_query_ids:
            spec_tmp = spectrum_ID1
            spectrum_ID1 = spectrum_ID2
            spectrum_ID2 = spec_tmp
        q_idxs_tmp = np.where(df_query.iloc[:,0].astype(str) == spectrum_ID1)[0]
        r_idxs_tmp = np.where(df_reference.iloc[:,0].astype(str) == spectrum_ID2)[0]
        q_spec = np.asarray(pd.concat([df_query.iloc[q_idxs_tmp,1], df_query.iloc[q_idxs_tmp,2]], axis=1).reset_index(drop=True))
        r_spec = np.asarray(pd.concat([df_reference.iloc[r_idxs_tmp,1], df_reference.iloc[r_idxs_tmp,2]], axis=1).reset_index(drop=True))

    #print(q_spec)
    #print(r_spec)
    q_spec = convert_spec(q_spec,mzs)
    r_spec = convert_spec(r_spec,mzs)
    

    # plot the untransformed spectra
    plt.subplot(2,1,1)

    # display warning message if either spectra have no non-zero ion fragments
    if np.max(q_spec[:,1]) == 0 or np.max(r_spec[:,1]) == 0:
        plt.text(0.5, 0.5, 'The query and/or reference spectrum has no non-zero intensities after transformations.\n Change transformation parameters.', ha='center', va='center', fontsize=7, color='black')
        plt.xticks([])
        plt.yticks([])
    else:
        q_spec_pre_trans = q_spec.copy()
        r_spec_pre_trans = r_spec.copy()
        q_spec_pre_trans[:,1] = q_spec_pre_trans[:,1].astype(float)
        r_spec_pre_trans[:,1] = r_spec_pre_trans[:,1].astype(float)

        # apply transformation to y-axis if relevant
        if y_axis_transformation == 'normalized':
            q_spec_pre_trans[:,1] = q_spec_pre_trans[:,1] / np.max(q_spec_pre_trans[:,1])
            r_spec_pre_trans[:,1] = r_spec_pre_trans[:,1] / np.max(r_spec_pre_trans[:,1])
        elif y_axis_transformation == 'log10':
            q_spec_pre_trans[:,1] = np.log10(q_spec_pre_trans[:,1]+1)
            r_spec_pre_trans[:,1] = np.log10(r_spec_pre_trans[:,1]+1)
        elif y_axis_transformation == 'sqrt':
            q_spec_pre_trans[:,1] = np.sqrt(q_spec_pre_trans[:,1])
            r_spec_pre_trans[:,1] = np.sqrt(r_spec_pre_trans[:,1])
        plt.vlines(x=q_spec_pre_trans[:,0], ymin=[0]*len(q_spec_pre_trans[:,0]), ymax=q_spec_pre_trans[:,1], linewidth=3, color='blue', label=f'Spectrum ID1: {spectrum_ID1}')
        plt.vlines(x=r_spec_pre_trans[:,0], ymin=[0]*len(r_spec_pre_trans[:,0]), ymax=-r_spec_pre_trans[:,1], linewidth=3, color='red', label=f'Spectrum ID2: {spectrum_ID2}')
        plt.xlabel('m/z',fontsize=8)
        plt.ylabel('Intensity', fontsize=8)
        plt.xticks(fontsize=8)
        plt.yticks(fontsize=8)
        plt.title('Untransformed Query and Reference Spectra', fontsize=12)

    for transformation in spectrum_preprocessing_order:
        if transformation == 'W': # weight factor transformation
            q_spec[:,1] = wf_transform(q_spec[:,0], q_spec[:,1], wf_mz, wf_intensity)
            r_spec[:,1] = wf_transform(r_spec[:,0], r_spec[:,1], wf_mz, wf_intensity)
        if transformation == 'L': # low-entropy transformation
            q_spec[:,1] = LE_transform(q_spec[:,1], LET_threshold, normalization_method)
            r_spec[:,1] = LE_transform(r_spec[:,1], LET_threshold, normalization_method)
        if transformation == 'N': # noise removal
            q_spec = remove_noise(q_spec, nr = noise_threshold)
            if high_quality_reference_library == False:
                r_spec = remove_noise(r_spec, nr = noise_threshold)
        if transformation == 'F': # filtering with respect to mz and/or intensity
            q_spec = filter_spec_gcms(q_spec, mz_min = mz_min, mz_max = mz_max, int_min = int_min, int_max = int_max)
            if high_quality_reference_library == False:
                r_spec = filter_spec_gcms(r_spec, mz_min = mz_min, mz_max = mz_max, int_min = int_min, int_max = int_max)

    # compute similarity score; if the spectra contain one point at most, their similarity is considered to be 0
    if q_spec.shape[0] > 1:
        if similarity_measure == 'cosine':
            similarity_score = S_cos(q_spec[:,1], r_spec[:,1])
        else:
            q_spec[:,1] = normalize(q_spec[:,1], method = normalization_method)
            r_spec[:,1] = normalize(r_spec[:,1], method = normalization_method)

            if similarity_measure == 'shannon':
                similarity_score = S_shannon(q_spec[:,1].astype('float'), r_spec[:,1].astype('float'))
            elif similarity_measure == 'renyi':
                similarity_score = S_renyi(q_spec[:,1], r_spec[:,1], q)
            elif similarity_measure == 'tsallis':
                similarity_score = S_tsallis(q_spec[:,1], r_spec[:,1], q)
    else:
        similarity_score = 0

 
    # plot the transformed spectra
    plt.subplot(2,1,2)

    # display warning message if either spectra are empty or have no non-zero intensity ion fragments
    if q_spec.shape[0] == 0 or r_spec.shape[0] == 0:
        plt.text(0.5, 0.5, 'The query and/or reference spectrum has no ion fragments left after transformations.\n Change transformation parameters.', ha='center', va='center', fontsize=7, color='black')
        plt.xticks([])
        plt.yticks([])
    elif np.max(q_spec[:,1]) == 0 or np.max(r_spec[:,1]) == 0:
        plt.text(0.5, 0.5, 'The query and/or reference spectrum has no non-zero intensities after transformations.\n Change transformation parameters.', ha='center', va='center', fontsize=7, color='black')
        plt.xticks([])
        plt.yticks([])
    else:
        # apply transformation to y-axis if relevant
        if y_axis_transformation == 'normalized':
            q_spec[:,1] = q_spec[:,1] / np.max(q_spec[:,1])
            r_spec[:,1] = r_spec[:,1] / np.max(r_spec[:,1])
        elif y_axis_transformation == 'log10':
            q_spec[:,1] = np.log10(q_spec[:,1]+1)
            r_spec[:,1] = np.log10(r_spec[:,1]+1)
        elif y_axis_transformation == 'sqrt':
            q_spec[:,1] = np.sqrt(q_spec[:,1])
            r_spec[:,1] = np.sqrt(r_spec[:,1])
        plt.vlines(x=mzs, ymin=[0]*len(mzs), ymax=q_spec[:,1], linewidth=3, color='blue')
        plt.vlines(x=mzs, ymin=[0]*len(mzs), ymax=-r_spec[:,1], linewidth=3, color='red')
        plt.xlabel('m/z', fontsize=8)
        plt.ylabel('Intensity', fontsize=8)
        plt.xticks(fontsize=8)
        plt.yticks(fontsize=8)
        plt.title(f'Transformed Query and Reference Spectra\n Similarity Score: {round(similarity_score,4)}', fontsize=12)



# adjust margins of figure
plt.subplots_adjust(top = 0.8, hspace = 0.7)

# include legend
plt.figlegend(loc = 'upper center')

# write figure to PDF
plt.savefig(path_output, format='pdf')


