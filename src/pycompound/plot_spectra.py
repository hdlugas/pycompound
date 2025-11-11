
from .processing import *
from .similarity_measures import *
import pandas as pd
from pathlib import Path
import sys
import matplotlib.pyplot as plt


def generate_plots_on_HRMS_data(query_data=None, reference_data=None, spectrum_ID1=None, spectrum_ID2=None, similarity_measure='cosine', weights={'Cosine':0.25,'Shannon':0.25,'Renyi':0.25,'Tsallis':0.25}, spectrum_preprocessing_order='FCNMWL', high_quality_reference_library=False, mz_min=0, mz_max=9999999, int_min=0, int_max=9999999, window_size_centroiding=0.5, window_size_matching=0.5, noise_threshold=0.0, wf_mz=0.0, wf_intensity=1.0, LET_threshold=0.0, entropy_dimension=1.1, y_axis_transformation='normalized', output_path=None, return_plot=False):
    '''
    plots two spectra against each other before and after preprocessing transformations for high-resolution mass spectrometry data

    --query_data: mgf, mzML, or csv file of query mass spectrum/spectra to be identified. If csv file, each row should correspond to a mass spectrum, the left-most column should contain an identifier, and each of the other columns should correspond to a single mass/charge ratio. Mandatory argument.
    --reference_data: mgf, mzML, or csv file of the reference mass spectra. If csv file, each row should correspond to a mass spectrum, the left-most column should contain in identifier (i.e. the CAS registry number or the compound name), and the remaining column should correspond to a single mass/charge ratio. Mandatory argument.
    --spectrum_ID1: ID of one spectrum to be plotted. Default is first spectrum in the query library. Optional argument.
    --spectrum_ID2: ID of another spectrum to be plotted. Default is first spectrum in the reference library. Optional argument.
    --similarity_measure: cosine, shannon, renyi, tsallis, mixture, jaccard, dice, 3w_jaccard, sokal_sneath, binary_cosine, mountford, mcconnaughey, driver_kroeber, simpson, braun_banquet, fager_mcgowan, kulczynski, intersection, hamming, hellinger. Default: cosine.
    --weights: dict of weights to give to each non-binary similarity measure (i.e. cosine, shannon, renyi, and tsallis) when the mixture similarity measure is specified. Default: 0.25 for each of the four non-binary similarity measures.
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
    --y_axis_transformation: transformation to apply to y-axis (i.e. intensity axis) of plots. Options: \'normalized\', \'none\', \'log10\', and \'sqrt\'. Default: normalized.')
    --output_path: path to output PDF file containing the plots of the spectra before and after preprocessing transformations. If no argument is passed, then the plots will be saved to the PDF ./spectrum1_{spectrum_ID1}_spectrum2_{spectrum_ID2}_plot.pdf in the current working directory.
    '''

    if query_data is None:
        print('\nError: No argument passed to the mandatory query_data. Please pass the path to the CSV file of the query data.')
        sys.exit()
    else:
        extension = query_data.rsplit('.',1)
        extension = extension[(len(extension)-1)]
        if extension == 'mgf' or extension == 'MGF' or extension == 'mzML' or extension == 'mzml' or extension == 'MZML' or extension == 'cdf' or extension == 'CDF':
            output_path_tmp = query_data[:-3] + 'csv'
            build_library_from_raw_data(input_path=query_data, output_path=output_path_tmp, is_reference=True)
            df_query = pd.read_csv(output_path_tmp)
        if extension == 'csv' or extension == 'CSV':
            df_query = pd.read_csv(query_data)
        unique_query_ids = df_query.iloc[:,0].unique().tolist()
        unique_query_ids = [str(tmp) for tmp in unique_query_ids]

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
        unique_reference_ids = df_reference.iloc[:,0].unique().tolist()
        unique_reference_ids = [str(tmp) for tmp in unique_reference_ids]


    if spectrum_ID1 is not None:
        spectrum_ID1 = str(spectrum_ID1)
    else:
        spectrum_ID1 = str(df_query.iloc[0,0])
        print('No argument passed to spectrum_ID1; using the first spectrum in query_data.')

    if spectrum_ID2 is not None:
        spectrum_ID2 = str(spectrum_ID2)
    else:
        spectrum_ID2 = str(df_reference.iloc[0,0])
        print('No argument passed to spectrum_ID2; using the first spectrum in reference_data.')

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

    if similarity_measure not in ['cosine','shannon','renyi','tsallis','mixture','jaccard','dice','3w_jaccard','sokal_sneath','binary_cosine','mountford','mcconnaughey','driver_kroeber','simpson','braun_banquet','fager_mcgowan','kulczynski','intersection','hamming','hellinger']:
        print('\nError: similarity_measure must be either cosine, shannon, renyi, tsallis, mixture, jaccard, dice, 3w_jaccard, sokal_sneath, binary_cosine, mountford, mcconnaughey, driver_kroeber, simpson, braun_banquet, fager_mcgowan, kulczynski, intersection, hamming, or hellinger.')
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

    normalization_method = 'standard' #consider including additional normalization methods to transform intensities into a probability distribution; softmax results in many numerical errors/warnings

    if y_axis_transformation not in ['normalized','none','log10','sqrt']:
        print('Error: y_axis_transformation must be either \'normalized\', \'none\', \'log10\', or \'sqrt\'.')
        sys.exit()

    if output_path is None:
        print(f'Warning: plots will be saved to the PDF ./spectrum1_{spectrum_ID1}_spectrum2_{spectrum_ID2}_plot.pdf in the current working directory.')
        output_path = f'{Path.cwd()}/spectrum1_{spectrum_ID1}_spectrum2_{spectrum_ID2}.pdf'


    if spectrum_ID1 in unique_query_ids and spectrum_ID2 in unique_query_ids:
        query_idx = unique_query_ids.index(spectrum_ID1)
        reference_idx = unique_query_ids.index(spectrum_ID2)
        q_idxs_tmp = np.where(df_query.iloc[:,0].astype(str) == unique_query_ids[query_idx])[0]
        r_idxs_tmp = np.where(df_query.iloc[:,0].astype(str) == unique_query_ids[reference_idx])[0]
        q_spec = np.asarray(pd.concat([df_query.iloc[q_idxs_tmp,1], df_query.iloc[q_idxs_tmp,2]], axis=1).reset_index(drop=True))
        r_spec = np.asarray(pd.concat([df_query.iloc[r_idxs_tmp,1], df_query.iloc[r_idxs_tmp,2]], axis=1).reset_index(drop=True))
    elif spectrum_ID1 in unique_reference_ids and spectrum_ID2 in unique_reference_ids:
        query_idx = unique_reference_ids.index(spectrum_ID1)
        reference_idx = unique_reference_ids.index(spectrum_ID2)
        q_idxs_tmp = np.where(df_reference.iloc[:,0].astype(str) == unique_reference_ids[query_idx])[0]
        r_idxs_tmp = np.where(df_reference.iloc[:,0].astype(str) == unique_reference_ids[reference_idx])[0]
        q_spec = np.asarray(pd.concat([df_reference.iloc[q_idxs_tmp,1], df_reference.iloc[q_idxs_tmp,2]], axis=1).reset_index(drop=True))
        r_spec = np.asarray(pd.concat([df_reference.iloc[r_idxs_tmp,1], df_reference.iloc[r_idxs_tmp,2]], axis=1).reset_index(drop=True))
    else:
        if spectrum_ID1 in unique_reference_ids and spectrum_ID2 in unique_query_ids:
            spec_tmp = spectrum_ID1
            spectrum_ID1 = spectrum_ID2
            spectrum_ID2 = spec_tmp
        query_idx = unique_query_ids.index(spectrum_ID1)
        reference_idx = unique_reference_ids.index(spectrum_ID2)
        q_idxs_tmp = np.where(df_query.iloc[:,0].astype(str) == unique_query_ids[query_idx])[0]
        r_idxs_tmp = np.where(df_reference.iloc[:,0].astype(str) == unique_reference_ids[reference_idx])[0]
        q_spec = np.asarray(pd.concat([df_query.iloc[q_idxs_tmp,1], df_query.iloc[q_idxs_tmp,2]], axis=1).reset_index(drop=True))
        r_spec = np.asarray(pd.concat([df_reference.iloc[r_idxs_tmp,1], df_reference.iloc[r_idxs_tmp,2]], axis=1).reset_index(drop=True))


    q_spec_pre_trans = q_spec.copy()
    r_spec_pre_trans = r_spec.copy()
    q_spec_pre_trans[:,1] = q_spec_pre_trans[:,1].astype(float)
    r_spec_pre_trans[:,1] = r_spec_pre_trans[:,1].astype(float)

    if y_axis_transformation == 'normalized':
        q_spec_pre_trans[:,1] = q_spec_pre_trans[:,1] / np.max(q_spec_pre_trans[:,1])
        r_spec_pre_trans[:,1] = r_spec_pre_trans[:,1] / np.max(r_spec_pre_trans[:,1])
        ylab = 'Normalized Intensity'
    elif y_axis_transformation == 'log10':
        q_spec_pre_trans[:,1] = np.log10(np.array(q_spec_pre_trans[:,1]+1,dtype=float))
        r_spec_pre_trans[:,1] = np.log10(np.array(r_spec_pre_trans[:,1]+1,dtype=float))
        ylab = 'log10(Intensity)'
    elif y_axis_transformation == 'sqrt':
        q_spec_pre_trans[:,1] = np.sqrt(np.array(q_spec_pre_trans[:,1],dtype=float))
        r_spec_pre_trans[:,1] = np.sqrt(np.array(r_spec_pre_trans[:,1],dtype=float))
        ylab = 'sqrt(Intensity)'
    else:
        ylab = 'Raw Intensity'

    fig, axes = plt.subplots(nrows=2, ncols=1)

    plt.subplot(2,1,1)
    plt.vlines(x=q_spec_pre_trans[:,0], ymin=[0]*q_spec_pre_trans.shape[0], ymax=q_spec_pre_trans[:,1], linewidth=3, color='blue', label=f'Spectrum ID 1: {spectrum_ID1}')
    plt.vlines(x=r_spec_pre_trans[:,0], ymin=[0]*r_spec_pre_trans.shape[0], ymax=-r_spec_pre_trans[:,1], linewidth=3, color='red', label=f'Spectrum ID 2: {spectrum_ID2}')
    plt.xlabel('m/z',fontsize=7)
    plt.ylabel(ylab, fontsize=7)
    plt.xticks(fontsize=7)
    plt.yticks(fontsize=7)
    plt.title('Untransformed Spectra', fontsize=10)

    mz_min_tmp_q = round(q_spec[:,0].min(),1)
    mz_min_tmp_r = round(r_spec[:,0].min(),1)
    int_min_tmp_q = round(q_spec[:,1].min(),1)
    int_min_tmp_r = round(r_spec[:,1].min(),1)
    mz_max_tmp_q = round(q_spec[:,0].max(),1)
    mz_max_tmp_r = round(r_spec[:,0].max(),1)
    int_max_tmp_q = round(q_spec[:,1].max(),1)
    int_max_tmp_r = round(r_spec[:,1].max(),1)
    mz_min_tmp = min([mz_min_tmp_q,mz_min_tmp_r])
    mz_max_tmp = min([mz_max_tmp_q,mz_max_tmp_r])
    int_min_tmp = min([int_min_tmp_q,int_min_tmp_r])
    int_max_tmp = max([int_max_tmp_q,int_max_tmp_r])

    is_matched = False
    for transformation in spectrum_preprocessing_order:
        if transformation == 'C' and q_spec.shape[0] > 1 and r_spec.shape[1] > 1:
            q_spec = centroid_spectrum(q_spec, window_size=window_size_centroiding) 
            r_spec = centroid_spectrum(r_spec, window_size=window_size_centroiding) 
        if transformation == 'M' and q_spec.shape[0] > 1 and r_spec.shape[1] > 1:
            m_spec = match_peaks_in_spectra(spec_a=q_spec, spec_b=r_spec, window_size=window_size_matching)
            q_spec = m_spec[:,0:2]
            r_spec = m_spec[:,[0,2]]
            is_matched = True
        if transformation == 'W' and q_spec.shape[0] > 1 and r_spec.shape[1] > 1:
            q_spec[:,1] = wf_transform(q_spec[:,0], q_spec[:,1], wf_mz, wf_intensity)
            r_spec[:,1] = wf_transform(r_spec[:,0], r_spec[:,1], wf_mz, wf_intensity)
        if transformation == 'L' and q_spec.shape[0] > 1 and r_spec.shape[1] > 1:
            q_spec[:,1] = LE_transform(q_spec[:,1], LET_threshold, normalization_method=normalization_method)
            r_spec[:,1] = LE_transform(r_spec[:,1], LET_threshold, normalization_method=normalization_method)
        if transformation == 'N' and q_spec.shape[0] > 1 and r_spec.shape[1] > 1:
            q_spec = remove_noise(q_spec, nr = noise_threshold)
            if high_quality_reference_library == False or high_quality_reference_library == 'False':
                r_spec = remove_noise(r_spec, nr = noise_threshold)
        if transformation == 'F' and q_spec.shape[0] > 1 and r_spec.shape[1] > 1:
            q_spec = filter_spec_lcms(q_spec, mz_min = mz_min, mz_max = mz_max, int_min = int_min, int_max = int_max, is_matched = is_matched)
            if high_quality_reference_library == False or high_quality_reference_library == 'False':
                r_spec = filter_spec_lcms(r_spec, mz_min = mz_min, mz_max = mz_max, int_min = int_min, int_max = int_max, is_matched = is_matched)

    q_ints = q_spec[:,1]
    r_ints = r_spec[:,1]

    if np.sum(q_ints) != 0 and np.sum(r_ints) != 0 and q_spec.shape[0] > 1 and r_spec.shape[1] > 1:
        similarity_score = get_similarity(similarity_measure, q_ints, r_ints, weights, entropy_dimension)
    else:
        similarity_score = 0

    plt.subplot(2,1,2)

    if q_spec.shape[0] > 1:
        if np.max(q_spec[:,1]) == 0 or np.max(r_spec[:,1]) == 0:
            plt.text(0.5, 0.5, 'The query and/or reference spectrum has no non-zero intensities after transformations.\n Change transformation parameters.', ha='center', va='center', fontsize=7, color='black')
            plt.xticks([])
            plt.yticks([])
        else:
            if y_axis_transformation == 'normalized':
                q_spec[:,1] = q_spec[:,1] / np.max(q_spec[:,1])
                r_spec[:,1] = r_spec[:,1] / np.max(r_spec[:,1])
                ylab='Normalized Intensity'
            elif y_axis_transformation == 'log10':
                q_spec[:,1] = np.log10(q_spec[:,1]+1)
                r_spec[:,1] = np.log10(r_spec[:,1]+1)
                ylab='log10(Intensity)'
            elif y_axis_transformation == 'sqrt':
                q_spec[:,1] = np.sqrt(q_spec[:,1])
                r_spec[:,1] = np.sqrt(r_spec[:,1])
                ylab='sqrt(Intensity)'
            else:
                ylab = 'Raw Intensity'
            plt.vlines(x=q_spec[:,0], ymin=[0]*q_spec.shape[0], ymax=q_spec[:,1], linewidth=3, color='blue')
            plt.vlines(x=r_spec[:,0], ymin=[0]*r_spec.shape[0], ymax=-r_spec[:,1], linewidth=3, color='red')
            plt.xlabel('m/z', fontsize=7)
            plt.ylabel(ylab, fontsize=7)
            plt.xticks(fontsize=7)
            plt.yticks(fontsize=7)
            plt.title(f'Transformed Spectra', fontsize=10)
    else:
        plt.text(0.5, 0.5, 'All points in the spectra were removed during preprocessing. \nChange the spectrum_preprocesing_order and/or change other spectrum-preprocessing parameters.', ha='center', va='center', fontsize=7, color='black')
        plt.xticks([])
        plt.yticks([])


    print('\n\n\n')
    print(high_quality_reference_library)
    print('\n\n\n')
    plt.subplots_adjust(top=0.8, hspace=0.92, bottom=0.3)
    plt.figlegend(loc = 'upper center')
    fig.text(0.05, 0.18, f'Similarity Measure: {similarity_measure.capitalize()}', fontsize=7)
    fig.text(0.05, 0.15, f'Similarity Score: {round(similarity_score,4)}', fontsize=7)
    fig.text(0.05, 0.12, f"Spectrum Preprocessing Order: {''.join(spectrum_preprocessing_order)}", fontsize=7)
    fig.text(0.05, 0.09, f'High Quality Reference Library: {str(high_quality_reference_library)}', fontsize=7)
    fig.text(0.05, 0.06, f'Window Size (Centroiding): {window_size_centroiding}', fontsize=7)
    fig.text(0.05, 0.03, f'Window Size (Matching): {window_size_matching}', fontsize=7)
    fig.text(0.45, 0.18, f'Raw-Scale M/Z Range: [{mz_min_tmp},{mz_max_tmp}]', fontsize=7)
    fig.text(0.45, 0.15, f'Raw-Scale Intensity Range: [{int_min_tmp},{int_max_tmp}]', fontsize=7)
    fig.text(0.45, 0.12, f'Noise Threshold: {noise_threshold}', fontsize=7)
    fig.text(0.45, 0.09, f'Weight Factors (m/z,intensity): ({wf_mz},{wf_intensity})', fontsize=7)
    fig.text(0.45, 0.06, f'Low-Entropy Threshold: {LET_threshold}', fontsize=7)
    if similarity_measure == 'mixture':
        fig.text(0.45, 0.03, f'Weights for mixture similarity: {weights}', fontsize=7)

    plt.savefig(output_path, format='pdf')

    if return_plot == True:
        return fig 




def generate_plots_on_NRMS_data(query_data=None, reference_data=None, spectrum_ID1=None, spectrum_ID2=None, similarity_measure='cosine', weights={'Cosine':0.25,'Shannon':0.25,'Renyi':0.25,'Tsallis':0.25}, spectrum_preprocessing_order='FNLW', high_quality_reference_library=False, mz_min=0, mz_max=9999999, int_min=0, int_max=9999999, noise_threshold=0.0, wf_mz=0.0, wf_intensity=1.0, LET_threshold=0.0, entropy_dimension=1.1, y_axis_transformation='normalized', output_path=None, return_plot=False):
    '''
    plots two spectra against each other before and after preprocessing transformations for high-resolution mass spectrometry data

    --query_data: cdf or csv file of query mass spectrum/spectra to be identified. If csv file, each row should correspond to a mass spectrum, the left-most column should contain an identifier, and each of the other columns should correspond to a single mass/charge ratio. Mandatory argument.
    --reference_data: cdf of csv file of the reference mass spectra. If csv file, each row should correspond to a mass spectrum, the left-most column should contain in identifier (i.e. the CAS registry number or the compound name), and the remaining column should correspond to a single mass/charge ratio. Mandatory argument.
    --similarity_measure: cosine, shannon, renyi, tsallis, mixture, jaccard, dice, 3w_jaccard, sokal_sneath, binary_cosine, mountford, mcconnaughey, driver_kroeber, simpson, braun_banquet, fager_mcgowan, kulczynski, intersection, hamming, hellinger. Default: cosine.
    --weights: dict of weights to give to each non-binary similarity measure (i.e. cosine, shannon, renyi, and tsallis) when the mixture similarity measure is specified. Default: 0.25 for each of the four non-binary similarity measures.
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
    --y_axis_transformation: transformation to apply to y-axis (i.e. intensity axis) of plots. Options: \'normalized\', \'none\', \'log10\', and \'sqrt\'. Default: normalized.')
    --output_path: path to output PDF file containing the plots of the spectra before and after preprocessing transformations. If no argument is passed, then the plots will be saved to the PDF ./spectrum1_{spectrum_ID1}_spectrum2_{spectrum_ID2}_plot.pdf in the current working directory.
    '''

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
            unique_reference_ids = df_reference.iloc[:,0].unique()


    if spectrum_ID1 is not None:
        spectrum_ID1 = str(spectrum_ID1)
    else:
        spectrum_ID1 = str(df_query.iloc[0,0])
        print('No argument passed to spectrum_ID1; using the first spectrum in query_data.')

    if spectrum_ID2 is not None:
        spectrum_ID2 = str(spectrum_ID2)
    else:
        spectrum_ID2 = str(df_reference.iloc[0,0])
        print('No argument passed to spectrum_ID2; using the first spectrum in reference_data.')

    if spectrum_preprocessing_order is not None:
        spectrum_preprocessing_order = list(spectrum_preprocessing_order)
    else:
        spectrum_preprocessing_order = ['F','N','W','L']
    if set(spectrum_preprocessing_order) - {'F','N','W','L'}:
        print(f'Error: spectrum_preprocessing_order must contain only \'F\', \'N\', \'W\', \'L\'.')
        sys.exit()

    if similarity_measure not in ['cosine','shannon','renyi','tsallis','mixture','jaccard','dice','3w_jaccard','sokal_sneath','binary_cosine','mountford','mcconnaughey','driver_kroeber','simpson','braun_banquet','fager_mcgowan','kulczynski','intersection','hamming','hellinger']:
        print('\nError: similarity_measure must be either cosine, shannon, renyi, tsallis, mixture, jaccard, dice, 3w_jaccard, sokal_sneath, binary_cosine, mountford, mcconnaughey, driver_kroeber, simpson, braun_banquet, fager_mcgowan, kulczynski, intersection, hamming, or hellinger.')
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

    normalization_method = 'standard' #consider including additional normalization methods to transform intensities into a probability distribution; softmax results in many numerical errors/warnings

    if y_axis_transformation not in ['normalized','none','log10','sqrt']:
        print('Error: y_axis_transformation must be either \'normalized\', \'none\', \'log10\', or \'sqrt\'.')
        sys.exit()

    if output_path is None:
        print(f'Warning: plots will be saved to the PDF ./spectrum1_{spectrum_ID1}_spectrum2_{spectrum_ID2}_plot.pdf in the current working directory.')
        output_path = f'{Path.cwd()}/spectrum1_{spectrum_ID1}_spectrum2_{spectrum_ID2}.pdf'

    min_mz = np.min([np.min(df_query.iloc[:,1]), np.min(df_reference.iloc[:,1])])
    max_mz = np.max([np.max(df_query.iloc[:,1]), np.max(df_reference.iloc[:,1])])
    mzs = np.linspace(min_mz,max_mz,(max_mz-min_mz+1))

    unique_query_ids = df_query.iloc[:,0].unique().tolist()
    unique_reference_ids = df_reference.iloc[:,0].unique().tolist()
    unique_query_ids = [str(ID) for ID in unique_query_ids]
    unique_reference_ids = [str(ID) for ID in unique_reference_ids]
    common_IDs = np.intersect1d([str(ID) for ID in unique_query_ids], [str(ID) for ID in unique_reference_ids])
    if len(common_IDs) > 0:
        print(f'Warning: the query and reference library have overlapping IDs: {common_IDs}')

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

    q_spec = convert_spec(q_spec,mzs)
    r_spec = convert_spec(r_spec,mzs)

    int_min_tmp_q = min(q_spec[q_spec[:,1].nonzero(),1][0])
    int_min_tmp_r = min(r_spec[r_spec[:,1].nonzero(),1][0])
    int_max_tmp_q = max(q_spec[q_spec[:,1].nonzero(),1][0])
    int_max_tmp_r = max(r_spec[r_spec[:,1].nonzero(),1][0])
    int_min_tmp = int(min([int_min_tmp_q,int_min_tmp_r]))
    int_max_tmp = int(max([int_max_tmp_q,int_max_tmp_r]))
    
    fig, axes = plt.subplots(nrows=2, ncols=1)

    plt.subplot(2,1,1)

    if np.max(q_spec[:,1]) == 0 or np.max(r_spec[:,1]) == 0:
        plt.text(0.5, 0.5, 'The query and/or reference spectrum has no non-zero intensities after transformations.\n Change transformation parameters.', ha='center', va='center', fontsize=7, color='black')
        plt.xticks([])
        plt.yticks([])
    else:
        q_spec_pre_trans = q_spec.copy()
        r_spec_pre_trans = r_spec.copy()
        q_spec_pre_trans[:,1] = q_spec_pre_trans[:,1].astype(float)
        r_spec_pre_trans[:,1] = r_spec_pre_trans[:,1].astype(float)

        if y_axis_transformation == 'normalized':
            q_spec_pre_trans[:,1] = q_spec_pre_trans[:,1] / np.max(q_spec_pre_trans[:,1])
            r_spec_pre_trans[:,1] = r_spec_pre_trans[:,1] / np.max(r_spec_pre_trans[:,1])
            ylab = 'Normalized Intensity'
        elif y_axis_transformation == 'log10':
            q_spec_pre_trans[:,1] = np.log10(q_spec_pre_trans[:,1]+1)
            r_spec_pre_trans[:,1] = np.log10(r_spec_pre_trans[:,1]+1)
            ylab = 'log10(Intensity)'
        elif y_axis_transformation == 'sqrt':
            q_spec_pre_trans[:,1] = np.sqrt(q_spec_pre_trans[:,1])
            r_spec_pre_trans[:,1] = np.sqrt(r_spec_pre_trans[:,1])
            ylab = 'sqrt(Intensity)'
        else:
            ylab = 'Raw Intensity'
        plt.vlines(x=q_spec_pre_trans[:,0], ymin=[0]*len(q_spec_pre_trans[:,0]), ymax=q_spec_pre_trans[:,1], linewidth=3, color='blue', label=f'Spectrum ID1: {spectrum_ID1}')
        plt.vlines(x=r_spec_pre_trans[:,0], ymin=[0]*len(r_spec_pre_trans[:,0]), ymax=-r_spec_pre_trans[:,1], linewidth=3, color='red', label=f'Spectrum ID2: {spectrum_ID2}')
        plt.xlabel('m/z',fontsize=7)
        plt.ylabel(ylab, fontsize=7)
        plt.xticks(fontsize=7)
        plt.yticks(fontsize=7)
        plt.title('Untransformed Query and Reference Spectra', fontsize=10)

    for transformation in spectrum_preprocessing_order:
        if transformation == 'W':
            q_spec[:,1] = wf_transform(q_spec[:,0], q_spec[:,1], wf_mz, wf_intensity)
            r_spec[:,1] = wf_transform(r_spec[:,0], r_spec[:,1], wf_mz, wf_intensity)
        if transformation == 'L':
            q_spec[:,1] = LE_transform(q_spec[:,1], LET_threshold, normalization_method)
            r_spec[:,1] = LE_transform(r_spec[:,1], LET_threshold, normalization_method)
        if transformation == 'N':
            q_spec = remove_noise(q_spec, nr = noise_threshold)
            if high_quality_reference_library == False or high_quality_reference_library == 'False':
                r_spec = remove_noise(r_spec, nr = noise_threshold)
        if transformation == 'F':
            q_spec = filter_spec_gcms(q_spec, mz_min = mz_min, mz_max = mz_max, int_min = int_min, int_max = int_max)
            if high_quality_reference_library == False or high_quality_reference_library == 'False':
                r_spec = filter_spec_gcms(r_spec, mz_min = mz_min, mz_max = mz_max, int_min = int_min, int_max = int_max)

    if q_spec.shape[0] > 1:
        similarity_score = get_similarity(similarity_measure, q_spec[:,1], r_spec[:,1], weights, entropy_dimension)
    else:
        similarity_score = 0

 
    plt.subplot(2,1,2)

    if q_spec.shape[0] == 0 or r_spec.shape[0] == 0:
        plt.text(0.5, 0.5, 'The query and/or reference spectrum has no ion fragments left after transformations.\n Change transformation parameters.', ha='center', va='center', fontsize=7, color='black')
        plt.xticks([])
        plt.yticks([])
    elif np.max(q_spec[:,1]) == 0 or np.max(r_spec[:,1]) == 0:
        plt.text(0.5, 0.5, 'The query and/or reference spectrum has no non-zero intensities after transformations.\n Change transformation parameters.', ha='center', va='center', fontsize=7, color='black')
        plt.xticks([])
        plt.yticks([])
    else:
        if y_axis_transformation == 'normalized':
            q_spec[:,1] = q_spec[:,1] / np.max(q_spec[:,1])
            r_spec[:,1] = r_spec[:,1] / np.max(r_spec[:,1])
            ylab='Normalized Intensity'
        elif y_axis_transformation == 'log10':
            q_spec[:,1] = np.log10(q_spec[:,1]+1)
            r_spec[:,1] = np.log10(r_spec[:,1]+1)
            ylab='log10(Intensity)'
        elif y_axis_transformation == 'sqrt':
            q_spec[:,1] = np.sqrt(q_spec[:,1])
            r_spec[:,1] = np.sqrt(r_spec[:,1])
            ylab='sqrt(Intensity)'
        else:
            ylab = 'Raw Intensity'
        plt.vlines(x=mzs, ymin=[0]*len(mzs), ymax=q_spec[:,1], linewidth=3, color='blue')
        plt.vlines(x=mzs, ymin=[0]*len(mzs), ymax=-r_spec[:,1], linewidth=3, color='red')
        plt.xlabel('m/z', fontsize=7)
        plt.ylabel(ylab, fontsize=7)
        plt.xticks(fontsize=7)
        plt.yticks(fontsize=7)
        plt.title(f'Transformed Query and Reference Spectra', fontsize=10)


    plt.subplots_adjust(top=0.8, hspace=0.92, bottom=0.3)
    plt.figlegend(loc = 'upper center')
    fig.text(0.05, 0.15, f'Similarity Measure: {similarity_measure.capitalize()}', fontsize=7)
    fig.text(0.05, 0.12, f'Similarity Score: {round(similarity_score,4)}', fontsize=7)
    fig.text(0.05, 0.09, f"Spectrum Preprocessing Order: {''.join(spectrum_preprocessing_order)}", fontsize=7)
    fig.text(0.05, 0.06, f'High Quality Reference Library: {str(high_quality_reference_library)}', fontsize=7)
    fig.text(0.05, 0.03, f'Raw-Scale M/Z Range: [{min_mz},{max_mz}]', fontsize=7)
    fig.text(0.45, 0.15, f'Raw-Scale Intensity Range: [{int_min_tmp},{int_max_tmp}]', fontsize=7)
    fig.text(0.45, 0.12, f'Noise Threshold: {noise_threshold}', fontsize=7)
    fig.text(0.45, 0.09, f'Weight Factors (m/z,intensity): ({wf_mz},{wf_intensity})', fontsize=7)
    fig.text(0.45, 0.06, f'Low-Entropy Threshold: {LET_threshold}', fontsize=7)
    if similarity_measure=='mixture':
        fig.text(0.45, 0.03, f'Weights for mixture similarity: {weights}', fontsize=7)
    plt.savefig(output_path, format='pdf')

    if return_plot == True:
        return fig

