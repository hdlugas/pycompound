
from pycompound.spec_lib_matching import tune_params_on_HRMS_data_grid
from pycompound.spec_lib_matching import tune_params_on_NRMS_data_grid
import argparse
import json
from pathlib import Path
import sys

parser = argparse.ArgumentParser()

parser.add_argument('--query_data', type=str, metavar='\b', help='CSV file of query mass spectrum/spectra to be identified. Each row should correspond to a mass spectrum, the left-most column should contain an identifier, and each of the other columns should correspond to a single mass/charge ratio. Mandatory argument.')
parser.add_argument('--reference_data', type=str, metavar='\b', help='CSV file of the reference mass spectra. Each row should correspond to a mass spectrum, the left-most column should contain in identifier (i.e. the CAS registry number or the compound name), and the remaining column should correspond to a single mass/charge ratio. Mandatory argument.')
parser.add_argument('--precursor_ion_mz_tolerance', type=str, metavar='\b', default=None, help='Precursor ion m/z tolerance (positive real number; only applicable to HRMS)). Default=None')
parser.add_argument('--ionization_mode', type=str, metavar='\b', default=None, help='Ionization mode (only applicable to HRMS). Options: \'Positive\', \'Negative\', or \'N/A\'.')
parser.add_argument('--adduct', type=str, metavar='\b', default='H', help='Adduct (only applicable to HRMS). Options: \'H\', \'NH3\', \'NH4\', \'OH\', \'Cl\', \'K\', \'Li\', \'Na\'. Default: \'H\'.')
parser.add_argument('--likely_reference_ids', type=str, metavar='\b', help='CSV file with one column containing the IDs of a subset of all compounds in the reference_data to be used in spectral library matching. Each ID in this file must be an ID in the reference library. Default: none (i.e. default is to use entire reference library)')
parser.add_argument('--similarity_measure', type=str, default='cosine', metavar='\b', help='Similarity measure: options are cosine, shannon, renyi, tsallis, mixture, jaccard, dice, 3w_jaccard, sokal_sneath, binary_cosine, mountford, mcconnaughey, driver_kroeber, simpson, braun_banquet, fager_mcgowan, kulczynski, intersection, hamming, or hellinger. Default: cosine.')
parser.add_argument('--weights', type=json.loads, default={'Cosine':0.25,'Shannon':0.25,'Renyi':0.25,'Tsallis':0.25}, metavar='\b', help='dict of weights to give to each non-binary similarity measure (i.e. cosine, shannon, renyi, and tsallis) when the mixture similarity measure is specified. Default: 0.25 for each of the four non-binary similarity measures.')
parser.add_argument('--chromatography_platform', type=str, metavar='\b', help='Chromatography platform: options are \'HRMS\' and \'NRMS\'. Mandatory argument.')
parser.add_argument('--spectrum_preprocessing_order', type=str, metavar='\b', default=None, help='The LC-MS/MS spectrum preprocessing transformations and the order in which they are to be applied. Note that these transformations are applied prior to computing similarity scores. Format must be a string with 2-6 characters chosen from C, F, M, N, L, W representing centroiding, filtering based on mass/charge and intensity values, matching, noise removal, low-entropy trannsformation, and weight-factor-transformation, respectively. For example, if \'WCM\' is passed, then each spectrum will undergo a weight factor transformation, then centroiding, and then matching. Note that if an argument is passed, then \'M\' must be contained in the argument, since matching is a required preprocessing step in spectral library matching of LC-MS/MS data. Furthermore, \'C\' must be performed before matching since centroiding can change the number of ion fragments in a given spectrum. Default: FCNMWL for HRMS, FNLW for NRMS')
parser.add_argument('--high_quality_reference_library', type=str, default='False', metavar='\b', help='True/False flag indicating whether the reference library is considered to be of high quality. If True, then the spectrum preprocessing transformations of filtering and noise removal are performed only on the query spectrum/spectra. If False, all spectrum preprocessing transformations specified will be applied to both the query and reference spectra. Default: False')
parser.add_argument('--mz_min', type=str, default='0', metavar='\b', help='Remove all peaks with mass/charge less than mz_min in each spectrum. Default: 0')
parser.add_argument('--mz_max', type=str, default='999999999999', metavar='\b', help='Remove all peaks with mass/charge greater than mz_max in each spectrum. Default: 999999999999')
parser.add_argument('--int_min', type=str, default='0', metavar='\b', help='Remove all peaks with intensity less than int_min in each spectrum. Default: 0')
parser.add_argument('--int_max', type=str, default='999999999999', metavar='\b', help='Remove all peaks with intensity greater than int_max in each spectrum. Default: 999999999999')
parser.add_argument('--window_size_centroiding', type=str, default='0.5', metavar='\b', help='Window size parameter used in centroiding a given spectrum. Only for HRMS. Default: 0.5')
parser.add_argument('--window_size_matching', type=str, default='0.5', metavar='\b', help='Window size parameter used in matching a query spectrum and a reference library spectrum. Only for HRMS. Default: 0.5')
parser.add_argument('--noise_threshold', type=str, default='0', metavar='\b', help='Ion fragments (i.e. points in a given mass spectrum) with intensity less than max(intensities)*noise_threshold are removed. Default: 0')
parser.add_argument('--wf_mz', type=str, default='0', metavar='\b', help='Mass/charge weight factor parameter. Default: 0.')
parser.add_argument('--wf_intensity', type=str, default='1', metavar='\b', help='Intensity weight factor parameter. Default: 1.')
parser.add_argument('--LET_threshold', type=str, default='0', metavar='\b', help='Low-entropy transformation threshold parameter. Spectra with Shannon entropy less than LET_threshold are transformed according to intensitiesNew=intensitiesOriginal^{(1+S)/(1+LET_threshold)}. Default: 0.')
parser.add_argument('--entropy_dimension', type=str, default='1.1', metavar='\b', help='Entropy dimension parameter. Must have positive value other than 1. When the entropy dimension is 1, then Renyi and Tsallis entropy are equivalent to Shannon entropy. Therefore, this parameter only applies to the renyi and tsallis similarity measures. This parameter will be ignored if similarity measure cosine or shannon is chosen. Default: 1.1.')
parser.add_argument('--output_path', type=str, default=f'{Path.cwd()}/output_tuning.txt', metavar='\b', help='Output TXT file containing one row for each parameter set used along with its corresponding accuracy. If no argument passed, then this TXT file is written to the current working directory with filename \'output_all_similarity_scores\'.txt.')

args = parser.parse_args()

if args.chromatography_platform == 'HRMS' and args.spectrum_preprocessing_order == None:
    spectrum_preprocessing_order = 'FCNMWL'
elif args.chromatography_platform == 'NRMS' and args.spectrum_preprocessing_order == None:
    spectrum_preprocessing_order = 'FNLW'
else:
    print('Error: chromatography_platform must be either \'HRMS\' or \'NRMS\'')
    sys.exit()


grid = {'similarity_measure':args.similarity_measure.split(','), 'weight':[args.weights], 'spectrum_preprocessing_order':spectrum_preprocessing_order.split(','), 'mz_min':args.mz_min.split(','), 'mz_max':args.mz_max.split(','), 'int_min':args.int_min.split(','), 'int_max':args.int_max.split(','), 'window_size_centroiding':args.window_size_centroiding.split(','), 'window_size_matching':args.window_size_matching.split(','), 'noise_threshold':args.noise_threshold.split(','), 'wf_mz':args.wf_mz.split(','), 'wf_int':args.wf_intensity.split(','), 'LET_threshold':args.LET_threshold.split(','), 'entropy_dimension':args.entropy_dimension.split(','), 'high_quality_reference_library':args.high_quality_reference_library.split(',')}

if args.chromatography_platform == 'HRMS':
    grid['mz_min'] = [float(x) for x in grid['mz_min']]
    grid['mz_max'] = [float(x) for x in grid['mz_max']]
elif args.chromatography_platform == 'NRMS':
    grid['mz_min'] = [int(x) for x in grid['mz_min']]
    grid['mz_max'] = [int(x) for x in grid['mz_max']]

grid['int_min'] = [float(x) for x in grid['int_min']]
grid['int_max'] = [float(x) for x in grid['int_max']]
grid['window_size_centroiding'] = [float(x) for x in grid['window_size_centroiding']]
grid['window_size_matching'] = [float(x) for x in grid['window_size_matching']]
grid['noise_threshold'] = [float(x) for x in grid['noise_threshold']]
grid['wf_mz'] = [float(x) for x in grid['wf_mz']]
grid['wf_int'] = [float(x) for x in grid['wf_int']]
grid['LET_threshold'] = [float(x) for x in grid['LET_threshold']]
grid['entropy_dimension'] = [float(x) for x in grid['entropy_dimension']]

if args.precursor_ion_mz_tolerance == None:
    precursor_ion_mz_tolerance_tmp = None
else:
    precursor_ion_mz_tolerance_tmp = float(args.precursor_ion_mz_tolerance)


if args.chromatography_platform == 'HRMS':
    tune_params_on_HRMS_data_grid(query_data=args.query_data,
                                  reference_data=args.reference_data,
                                  precursor_ion_mz_tolerance=precursor_ion_mz_tolerance_tmp,
                                  ionization_mode=args.ionization_mode,
                                  adduct=args.adduct,
                                  grid=grid,
                                  output_path=args.output_path)

if args.chromatography_platform == 'NRMS':
    tune_params_on_NRMS_data_grid(query_data=args.query_data,
    reference_data=args.reference_data,
    grid=grid,
    output_path=args.output_path)

