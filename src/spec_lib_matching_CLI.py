
from pycompound.spec_lib_matching import run_spec_lib_matching_on_HRMS_data 
from pycompound.spec_lib_matching import run_spec_lib_matching_on_NRMS_data 
from pathlib import Path
import pandas as pd
import argparse
import sys
import json


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
parser.add_argument('--print_id_results', type=str, default=False, metavar='\b', help='Flag that prints identification results if True. Default: False')
parser.add_argument('--output_identification', type=str, default=f'{Path.cwd()}/output_identification.txt', metavar='\b', help='Output TXT file containing the most-similar reference spectra for each query spectrum along with the corresponding similarity scores. Default is to save identification output in current working directory (i.e. same directory this script is contained in) with filename \'output_identification.txt\'.')
parser.add_argument('--output_similarity_scores', type=str, default=f'{Path.cwd()}/output_all_similarity_scores.txt', metavar='\b', help='Output TXT file containing similarity scores between all query spectrum/spectra and all reference spectra. Each row corresponds to a query spectrum, the left-most column contains the query spectrum/spectra identifier, and the remaining column contain the similarity scores with respect to all reference library spectra. If no argument passed, then this TXT file is written to the current working directory with filename \'output_all_similarity_scores\'.txt.')

args = parser.parse_args()


if args.chromatography_platform == 'HRMS':
    run_spec_lib_matching_on_HRMS_data(query_data=args.query_data,
                                       reference_data=args.reference_data,
                                       precursor_ion_mz_tolerance=args.precursor_ion_mz_tolerance,
                                       ionization_mode=args.ionization_mode,
                                       adduct=args.adduct,
                                       likely_reference_ids=args.likely_reference_ids,
                                       similarity_measure=args.similarity_measure,
                                       weights=args.weights,
                                       spectrum_preprocessing_order=args.spectrum_preprocessing_order,
                                       high_quality_reference_library=args.high_quality_reference_library,
                                       mz_min=args.mz_min,
                                       mz_max=args.mz_max,
                                       int_min=args.int_min,
                                       int_max=args.int_max,
                                       window_size_centroiding=args.window_size_centroiding,
                                       window_size_matching=args.window_size_matching,
                                       noise_threshold=args.noise_threshold,
                                       wf_mz=args.wf_mz,
                                       wf_intensity=args.wf_intensity,
                                       LET_threshold=args.LET_threshold,
                                       entropy_dimension=args.entropy_dimension,
                                       n_top_matches_to_save=args.n_top_matches_to_save,
                                       print_id_results=args.print_id_results,
                                       output_identification=args.output_identification,
                                       output_similarity_scores=args.output_similarity_scores)


if args.chromatography_platform == 'NRMS':
    run_spec_lib_matching_on_NRMS_data(query_data=args.query_data,
                                       reference_data=args.reference_data,
                                       likely_reference_ids=args.likely_reference_ids,
                                       similarity_measure=args.similarity_measure,
                                       weights=args.weights,
                                       spectrum_preprocessing_order=args.spectrum_preprocessing_order,
                                       high_quality_reference_library=args.high_quality_reference_library,
                                       mz_min=args.mz_min,
                                       mz_max=args.mz_max,
                                       int_min=args.int_min,
                                       int_max=args.int_max,
                                       noise_threshold=args.noise_threshold,
                                       wf_mz=args.wf_mz,
                                       wf_intensity=args.wf_intensity,
                                       LET_threshold=args.LET_threshold,
                                       entropy_dimension=args.entropy_dimension,
                                       n_top_matches_to_save=args.n_top_matches_to_save,
                                       print_id_results=args.print_id_results,
                                       output_identification=args.output_identification,
                                       output_similarity_scores=args.output_similarity_scores)

