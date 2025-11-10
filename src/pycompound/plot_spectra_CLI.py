
from pycompound.plot_spectra import generate_plots_on_HRMS_data
from pycompound.plot_spectra import generate_plots_on_NRMS_data
import pandas as pd
import argparse
import json
from pathlib import Path
import sys


parser = argparse.ArgumentParser()

parser.add_argument('--query_data', type=str, metavar='\b', help='CSV file of query mass spectrum/spectra to be identified. Each row should correspond to a mass spectrum, the left-most column should contain an identifier, and each of the other columns should correspond to a single mass/charge ratio. Mandatory argument.')
parser.add_argument('--reference_data', type=str, metavar='\b', help='CSV file of the reference mass spectra. Each row should correspond to a mass spectrum, the left-most column should contain in identifier (i.e. the CAS registry number or the compound name), and the remaining column should correspond to a single mass/charge ratio. Mandatory argument.')
parser.add_argument('--spectrum_ID1', type=str, metavar='\b', help='The identifier of the query spectrum to be plotted. Default: first query spectrum in query_data.')
parser.add_argument('--spectrum_ID2', type=str, metavar='\b', help='The identifier of the reference spectrum to be plotted. Default: first reference spectrum in reference_data.')
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
parser.add_argument('--y_axis_transformation', type=str, default='normalized', metavar='\b', help='Transformation to apply to y-axis (i.e. intensity axis) of plots. Options: \'normalized\', \'none\', \'log10\', and \'sqrt\'. Default: normalized.')
parser.add_argument('--output_path', type=str, metavar='\b', help='Output PDF file containing the plots of the query and reference spectra before and after preprocessing transformations. If no argument is passed, then the plots will be saved to the PDF ./query_spec_{query_spectrum_ID}_reference_spec_{reference_spectrum_ID}_plot.pdf in the current working directory.')

args = parser.parse_args()

if args.chromatography_platform == 'HRMS':
    spectrum_preprocessing_order = 'FCNMWL'
elif args.chromatography_platform == 'NRMS':
    spectrum_preprocessing_order = 'FNLW'
else:
    print('Error: chromatography_platform must be either \'HRMS\' or \'NRMS\'')
    sys.exit()


if args.chromatography_platform == 'HRMS':
    generate_plots_on_HRMS_data(query_data=args.query_data, reference_data=args.reference_data, spectrum_ID1=args.spectrum_ID1, spectrum_ID2=args.spectrum_ID2, similarity_measure=args.similarity_measure, weights=args.weights, spectrum_preprocessing_order=spectrum_preprocessing_order, high_quality_reference_library=args.high_quality_reference_library, mz_min=args.mz_min, mz_max=args.mz_max, int_min=args.int_min, int_max=args.int_max, window_size_centroiding=args.window_size_centroiding, window_size_matching=args.window_size_matching, noise_threshold=args.noise_threshold, wf_mz=args.wf_mz, wf_intensity=args.wf_intensity, LET_threshold=args.LET_threshold, entropy_dimension=args.entropy_dimension, y_axis_transformation=args.y_axis_transformation, output_path=args.output_path)
elif args.chromatography_platform == 'NRMS':
    generate_plots_on_NRMS_data(query_data=args.query_data, reference_data=args.reference_data, spectrum_ID1=args.spectrum_ID1, spectrum_ID2=args.spectrum_ID2, similarity_measure=args.similarity_measure, weights=args.weights, spectrum_preprocessing_order=spectrum_preprocessing_order, high_quality_reference_library=args.high_quality_reference_library, mz_min=args.mz_min, mz_max=args.mz_max, int_min=args.int_min, int_max=args.int_max, noise_threshold=args.noise_threshold, wf_mz=args.wf_mz, wf_intensity=args.wf_intensity, LET_threshold=args.LET_threshold, entropy_dimension=args.entropy_dimension, y_axis_transformation=args.y_axis_transformation, output_path=args.output_path)


