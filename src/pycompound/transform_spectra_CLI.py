
from pycompound.processing import transform_spectra
import pandas as pd
import argparse
import json
from pathlib import Path
import sys


def main() -> None:
    parser = argparse.ArgumentParser()

    parser.add_argument('--spectra_data', type=str, metavar='\b', help='file containing spectra mass spectrum/spectra to be transformed. Each row should correspond to a mass spectrum, the left-most column should contain an identifier, and each of the other columns should correspond to a single mass/charge ratio. Mandatory argument.')
    parser.add_argument('--spectrum_preprocessing_order', type=str, default='FNWL', metavar='\b', help='The spectrum preprocessing transformations and the order in which they are to be applied. Format must be a string with 2-4 characters chosen from F, N, L, W representing filtering based on mass/charge and intensity values, noise removal, low-entropy transformation, and weight-factor-transformation, respectively. For example, if \'WL\' is passed, then each spectrum will undergo a weight factor transformation followed by a low-entropy transformation. Default: FNWL')
    parser.add_argument('--mz_min', type=int, default=0, metavar='\b', help='Remove all peaks with mass/charge less than mz_min in each spectrum. Default: 0')
    parser.add_argument('--mz_max', type=int, default=999999999999, metavar='\b', help='Remove all peaks with mass/charge greater than mz_max in each spectrum. Default: 999999999999')
    parser.add_argument('--int_min', type=float, default=0, metavar='\b', help='Remove all peaks with intensity less than int_min in each spectrum. Default: 0')
    parser.add_argument('--int_max', type=float, default=999999999999, metavar='\b', help='Remove all peaks with intensity greater than int_max in each spectrum. Default: 999999999999')
    parser.add_argument('--window_size_centroiding', type=float, default=0.5, metavar='\b', help='Window size parameter used in centroiding a given spectrum. Only for HRMS. Default: 0.5')
    parser.add_argument('--window_size_matching', type=float, default=0.5, metavar='\b', help='Window size parameter used in matching a query spectrum and a reference library spectrum. Only for HRMS. Default: 0.5')
    parser.add_argument('--noise_threshold', type=float, default=0, metavar='\b', help='Ion fragments (i.e., points in a given mass spectrum) with intensity less than max(intensities)*noise_threshold are removed. Default: 0')
    parser.add_argument('--wf_mz', type=float, default=0, metavar='\b', help='Mass/charge weight factor parameter. Default: 0.')
    parser.add_argument('--wf_intensity', type=float, default=1, metavar='\b', help='Intensity weight factor parameter. Default: 1.')
    parser.add_argument('--LET_threshold', type=float, default=0, metavar='\b', help='Low-entropy transformation threshold parameter. Spectra with Shannon entropy less than LET_threshold are transformed according to intensitiesNew=intensitiesOriginal^{(1+S)/(1+LET_threshold)}. Default: 0.')
    parser.add_argument('--output_path', type=str, metavar='\b', help='Output (txt, mgf, msp) file containing the processed spectra. If no argument is passed, then the plots will be saved to the text file ./processed_spectra.txt in the current working directory.')

    args = parser.parse_args()
    if output_path is None:
        output_path = f'{Path.cwd()}/processed_spectra.txt'
        print(f'Warning: writing processed spectral data to {output_path}')
    output_extension = output_path.rsplit('.',1)
    output_extension = output_extension[(len(output_extension)-1)]
    if output_extension not in ['mgf', 'MGF', 'msp', 'MSP', 'txt', 'TXT']:
        print('Error: output_path must specify a txt, mgf, or msp file')
        sys.exit()


    transform_spectra(spectra_data=args.spectra_data, spectrum_preprocessing_order=spectrum_preprocessing_order, high_quality_reference_library=args.high_quality_reference_library, mz_min=args.mz_min, mz_max=args.mz_max, int_min=args.int_min, int_max=args.int_max, noise_threshold=args.noise_threshold, wf_mz=args.wf_mz, wf_intensity=args.wf_intensity, LET_threshold=args.LET_threshold, entropy_dimension=args.entropy_dimension, y_axis_transformation=args.y_axis_transformation, output_path=args.output_path)



if __name__ == "__main__":
    main()

