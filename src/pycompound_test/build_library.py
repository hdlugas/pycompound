
# this script has a function to extract the mass spectra from an mgf, mzML, or cdf file and write them in the necessary format for use in spectral library matching

import netCDF4 as nc
import pandas as pd
from pathlib import Path
from pyteomics import mgf
from pyteomics import mzml
import sys

def build_library_from_raw_data(input_path=None, output_path=None, is_reference=False):
    '''
    Converts mgf, mzML, or cdf file to the necessary format for spectral library matching.

    --input_path: Path to input file (must be mgf, mzML, or cdf file). Mandatory argument.
    --output_path: Path to output CSV file. Default: current working directory.
    --is_reference: Boolean flag indicating whether IDs of spectra should be written to output. Only pass true if building a reference library with known compound IDs. Only applicable to mgf files. Options: \'True\', \'False\'. Optional argument. Default: False.
    '''

    if input_path is None:
        print('Error: please specify input_path (i.e. the path to the input mgf, mzML, or cdf file). Mandatory argument.')
        sys.exit()

    if output_path is None:
        #print('Warning: no output_path specified, so library is written to {Path.cwd()}/build_library.csv')
        tmp = input_path.split('/')
        tmp = tmp[(len(tmp)-1)]
        basename = tmp.split('.')[0]
        output_path = f'{Path.cwd()}/{basename}.csv'
        print(f'Warning: no output_path specified, so library is written to {output_path}')

    if is_reference not in [True,False]:
        print('Error: is_reference must be either \'True\' or \'False\'.')
        sys.exit()

    # determine whether an mgf or a mzML file was passed to --input_path
    last_three_chars = input_path[(len(input_path)-3):len(input_path)]
    last_four_chars = input_path[(len(input_path)-4):len(input_path)]
    if last_three_chars == 'mgf' or last_three_chars == 'MGF':
        input_file_type = 'mgf'
    elif last_four_chars == 'mzML' or last_four_chars == 'mzml' or last_four_chars == 'MZML':
        input_file_type = 'mzML'
    elif last_three_chars == 'cdf' or last_three_chars == 'CDF':
        input_file_type = 'cdf'
    else:
        print('ERROR: either an \'mgf\', \'mzML\', or \'cdf\' file must be passed to --input_path')
        sys.exit()


    # obtain a list of spectra from the input file
    spectra = []

    if input_file_type == 'mgf':
        with mgf.read(input_path, index_by_scans = True) as reader:
            for spec in reader:
                spectra.append(spec)

    if input_file_type == 'mzML':
        with mzml.read(input_path) as reader:
            for spec in reader:
                spectra.append(spec)


    # extract the relevant information from each spectra (i.e m/z ratios and intensities)
    if input_file_type == 'mgf' or input_file_type == 'mzML':
        ids = []
        mzs = []
        ints = []
        for i in range(0,len(spectra)):
            for j in range(0,len(spectra[i]['m/z array'])):
                if input_file_type == 'mzML':
                    ids.append(f'ID_{i+1}')
                else:
                    if is_reference == False:
                        ids.append(f'ID_{i+1}')
                    elif is_reference == True:
                        ids.append(spectra[i]['params']['name'])
                mzs.append(spectra[i]['m/z array'][j])
                ints.append(spectra[i]['intensity array'][j])

    if input_file_type == 'cdf':
        dataset = nc.Dataset(input_path, 'r')
        all_mzs = dataset.variables['mass_values'][:]
        all_ints = dataset.variables['intensity_values'][:]
        scan_idxs = dataset.variables['scan_index'][:]
        dataset.close()

        ids = []
        mzs = []
        ints = []
        for i in range(0,(len(scan_idxs)-1)):
            if i % 1000 == 0:
                print(f'analyzed {i} out of {len(scan_idxs)} scans')
            s_idx = scan_idxs[i]
            e_idx = scan_idxs[i+1]

            mzs_tmp = all_mzs[s_idx:e_idx]
            ints_tmp = all_ints[s_idx:e_idx]

            for j in range(0,len(mzs_tmp)):
                ids.append(f'ID_{i+1}')
                mzs.append(mzs_tmp[j])
                ints.append(ints_tmp[j])


    # write CSV file of spectra for use in spectral library matching
    df = pd.DataFrame({'id':ids, 'mz_ratio':mzs, 'intensity':ints})
    df.to_csv(output_path, index=False)


