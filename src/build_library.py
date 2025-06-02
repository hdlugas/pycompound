
# this script extracts the mass spectra from an mgf, mzML, or cdf file and saves them in the necessary format for use in spectral library matching

# import libraries
import argparse
import netCDF4 as nc
import pandas as pd
from pathlib import Path
from pyteomics import mgf
from pyteomics import mzml
import sys

# create argparse object so command-line input can be imported
parser = argparse.ArgumentParser()

# import command-line arguments
parser.add_argument('--input_path', metavar='\b', help='Path to input file (must be mgf, mlMZ, or cdf file). Mandatory argument.')
parser.add_argument('--output_path', metavar='\b', help='Path to output CSV file. Default: current working directory.')
parser.add_argument('--is_reference', metavar='\b', help='Boolean flag indicating whether IDs of spectra should be written to output. Only pass True if building a reference library with known compound IDs. Only applicable to MGF files. Options: \'True\', \'False\'. Optional argument. Default: False.')

# parse the user-input arguments
args = parser.parse_args()

# import the path to the MGF file
if args.input_path is None:
    print(f'ERROR: the path to the input file must be passed to --input_path.')
    sys.exit()
else:
    input_path = args.input_path

# import the path to the output CSV file
if args.output_path is None:
    print(f'WARNING: no argument was passed to --output_path. The resulting CSV file will be written to the current working directory.')
    output_path = f'{Path.cwd()}/built_library.csv'
else:
    output_path = args.output_path

# import the flag indicating whether 
if args.is_reference is None:
    is_reference = False
else:
    is_reference = args.is_reference

if is_reference != 'True' and is_reference is not 'False':
    print(f'ERROR: is_reference must be either \'True\' or \'False\'')
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
    print(f'ERROR: either an \'mgf\', \'mzML\', or \'cdf\' file must be passed to --input_path')


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
                if is_reference == 'False':
                    ids.append(f'ID_{i+1}')
                elif is_reference == 'True':
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



