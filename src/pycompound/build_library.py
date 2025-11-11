
import netCDF4 as nc
import numpy as np
import pandas as pd
from pathlib import Path
from pyteomics import mgf
from pyteomics import mzml
import sys
import json

def build_library_from_raw_data(input_path=None, output_path=None, is_reference=False):
    '''
    Converts mgf, mzML, cdf, json, or msp file to the necessary format for spectral library matching.

    --input_path: Path to input file (must be mgf, mzML, cdf, json, or msp file). Mandatory argument.
    --output_path: Path to output TXT file. Default: current working directory.
    --is_reference: Boolean flag indicating whether IDs of spectra should be written to output. Only pass true if building a reference library with known compound IDs. Only applicable to mgf and msp files. Options: \'True\', \'False\'. Optional argument. Default: False.
    '''

    if input_path is None:
        print('Error: please specify input_path (i.e. the path to the input mgf, mzML, cdf, json, or msp file). Mandatory argument.')
        sys.exit()

    if output_path is None:
        tmp = input_path.split('/')
        tmp = tmp[(len(tmp)-1)]
        basename = tmp.split('.')[0]
        output_path = f'{Path.cwd()}/{basename}.csv'
        print(f'Warning: no output_path specified, so library is written to {output_path}')

    if is_reference not in [True,False]:
        print('Error: is_reference must be either \'True\' or \'False\'.')
        sys.exit()

    last_three_chars = input_path[(len(input_path)-3):len(input_path)]
    last_four_chars = input_path[(len(input_path)-4):len(input_path)]
    if last_three_chars == 'mgf' or last_three_chars == 'MGF':
        input_file_type = 'mgf'
    elif last_four_chars == 'mzML' or last_four_chars == 'mzml' or last_four_chars == 'MZML':
        input_file_type = 'mzML'
    elif last_four_chars == 'json' or last_four_chars == 'JSON':
        input_file_type = 'json'
    elif last_three_chars == 'cdf' or last_three_chars == 'CDF':
        input_file_type = 'cdf'
    elif last_three_chars == 'msp' or last_three_chars == 'MSP':
        input_file_type = 'msp'
    else:
        print('ERROR: either an \'mgf\', \'mzML\', \'cdf\', \'json\', or \'msp\' file must be passed to --input_path')
        sys.exit()


    spectra = []
    if input_file_type == 'mgf':
        #with mgf.read(input_path, index_by_scans = True) as reader:
        with mgf.read(input_path, use_index=False) as reader:
            for spec in reader:
                spectra.append(spec)
    if input_file_type == 'mzML':
        with mzml.read(input_path) as reader:
            for spec in reader:
                spectra.append(spec)


    if input_file_type == 'mgf' or input_file_type == 'mzML':
        ids = []
        mzs = []
        ints = []
        precursor_ion_mzs = []
        for i in range(0,len(spectra)):
            for j in range(0,len(spectra[i]['m/z array'])):
                if input_file_type == 'mzML':
                    if is_reference == False:
                        ids.append(f'ID_{i+1}')
                    else:
                        ids.append(spectra[i]['id'])
                elif input_file_type == 'mgf':
                    precursor_ion_mzs.append(spectra[i]['params']['pepmass'][0])
                    if is_reference == False:
                        ids.append(f'ID_{i+1}')
                    else:
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



    if input_file_type == "msp":
        ids = []
        mzs = []
        ints = []
        precursor_ion_mzs = []
        spectrum_id = None
        precursor_ion_mz = None
        with open(input_path, "r", encoding="utf-8", errors="ignore") as f:
            i = 0
            for line in f:
                line = line.strip()
                if not line:
                    continue

                if line.startswith("Name:"):
                    i += 1
                    if not is_reference:
                        spectrum_id = f"ID_{i}"
                    else:
                        spectrum_id = line.replace("Name:", "", 1).strip()

                elif line.startswith("PrecursorMZ:"):
                    try:
                        precursor_ion_mz = float(line.replace("PrecursorMZ:", "", 1).strip())
                    except ValueError:
                        precursor_ion_mz = None

                elif line[0].isdigit():
                    try:
                        mz, intensity = map(float, line.split()[:2])
                    except ValueError:
                        continue

                    if spectrum_id is None:
                        continue

                    ids.append(spectrum_id)
                    mzs.append(mz)
                    ints.append(intensity)
                    precursor_ion_mzs.append(precursor_ion_mz)



    if input_file_type == 'json':
        data = json.load(open(input_path))
        ids = []
        mzs = []
        ints = []
        precursor_ion_mzs = []
        for i in range(0,len(data)):
            spec_ID_tmp = data[i]['spectrum_id']
            tmp = data[i]['peaks_json']
            tmp = tmp[1:-1].split(",")
            tmp = [a.replace("[","") for a in tmp]
            tmp = [a.replace("]","") for a in tmp]
            mzs_tmp = tmp[0::2]
            ints_tmp = tmp[1::2]
            if is_reference == False:
                ids.extend([f'ID_{i+1}'] * len(mzs_tmp))
            elif is_reference == True:
                ids.extend([spec_ID_tmp] * len(mzs_tmp))
            mzs.extend(mzs_tmp)
            ints.extend(ints_tmp)
            precursor_ion_mzs.extend([data[i]['Precursor_MZ']] * len(mzs_tmp))


    if len(precursor_ion_mzs) > 0:
        df = pd.DataFrame({'id':ids, 'mz_ratio':mzs, 'intensity':ints, 'precursor_ion_mz':precursor_ion_mzs})
    else:
        df = pd.DataFrame({'id':ids, 'mz_ratio':mzs, 'intensity':ints})

    df.to_csv(output_path, index=False, sep='\t')


