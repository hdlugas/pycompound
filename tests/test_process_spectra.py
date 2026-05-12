
from pycompound.processing import transform_spectra
from pathlib import Path
import os

print('\n\ntest #1')
transform_spectra(spectra_data = f'{Path.cwd()}/data/lcms_query.txt',
                       output_path = f'{Path.cwd()}/data/example_processed_data1.mgf',
                       spectrum_preprocessing_order = 'WL',
                       wf_mz = 0.8,
                       wf_intensity = 1.1,
                       LET_threshold = 3)


print('\n\ntest #2')
transform_spectra(spectra_data = f'{Path.cwd()}/data/lcms_query.txt',
                       output_path = f'{Path.cwd()}/data/example_processed_data2.msp',
                       spectrum_preprocessing_order = 'WL',
                       wf_mz = 0.8,
                       wf_intensity = 1.1,
                       LET_threshold = 3)

print('\n\ntest #3')
transform_spectra(spectra_data = f'{Path.cwd()}/data/lcms_query.txt',
                       output_path = f'{Path.cwd()}/data/example_processed_data3.txt',
                       spectrum_preprocessing_order = 'WL',
                       wf_mz = 0.8,
                       wf_intensity = 1.1,
                       LET_threshold = 3)

print('\n\ntest #4')
transform_spectra(spectra_data = f'{Path.cwd()}/data/GNPS-SELLECKCHEM-FDA-PART1.mgf',
                       output_path = f'{Path.cwd()}/data/example_processed_data4.txt',
                       spectrum_preprocessing_order = 'FLNW',
                       wf_mz = 0.8,
                       wf_intensity = 1.1,
                       LET_threshold = 3,
                       mz_min = 50,
                       noise_threshold=0.1)

print('\n\ntest #5')
transform_spectra(spectra_data = f'{Path.cwd()}/data/MoNA-export-Human_Plasma_Quant.msp',
                       output_path = f'{Path.cwd()}/data/example_processed_data5.txt',
                       spectrum_preprocessing_order = 'FLNW',
                       wf_mz = 0.8,
                       wf_intensity = 1.1,
                       LET_threshold = 3,
                       mz_min = 50,
                       noise_threshold=0.1)

"""
print('\n\ntest #1')
transform_spectra_HRMS(spectra_data = f'{Path.cwd()}/data/lcms_query.txt',
                       output_path = f'{Path.cwd()}/data/example_processed_lcms_data.sam')
"""
