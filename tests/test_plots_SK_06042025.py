
from plot_spectra import generate_plots_on_HRMS_data
from plot_spectra import generate_plots_on_NRMS_data
from pathlib import Path
import os

print('\n\ntest #0:')
generate_plots_on_HRMS_data(
        query_data=f'{Path.cwd()}/../data/lcms_query_library.csv',
        reference_data=f'{Path.cwd()}/../data/lcms_reference_library.csv',
        output_path=f'{Path.cwd()}/../tests/plots/test_no_wf_normalized_y_axis_no_mz_zoom.pdf')

print('\n\ntest #1:')
generate_plots_on_HRMS_data(
        query_data=f'{Path.cwd()}/../data/lcms_query_library.csv',
        reference_data=f'{Path.cwd()}/../data/lcms_reference_library.csv',
        wf_mz=2,
        wf_intensity=0.5,
        output_path=f'{Path.cwd()}/../tests/plots/test_wf_normalized_y_axis_no_mz_zoom.pdf')

print('\n\ntest #2:')
generate_plots_on_HRMS_data(
        query_data=f'{Path.cwd()}/../data/lcms_query_library.csv',
        reference_data=f'{Path.cwd()}/../data/lcms_reference_library.csv',
        y_axis_transformation='log10',
        output_path=f'{Path.cwd()}/../tests/plots/test_no_wf_log10_y_axis_no_mz_zoom.pdf')

print('\n\ntest #3:')
generate_plots_on_HRMS_data(
        query_data=f'{Path.cwd()}/../data/lcms_query_library.csv',
        reference_data=f'{Path.cwd()}/../data/lcms_reference_library.csv',
        y_axis_transformation='sqrt',
        output_path=f'{Path.cwd()}/../tests/plots/test_no_wf_sqrt_y_axis_no_mz_zoom.pdf')

print('\n\ntest #4:')
generate_plots_on_HRMS_data(
        query_data=f'{Path.cwd()}/../data/lcms_query_library.csv',
        reference_data=f'{Path.cwd()}/../data/lcms_reference_library.csv',
        mz_min = 400,
        mz_max = 650,
        y_axis_transformation='sqrt',
        output_path=f'{Path.cwd()}/../tests/plots/test_no_wf_normalized_y_axis_mz_zoom.pdf')


