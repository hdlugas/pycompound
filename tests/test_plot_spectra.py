
from pycompound.plot_spectra import generate_plots_on_HRMS_data
from pycompound.plot_spectra import generate_plots_on_NRMS_data
from pathlib import Path
import os

os.makedirs(f'{Path.cwd()}/plots', exist_ok=True)

print('\n\ntest #1:')
generate_plots_on_HRMS_data(
        query_data=f'{Path.cwd()}/data/lcms_query_library.csv',
        reference_data=f'{Path.cwd()}/data/lcms_reference_library.csv',
        high_quality_reference_library=True,
        noise_threshold=0.1,
        mz_min=100,
        output_path=f'{Path.cwd()}/plots/test1.pdf')

print('\n\ntest #2:')
generate_plots_on_HRMS_data(
        query_data=f'{Path.cwd()}/data/lcms_query_library.csv',
        reference_data=f'{Path.cwd()}/data/lcms_reference_library.csv',
        noise_threshold=0.1,
        similarity_measure='shannon',
        output_path=f'{Path.cwd()}/plots/test2.pdf')

print('\n\ntest #3:')
generate_plots_on_HRMS_data(
        query_data=f'{Path.cwd()}/data/lcms_query_library.csv',
        reference_data=f'{Path.cwd()}/data/lcms_reference_library.csv',
        similarity_measure='renyi',
        entropy_dimension=1.2,
        output_path=f'{Path.cwd()}/plots/test3.pdf')

print('\n\ntest #4:')
generate_plots_on_HRMS_data(
        query_data=f'{Path.cwd()}/data/lcms_query_library.csv',
        reference_data=f'{Path.cwd()}/data/lcms_reference_library.csv',
        similarity_measure='tsallis',
        entropy_dimension=1.2,
        output_path=f'{Path.cwd()}/plots/test4.pdf')

print('\n\ntest #5:')
generate_plots_on_HRMS_data(
        query_data=f'{Path.cwd()}/data/lcms_query_library.csv',
        reference_data=f'{Path.cwd()}/data/lcms_reference_library.csv',
        similarity_measure='tsallis',
        entropy_dimension=1.2,
        output_path=f'{Path.cwd()}/plots/test5.pdf')

print('\n\ntest #6:')
generate_plots_on_HRMS_data(
        query_data=f'{Path.cwd()}/data/lcms_query_library.csv',
        reference_data=f'{Path.cwd()}/data/lcms_reference_library.csv',
        wf_intensity=0.8,
        wf_mz=1.1,
        output_path=f'{Path.cwd()}/plots/test6.pdf')

print('\n\ntest #7:')
generate_plots_on_HRMS_data(
        query_data=f'{Path.cwd()}/data/lcms_query_library.csv',
        reference_data=f'{Path.cwd()}/data/lcms_reference_library.csv',
        window_size_centroiding=0.1,
        output_path=f'{Path.cwd()}/plots/test7.pdf')

print('\n\ntest #8:')
generate_plots_on_HRMS_data(
        query_data=f'{Path.cwd()}/data/lcms_query_library.csv',
        reference_data=f'{Path.cwd()}/data/lcms_reference_library.csv',
        window_size_matching=0.25,
        output_path=f'{Path.cwd()}/plots/test8.pdf')

print('\n\ntest #9:')
generate_plots_on_HRMS_data(
        query_data=f'{Path.cwd()}/data/lcms_query_library.csv',
        reference_data=f'{Path.cwd()}/data/lcms_reference_library.csv',
        spectrum_preprocessing_order='WCM',
        wf_mz=0.8,
        wf_intensity=1.1,
        output_path=f'{Path.cwd()}/plots/test9.pdf')

print('\n\ntest #10:')
generate_plots_on_HRMS_data(
        query_data=f'{Path.cwd()}/data/lcms_query_library.csv',
        reference_data=f'{Path.cwd()}/data/lcms_reference_library.csv',
        LET_threshold=3,
        output_path=f'{Path.cwd()}/plots/test10.pdf')

print('\n\ntest #11:')
generate_plots_on_HRMS_data(
        query_data=f'{Path.cwd()}/data/lcms_query_library.csv',
        reference_data=f'{Path.cwd()}/data/lcms_reference_library.csv',
        spectrum_ID1 = 212,
        spectrum_ID2 = 100,
        LET_threshold=3,
        output_path=f'{Path.cwd()}/plots/test11.pdf')

print('\n\ntest #12:')
generate_plots_on_HRMS_data(
        query_data=f'{Path.cwd()}/data/lcms_query_library.csv',
        reference_data=f'{Path.cwd()}/data/lcms_reference_library.csv',
        spectrum_ID1 = 'Jamaicamide A M+H',
        spectrum_ID2 = 'Malyngamide J M+H',
        LET_threshold=3,
        output_path=f'{Path.cwd()}/plots/test12.pdf')

print('\n\ntest #13:')
generate_plots_on_HRMS_data(
        query_data=f'{Path.cwd()}/data/lcms_query_library.csv',
        reference_data=f'{Path.cwd()}/data/lcms_reference_library.csv',
        spectrum_ID1 = 'Jamaicamide A M+H',
        spectrum_ID2 = 'Jamaicamide A M+H',
        LET_threshold=3,
        output_path=f'{Path.cwd()}/plots/test13.pdf')

print('\n\ntest #14:')
generate_plots_on_NRMS_data(
        query_data=f'{Path.cwd()}/data/gcms_query_library.csv',
        reference_data=f'{Path.cwd()}/data/gcms_reference_library.csv',
        output_path=f'{Path.cwd()}/plots/test14.pdf')

print('\n\ntest #15:')
generate_plots_on_NRMS_data(
        query_data=f'{Path.cwd()}/data/gcms_query_library.csv',
        reference_data=f'{Path.cwd()}/data/gcms_reference_library.csv',
        spectrum_ID1 = 463514,
        spectrum_ID2 = 112312,
        output_path=f'{Path.cwd()}/plots/test15.pdf')

print('\n\ntest #17:')
generate_plots_on_NRMS_data(
        query_data=f'{Path.cwd()}/data/gcms_query_library.csv',
        reference_data=f'{Path.cwd()}/data/gcms_reference_library.csv',
        output_path=f'{Path.cwd()}/plots/test17.pdf')

print('\n\ntest #18:')
generate_plots_on_NRMS_data(
        query_data=f'{Path.cwd()}/data/gcms_query_library.csv',
        reference_data=f'{Path.cwd()}/data/gcms_reference_library.csv',
        y_axis_transformation='none',
        output_path=f'{Path.cwd()}/plots/test18.pdf')

print('\n\ntest #19:')
generate_plots_on_NRMS_data(
        query_data=f'{Path.cwd()}/data/gcms_query_library.csv',
        reference_data=f'{Path.cwd()}/data/gcms_reference_library.csv',
        y_axis_transformation='log10',
        output_path=f'{Path.cwd()}/plots/test19.pdf')

print('\n\ntest #20:')
generate_plots_on_NRMS_data(
        query_data=f'{Path.cwd()}/data/gcms_query_library.csv',
        reference_data=f'{Path.cwd()}/data/gcms_reference_library.csv',
        y_axis_transformation='sqrt',
        output_path=f'{Path.cwd()}/plots/test20.pdf')

print('\n\ntest #21:')
generate_plots_on_HRMS_data(
        query_data=f'{Path.cwd()}/data/lcms_query_library.csv',
        reference_data=f'{Path.cwd()}/data/lcms_reference_library.csv',
        output_path=f'{Path.cwd()}/plots/test_no_wf_normalized_y_axis_no_mz_zoom.pdf')

print('\n\ntest #22:')
generate_plots_on_HRMS_data(
        query_data=f'{Path.cwd()}/data/lcms_query_library.csv',
        reference_data=f'{Path.cwd()}/data/lcms_reference_library.csv',
        wf_mz=2,
        wf_intensity=0.5,
        output_path=f'{Path.cwd()}/plots/test_wf_normalized_y_axis_no_mz_zoom.pdf')

print('\n\ntest #23:')
generate_plots_on_HRMS_data(
        query_data=f'{Path.cwd()}/data/lcms_query_library.csv',
        reference_data=f'{Path.cwd()}/data/lcms_reference_library.csv',
        y_axis_transformation='log10',
        output_path=f'{Path.cwd()}/plots/test_no_wf_log10_y_axis_no_mz_zoom.pdf')

print('\n\ntest #24:')
generate_plots_on_HRMS_data(
        query_data=f'{Path.cwd()}/data/lcms_query_library.csv',
        reference_data=f'{Path.cwd()}/data/lcms_reference_library.csv',
        y_axis_transformation='sqrt',
        output_path=f'{Path.cwd()}/plots/test_no_wf_sqrt_y_axis_no_mz_zoom.pdf')

print('\n\ntest #25:')
generate_plots_on_HRMS_data(
        query_data=f'{Path.cwd()}/data/lcms_query_library.csv',
        reference_data=f'{Path.cwd()}/data/lcms_reference_library.csv',
        mz_min = 400,
        mz_max = 650,
        y_axis_transformation='sqrt',
        output_path=f'{Path.cwd()}/plots/test_no_wf_sqrt_y_axis_mz_zoom.pdf')

print('\n\ntest #26:')
generate_plots_on_HRMS_data(
        query_data=f'{Path.cwd()}/data/lcms_query_library.csv',
        reference_data=f'{Path.cwd()}/data/lcms_reference_library.csv',
        high_quality_reference_library=False,
        output_path=f'{Path.cwd()}/plots/test_HRMS.pdf')

print('\n\ntest #27:')
generate_plots_on_NRMS_data(
        query_data=f'{Path.cwd()}/data/gcms_query_library.csv',
        reference_data=f'{Path.cwd()}/data/gcms_reference_library.csv',
        high_quality_reference_library=False,
        output_path=f'{Path.cwd()}/plots/test_NRMS.pdf')

print('\n\ntest #28:')
generate_plots_on_HRMS_data(
        query_data=f'{Path.cwd()}/data/lcms_query_library.csv',
        reference_data=f'{Path.cwd()}/data/lcms_reference_library.csv',
        similarity_measure='jaccard',
        output_path=f'{Path.cwd()}/plots/test28.pdf')

print('\n\ntest #28:')
generate_plots_on_HRMS_data(
        query_data=f'{Path.cwd()}/data/lcms_query_library.csv',
        reference_data=f'{Path.cwd()}/data/lcms_reference_library.csv',
        similarity_measure='hamming',
        output_path=f'{Path.cwd()}/plots/test28.pdf')

print('\n\ntest #29:')
generate_plots_on_NRMS_data(
        query_data=f'{Path.cwd()}/data/gcms_query_library.csv',
        reference_data=f'{Path.cwd()}/data/gcms_reference_library.csv',
        similarity_measure='sokal_sneath',
        output_path=f'{Path.cwd()}/plots/test29.pdf')

print('\n\ntest #30:')
generate_plots_on_NRMS_data(
        query_data=f'{Path.cwd()}/data/gcms_query_library.csv',
        reference_data=f'{Path.cwd()}/data/gcms_reference_library.csv',
        similarity_measure='simpson',
        output_path=f'{Path.cwd()}/plots/test30.pdf')

print('\n\ntest #31:')
generate_plots_on_NRMS_data(
        query_data=f'{Path.cwd()}/data/gcms_query_library.csv',
        reference_data=f'{Path.cwd()}/data/gcms_reference_library.csv',
        similarity_measure='mixture',
        weights={'Cosine':0.5, 'Shannon':0.3, 'Renyi':0.1, 'Tsallis':0.1},
        output_path=f'{Path.cwd()}/plots/test31.pdf')

print('\n\ntest #32:')
generate_plots_on_HRMS_data(
        query_data=f'{Path.cwd()}/data/lcms_query_library.csv',
        reference_data=f'{Path.cwd()}/data/lcms_reference_library.csv',
        similarity_measure='mixture',
        weights={'Cosine':0.1, 'Shannon':0.2, 'Renyi':0.3, 'Tsallis':0.4},
        output_path=f'{Path.cwd()}/plots/test32.pdf')

