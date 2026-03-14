from pycompound.plot_spectra import generate_plots_on_NRMS_data
from pycompound.spec_lib_matching import run_spec_lib_matching_on_NRMS_data
from pycompound.spec_lib_matching import tune_params_on_NRMS_data_grid
from pycompound.spec_lib_matching import tune_params_DE
from pathlib import Path
import os

path_to_query1 = f'{Path.cwd()}/tests/data/gcms_query.txt'
path_to_query2 = f'{Path.cwd()}/tests/data/gcms_query_tuning.txt'
path_to_ref = f'{Path.cwd()}/tests/data/trimmed_gcms_reference_library.txt'

print("############ Starting Toy Examples ############")
print()
print(">>> Plot Spectra")

##### plot spectra #####
generate_plots_on_NRMS_data(
        query_data = path_to_query1,
        reference_data = path_to_ref,
        similarity_measure = 'cosine',
        spectrum_ID1 = 'ID_1',
        spectrum_ID2 = '463-51-4',
        output_path = f'{Path.cwd()}/python_package_plotting_example.pdf')

print()
print("Run Spectral Library Matching")

##### run spectral library matching #####
run_spec_lib_matching_on_NRMS_data(
        query_data = path_to_query2,
        reference_data = path_to_ref,
        similarity_measure = 'cosine',
        print_id_results = True)

print()
print("Tune Parameters via Grid Search")

##### tune parameters via exhaustive grid search #####
tune_params_on_NRMS_data_grid(
        query_data = path_to_query2,
        reference_data = path_to_ref,
        grid={'wf_mz':[0.0,2.0], 'wf_int':[1.0,2.0]},
        output_path=f'{Path.cwd()}/test_grid_tuning.txt')

print()
print("Tune Parameters via DE Optimization")

##### tune parameters via differential evolution optimization #####
tune_params_DE(
        query_data = path_to_query2,
        reference_data = path_to_ref,
        chromatography_platform = 'NRMS',
        similarity_measure = 'cosine',
        optimize_params = ["wf_mz","wf_int"],
        param_bounds = {"wf_mz":(0.0,5.0),"wf_int":(0.0,5.0)},
        maxiters = 10,
        de_workers = 1)

print()
print("############ Completed Toy Examples ############")