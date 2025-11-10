
from pycompound.spec_lib_matching import run_spec_lib_matching_on_HRMS_data 
from pycompound.spec_lib_matching import run_spec_lib_matching_on_NRMS_data 
from pathlib import Path
import os


'''
print('\n\ntest #1:')
run_spec_lib_matching_on_HRMS_data(query_data=f'{Path.cwd()}/data/lcms_query_library.txt',
                                   reference_data=f'{Path.cwd()}/data/lcms_reference_library.txt',
                                   similarity_measure='cosine')

print('\n\ntest #2:')
run_spec_lib_matching_on_HRMS_data(query_data=f'{Path.cwd()}/data/lcms_query_library.txt',
                                   reference_data=f'{Path.cwd()}/data/lcms_reference_library.txt',
                                   spectrum_preprocessing_order='MAB')

print('\n\ntest #3:')
run_spec_lib_matching_on_HRMS_data(query_data=f'{Path.cwd()}/data/lcms_query_library.txt',
                                   reference_data=f'{Path.cwd()}/data/lcms_reference_library.txt',
                                   window_size_centroiding='small')

print('\n\ntest #4:')
run_spec_lib_matching_on_HRMS_data(query_data=f'{Path.cwd()}/data/lcms_query_library.txt',
                                   reference_data=f'{Path.cwd()}/data/lcms_reference_library.txt',
                                   mz_min='hello')

print('\n\ntest #5:')
run_spec_lib_matching_on_HRMS_data(query_data=f'{Path.cwd()}/data/lcms_query_library.txt',
                                   reference_data=f'{Path.cwd()}/data/lcms_reference_library.txt',
                                   window_size_matching='big')

print('\n\ntest #6:')
run_spec_lib_matching_on_HRMS_data(query_data=f'{Path.cwd()}/data/lcms_query_library.txt',
                                   reference_data=f'{Path.cwd()}/data/lcms_reference_library.txt',
                                   noise_threshold='world')

print('\n\ntest #7:')
run_spec_lib_matching_on_HRMS_data(query_data=f'{Path.cwd()}/data/lcms_query_library.txt',
                                   reference_data=f'{Path.cwd()}/data/lcms_reference_library.txt',
                                   noise_threshold=-0.5)

print('\n\ntest #8:')
run_spec_lib_matching_on_HRMS_data(query_data=f'{Path.cwd()}/data/lcms_query_library.txt',
                                   reference_data=f'{Path.cwd()}/data/lcms_reference_library.txt',
                                   wf_intensity='a')

print('\n\ntest #9:')
run_spec_lib_matching_on_HRMS_data(query_data=f'{Path.cwd()}/data/lcms_query_library.txt',
                                   reference_data=f'{Path.cwd()}/data/lcms_reference_library.txt',
                                   entropy_dimension=-1)
print('\n\ntest #10:')
run_spec_lib_matching_on_HRMS_data(query_data=f'{Path.cwd()}/data/lcms_query_library.txt',
                                   reference_data=f'{Path.cwd()}/data/lcms_reference_library.txt',
                                   similarity_measure='renyi',
                                   wf_mz=2,
                                   wf_intensity=0.5,
                                   entropy_dimension=2,
                                   n_top_matches_to_save=3,
                                   print_id_results='asdf')

print('\n\ntest #12:')
run_spec_lib_matching_on_HRMS_data(query_data=f'{Path.cwd()}/data/lcms_query_library.txt',
                                   reference_data=f'{Path.cwd()}/data/lcms_reference_library.txt',
                                   similarity_measure='tsallis',
                                   wf_mz=2,
                                   wf_intensity=0.5,
                                   entropy_dimension=2,
                                   n_top_matches_to_save=3,
                                   print_id_results=True)

print('\n\ntest #13:')
run_spec_lib_matching_on_HRMS_data(query_data=f'{Path.cwd()}/data/lcms_query_library.txt',
                                   reference_data=f'{Path.cwd()}/data/lcms_reference_library.txt',
                                   similarity_measure='shannon',
                                   wf_intensity=0.5,
                                   n_top_matches_to_save=2,
                                   print_id_results=True)

print('\n\ntest #14:')
run_spec_lib_matching_on_HRMS_data(query_data=f'{Path.cwd()}/data/lcms_query_library.txt',
                                   reference_data=f'{Path.cwd()}/data/lcms_reference_library.txt',
                                   similarity_measure='renyi',
                                   wf_mz=2,
                                   wf_intensity=0.5,
                                   entropy_dimension=2,
                                   n_top_matches_to_save=3,
                                   print_id_results=True)


print('\n\ntest #15:')
run_spec_lib_matching_on_HRMS_data(query_data=f'{Path.cwd()}/data/lcms_query_library.txt',
                                   reference_data=f'{Path.cwd()}/data/lcms_reference_library.txt',
                                   similarity_measure='cosine',
                                   wf_mz=2,
                                   wf_intensity=0.5,
                                   n_top_matches_to_save=3)

print('\n\ntest #16:')
run_spec_lib_matching_on_HRMS_data(query_data=f'{Path.cwd()}/data/lcms_query_library.txt',
                                   reference_data=f'{Path.cwd()}/data/lcms_reference_library.txt',
                                   similarity_measure='cosine',
                                   spectrum_preprocessing_order='LFWNCM')

print('\n\ntest #17:')
run_spec_lib_matching_on_HRMS_data(query_data=f'{Path.cwd()}/data/lcms_query_library.txt',
                                   reference_data=f'{Path.cwd()}/data/lcms_reference_library.txt',
                                   similarity_measure='cosine',
                                   mz_min=100,
                                   mz_max=400,
                                   int_min=50,
                                   int_max=150000)

print('\n\ntest #18:')
run_spec_lib_matching_on_HRMS_data(query_data=f'{Path.cwd()}/data/lcms_query_library.txt',
                                   reference_data=f'{Path.cwd()}/data/lcms_reference_library.txt',
                                   window_size_centroiding=0.1,
                                   window_size_matching=0.05)

print('\n\ntest #19:')
run_spec_lib_matching_on_HRMS_data(query_data=f'{Path.cwd()}/data/lcms_query_library.txt',
                                   reference_data=f'{Path.cwd()}/data/lcms_reference_library.txt',
                                   noise_threshold=0.1)

print('\n\ntest #20:')
run_spec_lib_matching_on_HRMS_data(query_data=f'{Path.cwd()}/data/lcms_query_library.txt',
                                   reference_data=f'{Path.cwd()}/data/lcms_reference_library.txt',
                                   LET_threshold=3)

print('\n\ntest #21:')
run_spec_lib_matching_on_HRMS_data(query_data=f'{Path.cwd()}/data/lcms_query_library.txt',
                                   reference_data=f'{Path.cwd()}/data/lcms_reference_library.txt',
                                   spectrum_preprocessing_order='WMC')

print('\n\ntest #22:')
run_spec_lib_matching_on_HRMS_data(query_data=f'{Path.cwd()}/data/lcms_query_library.txt',
                                   reference_data=f'{Path.cwd()}/data/lcms_reference_library.txt',
                                   spectrum_preprocessing_order='ML')

print('\n\ntest #23:')
run_spec_lib_matching_on_HRMS_data(query_data=f'{Path.cwd()}/data/lcms_query_library.txt',
                                   reference_data=f'{Path.cwd()}/data/lcms_reference_library.txt',
                                   noise_threshold=0.1,
                                   spectrum_preprocessing_order='MNL')

print('\n\ntest #24:')
run_spec_lib_matching_on_NRMS_data(query_data=f'{Path.cwd()}/data/gcms_query_library_tuning.txt',
                                   reference_data=f'{Path.cwd()}/data/gcms_reference_library.txt',
                                   print_id_results=True)

print('\n\ntest #25:')
run_spec_lib_matching_on_NRMS_data(query_data=f'{Path.cwd()}/data/gcms_query_library.txt',
                                   reference_data=f'{Path.cwd()}/data/gcms_reference_library.txt',
                                   spectrum_preprocessing_order='MAB')

print('\n\ntest #26:')
run_spec_lib_matching_on_NRMS_data(query_data=f'{Path.cwd()}/data/gcms_query_library.txt',
                                   reference_data=f'{Path.cwd()}/data/gcms_reference_library.txt',
                                   mz_min='hello')

print('\n\ntest #27:')
run_spec_lib_matching_on_NRMS_data(query_data=f'{Path.cwd()}/data/gcms_query_library.txt',
                                   reference_data=f'{Path.cwd()}/data/gcms_reference_library.txt',
                                   noise_threshold='world')

print('\n\ntest #28:')
run_spec_lib_matching_on_NRMS_data(query_data=f'{Path.cwd()}/data/gcms_query_library.txt',
                                   reference_data=f'{Path.cwd()}/data/gcms_reference_library.txt',
                                   noise_threshold=-0.5)

print('\n\ntest #29:')
run_spec_lib_matching_on_NRMS_data(query_data=f'{Path.cwd()}/data/gcms_query_library.txt',
                                   reference_data=f'{Path.cwd()}/data/gcms_reference_library.txt',
                                   wf_intensity='a')

print('\n\ntest #30:')
run_spec_lib_matching_on_NRMS_data(query_data=f'{Path.cwd()}/data/gcms_query_library.txt',
                                   reference_data=f'{Path.cwd()}/data/gcms_reference_library.txt',
                                   entropy_dimension=-1)

print('\n\ntest #33:')
run_spec_lib_matching_on_NRMS_data(query_data=f'{Path.cwd()}/data/gcms_query_library.txt',
                                   reference_data=f'{Path.cwd()}/data/gcms_reference_library.txt',
                                   similarity_measure='shannon')

print('\n\ntest #33:')
run_spec_lib_matching_on_NRMS_data(query_data=f'{Path.cwd()}/data/gcms_query_library.txt',
                                   reference_data=f'{Path.cwd()}/data/gcms_reference_library.txt',
                                   spectrum_preprocessing_order='NF',
                                   mz_min=50,
                                   mz_max=600,
                                   noise_threshold=0.1)

print('\n\ntest #34:')
run_spec_lib_matching_on_NRMS_data(query_data=f'{Path.cwd()}/data/gcms_query_library.txt',
                                   reference_data=f'{Path.cwd()}/data/gcms_reference_library.txt',
                                   spectrum_preprocessing_order='FN',
                                   mz_min=50,
                                   mz_max=600,
                                   noise_threshold=0.1)

print('\n\ntest #35:')
run_spec_lib_matching_on_NRMS_data(query_data=f'{Path.cwd()}/data/gcms_query_library.txt',
                                   reference_data=f'{Path.cwd()}/data/gcms_reference_library.txt',
                                   spectrum_preprocessing_order='FNW',
                                   wf_mz=1.5,
                                   mz_min=50,
                                   mz_max=600,
                                   noise_threshold=0.1)

print('\n\ntest #36:')
run_spec_lib_matching_on_NRMS_data(query_data=f'{Path.cwd()}/data/gcms_query_library.txt',
                                   reference_data=f'{Path.cwd()}/data/gcms_reference_library.txt',
                                   spectrum_preprocessing_order='FNW',
                                   similarity_measure='renyi',
                                   entropy_dimension=1.5,
                                   wf_mz=1.5,
                                   mz_min=50,
                                   mz_max=600,
                                   noise_threshold=0.1)

print('\n\ntest #37:')
run_spec_lib_matching_on_NRMS_data(query_data=f'{Path.cwd()}/data/gcms_query_library.txt',
                                   reference_data=f'{Path.cwd()}/data/gcms_reference_library.txt',
                                   high_quality_reference_library=True,
                                   noise_threshold=0.1)

print('\n\ntest #38:')
run_spec_lib_matching_on_NRMS_data(query_data=f'{Path.cwd()}/data/gcms_query_library.txt',
                                   reference_data=f'{Path.cwd()}/data/gcms_reference_library.txt',
                                   high_quality_reference_library=False,
                                   noise_threshold=0.1)

print('\n\ntest #39:')
run_spec_lib_matching_on_NRMS_data(query_data=f'{Path.cwd()}/data/gcms_query_library.txt',
                                   reference_data=f'{Path.cwd()}/data/gcms_reference_library.txt',
                                   likely_reference_ids=f'{Path.cwd()}/data/likely_gcms_ids.txt')

print('\n\ntest #40:')
run_spec_lib_matching_on_HRMS_data(query_data=f'{Path.cwd()}/data/lcms_query_library.txt',
                                   reference_data=f'{Path.cwd()}/data/lcms_reference_library.txt',
                                   likely_reference_ids=f'{Path.cwd()}/data/likely_lcms_ids.txt')

#print('\n\ntest #41:')
#run_spec_lib_matching_on_NRMS_data(query_data=f'{Path.cwd()}/data/n_Alkane+76_mixtures_1_.cdf',
#                                   reference_data=f'{Path.cwd()}/data/gcms_reference_library.txt')

print('\n\ntest #42:')
run_spec_lib_matching_on_HRMS_data(query_data=f'{Path.cwd()}/data/lcms_query_library.txt',
                                   reference_data=f'{Path.cwd()}/data/lcms_reference_library.txt')

print('\n\ntest #43:')
run_spec_lib_matching_on_HRMS_data(query_data=f'{Path.cwd()}/data/lcms_query_library.txt',
                                   reference_data=[f'{Path.cwd()}/data/lcms_reference_library.txt',f'{Path.cwd()}/data/GNPS-SELLECKCHEM-FDA-PART1.mgf'])

print('\n\ntest #44:')
run_spec_lib_matching_on_NRMS_data(query_data=f'{Path.cwd()}/data/gcms_query_library.txt',
                                   reference_data=f'{Path.cwd()}/data/gcms_reference_library.txt')

print('\n\ntest #45:')
run_spec_lib_matching_on_NRMS_data(query_data=f'{Path.cwd()}/data/gcms_query_library.txt',
                                   reference_data=[f'{Path.cwd()}/data/gcms_reference_library.txt',f'{Path.cwd()}/data/gcms_reference_library.txt'])

print('\n\ntest #46:')
run_spec_lib_matching_on_HRMS_data(query_data=f'{Path.cwd()}/data/lcms_query_library.txt',
                                   reference_data=f'{Path.cwd()}/data/lcms_reference_library.txt',
                                   similarity_measure='jaccard')

print('\n\ntest #47:')
run_spec_lib_matching_on_NRMS_data(query_data=f'{Path.cwd()}/data/gcms_query_library.txt',
                                   reference_data=f'{Path.cwd()}/data/gcms_reference_library.txt',
                                   similarity_measure='hellinger')

print('\n\ntest #48:')
run_spec_lib_matching_on_NRMS_data(query_data=f'{Path.cwd()}/data/gcms_query_library.txt',
                                   reference_data=f'{Path.cwd()}/data/gcms_reference_library.txt',
                                   similarity_measure='mixture')

print('\n\ntest #49:')
run_spec_lib_matching_on_NRMS_data(query_data=f'{Path.cwd()}/data/gcms_query_library.txt',
                                   reference_data=f'{Path.cwd()}/data/gcms_reference_library.txt',
                                   similarity_measure='mixture',
                                   weights={'Cosine':0.7, 'Shannon':0.1, 'Renyi':0.1, 'Tsallis':0.1})

print('\n\ntest #50:')
run_spec_lib_matching_on_HRMS_data(query_data=f'{Path.cwd()}/data/1min.mzML',
                                   reference_data=f'{Path.cwd()}/data/lcms_reference_library.txt')

print('\n\ntest #51:')
run_spec_lib_matching_on_HRMS_data(query_data=f'{Path.cwd()}/data/GNPS-SELLECKCHEM-FDA-PART1.mgf',
                                   reference_data=f'{Path.cwd()}/data/full_GNPS_reference_library.txt',
                                   ionization_mode='Positive',
                                   adduct='H',
                                   precursor_ion_mz_tolerance=1.0,
                                   print_id_results=True)

print('\n\ntest #52:')
run_spec_lib_matching_on_HRMS_data(query_data=f'{Path.cwd()}/data/GNPS-NIH-SMALLMOLECULEPHARMACOLOGICALLYACTIVE.json',
                                   reference_data=f'{Path.cwd()}/data/full_GNPS_reference_library.txt',
                                   ionization_mode='Negative',
                                   adduct='H',
                                   precursor_ion_mz_tolerance=0.5,
                                   print_id_results=True)

print('\n\ntest #53:')
run_spec_lib_matching_on_HRMS_data(query_data=f'{Path.cwd()}/data/GNPS-NIH-SMALLMOLECULEPHARMACOLOGICALLYACTIVE.json',
                                   reference_data=f'{Path.cwd()}/data/full_GNPS_reference_library.txt',
                                   adduct='NH3',
                                   print_id_results=True)

print('\n\ntest #54:')
run_spec_lib_matching_on_HRMS_data(query_data=f'{Path.cwd()}/data/MoNA-export-Human_Plasma_Quant.msp',
                                   reference_data=f'{Path.cwd()}/data/full_GNPS_reference_library.txt',
                                   ionization_mode='Positive',
                                   adduct='H',
                                   precursor_ion_mz_tolerance=0.5,
                                   print_id_results=True)
'''

