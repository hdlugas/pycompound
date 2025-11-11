#!/bin/bash


echo $'\n\n\ntest #1'
python ../src/pycompound/tuning_CLI_grid.py \
        --query_data ${PWD}/data/lcms_query_library_tuning.txt \
        --reference_data ${PWD}/data/full_GNPS_reference_library.txt \
        --precursor_ion_mz_tolerance 0.1 \
        --ionization_mode Positive \
        --adduct H \
        --window_size_matching 0.5 \
        --chromatography_platform HRMS \
        --output_path ${PWD}/output_tuning_HRMS_CLI_1.txt

echo $'\n\n\ntest #2'
python ../src/pycompound/tuning_CLI_grid.py \
        --query_data ${PWD}/data/lcms_query_library_tuning.txt \
        --reference_data ${PWD}/data/full_GNPS_reference_library.txt \
        --precursor_ion_mz_tolerance 0.2 \
        --ionization_mode Positive \
        --adduct H \
        --window_size_matching 0.5,0.1,0.05 \
        --chromatography_platform HRMS \
        --output_path ${PWD}/output_tuning_HRMS_CLI_2.txt

echo $'\n\n\ntest #3'
python ../src/pycompound/tuning_CLI_grid.py \
        --query_data ${PWD}/data/lcms_query_library_tuning.txt \
        --reference_data ${PWD}/data/lcms_reference_library.txt \
        --window_size_matching 0.5 \
        --chromatography_platform HRMS \
        --output_path ${PWD}/output_tuning_HRMS_CLI_3.txt \

echo $'\n\n\ntest #4'
python ../src/pycompound/tuning_CLI_grid.py \
        --query_data ${PWD}/data/gcms_query_library_tuning.txt \
        --reference_data ${PWD}/data/gcms_reference_library.txt \
        --similarity_measure cosine,shannon,renyi \
        --wf_mz 0,2,3 \
        --noise_threshold 0,0.1 \
        --chromatography_platform NRMS \
        --output_path ${PWD}/output_tuning_NRMS_CLI.txt \

echo $'\n\n\ntest #5'
python ../src/pycompound/tuning_CLI_DE.py \
  --chromatography_platform NRMS \
  --query_data ${PWD}/data/gcms_query_library_tuning.txt \
  --reference_data ${PWD}/data/gcms_reference_library.txt \
  --similarity_measure cosine \
  --opt noise_threshold wf_mz \
  --bound noise_threshold=0.0:0.20 \
  --bound wf_mz=0.0:5.0 \
  --maxiter 3 \
  --seed 1 \
  --workers 4

echo $'\n\n\ntest #6'
python ../src/pycompound/tuning_CLI_DE.py \
  --chromatography_platform HRMS \
  --query_data ${PWD}/data/lcms_query_library_tuning.txt \
  --reference_data ${PWD}/data/full_GNPS_reference_library.txt \
  --precursor_ion_mz_tolerance 0.1 \
  --ionization_mode Positive \
  --adduct H \
  --similarity_measure cosine \
  --opt window_size_centroiding noise_threshold wf_mz \
  --bound window_size_centroiding=0.0:0.4 \
  --bound noise_threshold=0.0:0.20 \
  --bound wf_mz=0.0:5.0 \
  --maxiter 3 \
  --seed 1 \
  --workers 5


