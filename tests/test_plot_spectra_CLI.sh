#1/bin/bash

echo $'\n\n\ntest #0'
python ../src/pycompound/plot_spectra_CLI.py \
        --query_data ${PWD}/data/lcms_query_library.txt \
        --reference_data ${PWD}/data/trimmed_GNPS_reference_library.txt \
        --wf_mz 2 \
        --window_size_matching 0.3 \
        --chromatography_platform HRMS \
        --output_path ${PWD}/output_plotting_HRMS.pdf

echo $'\n\n\ntest #1'
python ../src/pycompound/plot_spectra_CLI.py \
        --query_data ${PWD}/data/gcms_query_library.txt \
        --reference_data ${PWD}/data/trimmed_gcms_reference_library.txt \
        --spectrum_ID1 463514 \
        --spectrum_ID2 112312 \
        --noise_threshold 0.1 \
        --LET_threshold 2 \
        --chromatography_platform NRMS \
        --output_path ${PWD}/output_plotting_NRMS_1.pdf

echo $'\n\n\ntest #2'
python ../src/pycompound/plot_spectra_CLI.py \
        --query_data ${PWD}/data/gcms_query_library.txt \
        --reference_data ${PWD}/data/trimmed_gcms_reference_library.txt \
        --similarity_measure mixture \
        --weights '{"Cosine":0.7,"Shannon":0.1,"Renyi":0.1,"Tsallis":0.1}' \
        --chromatography_platform NRMS \
        --output_path ${PWD}/output_plotting_NRMS_2.pdf

