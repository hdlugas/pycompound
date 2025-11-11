#!/bin/bash

cd ${PWD}/../src/pycompound

echo $'\n\n\ntest #0'
python spec_lib_matching_CLI.py \
        --query_data ../../tests/data/gcms_query_library.txt \
        --reference_data ../../tests/data/trimmed_gcms_reference_library.txt \
        --chromatography_platform NRMS \
        --output_identification ${PWD}/../../tests/output_identification_NRMS.txt \
        --output_similarity_scores ${PWD}/../../tests/output_similarity_scores_NRMS.txt


echo $'\n\n\ntest #1'
python spec_lib_matching_CLI.py \
        --query_data ../../tests/data/lcms_query_library.txt \
        --reference_data ../../tests/data/trimmed_GNPS_reference_library.txt \
        --chromatography_platform HRMS \
        --output_identification ${PWD}/../../tests/output_identification_HRMS_1.txt \
        --output_similarity_scores ${PWD}/../../tests/output_similarity_scores_HRMS_1.txt


echo $'\n\n\ntest #2'
python spec_lib_matching_CLI.py \
        --query_data ../../tests/data/lcms_query_library.txt \
        --reference_data ../../tests/data/trimmed_GNPS_reference_library.txt \
        --chromatography_platform HRMS \
        --similarity_measure mixture \
        --weights '{"Cosine":0.7,"Shannon":0.1,"Renyi":0.1,"Tsallis":0.1}' \
        --output_identification ${PWD}/../../tests/output_identification_HRMS_2.txt \
        --output_similarity_scores ${PWD}/../../tests/output_similarity_scores_HRMS_2.txt


echo $'\n\n\ntest #3'
python spec_lib_matching_CLI.py \
        --query_data ../../tests/data/lcms_query_library.txt \
        --reference_data ../../tests/data/trimmed_GNPS_reference_library.txt \
        --chromatography_platform HRMS \
        --precursor_ion_mz 0.8 \
        --ionization_mode Positive \
        --adduct H \
        --similarity_measure shannon \
        --output_identification ${PWD}/../../tests/output_identification_HRMS_3.txt \
        --output_similarity_scores ${PWD}/../../tests/output_similarity_scores_HRMS_3.txt

