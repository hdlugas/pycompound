#!/usr/bin/env bash
QUERY_PATH1=${PWD}/tests/data/gcms_query.txt
QUERY_PATH2=${PWD}/tests/data/gcms_query_tuning.txt
REF_PATH=${PWD}/tests/data/trimmed_gcms_reference_library.txt

echo "############ Starting Toy Examples ############"
echo ""
echo ">>> Plot Spectra"

##### plot spectra #####
python src/pycompound/plot_spectra_CLI.py \
        --query_data $QUERY_PATH1 \
        --reference_data $REF_PATH \
        --similarity_measure cosine \
        --chromatography_platform NRMS \
        --spectrum_ID1 "ID_1" \
        --spectrum_ID2 "463-51-4" \
        --output_path ${PWD}/CLI_plotting_example.pdf

echo ""
echo ">>> Run Spectral library Matching" 

##### run spectral library matching #####
python src/pycompound/spec_lib_matching_CLI.py \
        --query_data $QUERY_PATH2 \
        --reference_data $REF_PATH \
        --chromatography_platform NRMS \
        --output_identification ${PWD}/CLI_identification_output_example.txt \
        --output_similarity_scores ${PWD}/CLI_similarity_scores_output_example.txt

echo ""
echo ">>> Tune Parameters via Grid Search"

##### tune parameters via exhaustive grid search #####
python src/pycompound/tuning_CLI_grid.py \
        --query_data $QUERY_PATH2 \
        --reference_data $REF_PATH \
        --chromatography_platform NRMS \
        --wf_int 1,2 \
        --wf_mz 0,2 \
        --output_path ${PWD}/CLI_grid_tuning_output_example.txt

echo ""
echo ">>> Tune Parameters via DE Optimization"

##### tune parameters via differential evolution optimization #####
python src/pycompound/tuning_CLI_DE.py \
        --query_data $QUERY_PATH2 \
        --reference_data $REF_PATH \
        --chromatography_platform NRMS \
        --opt wf_mz wf_int \
        --bound wf_mz=0.0:5.0 \
        --bound wf_int=0.0:5.0 \
        --maxiter 10 \
        --workers 5

echo ""
echo "############ Completed Toy Examples ############"
