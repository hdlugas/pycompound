#!/bin/bash

echo $'\n\n\ntest #1'
pycompound_transform_spectra \
        --spectra_data ${PWD}/data/GNPS-SELLECKCHEM-FDA-PART1.mgf \
        --wf_mz 2 \
        --spectrum_preprocessing_order W \
        --output_path ${PWD}/data/example_transformed_spectra_CLI1.mgf

echo $'\n\n\ntest #2'
pycompound_transform_spectra \
        --spectra_data ${PWD}/data/lcms_query.txt \
        --LET_threshold 3 \
        --wf_mz 1.1 \
        --spectrum_preprocessing_order LW \
        --output_path ${PWD}/data/example_transformed_spectra_CLI2.txt

echo $'\n\n\ntest #3'
pycompound_transform_spectra \
        --spectra_data ${PWD}/data/MoNA-export-Human_Plasma_Quant.msp \
        --noise_threshold 0.075 \
        --spectrum_preprocessing_order N \
        --output_path ${PWD}/data/example_transformed_spectra_CLI3.msp

