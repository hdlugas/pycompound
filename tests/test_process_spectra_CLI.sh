#!/bin/bash

echo $'\n\n\ntest #1'
pycompound_transform_spectra \
        --spectra_data ${PWD}/data/lcms_query.txt \
        --wf_mz 2 \
        --spectra_preprocessing_transformation W \
        --output_path ${PWD}/example_transformed_spectra_CLI1.txt

