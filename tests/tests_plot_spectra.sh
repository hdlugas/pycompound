#!/bin/bash

cd ${PWD}/../src

: <<'COMMENT'

COMMENT

echo $'\n\n\n\n\ntest #0'
python plot_spectra.py \
  --query_data ../data/lcms_query_library.csv \
  --reference_data ../data/lcms_reference_library.csv \
  --similarity_measure tsallis \
  --chromatography_platform HRMS \
  --wf_mz 2 \
  --wf_intensity 0.5 \
  --normalization_method standard \
  --entropy_dimension 2 \
  --save_plots ../tests/test0.pdf

echo $'\n\n\n\n\ntest #1'
python plot_spectra.py \
  --query_data ../data/lcms_query_library.csv \
  --reference_data ../data/lcms_reference_library.csv \
  --similarity_measure tsallis \
  --chromatography_platform HRMS \
  --wf_mz 2 \
  --wf_intensity 0.5 \
  --normalization_method standard \
  --entropy_dimension 2 \
  --save_plots ../tests/test1.pdf

echo $'\n\n\n\n\ntest #2'
python plot_spectra.py \
  --query_data ../data/lcms_query_library.csv \
  --reference_data ../data/lcms_reference_library.csv \
  --similarity_measure tsallis \
  --chromatography_platform HRMS \
  --wf_mz 2 \
  --wf_intensity 0.5 \
  --normalization_method standard \
  --entropy_dimension 0.5 \
  --save_plots ../tests/test2.pdf

echo $'\n\n\n\n\ntest #3'
python plot_spectra.py \
  --query_data ../data/lcms_query_library.csv \
  --reference_data ../data/lcms_reference_library.csv \
  --similarity_measure tsallis \
  --chromatography_platform HRMS \
  --normalization_method standard \
  --entropy_dimension 2 \
  --save_plots ../tests/test3.pdf

echo $'\n\n\n\n\ntest #4'
python plot_spectra.py \
  --query_data ../data/lcms_query_library.csv \
  --reference_data ../data/lcms_reference_library.csv \
  --similarity_measure renyi \
  --chromatography_platform HRMS \
  --normalization_method standard \
  --entropy_dimension 1.1 \
  --save_plots ../tests/test4.pdf


echo $'\n\n\n\n\ntest #5'
python plot_spectra.py \
  --query_data ../data/lcms_query_library.csv \
  --reference_data ../data/lcms_reference_library.csv \
  --similarity_measure renyi \
  --chromatography_platform HRMS \
  --normalization_method standard \
  --entropy_dimension 0.9 \
  --save_plots ../tests/test5.pdf

echo $'\n\n\n\n\ntest #6'
python plot_spectra.py \
  --query_data ../data/lcms_query_library.csv \
  --reference_data ../data/lcms_reference_library.csv \
  --similarity_measure shannon \
  --chromatography_platform HRMS \
  --normalization_method standard \
  --save_plots ../tests/test6.pdf

echo $'\n\n\n\n\ntest #7'
python plot_spectra.py \
  --query_data ../data/lcms_query_library.csv \
  --reference_data ../data/lcms_reference_library.csv \
  --similarity_measure cosine \
  --chromatography_platform HRMS \
  --save_plots ../tests/test7.pdf

echo $'\n\n\n\n\ntest #8'
python plot_spectra.py \
  --query_data ../data/lcms_query_library.csv \
  --reference_data ../data/lcms_reference_library.csv \
  --spectrum_preprocessing_order LFWNCM \
  --similarity_measure cosine \
  --chromatography_platform HRMS \
  --wf_mz 0.5 \
  --save_plots ../tests/test8.pdf

echo $'\n\n\n\n\ntest #9'
python plot_spectra.py \
  --query_data ../data/lcms_query_library.csv \
  --reference_data ../data/lcms_reference_library.csv \
  --spectrum_preprocessing_order LFWNCM \
  --similarity_measure cosine \
  --chromatography_platform HRMS \
  --mz_min 200 \
  --mz_max 250 \
  --int_min 50 \
  --int_max 500 \
  --save_plots ../tests/test9.pdf

echo $'\n\n\n\n\ntest #10'
python plot_spectra.py \
  --query_data ../data/lcms_query_library.csv \
  --reference_data ../data/lcms_reference_library.csv \
  --spectrum_preprocessing_order FWCNLM \
  --similarity_measure cosine \
  --chromatography_platform HRMS \
  --mz_max 100 \
  --save_plots ../tests/test10.pdf

echo $'\n\n\n\n\ntest #11'
python plot_spectra.py \
  --query_data ../data/lcms_query_library.csv \
  --reference_data ../data/lcms_reference_library.csv \
  --spectrum_preprocessing_order FWCNLM \
  --similarity_measure cosine \
  --chromatography_platform HRMS \
  --int_max 300 \
  --save_plots ../tests/test11.pdf

echo $'\n\n\n\n\ntest #12'
python plot_spectra.py \
  --query_data ../data/lcms_query_library.csv \
  --reference_data ../data/lcms_reference_library.csv \
  --spectrum_preprocessing_order FWCNLM \
  --similarity_measure cosine \
  --chromatography_platform HRMS \
  --int_min 100 \
  --save_plots ../tests/test12.pdf

echo $'\n\n\n\n\ntest #13'
python plot_spectra.py \
  --query_data ../data/lcms_query_library.csv \
  --reference_data ../data/lcms_reference_library.csv \
  --spectrum_preprocessing_order FWCNLM \
  --similarity_measure cosine \
  --chromatography_platform HRMS \
  --window_size_centroiding 0.1 \
  --save_plots ../tests/test13.pdf

echo $'\n\n\n\n\ntest #14'
python plot_spectra.py \
  --query_data ../data/lcms_query_library.csv \
  --reference_data ../data/lcms_reference_library.csv \
  --spectrum_preprocessing_order FWCNLM \
  --similarity_measure cosine \
  --chromatography_platform HRMS \
  --window_size_matching 0.1 \
  --save_plots ../tests/test14.pdf

echo $'\n\n\n\n\ntest #15'
python plot_spectra.py \
  --query_data ../data/lcms_query_library.csv \
  --reference_data ../data/lcms_reference_library.csv \
  --spectrum_preprocessing_order FWCNLM \
  --similarity_measure cosine \
  --chromatography_platform HRMS \
  --noise_threshold 0.4 \
  --save_plots ../tests/test15.pdf

echo $'\n\n\n\n\ntest #16'
python plot_spectra.py \
  --query_data ../data/lcms_query_library.csv \
  --reference_data ../data/lcms_reference_library.csv \
  --spectrum_preprocessing_order FWCNLM \
  --similarity_measure cosine \
  --chromatography_platform HRMS \
  --LET_threshold 3 \
  --save_plots ../tests/test16.pdf

echo $'\n\n\n\n\ntest #17'
python plot_spectra.py \
  --query_data ../data/lcms_query_library.csv \
  --reference_data ../data/lcms_reference_library.csv \
  --spectrum_preprocessing_order CMWNLF \
  --similarity_measure cosine \
  --chromatography_platform HRMS \
  --LET_threshold 2 \
  --save_plots ../tests/test17.pdf

echo $'\n\n\n\n\ntest #18'
python plot_spectra.py \
  --query_data ../data/lcms_query_library.csv \
  --reference_data ../data/lcms_reference_library.csv \
  --spectrum_preprocessing_order CMWNLF \
  --similarity_measure cosine \
  --chromatography_platform HRMS \
  --window_size_centroiding 0.05 \
  --wf_mz 0.5 \
  --wf_int 1.3 \
  --LET_threshold 2 \
  --save_plots ../tests/test18.pdf

echo $'\n\n\n\n\ntest #19'
python plot_spectra.py \
  --query_data ../data/lcms_query_library.csv \
  --reference_data ../data/lcms_reference_library.csv \
  --spectrum_preprocessing_order WCNLMF \
  --similarity_measure shannon \
  --chromatography_platform HRMS \
  --wf_mz 0.5 \
  --wf_int 1.3 \
  --LET_threshold 2 \
  --save_plots ../tests/test19.pdf

echo $'\n\n\n\n\ntest #20'
python plot_spectra.py \
  --query_data ../data/lcms_query_library.csv \
  --reference_data ../data/lcms_reference_library.csv \
  --spectrum_preprocessing_order WCNLFM \
  --similarity_measure renyi \
  --chromatography_platform HRMS \
  --wf_mz 0.5 \
  --wf_int 1.3 \
  --LET_threshold 2 \
  --save_plots ../tests/test20.pdf

echo $'\n\n\n\n\ntest #21'
python plot_spectra.py \
  --query_data ../data/lcms_query_library.csv \
  --reference_data ../data/lcms_reference_library.csv \
  --spectrum_preprocessing_order LWCMNF \
  --similarity_measure tsallis \
  --chromatography_platform HRMS \
  --wf_mz 0.5 \
  --wf_int 1.3 \
  --LET_threshold 2 \
  --save_plots ../tests/test21.pdf

echo $'\n\n\n\n\ntest #22'
python plot_spectra.py \
  --query_data ../data/lcms_query_library.csv \
  --reference_data ../data/lcms_reference_library.csv \
  --spectrum_preprocessing_order CMWNLF \
  --similarity_measure tsallis \
  --chromatography_platform HRMS \
  --LET_threshold 2 \
  --save_plots ../tests/test22.pdf

echo $'\n\n\n\n\ntest #23'
python plot_spectra.py \
  --query_data ../data/lcms_query_library.csv \
  --reference_data ../data/lcms_reference_library.csv \
  --spectrum_preprocessing_order ML \
  --similarity_measure tsallis \
  --chromatography_platform HRMS \
  --LET_threshold 3 \
  --save_plots ../tests/test23.pdf

echo $'\n\n\n\n\ntest #24'
python plot_spectra.py \
  --query_data ../data/lcms_query_library.csv \
  --reference_data ../data/lcms_reference_library.csv \
  --spectrum_preprocessing_order MNL \
  --similarity_measure tsallis \
  --chromatography_platform HRMS \
  --noise_threshold 0.1 \
  --LET_threshold 3 \
  --save_plots ../tests/test24.pdf

echo $'\n\n\n\n\ntest #25'
python plot_spectra.py \
  --query_data ../data/lcms_query_library.csv \
  --reference_data ../data/lcms_reference_library.csv \
  --spectrum_preprocessing_order MW \
  --similarity_measure tsallis \
  --chromatography_platform HRMS \
  --wf_mz 0.55 \
  --wf_int 1.1 \
  --save_plots ../tests/test25.pdf

echo $'\n\n\n\n\ntest #26'
python plot_spectra.py \
  --query_data ../data/lcms_query_library.csv \
  --reference_data ../data/lcms_reference_library.csv \
  --spectrum_preprocessing_order WM \
  --similarity_measure tsallis \
  --chromatography_platform HRMS \
  --wf_mz 0.55 \
  --wf_int 1.1 \
  --save_plots ../tests/test26.pdf

echo $'\n\n\n\n\ntest #27'
python plot_spectra.py \
  --query_data ../data/lcms_query_library.csv \
  --reference_data ../data/lcms_reference_library.csv \
  --spectrum_preprocessing_order FWL \
  --similarity_measure tsallis \
  --chromatography_platform HRMS \
  --wf_mz 0.55 \
  --wf_int 1.1 \
  --LET_threshold 4 \
  --save_plots ../tests/test27.pdf

echo $'\n\n\n\n\ntest #28'
python plot_spectra.py \
  --query_data ../data/lcms_query_library.csv \
  --reference_data ../data/lcms_reference_library.csv \
  --spectrum_preprocessing_order MCWL \
  --similarity_measure tsallis \
  --chromatography_platform HRMS \
  --wf_mz 0.55 \
  --wf_int 1.1 \
  --LET_threshold 4 \
  --save_plots ../tests/test27.pdf

echo $'\n\n\n\n\ntest #29'
python plot_spectra.py \
  --query_data ../data/lcms_query_library.csv \
  --reference_data ../data/lcms_reference_library.csv \
  --similarity_measure tsallis \
  --chromatography_platform HRMS \
  --high_quality_reference_library True \
  --wf_mz 0.55 \
  --wf_int 1.1 \
  --LET_threshold 3 \
  --save_plots ../tests/test29.pdf

echo $'\n\n\n\n\ntest #30'
python plot_spectra.py \
  --query_data ../data/lcms_query_library.csv \
  --reference_data ../data/lcms_reference_library.csv \
  --spectrum_ID1 212 \
  --similarity_measure tsallis \
  --chromatography_platform HRMS \
  --high_quality_reference_library True \
  --wf_mz 0.55 \
  --wf_int 1.1 \
  --LET_threshold 3 \
  --save_plots ../tests/test30.pdf

echo $'\n\n\n\n\ntest #31'
python plot_spectra.py \
  --query_data ../data/lcms_query_library.csv \
  --reference_data ../data/lcms_reference_library.csv \
  --spectrum_ID2 "Malyngamide J M+H" \
  --similarity_measure tsallis \
  --chromatography_platform HRMS \
  --high_quality_reference_library True \
  --wf_mz 0.55 \
  --wf_int 1.1 \
  --LET_threshold 3 \
  --save_plots ../tests/test31.pdf

echo $'\n\n\n\n\ntest #32'
python plot_spectra.py \
  --query_data ../data/gcms_query_library.csv \
  --reference_data ../data/gcms_reference_library.csv \
  --chromatography_platform NRMS \
  --similarity_measure cosine \
  --normalization_method standard \
  --save_plots ../tests/test32.pdf

echo $'\n\n\n\n\ntest #33'
python plot_spectra.py \
  --query_data ../data/gcms_query_library.csv \
  --reference_data ../data/gcms_reference_library.csv \
  --chromatography_platform NRMS \
  --similarity_measure tsallis \
  --wf_mz 2 \
  --wf_intensity 0.5 \
  --normalization_method standard \
  --entropy_dimension 2 \
  --save_plots ../tests/test33.pdf

echo $'\n\n\n\n\ntest #34'
python plot_spectra.py \
  --query_data ../data/gcms_query_library.csv \
  --reference_data ../data/gcms_reference_library.csv \
  --chromatography_platform NRMS \
  --similarity_measure tsallis \
  --wf_mz 2 \
  --wf_intensity 0.5 \
  --normalization_method standard \
  --entropy_dimension 0.5 \
  --save_plots ../tests/test34.pdf

echo $'\n\n\n\n\ntest #35'
python plot_spectra.py \
  --query_data ../data/gcms_query_library.csv \
  --reference_data ../data/gcms_reference_library.csv \
  --chromatography_platform NRMS \
  --similarity_measure tsallis \
  --normalization_method standard \
  --entropy_dimension 2 \
  --save_plots ../tests/test35.pdf

echo $'\n\n\n\n\ntest #36'
python plot_spectra.py \
  --query_data ../data/gcms_query_library.csv \
  --reference_data ../data/gcms_reference_library.csv \
  --chromatography_platform NRMS \
  --similarity_measure renyi \
  --normalization_method standard \
  --entropy_dimension 1.1 \
  --save_plots ../tests/test36.pdf

echo $'\n\n\n\n\ntest #37'
python plot_spectra.py \
  --query_data ../data/gcms_query_library.csv \
  --reference_data ../data/gcms_reference_library.csv \
  --chromatography_platform NRMS \
  --similarity_measure renyi \
  --normalization_method standard \
  --entropy_dimension 0.9 \
  --save_plots ../tests/test37.pdf

echo $'\n\n\n\n\ntest #38'
python plot_spectra.py \
  --query_data ../data/gcms_query_library.csv \
  --reference_data ../data/gcms_reference_library.csv \
  --chromatography_platform NRMS \
  --similarity_measure shannon \
  --normalization_method standard \
  --save_plots ../tests/test38.pdf

echo $'\n\n\n\n\ntest #39'
python plot_spectra.py \
  --query_data ../data/gcms_query_library.csv \
  --reference_data ../data/gcms_reference_library.csv \
  --chromatography_platform NRMS \
  --similarity_measure cosine \
  --save_plots ../tests/test39.pdf

echo $'\n\n\n\n\ntest #40'
python plot_spectra.py \
  --query_data ../data/gcms_query_library.csv \
  --reference_data ../data/gcms_reference_library.csv \
  --chromatography_platform NRMS \
  --spectrum_preprocessing_order FNLW\
  --similarity_measure cosine \
  --wf_mz 0.5 \
  --save_plots ../tests/test40.pdf

echo $'\n\n\n\n\ntest #41'
python plot_spectra.py \
  --query_data ../data/gcms_query_library.csv \
  --reference_data ../data/gcms_reference_library.csv \
  --chromatography_platform NRMS \
  --spectrum_preprocessing_order FNW\
  --similarity_measure cosine \
  --mz_min 10 \
  --mz_max 450 \
  --int_min 10 \
  --int_max 1500 \
  --save_plots ../tests/test41.pdf

echo $'\n\n\n\n\ntest #42'
python plot_spectra.py \
  --query_data ../data/gcms_query_library.csv \
  --reference_data ../data/gcms_reference_library.csv \
  --chromatography_platform NRMS \
  --spectrum_preprocessing_order FNLW\
  --similarity_measure cosine \
  --mz_max 100 \
  --save_plots ../tests/test42.pdf

echo $'\n\n\n\n\ntest #43'
python plot_spectra.py \
  --query_data ../data/gcms_query_library.csv \
  --reference_data ../data/gcms_reference_library.csv \
  --chromatography_platform NRMS \
  --spectrum_preprocessing_order FNLW \
  --similarity_measure cosine \
  --int_max 300 \
  --save_plots ../tests/test43.pdf


echo $'\n\n\n\n\ntest #44'
python plot_spectra.py \
  --query_data ../data/gcms_query_library.csv \
  --reference_data ../data/gcms_reference_library.csv \
  --chromatography_platform NRMS \
  --similarity_measure cosine \
  --int_min 100 \
  --save_plots ../tests/test44.pdf

echo $'\n\n\n\n\ntest #45'
python plot_spectra.py \
  --query_data ../data/gcms_query_library.csv \
  --reference_data ../data/gcms_reference_library.csv \
  --chromatography_platform NRMS \
  --spectrum_preprocessing_order NFLW \
  --similarity_measure cosine \
  --noise_threshold 0.1 \
  --save_plots ../tests/test45.pdf

echo $'\n\n\n\n\ntest #46'
python plot_spectra.py \
  --query_data ../data/gcms_query_library.csv \
  --reference_data ../data/gcms_reference_library.csv \
  --chromatography_platform NRMS \
  --spectrum_preprocessing_order NFLW \
  --similarity_measure cosine \
  --noise_threshold 0.4 \
  --save_plots ../tests/test46.pdf

echo $'\n\n\n\n\ntest #47'
python plot_spectra.py \
  --query_data ../data/gcms_query_library.csv \
  --reference_data ../data/gcms_reference_library.csv \
  --chromatography_platform NRMS \
  --spectrum_preprocessing_order FNLW \
  --similarity_measure cosine \
  --LET_threshold 2 \
  --save_plots ../tests/test47.pdf

echo $'\n\n\n\n\ntest #48'
python plot_spectra.py \
  --query_data ../data/gcms_query_library.csv \
  --reference_data ../data/gcms_reference_library.csv \
  --chromatography_platform NRMS \
  --spectrum_preprocessing_order LWF \
  --similarity_measure cosine \
  --LET_threshold 2 \
  --save_plots ../tests/test48.pdf

echo $'\n\n\n\n\ntest #49'
python plot_spectra.py \
  --query_data ../data/gcms_query_library.csv \
  --reference_data ../data/gcms_reference_library.csv \
  --chromatography_platform NRMS \
  --spectrum_preprocessing_order LWF \
  --similarity_measure cosine \
  --LET_threshold 1.5 \
  --save_plots ../tests/test49.pdf

echo $'\n\n\n\n\ntest #50'
python plot_spectra.py \
  --query_data ../data/gcms_query_library.csv \
  --reference_data ../data/gcms_reference_library.csv \
  --spectrum_ID1 ID_2 \
  --chromatography_platform NRMS \
  --spectrum_preprocessing_order LWF \
  --similarity_measure cosine \
  --LET_threshold 1.5 \
  --save_plots ../tests/test50.pdf

echo $'\n\n\n\n\ntest #51'
python plot_spectra.py \
  --query_data ../data/gcms_query_library.csv \
  --reference_data ../data/gcms_reference_library.csv \
  --chromatography_platform NRMS \
  --spectrum_ID2 463514 \
  --spectrum_preprocessing_order LWF \
  --similarity_measure cosine \
  --LET_threshold 1.5 \
  --save_plots ../tests/test51.pdf

echo $'\n\n\n\n\ntest #52'
python plot_spectra.py \
  --query_data ../data/gcms_query_library.csv \
  --reference_data ../data/gcms_reference_library.csv \
  --chromatography_platform NRMS \
  --y_axis_transformation none \
  --save_plots ../tests/test52.pdf

echo $'\n\n\n\n\ntest #53'
python plot_spectra.py \
  --query_data ../data/gcms_query_library.csv \
  --reference_data ../data/gcms_reference_library.csv \
  --chromatography_platform NRMS \
  --y_axis_transformation log10 \
  --save_plots ../tests/test53.pdf

echo $'\n\n\n\n\ntest #54'
python plot_spectra.py \
  --query_data ../data/gcms_query_library.csv \
  --reference_data ../data/gcms_reference_library.csv \
  --chromatography_platform NRMS \
  --y_axis_transformation sqrt \
  --save_plots ../tests/test54.pdf

echo $'\n\n\n\n\ntest #55'
python plot_spectra.py \
  --query_data ../data/gcms_query_library.csv \
  --reference_data ../data/gcms_reference_library.csv \
  --spectrum_ID1 ID_2 \
  --spectrum_ID2 463514 \
  --chromatography_platform NRMS \
  --spectrum_preprocessing_order LWF \
  --save_plots ../tests/test55.pdf

echo $'\n\n\n\n\ntest #56'
python plot_spectra.py \
  --query_data ../data/gcms_query_library.csv \
  --reference_data ../data/gcms_reference_library.csv \
  --spectrum_ID1 ID_1 \
  --spectrum_ID2 ID_2 \
  --chromatography_platform NRMS \
  --spectrum_preprocessing_order LWF \
  --save_plots ../tests/test56.pdf

echo $'\n\n\n\n\ntest #57'
python plot_spectra.py \
  --query_data ../data/gcms_query_library.csv \
  --reference_data ../data/gcms_reference_library.csv \
  --spectrum_ID1 463514 \
  --spectrum_ID2 112312 \
  --chromatography_platform NRMS \
  --spectrum_preprocessing_order LWF \
  --save_plots ../tests/test57.pdf

echo $'\n\n\n\n\ntest #58'
python plot_spectra.py \
  --query_data ../data/lcms_query_library.csv \
  --reference_data ../data/lcms_reference_library.csv \
  --spectrum_ID1 212 \
  --spectrum_ID2 100 \
  --chromatography_platform HRMS \
  --save_plots ../tests/test58.pdf

echo $'\n\n\n\n\ntest #59'
python plot_spectra.py \
  --query_data ../data/lcms_query_library.csv \
  --reference_data ../data/lcms_reference_library.csv \
  --spectrum_ID1 "Jamaicamide A M+H" \
  --spectrum_ID2 "Malyngamide J M+H" \
  --chromatography_platform HRMS \
  --save_plots ../tests/test59.pdf

echo -e '\n\nFinished Testing.\n'


