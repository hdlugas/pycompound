@echo off

:: define paths
set QUERY_PATH1=%cd%\tests\data\gcms_query.txt
set QUERY_PATH2=%cd%\tests\data\gcms_query_tuning.txt
set REF_PATH=%cd%\tests\data\trimmed_gcms_reference_library.txt

echo "############ Starting Toy Examples ############"
echo ""
echo ">>> Plot Spectra"

:: plot spectra
python src\pycompound\plot_spectra_CLI.py ^
 --query_data %QUERY_PATH1% ^
 --reference_data %REF_PATH% ^
 --similarity_measure cosine ^
 --chromatography_platform NRMS ^
 --spectrum_ID1 ID_1 ^
 --spectrum_ID2 463-51-4 ^
 --output_path %cd%\CLI_plotting_example.pdf

echo ""
echo ">>> Run Spectral library Matching" 

:: run spectral library matching
python src\pycompound\spec_lib_matching_CLI.py ^
 --query_data %QUERY_PATH2% ^
 --reference_data %REF_PATH% ^
 --chromatography_platform NRMS ^
 --output_identification %cd%\CLI_identification_output_example.txt ^
 --output_similarity_scores %cd%\CLI_similarity_scores_output_example.txt

echo ""
echo ">>> Tune Parameters via Grid Search"

:: grid tuning
python src\pycompound\tuning_CLI_grid.py ^
 --query_data %QUERY_PATH2% ^
 --reference_data %REF_PATH% ^
 --chromatography_platform NRMS ^
 --wf_int 1,2 ^
 --wf_mz 0,2 ^
 --output_path %cd%\CLI_grid_tuning_output_example.txt

echo ""
echo ">>> Tune Parameters via DE Optimization

:: differential evolution tuning
python src\pycompound\tuning_CLI_DE.py ^
 --query_data %QUERY_PATH2% ^
 --reference_data %REF_PATH% ^
 --chromatography_platform NRMS ^
 --opt wf_mz wf_int ^
 --bound wf_mz=0.0:5.0 ^
 --bound wf_int=0.0:5.0 ^
 --maxiter 10 ^
 --workers 5

echo ""
echo "############ Completed Toy Examples ############"