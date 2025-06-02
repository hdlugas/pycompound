#!/bin/bash

cd ${PWD}/../src


echo $'\n\n\ntest #1'
python build_library.py \
  --input_path ../data/GNPS-SELLECKCHEM-FDA-PART1.mgf \
  --output_path ../data/lcms_library_from_mgf.csv

echo $'\n\n\ntest #2'
python build_library.py \
  --input_path ../data/1min.mzML \
  --output_path ../data/lcms_library_from_mzML.csv

echo $'\n\n\ntest #3'
python build_library.py \
  --input_path ../data/C01_MTBSTFA.cdf \
  --output_path ../data/gcms_library_from_cdf.csv

echo $'\n\n\ntest #4'
python build_library.py \
  --output_path ../data/lcms_library_from_mgf2.csv

echo $'\n\n\ntest #5'
python build_library.py \
  --output_path ../data/lcms_library_from_mzML2.csv

echo $'\n\n\ntest #6'
python build_library.py \
  --output_path ../data/gcms_library_from_cdf.csv

echo $'\n\n\ntest #7'
python build_library.py \
  --input_path ../data/GNPS-SELLECKCHEM-FDA-PART1.mgf \

echo $'\n\n\ntest #8'
python build_library.py \
  --input_path ../data/1min.mzML \

echo $'\n\n\ntest #9'
python build_library.py \
  --input_path ../data/C01_MTBSTFA.cdf \

echo $'\n\nTesting complete.\n'



