# PyCompound

A Python-based tool for spectral library matching, PyCompound is available as a Python package (pycompound) with a command-line interface (CLI) available and as a GUI application build with Python/Shiny. It performs spectral library matching to identify chemical compounds, offering a range of spectrum preprocessing transformations and similarity measures, including Cosine, three entropy-based similarity measures, and a plethora of binary similarity measures. PyCompound also includes functionality to tune parameters commonly used in a compound identification workflow given a query library of spectra with known ID. PyCompound supports both high-resolution mass spectrometry (HRMS) data (e.g., LC-MS/MS) and nominal-resolution mass spectrometry (NRMS) data (e.g., GC-MS). For the full documentation, see the GitHub repository https://github.com/hdlugas/pycompound.

# Installation

## Install directly from PyPI (recommended):
```
conda create -n pycompound_env python=3.12 -y
conda activate pycompound_env
pip install pycompound==0.1.11
```

## Source
Install directly from GitHub:
```
conda create -n pycompound_env -y python=3.12
conda activate pycompound_env
pip install git+https://github.com/hdlugas/pycompound.git
```

