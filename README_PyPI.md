<p align="center">
  <img src="images/PyCompound_logo.png" alt="PyCompound Logo" width="180"/>
</p>

PyCompound is a Python-based tool for spectral library matching designed to identify chemical compounds from mass spectrometry data. It is available in three formats: a Python package, a command-line interface (CLI), and a graphical user interface (GUI) built with Python/Shiny. PyCompound provides a flexible and extensible framework for spectral library matching and introduces several key features. These include entropy-based similarity measures such as Shannon, Tsallis, and the Rényi entropy similarity measure introduced here for the first time, as well as conventional similarity metrics, including cosine and binary similarity measures. PyCompound supports customizable preprocessing workflows that allow users to explicitly control the order of spectral preprocessing steps. In addition, PyCompound includes transformation parameter optimization using metaheuristic algorithms and supports the construction of user-defined mixture or composite similarity measures by combining two or more similarity metrics. PyCompound supports both high-resolution mass spectrometry (HRMS) data (e.g., LC-MS/MS) and nominal-resolution mass spectrometry (NRMS) data (e.g., GC-MS). For the full documentation, including toy examples, see the GitHub repository https://github.com/hdlugas/pycompound.

# Installation

## Install directly from PyPI (recommended):
```
conda create -n pycompound_env python=3.12 -y
conda activate pycompound_env
pip install pycompound==0.1.14
```

## Source
Install directly from GitHub:
```
conda create -n pycompound_env -y python=3.12
conda activate pycompound_env
pip install git+https://github.com/hdlugas/pycompound.git
```

## Shiny
The Shiny app is available at https://connect.posit.cloud/fy7392. Video tutorials of the PyCompound Shiny application are available on YouTube (https://www.youtube.com/@PyCompound).
