# PyCompound
A Python-based tool for spectral library matching, PyCompound is available as a Python package with a command-line interface (CLI) available and as a GUI application build with Python/Shiny. It performs spectral library matching to identify chemical compounds, offering a range of spectrum preprocessing transformations and similarity measures, including Cosine and three entropy-based similarity measures. PyCompound supports both high-resolution mass spectrometry (HRMS) data (e.g., LC-MS/MS) and nominal-resolution mass spectrometry (NRMS) data (e.g., GC-MS).

## Table of Contents
- [1. Install dependencies](#create-conda-env)
- [2. Functionality](#functionality)
   - [2.1 Spectrum Preprocessing Transformations](#spec-preprocessing-transformations)
   - [2.2 Similarity Measures](#similarity-measures)
- [3. Usage](#usage)
   - [3.1 Parameter descriptions](#param_descriptions)
   - [3.2 Obtain LC-MS/MS or GC-MS library from MGF, mzML, or cdf file](#process-data)
   - [3.3 Run spectral library matching](#run-spec-lib-matching)
   - [3.4 Tune parameters](#tuning)
   - [3.5 Plot a query spectrum against a reference spectrum before and after spectrum preprocessing transformations](#plotting)
   - [3.6 Shiny application](#shiny)
- [4. Bugs/Questions?](#bugs-questions)

<a name="create-conda-env"></a>
## 1. Install dependencies
PyCompound requires the Python dependencies Matplotlib, NumPy, Pandas, SciPy, Pyteomics, and netCDF4. Specifically, this software was validated with python=3.12.4, matplotlib=3.8.4, numpy=1.26.4, pandas=2.2.2, scipy=1.13.1, pyteomics=4.7.2, netCDF4=1.6.5, lxml=5.1.0, joblib=1.5.2, and shiny=1.4.0, although it may work with other versions of these tools. A user may consider creating a conda environment (see [https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html) for guidance on getting started with conda if you are unfamiliar). For a system with conda installed, one can create the environment pycompound_env, activate it, and install the necessary dependencies with:
```
conda create -n pycompound_env python=3.12
conda activate pycompound_env
pip install pycompound==0.1.7
```

<a name="functionality"></a>
## 2. Functionality

<a name="spec-preprocessing-transformations"></a>
## 2.1 Spectrum Preprocessing Transformations
The following spectrum preprocessing transformations are offered:

-   Filtering: Given user-defined parameters (mz_min,mz_max),
    (int_min,int_max) and spectrum $I$ with m/z values
    $(m_{1},m_{2},...,m_{n})$ and intensities $(x_{1},x_{2},...,x_{n})$,
    the transformed spectrum $I^{\star}$ consists of the peaks
    $(m_{i},x_{i})$ in $I$ such that mz_min $\leq m_{i}\leq$ mz_max and
    int_min $\leq x_{i}\leq$ int_max.
    
-   Weight Factor Transformation: Given a pair of user-defined weight
    factor parameters $(\text{a,b})$ and spectrum $I$ with m/z values
    $(m_{1},m_{2},...,m_{n})$ and intensities $(x_{1},x_{2},...,x_{n})$,
    the transformed spectrum $I^{\star}$ has the same m/z values as $I$
    and has intensities given by
    $I^{\star}:=(m_{1}^{\text{a}}\cdot x_{1}^{\text{b}},m_{2}^{\text{a}}\cdot x_{2}^{\text{b}},...,m_{n}^{\text{a}}\cdot x_{n}^{\text{b}})$.

-   Low-Entropy Transformation: Given a user-defined low-entropy
    threshold parameter $T$ and spectrum $I$ with intensities
    $(x_{1},x_{2},...,x_{n})$, $\sum_{i=1}^nx_i = 1$, and Shannon
    entropy $H_{Shannon}(I)=-\sum_{i=1}^{n}x_{i}\cdot ln(x_{i})$, the
    transformed spectrum intensities
    $I^{\star}=(x_{1}^{\star},x_{2}^{\star},...,x_{n}^{\star})$ are such
    that, for all $i\in\{1,2,...,n\}$, $x_{i}^{\star}=x_{i}$ if
    $H_{Shannon}(I)\geq T$ and
    $x_{i}^{\star}=x_{i}^{\frac{1+H_{Shannon}(I)}{1+T}}$ if
    $H_{Shannon}(I)<T$.

-   Centroiding (only applicable to HRMS data): Given a user-defined
    window-size parameter $w_{centroiding}$ and a spectrum $I$ with m/z
    values $(m_{1},m_{2},...,m_{n})$ and intensities
    $(x_{1},x_{2},...,x_{n})$, the transformed spectrum $I^{\star}$
    merges adjacent peaks $(m_{i},x_{i}),(m_{i+1},x_{i+1})$ into the
    peak
    $(\frac{m_{i}\cdot x_{i}+m_{i+1}\cdot x_{i+1}}{x_{i}+x_{i+1}},x_{i}+x_{i+1})$
    if $|m_{i}-m_{i+1}|< w_{centroiding}$ for
    $i\in\{1,2,...,n-1\}$. This centroiding procedure generalizes to
    more than two peaks whose m/z values are within a distance
    $w_{centroiding}$ of each other.

-   Noise Removal: Given a user-defined noise removal parameter $r$ and
    a spectrum $I$ with intensities $(x_{1},x_{2},...,x_{n})$, noise
    removal removes peaks from $I$ with
    $x_{j}< r\cdot\text{max}(\{x_{1},x_{2},...,x_{n}\})$ for
    $j\in\{1,2,...,n\}$.

-   Matching (only applicable to HRMS data): Given a user-defined
    window-size parameter $w_{matching}$ and two spectra $I$, $J$ with
    m/z ratios $(a_{1},a_{2},...,a_{n}), (b_{1},b_{2},...,b_{m})$ and
    intensities $(x_{1},x_{2},...,x_{n}), (y_{1},y_{2},...,y_{m})$,
    respectively, of which we would like to measure the similarity
    between, the matching procedure outputs two spectra
    $I^{\star},J^{\star}$ containing the same number of peaks with
    $I^{\star}$ and $J^{\star}$ having intensities and
    identical m/z ratios. Specifically, for a given peak $(a_{i},x_{i})$
    of $I$, if there are no peaks $(b_{j},y_{j})$ in $J$ with
    $|a_{i}-b_{j}|< w_{matching}$, then the peak $(a_{i},x_{i})$
    remains in $I^{\star}$ and the peak $(a_{i},0)$ is included in
    $J^{\star}$. If there is at least one peak $(b_{j},y_{j})$ with
    $|a_{i}-b_{j}|< w_{matching}$, then the peak $(a_{i},x_{i})$
    remains in $I^{\star}$ and the peak
    $(a_{i},\sum_{j\text{ such that }|a_{i}-b_{j}|< w_{matching}}b_{j})$
    is included in $J^{\star}$. This procedure is applied when
    transposing the roles of $I$ and $J$ as well.

<a name="similarity-measures"></a>
## 2.2 Similarity Measures
Given a pair of processed spectra intensities
$I=(a_{1},a_{2},...,a_{n}), J=(b_{1},b_{2},...,b_{n})\in\mathbb{R}^{n}$
with $0\leq a_{i},b_{i}\leq 1$ for all $i\in\{1,2,...,n\}$ and
$\sum_{i=1}^{n}a_{i}=\sum_{i=1}^{n}b_{i}=1$, PyCompound provides
functionality for computing the following similarity measures:

-   Cosine Similarity Measure:

```math
S_{Cosine}(I,J)=\frac{I\circ J}{|I|_{2}\cdot |J|_{2}}
```
where multiplication in the numerator refers to the dot product $I\circ J=a_{1}b_{1}+a_{2}b_{2}+...+a_{n}b_{n}$ of $I$ and $J$ and multiplication in the denominator refers to multiplication of the $L^{2}$-norms of $I$ and $J$, $\vert I\vert_{2}=\sqrt{a_{1}^{2}+a_{2}^{2}+...+a_{n}^{2}}, \vert J\vert_{2}=\sqrt{b_{1}^{2}+b_{2}^{2}+...+b_{n}^{2}}$.

-   Shannon Entropy Similarity Measure:

```math
S_{Shannon}(I,J) = 1-\frac{2\cdot H_{Shannon}\left(\frac{I+J}{2}\right) - H_{Shannon}(I)-H_{Shannon}(J)}{ln(4)},
```

```math
H_{Shannon}(I)=-\sum_{i=1}^{n}a_{i}\cdot ln(a_{i})
```

-    Tsallis Entropy Similarity Measure:

```math
S_{Tsallis}(I,J,q)=1-\frac{2\times H_{Tsallis}(I/2+J/2,q)-H_{Tsallis}(I,q)-H_{Tsallis}(J,q)}{N_{Tsallis}(I,J,q)},
```

```math
N_{Tsallis}(I,J,q):=\frac{\sum_{i=1}^{n}\left(2\left(\frac{a_{i}}{2}\right)^{q}+2\left(\frac{b_{i}}{2}\right)^{q}-a_{i}^{q}-b_{i}^{q}\right)}{1-q},
```

```math
H_{Tsallis}(I,q)=\frac{\left(\sum_{i=1}^{n}a_{i}^{q}\right)-1}{1-q},
```

```math
q\neq 1, \ q>0
```

-   Rényi Entropy Similarity Measure:

```math
S_{Renyi}(I,J,q)=1-\frac{2\times H_{Renyi}(I/2+J/2,q)-H_{Renyi}(I,q)-H_{Renyi}(J,q)}{N_{Renyi}(I,J,q)},
```

```math
N_{Renyi}(I,J,q):=\left(\frac{1}{1-q}\right)\left(2\times ln\left(\sum_{i}(a_{i}/2)^{q}+\sum_{j}(b_{j}/2)^{q}\right)-ln(\sum_{i}a_{i}^{q})-ln(\sum_{i}b_{i}^{q})\right),
```

```math
H_{Renyi}(I,q)=\frac{1}{1-q}ln(\sum_{i=1}^{n}a_{i}^{q}),
```

```math
q\neq 1, \ q>0
```

Additionally, the plethora of binary similarity measures considered in https://doi.org/10.3390/metabo12080694 are available along with a mixture similarity measure that is a weighted sum of the four non-binary similarity measures (i.e. Cosine, Shannon Entropy, Renyi, and Tsallis).

<a name="usage"></a>
## 3. Usage
PyCompound has three main capabilities:
1. Plotting a query spectrum vs. a reference spectrum before and after preprocessing transformations.
2. Running spectral library matching to identify compounds based on their mass spectrometry data
3. Tuning parameters to maximize accuracy given a query dataset with known compuond IDs (e.g. from targeted metabolomics experiments).

These tasks are implemented separately for the cases of (i) NRMS and (ii) HRMS data due to the different spectrum preprocessing transformations stemming from a different format in the mass to charge (m/z) ratios in NRMS vs HRMS data. Example scripts which implement these tasks can be found in the pycompound/tests directory.

<a name="param_descriptions"></a>
### 3.1 Parameter descriptions

For the function build_library_from_raw_data:
```
--input_path: Path to input file (must be either mgf, mzMZ, msp, cdf, or json file). Mandatory argument.

--output_path: Path to output text file. Default: current working directory.

--is_reference: Boolean flag indicating whether IDs of spectra should be written to output. Only pass True if building a library with known compound IDs. Only applicable to MGF files. Options: \'True\', \'False\'. Optional argument. Default: False.
```

Common parameters:
```
--query_data (mandatory argument):
  * HRMS case: mgf, mzML, msp, json, or txt file of query mass spectrum/spectra to be identified. If txt file, must have at least 3 columns with each row corresponding to a single ion fragment of a mass spectrum, one 'id' column containing an identifier, one 'mz_ratio' column corresponding to the mass to charge (m/z) ratios, and one 'intensity' column containing the intensities. For example, if spectrum A has 3 ion fragments, then there would be three rows in this text file corresponding to spectrum A. Optional columns for the text file are 'precursor_ion_mz', 'ionization_mode', and 'adduct'.
  * NRMS case: cdf or txt file of query mass spectrum/spectra to be identified. If txt file, same format as in HRMS case is required.

--reference_data (mandatory argument): Same format text file as query_data except of reference library spectra. We recommend using the reference libraries from our Zenodo database ([https://zenodo.org/records/12786324](https://zenodo.org/records/12786324); stored on Zenodo due to file size limitations on GitHub).

--precursor_ion_mz_tolerance (only applicable to HRMS): positive float representing a window size around each query spectrum's precursor ion mass:charge ratio in which candidate reference spectra must lie to be considered in compound identification. Default: None. 

--ionization_mode (only applicable to HRMS): Positive, Negative, or None. Default: None.

--adduct (only applicable to HRMS): Options: H, NH3, NH4, Na, K, N/A. Default: N/A.

--likely_reference_IDs: text file with one column containing the IDs of a subset of all compounds in the reference_data to be used in spectral library matching. Each ID in this file must be an ID in the reference library. Default: None (i.e. default is to use entire reference library)

--similarity_measure: cosine, shannon, renyi, tsallis, mixture, jaccard, dice, 3w_jaccard, sokal_sneath, binary_cosine, mountford, mcconnaughey, driver_kroeber, simpson, braun_banquet, fager_mcgowan, kulczynski, intersection, hamming, hellinger. Default: cosine.

--weights: dict of weights to give to each non-binary similarity measure (i.e. cosine, shannon, renyi, and tsallis) when the mixture similarity measure is specified. Default: 0.25 for each of the four non-binary similarity measures.

--spectrum_preprocessing_order: The spectrum preprocessing transformations and the order in which they are to be applied. Note that these transformations are applied prior to computing similarity scores. Format must be a string with 2-6 characters chosen from C, F, M, N, L, W representing centroiding, filtering based on mass/charge and intensity values, matching, noise removal, low-entropy trannsformation, and weight-factor-transformation, respectively. For example, if \'WCM\' is passed, then each spectrum will undergo a weight factor transformation, then centroiding, and then matching. Note that if an argument is passed, then \'M\' must be contained in the argument, since matching is a required preprocessing step in spectral library matching of HRMS data. Furthermore, \'C\' must be performed before matching since centroiding can change the number of ion fragments in a given spectrum. Note that C and M are not applicable to NRMS data. Default: FCNMWL for HRMS and FNLW for NRMS.')

--high_quality_reference_library: True/False flag indicating whether the reference library is considered to be of high quality. If True, then the spectrum preprocessing transformations of filtering and noise removal are performed only on the query spectrum/spectra. If False, all spectrum preprocessing transformations specified will be applied to both the query and reference spectra. Default: False')

--mz_min: Remove all peaks with mass/charge value less than mz_min in each spectrum. Default: 0

--mz_max: Remove all peaks with mass/charge value greater than mz_max in each spectrum. Default: 9999999

--int_min: Remove all peaks with intensity value less than int_min in each spectrum. Default: 0

--int_max: Remove all peaks with intensity value greater than int_max in each spectrum. Default: 9999999

--window_size_centroiding (only for HRMS): Window size parameter used in centroiding a given spectrum. Default: 0.5

--window_size_matching (only for HRMS): Window size parameter used in matching a query spectrum and a reference library spectrum. Default: 0.5

--noise_threshold: Ion fragments (i.e. points in a given mass spectrum) with intensity less than max(intensities)*noise_threshold are removed. Default: 0.0

--wf_mz: Mass/charge weight factor parameter. Default: 0.0

--wf_intensity: Intensity weight factor parameter. Default: 0.0

--LET_threshold: Low-entropy transformation threshold parameter. Spectra with Shannon entropy less than LET_threshold are transformed according to intensitiesNew=intensitiesOriginal^{(1+S)/(1+LET_threshold)}. Default: 0.0

--entropy_dimension: Entropy dimension parameter. Must have positive value other than 1. When the entropy dimension is 1, then Renyi and Tsallis entropy are equivalent to Shannon entropy. Therefore, this parameter only applies to the renyi and tsallis similarity measures. This parameter will be ignored if similarity measure cosine or shannon is chosen. Default: 1.1
```

Parameters specific to run_spec_lib_matching_on_HRMS_data and run_spec_lib_matching_on_NRMS_data:

```
--n_top_matches_to_save: The number of top matches to report. For example, if n_top_matches_to_save=5, then for each query spectrum, the five reference spectra with the largest similarity with the given query spectrum will be reported. Default: 1

--print_id_results: Flag that prints identification results if True. Default: False

--output_identification: Output text file containing the most-similar reference spectra for each query spectrum along with the corresponding similarity scores. Default is to save identification output in current working directory with filename 'output_identification.txt'.

--output_similarity_scores: Output text file containing similarity scores between all query spectrum/spectra and all reference spectra. Each row corresponds to a query spectrum, the left-most column contains the query spectrum/spectra identifier, and the remaining column contain the similarity scores with respect to all reference library spectra. If no argument passed, then this text file is written to the current working directory with filename output_all_similarity_scores.txt.
```

Parameters specific to tune_params_on_HRMS_data_grid and tune_params_on_NRMS_data_grid:
```
`` grid: dict object such as {'similarity_measure':['cosine','shannon'], 'spectrum_preprocessing_order':['FCNMWL'], 'mz_min':[0], 'mz_max':[9999999], 'int_min':[0], 'int_max':[99999999], 'window_size_centroiding':[0.5], 'window_size_matching':[0.5], 'noise_threshold':[0.0,0.1], 'wf_mz':[0.0], 'wf_int':[1.0], 'LET_threshold':[0.0], 'entropy_dimension':[1.1], 'high_quality_reference_library':[False]} with all possible combinations of parameters being utilized.

--output_path: path to output text file containing the accuracies for each possible combination of parameters. If no argument is passed, then the plots will be saved to ./tuning_param_output.txt in the current working directory.
```

Parameters specific to tune_params_DE:
```
-- optimize_params: list of continuous parameters (i.e. window_size_centroiding, window_size_matching, noise_threshold, wf_mz, wf_int, LET_threshold; window_size parameters only applicable to HRMS data) to optimize via differential evolution.

-- param_bounds: dict with keys being the parameters to optimize and values being a tuple of length 2 of the lower and upper bounds of acceptable parameter values. 

-- maxiters: maximum number of iterations of differential evolution.

-- de_workers: number of CPUs to utilize.
```

Parameters specific to generate_plots_on_HRMS_data and generate_plots_on_NRMS_data:
```
--spectrum_ID1: ID of one spectrum to be plotted. Default is first spectrum in the query library. Optional argument.

--spectrum_ID2: ID of another spectrum to be plotted. Default is first spectrum in the reference library. Optional argument.

--y_axis_transformation: transformation to apply to y-axis (i.e. intensity axis) of plots. Options: 'normalized', 'none', 'log10', and 'sqrt'. Default: 'normalized.')

--output_path: path to output PDF file containing the plots of the spectra before and after preprocessing transformations. If no argument is passed, then the plots will be saved to the PDF ./spectrum1_{spectrum_ID1}_spectrum2_{spectrum_ID2}_plot.pdf in the current working directory.
```


<a name="process-data"></a>
### 3.2 Obtain LC-MS/MS or GC-MS library from MGF, mzML, cdf, msp, or json file
To obtain a text file of LC-MS/MS spectra in the format necessary for spectral library matching from raw data in the form of an mgf, mzML, msp, json, or cdf file inside Python, one can run:
```
from pycompound.build_library import build_library_from_raw_data

build_library_from_raw_data(input_path='path_to_input_file', output_path='path_to_output_file', is_reference=False)
```

Since the other functionality provided by pycompound is capable of being directly run on mgf, mzML, msp, json, and cdf files, you may not need to directly build a library yourself. Some example mgf and json files one can use to build an LC-MS/MS library can be found from the Global Natural Products Social Molecular Networking (GNPS) databases here: [https://external.gnps2.org/gnpslibrary](https://external.gnps2.org/gnpslibrary). Some example mzML files one can use to build an LC-MS/MS library can be found in this repository: [https://github.com/HUPO-PSI/mzML](https://github.com/HUPO-PSI/mzML). Some example MSP files can be found here: [https://mona.fiehnlab.ucdavis.edu/downloads](https://mona.fiehnlab.ucdavis.edu/downloads). The mgf, mzML, msp, and json files provided in this repository are trimmed versions of files found in these referenced repositories. The script tests/test_build_libraries.py demonstrates this usage.

Full LC-MS/MS and GC-MS reference libraries are available at the Zenodo database ([https://zenodo.org/records/12786324](https://zenodo.org/records/12786324)). 

<a name="run-spec-lib-matching"></a>
### 3.3 Run spectral library matching
The files tests/test_spec_lib_matching.py, tests/test_spec_lib_matching_CLI, and tests/example_code_for_python_use.py demonstrate how some of the spectrum preprocessing functionality and similarity measures can be implemented either directly in Python or in the CLI wrapper. The two main functions - one for HRMS data and one for NRMS data - can be implemented as shown below inside Python:
```
from pycompound.spec_lib_matching import run_spec_lib_matching_on_HRMS_data
from pycompound.spec_lib_matching import run_spec_lib_matching_on_NRMS_data

run_spec_lib_matching_on_HRMS_data(
        query_data='path_to_query_library',
        reference_data='path_to_reference_library',
        likely_reference_IDs=None,
        similarity_measure='cosine',
        spectrum_preprocessing_order='FCNMWL',
        high_quality_reference_library=False,
        mz_min=0,
        mz_max=9999999,
        int_min=0,
        int_max=9999999,
        window_size_centroiding=0.5,
        window_size_matching=0.5,
        noise_threshold=0.0,
        wf_mz=0.0,
        wf_intensity=1.0,
        LET_threshold=0.0,
        entropy_dimension=1.1,
        n_top_matches_to_save=1,
        print_id_results=False,
        output_identification=None,
        output_similarity_scores=None)

run_spec_lib_matching_on_NRMS_data(
        query_data='path_to_query_library',
        reference_data='path_to_reference_library',
        likely_reference_IDs=None,
        similarity_measure='cosine',
        spectrum_preprocessing_order='FNLW',
        high_quality_reference_library=False,
        mz_min=0,
        mz_max=9999999,
        int_min=0,
        int_max=9999999,
        noise_threshold=0.0,
        wf_mz=0.0,
        wf_intensity=1.0,
        LET_threshold=0.0,
        entropy_dimension=1.1,
        n_top_matches_to_save=1,
        print_id_results=False,
        output_identification=None,
        output_similarity_scores=None)
```

To use the CLI version, one can run the following from the terminal:
```
python spec_lib_matching_CLI.py \
        --query_data ${PWD}/../tests/data/lcms_query_library.txt \
        --reference_data ${PWD}/../tests/data/full_GNPS_reference_library.txt \
        --chromatography_platform HRMS \
        --likely_reference_IDs None \
        --similarity_measure cosine \
        --spectrum_preprocessing_order FCNMWL \
        --high_quality_reference_library False \
        --mz_min 0 \
        --mz_max 9999999 \
        --int_min 0 \
        --int_max 9999999 \
        --window_size_centroiding 0.5 \
        --window_size_matching 0.5 \
        --noise_threshold 0.0 \
        --wf_mz 0.0 \
        --wf_intensity 1.0 \
        --LET_threshold 0.0 \
        --entropy_dimension 1.1 \
        --n_top_matches_to_save 1 \
        --print_id_results False \
        --output_identification ${PWD}/../tests/output_identification_HRMS.txt \
        --output_similarity_scores ${PWD}/../tests/output_similarity_scores_HRMS.txt

python spec_lib_matching_CLI.py \
        --query_data ${PWD}/../tests/data/lcms_query_library.txt \
        --reference_data ${PWD}/../tests/data/full_GNPS_reference_library.txt \
        --chromatography_platform NRMS \
        --likely_reference_IDs None \
        --similarity_measure cosine \
        --spectrum_preprocessing_order FCNMWL \
        --high_quality_reference_library False \
        --mz_min 0 \
        --mz_max 9999999 \
        --int_min 0 \
        --int_max 9999999 \
        --noise_threshold 0.0 \
        --wf_mz 0.0 \
        --wf_intensity 1.0 \
        --LET_threshold 0.0 \
        --entropy_dimension 1.1 \
        --n_top_matches_to_save 1 \
        --print_id_results False \
        --output_identification ${PWD}/../tests/output_identification_NRMS.txt \
        --output_similarity_scores ${PWD}/../tests/output_similarity_scores_NRMS.txt
```

For a user who may wish to incorporate our transformations and similarity measures directly in their python code similar to the example script tests/example_code_for_python_use.py, the available transformations and similarity measures are:
```
# Weight factor transformation
wf_transform(spec_mzs, spec_ints, wf_mz, wf_int)
"""
Perform weight factor transformation on a spectrum
Args:
   spec_mzs: 1d numpy array representing mass/charge values 
   spec_ints: 1d numpy array representing intensity values 
   wf_mz: float
   wf_int: float
Returns:
   np.ndarray: 1d numpy array of weight-factor-transformed spectrum intensities
"""

# Low-entropy transformation
LE_transform(intensity, thresh, normalization_method)
"""
Transforms spectrum's intensities if the Shannon entropy of the intensities is below some threshold
Args:
   intensity: 1d numpy array
   thresh: nonnegative float
   normalization_method: either 'standard' or 'softmax'
Returns:
   np.ndarray: 1d numpy array of transformed intensities
"""

# Filter HR-MS such as LC-MS/MS spectrum
filter_spec_lcms(spec, mz_min, mz_max, int_min, int_max, is_matched)
"""
Filter an MS/MS spectrum based on m/z and intensity values
Args:
   spec: N x 2 numpy array with first column being m/z and second column being intensity
   mz_min: minimum m/z value
   mz_max: maximum m/z value
   int_min: minimum intensity value
   int_max: maximum intensity value
   is_matched: flag to indicate whether the given spectrum has already been matched to another spectrum
Returns:
   np.ndarray: N x 2 numpy array with intensity of 0 put anywhere outside of the m/z and/or intensity bounds
"""

# Filter NR-MS such as GC-MS spectrum
filter_spec_gcms(spec, mz_min, mz_max, int_min, int_max)
"""
Filter an MS spectrum based on m/z and intensity values
Args:
   spec: N x 2 numpy array with first column being m/z and second column being intensity
   mz_min: minimum m/z value
   mz_max: maximum m/z value
   int_min: minimum intensity value
   int_max: maximum intensity value
Returns:
   np.ndarray: N x 2 numpy array with intensity of 0 put anywhere outside of the m/z and/or intensity bounds
"""

# Remove low-intensity noise
remove_noise(spec, nr)
"""
Remove low-intensity ion fragments
Args:
   spec: N x 2 numpy array with first column being m/z and second column being intensity
   nr: noise removal parameter; ion fragments with intensity less than max(intensity)*nr have intensity set to 0
Returns:
   np.ndarray: N x 2 numpy array
"""

# Centroid spectrum by merging close m/z peaks
centroid_spectrum(spec, window_size)
"""
Centroid a spectrum by merging ion fragments that are 'close' with respect to m/z value
Args:
   spec: N x 2 numpy array with the first column being mass/charge and the second column being intensity
   window_size: window-size parameter
Returns:
   np.ndarray: M x 2 numpy array with M <= N due to peaks being merged
"""

# Match peaks between two spectra
match_peaks_in_spectra(spec_a, spec_b, window_size)
"""
Align two spectra so that we obtain a list of intensity values from each spectrum of the same length
Args:
   spec_a: N x 2 numpy array with the first column being mass/charge and the second column being intensity
   spec_b: M x 2 numpy array with the first column being mass/charge and the second column being intensity
   window_size: window-size parameter
Returns:
   np.ndarray: K x 3 numpy array with first column being mass/charge, second column being matched intensities of spec_a, and third column being matched intensities of spec_b
"""

# Assign 0 to the intensities without m/z values
convert_spec(spec, mzs)
"""
Set intensity values to 0 where m/z values are missing
Args:
   spec: N x 2 dimensional numpy array
   mzs: length M list of entire span of mass/charge values considering both the query and reference libraries
Returns:
   np.ndarray: M x 2 dimensional numpy array
"""

# Cosine similarity
S_cos(ints_a, ints_b)
"""
Cosine similarity measure
Args:
   ints_a: 1d numpy array of intensities of a spectrum
   ints_b: 1d numpy array of intensities of a spectrum
Returns:
   float: float between 0 and 1 indicating the similarity of the two spectra
"""

# Shnnon entropy similarity
S_shannon(ints_a, ints_b)
"""
Shannon entropy similarity measure
Args:
   ints_a: 1d numpy array of intensities of a spectrum
   ints_b: 1d numpy array of intensities of a spectrum
Returns:
   float: float between 0 and 1 indicating the similarity of the two spectra
"""

# Renyi entropy similarity
S_renyi(ints_a, ints_b, q)
"""
Renyi entropy similarity measure
Args:
   ints_a: 1d numpy array of intensities of a spectrum
   ints_b: 1d numpy array of intensities of a spectrum
   q: positive float representing 'entropy dimension'
Returns:
   float: float between 0 and 1 indicating the similarity of the two spectra
"""

# Tsallis entropy similarity
S_tsallis(ints_a, ints_b, q)
"""
Tsallis entropy similarity measure
Args:
   ints_a: 1d numpy array of intensities of a spectrum
   ints_b: 1d numpy array of intensities of a spectrum
   q: positive float representing 'entropy dimension'
Returns:
   float: float between 0 and 1 indicating the similarity of the two spectra
"""
```


<a name="tuning"></a>
### 3.4 Tune parameters
Note that in order to tune parameters such as noise_threshold, LET_threshold etc., one must have a query library with compounds whose ground truth ID is known (e.g. from targeted metabolomics experiments). PyCompound offers two different methods of tuning parameters: one being an exhaustive grid search of pre-specified values, and the other being an optimization approach using differential evolution to optimize continuous parameters with respect to accuracy. The usage of the functions to tune parameters within Python is:
```
from pycompound.spec_lib_matching import tune_params_on_HRMS_data_grid
from pycompound.spec_lib_matching import tune_params_on_NRMS_data_grid
from pycompound.spec_lib_matching import tune_params_DE
from pathlib import Path

tune_params_on_HRMS_data_grid(
    query_data=f'{Path.cwd()}/tests/data/lcms_query_library_tuning.txt',
    reference_data=f'{Path.cwd()}/tests/data/full_GNPS_reference_library.txt',
    precursor_ion_mz_tolerance=0.5,
    ionization_mode='Positive',
    adduct='H',
    grid={'similarity_measure':['cosine'], 'spectrum_preprocessing_order':['FCNMWL'], 'mz_min':[0], 'mz_max':[9999999], 'int_min':[0], 'int_max':[99999999], 'window_size_centroiding':[0.5], 'window_size_matching':[0.1,0.5], 'noise_threshold':[0.0], 'wf_mz':[0.0], 'wf_int':[1.0], 'LET_threshold':[0.0], 'entropy_dimension':[1.1], 'high_quality_reference_library':[False]},
    output_path=f'{Path.cwd()}/tuning_param_output_HRMS.txt'
)

tune_params_on_NRMS_data_grid(
    query_data=f'{Path.cwd()}/tests/data/gcms_query_library_tuning.txt',
    reference_data=f'{Path.cwd()}/tests/data/gcms_reference_library.txt',
    grid={'similarity_measure':['cosine','shannon'], 'spectrum_preprocessing_order':['FNLW'], 'mz_min':[0], 'mz_max':[9999999], 'int_min':[0], 'int_max':[99999999], 'noise_threshold':[0.0,0.1], 'wf_mz':[0.0], 'wf_int':[1.0], 'LET_threshold':[0.0,3.0], 'entropy_dimension':[1.1], 'high_quality_reference_library':[False]},
    output_path=f'{Path.cwd()}/tuning_param_output_NRMS.txt'
)

tune_params_DE(
    query_data=f'{Path.cwd()}/tests/data/lcms_query_library_tuning.txt',
    reference_data=f'{Path.cwd()}/tests/data/full_GNPS_reference_library.txt',
    precursor_ion_mz_tolerance=0.1,
    ionization_mode='Positive',
    adduct='H',
    chromatography_platform='HRMS',
    similarity_measure='shannon',
    optimize_params=["wf_mz","wf_int"],
    param_bounds={"wf_mz":(0.0,5.0),"wf_int":(0.0,5.0)},
    default_params={"window_size_centroiding": 0.5, "window_size_matching":0.5, "noise_threshold":0.10, "wf_mz":0.0, "wf_int":1.0, "LET_threshold":0.0, "entropy_dimension":1.1},
    maxiters=2,
    de_workers=-1
)

tune_params_DE(
    query_data=f'{Path.cwd()}/tests/data/gcms_query_library_tuning.txt',
    reference_data=f'{Path.cwd()}/tests/data/gcms_reference_library.txt',
    chromatography_platform='NRMS',
    similarity_measure='renyi',
    optimize_params=["wf_mz","wf_int","LET_threshold","entropy_dimension"],
    param_bounds={"wf_mz":(0.0,5.0),"wf_int":(0.0,5.0),"LET_threshold":(0,5),"entropy_dimension":(1.01,3)},
    default_params={"noise_threshold":0.10, "wf_mz":0.0, "wf_int":1.0, "LET_threshold":0.0, "entropy_dimension":1.1},
    de_workers=-1
)
```

The CLI version can be run with:
```
python ../src/tuning_CLI_grid.py \
  --query_data ${PWD}/tests/data/lcms_query_library_tuning.txt \
  --reference_data ${PWD}/tests/data/full_GNPS_reference_library.txt \
  --precursor_ion_mz_tolerance 0.1 \
  --ionization_mode Positive \
  --adduct H \
  --chromatography_platform HRMS \
  --similarity_measure cosine \
  --spectrum_preprocessing_order FCNMWL \
  --high_quality_reference_library False \
  --mz_min 0 \
  --mz_max 9999999 \
  --int_min 0 \
  --int_max 9999999 \
  --window_size_centroiding 0.5 \
  --window_size_matching 0.1,0.5 \
  --noise_threshold 0.0 \
  --wf_mz 2,3 \
  --wf_intensity 1.0 \
  --LET_threshold 0.0 \
  --entropy_dimension 1.1 \
  --output_path ${PWD}/output_tuning_HRMS_grid.txt \

python ../src/pycompound/tuning_CLI_grid.py \
  --query_data ${PWD}/tests/data/gcms_query_library_tuning.txt \
  --reference_data ${PWD}/tests/data/gcms_reference_library.txt \
  --chromatography_platform NRMS \
  --similarity_measure cosine,shannon \
  --spectrum_preprocessing_order FCNMWL \
  --high_quality_reference_library False \
  --mz_min 0 \
  --mz_max 9999999 \
  --int_min 0 \
  --int_max 9999999 \
  --noise_threshold 0.0,0.1 \
  --wf_mz 0 \
  --wf_intensity 1.0 \
  --LET_threshold 0.0 \
  --entropy_dimension 1.1 \
  --output_path ${PWD}/output_tuning_NRMS_grid.txt \

python ../src/pycompound/tuning_CLI_DE.py \
  --chromatography_platform HRMS \
  --query_data ${PWD}/data/lcms_query_library_tuning.txt \
  --reference_data ${PWD}/data/full_GNPS_reference_library.txt \
  --precursor_ion_mz_tolerance 0.1 \
  --ionization_mode Positive \
  --adduct H \
  --similarity_measure cosine \
  --opt window_size_centroiding noise_threshold wf_mz \
  --bound window_size_centroiding=0.0:0.4 \
  --bound noise_threshold=0.0:0.20 \
  --bound wf_mz=0.0:5.0 \
  --maxiter 3 \
  --seed 1 \
  --workers 5

python ../src/pycompound/tuning_CLI_DE.py \
  --query_data ${PWD}/tests/data/gcms_query_library_tuning.txt \
  --reference_data ${PWD}/tests/data/gcms_reference_library.txt \
  --chromatography_platform NRMS \
  --similarity_measure cosine \
  --opt noise_threshold wf_mz \
  --bound noise_threshold=0.0:0.20 \
  --bound wf_mz=0.0:5.0 \
  --maxiter 3 \
  --seed 1 \
  --workers 4

```


<a name="plotting"></a>
### 3.5 Plot a query spectrum against a reference spectrum before and after spectrum preprocessing transformations
These functions plot a query spectrum against a reference spectrum, both before and after preprocessing. They support HRMS and NRMS data and can be used directly within Python with usage:
```
from pycompound.plot_spectra import generate_plots_on_HRMS_data
from pycompound.plot_spectra import generate_plots_on_NRMS_data

generate_plots_on_HRMS_data(
        query_data='path_to_query_library',
        reference_data='path_to_reference_data',
        spectrum_ID1=None,
        spectrum_ID2=None,
        similarity_measure='cosine',
        spectrum_preprocessing_order='FCNMWL',
        high_quality_reference_library=False,
        mz_min=0,
        mz_max=9999999,
        int_min=0,
        int_max=9999999,
        window_size_centroiding=0.5,
        window_size_matching=0.5,
        noise_threshold=0.0,
        wf_mz=0.0,
        wf_intensity=1.0,
        LET_threshold=0.0,
        entropy_dimension=1.1,
        y_axis_transformation='normalized',
        output_path=None
)

generate_plots_on_NRMS_data(
        query_data='path_to_query_library',
        reference_data='path_to_reference_data',
        spectrum_ID1=None,
        spectrum_ID2=None,
        similarity_measure='cosine',
        spectrum_preprocessing_order='FNLW',
        high_quality_reference_library=False,
        mz_min=0,
        mz_max=9999999,
        int_min=0,
        int_max=9999999,
        noise_threshold=0.0,
        wf_mz=0.0,
        wf_intensity=1.0,
        LET_threshold=0.0,
        entropy_dimension=1.1,
        y_axis_transformation='normalized',
        output_path=None
)
```

To use the command line version, one can run the following from the terminal:
```
python plot_spectra_CLI.py \
  --query_data ${PWD}/tests/data/lcms_query_library.txt \
  --reference_data ${PWD}/tests/data/full_GNPS_reference_library.txt \
  --spectrum_ID1 463514 \
  --spectrum_ID2 112312 \
  --chromatography_platform HRMS \
  --similarity_measure cosine \
  --spectrum_preprocessing_order FCNMWL \
  --high_quality_reference_library False \
  --mz_min 0 \
  --mz_max 9999999 \
  --int_min 0 \
  --int_max 9999999 \
  --window_size_centroiding 0.5 \
  --window_size_matching 0.5 \
  --noise_threshold 0.0 \
  --wf_mz 0.0 \
  --wf_intensity 1.0 \
  --LET_threshold 0.0 \
  --entropy_dimension 1.1 \
  --output_path ${PWD}/output_plotting_HRMS.pdf \

python plot_spectra_CLI.py \
  --query_data ${PWD}/data/gcms_query_library.txt \
  --reference_data ${PWD}/data/gcms_reference_library.txt \
  --spectrum_ID1 463514 \
  --spectrum_ID2 112312 \
  --chromatography_platform NRMS \
  --similarity_measure tsallis \
  --spectrum_preprocessing_order FCNMWL \
  --high_quality_reference_library False \
  --mz_min 0 \
  --mz_max 9999999 \
  --int_min 0 \
  --int_max 9999999 \
  --noise_threshold 0.0 \
  --wf_mz 0.0 \
  --wf_intensity 1.0 \
  --LET_threshold 0.0 \
  --entropy_dimension 1.1 \
  --output_path ${PWD}/output_plotting_NRMS.pdf \
```

An example of such a generated plot is seen below.

<br />

![image](https://github.com/user-attachments/assets/de22a402-1329-4bb3-a664-9423159264c8)

<br />

This plot compares two MS/MS spectra: Spectrum ID 1 (unknown, in blue) and Spectrum ID 2 (Hectochlorin M+H, in red). The top panel displays the untransformed spectra, while the bottom panel shows the transformed spectra following preprocessing steps. The footnote details are as follows:

-   Filtering: Given user-defined parameters (mz_min,mz_max),

-   Similarity Measure: Cosine -- The similarity measure used is cosine correlation.

-   Similarity Score: 0.9946 -- The cosine similarity score between the two transformed spectra.

-   Spectrum Preprocessing Order: FCNMWL -- The sequence of preprocessing steps applied: Filtering (F), Centroiding (C), Noise removal (N), Matching (M), Weight factor transformation (W), and Low-entropy transformation (L).

-   High Quality Reference Library: False -- Both query and reference spectra underwent the same preprocessing transformations.

-   Window Size (Centroiding): 0.5 -- A 0.5 Da window was used for centroiding peaks.

-   Window Size (Matching): 0.5 -- Peaks were aligned using a 0.5 Da m/z tolerance window.

-   Raw-Scale M/Z Range: [217.7, 628.8] -- The maximum and minimum of m/z values of peaks with non-zero intensities.

-   Raw-Scale Intensity Range: [3885.0, 5549140] -- The maximum and minimum of absolute non-zero intensity values of the raw spectra before normalization.

-   Noise Threshold: 0.0 -- No noise threshold was applied.

-   Weight Factors (m/z, intensity): (0.0, 1.0) -- Non-zero intensities were transformed using weights of 0.0 for m/z and 1.0 for intensity.

-   Low-Entropy Threshold: 0.0 -- No low-entropy transformation was applied.


<a name="shiny"></a>
### 3.6 Shiny application
PyCompound is also available as a Shiny application. The Shiny application offers the same functionality as the Python package and its CLI interface. Simply run the Python script src/pycompound_shiny.py with a command such as <shiny run --launch-browser pycompound_shiny.py> to launch the Shiny application. Alternatively, one can you the publicly available web version at [https://0199ee0c-c2ce-4fdc-5ade-623633df1622.share.connect.posit.cloud/](https://0199ee0c-c2ce-4fdc-5ade-623633df1622.share.connect.posit.cloud/). If you plan to perform some heavy computations such as parameter tuning on large datasets, we recommend either using the Python package, its CLI wrapper, or running the Shiny app on your local machine to take advantage of multithreading (which isn't offered on the POSIT-hosted Shiny app).


<a name="bugs-questions"></a>
## 4. Bugs/Questions?
If you notice any bugs in this software or have any questions, please create a new issue in this repository.

