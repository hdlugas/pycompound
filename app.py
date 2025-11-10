
from shiny import App, ui, reactive, render, req
from shiny.types import SilentException
from pathlib import Path
from contextlib import redirect_stdout, redirect_stderr
import contextlib
import subprocess
import traceback
import asyncio
import io
import os
import sys
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import netCDF4 as nc
from pyteomics import mgf, mzml
import ast
from numbers import Real
import logging
from scipy.optimize import differential_evolution
import scipy
import scipy.stats
from itertools import product
import json
import re
import urllib.parse
import urllib.request
import matplotlib

matplotlib.rcParams['svg.fonttype'] = 'none'

_LOG_QUEUE: asyncio.Queue[str] = asyncio.Queue()

_ADDUCT_PAT = re.compile(r"\s*(?:\[(M[^\]]+)\]|(M[+-][A-Za-z0-9]+)\+?)\s*$", re.IGNORECASE)

def start_log_consumer():
    if getattr(start_log_consumer, "_started", False):
        return
    start_log_consumer._started = True

    async def _consume():
        while True:
            s = await _LOG_QUEUE.get()
            match_log_rv.set(match_log_rv.get() + s)
            await reactive.flush()

    asyncio.create_task(_consume())


def start_log_consumer():
    if getattr(start_log_consumer, "_started", False):
        return
    start_log_consumer._started = True

    async def _consume():
        while True:
            s = await _LOG_QUEUE.get()
            match_log_rv.set(match_log_rv.get() + s)
            await reactive.flush()

    asyncio.create_task(_consume())



def _strip_adduct(name: str) -> str:
    return _ADDUCT_PAT.sub("", name).strip()

def get_pubchem_url(query: str) -> str:
    base_name = _strip_adduct(query)
    endpoint = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/" + urllib.parse.quote(base_name) + "/cids/TXT")
    try:
        with urllib.request.urlopen(endpoint, timeout=10) as r:
            txt = r.read().decode("utf-8").strip()
        cid = txt.splitlines()[0].strip()
        if cid.isdigit():
            return f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"
    except Exception:
        pass
    q = urllib.parse.quote(base_name)
    return f"https://pubchem.ncbi.nlm.nih.gov/#query={q}"



def build_library_from_raw_data(input_path=None, output_path=None, is_reference=False):
    if input_path is None:
        print('Error: please specify input_path (i.e. the path to the input mgf, mzML, cdf, json, or msp file). Mandatory argument.')
        sys.exit()

    if output_path is None:
        tmp = input_path.split('/')
        tmp = tmp[(len(tmp)-1)]
        basename = tmp.split('.')[0]
        output_path = f'{Path.cwd()}/{basename}.csv'
        print(f'Warning: no output_path specified, so library is written to {output_path}')

    if is_reference not in [True,False]:
        print('Error: is_reference must be either \'True\' or \'False\'.')
        sys.exit()

    last_three_chars = input_path[(len(input_path)-3):len(input_path)]
    last_four_chars = input_path[(len(input_path)-4):len(input_path)]
    if last_three_chars == 'mgf' or last_three_chars == 'MGF':
        input_file_type = 'mgf'
    elif last_four_chars == 'mzML' or last_four_chars == 'mzml' or last_four_chars == 'MZML':
        input_file_type = 'mzML'
    elif last_four_chars == 'json' or last_four_chars == 'JSON':
        input_file_type = 'json'
    elif last_three_chars == 'cdf' or last_three_chars == 'CDF':
        input_file_type = 'cdf'
    elif last_three_chars == 'msp' or last_three_chars == 'MSP':
        input_file_type = 'msp'
    else:
        print('ERROR: either an \'mgf\', \'mzML\', \'cdf\', \'json\', or \'msp\' file must be passed to --input_path')
        sys.exit()



def generate_plots_on_HRMS_data(query_data=None, reference_data=None, precursor_ion_mz=None, precursor_ion_mz_tolerance=None, ionization_mode=None, collision_energy=None, spectrum_ID1=None, spectrum_ID2=None, print_url_spectrum1='No', print_url_spectrum2='No', similarity_measure='cosine', weights={'Cosine':0.25,'Shannon':0.25,'Renyi':0.25,'Tsallis':0.25}, spectrum_preprocessing_order='FCNMWL', high_quality_reference_library=False, mz_min=0, mz_max=9999999, int_min=0, int_max=9999999, window_size_centroiding=0.5, window_size_matching=0.5, noise_threshold=0.0, wf_mz=0.0, wf_intensity=1.0, LET_threshold=0.0, entropy_dimension=1.1, y_axis_transformation='normalized', output_path=None, return_plot=False):

    if query_data is None:
        print('\nError: No argument passed to the mandatory query_data. Please pass the path to the TXT file of the query data.')
        sys.exit()
    else:
        extension = query_data.rsplit('.',1)
        extension = extension[(len(extension)-1)]
        if extension == 'mgf' or extension == 'MGF' or extension == 'mzML' or extension == 'mzml' or extension == 'MZML' or extension == 'cdf' or extension == 'CDF' or extension == 'msp' or extension == 'MSP' or extension == 'json' or extension == 'JSON':
            output_path_tmp = query_data[:-3] + 'txt'
            #build_library_from_raw_data(input_path=query_data, output_path=output_path_tmp, is_reference=True)
            build_library_from_raw_data(input_path=query_data, output_path=output_path_tmp, is_reference=False)
            df_query = pd.read_csv(output_path_tmp, sep='\t')
        if extension == 'txt' or extension == 'TXT':
            df_query = pd.read_csv(query_data, sep='\t')
        unique_query_ids = df_query['id'].unique().tolist()
        unique_query_ids = [str(tmp) for tmp in unique_query_ids]

    if reference_data is None:
        print('\nError: No argument passed to the mandatory reference_data. Please pass the path to the TXT file of the reference data.')
        sys.exit()
    else:
        extension = reference_data.rsplit('.',1)
        extension = extension[(len(extension)-1)]
        if extension == 'mgf' or extension == 'MGF' or extension == 'mzML' or extension == 'mzml' or extension == 'MZML' or extension == 'cdf' or extension == 'CDF' or extension == 'msp' or extension == 'MSP' or extension == 'json' or extension == 'JSON':
            output_path_tmp = reference_data[:-3] + 'txt'
            build_library_from_raw_data(input_path=reference_data, output_path=output_path_tmp, is_reference=True)
            df_reference = pd.read_csv(output_path_tmp, sep='\t')
        if extension == 'txt' or extension == 'TXT':
            df_reference = pd.read_csv(reference_data, sep='\t')
            cols_tmp = df_reference.columns.tolist()
            if 'precursor_ion_mz' in cols_tmp and 'ionization_mode' in cols_tmp and 'collision_energy' in cols_tmp:
                if precursor_ion_mz is not None and precursor_ion_mz_tolerance is not None:
                    df_reference = df_reference.loc[(df_reference['precursor_ion_mz'] > (precursor_ion_mz-precursor_ion_mz_tolerance) & df_reference['precursor_ion_mz'] < (precursor_ion_mz+precursor_ion_mz_tolerance))]
                if ionization_mode is not None:
                    df_reference = df_reference.loc[df_reference['ionization_mode'==ionization_mode]]
                if collision_energy is not None:
                    df_reference = df_reference.loc[df_reference['collision_energy'==collision_energy]]
                df_reference = df_reference.drop(columns=['precursor_ion_mz','ionization_mode','collision_energy'])
        unique_reference_ids = df_reference['id'].unique().tolist()
        unique_reference_ids = [str(tmp) for tmp in unique_reference_ids]

    if spectrum_ID1 is not None:
        spectrum_ID1 = str(spectrum_ID1)
    else:
        spectrum_ID1 = str(df_query['id'].iloc[0])
        print('No argument passed to spectrum_ID1; using the first spectrum in query_data.')

    if spectrum_ID2 is not None:
        spectrum_ID2 = str(spectrum_ID2)
    else:
        spectrum_ID2 = str(df_reference['id'].iloc[0])
        print('No argument passed to spectrum_ID2; using the first spectrum in reference_data.')

    if spectrum_preprocessing_order is not None:
        spectrum_preprocessing_order = list(spectrum_preprocessing_order)
    else:
        spectrum_preprocessing_order = ['F', 'C', 'N', 'M', 'W', 'L']
    if 'M' not in spectrum_preprocessing_order:
        print(f'Error: \'M\' must be a character in spectrum_preprocessing_order.')
        sys.exit()
    if 'C' in spectrum_preprocessing_order:
        if spectrum_preprocessing_order.index('C') > spectrum_preprocessing_order.index('M'):
            print(f'Error: \'C\' must come before \'M\' in spectrum_preprocessing_order.')
            sys.exit()
    if set(spectrum_preprocessing_order) - {'F','C','N','M','W','L'}:
        print(f'Error: spectrum_preprocessing_order must contain only \'C\', \'F\', \'M\', \'N\', \'L\', \'W\'.')
        sys.exit()

    if similarity_measure not in ['cosine','shannon','renyi','tsallis','mixture','jaccard','dice','3w_jaccard','sokal_sneath','binary_cosine','mountford','mcconnaughey','driver_kroeber','simpson','braun_banquet','fager_mcgowan','kulczynski','intersection','hamming','hellinger']:
        print('\nError: similarity_measure must be either cosine, shannon, renyi, tsallis, mixture, jaccard, dice, 3w_jaccard, sokal_sneath, binary_cosine, mountford, mcconnaughey, driver_kroeber, simpson, braun_banquet, fager_mcgowan, kulczynski, intersection, hamming, or hellinger.')
        sys.exit()

    if isinstance(int_min,int) is True:
        int_min = float(int_min)
    if isinstance(int_max,int) is True:
        int_max = float(int_max)
    if isinstance(mz_min,int) is False or isinstance(mz_max,int) is False or isinstance(int_min,float) is False or isinstance(int_max,float) is False:
        print('Error: mz_min must be a non-negative integer, mz_max must be a positive integer, int_min must be a non-negative float, and int_max must be a positive float')
        sys.exit()
    if mz_min < 0:
        print('\nError: mz_min should be a non-negative integer')
        sys.exit()
    if mz_max <= 0:
        print('\nError: mz_max should be a positive integer')
        sys.exit()
    if int_min < 0:
        print('\nError: int_min should be a non-negative float')
        sys.exit()
    if int_max <= 0:
        print('\nError: int_max should be a positive float')
        sys.exit()

    if isinstance(window_size_centroiding,float) is False or window_size_centroiding <= 0.0:
        print('Error: window_size_centroiding must be a positive float.')
        sys.exit()
    if isinstance(window_size_matching,float) is False or window_size_matching<= 0.0:
        print('Error: window_size_matching must be a positive float.')
        sys.exit()

    if isinstance(noise_threshold,int) is True:
        noise_threshold = float(noise_threshold)
    if isinstance(noise_threshold,float) is False or noise_threshold < 0:
        print('Error: noise_threshold must be a positive float.')
        sys.exit()

    if isinstance(wf_intensity,int) is True:
        wf_intensity = float(wf_intensity)
    if isinstance(wf_mz,int) is True:
        wf_mz = float(wf_mz)
    if isinstance(wf_intensity,float) is False or isinstance(wf_mz,float) is False:
        print('Error: wf_mz and wf_intensity must be integers or floats')
        sys.exit()

    if entropy_dimension <= 0:
        print('\nError: entropy_dimension should be a positive float')
        sys.exit()
    else:
        q = entropy_dimension

    normalization_method = 'standard' #consider including additional normalization methods to transform intensities into a probability distribution; softmax results in many numerical errors/warnings

    if y_axis_transformation not in ['normalized','none','log10','sqrt']:
        print('Error: y_axis_transformation must be either \'normalized\', \'none\', \'log10\', or \'sqrt\'.')
        sys.exit()

    if output_path is None:
        print(f'Warning: plots will be saved to the svg ./spectrum1_{spectrum_ID1}_spectrum2_{spectrum_ID2}.svg in the current working directory.')
        output_path = f'{Path.cwd()}/spectrum1_{spectrum_ID1}_spectrum2_{spectrum_ID2}.svg'


    if spectrum_ID1 in unique_query_ids and spectrum_ID2 in unique_query_ids:
        query_idx = unique_query_ids.index(spectrum_ID1)
        reference_idx = unique_query_ids.index(spectrum_ID2)
        q_idxs_tmp = np.where(df_query.iloc[:,0].astype(str) == unique_query_ids[query_idx])[0]
        r_idxs_tmp = np.where(df_query.iloc[:,0].astype(str) == unique_query_ids[reference_idx])[0]
        q_spec = np.asarray(pd.concat([df_query.iloc[q_idxs_tmp,1], df_query.iloc[q_idxs_tmp,2]], axis=1).reset_index(drop=True))
        r_spec = np.asarray(pd.concat([df_query.iloc[r_idxs_tmp,1], df_query.iloc[r_idxs_tmp,2]], axis=1).reset_index(drop=True))
    elif spectrum_ID1 in unique_reference_ids and spectrum_ID2 in unique_reference_ids:
        query_idx = unique_reference_ids.index(spectrum_ID1)
        reference_idx = unique_reference_ids.index(spectrum_ID2)
        q_idxs_tmp = np.where(df_reference.iloc[:,0].astype(str) == unique_reference_ids[query_idx])[0]
        r_idxs_tmp = np.where(df_reference.iloc[:,0].astype(str) == unique_reference_ids[reference_idx])[0]
        q_spec = np.asarray(pd.concat([df_reference.iloc[q_idxs_tmp,1], df_reference.iloc[q_idxs_tmp,2]], axis=1).reset_index(drop=True))
        r_spec = np.asarray(pd.concat([df_reference.iloc[r_idxs_tmp,1], df_reference.iloc[r_idxs_tmp,2]], axis=1).reset_index(drop=True))
    else:
        if spectrum_ID1 in unique_reference_ids and spectrum_ID2 in unique_query_ids:
            spec_tmp = spectrum_ID1
            spectrum_ID1 = spectrum_ID2
            spectrum_ID2 = spec_tmp
        query_idx = unique_query_ids.index(spectrum_ID1)
        reference_idx = unique_reference_ids.index(spectrum_ID2)
        q_idxs_tmp = np.where(df_query['id'].astype(str) == unique_query_ids[query_idx])[0]
        r_idxs_tmp = np.where(df_reference['id'].astype(str) == unique_reference_ids[reference_idx])[0]
        q_spec = np.asarray(pd.concat([df_query['mz_ratio'].iloc[q_idxs_tmp], df_query['intensity'].iloc[q_idxs_tmp]], axis=1).reset_index(drop=True))
        r_spec = np.asarray(pd.concat([df_reference['mz_ratio'].iloc[r_idxs_tmp], df_reference['intensity'].iloc[r_idxs_tmp]], axis=1).reset_index(drop=True))


    q_spec_pre_trans = q_spec.copy()
    r_spec_pre_trans = r_spec.copy()
    q_spec_pre_trans[:,1] = q_spec_pre_trans[:,1].astype(float)
    r_spec_pre_trans[:,1] = r_spec_pre_trans[:,1].astype(float)

    if y_axis_transformation == 'normalized':
        q_spec_pre_trans[:,1] = q_spec_pre_trans[:,1] / np.max(q_spec_pre_trans[:,1])
        r_spec_pre_trans[:,1] = r_spec_pre_trans[:,1] / np.max(r_spec_pre_trans[:,1])
        ylab = 'Normalized Intensity'
    elif y_axis_transformation == 'log10':
        q_spec_pre_trans[:,1] = np.log10(np.array(q_spec_pre_trans[:,1]+1,dtype=float))
        r_spec_pre_trans[:,1] = np.log10(np.array(r_spec_pre_trans[:,1]+1,dtype=float))
        ylab = 'log10(Intensity)'
    elif y_axis_transformation == 'sqrt':
        q_spec_pre_trans[:,1] = np.sqrt(np.array(q_spec_pre_trans[:,1],dtype=float))
        r_spec_pre_trans[:,1] = np.sqrt(np.array(r_spec_pre_trans[:,1],dtype=float))
        ylab = 'sqrt(Intensity)'
    else:
        ylab = 'Raw Intensity'

    fig, axes = plt.subplots(nrows=2, ncols=1)

    plt.subplot(2,1,1)
    plt.vlines(x=q_spec_pre_trans[:,0], ymin=[0]*q_spec_pre_trans.shape[0], ymax=q_spec_pre_trans[:,1], linewidth=3, color='blue', label=f'Spectrum ID 1: {spectrum_ID1}')
    plt.vlines(x=r_spec_pre_trans[:,0], ymin=[0]*r_spec_pre_trans.shape[0], ymax=-r_spec_pre_trans[:,1], linewidth=3, color='red', label=f'Spectrum ID 2: {spectrum_ID2}')
    plt.xlabel('m/z',fontsize=7)
    plt.ylabel(ylab, fontsize=7)
    plt.xticks(fontsize=7)
    plt.yticks(fontsize=7)
    plt.title('Untransformed Spectra', fontsize=10)

    mz_min_tmp_q = round(q_spec[:,0].min(),1)
    mz_min_tmp_r = round(r_spec[:,0].min(),1)
    int_min_tmp_q = round(q_spec[:,1].min(),1)
    int_min_tmp_r = round(r_spec[:,1].min(),1)
    mz_max_tmp_q = round(q_spec[:,0].max(),1)
    mz_max_tmp_r = round(r_spec[:,0].max(),1)
    int_max_tmp_q = round(q_spec[:,1].max(),1)
    int_max_tmp_r = round(r_spec[:,1].max(),1)
    mz_min_tmp = min([mz_min_tmp_q,mz_min_tmp_r])
    mz_max_tmp = min([mz_max_tmp_q,mz_max_tmp_r])
    int_min_tmp = min([int_min_tmp_q,int_min_tmp_r])
    int_max_tmp = max([int_max_tmp_q,int_max_tmp_r])

    is_matched = False
    for transformation in spectrum_preprocessing_order:
        if transformation == 'C' and q_spec.shape[0] > 1 and r_spec.shape[1] > 1:
            q_spec = centroid_spectrum(q_spec, window_size=window_size_centroiding) 
            r_spec = centroid_spectrum(r_spec, window_size=window_size_centroiding) 
        if transformation == 'M' and q_spec.shape[0] > 1 and r_spec.shape[1] > 1:
            m_spec = match_peaks_in_spectra(spec_a=q_spec, spec_b=r_spec, window_size=window_size_matching)
            q_spec = m_spec[:,0:2]
            r_spec = m_spec[:,[0,2]]
            is_matched = True
        if transformation == 'W' and q_spec.shape[0] > 1 and r_spec.shape[1] > 1:
            q_spec[:,1] = wf_transform(q_spec[:,0], q_spec[:,1], wf_mz, wf_intensity)
            r_spec[:,1] = wf_transform(r_spec[:,0], r_spec[:,1], wf_mz, wf_intensity)
        if transformation == 'L' and q_spec.shape[0] > 1 and r_spec.shape[1] > 1:
            q_spec[:,1] = LE_transform(q_spec[:,1], LET_threshold, normalization_method=normalization_method)
            r_spec[:,1] = LE_transform(r_spec[:,1], LET_threshold, normalization_method=normalization_method)
        if transformation == 'N' and q_spec.shape[0] > 1 and r_spec.shape[1] > 1:
            q_spec = remove_noise(q_spec, nr = noise_threshold)
            if high_quality_reference_library == False or high_quality_reference_library == 'False':
                r_spec = remove_noise(r_spec, nr = noise_threshold)
        if transformation == 'F' and q_spec.shape[0] > 1 and r_spec.shape[1] > 1:
            q_spec = filter_spec_lcms(q_spec, mz_min = mz_min, mz_max = mz_max, int_min = int_min, int_max = int_max, is_matched = is_matched)
            if high_quality_reference_library == False or high_quality_reference_library == 'False':
                r_spec = filter_spec_lcms(r_spec, mz_min = mz_min, mz_max = mz_max, int_min = int_min, int_max = int_max, is_matched = is_matched)

    q_ints = q_spec[:,1]
    r_ints = r_spec[:,1]

    if np.sum(q_ints) != 0 and np.sum(r_ints) != 0 and q_spec.shape[0] > 1 and r_spec.shape[1] > 1:
        similarity_score = get_similarity(similarity_measure, q_ints, r_ints, weights, entropy_dimension)
    else:
        similarity_score = 0

    plt.subplot(2,1,2)

    if q_spec.shape[0] > 1:
        if np.max(q_spec[:,1]) == 0 or np.max(r_spec[:,1]) == 0:
            plt.text(0.5, 0.5, 'The query and/or reference spectrum has no non-zero intensities after transformations.\n Change transformation parameters.', ha='center', va='center', fontsize=7, color='black')
            plt.xticks([])
            plt.yticks([])
        else:
            if y_axis_transformation == 'normalized':
                q_spec[:,1] = q_spec[:,1] / np.max(q_spec[:,1])
                r_spec[:,1] = r_spec[:,1] / np.max(r_spec[:,1])
                ylab='Normalized Intensity'
            elif y_axis_transformation == 'log10':
                q_spec[:,1] = np.log10(q_spec[:,1]+1)
                r_spec[:,1] = np.log10(r_spec[:,1]+1)
                ylab='log10(Intensity)'
            elif y_axis_transformation == 'sqrt':
                q_spec[:,1] = np.sqrt(q_spec[:,1])
                r_spec[:,1] = np.sqrt(r_spec[:,1])
                ylab='sqrt(Intensity)'
            else:
                ylab = 'Raw Intensity'
            plt.vlines(x=q_spec[:,0], ymin=[0]*q_spec.shape[0], ymax=q_spec[:,1], linewidth=3, color='blue')
            plt.vlines(x=r_spec[:,0], ymin=[0]*r_spec.shape[0], ymax=-r_spec[:,1], linewidth=3, color='red')
            plt.xlabel('m/z', fontsize=7)
            plt.ylabel(ylab, fontsize=7)
            plt.xticks(fontsize=7)
            plt.yticks(fontsize=7)
            plt.title(f'Transformed Spectra', fontsize=10)
    else:
        plt.text(0.5, 0.5, 'All points in the spectra were removed during preprocessing. \nChange the spectrum_preprocesing_order and/or change other spectrum-preprocessing parameters.', ha='center', va='center', fontsize=7, color='black')
        plt.xticks([])
        plt.yticks([])

    plt.subplots_adjust(top=0.8, hspace=0.92, bottom=0.3)
    plt.figlegend(loc='upper center')

    fig.text(0.05, 0.20, f'Similarity Measure: {similarity_measure.capitalize()}', fontsize=7)
    fig.text(0.05, 0.17, f'Similarity Score: {round(similarity_score, 4)}', fontsize=7)
    fig.text(0.05, 0.14, f"Spectrum Preprocessing Order: {''.join(spectrum_preprocessing_order)}", fontsize=7)
    fig.text(0.05, 0.11, f'High Quality Reference Library: {str(high_quality_reference_library)}', fontsize=7)
    fig.text(0.05, 0.08, f'Window Size (Centroiding): {window_size_centroiding}', fontsize=7)
    fig.text(0.05, 0.05, f'Window Size (Matching): {window_size_matching}', fontsize=7)
    if similarity_measure == 'mixture':
        fig.text(0.05, 0.02, f'Weights for mixture similarity: {weights}', fontsize=7)

    fig.text(0.40, 0.20, f'Raw-Scale M/Z Range: [{mz_min_tmp},{mz_max_tmp}]', fontsize=7)
    fig.text(0.40, 0.17, f'Raw-Scale Intensity Range: [{int_min_tmp},{int_max_tmp}]', fontsize=7)
    fig.text(0.40, 0.14, f'Noise Threshold: {noise_threshold}', fontsize=7)
    fig.text(0.40, 0.11, f'Weight Factors (m/z,intensity): ({wf_mz},{wf_intensity})', fontsize=7)
    fig.text(0.40, 0.08, f'Low-Entropy Threshold: {LET_threshold}', fontsize=7)

    if print_url_spectrum1 == 'Yes' and print_url_spectrum2 == 'Yes':
        url_tmp1 = get_pubchem_url(query=spectrum_ID1)
        url_tmp2 = get_pubchem_url(query=spectrum_ID2)
        t1 = fig.text(0.40, 0.05, f'PubChem URL for {spectrum_ID1}: {url_tmp1}', fontsize=7)
        t2 = fig.text(0.40, 0.02, f'PubChem URL for {spectrum_ID2}: {url_tmp2}', fontsize=7)
        t1.set_url(url_tmp1)
        t2.set_url(url_tmp2)

    if print_url_spectrum1 == 'Yes' and print_url_spectrum2 == 'No':
        url_tmp1 = get_pubchem_url(query=spectrum_ID1)
        t1 = fig.text(0.40, 0.05, f'PubChem URL for {spectrum_ID1}: {url_tmp1}', fontsize=7)
        t1.set_url(url_tmp1)

    if print_url_spectrum1 == 'No' and print_url_spectrum2 == 'Yes':
        url_tmp2 = get_pubchem_url(query=spectrum_ID2)
        t2 = fig.text(0.40, 0.05, f'PubChem URL for {spectrum_ID2}: {url_tmp2}', fontsize=7)
        t2.set_url(url_tmp2)

    fig.savefig(output_path, format='svg')

    if return_plot == True:
        return fig




def generate_plots_on_NRMS_data(query_data=None, reference_data=None, spectrum_ID1=None, spectrum_ID2=None, print_url_spectrum1='No', print_url_spectrum2='No', similarity_measure='cosine', weights={'Cosine':0.25,'Shannon':0.25,'Renyi':0.25,'Tsallis':0.25}, spectrum_preprocessing_order='FNLW', high_quality_reference_library=False, mz_min=0, mz_max=9999999, int_min=0, int_max=9999999, noise_threshold=0.0, wf_mz=0.0, wf_intensity=1.0, LET_threshold=0.0, entropy_dimension=1.1, y_axis_transformation='normalized', output_path=None, return_plot=False):

    if query_data is None:
        print('\nError: No argument passed to the mandatory query_data. Please pass the path to the TXT file of the query data.')
        sys.exit()
    else:
        extension = query_data.rsplit('.',1)
        extension = extension[(len(extension)-1)]
        if extension == 'mgf' or extension == 'MGF' or extension == 'mzML' or extension == 'mzml' or extension == 'MZML' or extension == 'cdf' or extension == 'CDF' or extension == 'msp' or extension == 'MSP' or extension == 'json' or extension == 'JSON':
            output_path_tmp = query_data[:-3] + 'txt'
            build_library_from_raw_data(input_path=query_data, output_path=output_path_tmp, is_reference=False)
            df_query = pd.read_csv(output_path_tmp, sep='\t')
        if extension == 'txt' or extension == 'TXT':
            df_query = pd.read_csv(query_data, sep='\t')
        unique_query_ids = df_query['id'].unique()

    if reference_data is None:
        print('\nError: No argument passed to the mandatory reference_data. Please pass the path to the TXT file of the reference data.')
        sys.exit()
    else:
        extension = reference_data.rsplit('.',1)
        extension = extension[(len(extension)-1)]
        if extension == 'mgf' or extension == 'MGF' or extension == 'mzML' or extension == 'mzml' or extension == 'MZML' or extension == 'cdf' or extension == 'CDF' or extension == 'msp' or extension == 'MSP' or extension == 'json' or extension == 'JSON':
            output_path_tmp = reference_data[:-3] + 'txt'
            build_library_from_raw_data(input_path=reference_data, output_path=output_path_tmp, is_reference=True)
            df_reference = pd.read_csv(output_path_tmp, sep='\t')
        if extension == 'txt' or extension == 'TXT':
            df_reference = pd.read_csv(reference_data, sep='\t')
            unique_reference_ids = df_reference['id'].unique()


    if spectrum_ID1 is not None:
        spectrum_ID1 = str(spectrum_ID1)
    else:
        spectrum_ID1 = str(df_query.iloc[0,0])
        print('No argument passed to spectrum_ID1; using the first spectrum in query_data.')

    if spectrum_ID2 is not None:
        spectrum_ID2 = str(spectrum_ID2)
    else:
        spectrum_ID2 = str(df_reference.iloc[0,0])
        print('No argument passed to spectrum_ID2; using the first spectrum in reference_data.')

    if spectrum_preprocessing_order is not None:
        spectrum_preprocessing_order = list(spectrum_preprocessing_order)
    else:
        spectrum_preprocessing_order = ['F','N','W','L']
    if set(spectrum_preprocessing_order) - {'F','N','W','L'}:
        print(f'Error: spectrum_preprocessing_order must contain only \'F\', \'N\', \'W\', \'L\'.')
        sys.exit()

    if similarity_measure not in ['cosine','shannon','renyi','tsallis','mixture','jaccard','dice','3w_jaccard','sokal_sneath','binary_cosine','mountford','mcconnaughey','driver_kroeber','simpson','braun_banquet','fager_mcgowan','kulczynski','intersection','hamming','hellinger']:
        print('\nError: similarity_measure must be either cosine, shannon, renyi, tsallis, mixture, jaccard, dice, 3w_jaccard, sokal_sneath, binary_cosine, mountford, mcconnaughey, driver_kroeber, simpson, braun_banquet, fager_mcgowan, kulczynski, intersection, hamming, or hellinger.')
        sys.exit()

    if isinstance(int_min,int) is True:
        int_min = float(int_min)
    if isinstance(int_max,int) is True:
        int_max = float(int_max)
    if isinstance(mz_min,int) is False or isinstance(mz_max,int) is False or isinstance(int_min,float) is False or isinstance(int_max,float) is False:
        print('Error: mz_min must be a non-negative integer, mz_max must be a positive integer, int_min must be a non-negative float, and int_max must be a positive float')
        sys.exit()
    if mz_min < 0:
        print('\nError: mz_min should be a non-negative integer')
        sys.exit()
    if mz_max <= 0:
        print('\nError: mz_max should be a positive integer')
        sys.exit()
    if int_min < 0:
        print('\nError: int_min should be a non-negative float')
        sys.exit()
    if int_max <= 0:
        print('\nError: int_max should be a positive float')
        sys.exit()

    if isinstance(noise_threshold,int) is True:
        noise_threshold = float(noise_threshold)
    if isinstance(noise_threshold,float) is False or noise_threshold < 0:
        print('Error: noise_threshold must be a positive float.')
        sys.exit()

    if isinstance(wf_intensity,int) is True:
        wf_intensity = float(wf_intensity)
    if isinstance(wf_mz,int) is True:
        wf_mz = float(wf_mz)
    if isinstance(wf_intensity,float) is False or isinstance(wf_mz,float) is False:
        print('Error: wf_mz and wf_intensity must be integers or floats')
        sys.exit()

    if entropy_dimension <= 0:
        print('\nError: entropy_dimension should be a positive float')
        sys.exit()
    else:
        q = entropy_dimension

    normalization_method = 'standard' #consider including additional normalization methods to transform intensities into a probability distribution; softmax results in many numerical errors/warnings

    if y_axis_transformation not in ['normalized','none','log10','sqrt']:
        print('Error: y_axis_transformation must be either \'normalized\', \'none\', \'log10\', or \'sqrt\'.')
        sys.exit()

    if output_path is None:
        print(f'Warning: plots will be saved to the svg ./spectrum1_{spectrum_ID1}_spectrum2_{spectrum_ID2}.svg in the current working directory.')
        output_path = f'{Path.cwd()}/spectrum1_{spectrum_ID1}_spectrum2_{spectrum_ID2}.svg'

    min_mz = np.min([df_query['mz_ratio'].min(), df_reference['mz_ratio'].min()])
    max_mz = np.max([df_query['mz_ratio'].max(), df_reference['mz_ratio'].max()])
    mzs = np.linspace(min_mz,max_mz,(max_mz-min_mz+1))

    unique_query_ids = df_query['id'].unique().tolist()
    unique_reference_ids = df_reference['id'].unique().tolist()
    unique_query_ids = [str(ID) for ID in unique_query_ids]
    unique_reference_ids = [str(ID) for ID in unique_reference_ids]
    common_IDs = np.intersect1d([str(ID) for ID in unique_query_ids], [str(ID) for ID in unique_reference_ids])
    if len(common_IDs) > 0:
        print(f'Warning: the query and reference library have overlapping IDs: {common_IDs}')

    if spectrum_ID1 in unique_query_ids and spectrum_ID2 in unique_query_ids:
        q_idxs_tmp = np.where(df_query.iloc[:,0].astype(str) == spectrum_ID1)[0]
        r_idxs_tmp = np.where(df_query.iloc[:,0].astype(str) == spectrum_ID2)[0]
        q_spec = np.asarray(pd.concat([df_query.iloc[q_idxs_tmp,1], df_query.iloc[q_idxs_tmp,2]], axis=1).reset_index(drop=True))
        r_spec = np.asarray(pd.concat([df_query.iloc[r_idxs_tmp,1], df_query.iloc[r_idxs_tmp,2]], axis=1).reset_index(drop=True))
    elif spectrum_ID1 in unique_reference_ids and spectrum_ID2 in unique_reference_ids:
        q_idxs_tmp = np.where(df_reference.iloc[:,0].astype(str) == spectrum_ID1)[0]
        r_idxs_tmp = np.where(df_reference.iloc[:,0].astype(str) == spectrum_ID2)[0]
        q_spec = np.asarray(pd.concat([df_reference.iloc[q_idxs_tmp,1], df_reference.iloc[q_idxs_tmp,2]], axis=1).reset_index(drop=True))
        r_spec = np.asarray(pd.concat([df_reference.iloc[r_idxs_tmp,1], df_reference.iloc[r_idxs_tmp,2]], axis=1).reset_index(drop=True))
    else:
        if spectrum_ID1 in unique_reference_ids and spectrum_ID2 in unique_query_ids:
            spec_tmp = spectrum_ID1
            spectrum_ID1 = spectrum_ID2
            spectrum_ID2 = spec_tmp
        q_idxs_tmp = np.where(df_query['id'].astype(str) == spectrum_ID1)[0]
        r_idxs_tmp = np.where(df_reference['id'].astype(str) == spectrum_ID2)[0]
        q_spec = np.asarray(pd.concat([df_query['mz_ratio'].iloc[q_idxs_tmp], df_query['intensity'].iloc[q_idxs_tmp]], axis=1).reset_index(drop=True))
        r_spec = np.asarray(pd.concat([df_reference['mz_ratio'].iloc[r_idxs_tmp], df_reference['intensity'].iloc[r_idxs_tmp]], axis=1).reset_index(drop=True))

    q_spec = convert_spec(q_spec,mzs)
    r_spec = convert_spec(r_spec,mzs)

    int_min_tmp_q = min(q_spec[q_spec[:,1].nonzero(),1][0])
    int_min_tmp_r = min(r_spec[r_spec[:,1].nonzero(),1][0])
    int_max_tmp_q = max(q_spec[q_spec[:,1].nonzero(),1][0])
    int_max_tmp_r = max(r_spec[r_spec[:,1].nonzero(),1][0])
    int_min_tmp = int(min([int_min_tmp_q,int_min_tmp_r]))
    int_max_tmp = int(max([int_max_tmp_q,int_max_tmp_r]))
    
    fig, axes = plt.subplots(nrows=2, ncols=1)

    plt.subplot(2,1,1)

    if np.max(q_spec[:,1]) == 0 or np.max(r_spec[:,1]) == 0:
        plt.text(0.5, 0.5, 'The query and/or reference spectrum has no non-zero intensities after transformations.\n Change transformation parameters.', ha='center', va='center', fontsize=7, color='black')
        plt.xticks([])
        plt.yticks([])
    else:
        q_spec_pre_trans = q_spec.copy()
        r_spec_pre_trans = r_spec.copy()
        q_spec_pre_trans[:,1] = q_spec_pre_trans[:,1].astype(float)
        r_spec_pre_trans[:,1] = r_spec_pre_trans[:,1].astype(float)

        if y_axis_transformation == 'normalized':
            q_spec_pre_trans[:,1] = q_spec_pre_trans[:,1] / np.max(q_spec_pre_trans[:,1])
            r_spec_pre_trans[:,1] = r_spec_pre_trans[:,1] / np.max(r_spec_pre_trans[:,1])
            ylab = 'Normalized Intensity'
        elif y_axis_transformation == 'log10':
            q_spec_pre_trans[:,1] = np.log10(q_spec_pre_trans[:,1]+1)
            r_spec_pre_trans[:,1] = np.log10(r_spec_pre_trans[:,1]+1)
            ylab = 'log10(Intensity)'
        elif y_axis_transformation == 'sqrt':
            q_spec_pre_trans[:,1] = np.sqrt(q_spec_pre_trans[:,1])
            r_spec_pre_trans[:,1] = np.sqrt(r_spec_pre_trans[:,1])
            ylab = 'sqrt(Intensity)'
        else:
            ylab = 'Raw Intensity'
        plt.vlines(x=q_spec_pre_trans[:,0], ymin=[0]*len(q_spec_pre_trans[:,0]), ymax=q_spec_pre_trans[:,1], linewidth=3, color='blue', label=f'Spectrum ID1: {spectrum_ID1}')
        plt.vlines(x=r_spec_pre_trans[:,0], ymin=[0]*len(r_spec_pre_trans[:,0]), ymax=-r_spec_pre_trans[:,1], linewidth=3, color='red', label=f'Spectrum ID2: {spectrum_ID2}')
        plt.xlabel('m/z',fontsize=7)
        plt.ylabel(ylab, fontsize=7)
        plt.xticks(fontsize=7)
        plt.yticks(fontsize=7)
        plt.title('Untransformed Query and Reference Spectra', fontsize=10)

    for transformation in spectrum_preprocessing_order:
        if transformation == 'W':
            q_spec[:,1] = wf_transform(q_spec[:,0], q_spec[:,1], wf_mz, wf_intensity)
            r_spec[:,1] = wf_transform(r_spec[:,0], r_spec[:,1], wf_mz, wf_intensity)
        if transformation == 'L':
            q_spec[:,1] = LE_transform(q_spec[:,1], LET_threshold, normalization_method)
            r_spec[:,1] = LE_transform(r_spec[:,1], LET_threshold, normalization_method)
        if transformation == 'N':
            q_spec = remove_noise(q_spec, nr = noise_threshold)
            if high_quality_reference_library == False or high_quality_reference_library == 'False':
                r_spec = remove_noise(r_spec, nr = noise_threshold)
        if transformation == 'F':
            q_spec = filter_spec_gcms(q_spec, mz_min = mz_min, mz_max = mz_max, int_min = int_min, int_max = int_max)
            if high_quality_reference_library == False or high_quality_reference_library == 'False':
                r_spec = filter_spec_gcms(r_spec, mz_min = mz_min, mz_max = mz_max, int_min = int_min, int_max = int_max)

    if q_spec.shape[0] > 1:
        similarity_score = get_similarity(similarity_measure, q_spec[:,1], r_spec[:,1], weights, entropy_dimension)
    else:
        similarity_score = 0

 
    plt.subplot(2,1,2)

    if q_spec.shape[0] == 0 or r_spec.shape[0] == 0:
        plt.text(0.5, 0.5, 'The query and/or reference spectrum has no ion fragments left after transformations.\n Change transformation parameters.', ha='center', va='center', fontsize=7, color='black')
        plt.xticks([])
        plt.yticks([])
    elif np.max(q_spec[:,1]) == 0 or np.max(r_spec[:,1]) == 0:
        plt.text(0.5, 0.5, 'The query and/or reference spectrum has no non-zero intensities after transformations.\n Change transformation parameters.', ha='center', va='center', fontsize=7, color='black')
        plt.xticks([])
        plt.yticks([])
    else:
        if y_axis_transformation == 'normalized':
            q_spec[:,1] = q_spec[:,1] / np.max(q_spec[:,1])
            r_spec[:,1] = r_spec[:,1] / np.max(r_spec[:,1])
            ylab='Normalized Intensity'
        elif y_axis_transformation == 'log10':
            q_spec[:,1] = np.log10(q_spec[:,1]+1)
            r_spec[:,1] = np.log10(r_spec[:,1]+1)
            ylab='log10(Intensity)'
        elif y_axis_transformation == 'sqrt':
            q_spec[:,1] = np.sqrt(q_spec[:,1])
            r_spec[:,1] = np.sqrt(r_spec[:,1])
            ylab='sqrt(Intensity)'
        else:
            ylab = 'Raw Intensity'
        plt.vlines(x=mzs, ymin=[0]*len(mzs), ymax=q_spec[:,1], linewidth=3, color='blue')
        plt.vlines(x=mzs, ymin=[0]*len(mzs), ymax=-r_spec[:,1], linewidth=3, color='red')
        plt.xlabel('m/z', fontsize=7)
        plt.ylabel(ylab, fontsize=7)
        plt.xticks(fontsize=7)
        plt.yticks(fontsize=7)
        plt.title(f'Transformed Query and Reference Spectra', fontsize=10)

    plt.subplots_adjust(top=0.8, hspace=0.92, bottom=0.3)
    plt.figlegend(loc='upper center')

    fig.text(0.05, 0.20, f'Similarity Measure: {similarity_measure.capitalize()}', fontsize=7)
    fig.text(0.05, 0.17, f'Similarity Score: {round(similarity_score, 4)}', fontsize=7)
    fig.text(0.05, 0.14, f"Spectrum Preprocessing Order: {''.join(spectrum_preprocessing_order)}", fontsize=7)
    fig.text(0.05, 0.11, f'High Quality Reference Library: {str(high_quality_reference_library)}', fontsize=7)
    fig.text(0.05, 0.08, f'Weight Factors (m/z,intensity): ({wf_mz},{wf_intensity})', fontsize=7)
    if similarity_measure == 'mixture':
        fig.text(0.05, 0.05, f'Weights for mixture similarity: {weights}', fontsize=7)

    fig.text(0.40, 0.20, f'Raw-Scale M/Z Range: [{min_mz},{max_mz}]', fontsize=7)
    fig.text(0.40, 0.17, f'Raw-Scale Intensity Range: [{int_min_tmp},{int_max_tmp}]', fontsize=7)
    fig.text(0.40, 0.14, f'Noise Threshold: {noise_threshold}', fontsize=7)
    fig.text(0.40, 0.11, f'Low-Entropy Threshold: {LET_threshold}', fontsize=7)

    if print_url_spectrum1 == 'Yes' and print_url_spectrum2 == 'Yes':
        url_tmp1 = get_pubchem_url(query=spectrum_ID1)
        url_tmp2 = get_pubchem_url(query=spectrum_ID2)
        t1 = fig.text(0.40, 0.08, f'PubChem URL for {spectrum_ID1}: {url_tmp1}', fontsize=7)
        t2 = fig.text(0.40, 0.05, f'PubChem URL for {spectrum_ID2}: {url_tmp2}', fontsize=7)
        t1.set_url(url_tmp1)
        t2.set_url(url_tmp2)

    if print_url_spectrum1 == 'Yes' and print_url_spectrum2 == 'No':
        url_tmp1 = get_pubchem_url(query=spectrum_ID1)
        t1 = fig.text(0.40, 0.08, f'PubChem URL for {spectrum_ID1}: {url_tmp1}', fontsize=7)
        t1.set_url(url_tmp1)

    if print_url_spectrum1 == 'No' and print_url_spectrum2 == 'Yes':
        url_tmp2 = get_pubchem_url(query=spectrum_ID2)
        t2 = fig.text(0.40, 0.08, f'PubChem URL for {spectrum_ID2}: {url_tmp2}', fontsize=7)
        t2.set_url(url_tmp2)

    fig.savefig(output_path, format='svg')

    if return_plot == True:
        return fig


def wf_transform(spec_mzs, spec_ints, wf_mz, wf_int):
    spec_ints = np.power(spec_mzs, wf_mz) * np.power(spec_ints, wf_int)
    return(spec_ints)


def LE_transform(intensity, thresh, normalization_method):
    intensity_tmp = normalize(intensity, method=normalization_method)
    if np.sum(intensity_tmp) > 0:
        S = scipy.stats.entropy(intensity_tmp.astype('float'))
        if S > 0 and S < thresh:
            w = (1 + S) / (1 + thresh) 
            intensity = np.power(intensity_tmp, w)
    else:
        intensity = np.zeros(len(intensity))
    return intensity 


def normalize(intensities,method='standard'):
    if np.sum(intensities) > 0:
        if method == 'softmax':
            if np.any(intensities > 700):
                print("Warning: some intensities are too large to exponentiate. Applying standard normalization.")
                intensities /= np.sum(intensities)
            else:
                intensities2 = np.exp(intensities)
                if np.isinf(intensities2).sum() == 0:
                    intensities = intensities / np.sum(intensities2)
        elif method == 'standard':
            intensities /= np.sum(intensities)
    return(intensities)


def filter_spec_lcms(spec, mz_min = 0, mz_max = 999999999999, int_min = 0, int_max = 999999999999, is_matched = False):
    if is_matched == False:
        spec = spec[spec[:,0] >= mz_min]
        spec = spec[spec[:,0] <= mz_max]
        spec = spec[spec[:,1] >= int_min]
        spec = spec[spec[:,1] <= int_max]
    else:
        spec = spec[spec[:,0] >= mz_min]
        spec = spec[spec[:,0] <= mz_max]
        spec[spec[:,1] >= int_min] = 0
        spec[spec[:,1] <= int_max] = 0
    return(spec)


def filter_spec_gcms(spec, mz_min = 0, mz_max = 999999999999, int_min = 0, int_max = 999999999999):
    spec[np.where(spec[:,0] < mz_min)[0],1] = 0
    spec[np.where(spec[:,0] > mz_max)[0],1] = 0
    spec[np.where(spec[:,1] < int_min)[0],1] = 0
    spec[np.where(spec[:,1] > int_max)[0],1] = 0
    return(spec)


def remove_noise(spec, nr):
    if spec.shape[0] > 1:
        if nr is not None:
            spec[np.where(spec[:,1] < np.max(spec[:,1]) * nr)[0]] = 0

    return(spec)


def centroid_spectrum(spec, window_size):
    spec = spec[np.argsort(spec[:,0])]

    mz_array = spec[:, 0]
    need_centroid = 0
    if mz_array.shape[0] > 1:
        mz_delta = mz_array[1:] - mz_array[:-1]
        if np.min(mz_delta) <= window_size:
            need_centroid = 1

    if need_centroid:
        intensity_order = np.argsort(-spec[:, 1])
        spec_new = []
        for i in intensity_order:
            mz_delta_allowed = window_size

            if spec[i, 1] > 0:
                i_left = i - 1
                while i_left >= 0:
                    mz_delta_left = spec[i, 0] - spec[i_left, 0]
                    if mz_delta_left <= mz_delta_allowed:
                        i_left -= 1
                    else:
                        break
                i_left += 1

                i_right = i + 1
                while i_right < spec.shape[0]:
                    mz_delta_right = spec[i_right, 0] - spec[i, 0]
                    if mz_delta_right <= mz_delta_allowed:
                        i_right += 1
                    else:
                        break

                intensity_sum = np.sum(spec[i_left:i_right, 1])
                intensity_weighted_sum = np.sum(spec[i_left:i_right, 0] * spec[i_left:i_right, 1])

                spec_new.append([intensity_weighted_sum / intensity_sum, intensity_sum])
                spec[i_left:i_right, 1] = 0

        spec_new = np.array(spec_new)
        spec_new = spec_new[np.argsort(spec_new[:, 0])]
        if spec_new.shape[0] > 1:
            spec_new = spec_new[np.argsort(spec_new[:, 0])]
            return spec_new
        else:
            return np.array([[0,0]])
    else:
        return spec



def match_peaks_in_spectra(spec_a, spec_b, window_size):
    a = 0
    b = 0

    spec_merged = []
    peak_b_int = 0.
    while a < spec_a.shape[0] and b < spec_b.shape[0]:
        mass_delta = spec_a[a, 0] - spec_b[b, 0]
        
        if mass_delta < -window_size:
            spec_merged.append([spec_a[a, 0], spec_a[a, 1], peak_b_int])
            peak_b_int = 0.
            a += 1
        elif mass_delta > window_size:
            spec_merged.append([spec_b[b, 0], 0., spec_b[b, 1]])
            b += 1
        else:
            peak_b_int += spec_b[b, 1]
            b += 1

    if peak_b_int > 0.:
        spec_merged.append([spec_a[a, 0], spec_a[a, 1], peak_b_int])
        peak_b_int = 0.
        a += 1

    if b < spec_b.shape[0]:
        spec_merged += [[x[0], 0., x[1]] for x in spec_b[b:]]

    if a < spec_a.shape[0]:
        spec_merged += [[x[0], x[1], 0.] for x in spec_a[a:]]

    if spec_merged:
        spec_merged = np.array(spec_merged, dtype=np.float64)
    else:
        spec_merged = np.array([[0., 0., 0.]], dtype=np.float64)
    return spec_merged



def convert_spec(spec, mzs):
    ints_tmp = []
    for i in range(0,len(mzs)):
        if mzs[i] in spec[:,0]:
            int_tmp = spec[np.where(spec[:,0] == mzs[i])[0][0],1]
        else:
            int_tmp = 0
        ints_tmp.append(int_tmp)
    out = np.transpose(np.array([mzs,ints_tmp]))
    return out


def get_reference_df(reference_data, likely_reference_IDs=None):
    extension = reference_data.rsplit('.',1)
    extension = extension[(len(extension)-1)]
    if extension == 'mgf' or extension == 'MGF' or extension == 'mzML' or extension == 'mzml' or extension == 'MZML' or extension == 'cdf' or extension == 'CDF' or extension == 'msp' or extension == 'MSP' or extension == 'json' or extension == 'JSON':
        output_path_tmp = reference_data[:-3] + 'txt'
        build_library_from_raw_data(input_path=reference_data, output_path=output_path_tmp, is_reference=True)
        df_reference = pd.read_csv(output_path_tmp, sep='\t')
    if extension == 'txt' or extension == 'TXT':
        df_reference = pd.read_csv(reference_data, sep='\t')
    if likely_reference_IDs is not None:
        likely_reference_IDs = pd.read_csv(likely_reference_IDs, header=None)
        df_reference = df_reference.loc[df_reference.iloc[:,0].isin(likely_reference_IDs.iloc[:,0].tolist())]
    return df_reference



def S_cos(ints_a, ints_b):
    if np.sum(ints_a) == 0 or np.sum(ints_b) == 0:
        return(0)
    else:
        return np.dot(ints_a,ints_b) / (np.sqrt(sum(np.power(ints_a,2))) * np.sqrt(sum(np.power(ints_b,2))))


def ent_renyi(ints, q):
    return np.log(sum(np.power(ints,q))) / (1-q)


def ent_tsallis(ints, q):
    return (sum(np.power(ints,q))-1) / (1-q)


def S_shannon(ints_a, ints_b):
    ent_a = scipy.stats.entropy(ints_a)
    ent_b = scipy.stats.entropy(ints_b)
    ent_ab = scipy.stats.entropy(ints_a + ints_b)
    return(1 - (2 * ent_ab - ent_a - ent_b)/np.log(4))


def S_renyi(ints_a, ints_b, q):
    if q == 1:
        print('Warning: the Renyi Entropy Similarity Measure is equivalent to the Shannon Entropy Similarity Measure when the entropy dimension is 1')
        return S_shannon(ints_a, ints_b)
    else:
        ent_a = ent_renyi(ints_a, q)
        ent_b = ent_renyi(ints_b, q)
        ent_merg = ent_renyi(ints_a/2 + ints_b/2, q)
        N = (1/(1-q)) * (2*np.log(np.sum(np.power(ints_a/2,q))+np.sum(np.power(ints_b/2,q))) - np.log(np.sum(np.power(ints_a,q))) - np.log(np.sum(np.power(ints_b,q))))
        return 1 - (2 * ent_merg - ent_a - ent_b) / N


def S_tsallis(ints_a, ints_b, q):
    if q == 1:
        print('Warning: the Tsallis Entropy Similarity Measure is equivalent to the Shannon Entropy Similarity Measure when the entropy dimension is 1')
        return S_shannon(ints_a, ints_b)
    else:
        ent_a = ent_tsallis(ints_a, q)
        ent_b = ent_tsallis(ints_b, q)
        ent_merg = ent_tsallis(ints_a/2 + ints_b/2, q)
        N = np.sum(2*np.power(ints_a/2,q)+2*np.power(ints_b/2,q)-np.power(ints_a,q)-np.power(ints_b,q)) / (1-q)
        return 1 - (2 * ent_merg - ent_a - ent_b) / N

def S_mixture(ints_a, ints_b, weights={'Cosine':0.25, 'Shannon':0.25, 'Renyi':0.25, 'Tsallis':0.25}, q=1.1):
    if set(weights.keys()).issubset(set(['Cosine','Shannon','Renyi','Tsallis'])) is False:
        print('Error: the keys to the weight parameter dict of the function S_mixture must be one of the four: Cosine, Shannon, Renyi, Tsallis')
        sys.exit()

    similarity = 0
    for key, value in weights.items():
        if key == 'Cosine':
            similarity += value * S_cos(ints_a,ints_b)
        if key == 'Shannon':
            similarity += value * S_shannon(ints_a,ints_b)
        if key == 'Renyi':
            similarity += value * S_renyi(ints_a,ints_b,q)
        if key == 'Tsallis':
            similarity += value * S_tsallis(ints_a,ints_b,q)
    return similarity


def get_contingency_entries(ints_a, ints_b):
    a = 0
    b = 0
    c = 0

    for x, y in zip(ints_a, ints_b):
        if x != 0 and y != 0:
            c += 1
        elif x != 0 and y == 0:
            a += 1
        elif x == 0  and y != 0:
            b += 1
    return [a,b,c]


def S_jaccard(ints_a, ints_b):
    tmp = get_contingency_entries(ints_a, ints_b)
    a = tmp[0]
    b = tmp[1]
    c = tmp[2]
    denom = a + b + c
    if denom == 0:
        similarity = 0
    else:
        similarity = c / (a + b + c)
    return similarity


def S_dice(ints_a, ints_b):
    tmp = get_contingency_entries(ints_a, ints_b)
    a = tmp[0]
    b = tmp[1]
    c = tmp[2]
    denom = a + b + 2 * c
    if denom == 0:
        similarity = 0
    else:
        similarity = 2 * c / denom
    return similarity


def S_3w_jaccard(ints_a, ints_b):
    tmp = get_contingency_entries(ints_a, ints_b)
    a = tmp[0]
    b = tmp[1]
    c = tmp[2]
    denom = a + b + 3 * c
    if denom == 0:
        similarity = 0
    else:
        similarity = 3 * c / denom
    return similarity


def S_sokal_sneath(ints_a, ints_b):
    tmp = get_contingency_entries(ints_a, ints_b)
    a = tmp[0]
    b = tmp[1]
    c = tmp[2]
    denom = 2 * a + 2 * b + c
    if denom == 0:
        similarity = 0
    else:
        similarity = c / denom
    return similarity


def S_binary_cosine(ints_a, ints_b):
    tmp = get_contingency_entries(ints_a, ints_b)
    a = tmp[0]
    b = tmp[1]
    c = tmp[2]
    denom = np.sqrt((a + c) * (b + c))
    if denom == 0:
        similarity = 0
    else:
        similarity = c / denom
    return similarity


def S_mountford(ints_a, ints_b):
    tmp = get_contingency_entries(ints_a, ints_b)
    a = tmp[0]
    b = tmp[1]
    c = tmp[2]
    denom = c * (a + b) + 2 * a * b
    if denom == 0:
        similarity = 1
    else:
        similarity = 2 * c / denom
    return similarity


def S_mcconnaughey(ints_a, ints_b):
    tmp = get_contingency_entries(ints_a, ints_b)
    a = tmp[0]
    b = tmp[1]
    c = tmp[2]
    denom = (a + c) * (b + c)
    if denom == 0:
        similarity = 0
    else:
        similarity = (c**2 - a * b) / denom
    return similarity


def S_driver_kroeber(ints_a, ints_b):
    tmp = get_contingency_entries(ints_a, ints_b)
    a = tmp[0]
    b = tmp[1]
    c = tmp[2]
    denom = 2 * (a + c) * (b + c)
    if denom == 0:
        similarity = 0
    else:
        similarity = c * (a + b + 2 * c) / denom
    return similarity


def S_simpson(ints_a, ints_b):
    tmp = get_contingency_entries(ints_a, ints_b)
    a = tmp[0]
    b = tmp[1]
    c = tmp[2]
    denom = min(a + c, b + c)
    if denom == 0:
        similarity = 0
    else:
        similarity = c / denom
    return similarity


def S_braun_banquet(ints_a, ints_b):
    tmp = get_contingency_entries(ints_a, ints_b)
    a = tmp[0]
    b = tmp[1]
    c = tmp[2]
    denom = max(a + c, b + c)
    if denom == 0:
        similarity = 0
    else:
        similarity = c / denom
    return similarity


def S_fager_mcgowan(ints_a, ints_b):
    tmp = get_contingency_entries(ints_a, ints_b)
    a = tmp[0]
    b = tmp[1]
    c = tmp[2]
    denom1 = np.sqrt((a + c) * (b + c))
    denom2 = 2 * np.sqrt(max(a + c, b + c))
    if denom1 == 0 or denom2 == 0:
        similarity = 0
    else:
        similarity = c / denom1 - 1 / denom2
    return similarity


def S_kulczynski(ints_a, ints_b):
    tmp = get_contingency_entries(ints_a, ints_b)
    a = tmp[0]
    b = tmp[1]
    c = tmp[2]
    denom = a + b
    if denom == 0:
        similarity = 1
    else:
        similarity = c / denom
    return similarity


def S_intersection(ints_a, ints_b):
    tmp = get_contingency_entries(ints_a, ints_b)
    c = tmp[2]
    return c


def S_hamming(ints_a, ints_b):
    tmp = get_contingency_entries(ints_a, ints_b)
    a = tmp[0]
    b = tmp[1]
    denom = a + b
    if denom == 0:
        similarity = 1
    else:
        similarity = 1 / denom
    return similarity


def S_hellinger(ints_a, ints_b):
    tmp = get_contingency_entries(ints_a, ints_b)
    a = tmp[0]
    b = tmp[1]
    c = tmp[2]
    similarity = 1 - np.sqrt((1 - c / np.sqrt((a + c) * (b + c))))
    return similarity


def get_similarity(similarity_measure, q_ints, r_ints, weights, q):

    if similarity_measure == 'cosine':
        similarity = S_cos(q_ints, r_ints)

    elif similarity_measure in ['shannon', 'renyi', 'tsallis']:
            q_ints = normalize(q_ints, method = 'standard')
            r_ints = normalize(r_ints, method = 'standard')
            if similarity_measure == 'shannon':
                similarity = S_shannon(q_ints, r_ints)
            elif similarity_measure == 'renyi':
                similarity = S_renyi(q_ints, r_ints, q)
            elif similarity_measure == 'tsallis':
                similarity = S_tsallis(q_ints, r_ints, q)

    elif similarity_measure == 'mixture':
        similarity = S_mixture(q_ints, r_ints, weights, q)

    elif similarity_measure == 'jaccard':
        similarity = S_jaccard(q_ints, r_ints)

    elif similarity_measure == 'dice':
        similarity = S_dice(q_ints, r_ints)

    elif similarity_measure == '3w_jaccard':
        similarity = S_3w_jaccard(q_ints, r_ints)

    elif similarity_measure == 'sokal_sneath':
        similarity = S_sokal_sneath(q_ints, r_ints)

    elif similarity_measure == 'binary_cosine':
        similarity = S_binary_cosine(q_ints, r_ints)

    elif similarity_measure == 'mountford':
        similarity = S_mountford(q_ints, r_ints)

    elif similarity_measure == 'mcconnaughey':
        similarity = S_mcconnaughey(q_ints, r_ints)

    elif similarity_measure == 'driver_kroeber':
        similarity = S_driver_kroeber(q_ints, r_ints)

    elif similarity_measure == 'simpson':
        similarity = S_simpson(q_ints, r_ints)

    elif similarity_measure == 'braun_banquet':
        similarity = S_braun_banquet(q_ints, r_ints)

    elif similarity_measure == 'fager_mcgowan':
        similarity = S_fager_mcgowan(q_ints, r_ints)

    elif similarity_measure == 'kulczynski':
        similarity = S_kulczynski(q_ints, r_ints)

    elif similarity_measure == 'intersection':
        similarity = S_intersection(q_ints, r_ints)

    elif similarity_measure == 'hamming':
        similarity = S_hamming(q_ints, r_ints)

    elif similarity_measure == 'hellinger':
        similarity = S_hellinger(q_ints, r_ints)

    return similarity


def _vector_to_full_params(X, default_params, optimize_params):
    params = default_params.copy()
    for name, val in zip(optimize_params, X):
        params[name] = float(val)
    return params


def objective_function_HRMS(X, ctx):
    p = _vector_to_full_params(X, ctx["default_params"], ctx["optimize_params"])
    acc = get_acc_HRMS(
        ctx["df_query"], ctx["df_reference"],
        ctx["precursor_ion_mz_tolerance"], ctx["ionization_mode"], ctx["adduct"],
        ctx["similarity_measure"], ctx["weights"], ctx["spectrum_preprocessing_order"],
        ctx["mz_min"], ctx["mz_max"], ctx["int_min"], ctx["int_max"],
        p["window_size_centroiding"], p["window_size_matching"], p["noise_threshold"],
        p["wf_mz"], p["wf_int"], p["LET_threshold"],
        p["entropy_dimension"],
        ctx["high_quality_reference_library"],
        verbose=False
    )
    print(f"\nparams({ctx['optimize_params']}) = {np.array(X)}\naccuracy: {acc*100}%")
    return 1.0 - acc

def objective_function_NRMS(X, ctx):
    p = _vector_to_full_params(X, ctx["default_params"], ctx["optimize_params"])
    acc = get_acc_NRMS(
        ctx["df_query"], ctx["df_reference"],
        ctx["unique_query_ids"], ctx["unique_reference_ids"],
        ctx["similarity_measure"], ctx["weights"], ctx["spectrum_preprocessing_order"],
        ctx["mz_min"], ctx["mz_max"], ctx["int_min"], ctx["int_max"],
        p["noise_threshold"], p["wf_mz"], p["wf_int"], p["LET_threshold"], p["entropy_dimension"],
        ctx["high_quality_reference_library"],
        verbose=False
    )
    print(f"\nparams({ctx['optimize_params']}) = {np.array(X)}\naccuracy: {acc*100}%")
    return 1.0 - acc



def tune_params_DE(query_data=None, reference_data=None, precursor_ion_mz_tolerance=None, ionization_mode=None, adduct=None, chromatography_platform='HRMS', similarity_measure='cosine', weights=None, spectrum_preprocessing_order='CNMWL', mz_min=0, mz_max=999999999, int_min=0, int_max=999999999, high_quality_reference_library=False, optimize_params=["window_size_centroiding","window_size_matching","noise_threshold","wf_mz","wf_int","LET_threshold","entropy_dimension"], param_bounds={"window_size_centroiding":(0.0,0.5),"window_size_matching":(0.0,0.5),"noise_threshold":(0.0,0.25),"wf_mz":(0.0,5.0),"wf_int":(0.0,5.0),"LET_threshold":(0.0,5.0),"entropy_dimension":(1.0,3.0)}, default_params={"window_size_centroiding": 0.5, "window_size_matching":0.5, "noise_threshold":0.10, "wf_mz":0.0, "wf_int":1.0, "LET_threshold":0.0, "entropy_dimension":1.1}, maxiters=3, de_workers=1):

    if query_data is None:
        print('\nError: No argument passed to the mandatory query_data. Please pass the path to the TXT file of the query data.')
        sys.exit()
    else:
        extension = query_data.rsplit('.',1)
        extension = extension[(len(extension)-1)]
        if extension == 'mgf' or extension == 'MGF' or extension == 'mzML' or extension == 'mzml' or extension == 'MZML' or extension == 'cdf' or extension == 'CDF' or extension == 'msp' or extension == 'MSP' or extension == 'json' or extension == 'JSON':
            output_path_tmp = query_data[:-3] + 'txt'
            build_library_from_raw_data(input_path=query_data, output_path=output_path_tmp, is_reference=False)
            df_query = pd.read_csv(output_path_tmp, sep='\t')
        if extension == 'txt' or extension == 'TXT':
            df_query = pd.read_csv(query_data, sep='\t')
        unique_query_ids = df_query.iloc[:,0].unique()

    if reference_data is None:
        print('\nError: No argument passed to the mandatory reference_data. Please pass the path to the TXT file of the reference data.')
        sys.exit()
    else:
        if isinstance(reference_data,str):
            df_reference = get_reference_df(reference_data=reference_data)
            unique_reference_ids = df_reference.iloc[:,0].unique()
        else:
            dfs = []
            unique_reference_ids = []
            for f in reference_data:
                tmp = get_reference_df(reference_data=f)
                dfs.append(tmp)
                unique_reference_ids.extend(tmp.iloc[:,0].unique())
            df_reference = pd.concat(dfs, axis=0, ignore_index=True)

    if 'ionization_mode' in df_reference.columns.tolist() and ionization_mode != None and ionization_mode != 'N/A':
        df_reference = df_reference.loc[df_reference['ionization_mode']==ionization_mode]
    if 'adduct' in df_reference.columns.tolist() and adduct != None and adduct != 'N/A':
        df_reference = df_reference.loc[df_reference['adduct']==adduct]

    unique_query_ids = df_query['id'].unique().tolist()
    unique_reference_ids = df_reference['id'].unique().tolist()

    ctx = dict(
        df_query=df_query,
        df_reference=df_reference,
        precursor_ion_mz_tolerance=precursor_ion_mz_tolerance,
        ionization_mode=ionization_mode,
        adduct=adduct,
        similarity_measure=similarity_measure,
        weights=weights,
        spectrum_preprocessing_order=spectrum_preprocessing_order,
        mz_min=mz_min, mz_max=mz_max, int_min=int_min, int_max=int_max,
        high_quality_reference_library=high_quality_reference_library,
        default_params=default_params,
        optimize_params=optimize_params,
    )

    bounds = [param_bounds[p] for p in optimize_params]

    if chromatography_platform == 'HRMS':
        result = differential_evolution(objective_function_HRMS, bounds=bounds, args=(ctx,), maxiter=maxiters, tol=0.0, workers=de_workers, seed=1)
    else:
        result = differential_evolution(objective_function_NRMS, bounds=bounds, args=(ctx,), maxiter=maxiters, tol=0.0, workers=de_workers, seed=1)

    best_full_params = _vector_to_full_params(result.x, default_params, optimize_params)
    best_acc = 100.0 - (result.fun * 100.0)

    print("\n=== Differential Evolution Result ===")
    print(f"Optimized over: {optimize_params}")
    print("Best values (selected params):")
    for name in optimize_params:
        print(f"  {name}: {best_full_params[name]}")
    print("\nFull parameter set used in final evaluation:")
    for k, v in best_full_params.items():
        print(f"  {k}: {v}")
    print(f"\nBest accuracy: {best_acc:.3f}%")
    _log(f"best = {result.x}, acc={100*(1-result.fun):.3f}%")


default_HRMS_grid = {'similarity_measure':['cosine'], 'weight':[{'Cosine':0.25,'Shannon':0.25,'Renyi':0.25,'Tsallis':0.25}], 'spectrum_preprocessing_order':['FCNMWL'], 'mz_min':[0], 'mz_max':[9999999], 'int_min':[0], 'int_max':[99999999], 'window_size_centroiding':[0.5], 'window_size_matching':[0.5], 'noise_threshold':[0.0], 'wf_mz':[0.0], 'wf_int':[1.0], 'LET_threshold':[0.0], 'entropy_dimension':[1.1], 'high_quality_reference_library':[False]}
default_NRMS_grid = {'similarity_measure':['cosine'], 'weight':[{'Cosine':0.25,'Shannon':0.25,'Renyi':0.25,'Tsallis':0.25}], 'spectrum_preprocessing_order':['FCNMWL'], 'mz_min':[0], 'mz_max':[9999999], 'int_min':[0], 'int_max':[99999999], 'noise_threshold':[0.0], 'wf_mz':[0.0], 'wf_int':[1.0], 'LET_threshold':[0.0], 'entropy_dimension':[1.1], 'high_quality_reference_library':[False]}


def _eval_one_HRMS(df_query, df_reference,
              precursor_ion_mz_tolerance_tmp, ionization_mode_tmp, adduct_tmp,
              similarity_measure_tmp, weight,
              spectrum_preprocessing_order_tmp, mz_min_tmp, mz_max_tmp,
              int_min_tmp, int_max_tmp, noise_threshold_tmp,
              window_size_centroiding_tmp, window_size_matching_tmp,
              wf_mz_tmp, wf_int_tmp, LET_threshold_tmp,
              entropy_dimension_tmp, high_quality_reference_library_tmp):

    acc = get_acc_HRMS(
        df_query=df_query, df_reference=df_reference,
        precursor_ion_mz_tolerance=precursor_ion_mz_tolerance_tmp,
        ionization_mode=ionization_mode_tmp, adduct=adduct_tmp,
        similarity_measure=similarity_measure_tmp, weights=weight,
        spectrum_preprocessing_order=spectrum_preprocessing_order_tmp,
        mz_min=mz_min_tmp, mz_max=mz_max_tmp,
        int_min=int_min_tmp, int_max=int_max_tmp,
        window_size_centroiding=window_size_centroiding_tmp,
        window_size_matching=window_size_matching_tmp,
        noise_threshold=noise_threshold_tmp,
        wf_mz=wf_mz_tmp, wf_int=wf_int_tmp,
        LET_threshold=LET_threshold_tmp,
        entropy_dimension=entropy_dimension_tmp,
        high_quality_reference_library=high_quality_reference_library_tmp,
        verbose=False
    )

    return (
        acc, similarity_measure_tmp, json.dumps(weight), spectrum_preprocessing_order_tmp,
        mz_min_tmp, mz_max_tmp, int_min_tmp, int_max_tmp,
        noise_threshold_tmp, window_size_centroiding_tmp, window_size_matching_tmp,
        wf_mz_tmp, wf_int_tmp, LET_threshold_tmp, entropy_dimension_tmp,
        high_quality_reference_library_tmp
    )


def _eval_one_NRMS(df_query, df_reference, unique_query_ids, unique_reference_ids,
              similarity_measure_tmp, weight,
              spectrum_preprocessing_order_tmp, mz_min_tmp, mz_max_tmp,
              int_min_tmp, int_max_tmp, noise_threshold_tmp,
              wf_mz_tmp, wf_int_tmp, LET_threshold_tmp,
              entropy_dimension_tmp, high_quality_reference_library_tmp):

    acc = get_acc_NRMS(
        df_query=df_query, df_reference=df_reference,
        unique_query_ids=unique_query_ids, unique_reference_ids=unique_reference_ids,
        similarity_measure=similarity_measure_tmp, weights=weight,
        spectrum_preprocessing_order=spectrum_preprocessing_order_tmp,
        mz_min=mz_min_tmp, mz_max=mz_max_tmp,
        int_min=int_min_tmp, int_max=int_max_tmp,
        noise_threshold=noise_threshold_tmp,
        wf_mz=wf_mz_tmp, wf_int=wf_int_tmp,
        LET_threshold=LET_threshold_tmp,
        entropy_dimension=entropy_dimension_tmp,
        high_quality_reference_library=high_quality_reference_library_tmp,
    )

    return (
        acc, similarity_measure_tmp, json.dumps(weight), spectrum_preprocessing_order_tmp,
        mz_min_tmp, mz_max_tmp, int_min_tmp, int_max_tmp, noise_threshold_tmp, 
        wf_mz_tmp, wf_int_tmp, LET_threshold_tmp, entropy_dimension_tmp, high_quality_reference_library_tmp
    )




def tune_params_on_HRMS_data_grid_shiny(query_data=None, reference_data=None, precursor_ion_mz_tolerance=None, ionization_mode=None, adduct=None, grid=None, output_path=None, return_output=False):
    local_grid = {**default_HRMS_grid, **(grid or {})}
    for key, value in local_grid.items():
        globals()[key] = value

    if query_data is None:
        print('\nError: No argument passed to the mandatory query_data. Please pass the path to the data file.')
        sys.exit()
    else:
        extension = query_data.rsplit('.', 1)[-1]
        if extension in ('mgf','MGF','mzML','mzml','MZML','cdf','CDF'):
            output_path_tmp = query_data[:-3] + 'txt'
            build_library_from_raw_data(input_path=query_data, output_path=output_path_tmp, is_reference=False)
            df_query = pd.read_csv(output_path_tmp, sep='\t')
        elif extension in ('txt','TXT'):
            df_query = pd.read_csv(query_data, sep='\t')
        else:
            print(f'\nError: Unsupported query_data extension: {extension}')
            sys.exit()
        unique_query_ids = df_query.iloc[:, 0].unique()

    if reference_data is None:
        print('\nError: No argument passed to the mandatory reference_data. Please pass the path to the data file(s).')
        sys.exit()
    else:
        if isinstance(reference_data, str):
            df_reference = get_reference_df(reference_data=reference_data)
            unique_reference_ids = df_reference.iloc[:, 0].unique()
        else:
            dfs = []
            unique_reference_ids = []
            for f in reference_data:
                tmp = get_reference_df(reference_data=f)
                dfs.append(tmp)
                unique_reference_ids.extend(tmp.iloc[:, 0].unique())
            df_reference = pd.concat(dfs, axis=0, ignore_index=True)

    print(f'\nNote that there are {len(unique_query_ids)} unique query spectra, '
          f'{len(unique_reference_ids)} unique reference spectra, and '
          f'{len(set(unique_query_ids) & set(unique_reference_ids))} of the query and reference spectra IDs are in common.\n')

    if 'ionization_mode' in df_reference.columns.tolist() and ionization_mode != None and ionization_mode != 'N/A':
        df_reference = df_reference.loc[df_reference['ionization_mode']==ionization_mode]
    if 'adduct' in df_reference.columns.tolist() and adduct != None and adduct != 'N/A':
        df_reference = df_reference.loc[df_reference['adduct']==adduct]

    if output_path is None:
        output_path = f'{Path.cwd()}/tuning_param_output.txt'
        print(f'Warning: since output_path=None, the output will be written to the current working directory: {output_path}')

    param_grid = product(
        similarity_measure, weight, spectrum_preprocessing_order, mz_min, mz_max, int_min, int_max,
        noise_threshold, window_size_centroiding, window_size_matching, wf_mz, wf_int, LET_threshold,
        entropy_dimension, high_quality_reference_library
    )

    results = []
    total = (
        len(similarity_measure) * len(weight) * len(spectrum_preprocessing_order) * len(mz_min) * len(mz_max) *
        len(int_min) * len(int_max) * len(noise_threshold) * len(window_size_centroiding) *
        len(window_size_matching) * len(wf_mz) * len(wf_int) * len(LET_threshold) *
        len(entropy_dimension) * len(high_quality_reference_library)
    )
    done = 0
    for params in param_grid:
        res = _eval_one_HRMS(df_query, df_reference, precursor_ion_mz_tolerance, ionization_mode, adduct, *params)
        results.append(res)
        done += 1
        print(f'Completed {done}/{total} grid combinations.\n', flush=True)

    df_out = pd.DataFrame(results, columns=[
        'ACC','SIMILARITY.MEASURE','WEIGHT','SPECTRUM.PROCESSING.ORDER','MZ.MIN','MZ.MAX',
        'INT.MIN','INT.MAX','NOISE.THRESHOLD','WINDOW.SIZE.CENTROIDING','WINDOW.SIZE.MATCHING',
        'WF.MZ','WF.INT','LET.THRESHOLD','ENTROPY.DIMENSION','HIGH.QUALITY.REFERENCE.LIBRARY'
    ])

    if 'WEIGHT' in df_out.columns:
        df_out['WEIGHT'] = (
            df_out['WEIGHT'].astype(str)
                .str.replace("\"","",regex=False)
                .str.replace("{","",regex=False)
                .str.replace("}","",regex=False)
                .str.replace(":","",regex=False)
                .str.replace("Cosine","",regex=False)
                .str.replace("Shannon","",regex=False)
                .str.replace("Renyi","",regex=False)
                .str.replace("Tsallis","",regex=False)
                .str.replace(" ","",regex=False)
        )

    if return_output:
        return df_out
    else:
        df_out.to_csv(output_path, index=False, sep='\t', quoting=csv.QUOTE_NONE)
        print(f'Wrote results to {output_path}')



def tune_params_on_NRMS_data_grid(query_data=None, reference_data=None, grid=None, output_path=None, return_output=False):
    grid = {**default_NRMS_grid, **(grid or {})}
    for key, value in grid.items():
        globals()[key] = value

    if query_data is None:
        print('\nError: No argument passed to the mandatory query_data. Please pass the path to the TXT file of the query data.')
        sys.exit()
    else:
        extension = query_data.rsplit('.',1)
        extension = extension[(len(extension)-1)]
        if extension == 'mgf' or extension == 'MGF' or extension == 'mzML' or extension == 'mzml' or extension == 'MZML' or extension == 'cdf' or extension == 'CDF' or extension == 'msp' or extension == 'MSP' or extension == 'json' or extension == 'JSON':
            output_path_tmp = query_data[:-3] + 'txt'
            build_library_from_raw_data(input_path=query_data, output_path=output_path_tmp, is_reference=False)
            df_query = pd.read_csv(output_path_tmp, sep='\t')
        if extension == 'txt' or extension == 'TXT':
            df_query = pd.read_csv(query_data, sep='\t')
        unique_query_ids = df_query.iloc[:,0].unique()

    if reference_data is None:
        print('\nError: No argument passed to the mandatory reference_data. Please pass the path to the TXT file of the reference data.')
        sys.exit()
    else:
        if isinstance(reference_data,str):
            df_reference = get_reference_df(reference_data=reference_data)
            unique_reference_ids = df_reference.iloc[:,0].unique()
        else:
            dfs = []
            unique_reference_ids = []
            for f in reference_data:
                tmp = get_reference_df(reference_data=f)
                dfs.append(tmp)
                unique_reference_ids.extend(tmp.iloc[:,0].unique())
            df_reference = pd.concat(dfs, axis=0, ignore_index=True)

    print(f'\nNote that there are {len(unique_query_ids)} unique query spectra, {len(unique_reference_ids)} unique reference spectra, and {len(set(unique_query_ids) & set(unique_reference_ids))} of the query and reference spectra IDs are in common.\n')

    if output_path is None:
        output_path = f'{Path.cwd()}/tuning_param_output.txt'
        print(f'Warning: since output_path=None, the output will be written to the current working directory: {output_path}')

    param_grid = product(similarity_measure, weight, spectrum_preprocessing_order, mz_min, mz_max, int_min, int_max,
                         noise_threshold, wf_mz, wf_int, LET_threshold, entropy_dimension, high_quality_reference_library)
    results = Parallel(n_jobs=-1, verbose=10)(delayed(_eval_one_NRMS)(df_query, df_reference, unique_query_ids, unique_reference_ids, *params) for params in param_grid)

    df_out = pd.DataFrame(results, columns=[
        'ACC','SIMILARITY.MEASURE','WEIGHT','SPECTRUM.PROCESSING.ORDER', 'MZ.MIN','MZ.MAX','INT.MIN','INT.MAX',
        'NOISE.THRESHOLD','WF.MZ','WF.INT','LET.THRESHOLD','ENTROPY.DIMENSION', 'HIGH.QUALITY.REFERENCE.LIBRARY'
    ])
    df_out['WEIGHT'] = df_out['WEIGHT'].str.replace("\"","",regex=False)
    df_out['WEIGHT'] = df_out['WEIGHT'].str.replace("{","",regex=False)
    df_out['WEIGHT'] = df_out['WEIGHT'].str.replace("}","",regex=False)
    df_out['WEIGHT'] = df_out['WEIGHT'].str.replace(":","",regex=False)
    df_out['WEIGHT'] = df_out['WEIGHT'].str.replace("Cosine","",regex=False)
    df_out['WEIGHT'] = df_out['WEIGHT'].str.replace("Shannon","",regex=False)
    df_out['WEIGHT'] = df_out['WEIGHT'].str.replace("Renyi","",regex=False)
    df_out['WEIGHT'] = df_out['WEIGHT'].str.replace("Tsallis","",regex=False)
    df_out['WEIGHT'] = df_out['WEIGHT'].str.replace(" ","",regex=False)
    if return_output is False:
        df_out.to_csv(output_path, index=False, sep='\t', quoting=csv.QUOTE_NONE)
    else:
        return df_out



def tune_params_on_NRMS_data_grid_shiny(query_data=None, reference_data=None, grid=None, output_path=None, return_output=False):
    local_grid = {**default_NRMS_grid, **(grid or {})}
    for key, value in local_grid.items():
        globals()[key] = value

    if query_data is None:
        print('\nError: No argument passed to the mandatory query_data. Please pass the path to the data file.')
        sys.exit()
    else:
        extension = query_data.rsplit('.', 1)[-1]
        if extension in ('mgf','MGF','mzML','mzml','MZML','cdf','CDF'):
            output_path_tmp = query_data[:-3] + 'txt'
            build_library_from_raw_data(input_path=query_data, output_path=output_path_tmp, is_reference=False)
            df_query = pd.read_csv(output_path_tmp, sep='\t')
        elif extension in ('txt','TXT'):
            df_query = pd.read_csv(query_data, sep='\t')
        else:
            print(f'\nError: Unsupported query_data extension: {extension}')
            sys.exit()
        unique_query_ids = df_query.iloc[:, 0].unique()

    if reference_data is None:
        print('\nError: No argument passed to the mandatory reference_data. Please pass the path to the data file(s).')
        sys.exit()
    else:
        if isinstance(reference_data, str):
            df_reference = get_reference_df(reference_data=reference_data)
            unique_reference_ids = df_reference.iloc[:, 0].unique()
        else:
            dfs = []
            unique_reference_ids = []
            for f in reference_data:
                tmp = get_reference_df(reference_data=f)
                dfs.append(tmp)
                unique_reference_ids.extend(tmp.iloc[:, 0].unique())
            df_reference = pd.concat(dfs, axis=0, ignore_index=True)

    print(f'\nNote that there are {len(unique_query_ids)} unique query spectra, '
          f'{len(unique_reference_ids)} unique reference spectra, and '
          f'{len(set(unique_query_ids) & set(unique_reference_ids))} of the query and reference spectra IDs are in common.\n')

    if output_path is None:
        output_path = f'{Path.cwd()}/tuning_param_output.txt'
        print(f'Warning: since output_path=None, the output will be written to the current working directory: {output_path}')

    param_grid = product(
        similarity_measure, weight, spectrum_preprocessing_order, mz_min, mz_max, int_min, int_max,
        noise_threshold, wf_mz, wf_int, LET_threshold,
        entropy_dimension, high_quality_reference_library
    )

    results = []
    total = (
        len(similarity_measure) * len(weight) * len(spectrum_preprocessing_order) * len(mz_min) * len(mz_max) * len(int_min) *
        len(int_max) * len(noise_threshold) * len(wf_mz) * len(wf_int) * len(LET_threshold) * len(entropy_dimension) * len(high_quality_reference_library)
    )
    done = 0
    for params in param_grid:
        res = _eval_one_NRMS(df_query, df_reference, unique_query_ids, unique_reference_ids, *params)
        results.append(res)
        done += 1
        print(f'Completed {done}/{total} grid combinations.\n', flush=True)

    df_out = pd.DataFrame(results, columns=[
        'ACC','SIMILARITY.MEASURE','WEIGHT','SPECTRUM.PROCESSING.ORDER','MZ.MIN','MZ.MAX',
        'INT.MIN','INT.MAX','NOISE.THRESHOLD','WF.MZ','WF.INT','LET.THRESHOLD','ENTROPY.DIMENSION','HIGH.QUALITY.REFERENCE.LIBRARY'
    ])

    if 'WEIGHT' in df_out.columns:
        df_out['WEIGHT'] = (
            df_out['WEIGHT'].astype(str)
                .str.replace("\"","",regex=False)
                .str.replace("{","",regex=False)
                .str.replace("}","",regex=False)
                .str.replace(":","",regex=False)
                .str.replace("Cosine","",regex=False)
                .str.replace("Shannon","",regex=False)
                .str.replace("Renyi","",regex=False)
                .str.replace("Tsallis","",regex=False)
                .str.replace(" ","",regex=False)
        )

    if return_output:
        return df_out
    else:
        df_out.to_csv(output_path, index=False, sep='\t', quoting=csv.QUOTE_NONE)
        print(f'Wrote results to {output_path}')




def get_acc_HRMS(df_query, df_reference, precursor_ion_mz_tolerance, ionization_mode, adduct, similarity_measure, weights, spectrum_preprocessing_order, mz_min, mz_max, int_min, int_max, window_size_centroiding, window_size_matching, noise_threshold, wf_mz, wf_int, LET_threshold, entropy_dimension, high_quality_reference_library, verbose=True):
    n_top_matches_to_save = 1
    unique_reference_ids = df_reference['id'].dropna().astype(str).unique().tolist()
    unique_query_ids = df_query['id'].dropna().astype(str).unique().tolist()
    all_similarity_rows = []

    for query_idx, qid in enumerate(unique_query_ids):
        if verbose:
            print(f'query spectrum #{query_idx} is being identified')

        q_mask = (df_query['id'] == qid)
        q_idxs = np.where(q_mask)[0]
        if q_idxs.size == 0:
            all_similarity_rows.append([0.0]*len(unique_reference_ids))
            continue

        q_spec_base = np.asarray(pd.concat([df_query['mz_ratio'].iloc[q_idxs], df_query['intensity'].iloc[q_idxs]], axis=1).reset_index(drop=True))

        if 'precursor_ion_mz' in df_query.columns and 'precursor_ion_mz' in df_reference.columns and precursor_ion_mz_tolerance is not None:
            precursor = float(df_query['precursor_ion_mz'].iloc[q_idxs[0]])
            df_reference_tmp = df_reference.loc[df_reference['precursor_ion_mz'].between(precursor - precursor_ion_mz_tolerance, precursor + precursor_ion_mz_tolerance, inclusive='both'), ['id', 'mz_ratio', 'intensity']].copy()
        else:
            df_reference_tmp = df_reference[['id','mz_ratio','intensity']].copy()

        if df_reference_tmp.empty:
            all_similarity_rows.append([0.0]*len(unique_reference_ids))
            continue

        ref_groups = dict(tuple(df_reference_tmp.groupby('id', sort=False)))

        similarity_by_ref = {}

        for ref_id, r_df in ref_groups.items():
            q_spec = q_spec_base.copy()
            r_spec = np.asarray(pd.concat([r_df['mz_ratio'], r_df['intensity']], axis=1).reset_index(drop=True))

            is_matched = False
            for transformation in spectrum_preprocessing_order:
                if np.isinf(q_spec[:, 1]).any():
                    q_spec[:, 1] = 0.0
                if np.isinf(r_spec[:, 1]).any():
                    r_spec[:, 1] = 0.0

                if transformation == 'C' and q_spec.shape[0] > 1 and r_spec.shape[0] > 1:
                    q_spec = centroid_spectrum(q_spec, window_size=window_size_centroiding)
                    r_spec = centroid_spectrum(r_spec, window_size=window_size_centroiding)

                if transformation == 'M' and q_spec.shape[0] > 1 and r_spec.shape[0] > 1:
                    m_spec = match_peaks_in_spectra(
                        spec_a=q_spec, spec_b=r_spec, window_size=window_size_matching
                    )
                    if m_spec.size == 0:
                        q_spec = np.empty((0,2))
                        r_spec = np.empty((0,2))
                    else:
                        q_spec = m_spec[:, 0:2]
                        r_spec = m_spec[:, [0, 2]]
                    is_matched = True

                if transformation == 'W' and q_spec.shape[0] > 1 and r_spec.shape[0] > 1:
                    q_spec[:, 1] = wf_transform(q_spec[:, 0], q_spec[:, 1], wf_mz, wf_int)
                    r_spec[:, 1] = wf_transform(r_spec[:, 0], r_spec[:, 1], wf_mz, wf_int)

                if transformation == 'L' and q_spec.shape[0] > 1 and r_spec.shape[0] > 1:
                    q_spec[:, 1] = LE_transform(q_spec[:, 1], LET_threshold, normalization_method='standard')
                    r_spec[:, 1] = LE_transform(r_spec[:, 1], LET_threshold, normalization_method='standard')

                if transformation == 'N' and q_spec.shape[0] > 1 and r_spec.shape[0] > 1:
                    q_spec = remove_noise(q_spec, nr=noise_threshold)
                    if not high_quality_reference_library:
                        r_spec = remove_noise(r_spec, nr=noise_threshold)

                if transformation == 'F' and q_spec.shape[0] > 1 and r_spec.shape[0] > 1:
                    q_spec = filter_spec_lcms(
                        q_spec, mz_min=mz_min, mz_max=mz_max, int_min=int_min, int_max=int_max, is_matched=is_matched
                    )
                    if not high_quality_reference_library:
                        r_spec = filter_spec_lcms(
                            r_spec, mz_min=mz_min, mz_max=mz_max, int_min=int_min, int_max=int_max, is_matched=is_matched
                        )

            if q_spec.shape[0] > 1 and r_spec.shape[0] > 1:
                q_ints = q_spec[:, 1]
                r_ints = r_spec[:, 1]
                if np.sum(q_ints) != 0 and np.sum(r_ints) != 0:
                    sim = get_similarity(similarity_measure, q_ints, r_ints, weights, entropy_dimension)
                else:
                    sim = 0.0
            else:
                sim = 0.0

            similarity_by_ref[str(ref_id)] = float(sim)

        row = [similarity_by_ref.get(ref_id, 0.0) for ref_id in unique_reference_ids]
        all_similarity_rows.append(row)

    df_scores = pd.DataFrame(all_similarity_rows, index=unique_query_ids, columns=unique_reference_ids)
    df_scores.index.name = 'QUERY.SPECTRUM.ID'

    top_idx = df_scores.values.argmax(axis=1)
    top_scores = df_scores.values[np.arange(df_scores.shape[0]), top_idx]
    top_ids = [df_scores.columns[i] for i in top_idx]

    df_tmp = pd.DataFrame({'TRUE.ID': df_scores.index.to_list(), 'PREDICTED.ID': top_ids, 'SCORE': top_scores})
    if verbose:
        print(df_tmp)

    acc = (df_tmp['TRUE.ID'] == df_tmp['PREDICTED.ID']).mean()
    return acc



def get_acc_NRMS(df_query, df_reference, unique_query_ids, unique_reference_ids, similarity_measure, weights, spectrum_preprocessing_order, mz_min, mz_max, int_min, int_max, noise_threshold, wf_mz, wf_int, LET_threshold, entropy_dimension, high_quality_reference_library, verbose=True):

    n_top_matches_to_save = 1

    min_mz = int(np.min([np.min(df_query.iloc[:,1]), np.min(df_reference.iloc[:,1])]))
    max_mz = int(np.max([np.max(df_query.iloc[:,1]), np.max(df_reference.iloc[:,1])]))
    mzs = np.linspace(min_mz,max_mz,(max_mz-min_mz+1))

    all_similarity_scores =  []
    for query_idx in range(0,len(unique_query_ids)):
        q_idxs_tmp = np.where(df_query.iloc[:,0] == unique_query_ids[query_idx])[0]
        q_spec_tmp = np.asarray(pd.concat([df_query.iloc[q_idxs_tmp,1], df_query.iloc[q_idxs_tmp,2]], axis=1).reset_index(drop=True))
        q_spec_tmp = convert_spec(q_spec_tmp,mzs)

        similarity_scores = []
        for ref_idx in range(0,len(unique_reference_ids)):
            q_spec = q_spec_tmp
            if verbose is True and ref_idx % 1000 == 0:
                print(f'Query spectrum #{query_idx} has had its similarity with {ref_idx} reference library spectra computed')
            r_idxs_tmp = np.where(df_reference.iloc[:,0] == unique_reference_ids[ref_idx])[0]
            r_spec_tmp = np.asarray(pd.concat([df_reference.iloc[r_idxs_tmp,1], df_reference.iloc[r_idxs_tmp,2]], axis=1).reset_index(drop=True))
            r_spec = convert_spec(r_spec_tmp,mzs)

            for transformation in spectrum_preprocessing_order:
                if np.isinf(q_spec[:,1]).sum() > 0:
                    q_spec[:,1] = np.zeros(q_spec.shape[0])
                if np.isinf(r_spec[:,1]).sum() > 0:
                    r_spec[:,1] = np.zeros(r_spec.shape[0])
                if transformation == 'W':
                    q_spec[:,1] = wf_transform(q_spec[:,0], q_spec[:,1], wf_mz, wf_int)
                    r_spec[:,1] = wf_transform(r_spec[:,0], r_spec[:,1], wf_mz, wf_int)
                if transformation == 'L':
                    q_spec[:,1] = LE_transform(q_spec[:,1], LET_threshold, normalization_method='standard')
                    r_spec[:,1] = LE_transform(r_spec[:,1], LET_threshold, normalization_method='standard')
                if transformation == 'N':
                    q_spec = remove_noise(q_spec, nr = noise_threshold)
                    if high_quality_reference_library == False:
                        r_spec = remove_noise(r_spec, nr = noise_threshold)
                if transformation == 'F':
                    q_spec = filter_spec_gcms(q_spec, mz_min = mz_min, mz_max = mz_max, int_min = int_min, int_max = int_max)
                    if high_quality_reference_library == False:
                        r_spec = filter_spec_gcms(r_spec, mz_min = mz_min, mz_max = mz_max, int_min = int_min, int_max = int_max)

            q_ints = q_spec[:,1]
            r_ints = r_spec[:,1]

            if np.sum(q_ints) != 0 and np.sum(r_ints) != 0:
                similarity_score = get_similarity(similarity_measure, q_spec[:,1], r_spec[:,1], weights, entropy_dimension)
            else:
                similarity_score = 0

            similarity_scores.append(similarity_score)
        all_similarity_scores.append(similarity_scores)

    df_scores = pd.DataFrame(all_similarity_scores, columns = unique_reference_ids)
    df_scores.index = unique_query_ids
    df_scores.index.names = ['QUERY.SPECTRUM.ID']

    preds = []
    scores = []
    for i in range(0, df_scores.shape[0]):
        df_scores_tmp = df_scores
        preds_tmp = []
        scores_tmp = []
        for j in range(0, n_top_matches_to_save):
            top_ref_specs_tmp = df_scores_tmp.iloc[i,np.where(df_scores_tmp.iloc[i,:] == np.max(df_scores_tmp.iloc[i,:]))[0]]
            cols_to_keep = np.where(df_scores_tmp.iloc[i,:] != np.max(df_scores_tmp.iloc[i,:]))[0]
            df_scores_tmp = df_scores_tmp.iloc[:,cols_to_keep]

            preds_tmp.append(';'.join(map(str,top_ref_specs_tmp.index.to_list())))
            if len(top_ref_specs_tmp.values) == 0:
                scores_tmp.append(0)
            else:
                scores_tmp.append(top_ref_specs_tmp.values[0])
        preds.append(preds_tmp)
        scores.append(scores_tmp)

    preds = np.array(preds)
    scores = np.array(scores)
    out = np.c_[unique_query_ids,preds,scores]
    df_tmp = pd.DataFrame(out, columns=['TRUE.ID','PREDICTED.ID','SCORE'])
    acc = (df_tmp['TRUE.ID']==df_tmp['PREDICTED.ID']).mean()
    return acc



def run_spec_lib_matching_on_HRMS_data_shiny(query_data=None, reference_data=None, precursor_ion_mz_tolerance=None, ionization_mode=None, adduct=None, likely_reference_ids=None, similarity_measure='cosine', weights={'Cosine':0.25,'Shannon':0.25,'Renyi':0.25,'Tsallis':0.25}, spectrum_preprocessing_order='FCNMWL', high_quality_reference_library=False, mz_min=0, mz_max=9999999, int_min=0, int_max=9999999, window_size_centroiding=0.5, window_size_matching=0.5, noise_threshold=0.0, wf_mz=0.0, wf_intensity=1.0, LET_threshold=0.0, entropy_dimension=1.1, n_top_matches_to_save=1, print_id_results=False, output_identification=None, output_similarity_scores=None, return_ID_output=False, verbose=True):
    if query_data is None:
        print('\nError: No argument passed to the mandatory query_data. Please pass the path to the CSV file of the query data.')
        sys.exit()
    else:
        extension = query_data.rsplit('.',1)
        extension = extension[(len(extension)-1)]
        if extension == 'mgf' or extension == 'MGF' or extension == 'mzML' or extension == 'mzml' or extension == 'MZML' or extension == 'cdf' or extension == 'CDF' or extension == 'json' or extension == 'JSON':
            output_path_tmp = query_data[:-3] + 'txt'
            build_library_from_raw_data(input_path=query_data, output_path=output_path_tmp, is_reference=False)
            #build_library_from_raw_data(input_path=query_data, output_path=output_path_tmp, is_reference=True)
            df_query = pd.read_csv(output_path_tmp, sep='\t')
        if extension == 'txt' or extension == 'TXT':
            df_query = pd.read_csv(query_data, sep='\t')
        unique_query_ids = df_query['id'].unique()

    if reference_data is None:
        print('\nError: No argument passed to the mandatory reference_data. Please pass the path to the reference data.')
        sys.exit()
    else:
        if isinstance(reference_data,str):
            df_reference = get_reference_df(reference_data,likely_reference_ids)
        else:
            dfs = []
            for f in reference_data:
                tmp = get_reference_df(f,likely_reference_ids)
                dfs.append(tmp)
            df_reference = pd.concat(dfs, axis=0, ignore_index=True)

    if 'ionization_mode' in df_reference.columns.tolist() and ionization_mode != 'N/A':
        df_reference = df_reference.loc[df_reference['ionization_mode']==ionization_mode]
    if 'adduct' in df_reference.columns.tolist() and adduct != 'N/A':
        df_reference = df_reference.loc[df_reference['adduct']==adduct]

    if spectrum_preprocessing_order is not None:
        spectrum_preprocessing_order = list(spectrum_preprocessing_order)
    else:
        spectrum_preprocessing_order = ['F', 'C', 'N', 'M', 'W', 'L']
    if 'M' not in spectrum_preprocessing_order:
        print(f'Error: \'M\' must be a character in spectrum_preprocessing_order.')
        sys.exit()
    if 'C' in spectrum_preprocessing_order:
        if spectrum_preprocessing_order.index('C') > spectrum_preprocessing_order.index('M'):
            print(f'Error: \'C\' must come before \'M\' in spectrum_preprocessing_order.')
            sys.exit()
    if set(spectrum_preprocessing_order) - {'F','C','N','M','W','L'}:
        print(f'Error: spectrum_preprocessing_order must contain only \'C\', \'F\', \'M\', \'N\', \'L\', \'W\'.')
        sys.exit()


    if similarity_measure not in ['cosine','shannon','renyi','tsallis','mixture','jaccard','dice','3w_jaccard','sokal_sneath','binary_cosine','mountford','mcconnaughey','driver_kroeber','simpson','braun_banquet','fager_mcgowan','kuldzynski','intersection','hamming','hellinger']:
        print('\nError: similarity_measure must be either cosine, shannon, renyi, tsallis, mixture, jaccard, dice, 3w_jaccard, sokal_sneath, binary_cosine, mountford, mcconnaughey, driver_kroeber, simpson, braun_banquet, fager_mcgowan, kulczynski, intersection, hamming, or hellinger')
        sys.exit()

    if isinstance(int_min,int) is True:
        int_min = float(int_min)
    if isinstance(int_max,int) is True:
        int_max = float(int_max)
    if isinstance(mz_min,int) is False or isinstance(mz_max,int) is False or isinstance(int_min,float) is False or isinstance(int_max,float) is False:
        print('Error: mz_min must be a non-negative integer, mz_max must be a positive integer, int_min must be a non-negative float, and int_max must be a positive float')
        sys.exit()
    if mz_min < 0:
        print('\nError: mz_min should be a non-negative integer')
        sys.exit()
    if mz_max <= 0:
        print('\nError: mz_max should be a positive integer')
        sys.exit()
    if int_min < 0:
        print('\nError: int_min should be a non-negative float')
        sys.exit()
    if int_max <= 0:
        print('\nError: int_max should be a positive float')
        sys.exit()

    if isinstance(window_size_centroiding,float) is False or window_size_centroiding <= 0.0:
        print('Error: window_size_centroiding must be a positive float.')
        sys.exit()
    if isinstance(window_size_matching,float) is False or window_size_matching<= 0.0:
        print('Error: window_size_matching must be a positive float.')
        sys.exit()

    if isinstance(noise_threshold,int) is True:
        noise_threshold = float(noise_threshold)
    if isinstance(noise_threshold,float) is False or noise_threshold < 0:
        print('Error: noise_threshold must be a positive float.')
        sys.exit()

    if isinstance(wf_intensity,int) is True:
        wf_intensity = float(wf_intensity)
    if isinstance(wf_mz,int) is True:
        wf_mz = float(wf_mz)
    if isinstance(wf_intensity,float) is False or isinstance(wf_mz,float) is False:
        print('Error: wf_mz and wf_intensity must be integers or floats')
        sys.exit()

    if entropy_dimension <= 0:
        print('\nError: entropy_dimension should be a positive float')
        sys.exit()
    else:
        q = entropy_dimension

    normalization_method = 'standard'

    if n_top_matches_to_save <= 0 or isinstance(n_top_matches_to_save,int)==False:
        print('\nError: n_top_matches_to_save should be a positive integer')
        sys.exit()

    if isinstance(print_id_results,bool)==False:
        print('\nError: print_id_results must be either True or False')
        sys.exit()
    
    if output_identification is None:
        output_identification = f'{Path.cwd()}/output_identification.txt'
        print(f'Warning: writing identification output to {output_identification}')

    if output_similarity_scores is None:
        output_similarity_scores = f'{Path.cwd()}/output_all_similarity_scores.txt'
        print(f'Warning: writing similarity scores to {output_similarity_scores}')


    unique_reference_ids = df_reference['id'].unique().tolist()
    all_similarity_scores = []

    for query_idx in range(len(unique_query_ids)):
        if verbose:
            print(f'query spectrum #{query_idx} is being identified')

        q_mask = (df_query['id'] == unique_query_ids[query_idx])
        q_idxs_tmp = np.where(q_mask)[0]
        q_spec_tmp = np.asarray(pd.concat([df_query['mz_ratio'].iloc[q_idxs_tmp], df_query['intensity'].iloc[q_idxs_tmp]], axis=1).reset_index(drop=True))

        if 'precursor_ion_mz' in df_query.columns.tolist() and 'precursor_ion_mz' in df_reference.columns.tolist() and precursor_ion_mz_tolerance != None:
            precursor_ion_mz_tmp = df_query['precursor_ion_mz'].iloc[q_idxs_tmp[0]]
            df_reference_tmp = df_reference.loc[df_reference['precursor_ion_mz'].between(precursor_ion_mz_tmp-precursor_ion_mz_tolerance, precursor_ion_mz_tmp+precursor_ion_mz_tolerance, inclusive='both'),['id','mz_ratio','intensity']].copy()
        else:
            df_reference_tmp = df_reference.copy()

        ref_groups = dict(tuple(df_reference_tmp.groupby('id', sort=False)))
        unique_reference_ids_tmp = list(ref_groups.keys())

        similarity_by_ref = {}
        for ref_id in unique_reference_ids_tmp:
            q_spec = q_spec_tmp.copy()
            r_df = ref_groups[ref_id]
            r_spec = np.asarray(pd.concat([r_df['mz_ratio'], r_df['intensity']], axis=1).reset_index(drop=True))

            is_matched = False

            for transformation in spectrum_preprocessing_order:
                if np.isinf(q_spec[:, 1]).sum() > 0:
                    q_spec[:, 1] = np.zeros(q_spec.shape[0])
                if np.isinf(r_spec[:, 1]).sum() > 0:
                    r_spec[:, 1] = np.zeros(r_spec.shape[0])

                if transformation == 'C' and q_spec.shape[0] > 1 and r_spec.shape[0] > 1:
                    q_spec = centroid_spectrum(q_spec, window_size=window_size_centroiding)
                    r_spec = centroid_spectrum(r_spec, window_size=window_size_centroiding)

                if transformation == 'M' and q_spec.shape[0] > 1 and r_spec.shape[0] > 1:
                    m_spec = match_peaks_in_spectra(spec_a=q_spec, spec_b=r_spec, window_size=window_size_matching)
                    q_spec = m_spec[:, 0:2]
                    r_spec = m_spec[:, [0, 2]]
                    is_matched = True

                if transformation == 'W' and q_spec.shape[0] > 1 and r_spec.shape[0] > 1:
                    q_spec[:, 1] = wf_transform(q_spec[:, 0], q_spec[:, 1], wf_mz, wf_intensity)
                    r_spec[:, 1] = wf_transform(r_spec[:, 0], r_spec[:, 1], wf_mz, wf_intensity)

                if transformation == 'L' and q_spec.shape[0] > 1 and r_spec.shape[0] > 1:
                    q_spec[:, 1] = LE_transform(q_spec[:, 1], LET_threshold, normalization_method=normalization_method)
                    r_spec[:, 1] = LE_transform(r_spec[:, 1], LET_threshold, normalization_method=normalization_method)

                if transformation == 'N' and q_spec.shape[0] > 1 and r_spec.shape[0] > 1:
                    q_spec = remove_noise(q_spec, nr=noise_threshold)
                    if not high_quality_reference_library:
                        r_spec = remove_noise(r_spec, nr=noise_threshold)

                if transformation == 'F' and q_spec.shape[0] > 1 and r_spec.shape[0] > 1:
                    q_spec = filter_spec_lcms(
                        q_spec, mz_min=mz_min, mz_max=mz_max, int_min=int_min, int_max=int_max, is_matched=is_matched
                    )
                    if not high_quality_reference_library:
                        r_spec = filter_spec_lcms(
                            r_spec, mz_min=mz_min, mz_max=mz_max, int_min=int_min, int_max=int_max, is_matched=is_matched
                        )

            q_ints = q_spec[:, 1]
            r_ints = r_spec[:, 1]

            if np.sum(q_ints) != 0 and np.sum(r_ints) != 0 and q_spec.shape[0] > 1 and r_spec.shape[0] > 1:
                sim = get_similarity(similarity_measure, q_ints, r_ints, weights, entropy_dimension)
            else:
                sim = 0.0

            similarity_by_ref[ref_id] = sim

        row_scores = [similarity_by_ref.get(ref_id, 0.0) for ref_id in unique_reference_ids]
        all_similarity_scores.append(row_scores)

    df_scores = pd.DataFrame(all_similarity_scores, index=unique_query_ids, columns=unique_reference_ids)
    df_scores.index = unique_query_ids
    df_scores.index.names = ['QUERY.SPECTRUM.ID']


    preds = []
    scores = []
    for i in range(0, df_scores.shape[0]):
        df_scores_tmp = df_scores
        preds_tmp = []
        scores_tmp = []
        for j in range(0, n_top_matches_to_save):
            top_ref_specs_tmp = df_scores_tmp.iloc[i,np.where(df_scores_tmp.iloc[i,:] == np.max(df_scores_tmp.iloc[i,:]))[0]]
            cols_to_keep = np.where(df_scores_tmp.iloc[i,:] != np.max(df_scores_tmp.iloc[i,:]))[0]
            df_scores_tmp = df_scores_tmp.iloc[:,cols_to_keep]

            preds_tmp.append(';'.join(map(str,top_ref_specs_tmp.index.to_list())))
            if len(top_ref_specs_tmp.values) == 0:
                scores_tmp.append(0)
            else:
                scores_tmp.append(top_ref_specs_tmp.values[0])
        preds.append(preds_tmp)
        scores.append(scores_tmp)

    preds = np.array(preds)
    scores = np.array(scores)
    out = np.c_[preds,scores]

    cnames_preds = []
    cnames_scores = []
    for i in range(0,n_top_matches_to_save):
        cnames_preds.append(f'RANK.{i+1}.PRED')
        cnames_scores.append(f'RANK.{i+1}.SIMILARITY.SCORE')

    df_top_ref_specs = pd.DataFrame(out, columns = [*cnames_preds, *cnames_scores])
    df_top_ref_specs.index = unique_query_ids
    df_top_ref_specs.index.names = ['QUERY.SPECTRUM.ID']

    df_scores.columns = ['Reference Spectrum ID: ' + col for col in  list(map(str,df_scores.columns.tolist()))]

    if print_id_results == True:
        print(df_top_ref_specs.to_string())

    if return_ID_output is False:
        df_top_ref_specs.to_csv(output_identification, sep='\t')
        df_scores.to_csv(output_similarity_scores, sep='\t')
    else:
        return df_top_ref_specs




def run_spec_lib_matching_on_NRMS_data_shiny(query_data=None, reference_data=None, likely_reference_ids=None, spectrum_preprocessing_order='FNLW', similarity_measure='cosine', weights={'Cosine':0.25,'Shannon':0.25,'Renyi':0.25,'Tsallis':0.25}, high_quality_reference_library=False, mz_min=0, mz_max=9999999, int_min=0, int_max=9999999, noise_threshold=0.0, wf_mz=0.0, wf_intensity=1.0, LET_threshold=0.0, entropy_dimension=1.1, n_top_matches_to_save=1, print_id_results=False, output_identification=None, output_similarity_scores=None, return_ID_output=False, verbose=True):
    if query_data is None:
        print('\nError: No argument passed to the mandatory query_data. Please pass the path to the TXT file of the query data.')
        sys.exit()
    else:
        extension = query_data.rsplit('.',1)
        extension = extension[(len(extension)-1)]
        if extension == 'mgf' or extension == 'MGF' or extension == 'mzML' or extension == 'mzml' or extension == 'MZML' or extension == 'cdf' or extension == 'CDF' or extension == 'msp' or extension == 'MSP' or extension == 'json' or extension == 'JSON':
            output_path_tmp = query_data[:-3] + 'txt'
            build_library_from_raw_data(input_path=query_data, output_path=output_path_tmp, is_reference=False)
            df_query = pd.read_csv(output_path_tmp, sep='\t')
        if extension == 'txt' or extension == 'TXT':
            df_query = pd.read_csv(query_data, sep='\t')
        unique_query_ids = df_query.iloc[:,0].unique()

    if reference_data is None:
        print('\nError: No argument passed to the mandatory reference_data. Please pass the path to the TXT file of the reference data.')
        sys.exit()
    else:
        if isinstance(reference_data,str):
            df_reference = get_reference_df(reference_data,likely_reference_ids)
            unique_reference_ids = df_reference.iloc[:,0].unique()
        else:
            dfs = []
            unique_reference_ids = []
            for f in reference_data:
                tmp = get_reference_df(f,likely_reference_ids)
                dfs.append(tmp)
                unique_reference_ids.extend(tmp.iloc[:,0].unique())
            df_reference = pd.concat(dfs, axis=0, ignore_index=True)


    if spectrum_preprocessing_order is not None:
        spectrum_preprocessing_order = list(spectrum_preprocessing_order)
    else:
        spectrum_preprocessing_order = ['F','N','W','L']
    if set(spectrum_preprocessing_order) - {'F','N','W','L'}:
        print(f'Error: spectrum_preprocessing_order must contain only \'F\', \'N\', \'W\', \'L\'.')
        sys.exit()

    if similarity_measure not in ['cosine','shannon','renyi','tsallis','mixture','jaccard','dice','3w_jaccard','sokal_sneath','binary_cosine','mountford','mcconnaughey','driver_kroeber','simpson','braun_banquet','fager_mcgowan','kuldzynski','intersection','hamming','hellinger']:
        print('\nError: similarity_measure must be either cosine, shannon, renyi, tsallis, mixture, jaccard, dice, 3w_jaccard, sokal_sneath, binary_cosine, mountford, mcconnaughey, driver_kroeber, simpson, braun_banquet, fager_mcgowan, kulczynski, intersection, hamming, or hellinger')
        sys.exit()

    if isinstance(int_min,int) is True:
        int_min = float(int_min)
    if isinstance(int_max,int) is True:
        int_max = float(int_max)
    if isinstance(mz_min,int) is False or isinstance(mz_max,int) is False or isinstance(int_min,float) is False or isinstance(int_max,float) is False:
        print('Error: mz_min must be a non-negative integer, mz_max must be a positive integer, int_min must be a non-negative float, and int_max must be a positive float')
        sys.exit()
    if mz_min < 0:
        print('\nError: mz_min should be a non-negative integer')
        sys.exit()
    if mz_max <= 0:
        print('\nError: mz_max should be a positive integer')
        sys.exit()
    if int_min < 0:
        print('\nError: int_min should be a non-negative float')
        sys.exit()
    if int_max <= 0:
        print('\nError: int_max should be a positive float')
        sys.exit()

    if isinstance(noise_threshold,int) is True:
        noise_threshold = float(noise_threshold)
    if isinstance(noise_threshold,float) is False or noise_threshold < 0:
        print('Error: noise_threshold must be a positive float.')
        sys.exit()

    if isinstance(wf_intensity,int) is True:
        wf_intensity = float(wf_intensity)
    if isinstance(wf_mz,int) is True:
        wf_mz = float(wf_mz)
    if isinstance(wf_intensity,float) is False or isinstance(wf_mz,float) is False:
        print('Error: wf_mz and wf_intensity must be integers or floats')
        sys.exit()

    if entropy_dimension <= 0:
        print('\nError: entropy_dimension should be a positive float')
        sys.exit()
    else:
        q = entropy_dimension

    normalization_method = 'standard' 

    if n_top_matches_to_save <= 0 or isinstance(n_top_matches_to_save,int)==False:
        print('\nError: n_top_matches_to_save should be a positive integer')
        sys.exit()

    if isinstance(print_id_results,bool)==False:
        print('\nError: print_id_results must be either True or False')
        sys.exit()
    
    if output_identification is None:
        output_identification = f'{Path.cwd()}/output_identification.txt'
        print(f'Warning: writing identification output to {output_identification}')

    if output_similarity_scores is None:
        output_similarity_scores = f'{Path.cwd()}/output_all_similarity_scores.txt'
        print(f'Warning: writing similarity scores to {output_similarity_scores}')



    min_mz = int(np.min([np.min(df_query.iloc[:,1]), np.min(df_reference.iloc[:,1])]))
    max_mz = int(np.max([np.max(df_query.iloc[:,1]), np.max(df_reference.iloc[:,1])]))
    mzs = np.linspace(min_mz,max_mz,(max_mz-min_mz+1))

    all_similarity_scores =  []
    for query_idx in range(0,len(unique_query_ids)):
        q_idxs_tmp = np.where(df_query.iloc[:,0] == unique_query_ids[query_idx])[0]
        q_spec_tmp = np.asarray(pd.concat([df_query.iloc[q_idxs_tmp,1], df_query.iloc[q_idxs_tmp,2]], axis=1).reset_index(drop=True))
        q_spec_tmp = convert_spec(q_spec_tmp,mzs)

        similarity_scores = []
        for ref_idx in range(0,len(unique_reference_ids)):
            if verbose is True and ref_idx % 1000 == 0:
                print(f'Query spectrum #{query_idx} has had its similarity with {ref_idx} reference library spectra computed')
            q_spec = q_spec_tmp
            r_idxs_tmp = np.where(df_reference.iloc[:,0] == unique_reference_ids[ref_idx])[0]
            r_spec_tmp = np.asarray(pd.concat([df_reference.iloc[r_idxs_tmp,1], df_reference.iloc[r_idxs_tmp,2]], axis=1).reset_index(drop=True))
            r_spec = convert_spec(r_spec_tmp,mzs)

            for transformation in spectrum_preprocessing_order:
                if np.isinf(q_spec[:,1]).sum() > 0:
                    q_spec[:,1] = np.zeros(q_spec.shape[0])
                if np.isinf(r_spec[:,1]).sum() > 0:
                    r_spec[:,1] = np.zeros(r_spec.shape[0])
                if transformation == 'W': 
                    q_spec[:,1] = wf_transform(q_spec[:,0], q_spec[:,1], wf_mz, wf_intensity)
                    r_spec[:,1] = wf_transform(r_spec[:,0], r_spec[:,1], wf_mz, wf_intensity)
                if transformation == 'L': 
                    q_spec[:,1] = LE_transform(q_spec[:,1], LET_threshold, normalization_method=normalization_method)
                    r_spec[:,1] = LE_transform(r_spec[:,1], LET_threshold, normalization_method=normalization_method)
                if transformation == 'N': 
                    q_spec = remove_noise(q_spec, nr = noise_threshold)
                    if high_quality_reference_library == False:
                        r_spec = remove_noise(r_spec, nr = noise_threshold)
                if transformation == 'F':
                    q_spec = filter_spec_gcms(q_spec, mz_min = mz_min, mz_max = mz_max, int_min = int_min, int_max = int_max)
                    if high_quality_reference_library == False:
                        r_spec = filter_spec_gcms(r_spec, mz_min = mz_min, mz_max = mz_max, int_min = int_min, int_max = int_max)

            q_ints = q_spec[:,1]
            r_ints = r_spec[:,1]

            if np.sum(q_ints) != 0 and np.sum(r_ints) != 0:
                similarity_score = get_similarity(similarity_measure, q_spec[:,1], r_spec[:,1], weights, entropy_dimension)
            else:
                similarity_score = 0

            similarity_scores.append(similarity_score)
        all_similarity_scores.append(similarity_scores)

    df_scores = pd.DataFrame(all_similarity_scores, columns = unique_reference_ids)
    df_scores.index = unique_query_ids
    df_scores.index.names = ['QUERY.SPECTRUM.ID']

    preds = []
    scores = []
    for i in range(0, df_scores.shape[0]):
        df_scores_tmp = df_scores
        preds_tmp = []
        scores_tmp = []
        for j in range(0, n_top_matches_to_save):
            top_ref_specs_tmp = df_scores_tmp.iloc[i,np.where(df_scores_tmp.iloc[i,:] == np.max(df_scores_tmp.iloc[i,:]))[0]]
            cols_to_keep = np.where(df_scores_tmp.iloc[i,:] != np.max(df_scores_tmp.iloc[i,:]))[0]
            df_scores_tmp = df_scores_tmp.iloc[:,cols_to_keep]

            preds_tmp.append(';'.join(map(str,top_ref_specs_tmp.index.to_list())))
            if len(top_ref_specs_tmp.values) == 0:
                scores_tmp.append(0)
            else:
                scores_tmp.append(top_ref_specs_tmp.values[0])
        preds.append(preds_tmp)
        scores.append(scores_tmp)

    preds = np.array(preds)
    scores = np.array(scores)
    out = np.c_[preds,scores]

    cnames_preds = []
    cnames_scores = []
    for i in range(0,n_top_matches_to_save):
        cnames_preds.append(f'RANK.{i+1}.PRED')
        cnames_scores.append(f'RANK.{i+1}.SIMILARITY.SCORE')

    df_top_ref_specs = pd.DataFrame(out, columns = [*cnames_preds, *cnames_scores])
    df_top_ref_specs.index = unique_query_ids
    df_top_ref_specs.index.names = ['QUERY.SPECTRUM.ID']

    if print_id_results == True:
        print(df_top_ref_specs.to_string())

    df_scores.columns = ['Reference Spectrum ID: ' + col for col in  list(map(str,df_scores.columns.tolist()))]

    if return_ID_output is False:
        df_top_ref_specs.to_csv(output_identification, sep='\t')
        df_scores.columns = ['Reference Spectrum ID: ' + col for col in  list(map(str,df_scores.columns.tolist()))]
        df_scores.to_csv(output_similarity_scores, sep='\t')
    else:
        return df_top_ref_specs


class _UIWriter:
    def __init__(self, loop, q: asyncio.Queue[str]):
        self._loop = loop
        self._q = q
    def write(self, s: str):
        if s:
            self._loop.call_soon_threadsafe(self._q.put_nowait, s)
        return len(s)
    def flush(self):
        pass


def attach_logging_to_writer(writer):
    handler = logging.StreamHandler(writer)
    handler.setLevel(logging.INFO)
    root = logging.getLogger()
    root.addHandler(handler)
    root.setLevel(logging.INFO)
    return handler, root



def _run_with_redirects(fn, writer, *args, **kwargs):
    with redirect_stdout(writer), redirect_stderr(writer):
        return fn(*args, **kwargs)


def strip_text(s):
    return [x.strip() for x in s.strip('[]').split(',') if x.strip()]


def strip_numeric(s):
    return [float(x.strip()) for x in s.strip('[]').split(',') if x.strip()]


def strip_weights(s):
    obj = ast.literal_eval(s) if isinstance(s, (str, bytes)) else s
    keys = ['Cosine', 'Shannon', 'Renyi', 'Tsallis']

    if isinstance(obj, (list, tuple)):
        if len(obj) == 4 and all(isinstance(x, Real) for x in obj):
            tuples = [obj]
        else:
            tuples = list(obj)
    else:
        raise ValueError(f"Expected a 4-tuple or a sequence of 4-tuples, got {type(obj).__name__}")

    out = []
    for t in tuples:
        if not (isinstance(t, (list, tuple)) and len(t) == 4):
            raise ValueError(f"Each item must be a 4-tuple, got: {t!r}")
        out.append(dict(zip(keys, t)))
    return out


def build_library(input_path=None, output_path=None):
    last_three_chars = input_path[(len(input_path)-3):len(input_path)]
    last_four_chars = input_path[(len(input_path)-4):len(input_path)]
    if last_three_chars == 'txt' or last_three_chars == 'TXT':
        return pd.read_csv(input_path, sep='\t')
    else:
        if last_three_chars == 'mgf' or last_three_chars == 'MGF':
            input_file_type = 'mgf'
        elif last_four_chars == 'mzML' or last_four_chars == 'mzml' or last_four_chars == 'MZML':
            input_file_type = 'mzML'
        elif last_four_chars == 'json' or last_four_chars == 'JSON':
            input_file_type = 'json'
        elif last_three_chars == 'cdf' or last_three_chars == 'CDF':
            input_file_type = 'cdf'
        elif last_three_chars == 'msp' or last_three_chars == 'MSP':
            input_file_type = 'msp'
        else:
            print('ERROR: either an \'mgf\', \'mzML\', \'cdf\', \'msp\', \'json\', or \'txt\' file must be passed to --input_path')
            sys.exit()

        spectra = []
        if input_file_type == 'mgf':
            with mgf.read(input_path, index_by_scans = True) as reader:
                for spec in reader:
                    spectra.append(spec)
        if input_file_type == 'mzML':
            with mzml.read(input_path) as reader:
                for spec in reader:
                    spectra.append(spec)

        if input_file_type == 'mgf' or input_file_type == 'mzML':
            ids = []
            mzs = []
            ints = []
            for i in range(0,len(spectra)):
                for j in range(0,len(spectra[i]['m/z array'])):
                    if input_file_type == 'mzML':
                        ids.append(f'ID_{i+1}')
                    else:
                        ids.append(spectra[i]['params']['name'])
                    mzs.append(spectra[i]['m/z array'][j])
                    ints.append(spectra[i]['intensity array'][j])

        if input_file_type == 'cdf':
            dataset = nc.Dataset(input_path, 'r')
            all_mzs = dataset.variables['mass_values'][:]
            all_ints = dataset.variables['intensity_values'][:]
            scan_idxs = dataset.variables['scan_index'][:]
            dataset.close()

            ids = []
            mzs = []
            ints = []
            for i in range(0,(len(scan_idxs)-1)):
                if i % 1000 == 0:
                    print(f'analyzed {i} out of {len(scan_idxs)} scans')
                s_idx = scan_idxs[i]
                e_idx = scan_idxs[i+1]

                mzs_tmp = all_mzs[s_idx:e_idx]
                ints_tmp = all_ints[s_idx:e_idx]

                for j in range(0,len(mzs_tmp)):
                    ids.append(f'ID_{i+1}')
                    mzs.append(mzs_tmp[j])
                    ints.append(ints_tmp[j])

        if input_file_type == 'msp':
            ids = []
            mzs = []
            ints = []
            with open(input_path, 'r') as f:
                i = 0
                for line in f:
                    line = line.strip()
                    if line.startswith('Name:'):
                        i += 1
                        spectrum_id = line.replace('Name: ','')
                    elif line and line[0].isdigit():
                        try:
                            mz, intensity = map(float, line.split()[:2])
                            ids.append(spectrum_id)
                            mzs.append(mz)
                            ints.append(intensity)
                        except ValueError:
                            continue

        if input_file_type == 'json':
            data = json.load(open(input_path))
            ids = []
            mzs = []
            ints = []
            for i in range(0,len(data)):
                spec_ID_tmp = data[i]['spectrum_id']
                tmp = data[i]['peaks_json']
                tmp = tmp[1:-1].split(",")
                tmp = [a.replace("[","") for a in tmp]
                tmp = [a.replace("]","") for a in tmp]
                mzs_tmp = tmp[0::2]
                ints_tmp = tmp[1::2]
                ids.extend([spec_ID_tmp] * len(mzs_tmp))
                mzs.extend(mzs_tmp)
                ints.extend(ints_tmp)

        df = pd.DataFrame({'id':ids, 'mz_ratio':mzs, 'intensity':ints})
        return df



def extract_first_column_ids(file_path: str, max_ids: int = 20000):
    suffix = Path(file_path).suffix.lower()

    if suffix == ".txt":
        df = pd.read_csv(file_path, sep='\t')
        if 'id' in df.columns.tolist():
            ids = df['id'].astype(str).dropna()
        else:
            ids = df.iloc[:, 0].astype(str).dropna()
        ids = [x for x in ids if x.strip() != ""]
        seen = set()
        uniq = []
        for x in ids:
            if x not in seen:
                uniq.append(x)
                seen.add(x)
        return uniq[:max_ids]

    ids = []
    try:
        with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                ls = line.strip()
                if ls.startswith("TITLE="):
                    ids.append(ls.split("=", 1)[1].strip())
                elif ls.lower().startswith("name:"):
                    ids.append(ls.split(":", 1)[1].strip())
                if len(ids) >= max_ids:
                    break
    except Exception:
        pass

    if ids:
        seen = set()
        uniq = []
        for x in ids:
            if x not in seen:
                uniq.append(x)
                seen.add(x)
        return uniq
    return []


def _open_plot_window(session, svg_bytes: bytes, title: str = "plot.svg"):
    """Send SVG bytes to browser and open in a new window as a data URL."""
    b64 = base64.b64encode(svg_bytes).decode("ascii")
    data_url = f"data:image/svg;base64,{b64}"
    session.send_custom_message("open-plot-window", {"svg": data_url, "title": title})


def plot_spectra_ui(platform: str):
    base_inputs = [
        ui.input_file("query_data", "Upload query dataset (mgf, mzML, cdf, msp, or txt):"),
        ui.input_file("reference_data", "Upload reference dataset (mgf, mzML, cdf, msp, or txt):"),
        ui.input_selectize(
            "spectrum_ID1",
            "Select spectrum ID 1 (default is the first spectrum in the library):",
            choices=[],
            multiple=False,
            options={"placeholder": "Upload a library..."},
        ),
        ui.input_selectize(
            "spectrum_ID2",
            "Select spectrum ID 2 (default is the first spectrum in the library):",
            choices=[],
            multiple=False,
            options={"placeholder": "Upload a library..."},
        ),
        ui.input_select('print_url_spectrum1', 'Print PubChem URL for spectrum 1:', ['No', 'Yes']),
        ui.input_select('print_url_spectrum2', 'Print PubChem URL for spectrum 2:', ['No', 'Yes']),
        ui.input_select("similarity_measure", "Select similarity measure:", ["cosine","shannon","renyi","tsallis","mixture","jaccard","dice","3w_jaccard","sokal_sneath","binary_cosine","mountford","mcconnaughey","driver_kroeber","simpson","braun_banquet","fager_mcgowan","kulczynski","intersection","hamming","hellinger"]),
        ui.input_text('weights', 'Weights for mixture similarity measure (cosine, shannon, renyi, tsallis):', '0.25, 0.25, 0.25, 0.25'),
        ui.input_select(
            "high_quality_reference_library",
            "Indicate whether the reference library is considered high quality. If True, filtering and noise removal are only applied to the query spectra.",
            [False, True],
        ),
    ]

    if platform == "HRMS":
        extra_inputs = [
            ui.input_text("spectrum_preprocessing_order", "Sequence of characters for preprocessing order (C (centroiding), F (filtering), M (matching), N (noise removal), L (low-entropy transformation), W (weight factor transformation)). M must be included, C before M if used.", "FCNMWL",),
            ui.input_numeric("window_size_centroiding", "Centroiding window-size:", 0.5),
            ui.input_numeric("window_size_matching", "Matching window-size:", 0.5),
        ]
    else:
        extra_inputs = [
            ui.input_text("spectrum_preprocessing_order", "Sequence of characters for preprocessing order (F (filtering), N (noise removal), L (low-entropy transformation), W (weight factor transformation)).", "FNLW",)
        ]

    numeric_inputs = [
        ui.input_numeric("mz_min", "Minimum m/z for filtering:", 0),
        ui.input_numeric("mz_max", "Maximum m/z for filtering:", 99999999),
        ui.input_numeric("int_min", "Minimum intensity for filtering:", 0),
        ui.input_numeric("int_max", "Maximum intensity for filtering:", 999999999),
        ui.input_numeric("noise_threshold", "Noise removal threshold:", 0.0),
        ui.input_numeric("wf_mz", "Mass/charge weight factor:", 0.0),
        ui.input_numeric("wf_int", "Intensity weight factor:", 1.0),
        ui.input_numeric("LET_threshold", "Low-entropy threshold:", 0.0),
        ui.input_numeric("entropy_dimension", "Entropy dimension (Renyi/Tsallis only):", 1.1),
    ]

    select_input = ui.input_select("y_axis_transformation", "Transformation to apply to intensity axis:", ["normalized", "none", "log10", "sqrt"])

    run_button_plot_spectra = ui.download_button("run_btn_plot_spectra", "Run", style="font-size:16px; padding:15px 30px; width:200px; height:80px")
    back_button = ui.input_action_button("back", "Back to main menu", style="font-size:16px; padding:15px 30px; width:200px; height:80px")

    if platform == "HRMS":
        inputs_columns = ui.layout_columns(
            ui.div(base_inputs[0:6], style="display:flex; flex-direction:column; gap:10px;"),
            ui.div([base_inputs[6:9], extra_inputs[0]], style="display:flex; flex-direction:column; gap:10px;"),
            ui.div(extra_inputs[1:3], numeric_inputs[0:3], style="display:flex; flex-direction:column; gap:10px;"),
            ui.div([numeric_inputs[3:10], select_input], style="display:flex; flex-direction:column; gap:10px;"),
            col_widths=(3,3,3,3),
        )
    elif platform == "NRMS":
        inputs_columns = ui.layout_columns(
            ui.div(base_inputs[0:6], style="display:flex; flex-direction:column; gap:10px;"),
            ui.div([base_inputs[6:9], *extra_inputs], style="display:flex; flex-direction:column; gap:10px;"),
            ui.div(numeric_inputs[0:5], style="display:flex; flex-direction:column; gap:10px;"),
            ui.div([numeric_inputs[5:10], select_input], style="display:flex; flex-direction:column; gap:10px;"),
            col_widths=(3,3,3,3),
        )

    return ui.div(
        ui.TagList(
            ui.h2("Plot Spectra"),
            inputs_columns,
            run_button_plot_spectra,
            back_button,
            ui.div(ui.output_text("plot_query_status"), style="margin-top:8px; font-size:14px"),
            ui.div(ui.output_text("plot_reference_status"), style="margin-top:8px; font-size:14px")
        ),
    )



def run_spec_lib_matching_ui(platform: str):
    base_inputs = [
        ui.input_file("query_data", "Upload query dataset (mgf, mzML, cdf, msp, or txt):"),
        ui.input_file("reference_data", "Upload reference dataset (mgf, mzML, cdf, msp, or txt):"),
        ui.input_select("similarity_measure", "Select similarity measure:", ["cosine","shannon","renyi","tsallis","mixture","jaccard","dice","3w_jaccard","sokal_sneath","binary_cosine","mountford","mcconnaughey","driver_kroeber","simpson","braun_banquet","fager_mcgowan","kulczynski","intersection","hamming","hellinger"]),
        ui.input_text('weights', 'Weights for mixture similarity measure (cosine, shannon, renyi, tsallis):', '0.25, 0.25, 0.25, 0.25'),
        ui.input_file('compound_ID_output_file', 'Upload output from spectral library matching to plot top matches (optional)'),
        ui.input_selectize("q_spec", "Select query spectrum (only applicable for plotting; default is the first spectrum in the compound ID output):", choices=[], multiple=False, options={"placeholder": "Upload compound ID output..."}),
        ui.input_selectize("r_spec", "Select reference spectrum (only applicable for plotting; default is the rank 1 reference spectrum):", choices=[], multiple=False, options={"placeholder": "Upload compound ID output..."}),
        ui.input_select('print_url_spectrum1', 'Print PubChem URL for query spectrum (only applicable for plotting):', ['No', 'Yes']),
        ui.input_select('print_url_spectrum2', 'Print PubChem URL for reference spectrum (only applicable for plotting):', ['No', 'Yes']),
        ui.input_select("high_quality_reference_library", "Indicate whether the reference library is considered high quality. If True, filtering and noise removal are only applied to the query spectra.", [False, True])
    ]

    if platform == "HRMS":
        extra_inputs = [
            ui.input_numeric("precursor_ion_mz_tolerance", "Precursor ion mass tolerance (leave blank if not applicable):", None),
            ui.input_select("ionization_mode", "Ionization mode:", ['Positive','Negative','N/A'], selected='N/A'),
            ui.input_select("adduct", "Adduct:", ['H','NH3','NH4','Na','K','N/A'], selected='N/A'),
            ui.input_text("spectrum_preprocessing_order","Sequence of characters for preprocessing order (C (centroiding), F (filtering), M (matching), N (noise removal), L (low-entropy transformation), W (weight factor transformation)). M must be included, C before M if used.","FCNMWL"),
            ui.input_numeric("window_size_centroiding", "Centroiding window-size:", 0.5),
            ui.input_numeric("window_size_matching", "Matching window-size:", 0.5),
        ]
    else:
        extra_inputs = [ui.input_text("spectrum_preprocessing_order","Sequence of characters for preprocessing order (F (filtering), N (noise removal), L (low-entropy transformation), W (weight factor transformation)).","FNLW")]

    numeric_inputs = [
        ui.input_numeric("mz_min", "Minimum m/z for filtering:", 0),
        ui.input_numeric("mz_max", "Maximum m/z for filtering:", 99999999),
        ui.input_numeric("int_min", "Minimum intensity for filtering:", 0),
        ui.input_numeric("int_max", "Maximum intensity for filtering:", 999999999),
        ui.input_numeric("noise_threshold", "Noise removal threshold:", 0.0),
        ui.input_numeric("wf_mz", "Mass/charge weight factor:", 0.0),
        ui.input_numeric("wf_int", "Intensity weight factor:", 1.0),
        ui.input_numeric("LET_threshold", "Low-entropy threshold:", 0.0),
        ui.input_numeric("entropy_dimension", "Entropy dimension (Renyi/Tsallis only):", 1.1),
        ui.input_numeric("n_top_matches_to_save", "Number of top matches to save:", 3),
    ]


    run_button_spec_lib_matching = ui.download_button("run_btn_spec_lib_matching", "Run Spectral Library Matching", style="font-size:16px; padding:15px 30px; width:200px; height:80px")
    run_button_plot_spectra_within_spec_lib_matching = ui.download_button("run_btn_plot_spectra_within_spec_lib_matching", "Plot Spectra", style="font-size:16px; padding:15px 30px; width:200px; height:80px")
    back_button = ui.input_action_button("back", "Back to main menu", style="font-size:16px; padding:15px 30px; width:200px; height:80px")

    if platform == "HRMS":
        inputs_columns = ui.layout_columns(
            ui.div([base_inputs[0:2], extra_inputs[0:3], base_inputs[2:4]], style="display:flex; flex-direction:column; gap:10px;"),
            ui.div([base_inputs[4:10]], style="display:flex; flex-direction:column; gap:10px;"),
            ui.div([extra_inputs[3:6], numeric_inputs[0:3]], style="display:flex; flex-direction:column; gap:10px;"),
            ui.div(numeric_inputs[3:10], style="display:flex; flex-direction:column; gap:10px;"),
            col_widths=(3,3,3,3)
        )
    elif platform == "NRMS":
        inputs_columns = ui.layout_columns(
            ui.div(base_inputs[0:6], style="display:flex; flex-direction:column; gap:10px;"),
            ui.div([base_inputs[6:10], *extra_inputs], style="display:flex; flex-direction:column; gap:10px;"),
            ui.div(numeric_inputs[0:5], style="display:flex; flex-direction:column; gap:10px;"),
            ui.div(numeric_inputs[5:10], style="display:flex; flex-direction:column; gap:10px;"),
            col_widths=(3,3,3,3)
        )

    log_panel = ui.card(
        ui.card_header("Identification log"),
        ui.output_text_verbatim("match_log"),
        style="max-height:300px; overflow:auto"
    )

    return ui.div(
        ui.TagList(
            ui.h2("Run Spectral Library Matching"),
            inputs_columns,
            run_button_spec_lib_matching,
            run_button_plot_spectra_within_spec_lib_matching,
            back_button,
            log_panel
        ),
    )



def run_parameter_tuning_grid_ui(platform: str):
    base_inputs = [
        ui.input_file("query_data", "Upload query dataset (mgf, mzML, cdf, msp, or txt):"),
        ui.input_file("reference_data", "Upload reference dataset (mgf, mzML, cdf, msp, or txt):"),
        ui.input_selectize("similarity_measure", "Select similarity measure(s):", ["cosine","shannon","renyi","tsallis","mixture","jaccard","dice","3w_jaccard","sokal_sneath","binary_cosine","mountford","mcconnaughey","driver_kroeber","simpson","braun_banquet","fager_mcgowan","kulczynski","intersection","hamming","hellinger"], multiple=True, selected='cosine'),
        ui.input_text('weights', 'Weights for mixture similarity measure (cosine, shannon, renyi, tsallis):', '((0.25, 0.25, 0.25, 0.25))'),
        ui.input_text("high_quality_reference_library", "Indicate whether the reference library is considered high quality. If True, filtering and noise removal are only applied to the query spectra.", '[True]')
    ]

    if platform == "HRMS":
        extra_inputs = [
            ui.input_numeric("precursor_ion_mz_tolerance", "Precursor ion mass tolerance (leave blank if not applicable):", None),
            ui.input_select("ionization_mode", "Ionization mode:", ['Positive','Negative','N/A'], selected='N/A'),
            ui.input_select("adduct", "Adduct:", ['H','NH3','NH4','Na','K','N/A'], selected='N/A'),
            ui.input_text("spectrum_preprocessing_order", "Sequence of characters for preprocessing order (C (centroiding), F (filtering), M (matching), N (noise removal), L (low-entropy transformation), W (weight factor transformation)). M must be included, C before M if used.", "[FCNMWL,CWM]"),
            ui.input_text("window_size_centroiding", "Centroiding window-size:", "[0.5]"),
            ui.input_text("window_size_matching", "Matching window-size:", "[0.1,0.5]"),
        ]
    else:
        extra_inputs = [
            ui.input_text(
                "spectrum_preprocessing_order",
                "Sequence of characters for preprocessing order (F (filtering), N (noise removal), L (low-entropy transformation), W (weight factor transformation)).",
                "[FNLW,WNL]",
            )
        ]

    numeric_inputs = [
        ui.input_text("mz_min", "Minimum m/z for filtering:", '[0]'),
        ui.input_text("mz_max", "Maximum m/z for filtering:", '[99999999]'),
        ui.input_text("int_min", "Minimum intensity for filtering:", '[0]'),
        ui.input_text("int_max", "Maximum intensity for filtering:", '[999999999]'),
        ui.input_text("noise_threshold", "Noise removal threshold:", '[0.0]'),
        ui.input_text("wf_mz", "Mass/charge weight factor:", '[0.0]'),
        ui.input_text("wf_int", "Intensity weight factor:", '[1.0]'),
        ui.input_text("LET_threshold", "Low-entropy threshold:", '[0.0]'),
        ui.input_text("entropy_dimension", "Entropy dimension (Renyi/Tsallis only):", '[1.1]')
    ]


    run_button_parameter_tuning_grid = ui.download_button("run_btn_parameter_tuning_grid", "Tune parameters (grid search)", style="font-size:16px; padding:15px 30px; width:200px; height:80px")
    back_button = ui.input_action_button("back", "Back to main menu", style="font-size:16px; padding:15px 30px; width:200px; height:80px")

    if platform == "HRMS":
        inputs_columns = ui.layout_columns(
            ui.div(base_inputs[0:6], style="display:flex; flex-direction:column; gap:10px;"),
            ui.div([base_inputs[6:7], *extra_inputs], style="display:flex; flex-direction:column; gap:10px;"),
            ui.div(numeric_inputs[0:5], style="display:flex; flex-direction:column; gap:10px;"),
            ui.div(numeric_inputs[5:9], style="display:flex; flex-direction:column; gap:10px;"),
            col_widths=(3, 3, 3, 3),
        )
    elif platform == "NRMS":
        inputs_columns = ui.layout_columns(
            ui.div(base_inputs[0:6], style="display:flex; flex-direction:column; gap:10px;"),
            ui.div([base_inputs[6:7], *extra_inputs], style="display:flex; flex-direction:column; gap:10px;"),
            ui.div(numeric_inputs[0:5], style="display:flex; flex-direction:column; gap:10px;"),
            ui.div(numeric_inputs[5:9], style="display:flex; flex-direction:column; gap:10px;"),
            col_widths=(3, 3, 3, 3),
        )

    log_panel = ui.card(
        ui.card_header("Identification log"),
        ui.output_text_verbatim("match_log"),
        style="max-height:300px; overflow:auto"
    )

    return ui.div(
        ui.TagList(
            ui.h2("Tune parameters (grid search)"),
            inputs_columns,
            run_button_parameter_tuning_grid,
            back_button,
            log_panel
        ),
    )



PARAMS_HRMS = {
    "window_size_centroiding": (0.0, 0.5),
    "window_size_matching":    (0.0, 0.5),
    "noise_threshold":         (0.0, 0.25),
    "wf_mz":                   (0.0, 5.0),
    "wf_int":                  (0.0, 5.0),
    "LET_threshold":           (0.0, 5.0),
    "entropy_dimension":       (1.0, 3.0)
}

PARAMS_NRMS = {
    "noise_threshold":         (0.0, 0.25),
    "wf_mz":                   (0.0, 5.0),
    "wf_int":                  (0.0, 5.0),
    "LET_threshold":           (0.0, 5.0),
    "entropy_dimension":       (1.0, 3.0)
}


def run_parameter_tuning_DE_ui(platform: str):
    # Pick param set per platform
    if platform == "HRMS":
        PARAMS = PARAMS_HRMS
    else:
        PARAMS = PARAMS_NRMS

    base_inputs = [
        ui.input_file("query_data", "Upload query dataset (mgf, mzML, cdf, msp, or txt):"),
        ui.input_file("reference_data", "Upload reference dataset (mgf, mzML, cdf, msp, or txt):"),
        ui.input_select("similarity_measure", "Select similarity measure:", ["cosine","shannon","renyi","tsallis","mixture","jaccard","dice","3w_jaccard","sokal_sneath","binary_cosine","mountford","mcconnaughey","driver_kroeber","simpson","braun_banquet","fager_mcgowan","kulczynski","intersection","hamming","hellinger"]),
        ui.input_text("weights", "Weights for mixture similarity measure (cosine, shannon, renyi, tsallis):", "0.25, 0.25, 0.25, 0.25"),
        ui.input_select("high_quality_reference_library", "Indicate whether the reference library is considered high quality. If True, filtering and noise removal are only applied to the query spectra.", [False, True])]

    if platform == "HRMS":
        extra_inputs = [
            ui.input_numeric("precursor_ion_mz_tolerance", "Precursor ion mass tolerance (leave blank if not applicable):", None),
            ui.input_select("ionization_mode", "Ionization mode:", ['Positive','Negative','N/A'], selected='N/A'),
            ui.input_select("adduct", "Adduct:", ['H','NH3','NH4','Na','K','N/A'], selected='N/A'),
            ui.input_text("spectrum_preprocessing_order", "Sequence of characters for preprocessing order (C (centroiding), F (filtering), M (matching), N (noise removal), L (low-entropy transformation), W (weight factor transformation)). M must be included, C before M if used.", "FCNMWL"),
            ui.input_numeric("window_size_centroiding", "Centroiding window-size:", 0.5),
            ui.input_numeric("window_size_matching", "Matching window-size:", 0.5),
        ]
    else:
        extra_inputs = [ui.input_text("spectrum_preprocessing_order", "Sequence of characters for preprocessing order (F (filtering), N (noise removal), L (low-entropy transformation), W (weight factor transformation)).", "FNLW")]

    numeric_inputs = [
        ui.input_numeric("mz_min", "Minimum m/z for filtering:", 0),
        ui.input_numeric("mz_max", "Maximum m/z for filtering:", 99_999_999),
        ui.input_numeric("int_min", "Minimum intensity for filtering:", 0),
        ui.input_numeric("int_max", "Maximum intensity for filtering:", 999_999_999),
        ui.input_numeric("noise_threshold", "Noise removal threshold:", 0.0),
        ui.input_numeric("wf_mz", "Mass/charge weight factor:", 0.0),
        ui.input_numeric("wf_int", "Intensity weight factor:", 1.0),
        ui.input_numeric("LET_threshold", "Low-entropy threshold:", 0.0),
        ui.input_numeric("entropy_dimension", "Entropy dimension (Renyi/Tsallis only):", 1.1),
        ui.input_numeric("max_iterations", "Maximum number of iterations:", 5),
    ]

    run_button_parameter_tuning_DE = ui.input_action_button("run_btn_parameter_tuning_DE", "Tune parameters (differential evolution optimization)", style="font-size:16px; padding:15px 30px; width:300px; height:100px")
    back_button = ui.input_action_button("back", "Back to main menu", style="font-size:16px; padding:15px 30px; width:300px; height:100px")

    if platform == "HRMS":
        inputs_columns = ui.layout_columns(
            ui.div(*base_inputs, style="display:flex; flex-direction:column; gap:10px;"),
            ui.div(*extra_inputs, style="display:flex; flex-direction:column; gap:10px;"),
            ui.div(*numeric_inputs[0:5], style="display:flex; flex-direction:column; gap:10px;"),
            ui.div(*numeric_inputs[5:11], style="display:flex; flex-direction:column; gap:10px;"),
            col_widths=(3, 3, 3, 3),
        )
    else:
        inputs_columns = ui.layout_columns(
            ui.div(*base_inputs, style="display:flex; flex-direction:column; gap:10px;"),
            ui.div(*extra_inputs, style="display:flex; flex-direction:column; gap:10px;"),
            ui.div(*numeric_inputs[0:5], style="display:flex; flex-direction:column; gap:10px;"),
            ui.div(*numeric_inputs[5:11], style="display:flex; flex-direction:column; gap:10px;"),
            col_widths=(3, 3, 3, 3),
        )

    return ui.page_fillable(
        ui.layout_sidebar(
            ui.sidebar(
                ui.h3("Select continuous parameters to optimize"),
                ui.input_checkbox_group("params", None, choices=list(PARAMS.keys()), selected=["noise_threshold", "LET_threshold"]),
                ui.hr(),
                ui.h4("Bounds for selected parameters"),
                ui.output_ui("bounds_inputs"),
                width=360,
            ),
            ui.div(
                ui.h2("Tune parameters (differential evolution optimization)"),
                inputs_columns,
                ui.div(run_button_parameter_tuning_DE, back_button, style=("display:flex; flex-direction:row; gap:12px; align-items:center; flex-wrap:wrap;")),
                ui.br(),
                ui.card(
                    ui.card_header("Live log"),
                    ui.output_text_verbatim("run_log"),
                ),
                style="display:flex; flex-direction:column; gap:16px;",
            ),
        )
    )



app_ui = ui.page_fluid(
    ui.head_content(ui.tags.link(rel="icon", href="emblem.png")),
    ui.div(ui.output_image("image"), style=("display:block; margin:20px auto; max-width:320px; height:auto; text-align:center")),
    ui.output_ui("main_ui"),
    ui.output_text("status_output"),
)




def server(input, output, session):

    current_page = reactive.Value("main_menu")
    
    plot_clicks = reactive.Value(0)
    match_clicks = reactive.Value(0)
    back_clicks = reactive.Value(0)

    run_status_plot_spectra = reactive.Value("")
    run_status_spec_lib_matching = reactive.Value("")
    run_status_plot_spectra_within_spec_lib_matching = reactive.Value("")
    run_status_parameter_tuning_grid = reactive.Value("")
    run_status_parameter_tuning_DE = reactive.Value("")
    is_tuning_grid_running = reactive.Value(False)
    is_tuning_DE_running = reactive.Value(False)
    match_log_rv = reactive.Value("")
    is_matching_rv = reactive.Value(False)
    is_any_job_running = reactive.Value(False)
    latest_txt_path_rv = reactive.Value("")
    latest_df_rv = reactive.Value(None)
    is_running_rv = reactive.Value(False)

    query_ids_rv = reactive.Value([])
    query_file_path_rv = reactive.Value(None)
    query_result_rv = reactive.Value(None)
    query_status_rv = reactive.Value("")
    reference_ids_rv = reactive.Value([])
    reference_file_path_rv = reactive.Value(None)
    reference_result_rv = reactive.Value(None)
    reference_status_rv = reactive.Value("")

    converted_query_path_rv = reactive.Value(None)
    converted_reference_path_rv = reactive.Value(None)

    df_rv = reactive.Value(None)


    def _discover_rank_cols(df: pd.DataFrame):
        pred_pat = re.compile(r"^RANK\.(\d+)\.PRED$")
        score_pat = re.compile(r"^RANK\.(\d+)\.SIMILARITY\.SCORE$")
        pred_map, score_map = {}, {}
        for c in df.columns:
            m = pred_pat.match(c)
            if m: pred_map[int(m.group(1))] = c
            m = score_pat.match(c)
            if m: score_map[int(m.group(1))] = c
        return [(k, pred_map[k], score_map.get(k)) for k in sorted(pred_map)]


    def _rank_choices_for_query(df: pd.DataFrame, qid: str):
        sub = df.loc[df["QUERY.SPECTRUM.ID"].astype(str) == str(qid)]
        if sub.empty:
            return {}, None
        row = sub.iloc[0]
        rank_cols = _discover_rank_cols(df)
        if not rank_cols:
            return {}, None

        choices = {}
        default_value = None
        for (k, pred_col, score_col) in rank_cols:
            pred = row.get(pred_col, None)
            if pd.isna(pred):
                continue
            pred = str(pred)
            score = row.get(score_col, None) if score_col else None
            score_str = f"{float(score):.6f}" if (score is not None and pd.notna(score)) else "NA"
            label = f"Rank {k} — {score_str} — {pred}"
            choices[label] = pred                 # values are plain names
            if k == 1:
                default_value = pred              # default = Rank 1 name

        if default_value is None and choices:
            default_value = next(iter(choices.values()))
        return choices, default_value


    @reactive.effect
    @reactive.event(input.compound_ID_output_file)
    async def _populate_ids_from_compound_ID_output_upload():
        files = input.compound_ID_output_file()
        if not files:
            return

        in_path = Path(files[0]["datapath"])
        try:
            query_status_rv.set(f"Reading table from: {in_path.name} …")
            await reactive.flush()

            df = await asyncio.to_thread(pd.read_csv, in_path, sep="\t", header=0)

            if "QUERY.SPECTRUM.ID" not in df.columns:
                raise ValueError("Missing required column: QUERY.SPECTRUM.ID")
            if not _discover_rank_cols(df):
                raise ValueError("No columns matching RANK.<k>.PRED found.")

            df_rv.set(df)

            ids = df["QUERY.SPECTRUM.ID"].astype(str).tolist()
            unique_ids_in_order = list(dict.fromkeys(ids))

            choices_dict, default_rank_value = _rank_choices_for_query(df, ids[0])
            choices_values = [str(v).strip() for v in choices_dict.values()]
            default_rank_value = str(default_rank_value).strip() if default_rank_value is not None else None

            ui.update_selectize("q_spec", choices=unique_ids_in_order, selected=ids[0])
            await reactive.flush()

            ui.update_selectize("r_spec", choices=choices_values, selected=choices_values[0])
            await reactive.flush()

        except Exception as e:
            query_status_rv.set(f"❌ Failed: {e}")
            await reactive.flush()
            raise


    @reactive.effect
    @reactive.event(input.q_spec)
    async def _update_rank_choices_on_compound_ID_change():
        df = df_rv.get()
        if df is None:
            return
        qid = input.q_spec()
        if not qid:
            return

        choices, default_rank_value = _rank_choices_for_query(df, qid)
        choices = list(choices.values())
        ui.update_selectize('r_spec', choices=choices, selected=default_rank_value)
        await reactive.flush()



    @output
    @render.ui
    def bounds_inputs():
        selected = input.params()
        if not selected:
            return ui.div(ui.em("Select one or more parameters above."))

        if input.chromatography_platform() == 'HRMS':
            PARAMS = PARAMS_HRMS
        else:
            PARAMS = PARAMS_NRMS
        blocks = []
        for name in selected:
            lo, hi = PARAMS.get(name, (0.0, 1.0))
            blocks.append(
                ui.card(
                    ui.card_header(name),
                    ui.layout_columns(
                        ui.input_numeric(f"min_{name}", "Lower", lo, step=0.001),
                        ui.input_numeric(f"max_{name}", "Upper", hi, step=0.001),
                    )
                )
            )
        return ui.div(*blocks)

    def _read_bounds_dict():
        selected = input.params()
        out = {}
        for name in selected:
            lo_default, hi_default = PARAMS.get(name, (0.0, 1.0))
            lo_id = f"min_{name}"
            hi_id = f"max_{name}"

            lo_val = input[lo_id]() if lo_id in input else lo_default
            hi_val = input[hi_id]() if hi_id in input else hi_default

            out[name] = (float(lo_val), float(hi_val))
        return out

    def _read_bounds():
        opt_params = input.params()
        bounds_dict = {}
        if input.chromatography_platform() == 'HRMS':
            PARAMS = PARAMS_HRMS
        else:
            PARAMS = PARAMS_NRMS

        for p in opt_params:
            lo_id, hi_id = f"min_{p}", f"max_{p}"
            lo_default, hi_default = PARAMS.get(p, (0.0, 1.0))
            lo = input[lo_id]() if lo_id in input else lo_default
            hi = input[hi_id]() if hi_id in input else hi_default
            if lo > hi:
                lo, hi = hi, lo
            bounds_dict[p] = (float(lo), float(hi))

        bounds_list = [bounds_dict[p] for p in opt_params]
        return opt_params, bounds_dict, bounds_list

    def _reset_plot_spectra_state():
        query_status_rv.set("")
        reference_status_rv.set("")
        query_ids_rv.set([])
        reference_ids_rv.set([])
        query_file_path_rv.set(None)
        reference_file_path_rv.set(None)
        query_result_rv.set(None)
        reference_result_rv.set(None)
        converted_query_path_rv.set(None)
        converted_reference_path_rv.set(None)
        try:
            ui.update_selectize("spectrum_ID1", choices=[], selected=None)
            ui.update_selectize("spectrum_ID2", choices=[], selected=None)
        except Exception:
            pass


    def _reset_spec_lib_matching_state():
        match_log_rv.set("")
        is_matching_rv.set(False)
        is_any_job_running.set(False)
        try:
            ui.update_selectize("spectrum_ID1", choices=[], selected=None)
            ui.update_selectize("spectrum_ID2", choices=[], selected=None)
        except Exception:
            pass


    def _reset_parameter_tuning_state():
        match_log_rv.set("")
        is_tuning_grid_running.set(False)
        is_tuning_DE_running.set(False)
        is_any_job_running.set(False)


    @reactive.effect
    @reactive.event(input.back)
    def _clear_on_back_from_pages():
        page = current_page()
        if page == "plot_spectra":
            _reset_plot_spectra_state()
        elif page == "run_spec_lib_matching":
            _reset_spec_lib_matching_state()
        elif page == "run_parameter_tuning_grid":
            _reset_parameter_tuning_state()
        elif page == "run_parameter_tuning_DE":
            _reset_parameter_tuning_state()

    @reactive.effect
    def _clear_on_enter_pages():
        page = current_page()
        if page == "plot_spectra":
            _reset_plot_spectra_state()
        elif page == "run_spec_lib_matching":
            _reset_spec_lib_matching_state()
        elif page == "run_parameter_tuning_grid":
            _reset_parameter_tuning_state()
        elif page == "run_parameter_tuning_DE":
            _reset_parameter_tuning_state()


    def _drain_queue_nowait(q: asyncio.Queue) -> list[str]:
        out = []
        try:
            while True:
                out.append(q.get_nowait())
        except asyncio.QueueEmpty:
            pass
        return out


    class ReactiveWriter(io.TextIOBase):
        def __init__(self, loop: asyncio.AbstractEventLoop):
            self._loop = loop
        def write(self, s: str):
            if not s:
                return 0
            self._loop.call_soon_threadsafe(_LOG_QUEUE.put_nowait, s)
            return len(s)
        def flush(self):
            pass

    def _run_with_redirects(func, writer: ReactiveWriter, **kwargs):
        with contextlib.redirect_stdout(writer), contextlib.redirect_stderr(writer):
            return func(**kwargs)



    @reactive.effect
    async def _pump_logs():
        if not (is_any_job_running.get() or is_tuning_grid_running.get() or is_tuning_DE_running.get() or is_matching_rv.get()):
            return
        reactive.invalidate_later(0.05)
        msgs = _drain_queue_nowait(_LOG_QUEUE)
        if msgs:
            match_log_rv.set(match_log_rv.get() + "".join(msgs))
            await reactive.flush()


    def process_database(file_path: str):
        suffix = Path(file_path).suffix.lower()
        return {"path": file_path, "suffix": suffix}

    @render.text
    def plot_query_status():
        return query_status_rv.get() or ""

    @render.text
    def plot_reference_status():
        return reference_status_rv.get() or ""


    @reactive.effect
    @reactive.event(input.query_data)
    async def _on_query_upload():
        files = input.query_data()
        req(files and len(files) > 0)

        file_path = files[0]["datapath"]
        query_file_path_rv.set(file_path)

        query_status_rv.set(f"Processing query database: {Path(file_path).name} …")
        await reactive.flush()

        try:
            result = await asyncio.to_thread(process_database, file_path)
            query_result_rv.set(result)
            query_status_rv.set("✅ Query database processed.")
            await reactive.flush()
        except Exception as e:
            query_status_rv.set(f"❌ Failed to process query database: {e}")
            await reactive.flush()


    @reactive.effect
    @reactive.event(input.reference_data)
    async def _on_reference_upload():
        files = input.reference_data()
        req(files and len(files) > 0)

        file_path = files[0]["datapath"]
        reference_file_path_rv.set(file_path)

        reference_status_rv.set(f"Processing reference database: {Path(file_path).name} …")
        await reactive.flush()

        try:
            result = await asyncio.to_thread(process_database, file_path)
            reference_result_rv.set(result)
            reference_status_rv.set("✅ Reference database processed.")
            await reactive.flush()
        except Exception as e:
            reference_status_rv.set(f"❌ Failed to process reference database: {e}")
            await reactive.flush()


    @render.text
    def match_log():
        return match_log_rv.get()


    @reactive.Effect
    def _():
        if input.plot_spectra() > plot_clicks.get():
            current_page.set("plot_spectra")
            plot_clicks.set(input.plot_spectra())
        elif input.run_spec_lib_matching() > match_clicks.get():
            current_page.set("run_spec_lib_matching")
            match_clicks.set(input.run_spec_lib_matching())
        elif input.run_parameter_tuning_grid() > match_clicks.get():
            current_page.set("run_parameter_tuning_grid")
            match_clicks.set(input.run_parameter_tuning_grid())
        elif input.run_parameter_tuning_DE() > match_clicks.get():
            current_page.set("run_parameter_tuning_DE")
            match_clicks.set(input.run_parameter_tuning_DE())
        elif hasattr(input, "back") and input.back() > back_clicks.get():
            current_page.set("main_menu")
            back_clicks.set(input.back())


    @render.image
    def image():
        dir = Path(__file__).resolve().parent
        img: ImgData = {"src": str(dir / "www/emblem.png"), "width": "250px", "height": "250px"}
        return img

    @output
    @render.ui
    def main_ui():
        if current_page() == "main_menu":
            return ui.page_fluid(
                ui.h2("Main Menu"),
                ui.div("Overview:", style="text-align:left; font-size:24px; font-weight:bold"),
                ui.div("PyCompound is a Python-based tool designed for performing spectral library matching on either high-resolution mass spectrometry data (HRMS) or low-resolution mass spectrometry data (NRMS). PyCompound offers a range of spectrum preprocessing transformations and similarity measures. These spectrum preprocessing transformations include filtering on mass/charge and/or intensity values, weight factor transformation, low-entropy transformation, centroiding, noise removal, and matching. The available similarity measures include the canonical Cosine similarity measure, three entropy-based similarity measures, and a variety of binary similarity measures: Jaccard, Dice, 3W-Jaccard, Sokal-Sneath, Binary Cosine, Mountford, McConnaughey, Driver-Kroeber, Simpson, Braun-Banquet, Fager-McGowan, Kulczynski, Intersection, Hamming, and Hellinger.", style="margin-top:10px; text-align:left; font-size:16px; font-weight:500"),
                ui.div("Select options:", style="margin-top:30px; text-align:left; font-size:24px; font-weight:bold"),
                ui.div(ui.input_radio_buttons("chromatography_platform", "Specify chromatography platform:", ["HRMS","NRMS"]), style="font-size:18px; margin-top:10px; max-width:none"),
                ui.input_action_button("plot_spectra", "Plot two spectra before and after preprocessing transformations.", style="font-size:18px; padding:20px 40px; width:550px; height:100px; margin-top:10px; margin-right:50px"),
                ui.input_action_button("run_spec_lib_matching", "Run spectral library matching to perform compound identification on a query library of spectra.", style="font-size:18px; padding:20px 40px; width:550px; height:100px; margin-top:10px; margin-right:50px"),
                ui.input_action_button("run_parameter_tuning_grid", "Grid search: Tune parameters to maximize accuracy of compound identification given a query library with known spectrum IDs.", style="font-size:18px; padding:20px 40px; width:450px; height:120px; margin-top:10px; margin-right:50px"),
                ui.input_action_button("run_parameter_tuning_DE", "Differential evolution optimization: Tune parameters to maximize accuracy of compound identification given a query library with known spectrum IDs.", style="font-size:18px; padding:20px 40px; width:500px; height:150px; margin-top:10px; margin-right:50px"),
                ui.div(
                    "References:",
                    style="margin-top:35px; text-align:left; font-size:24px; font-weight:bold"
                ),
                ui.div(
                    "If Shannon Entropy similarity measure, low-entropy transformation, or centroiding are used:",
                    style="margin-top:10px; text-align:left; font-size:14px; font-weight:500"
                ),
                ui.div(
                    ui.HTML(
                        'Li, Y., Kind, T., Folz, J. et al. (2021) Spectral entropy outperforms MS/MS dot product similarity for small-molecule compound identification. Nat Methods, 18 1524–1531. <a href="https://doi.org/10.1038/s41592-021-01331-z" target="_blank">https://doi.org/10.1038/s41592-021-01331-z</a>.'
                    ),
                    style="text-align:left; font-size:14px; font-weight:500"
                ),
                ui.div(
                    "If Tsallis Entropy similarity measure or series of preprocessing transformations are used:",
                    style="margin-top:10px; text-align:left; font-size:14px; font-weight:500"
                ),
                ui.div(
                    ui.HTML(
                        'Dlugas, H., Zhang, X., Kim, S. (2025) Comparative analysis of continuous similarity measures for compound identification in mass spectrometry-based metabolomics. Chemometrics and Intelligent Laboratory Systems, 263, 105417. <a href="https://doi.org/10.1016/j.chemolab.2025.105417", target="_blank">https://doi.org/10.1016/j.chemolab.2025.105417</a>.'
                    ),
                    style="text-align:left; font-size:14px; font-weight:500"
                ),
                ui.div(
                    "If binary similarity measures are used:",
                    style="margin-top:10px; text-align:left; font-size:14px; font-weight:500"
                ),
                ui.div(
                    ui.HTML(
                        'Kim, S., Kato, I., & Zhang, X. (2022). Comparative Analysis of Binary Similarity Measures for Compound Identification in Mass Spectrometry-Based Metabolomics. Metabolites, 12(8), 694. <a href="https://doi.org/10.3390/metabo12080694" target="_blank">https://doi.org/10.3390/metabo12080694</a>.'
                    ),
                    style="text-align:left; font-size:14px; font-weight:500"
                ),

                ui.div(
                    "If weight factor transformation is used:",
                    style="margin-top:10px; text-align:left; font-size:14px; font-weight:500"
                ),
                ui.div(
                    ui.HTML(
                        'Kim, S., Koo, I., Wei, X., & Zhang, X. (2012). A method of finding optimal weight factors for compound identification in gas chromatography-mass spectrometry. Bioinformatics, 28(8), 1158-1163. <a href="https://doi.org/10.1093/bioinformatics/bts083" target="_blank">https://doi.org/10.1093/bioinformatics/bts083</a>.'
                    ),
                    style="margin-bottom:40px; text-align:left; font-size:14px; font-weight:500"
                ),
            )
        elif current_page() == "plot_spectra":
            return plot_spectra_ui(input.chromatography_platform())
        elif current_page() == "run_spec_lib_matching":
            return run_spec_lib_matching_ui(input.chromatography_platform())
        elif current_page() == "run_parameter_tuning_grid":
            return run_parameter_tuning_grid_ui(input.chromatography_platform())
        elif current_page() == "run_parameter_tuning_DE":
            return run_parameter_tuning_DE_ui(input.chromatography_platform())



    @reactive.effect
    @reactive.event(input.query_data)
    async def _populate_ids_from_query_upload():
        files = input.query_data()
        if not files:
            return

        in_path = Path(files[0]["datapath"])
        suffix = in_path.suffix.lower()

        try:
            if suffix == ".txt":
                txt_path = in_path
                converted_query_path_rv.set(str(txt_path))
            else:
                query_status_rv.set(f"Converting {in_path.name} →  TXT…")
                await reactive.flush()

                tmp_txt_path = in_path.with_suffix(".converted.txt")

                out_obj = await asyncio.to_thread(build_library, str(in_path), str(tmp_txt_path))

                if isinstance(out_obj, (str, os.PathLike, Path)):
                    txt_path = Path(out_obj)
                elif isinstance(out_obj, pd.DataFrame):
                    out_obj.to_csv(tmp_txt_path, index=False, sep='\t')
                    txt_path = tmp_txt_path
                else:
                    raise TypeError(f"build_library returned unsupported type: {type(out_obj)}")

                converted_query_path_rv.set(str(txt_path))

            query_status_rv.set(f"Reading IDs from: {txt_path.name} …")
            await reactive.flush()

            ids = await asyncio.to_thread(extract_first_column_ids, str(txt_path))
            query_ids_rv.set(ids)

            ui.update_selectize("spectrum_ID1", choices=ids, selected=(ids[0] if ids else None))

            query_status_rv.set(f"✅ Loaded {len(ids)} IDs from {txt_path.name}" if ids else f"⚠️ No IDs found in {txt_path.name}")
            await reactive.flush()

        except Exception as e:
            query_status_rv.set(f"❌ Failed: {e}")
            await reactive.flush()
            raise


    @reactive.effect
    @reactive.event(input.reference_data)
    async def _populate_ids_from_reference_upload():
        files = input.reference_data()
        if not files:
            return

        in_path = Path(files[0]["datapath"])
        suffix = in_path.suffix.lower()

        try:
            if suffix == ".txt":
                txt_path = in_path
                converted_reference_path_rv.set(str(txt_path))
            else:
                reference_status_rv.set(f"Converting {in_path.name} →  TXT…")
                await reactive.flush()

                tmp_txt_path = in_path.with_suffix(".converted.txt")

                out_obj = await asyncio.to_thread(build_library, str(in_path), str(tmp_txt_path))

                if isinstance(out_obj, (str, os.PathLike, Path)):
                    txt_path = Path(out_obj)
                elif isinstance(out_obj, pd.DataFrame):
                    out_obj.to_csv(tmp_txt_path, index=False, sep='\t')
                    txt_path = tmp_txt_path
                else:
                    raise TypeError(f"build_library returned unsupported type: {type(out_obj)}")

                converted_reference_path_rv.set(str(txt_path))

            reference_status_rv.set(f"Reading IDs from: {txt_path.name} …")
            await reactive.flush()

            ids = await asyncio.to_thread(extract_first_column_ids, str(txt_path))
            reference_ids_rv.set(ids)

            ui.update_selectize("spectrum_ID2", choices=ids, selected=(ids[0] if ids else None))

            reference_status_rv.set(
                f"✅ Loaded {len(ids)} IDs from {txt_path.name}" if ids else f"⚠️ No IDs found in {txt_path.name}"
            )
            await reactive.flush()

        except Exception as e:
            reference_status_rv.set(f"❌ Failed: {e}")
            await reactive.flush()
            raise


    @render.download(filename=lambda: f"plot.svg")
    def run_btn_plot_spectra():
        spectrum_ID1 = input.spectrum_ID1() or None
        spectrum_ID2 = input.spectrum_ID2() or None

        weights = [float(weight.strip()) for weight in input.weights().split(",") if weight.strip()]
        weights = {'Cosine':weights[0], 'Shannon':weights[1], 'Renyi':weights[2], 'Tsallis':weights[3]}

        high_quality_reference_library_tmp2 = False
        if input.high_quality_reference_library() != 'False':
            high_quality_reference_library_tmp2 = True

        if input.chromatography_platform() == "HRMS":
            fig = generate_plots_on_HRMS_data(query_data=input.query_data()[0]['datapath'], reference_data=input.reference_data()[0]['datapath'], spectrum_ID1=spectrum_ID1, spectrum_ID2=spectrum_ID2, print_url_spectrum1=input.print_url_spectrum1(), print_url_spectrum2=input.print_url_spectrum2(), similarity_measure=input.similarity_measure(), weights=weights, spectrum_preprocessing_order=input.spectrum_preprocessing_order(), high_quality_reference_library=high_quality_reference_library_tmp2, mz_min=input.mz_min(), mz_max=input.mz_max(), int_min=input.int_min(), int_max=input.int_max(), window_size_centroiding=input.window_size_centroiding(), window_size_matching=input.window_size_matching(), noise_threshold=input.noise_threshold(), wf_mz=input.wf_mz(), wf_intensity=input.wf_int(), LET_threshold=input.LET_threshold(), entropy_dimension=input.entropy_dimension(), y_axis_transformation=input.y_axis_transformation(), return_plot=True)
            plt.show()
        elif input.chromatography_platform() == "NRMS":
            fig = generate_plots_on_NRMS_data(query_data=input.query_data()[0]['datapath'], reference_data=input.reference_data()[0]['datapath'], spectrum_ID1=spectrum_ID1, spectrum_ID2=spectrum_ID2, print_url_spectrum1=input.print_url_spectrum1(), print_url_spectrum2=input.print_url_spectrum2(), similarity_measure=input.similarity_measure(), spectrum_preprocessing_order=input.spectrum_preprocessing_order(), high_quality_reference_library=high_quality_reference_library_tmp2, mz_min=input.mz_min(), mz_max=input.mz_max(), int_min=input.int_min(), int_max=input.int_max(), noise_threshold=input.noise_threshold(), wf_mz=input.wf_mz(), wf_intensity=input.wf_int(), LET_threshold=input.LET_threshold(), entropy_dimension=input.entropy_dimension(), y_axis_transformation=input.y_axis_transformation(), return_plot=True)
            plt.show()
        with io.BytesIO() as buf:
            fig.savefig(buf, format="svg", dpi=150, bbox_inches="tight")
            plt.close()
            yield buf.getvalue()




    @render.download(filename="identification_output.txt")
    async def run_btn_spec_lib_matching():
        match_log_rv.set("Running identification...\n")
        await reactive.flush()

        hq = input.high_quality_reference_library()
        if isinstance(hq, str):
            hq = hq.lower() == "true"
        elif isinstance(hq, (int, float)):
            hq = bool(hq)

        weights = [float(weight.strip()) for weight in input.weights().split(",") if weight.strip()]
        weights = {'Cosine': weights[0], 'Shannon': weights[1], 'Renyi': weights[2], 'Tsallis': weights[3]}

        common_kwargs = dict(
            query_data=input.query_data()[0]["datapath"],
            reference_data=input.reference_data()[0]["datapath"],
            likely_reference_ids=None,
            similarity_measure=input.similarity_measure(),
            weights=weights,
            spectrum_preprocessing_order=input.spectrum_preprocessing_order(),
            high_quality_reference_library=hq,
            mz_min=input.mz_min(), mz_max=input.mz_max(),
            int_min=input.int_min(), int_max=input.int_max(),
            noise_threshold=input.noise_threshold(),
            wf_mz=input.wf_mz(), wf_intensity=input.wf_int(),
            LET_threshold=input.LET_threshold(), entropy_dimension=input.entropy_dimension(),
            n_top_matches_to_save=input.n_top_matches_to_save(),
            print_id_results=True,
            output_identification=str(Path.cwd() / "identification_output.txt"),
            output_similarity_scores=str(Path.cwd() / "similarity_scores.txt"),
            return_ID_output=True,
        )

        # --- streaming setup (same pattern as your DE block) ---
        loop = asyncio.get_running_loop()
        q: asyncio.Queue[str | None] = asyncio.Queue()

        class UIWriter(io.TextIOBase):
            def write(self, s: str):
                if s:
                    loop.call_soon_threadsafe(q.put_nowait, s)
                return len(s)
            def flush(self): pass

        async def _drain():
            while True:
                msg = await q.get()
                if msg is None:
                    break
                match_log_rv.set(match_log_rv.get() + msg)
                await reactive.flush()

        drain_task = asyncio.create_task(_drain())
        writer = UIWriter()

        # --- worker wrappers that install redirects INSIDE the thread ---
        def _run_hrms():
            with redirect_stdout(writer), redirect_stderr(writer):
                # optional heartbeat
                print(">> Starting HRMS identification ...", flush=True)
                return run_spec_lib_matching_on_HRMS_data_shiny(
                    precursor_ion_mz_tolerance=input.precursor_ion_mz_tolerance(),
                    ionization_mode=input.ionization_mode(),
                    adduct=input.adduct(),
                    window_size_centroiding=input.window_size_centroiding(),
                    window_size_matching=input.window_size_matching(),
                    **common_kwargs
                )

        def _run_nrms():
            with redirect_stdout(writer), redirect_stderr(writer):
                print(">> Starting NRMS identification ...", flush=True)
                return run_spec_lib_matching_on_NRMS_data_shiny(**common_kwargs)

        # --- run in worker thread and stream output live ---
        try:
            if input.chromatography_platform() == "HRMS":
                df_out = await asyncio.to_thread(_run_hrms)
            else:
                df_out = await asyncio.to_thread(_run_nrms)

            match_log_rv.set(match_log_rv.get() + "\n✅ Identification finished.\n")
            await reactive.flush()

        except Exception as e:
            import traceback
            tb = "".join(traceback.format_exception(type(e), e, e.__traceback__))
            match_log_rv.set(match_log_rv.get() + f"\n❌ {type(e).__name__}: {e}\n{tb}\n")
            await reactive.flush()
            # make sure to stop the drainer before re-raising
            await q.put(None); await drain_task
            raise

        finally:
            await q.put(None)
            await drain_task

        yield df_out.to_csv(index=True, sep="\t")




    @render.download(filename="plot.svg")
    def run_btn_plot_spectra_within_spec_lib_matching():
        req(input.query_data(), input.reference_data())

        spectrum_ID1 = input.q_spec() or None
        spectrum_ID2 = input.r_spec() or None

        hq = input.high_quality_reference_library()
        if isinstance(hq, str):
            hq = hq.lower() == "true"
        elif isinstance(hq, (int, float)):
            hq = bool(hq)

        weights = [float(weight.strip()) for weight in input.weights().split(",") if weight.strip()]
        weights = {'Cosine':weights[0], 'Shannon':weights[1], 'Renyi':weights[2], 'Tsallis':weights[3]}

        common = dict(
            query_data=input.query_data()[0]['datapath'],
            reference_data=input.reference_data()[0]['datapath'],
            spectrum_ID1=spectrum_ID1,
            spectrum_ID2=spectrum_ID2,
            print_url_spectrum1=input.print_url_spectrum1(),
            print_url_spectrum2=input.print_url_spectrum2(), 
            similarity_measure=input.similarity_measure(),
            weights=weights,
            spectrum_preprocessing_order=input.spectrum_preprocessing_order(),
            high_quality_reference_library=hq,
            mz_min=input.mz_min(), mz_max=input.mz_max(),
            int_min=input.int_min(), int_max=input.int_max(),
            noise_threshold=input.noise_threshold(),
            wf_mz=input.wf_mz(), wf_intensity=input.wf_int(),
            LET_threshold=input.LET_threshold(), entropy_dimension=input.entropy_dimension(),
            y_axis_transformation="normalized",
            return_plot=True
        )

        if input.chromatography_platform() == "HRMS":
            fig = generate_plots_on_HRMS_data(
                window_size_centroiding=input.window_size_centroiding(),
                window_size_matching=input.window_size_matching(),
                **common
            )
            plt.show()
        else:
            fig = generate_plots_on_NRMS_data(**common)
            plt.show()

        with io.BytesIO() as buf:
            fig.savefig(buf, format="svg", dpi=150, bbox_inches="tight")
            plt.close()
            yield buf.getvalue()


    @render.download(filename="parameter_tuning_grid_output.txt")
    async def run_btn_parameter_tuning_grid():
        is_any_job_running.set(True)
        is_tuning_grid_running.set(True)
        match_log_rv.set("Running grid search of all parameters specified...\n")
        await reactive.flush()

        similarity_measure_tmp = list(input.similarity_measure())
        high_quality_reference_library_tmp = [x.strip().lower() == "true" for x in input.high_quality_reference_library().strip().strip("[]").split(",") if x.strip()]
        spectrum_preprocessing_order_tmp = strip_text(input.spectrum_preprocessing_order())
        mz_min_tmp = strip_numeric(input.mz_min())
        mz_max_tmp = strip_numeric(input.mz_max())
        int_min_tmp = strip_numeric(input.int_min())
        int_max_tmp = strip_numeric(input.int_max())
        noise_threshold_tmp = strip_numeric(input.noise_threshold())
        wf_mz_tmp = strip_numeric(input.wf_mz())
        wf_int_tmp = strip_numeric(input.wf_int())
        LET_threshold_tmp = strip_numeric(input.LET_threshold())
        entropy_dimension_tmp = strip_numeric(input.entropy_dimension())
        weights_tmp = strip_weights(input.weights())

        common_kwargs = dict(
            query_data=input.query_data()[0]["datapath"],
            reference_data=input.reference_data()[0]["datapath"],
            output_path=str(Path.cwd() / "parameter_tuning_grid_output.txt"),
            return_output=True,
        )

        loop = asyncio.get_running_loop()
        rw = ReactiveWriter(loop)

        try:
            if input.chromatography_platform() == "HRMS":
                precursor_ion_mz_tolerance = float(input.precursor_ion_mz_tolerance())
                ionization_mode = str(input.ionization_mode())
                adduct = str(input.adduct())
                window_size_centroiding_tmp = strip_numeric(input.window_size_centroiding())
                window_size_matching_tmp = strip_numeric(input.window_size_matching())
                grid = {
                    'similarity_measure': similarity_measure_tmp,
                    'weight': weights_tmp,
                    'spectrum_preprocessing_order': spectrum_preprocessing_order_tmp,
                    'mz_min': mz_min_tmp,
                    'mz_max': mz_max_tmp,
                    'int_min': int_min_tmp,
                    'int_max': int_max_tmp,
                    'noise_threshold': noise_threshold_tmp,
                    'wf_mz': wf_mz_tmp,
                    'wf_int': wf_int_tmp,
                    'LET_threshold': LET_threshold_tmp,
                    'entropy_dimension': entropy_dimension_tmp,
                    'high_quality_reference_library': high_quality_reference_library_tmp,
                    'window_size_centroiding': window_size_centroiding_tmp,
                    'window_size_matching': window_size_matching_tmp,
                }
                df_out = await asyncio.to_thread(_run_with_redirects, tune_params_on_HRMS_data_grid_shiny, rw, **common_kwargs, grid=grid, precursor_ion_mz_tolerance=precursor_ion_mz_tolerance, ionization_mode=ionization_mode, adduct=adduct)
            else:
                grid = {
                    'similarity_measure': similarity_measure_tmp,
                    'weight': weights_tmp,
                    'spectrum_preprocessing_order': spectrum_preprocessing_order_tmp,
                    'mz_min': mz_min_tmp,
                    'mz_max': mz_max_tmp,
                    'int_min': int_min_tmp,
                    'int_max': int_max_tmp,
                    'noise_threshold': noise_threshold_tmp,
                    'wf_mz': wf_mz_tmp,
                    'wf_int': wf_int_tmp,
                    'LET_threshold': LET_threshold_tmp,
                    'entropy_dimension': entropy_dimension_tmp,
                    'high_quality_reference_library': high_quality_reference_library_tmp,
                }
                df_out = await asyncio.to_thread(_run_with_redirects, tune_params_on_NRMS_data_grid_shiny, rw, **common_kwargs, grid=grid)

            match_log_rv.set(match_log_rv.get() + "\n✅ Parameter tuning finished.\n")
        except Exception as e:
            match_log_rv.set(match_log_rv.get() + f"\n❌ Error: {e}\n")
            raise
        finally:
            is_tuning_grid_running.set(False)
            is_any_job_running.set(False)
            await reactive.flush()

        yield df_out.to_csv(index=False, sep='\t').encode("utf-8")



    @reactive.effect
    @reactive.event(input.run_btn_parameter_tuning_DE)
    async def run_btn_parameter_tuning_DE():
        match_log_rv.set("Tuning specified continuous parameters using differential evolution...\n")
        is_any_job_running.set(True)
        is_tuning_DE_running.set(True)
        await reactive.flush()

        def _safe_float(v, default):
            try:
                if v is None:
                    return default
                return float(v)
            except Exception:
                return default

        def _iget(id, default=None):
            if id in input:
                try:
                    return input[id]()
                except SilentException:
                    return default
            return default

        loop = asyncio.get_running_loop()
        q: asyncio.Queue[str | None] = asyncio.Queue()

        class UIWriter(io.TextIOBase):
            def write(self, s: str):
                if s:
                    loop.call_soon_threadsafe(q.put_nowait, s)
                return len(s)
            def flush(self): pass

        async def _drain():
            while True:
                msg = await q.get()
                if msg is None:
                    break
                match_log_rv.set(match_log_rv.get() + msg)
                await reactive.flush()

        drain_task = asyncio.create_task(_drain())
        writer = UIWriter()

        try:
            qfile = _iget("query_data")[0]["datapath"]
            rfile = _iget("reference_data")[0]["datapath"]

            platform = _iget("chromatography_platform", "HRMS")
            sim = _iget("similarity_measure", "cosine")
            spro = _iget("spectrum_preprocessing_order", "FCNMWL")

            hq_raw = _iget("high_quality_reference_library", False)
            if isinstance(hq_raw, str):
                hq = hq_raw.lower() == "true"
            else:
                hq = bool(hq_raw)

            mz_min = _safe_float(_iget("mz_min", 0.0), 0.0)
            mz_max = _safe_float(_iget("mz_max", 99_999_999.0), 99_999_999.0)
            int_min = _safe_float(_iget("int_min", 0.0), 0.0)
            int_max = _safe_float(_iget("int_max", 999_999_999.0), 999_999_999.0)

            w_text = _iget("weights", "") or ""
            w_list = [float(w.strip()) for w in w_text.split(",") if w.strip()]
            w_list = (w_list + [0.0, 0.0, 0.0, 0.0])[:4]
            weights = {"Cosine": w_list[0], "Shannon": w_list[1], "Renyi": w_list[2], "Tsallis": w_list[3]}

            opt_params = tuple(_iget("params", ()) or ())
            bounds_dict = {}
            param_defaults = PARAMS_HRMS if platform == "HRMS" else PARAMS_NRMS
            for p in opt_params:
                lo = _safe_float(_iget(f"min_{p}", param_defaults.get(p, (0.0, 1.0))[0]),
                                 param_defaults.get(p, (0.0, 1.0))[0])
                hi = _safe_float(_iget(f"max_{p}", param_defaults.get(p, (0.0, 1.0))[1]),
                                 param_defaults.get(p, (0.0, 1.0))[1])
                if lo > hi:
                    lo, hi = hi, lo
                bounds_dict[p] = (lo, hi)

            defaults = {
                "window_size_centroiding": _safe_float(_iget("window_size_centroiding", 0.5), 0.5),
                "window_size_matching":    _safe_float(_iget("window_size_matching",    0.5), 0.5),
                "noise_threshold":         _safe_float(_iget("noise_threshold",         0.0), 0.0),
                "wf_mz":                   _safe_float(_iget("wf_mz",                   0.0), 0.0),
                "wf_int":                  _safe_float(_iget("wf_int",                  1.0), 1.0),
                "LET_threshold":           _safe_float(_iget("LET_threshold",           0.0), 0.0),
                "entropy_dimension":       _safe_float(_iget("entropy_dimension",       1.1), 1.1),
            }
            if platform == "NRMS":
                defaults.pop("window_size_centroiding", None)
                defaults.pop("window_size_matching", None)

        except Exception as e:
            import traceback
            tb = "".join(traceback.format_exception(type(e), e, e.__traceback__))
            match_log_rv.set(match_log_rv.get() + f"\n❌ Input snapshot failed:\n{tb}\n")
            is_tuning_DE_running.set(False); is_any_job_running.set(False)
            await q.put(None); await drain_task; await reactive.flush()
            return

        def _run():
            with redirect_stdout(writer), redirect_stderr(writer):
                return tune_params_DE(
                    query_data=qfile,
                    reference_data=rfile,
                    precursor_ion_mz_tolerance=float(input.precursor_ion_mz_tolerance()),
                    ionization_mode=input.ionization_mode(),
                    adduct=input.adduct(),
                    chromatography_platform=input.chromatography_platform(),
                    similarity_measure=sim,
                    weights=weights,
                    spectrum_preprocessing_order=spro,
                    mz_min=mz_min, mz_max=mz_max,
                    int_min=int_min, int_max=int_max,
                    high_quality_reference_library=hq,
                    optimize_params=list(opt_params),
                    param_bounds=bounds_dict,
                    default_params=defaults,
                    de_workers=1,
                    maxiters=input.max_iterations()
                )

        try:
            _ = await asyncio.to_thread(_run)
            match_log_rv.set(match_log_rv.get() + "\n✅ Differential evolution finished.\n")
        except Exception as e:
            import traceback
            tb = "".join(traceback.format_exception(type(e), e, e.__traceback__))
            match_log_rv.set(match_log_rv.get() + f"\n❌ {type(e).__name__}: {e}\n{tb}\n")
        finally:
            await q.put(None)
            await drain_task
            is_tuning_DE_running.set(False)
            is_any_job_running.set(False)
            await reactive.flush()


    @reactive.effect
    async def _pump_reactive_writer_logs():
        if not is_tuning_grid_running.get():
            return

        reactive.invalidate_later(0.1)
        msgs = _drain_queue_nowait(_LOG_QUEUE)
        if msgs:
            match_log_rv.set(match_log_rv.get() + "".join(msgs))
            await reactive.flush()


    @render.text
    def status_output():
        return run_status_plot_spectra.get()
        return run_status_spec_lib_matching.get()
        return run_status_parameter_tuning_grid.get()
        return run_status_parameter_tuning_DE.get()

    @output
    @render.text
    def run_log():
        return match_log_rv.get()


app = App(app_ui, server)


