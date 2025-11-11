
from pycompound.build_library import build_library_from_raw_data
from .processing import *
from .similarity_measures import *
import pandas as pd
from pathlib import Path
import json
from itertools import product
from joblib import Parallel, delayed
import csv
import sys, csv
from scipy.optimize import differential_evolution


def _vector_to_full_params(X, default_params, optimize_params):
    params = default_params.copy()
    for name, val in zip(optimize_params, X):
        params[name] = float(val)
    return params


def objective_function_HRMS(X, ctx):
    p = _vector_to_full_params(X, ctx["default_params"], ctx["optimize_params"])
    acc = get_acc_HRMS(
        ctx["df_query"],
        ctx["df_reference"],
        ctx["precursor_ion_mz_tolerance"],
        ctx["ionization_mode"], ctx["adduct"],
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
        ctx["df_query"], ctx["df_reference"], ctx['unique_query_ids'], ctx['unique_reference_ids'],
        ctx["similarity_measure"], ctx["weights"], ctx["spectrum_preprocessing_order"],
        ctx["mz_min"], ctx["mz_max"], ctx["int_min"], ctx["int_max"],
        p["noise_threshold"], p["wf_mz"], p["wf_int"], p["LET_threshold"], p["entropy_dimension"],
        ctx["high_quality_reference_library"],
        verbose=False
    )
    print(f"\nparams({ctx['optimize_params']}) = {np.array(X)}\naccuracy: {acc*100}%")
    return 1.0 - acc




def tune_params_DE(query_data=None, reference_data=None, chromatography_platform='HRMS', precursor_ion_mz_tolerance=None, ionization_mode=None, adduct=None, similarity_measure='cosine', weights=None, spectrum_preprocessing_order='CNMWL', mz_min=0, mz_max=999999999, int_min=0, int_max=999999999, high_quality_reference_library=False, optimize_params=["window_size_centroiding","window_size_matching","noise_threshold","wf_mz","wf_int","LET_threshold","entropy_dimension"], param_bounds={"window_size_centroiding":(0.0,0.5),"window_size_matching":(0.0,0.5),"noise_threshold":(0.0,0.25),"wf_mz":(0.0,5.0),"wf_int":(0.0,5.0),"LET_threshold":(0.0,5.0),"entropy_dimension":(1.0,3.0)}, default_params={"window_size_centroiding": 0.5, "window_size_matching":0.5, "noise_threshold":0.10, "wf_mz":0.0, "wf_int":1.0, "LET_threshold":0.0, "entropy_dimension":1.1}, maxiters=3, de_workers=1):

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

    if reference_data is None:
        print('\nError: No argument passed to the mandatory reference_data. Please pass the path to the TXT file of the reference data.')
        sys.exit()
    else:
        if isinstance(reference_data,str):
            df_reference = get_reference_df(reference_data=reference_data)
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
        unique_query_ids=unique_query_ids,
        unique_reference_ids=unique_reference_ids,
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
        result = differential_evolution(objective_function_HRMS, bounds=bounds, args=(ctx,), maxiter=maxiters, tol=0.0, workers=de_workers, seed=1, updating='deferred' if de_workers!=1 else 'immediate')
    else:
        result = differential_evolution(objective_function_NRMS, bounds=bounds, args=(ctx,), maxiter=maxiters, tol=0.0, workers=de_workers, seed=1, updating='deferred' if de_workers!=1 else 'immediate')

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
        verbose=False
    )

    return (
        acc, similarity_measure_tmp, json.dumps(weight), spectrum_preprocessing_order_tmp,
        mz_min_tmp, mz_max_tmp, int_min_tmp, int_max_tmp, noise_threshold_tmp, 
        wf_mz_tmp, wf_int_tmp, LET_threshold_tmp, entropy_dimension_tmp, high_quality_reference_library_tmp
    )



def tune_params_on_HRMS_data_grid(query_data=None, reference_data=None, precursor_ion_mz_tolerance=None, ionization_mode=None, adduct=None, grid=None, output_path=None, return_output=False):
    grid = {**default_HRMS_grid, **(grid or {})}
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
        unique_query_ids = df_query['id'].unique()

    if reference_data is None:
        print('\nError: No argument passed to the mandatory reference_data. Please pass the path to the TXT file of the reference data.')
        sys.exit()
    else:
        if isinstance(reference_data,str):
            df_reference = get_reference_df(reference_data=reference_data)
            unique_reference_ids = df_reference['id'].unique()
        else:
            dfs = []
            unique_reference_ids = []
            for f in reference_data:
                tmp = get_reference_df(reference_data=f)
                dfs.append(tmp)
                unique_reference_ids.extend(tmp['id'].unique())
            df_reference = pd.concat(dfs, axis=0, ignore_index=True)

    if 'ionization_mode' in df_reference.columns.tolist() and ionization_mode != 'N/A' and ionization_mode != None:
        df_reference = df_reference.loc[df_reference['ionization_mode']==ionization_mode].copy()
    if 'adduct' in df_reference.columns.tolist() and adduct != 'N/A' and adduct != None:
        df_reference = df_reference.loc[df_reference['adduct']==adduct].copy()
    unique_reference_ids_tmp2 = df_reference['id'].unique()

    print(f'\nNote that there are {len(unique_query_ids)} unique query spectra, {len(unique_reference_ids)} unique reference spectra, and {len(set(unique_query_ids) & set(unique_reference_ids_tmp2))} of the query and reference spectra IDs are in common.\n')

    if output_path is None:
        output_path = f'{Path.cwd()}/tuning_param_output.txt'
        print(f'Warning: since output_path=None, the output will be written to the current working directory: {output_path}')

    param_grid = product(similarity_measure, weight, spectrum_preprocessing_order, mz_min, mz_max, int_min, int_max, noise_threshold,
                         window_size_centroiding, window_size_matching, wf_mz, wf_int, LET_threshold, entropy_dimension, high_quality_reference_library)
    results = Parallel(n_jobs=-1, verbose=10)(delayed(_eval_one_HRMS)(df_query, df_reference, precursor_ion_mz_tolerance, ionization_mode, adduct,  *params) for params in param_grid)

    df_out = pd.DataFrame(results, columns=[
        'ACC','SIMILARITY.MEASURE','WEIGHT','SPECTRUM.PROCESSING.ORDER', 'MZ.MIN','MZ.MAX','INT.MIN','INT.MAX','NOISE.THRESHOLD',
        'WINDOW.SIZE.CENTROIDING','WINDOW.SIZE.MATCHING', 'WF.MZ','WF.INT','LET.THRESHOLD','ENTROPY.DIMENSION', 'HIGH.QUALITY.REFERENCE.LIBRARY'
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
    df_out.to_csv(output_path, index=False, sep='\t', quoting=csv.QUOTE_NONE)

    if return_output is False:
        df_out.to_csv(output_path, index=False, sep='\t', quoting=csv.QUOTE_NONE)
    else:
        return df_out



def tune_params_on_NRMS_data_grid(query_data=None, reference_data=None, grid=None, output_path=None, return_output=False):
    grid = {**default_NRMS_grid, **(grid or {})}
    for key, value in grid.items():
        globals()[key] = value

    if query_data is None:
        print('\nError: No argument passed to the mandatory query_data. Please pass the path to the CSV file of the query data.')
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
        print('\nError: No argument passed to the mandatory reference_data. Please pass the path to the CSV file of the reference data.')
        sys.exit()
    else:
        if isinstance(reference_data,str):
            df_reference = get_reference_df(reference_data=reference_data)
            unique_reference_ids = df_reference['id'].unique()
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

    df_out = pd.DataFrame(results, columns=['ACC','SIMILARITY.MEASURE','WEIGHT','SPECTRUM.PROCESSING.ORDER', 'MZ.MIN','MZ.MAX','INT.MIN','INT.MAX',
                                            'NOISE.THRESHOLD','WF.MZ','WF.INT','LET.THRESHOLD','ENTROPY.DIMENSION', 'HIGH.QUALITY.REFERENCE.LIBRARY'])
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



def run_spec_lib_matching_on_HRMS_data(query_data=None, reference_data=None, precursor_ion_mz_tolerance=None, ionization_mode=None, adduct=None, likely_reference_ids=None, similarity_measure='cosine', weights={'Cosine':0.25,'Shannon':0.25,'Renyi':0.25,'Tsallis':0.25}, spectrum_preprocessing_order='FCNMWL', high_quality_reference_library=False, mz_min=0, mz_max=9999999, int_min=0, int_max=9999999, window_size_centroiding=0.5, window_size_matching=0.5, noise_threshold=0.0, wf_mz=0.0, wf_intensity=1.0, LET_threshold=0.0, entropy_dimension=1.1, n_top_matches_to_save=1, print_id_results=False, output_identification=None, output_similarity_scores=None, return_ID_output=False, verbose=True):
    if query_data is None:
        print('\nError: No argument passed to the mandatory query_data. Please pass the path to the CSV file of the query data.')
        sys.exit()
    else:
        extension = query_data.rsplit('.',1)
        extension = extension[(len(extension)-1)]
        if extension == 'mgf' or extension == 'MGF' or extension == 'mzML' or extension == 'mzml' or extension == 'MZML' or extension == 'cdf' or extension == 'CDF' or extension == 'json' or extension == 'JSON' or extension == 'msp' or extension == 'MSP':
            output_path_tmp = query_data[:-3] + 'txt'
            build_library_from_raw_data(input_path=query_data, output_path=output_path_tmp, is_reference=False)
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

    if 'ionization_mode' in df_reference.columns.tolist() and ionization_mode != 'N/A' and ionization_mode != None:
        df_reference = df_reference.loc[df_reference['ionization_mode']==ionization_mode]
    if 'adduct' in df_reference.columns.tolist() and adduct != 'N/A' and adduct != None:
        df_reference = df_reference.loc[df_reference['adduct']==adduct]

    print(df_reference.loc[df_reference['id']=='Hectochlorin M+H'])

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
            #print('\nhere!!!!!!!!!!!!!!!')
            #print(r_spec)

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





def run_spec_lib_matching_on_NRMS_data(query_data=None, reference_data=None, likely_reference_ids=None, spectrum_preprocessing_order='FNLW', similarity_measure='cosine', weights={'Cosine':0.25,'Shannon':0.25,'Renyi':0.25,'Tsallis':0.25}, high_quality_reference_library=False, mz_min=0, mz_max=9999999, int_min=0, int_max=9999999, noise_threshold=0.0, wf_mz=0.0, wf_intensity=1.0, LET_threshold=0.0, entropy_dimension=1.1, n_top_matches_to_save=1, print_id_results=False, output_identification=None, output_similarity_scores=None, return_ID_output=False, verbose=True):
    if query_data is None:
        print('\nError: No argument passed to the mandatory query_data. Please pass the path to the CSV file of the query data.')
        sys.exit()
    else:
        extension = query_data.rsplit('.',1)
        extension = extension[(len(extension)-1)]
        if extension == 'mgf' or extension == 'MGF' or extension == 'mzML' or extension == 'mzml' or extension == 'MZML' or extension == 'cdf' or extension == 'CDF':
            output_path_tmp = query_data[:-3] + 'txt'
            build_library_from_raw_data(input_path=query_data, output_path=output_path_tmp, is_reference=False)
            df_query = pd.read_csv(output_path_tmp, sep='\t')
        if extension == 'txt' or extension == 'TXT':
            df_query = pd.read_csv(query_data, sep='\t')
        unique_query_ids = df_query.iloc[:,0].unique()

    if reference_data is None:
        print('\nError: No argument passed to the mandatory reference_data. Please pass the path to the CSV file of the reference data.')
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

