
#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Any
import numpy as np
import pandas as pd
from scipy.optimize import differential_evolution
from pycompound.spec_lib_matching import get_acc_HRMS, get_acc_NRMS


ALL_PARAMS = [
    "window_size_centroiding",
    "window_size_matching",
    "noise_threshold",
    "wf_mz",
    "wf_int",
    "LET_threshold",
    "entropy_dimension",
]

SUGGESTED_BOUNDS = {
    "window_size_centroiding": (0.0, 0.5),
    "window_size_matching": (0.0, 0.5),
    "noise_threshold": (0.0, 0.25),
    "wf_mz": (0.0, 5.0),
    "wf_int": (0.0, 5.0),
    "LET_threshold": (0.0, 5.0),
    "entropy_dimension": (1.0, 3.0),
}

DEFAULT_PARAMS = {
    "window_size_centroiding": 0.5,
    "window_size_matching": 0.5,
    "noise_threshold": 0.10,
    "wf_mz": 0.0,
    "wf_int": 1.0,
    "LET_threshold": 0.0,
    "entropy_dimension": 1.1,
}


def parse_bound(s: str) -> Tuple[str, Tuple[float, float]]:
    """Parse a command-line bound specification.

    Args:
        s: Bound string in the form "name=min:max".

    Returns:
        Tuple containing the parameter name and its lower and upper bounds.

    Raises:
        argparse.ArgumentTypeError: If the bound string is malformed.
    """
    if "=" not in s or ":" not in s:
        raise argparse.ArgumentTypeError(f"Bad --bound format '{s}'. Use name=min:max")

    name, rng = s.split("=", 1)
    lo, hi = rng.split(":", 1)

    try:
        lo_f, hi_f = float(lo), float(hi)
    except ValueError as e:
        raise argparse.ArgumentTypeError(f"Non-numeric bound in '{s}': {e}")

    if lo_f > hi_f:
        raise argparse.ArgumentTypeError(f"Lower bound > upper bound in '{s}'")

    return name.strip(), (lo_f, hi_f)


def parse_default(s: str) -> Tuple[str, float]:
    """Parse a command-line default parameter specification.

    Args:
        s: Default parameter string in the form "name=value".

    Returns:
        Tuple containing the parameter name and default value.

    Raises:
        argparse.ArgumentTypeError: If the default string is malformed.
    """
    if "=" not in s:
        raise argparse.ArgumentTypeError(f"Bad --default format '{s}'. Use name=value")

    name, val = s.split("=", 1)

    try:
        v = float(val)
    except ValueError as e:
        raise argparse.ArgumentTypeError(f"Non-numeric default in '{s}': {e}")

    return name.strip(), v


def _vector_to_full_params(
    x: np.ndarray,
    default_params: Dict[str, float],
    optimize_params: List[str],
) -> Dict[str, float]:
    """Convert an optimization vector into a full parameter dictionary.

    Args:
        x: Vector of optimized parameter values.
        default_params: Dictionary containing the full default parameter set.
        optimize_params: Parameter names corresponding to values in x.

    Returns:
        A copy of default_params with optimized values inserted.
    """
    params = dict(default_params)

    for name, val in zip(optimize_params, x):
        params[name] = float(val)

    return params


def objective_HRMS(x: np.ndarray, ctx: Dict[str, Any]) -> float:
    """Evaluate the HRMS differential evolution objective function.

    Args:
        x: Vector of candidate parameter values.
        ctx: Context dictionary containing data, fixed settings, and optimization metadata.

    Returns:
        Objective value equal to one minus HRMS identification accuracy.
    """
    p = _vector_to_full_params(
        x=x,
        default_params=ctx["default_params"],
        optimize_params=ctx["optimize_params"],
    )

    acc = get_acc_HRMS(
        ctx["df_query"],
        ctx["df_reference"],
        ctx["precursor_ion_mz_tolerance"],
        ctx["ionization_mode"],
        ctx["adduct"],
        ctx["similarity_measure"],
        ctx["weights"],
        ctx["spectrum_preprocessing_order"],
        ctx["mz_min"],
        ctx["mz_max"],
        ctx["int_min"],
        ctx["int_max"],
        p["window_size_centroiding"],
        p["window_size_matching"],
        p["noise_threshold"],
        p["wf_mz"],
        p["wf_int"],
        p["LET_threshold"],
        p["entropy_dimension"],
        ctx["high_quality_reference_library"],
        verbose=False,
    )

    print(f"\n{ctx['optimize_params']} = {np.array(x)}\naccuracy: {acc * 100}%")

    return 1.0 - acc


def objective_NRMS(x: np.ndarray, ctx: Dict[str, Any]) -> float:
    """Evaluate the NRMS differential evolution objective function.

    Args:
        x: Vector of candidate parameter values.
        ctx: Context dictionary containing data, fixed settings, and optimization metadata.

    Returns:
        Objective value equal to one minus NRMS identification accuracy.
    """
    p = _vector_to_full_params(
        x=x,
        default_params=ctx["default_params"],
        optimize_params=ctx["optimize_params"],
    )

    acc = get_acc_NRMS(
        ctx["df_query"],
        ctx["df_reference"],
        ctx["uq"],
        ctx["ur"],
        ctx["similarity_measure"],
        ctx["weights"],
        ctx["spectrum_preprocessing_order"],
        ctx["mz_min"],
        ctx["mz_max"],
        ctx["int_min"],
        ctx["int_max"],
        p["noise_threshold"],
        p["wf_mz"],
        p["wf_int"],
        p["LET_threshold"],
        p["entropy_dimension"],
        ctx["high_quality_reference_library"],
    )

    print(f"\n{ctx['optimize_params']} = {np.array(x)}\naccuracy: {acc * 100}%")

    return 1.0 - acc


def main() -> None:
    """Run the differential evolution tuning command-line interface.

    Args:
        None.

    Returns:
        None.
    """
    parser = argparse.ArgumentParser(
        description="Parameter tuning via Differential Evolution for HRMS/NRMS using pycompound."
    )

    parser.add_argument(
        "--chromatography_platform",
        choices=["HRMS", "NRMS"],
        default="HRMS",
        help="Chromatography platform.",
    )
    parser.add_argument(
        "--query_data",
        required=True,
        help="Path to query TXT file. Must contain an 'id' column.",
    )
    parser.add_argument(
        "--reference_data",
        required=True,
        nargs="+",
        help="Path(s) to reference TXT file(s). Must contain an 'id' column.",
    )
    parser.add_argument(
        "--precursor_ion_mz_tolerance",
        type=float,
        default=None,
        help="Precursor ion m/z tolerance. Only applicable to HRMS. Default: None.",
    )
    parser.add_argument(
        "--ionization_mode",
        choices=["Positive", "Negative", None],
        default=None,
        help="Ionization mode. Only applicable to HRMS.",
    )
    parser.add_argument(
        "--adduct",
        choices=["H", "NH3", "NH4", "OH", "Cl", "K", "Li", "Na", None],
        default=None,
        help="Adduct. Only applicable to HRMS.",
    )
    parser.add_argument(
        "--similarity_measure",
        default="cosine",
        choices=["cosine", "shannon", "renyi", "tsallis"],
        help="Similarity measure.",
    )
    parser.add_argument(
        "--weights",
        default="",
        help="Weights specification. Empty string means None.",
    )
    parser.add_argument(
        "--spectrum_preprocessing_order",
        default="CNMWL",
        help="Spectrum preprocessing order string.",
    )
    parser.add_argument("--mz-min", type=float, default=0.0)
    parser.add_argument("--mz-max", type=float, default=999_999_999.0)
    parser.add_argument("--int-min", type=float, default=0.0)
    parser.add_argument("--int-max", type=float, default=999_999_999.0)
    parser.add_argument(
        "--hq-ref-lib",
        action="store_true",
        help="Use high-quality reference library flag.",
    )
    parser.add_argument(
        "--opt",
        nargs="+",
        default=["window_size_centroiding", "noise_threshold", "wf_mz", "wf_int"],
        help=f"Parameters to optimize. Must be a subset of {ALL_PARAMS}.",
    )
    parser.add_argument(
        "--bound",
        action="append",
        default=[],
        type=parse_bound,
        help="Bound specification in the form 'name=min:max'. Repeatable.",
    )
    parser.add_argument(
        "--default",
        dest="defaults",
        action="append",
        default=[],
        type=parse_default,
        help="Override a default parameter using 'name=value'. Repeatable.",
    )
    parser.add_argument("--maxiter", type=int, default=15)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument(
        "--workers",
        type=int,
        default=-1,
        help="Use -1 for all cores; use 1 to disable parallelism.",
    )

    args = parser.parse_args()

    unknown = [x for x in args.opt if x not in ALL_PARAMS]
    if unknown:
        sys.exit(f"Error: unknown --opt params: {unknown}")

    qpath = Path(args.query_data)
    if not qpath.exists():
        sys.exit(f"Query TXT not found: {qpath}")

    df_query = pd.read_csv(qpath, sep="\t")
    if "id" not in df_query.columns:
        sys.exit("Query TXT must contain an 'id' column.")

    ref_paths = [Path(pth) for pth in args.reference_data]
    for r in ref_paths:
        if not r.exists():
            sys.exit(f"Reference TXT not found: {r}")

    df_reference = pd.concat(
        [pd.read_csv(r, sep="\t") for r in ref_paths],
        axis=0,
        ignore_index=True,
    )

    if "id" not in df_reference.columns:
        sys.exit("Reference TXT must contain an 'id' column.")

    uq = df_query["id"].unique().tolist()
    ur = df_reference["id"].unique().tolist()

    default_params = dict(DEFAULT_PARAMS)
    for name, val in args.defaults:
        if name not in DEFAULT_PARAMS:
            sys.exit(
                f"--default refers to unknown parameter '{name}'. "
                f"Allowed: {list(DEFAULT_PARAMS)}"
            )
        default_params[name] = val

    param_bounds: Dict[str, Tuple[float, float]] = dict(SUGGESTED_BOUNDS)
    for name, (lo, hi) in args.bound:
        if name not in SUGGESTED_BOUNDS:
            sys.exit(
                f"--bound refers to unknown parameter '{name}'. "
                f"Allowed: {list(SUGGESTED_BOUNDS)}"
            )
        param_bounds[name] = (lo, hi)

    bounds = [param_bounds[p] for p in args.opt]

    ctx: Dict[str, Any] = {
        "df_query": df_query,
        "df_reference": df_reference,
        "precursor_ion_mz_tolerance": args.precursor_ion_mz_tolerance,
        "ionization_mode": args.ionization_mode,
        "uq": uq,
        "ur": ur,
        "adduct": args.adduct,
        "similarity_measure": args.similarity_measure,
        "weights": None if args.weights.strip() == "" else args.weights,
        "spectrum_preprocessing_order": args.spectrum_preprocessing_order,
        "mz_min": float(args.mz_min),
        "mz_max": float(args.mz_max),
        "int_min": float(args.int_min),
        "int_max": float(args.int_max),
        "high_quality_reference_library": bool(args.hq_ref_lib),
        "default_params": default_params,
        "optimize_params": args.opt,
    }

    history_acc: List[float] = []

    def callback(xk: np.ndarray, convergence: float) -> None:
        """Record accuracy at each differential evolution callback.

        Args:
            xk: Current best parameter vector.
            convergence: Current convergence metric from differential evolution.

        Returns:
            None.
        """
        if args.chromatography_platform == "HRMS":
            acc_pct = (1.0 - objective_HRMS(xk, ctx)) * 100.0
        else:
            acc_pct = (1.0 - objective_NRMS(xk, ctx)) * 100.0

        history_acc.append(acc_pct)

    objective = objective_HRMS if args.chromatography_platform == "HRMS" else objective_NRMS

    result = differential_evolution(
        objective,
        bounds=bounds,
        args=(ctx,),
        maxiter=int(args.maxiter),
        tol=0.0,
        seed=int(args.seed),
        workers=int(args.workers),
        callback=callback,
        updating="deferred" if int(args.workers) != 1 else "immediate",
    )

    best_params = _vector_to_full_params(
        x=result.x,
        default_params=default_params,
        optimize_params=args.opt,
    )
    best_acc_pct = (1.0 - result.fun) * 100.0

    print("\n=== Differential Evolution Result ===")
    print(f"Mode: {args.chromatography_platform}")
    print(f"Optimized over: {args.opt}")
    print("Best values for selected parameters:")

    for name in args.opt:
        print(f"  {name}: {best_params[name]}")

    print("\nFull parameter set used in final evaluation:")

    for k in ALL_PARAMS:
        print(f"  {k}: {best_params[k]}")

    print(f"\nBest accuracy: {best_acc_pct:.3f}%")
    print(
        f"DE raw: success={result.success}, "
        f"nfev={result.nfev}, "
        f"nit={result.nit}, "
        f"message='{result.message}'"
    )


if __name__ == "__main__":
    main()

