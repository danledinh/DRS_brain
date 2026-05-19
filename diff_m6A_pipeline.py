# Differential m6A Analysis 

import os
import pandas as pd
import numpy as np
import statsmodels.api as sm
from tqdm import tqdm
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Differential m6A Analysis Pipeline")
    parser.add_argument('--input_dir', type=str, default="./input_data", help='Directory containing input files')
    parser.add_argument('--output_file', type=str, required=True, help='Path for output results TSV file')
    parser.add_argument('--file_counts_a', type=str, default="CDS_m6A_crossed_a_matrix.tsv", help='Filename for m6A counts (a)')
    parser.add_argument('--file_counts_A', type=str, default="CDS_m6A_crossed_A_matrix.tsv", help='Filename for canonical A counts (A)')
    parser.add_argument('--min_valid_A', type=float, default=1e4, help='Minimum valid A count per gene (filter)')
    parser.add_argument('--groupA', type=str, nargs='+', required=True, help='List of sample column names for group A (space separated)')
    parser.add_argument('--groupB', type=str, nargs='+', required=True, help='List of sample column names for group B (space separated)')
    return parser.parse_args()

# --- Data Loading and Filtering Functions ---
def load_and_filter_counts(file_a, file_A, min_valid_A=1e4):
    """
    Load count matrices, filter by minimum coverage, and align samples/genes.
    Returns: (a_mat_piv, A_mat_piv)
    """
    mat_a = pd.read_csv(file_a, sep="\t")
    mat_A = pd.read_csv(file_A, sep="\t")

    # Pivot to matrix format
    a_mat_piv = mat_a.pivot(index="name", columns="prefix", values="count_a").fillna(0).astype(int)
    A_mat_piv = mat_A.pivot(index="name", columns="prefix", values="count_A").fillna(0).astype(int)

    # Align genes and samples
    valid_idx = a_mat_piv.index.intersection(A_mat_piv.index)
    valid_col = a_mat_piv.columns.intersection(A_mat_piv.columns)
    a_mat_piv = a_mat_piv.loc[valid_idx, valid_col]
    A_mat_piv = A_mat_piv.loc[valid_idx, valid_col]

    # Filter by minimum valid A
    mask = (A_mat_piv > min_valid_A).all(axis=1)
    mask = mask.reindex(a_mat_piv.index, fill_value=False)
    a_mat_piv = a_mat_piv.loc[mask]
    A_mat_piv = A_mat_piv.loc[mask]

    return a_mat_piv, A_mat_piv

# --- Differential Analysis Function ---

def run_differential_analysis(a_mat_piv, A_mat_piv, groupA, groupB, output_dir, n_cap=1e4):
    """
    Run binomial GLM for each gene to compare m6A levels between two user-defined groups.
    groupA, groupB: lists of sample column names for each group
    """
    sample_names = groupA + groupB
    obs_cond = np.array([0]*len(groupA) + [1]*len(groupB))
    n_samples = len(sample_names)

    def _fit_gene(gn):
        k_raw = a_mat_piv.loc[gn, sample_names].values.astype(float)
        n_raw = (a_mat_piv + A_mat_piv).loc[gn, sample_names].values.astype(float)
        scale = np.minimum(1.0, n_cap / n_raw)
        k = np.round(k_raw * scale).astype(int)
        n = np.round(n_raw * scale).astype(int)
        try:
            endog = np.column_stack([k, n - k])
            exog  = np.column_stack([np.ones(n_samples), obs_cond])
            result = sm.GLM(endog, exog, family=Binomial()).fit(scale="dev", disp=0)
            log_or = result.params[1]
            p_val  = result.pvalues[1]
            obs_rate    = k_raw / n_raw * 100
            mu_low_obs  = obs_rate[:len(groupA)].mean()
            mu_high_obs = obs_rate[len(groupA):].mean()
            log2fc      = np.log2((mu_high_obs + 1e-6) / (mu_low_obs + 1e-6))
            return {'gene': gn, 'Log_Odds_Ratio': log_or,
                    'log2FC': log2fc,
                    'mu_low': mu_low_obs, 'mu_high': mu_high_obs,
                    'delta_m6a': mu_high_obs - mu_low_obs,
                    'P_Value_GLM': p_val, 'Converged': result.converged}
        except Exception:
            return {'gene': gn, 'Log_Odds_Ratio': np.nan, 'log2FC': np.nan,
                    'mu_low': np.nan, 'mu_high': np.nan, 'delta_m6a': np.nan,
                    'P_Value_GLM': np.nan, 'Converged': False}

    # Sequential (non-parallel) for-loop for clarity
    results_list = []
    for gn in tqdm(a_mat_piv.index, desc="Analyzing genes"):
        results_list.append(_fit_gene(gn))
    results_df = pd.DataFrame(results_list)
    return results_df

# --- Main Pipeline ---


def main():
    args = parse_args()
    file_a = os.path.join(args.input_dir, args.file_counts_a)
    file_A = os.path.join(args.input_dir, args.file_counts_A)
    a_mat_piv, A_mat_piv = load_and_filter_counts(
        file_a, file_A,
        min_valid_A=args.min_valid_A
    )
    # Validate user-provided group columns
    missing_A = [col for col in args.groupA if col not in a_mat_piv.columns]
    missing_B = [col for col in args.groupB if col not in a_mat_piv.columns]
    if missing_A or missing_B:
        raise ValueError(f"Missing columns in data: groupA missing {missing_A}, groupB missing {missing_B}")
    print(f"Loaded {a_mat_piv.shape[0]} genes and {a_mat_piv.shape[1]} samples.")
    print(f"GroupA: {args.groupA}")
    print(f"GroupB: {args.groupB}")
    results_df = run_differential_analysis(a_mat_piv, A_mat_piv, args.groupA, args.groupB, n_cap=1e4)
    results_df.to_csv(args.output_file, sep="\t", index=False)
    print("Analysis complete. Results saved to:", args.output_file)

if __name__ == "__main__":
    main()
