#!/usr/bin/env python3
"""
statistical_analysis.py
Statistical analysis for Rosetta ddG results

Usage:
    python statistical_analysis.py results.csv --pathogenic pathogenic_variants.txt --benign benign_variants.txt
"""

import argparse
import numpy as np
import pandas as pd
from scipy import stats


def load_data(csv_file):
    """Load ddG results from CSV."""
    return pd.read_csv(csv_file)


def compare_groups(pathogenic_ddg, benign_ddg):
    """
    Compare ddG distributions between pathogenic and benign variants.

    Returns dict with statistical test results.
    """
    pathogenic_ddg = np.array(pathogenic_ddg)
    benign_ddg = np.array(benign_ddg)

    n1, n2 = len(pathogenic_ddg), len(benign_ddg)

    # Welch's t-test
    t_stat, t_pvalue = stats.ttest_ind(pathogenic_ddg, benign_ddg, equal_var=False)

    # Mann-Whitney U test
    u_stat, u_pvalue = stats.mannwhitneyu(
        pathogenic_ddg, benign_ddg, alternative='two-sided'
    )

    # Cohen's d effect size
    var1 = np.var(pathogenic_ddg, ddof=1)
    var2 = np.var(benign_ddg, ddof=1)
    pooled_sd = np.sqrt(((n1-1)*var1 + (n2-1)*var2) / (n1+n2-2))
    cohens_d = (np.mean(pathogenic_ddg) - np.mean(benign_ddg)) / pooled_sd if pooled_sd > 0 else 0

    # ROC-AUC (requires sklearn)
    try:
        from sklearn.metrics import roc_curve, auc
        y_true = np.concatenate([np.ones(n1), np.zeros(n2)])
        y_score = np.concatenate([pathogenic_ddg, benign_ddg])
        fpr, tpr, thresholds = roc_curve(y_true, y_score)
        roc_auc = auc(fpr, tpr)

        # Optimal threshold (Youden's J)
        j_scores = tpr - fpr
        optimal_idx = np.argmax(j_scores)
        optimal_threshold = thresholds[optimal_idx]
    except ImportError:
        roc_auc = None
        optimal_threshold = None

    return {
        'pathogenic_n': n1,
        'pathogenic_mean': np.mean(pathogenic_ddg),
        'pathogenic_sd': np.std(pathogenic_ddg, ddof=1),
        'pathogenic_median': np.median(pathogenic_ddg),
        'benign_n': n2,
        'benign_mean': np.mean(benign_ddg),
        'benign_sd': np.std(benign_ddg, ddof=1),
        'benign_median': np.median(benign_ddg),
        't_statistic': t_stat,
        't_pvalue': t_pvalue,
        'mann_whitney_u': u_stat,
        'mann_whitney_pvalue': u_pvalue,
        'cohens_d': cohens_d,
        'roc_auc': roc_auc,
        'optimal_threshold': optimal_threshold
    }


def print_comparison(results):
    """Print comparison results."""
    print("=" * 70)
    print("STATISTICAL COMPARISON: Pathogenic vs Benign")
    print("=" * 70)
    print()

    print("GROUP STATISTICS:")
    print(f"  Pathogenic (n={results['pathogenic_n']}):")
    print(f"    Mean +/- SD: {results['pathogenic_mean']:.2f} +/- {results['pathogenic_sd']:.2f} kcal/mol")
    print(f"    Median:      {results['pathogenic_median']:.2f} kcal/mol")
    print()
    print(f"  Benign (n={results['benign_n']}):")
    print(f"    Mean +/- SD: {results['benign_mean']:.2f} +/- {results['benign_sd']:.2f} kcal/mol")
    print(f"    Median:      {results['benign_median']:.2f} kcal/mol")
    print()

    print("STATISTICAL TESTS:")
    print(f"  Welch's t-test:    t = {results['t_statistic']:.3f}, p = {results['t_pvalue']:.2e}")
    print(f"  Mann-Whitney U:    U = {results['mann_whitney_u']:.0f}, p = {results['mann_whitney_pvalue']:.2e}")
    print()

    print("EFFECT SIZE:")
    print(f"  Cohen's d: {results['cohens_d']:.3f}")
    if abs(results['cohens_d']) < 0.2:
        effect_interp = "negligible"
    elif abs(results['cohens_d']) < 0.5:
        effect_interp = "small"
    elif abs(results['cohens_d']) < 0.8:
        effect_interp = "medium"
    else:
        effect_interp = "large"
    print(f"  Interpretation: {effect_interp}")
    print()

    if results['roc_auc'] is not None:
        print("CLASSIFICATION PERFORMANCE:")
        print(f"  ROC-AUC: {results['roc_auc']:.3f}")
        print(f"  Optimal threshold: {results['optimal_threshold']:.2f} kcal/mol")

    print("=" * 70)


def aggregate_by_variant(df, groupby_col='variant'):
    """Aggregate results across multiple seeds/replicates."""
    agg = df.groupby(groupby_col).agg({
        'ddg_mean': ['mean', 'std', 'count']
    }).reset_index()

    agg.columns = [groupby_col, 'ddg_aggregate_mean', 'ddg_between_seed_sd', 'n_seeds']
    return agg


def main():
    parser = argparse.ArgumentParser(description='Statistical analysis of Rosetta ddG results')
    parser.add_argument('csv_file', help='CSV file with ddG results')
    parser.add_argument('--pathogenic', help='File listing pathogenic variants (one per line)')
    parser.add_argument('--benign', help='File listing benign variants (one per line)')
    parser.add_argument('--classification-col', default=None, help='Column name for classification')

    args = parser.parse_args()

    # Load data
    df = load_data(args.csv_file)
    print(f"Loaded {len(df)} records from {args.csv_file}")
    print()

    # Determine classification
    if args.pathogenic and args.benign:
        # Load from files
        with open(args.pathogenic) as f:
            pathogenic_variants = set(line.strip() for line in f if line.strip())
        with open(args.benign) as f:
            benign_variants = set(line.strip() for line in f if line.strip())

        pathogenic_ddg = df[df['variant'].isin(pathogenic_variants)]['ddg_mean'].values
        benign_ddg = df[df['variant'].isin(benign_variants)]['ddg_mean'].values

    elif args.classification_col and args.classification_col in df.columns:
        # Use existing column
        pathogenic_ddg = df[df[args.classification_col] == 'pathogenic']['ddg_mean'].values
        benign_ddg = df[df[args.classification_col] == 'benign']['ddg_mean'].values

    elif 'classification' in df.columns:
        # Default classification column
        pathogenic_ddg = df[df['classification'] == 'pathogenic']['ddg_mean'].values
        benign_ddg = df[df['classification'] == 'benign']['ddg_mean'].values

    else:
        print("ERROR: No classification information provided")
        print("Use --pathogenic and --benign to specify variant lists")
        print("Or ensure CSV has a 'classification' column")
        return

    if len(pathogenic_ddg) == 0 or len(benign_ddg) == 0:
        print("ERROR: One or both groups have no data")
        print(f"  Pathogenic: {len(pathogenic_ddg)} variants")
        print(f"  Benign: {len(benign_ddg)} variants")
        return

    # Run comparison
    results = compare_groups(pathogenic_ddg, benign_ddg)
    print_comparison(results)


if __name__ == "__main__":
    main()
