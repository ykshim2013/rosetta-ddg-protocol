#!/usr/bin/env python3
"""
parse_ddg_results.py
Parse Rosetta cartesian_ddG output files

Usage:
    python parse_ddg_results.py mutfile.ddg
    python parse_ddg_results.py results_directory/
"""

import re
import sys
import numpy as np
from pathlib import Path


def parse_rosetta_ddg(ddg_file):
    """
    Parse Rosetta cartesian_ddg output file.

    Args:
        ddg_file: Path to mutfile.ddg or ddg_predictions.out

    Returns:
        dict: Statistics for WT and mutant energies, plus ddG
    """
    wt_energies = []
    mut_energies = []
    mut_name = None

    with open(ddg_file, 'r') as f:
        for line in f:
            # Skip non-data lines
            if not line.startswith('COMPLEX:'):
                continue

            parts = line.split()
            if len(parts) < 4:
                continue

            # Extract total energy (4th field)
            try:
                energy = float(parts[3])
            except (ValueError, IndexError):
                continue

            # Classify as WT or mutant
            if 'WT:' in line or 'WT_:' in line:
                wt_energies.append(energy)
            elif 'MUT_' in line or 'MUT:' in line:
                mut_energies.append(energy)
                # Extract mutation name
                if mut_name is None:
                    match = re.search(r'MUT_?(\d*[A-Z]{3}|\w+):', line)
                    if match:
                        mut_name = match.group(1)

    if not wt_energies or not mut_energies:
        return None

    # Calculate ddG for each iteration (paired comparison)
    n_pairs = min(len(wt_energies), len(mut_energies))
    ddg_values = np.array(mut_energies[:n_pairs]) - np.array(wt_energies[:n_pairs])

    return {
        'mutation': mut_name,
        'n_iterations': n_pairs,
        'wt_mean': np.mean(wt_energies[:n_pairs]),
        'wt_sd': np.std(wt_energies[:n_pairs], ddof=1) if n_pairs > 1 else 0.0,
        'mut_mean': np.mean(mut_energies[:n_pairs]),
        'mut_sd': np.std(mut_energies[:n_pairs], ddof=1) if n_pairs > 1 else 0.0,
        'ddg_mean': np.mean(ddg_values),
        'ddg_sd': np.std(ddg_values, ddof=1) if n_pairs > 1 else 0.0,
        'ddg_median': np.median(ddg_values),
        'ddg_iqr': np.percentile(ddg_values, 75) - np.percentile(ddg_values, 25),
        'ddg_min': np.min(ddg_values),
        'ddg_max': np.max(ddg_values),
        'ddg_values': ddg_values.tolist()
    }


def interpret_ddg(ddg_mean):
    """Interpret ddG value."""
    if ddg_mean > 2.0:
        return "DESTABILIZING (likely pathogenic)"
    elif ddg_mean > 1.0:
        return "MILDLY DESTABILIZING"
    elif ddg_mean < -2.0:
        return "STABILIZING (likely benign)"
    elif ddg_mean < -1.0:
        return "MILDLY STABILIZING"
    else:
        return "NEUTRAL (VUS)"


def print_results(results, filepath=None):
    """Print formatted results."""
    print("=" * 60)
    print("ROSETTA CARTESIAN DDG RESULTS")
    print("=" * 60)

    if filepath:
        print(f"File: {filepath}")
    if results['mutation']:
        print(f"Mutation: {results['mutation']}")
    print(f"Iterations: {results['n_iterations']}")
    print()

    print("Energy Values (REU):")
    print(f"  Wildtype: {results['wt_mean']:.2f} +/- {results['wt_sd']:.2f}")
    print(f"  Mutant:   {results['mut_mean']:.2f} +/- {results['mut_sd']:.2f}")
    print()

    print("Stability Change (ddG, kcal/mol):")
    print(f"  Mean +/- SD: {results['ddg_mean']:.2f} +/- {results['ddg_sd']:.2f}")
    print(f"  Median:      {results['ddg_median']:.2f}")
    print(f"  Range:       [{results['ddg_min']:.2f}, {results['ddg_max']:.2f}]")
    print()

    interpretation = interpret_ddg(results['ddg_mean'])
    print(f"Interpretation: {interpretation}")
    print("=" * 60)


def process_directory(dir_path):
    """Process all ddG files in a directory."""
    results_list = []

    # Find all .ddg files
    ddg_files = list(Path(dir_path).glob('**/*.ddg'))

    if not ddg_files:
        print(f"No .ddg files found in {dir_path}")
        return

    print(f"Found {len(ddg_files)} .ddg files")
    print()

    for ddg_file in sorted(ddg_files):
        results = parse_rosetta_ddg(ddg_file)
        if results:
            # Get variant name from directory
            variant = ddg_file.parent.name
            results['variant'] = variant
            results['file'] = str(ddg_file)
            results_list.append(results)

    if not results_list:
        print("No valid results found")
        return

    # Print summary table
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"{'Variant':<15} {'ddG Mean':>10} {'ddG SD':>10} {'N':>5} {'Interpretation':<25}")
    print("-" * 80)

    for r in results_list:
        interp = interpret_ddg(r['ddg_mean'])[:25]
        print(f"{r['variant']:<15} {r['ddg_mean']:>10.2f} {r['ddg_sd']:>10.2f} {r['n_iterations']:>5} {interp:<25}")

    print("=" * 80)

    # Save to CSV
    csv_file = Path(dir_path) / "ddg_summary.csv"
    with open(csv_file, 'w') as f:
        f.write("variant,mutation,ddg_mean,ddg_sd,ddg_median,n_iterations,interpretation\n")
        for r in results_list:
            interp = interpret_ddg(r['ddg_mean'])
            f.write(f"{r['variant']},{r.get('mutation','')},{r['ddg_mean']:.3f},"
                    f"{r['ddg_sd']:.3f},{r['ddg_median']:.3f},{r['n_iterations']},{interp}\n")

    print(f"\nSaved summary to: {csv_file}")


def main():
    if len(sys.argv) < 2:
        print("Usage:")
        print("  python parse_ddg_results.py mutfile.ddg      # Single file")
        print("  python parse_ddg_results.py results_dir/     # Directory")
        sys.exit(1)

    path = Path(sys.argv[1])

    if path.is_dir():
        process_directory(path)
    elif path.is_file():
        results = parse_rosetta_ddg(path)
        if results:
            print_results(results, filepath=path)
        else:
            print(f"ERROR: Could not parse {path}")
            sys.exit(1)
    else:
        print(f"ERROR: Path not found: {path}")
        sys.exit(1)


if __name__ == "__main__":
    main()
