# Example Files

This directory contains example files to help you get started.

## Files

### example_mutfile
Standard Rosetta mutation file format for a single K135R mutation.

**Format:**
```
total <number_of_mutations>
<number_of_residues_in_this_mutation>
<WT_residue> <position> <MUT_residue>
```

### variant_list.txt
Example variant list for batch processing with `03_batch_processing.sh`.
One variant per line in format: `{WT}{position}{MUT}` (e.g., K135R).

### example_output.ddg
Example output from Rosetta cartesian_ddG showing the expected format.
Use this to understand how to parse results with `parse_ddg_results.py`.

## Usage

```bash
# Single variant
./scripts/02_cartesian_ddg.sh relaxed_structure.pdb examples/example_mutfile output/

# Batch processing
./scripts/03_batch_processing.sh relaxed_structure.pdb examples/variant_list.txt results/
```
