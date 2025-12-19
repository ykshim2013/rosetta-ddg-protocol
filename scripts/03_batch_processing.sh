#!/bin/bash
#
# 03_batch_processing.sh
# Process multiple variants in batch
#
# Usage: ./03_batch_processing.sh relaxed.pdb variant_list.txt [output_dir]
#

set -e

# ============================================================
# CONFIGURATION
# ============================================================
ROSETTA_BIN="${ROSETTA_BIN:-/path/to/rosetta/source/bin}"
PLATFORM="${PLATFORM:-macosclangrelease}"
N_ITERATIONS="${N_ITERATIONS:-10}"
N_PARALLEL="${N_PARALLEL:-1}"  # Parallel jobs (set >1 for multi-core)

# ============================================================
# INPUT VALIDATION
# ============================================================
INPUT_PDB=$1
VARIANT_LIST=$2
OUTPUT_BASE=${3:-"rosetta_results"}

if [ -z "$INPUT_PDB" ] || [ -z "$VARIANT_LIST" ]; then
    echo "Usage: $0 <relaxed.pdb> <variant_list.txt> [output_dir]"
    echo ""
    echo "variant_list.txt format (one variant per line):"
    echo "  K135R"
    echo "  D175G"
    echo "  R138C"
    echo ""
    echo "Environment variables:"
    echo "  ROSETTA_BIN: Path to Rosetta bin directory"
    echo "  N_ITERATIONS: Iterations per variant (default: 10)"
    echo "  N_PARALLEL: Parallel jobs (default: 1)"
    exit 1
fi

if [ ! -f "$INPUT_PDB" ]; then
    echo "ERROR: Input structure not found: $INPUT_PDB"
    exit 1
fi

if [ ! -f "$VARIANT_LIST" ]; then
    echo "ERROR: Variant list not found: $VARIANT_LIST"
    exit 1
fi

CARTESIAN_DDG="${ROSETTA_BIN}/cartesian_ddg.${PLATFORM}"

if [ ! -f "$CARTESIAN_DDG" ]; then
    echo "ERROR: Rosetta cartesian_ddg not found: $CARTESIAN_DDG"
    exit 1
fi

# ============================================================
# SETUP
# ============================================================
mkdir -p "$OUTPUT_BASE"

TOTAL_VARIANTS=$(wc -l < "$VARIANT_LIST" | tr -d ' ')

echo "========================================================"
echo "BATCH ROSETTA CARTESIAN DDG"
echo "========================================================"
echo "Input structure: $INPUT_PDB"
echo "Variant list:    $VARIANT_LIST"
echo "Total variants:  $TOTAL_VARIANTS"
echo "Iterations:      $N_ITERATIONS"
echo "Parallel jobs:   $N_PARALLEL"
echo "Output:          $OUTPUT_BASE/"
echo "========================================================"
echo ""

# ============================================================
# PROCESS FUNCTION
# ============================================================
process_variant() {
    local variant=$1
    local input_pdb=$2
    local output_base=$3
    local iterations=$4
    local cartesian_ddg=$5

    # Parse variant (e.g., K135R -> K, 135, R)
    local wt_res=${variant:0:1}
    local position=$(echo "$variant" | sed 's/^[A-Z]//' | sed 's/[A-Z]$//')
    local mut_res=${variant: -1}

    # Create variant directory
    local variant_dir="${output_base}/${variant}"
    mkdir -p "$variant_dir"

    # Create mutation file
    cat > "${variant_dir}/mutfile" << EOF
total 1
1
$wt_res $position $mut_res
EOF

    # Run cartesian_ddg
    $cartesian_ddg \
      -s "$input_pdb" \
      -ddg:mut_file "${variant_dir}/mutfile" \
      -ddg:iterations $iterations \
      -fa_max_dis 9.0 \
      -ddg:cartesian \
      -score:weights ref2015_cart \
      -ddg:dump_pdbs false \
      -ddg:suppress_checkpointing true \
      -out:path:all "$variant_dir" \
      > "${variant_dir}/rosetta.log" 2>&1

    if [ $? -eq 0 ]; then
        echo "Completed: $variant"
    else
        echo "FAILED: $variant (check ${variant_dir}/rosetta.log)"
    fi
}

export -f process_variant

# ============================================================
# RUN BATCH
# ============================================================
START_TIME=$(date +%s)

if [ "$N_PARALLEL" -gt 1 ]; then
    # Parallel execution
    cat "$VARIANT_LIST" | xargs -P $N_PARALLEL -I {} \
        bash -c "process_variant {} '$INPUT_PDB' '$OUTPUT_BASE' '$N_ITERATIONS' '$CARTESIAN_DDG'"
else
    # Sequential execution with progress
    CURRENT=0
    while read variant; do
        CURRENT=$((CURRENT + 1))
        echo "[$CURRENT/$TOTAL_VARIANTS] Processing: $variant"
        process_variant "$variant" "$INPUT_PDB" "$OUTPUT_BASE" "$N_ITERATIONS" "$CARTESIAN_DDG"
    done < "$VARIANT_LIST"
fi

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
HOURS=$((ELAPSED / 3600))
MINUTES=$(((ELAPSED % 3600) / 60))

# ============================================================
# SUMMARY
# ============================================================
echo ""
echo "========================================================"
echo "BATCH PROCESSING COMPLETE"
echo "========================================================"
echo "Total time: ${HOURS}h ${MINUTES}m"
echo "Results in: $OUTPUT_BASE/"
echo ""
echo "Next step: Aggregate results"
echo "  python scripts/parse_ddg_results.py $OUTPUT_BASE"
echo "========================================================"
