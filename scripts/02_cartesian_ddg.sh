#!/bin/bash
#
# 02_cartesian_ddg.sh
# Rosetta cartesian_ddG calculation
#
# Usage: ./02_cartesian_ddg.sh relaxed.pdb mutfile [output_dir]
#

set -e

# ============================================================
# CONFIGURATION - Update ROSETTA_BIN to your installation
# ============================================================
ROSETTA_BIN="${ROSETTA_BIN:-/path/to/rosetta/source/bin}"
PLATFORM="${PLATFORM:-macosclangrelease}"  # or linuxgccrelease
N_ITERATIONS="${N_ITERATIONS:-10}"

# ============================================================
# INPUT VALIDATION
# ============================================================
INPUT_PDB=$1
MUTFILE=$2
OUTPUT_DIR=${3:-.}

if [ -z "$INPUT_PDB" ] || [ -z "$MUTFILE" ]; then
    echo "Usage: $0 <relaxed.pdb> <mutfile> [output_dir]"
    echo ""
    echo "Environment variables:"
    echo "  ROSETTA_BIN: Path to Rosetta bin directory"
    echo "  PLATFORM: Binary suffix (default: macosclangrelease)"
    echo "  N_ITERATIONS: Number of ddG iterations (default: 10)"
    exit 1
fi

if [ ! -f "$INPUT_PDB" ]; then
    echo "ERROR: Input structure not found: $INPUT_PDB"
    exit 1
fi

if [ ! -f "$MUTFILE" ]; then
    echo "ERROR: Mutation file not found: $MUTFILE"
    exit 1
fi

CARTESIAN_DDG="${ROSETTA_BIN}/cartesian_ddg.${PLATFORM}"

if [ ! -f "$CARTESIAN_DDG" ]; then
    echo "ERROR: Rosetta cartesian_ddg not found: $CARTESIAN_DDG"
    echo ""
    echo "Please set ROSETTA_BIN to your Rosetta installation:"
    echo "  export ROSETTA_BIN=/path/to/rosetta/source/bin"
    exit 1
fi

# Create output directory if needed
mkdir -p "$OUTPUT_DIR"

# ============================================================
# RUN CARTESIAN DDG
# ============================================================
echo "========================================================"
echo "ROSETTA CARTESIAN DDG"
echo "========================================================"
echo "Input:      $INPUT_PDB"
echo "Mutfile:    $MUTFILE"
echo "Output:     $OUTPUT_DIR"
echo "Iterations: $N_ITERATIONS"
echo "Binary:     $CARTESIAN_DDG"
echo "========================================================"
echo ""
echo "Mutation file contents:"
cat "$MUTFILE"
echo ""
echo "Starting cartesian_ddG..."
echo ""

$CARTESIAN_DDG \
  -s "$INPUT_PDB" \
  -ddg:mut_file "$MUTFILE" \
  -ddg:iterations "$N_ITERATIONS" \
  -fa_max_dis 9.0 \
  -ddg:cartesian \
  -score:weights ref2015_cart \
  -ddg:dump_pdbs false \
  -ddg:suppress_checkpointing true \
  -out:path:all "$OUTPUT_DIR"

# ============================================================
# CHECK OUTPUT
# ============================================================
DDG_FILE="${OUTPUT_DIR}/$(basename ${MUTFILE}).ddg"

if [ -f "$DDG_FILE" ]; then
    echo ""
    echo "========================================================"
    echo "CARTESIAN DDG COMPLETE"
    echo "========================================================"
    echo "Output: $DDG_FILE"
    echo ""
    echo "Next step: Parse results"
    echo "  python scripts/parse_ddg_results.py $DDG_FILE"
    echo "========================================================"
else
    # Try alternative filename
    DDG_FILE="${OUTPUT_DIR}/mutfile.ddg"
    if [ -f "$DDG_FILE" ]; then
        echo ""
        echo "========================================================"
        echo "CARTESIAN DDG COMPLETE"
        echo "========================================================"
        echo "Output: $DDG_FILE"
        echo "========================================================"
    else
        echo ""
        echo "ERROR: Output file not found"
        echo "Check for errors above."
        ls -la ${OUTPUT_DIR}/*.ddg 2>/dev/null || echo "No .ddg files found."
        exit 1
    fi
fi
