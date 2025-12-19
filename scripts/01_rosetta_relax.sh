#!/bin/bash
#
# 01_rosetta_relax.sh
# Rosetta constrained relaxation for structure preparation
#
# Usage: ./01_rosetta_relax.sh input.pdb [output_prefix]
#

set -e

# ============================================================
# CONFIGURATION - Update ROSETTA_BIN to your installation
# ============================================================
ROSETTA_BIN="${ROSETTA_BIN:-/path/to/rosetta/source/bin}"
PLATFORM="${PLATFORM:-macosclangrelease}"  # or linuxgccrelease

# ============================================================
# INPUT VALIDATION
# ============================================================
INPUT_PDB=$1
OUTPUT_PREFIX=${2:-"relaxed_"}

if [ -z "$INPUT_PDB" ]; then
    echo "Usage: $0 <input.pdb> [output_prefix]"
    echo ""
    echo "Environment variables:"
    echo "  ROSETTA_BIN: Path to Rosetta bin directory"
    echo "  PLATFORM: Binary suffix (default: macosclangrelease)"
    exit 1
fi

if [ ! -f "$INPUT_PDB" ]; then
    echo "ERROR: Input file not found: $INPUT_PDB"
    exit 1
fi

RELAX="${ROSETTA_BIN}/relax.${PLATFORM}"

if [ ! -f "$RELAX" ]; then
    echo "ERROR: Rosetta relax not found: $RELAX"
    echo ""
    echo "Please set ROSETTA_BIN to your Rosetta installation:"
    echo "  export ROSETTA_BIN=/path/to/rosetta/source/bin"
    exit 1
fi

# ============================================================
# RUN ROSETTA RELAX
# ============================================================
echo "========================================================"
echo "ROSETTA RELAX - Structure Preparation"
echo "========================================================"
echo "Input:  $INPUT_PDB"
echo "Output: ${OUTPUT_PREFIX}*"
echo "Binary: $RELAX"
echo "========================================================"
echo ""
echo "Starting relaxation..."
echo "Expected time: 30-60 minutes"
echo ""

$RELAX \
  -s "$INPUT_PDB" \
  -relax:constrain_relax_to_start_coords \
  -relax:coord_constrain_sidechains \
  -relax:ramp_constraints false \
  -ex1 \
  -ex2 \
  -use_input_sc \
  -nstruct 1 \
  -out:prefix "${OUTPUT_PREFIX}" \
  -out:pdb

# ============================================================
# CHECK OUTPUT
# ============================================================
OUTPUT_FILE="${OUTPUT_PREFIX}$(basename ${INPUT_PDB%.pdb})_0001.pdb"

if [ -f "$OUTPUT_FILE" ]; then
    echo ""
    echo "========================================================"
    echo "RELAXATION COMPLETE"
    echo "========================================================"
    echo "Output: $OUTPUT_FILE"
    echo ""
    echo "Next step: Run cartesian_ddG"
    echo "  ./02_cartesian_ddg.sh $OUTPUT_FILE mutfile"
    echo "========================================================"
else
    echo ""
    echo "ERROR: Expected output file not found: $OUTPUT_FILE"
    echo "Check for errors above."
    ls -la ${OUTPUT_PREFIX}*.pdb 2>/dev/null || echo "No output PDB files found."
    exit 1
fi
