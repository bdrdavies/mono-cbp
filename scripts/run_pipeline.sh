#!/bin/bash
#
# Run the complete mono-cbp pipeline
#
# Usage: ./run_pipeline.sh [catalogue] [data_dir]
#

# Default paths
CATALOGUE="${1:-catalogues/TEBC_morph_05_P_7_ADJUSTED.csv}"
DATA_DIR="${2:-data}"
SECTOR_TIMES="catalogues/sector_times.csv"
OUTPUT_DIR="results"
CONFIG_FILE="mono_cbp/config_example.json"
PLOT_DIR="${OUTPUT_DIR}/plots"

echo "=========================================="
echo "mono-cbp Pipeline"
echo "=========================================="
echo "Catalogue: $CATALOGUE"
echo "Data directory: $DATA_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "=========================================="
echo

# Check if files exist
if [ ! -f "$CATALOGUE" ]; then
    echo "Error: Catalogue file not found: $CATALOGUE"
    exit 1
fi

if [ ! -d "$DATA_DIR" ]; then
    echo "Error: Data directory not found: $DATA_DIR"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Build command
CMD="mono-cbp run --catalogue $CATALOGUE --data-dir $DATA_DIR --output-dir $OUTPUT_DIR --tebc --plot-dir $PLOT_DIR"

# Add optional parameters
if [ -f "$SECTOR_TIMES" ]; then
    CMD="$CMD --sector-times $SECTOR_TIMES"
fi

if [ -f "$CONFIG_FILE" ]; then
    CMD="$CMD --config $CONFIG_FILE"
fi

# Run pipeline
echo "Running: $CMD"
echo
$CMD

# Check exit status
if [ $? -eq 0 ]; then
    echo
    echo "=========================================="
    echo "Pipeline completed successfully!"
    echo "Results saved to: $OUTPUT_DIR"
    echo "=========================================="
else
    echo
    echo "=========================================="
    echo "Pipeline failed!"
    echo "=========================================="
    exit 1
fi
