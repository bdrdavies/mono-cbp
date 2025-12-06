#!/bin/bash
#
# Run injection-retrieval testing
#
# Usage: ./inject_retrieve.sh [data_dir] [output_file]
#

DATA_DIR="${1:-data}"
MODELS_FILE="catalogues/transit_models.npz"
CATALOGUE="catalogues/TEBC_morph_05_P_7_ADJUSTED.csv"
OUTPUT_FILE="${2:-results/injection_results.csv}"
CONFIG_FILE="mono_cbp/config_example.json"

echo "=========================================="
echo "Injection-Retrieval Testing"
echo "=========================================="

# Check files
if [ ! -d "$DATA_DIR" ]; then
    echo "Error: Data directory not found: $DATA_DIR"
    exit 1
fi

if [ ! -f "$MODELS_FILE" ]; then
    echo "Error: Transit models file not found: $MODELS_FILE"
    echo "Create models using: python -c 'from mono_cbp.utils import create_transit_models'"
    exit 1
fi

if [ ! -f "$CATALOGUE" ]; then
    echo "Error: Eclipse parameters file not found: $CATALOGUE"
    exit 1
fi

# Create output directory
mkdir -p "$(dirname $OUTPUT_FILE)"

# Build command
CMD="mono-cbp inject-retrieve \
    --models $MODELS_FILE \
    --data-dir $DATA_DIR \
    --catalogue $CATALOGUE \
    --output $OUTPUT_FILE \
    --tebc"

# Add optional config
if [ -f "$CONFIG_FILE" ]; then
    CMD="$CMD --config $CONFIG_FILE"
fi

# Run injection-retrieval
$CMD

if [ $? -eq 0 ]; then
    echo
    echo "=========================================="
    echo "Injection-retrieval complete!"
    echo "Results saved to: $OUTPUT_FILE"
    echo "=========================================="
else
    echo "Error: Injection-retrieval failed"
    exit 1
fi
