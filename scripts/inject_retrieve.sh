#!/bin/bash
#
# Run injection-retrieval testing
#
# Usage: ./inject_retrieve.sh [data_dir] [output_file]
#

DATA_DIR="${1:-data/lightcurves/masked}"
MODELS_FILE="data/transit_models.npz"
ECLIPSE_PARAMS="data/catalogues/eclipse_params.csv"
OUTPUT_FILE="${2:-results/injection_results.csv}"
N_INJECTIONS=100

echo "=========================================="
echo "Injection-Retrieval Testing"
echo "=========================================="
echo "Data directory: $DATA_DIR"
echo "Transit models: $MODELS_FILE"
echo "Output file: $OUTPUT_FILE"
echo "Number of injections: $N_INJECTIONS"
echo "=========================================="
echo

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

if [ ! -f "$ECLIPSE_PARAMS" ]; then
    echo "Error: Eclipse parameters file not found: $ECLIPSE_PARAMS"
    exit 1
fi

# Create output directory
mkdir -p "$(dirname $OUTPUT_FILE)"

# Run injection-retrieval
mono-cbp inject-retrieve \
    --models "$MODELS_FILE" \
    --data-dir "$DATA_DIR" \
    --eclipse-params "$ECLIPSE_PARAMS" \
    --output "$OUTPUT_FILE" \
    --n-injections "$N_INJECTIONS"

if [ $? -eq 0 ]; then
    echo
    echo "=========================================="
    echo "Injection-retrieval complete!"
    echo "Results saved to: $OUTPUT_FILE"
    echo "=========================================="

    # Calculate recovery rate
    if [ -f "$OUTPUT_FILE" ]; then
        N_TOTAL=$(tail -n +2 "$OUTPUT_FILE" | wc -l)
        N_RECOVERED=$(tail -n +2 "$OUTPUT_FILE" | cut -d',' -f7 | grep -c "True")
        RECOVERY_RATE=$(echo "scale=2; 100 * $N_RECOVERED / $N_TOTAL" | bc)
        echo
        echo "Recovery rate: $RECOVERY_RATE% ($N_RECOVERED/$N_TOTAL)"
    fi
else
    echo "Error: Injection-retrieval failed"
    exit 1
fi
