#!/bin/bash
#
# Find transit events in masked light curves
#
# Usage: ./find_transits.sh [data_dir] [output_file]
#

DATA_DIR="${1:-data}"
ECLIPSE_PARAMS="catalogues/TEBC_morph_05_P_7.csv"
SECTOR_TIMES="catalogues/sector_times.csv"
OUTPUT_FILE="${2:-results/transit_events.txt}"
PLOT_DIR="results/vetting_plots"
SNIPPETS_DIR="results/event_snippets"

# Configuration
MAD_THRESHOLD=3.0
METHOD="cb"  # cb = cosine+biweight, cp = cosine+pspline

echo "=========================================="
echo "Transit Finding"
echo "=========================================="
echo "Data directory: $DATA_DIR"
echo "Eclipse params: $ECLIPSE_PARAMS"
echo "Output file: $OUTPUT_FILE"
echo "MAD threshold: $MAD_THRESHOLD"
echo "Method: $METHOD"
echo "=========================================="
echo

# Check files
if [ ! -d "$DATA_DIR" ]; then
    echo "Error: Data directory not found: $DATA_DIR"
    exit 1
fi

if [ ! -f "$ECLIPSE_PARAMS" ]; then
    echo "Error: Eclipse parameters file not found: $ECLIPSE_PARAMS"
    exit 1
fi

# Create output directories
mkdir -p "$(dirname $OUTPUT_FILE)"
mkdir -p "$PLOT_DIR"
mkdir -p "$SNIPPETS_DIR"

# Build command
CMD="mono-cbp find-transits \
    --catalogue $ECLIPSE_PARAMS \
    --data-dir $DATA_DIR \
    --output $OUTPUT_FILE \
    --threshold $MAD_THRESHOLD \
    --method $METHOD \
    --tebc"

# Add optional parameters
if [ -f "$SECTOR_TIMES" ]; then
    CMD="$CMD --sector-times $SECTOR_TIMES"
fi

# Run transit finding
echo "Running: $CMD"
echo
$CMD

if [ $? -eq 0 ]; then
    echo
    echo "=========================================="
    echo "Transit finding complete!"
    echo "Events saved to: $OUTPUT_FILE"
    echo "=========================================="

    # Count events
    if [ -f "$OUTPUT_FILE" ]; then
        N_EVENTS=$(tail -n +2 "$OUTPUT_FILE" | wc -l)
        echo "Total events detected: $N_EVENTS"
    fi
else
    echo "Error: Transit finding failed"
    exit 1
fi
