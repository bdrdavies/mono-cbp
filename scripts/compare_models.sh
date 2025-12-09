#!/bin/bash
#
# Run Bayesian model comparison on detected events
#
# Usage: ./compare_models.sh [event_dir] [output_file]
#

EVENT_DIR="${1:-results/event_snippets}"
OUTPUT_FILE="${2:-results/classifications.csv}"
CONFIG_FILE="examples/config_example.json"

echo "=========================================="
echo "Model Comparison"
echo "=========================================="
echo "Event directory: $EVENT_DIR"
echo "Output file: $OUTPUT_FILE"
echo "=========================================="
echo

# Check if event directory exists
if [ ! -d "$EVENT_DIR" ]; then
    echo "Error: Event directory not found: $EVENT_DIR"
    exit 1
fi

# Count events
N_EVENTS=$(ls -1 "$EVENT_DIR"/*.npz 2>/dev/null | wc -l)
if [ "$N_EVENTS" -eq 0 ]; then
    echo "Error: No event files found in $EVENT_DIR"
    exit 1
fi

echo "Found $N_EVENTS events to process"
echo

# Create output directory
mkdir -p "$(dirname $OUTPUT_FILE)"

# Build command
CMD="mono-cbp compare-models \
    --event-dir $EVENT_DIR \
    --output $OUTPUT_FILE"

# Add optional config
if [ -f "$CONFIG_FILE" ]; then
    CMD="$CMD --config $CONFIG_FILE"
fi

# Run model comparison
$CMD

if [ $? -eq 0 ]; then
    echo
    echo "=========================================="
    echo "Model comparison complete!"
    echo "Classifications saved to: $OUTPUT_FILE"
    echo "=========================================="

    # Display individual classifications
    if [ -f "$OUTPUT_FILE" ]; then
        echo
        echo "Classification Summary:"
        # Extract filename (remove .npz extension) and classification
        tail -n +2 "$OUTPUT_FILE" | awk -F',' '{
            gsub(/\.npz$/, "", $1);
            printf "   %s: %s\n", $1, $2
        }'

        # Show count summary
        echo
        echo "Classification Counts:"
        tail -n +2 "$OUTPUT_FILE" | cut -d',' -f2 | sort | uniq -c | sort -rn
    fi
else
    echo "Error: Model comparison failed"
    exit 1
fi
