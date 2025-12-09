#!/bin/bash
#
# Mask eclipses in eclipsing binary light curves
# Note: Modifies files in-place by adding eclipse_mask column
#
# Usage: ./mask_eclipses.sh [catalogue] [data_dir]
#

CATALOGUE="${1:-catalogues/TEBC_morph_05_P_7.csv}"
DATA_DIR="${2:-data}"
CONFIG_FILE="examples/config_example.json"

echo "=========================================="
echo "Eclipse Masking"
echo "=========================================="
echo "Catalogue: $CATALOGUE"
echo "Data directory: $DATA_DIR"
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

# Build command
CMD="mono-cbp mask-eclipses \
    --catalogue $CATALOGUE \
    --data-dir $DATA_DIR \
    --tebc"

# Add optional config
if [ -f "$CONFIG_FILE" ]; then
    CMD="$CMD --config $CONFIG_FILE"
fi

# Run masking (modifies files in-place)
$CMD

if [ $? -eq 0 ]; then
    echo
    echo "=========================================="
    echo "Eclipse masking complete!"
    echo "Files updated in-place in: $DATA_DIR"
    echo "=========================================="
else
    echo "Error: Eclipse masking failed"
    exit 1
fi
