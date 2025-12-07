## Command-Line Interface

The package provides both a dedicated CLI tool and convenient shell scripts for batch processing.

### Using the CLI Tool

After installation, the `mono-cbp` command is available system-wide:

```bash
# Run the complete pipeline
mono-cbp run \
  --catalogue catalogues/TEBC_morph_05_P_7.csv \
  --data-dir ./test_data \
  --output-dir ./results \
  --sector-times catalogues/sector_times.csv

# Mask eclipses only (modifies files in-place)
mono-cbp mask-eclipses \
  --catalogue catalogues/TEBC_morph_05_P_7.csv \
  --data-dir ./test_data

# Find transits only
mono-cbp find-transits \
  --catalogue catalogues/TEBC_morph_05_P_7.csv \
  --data-dir ./test_data \
  --output transit_events.txt \
  --plot-dir ./plots \
  --threshold 4.0 \
  --method cb

# Compare models for vetting
mono-cbp compare-models \
  --event-dir ./event_snippets \
  --output classifications.csv \
  --output-dir ./results

# Run injection-retrieval analysis
mono-cbp inject-retrieve \
  --models catalogues/transit_models.npz \
  --catalogue catalogues/TEBC_morph_05_P_7.csv \
  --data-dir ./test_data \
  --output inj-ret_results.csv \
  --n-injections 100
```

**Available commands:**
- `run` - Execute the complete pipeline (masking → finding → vetting)
- `mask-eclipses` - Mask primary and secondary eclipses only
- `find-transits` - Detect transit events only
- `compare-models` - Perform Bayesian model comparison on events
- `inject-retrieve` - Run injection-retrieval completeness analysis

**Common options:**
- `--config PATH` - Path to JSON configuration file (overrides defaults)
- `--tebc` - Use TEBC catalogue format with `*_2g` and `*_pf` columns
- `--plot-dir PATH` - Directory for saving diagnostic plots
- `--output-dir PATH` - Directory for saving results

Run any command with `--help` for detailed options:
```bash
mono-cbp find-transits --help
```

### Using Shell Scripts

Shell scripts in the `scripts/` directory provide alternative command-line interfaces with preset defaults:

```bash
# Run the complete pipeline
./scripts/run_pipeline.sh [catalogue] [data_dir]

# Eclipse masking only
./scripts/mask_eclipses.sh

# Transit finding only
./scripts/find_transits.sh

# Model comparison on detected events
./scripts/compare_models.sh

# Injection-retrieval testing
./scripts/inject_retrieve.sh
```

These scripts use default paths configured in the script files and can be edited for your specific project structure.