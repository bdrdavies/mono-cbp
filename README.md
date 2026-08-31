# mono-cbp: Search for Monotransits of Circumbinary Planets

A Python package for detecting circumbinary planets in TESS eclipsing binary light curves through the identification of single transit events ("monotransits").

[![Python Version](https://img.shields.io/badge/python-3.9%2B-blue)](https://www.python.org/downloads/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## Overview

`mono-cbp` is a pipeline designed to systematically search for circumbinary planets by detecting individual transits in eclipsing binary systems. The pipeline automates masking stellar eclipses, threshold crossing event (TCE) detection, Bayesian vetting, and completeness analysis, and is capable of processing large catalogues of eclipsing binaries quite quickly.

Circumbinary planets (CBPs) are planets that orbit both components of a stellar binary. Currently, only 14 transiting CBPs have been discovered, while thousands of transiting planets have been discovered orbiting around single stars. This discrepancy is partially due to observational bias, but also due to their physical properties; CBPs tend to have long periods, due to orbital stability constraints, and their orbital planes precess such that a CBP can oscillate between transiting and non-transiting configurations on timescales of decades.

For transit surveys that do not continuously observe a given field for multiple-year durations (such as *Kepler*), the probability of observing multiple transits of long-period planets (such as CBPs) drops drastically. Therefore, such planets manifest themselves as single-transit events, which can then be followed-up by other observatories for confirmation. *TESS* observes almost the entire sky, providing a wealth of photometric data for thousands of potential CBP hosts. `mono-cbp` allows one to find, in a systematic and efficient manner, candidate CBPs from *TESS* observations, and can be applied to other photometric datasets. 

## Key Features

- **Eclipse Masking**:  Mask primary and secondary eclipses in eclipsing binary light curves using eclipse positions and widths and binary ephemeris provided by an input catalogue
- **Transit Detection**: Removes unwanted trends from the input light curves and performs single-event detection using the by identifying TCEs (see [Hawthorn et al. 2024](https://academic.oup.com/mnras/article/528/2/1841/7589620?login=false))
- **Bayesian Model Comparison**: Event classification to discern transit-like events and systematics/detrending artefacts
- **Injection-Retrieval Testing**: Completeness analysis through synthetic transit injection and recovery statistics
- **Modular Architecture**: Use individual components independently or run the complete integrated pipeline
- **Configuration-Driven**: Customise parameters via Python dictionaries or JSON files without modifying code
- **Command-Line Interface**: Shell scripts and CLI subcommands for batch processing and reproducibility

## Installation

### Requirements

- Python 3.9 or higher (tested against 3.9, 3.10, and 3.13)

### From PyPI (Recommended)

The easiest way to install `mono-cbp` is from PyPI:

```bash
pip install mono-cbp
```

It is advisable to install `mono-cbp` into a Python environment using your favourite package manager, e.g. for `conda`:

```bash
conda create --name mono-cbp python=3.9
conda activate mono-cbp
pip install mono-cbp
```

This installs the package and creates the `mono-cbp` command-line tool.

### From Source

For development or to use the latest unreleased features:

```bash
git clone https://github.com/bdrdavies/mono-cbp.git
cd mono-cbp
pip install -e .
```

### Verify Installation

To check that the installation has been successful:

```bash
python -c "import mono_cbp; print(mono_cbp.__version__)"
mono-cbp --help
```

### Dependencies

All dependencies are automatically installed when you install `mono-cbp`.

See [pyproject.toml](pyproject.toml) for the complete dependency list and version constraints.

### Troubleshooting Installation

If you encounter issues:

- **Python version**: The package has been tested most thoroughly with Python 3.9.
- **Dependency conflicts**: If you have conflicts with existing packages, create a fresh environment
- **Import errors**: If you see import errors, ensure all dependencies installed correctly by checking `pip show mono-cbp` and comparing against [pyproject.toml](pyproject.toml)

## Examples & Tutorials

There are a series of Jupyter notebooks in the `examples/` directory to demonstrate how to use the package in your own code:

1. **[00_download_light_curves.ipynb](examples/00_download_light_curves.ipynb)** - Download TESS light curves in the `mono-cbp` format using [lightkurve](https://lightkurve.github.io/lightkurve/)
2. **[01_complete_pipeline.ipynb](examples/01_complete_pipeline.ipynb)** - End-to-end processing on sample data
3. **[02_eclipse_masking.ipynb](examples/02_eclipse_masking.ipynb)** - Eclipse masking demo
4. **[03_transit_finding.ipynb](examples/03_transit_finding.ipynb)** - TCE detection example
5. **[04_model_comparison.ipynb](examples/04_model_comparison.ipynb)** - Bayesian model comparison example
6. **[05_injection_retrieval.ipynb](examples/05_injection_retrieval.ipynb)** - Completeness testing

## Documentation

Documentation is available in the `docs/` directory:

- **[docs/quickstart.md](docs/quickstart.md)** - Quickstart guide
- **[docs/data_formats.md](docs/data_formats.md)** - Input and output data format specifications
- **[docs/configuration.md](docs/configuration.md)** - Configuration system reference
- **[docs/api_reference.md](docs/api_reference.md)** - API documentation

## Support & Contact

For questions, issues, or feature requests:
- **Issues:** Open an issue on [GitHub Issues](https://github.com/bdrdavies/mono-cbp/issues)
- **Documentation:** Review the [full documentation](docs/)
- **Email:** ben.d.r.davies@warwick.ac.uk

## Credit

If you make use of this code, please cite the following paper:
[B. D. R. Davies et al. "Finding Circumbinary Planets: A Semi-Automated Transit Search of TESS Eclipsing Binaries", MNRAS 548, 4 (June 2026)](https://ui.adsabs.harvard.edu/abs/2026MNRAS.548ag743D/abstract)

```
@ARTICLE{Davies2026,
       author = {{Davies}, Benjamin D.~R. and {Brown}, David J.~A. and {Gill}, Samuel and {French}, Jenni R.},
        title = "{Finding circumbinary planets: a semi-automated transit search of TESS eclipsing binaries}",
      journal = {\mnras},
     keywords = {exoplanets, planets and satellites: detection, binaries: eclipsing, planetary systems, software: data analysis, Earth and Planetary Astrophysics, Instrumentation and Methods for Astrophysics, Solar and Stellar Astrophysics},
         year = 2026,
        month = jun,
       volume = {548},
       number = {4},
          eid = {stag743},
        pages = {stag743},
          doi = {10.1093/mnras/stag743},
archivePrefix = {arXiv},
       eprint = {2604.09435},
 primaryClass = {astro-ph.EP},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2026MNRAS.548ag743D},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```