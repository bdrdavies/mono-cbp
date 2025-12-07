# Installation

## Requirements

- Python 3.8 or higher
- pip or conda package manager

## Install From Source

Clone the repository and install in development mode:

```bash
git clone https://github.com/bdrdavies/mono-cbp.git
cd mono-cbp
pip install -e .
```

This allows you to make changes to the source code and have them immediately reflected.

## Verification

Verify your installation:

```bash
python -c "import mono_cbp; print(mono_cbp.__version__)"
```

## Updating

To update to the latest version:

```bash
cd mono-cbp
git pull
pip install -e . --upgrade
```