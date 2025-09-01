[![MIT License](https://img.shields.io/badge/License-MIT-green.svg)](https://github.com/glsalierno/pubchem-toxinfo-cas-retriever/blob/main/LICENSE)
[![Python 3.x](https://img.shields.io/badge/python-3.x-blue.svg)](https://www.python.org/downloads/)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2014b%2B-orange.svg)](https://www.mathworks.com/products/matlab.html)
# CAS to SMILES and HSP Retrieval

This repository contains scripts to retrieve SMILES and other chemical properties from PubChem using a CAS number, then compute Hansen Solubility Parameters (HSP) using HSPiP software.

## Files
- `cas2smiles2HSP.m`: MATLAB script to process CAS lists and integrate with Python and HSPiP.
- `get_smiles_InChI_IUPAC_props.py`: Python script to fetch data from PubChem.
- `HSPiP_CLI_v7.py`: Pure Python script for batch processing SMILES with HSPiP CLI (retries, validation via RDKit).

## Requirements
- **Python 3.x**: With `requests` library (`pip install requests`). For `HSPiP_CLI_v7.py`: add `pandas`, `numpy`, `rdkit` (via conda), `pyperclip`, `argparse`.
- **MATLAB**: Base installation; no additional toolboxes needed.
- **HSPiP Software**: Installed with CLI license enabled (available from official HSPiP providers; CLI mode requires specific licensing). Update paths in scripts to your local HSPiP folder.
- A .mat file with CAS numbers (example structure provided in code comments).

## Usage
1. Place all files in the same directory.
2. Update paths in scripts as needed (e.g., HSPiP installation).
3. For MATLAB workflow: Run `cas2smiles2HSP.m`.
4. For pure Python HSPiP processing: `python HSPiP_CLI_v7.py <smiles_file.csv or .txt>` (after obtaining SMILES from PubChem).

## Notes
- This is an anonymized version for sharing; no personal data included.
- Ensure compliance with PubChem API terms and HSPiP licensing.
- For questions, open an issue.

Author: glsalierno  
Date: September 2025  
GitHub: [glsalierno](https://github.com/glsalierno)
