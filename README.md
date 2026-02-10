[![DOI](https://zenodo.org/badge/1048602949.svg)](https://doi.org/10.5281/zenodo.18475398)
[![MIT License](https://img.shields.io/badge/License-MIT-green.svg)](https://github.com/glsalierno/pubchem-toxinfo-cas-retriever/blob/main/LICENSE)
[![Python 3.x](https://img.shields.io/badge/python-3.x-blue.svg)](https://www.python.org/downloads/)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2014b%2B-orange.svg)](https://www.mathworks.com/products/matlab.html)
# CAS to SMILES and HSP Retrieval

This repository provides scripts to retrieve SMILES, InChI, IUPAC names, and other chemical properties from PubChem using CAS numbers, then compute Hansen Solubility Parameters (HSP) using the HSPiP software. It supports integration between MATLAB and Python for efficient chemical data processing.

## About the Repository
- **Purpose**: Automate the conversion of CAS numbers to chemical identifiers (via PubChem API) and calculate HSP values (via HSPiP CLI).
- **Key Features**: Batch processing, retries for failed queries, validation with RDKit, and MATLAB-Python bridging.
- **DOI**: [10.5281/zenodo.18475390](https://doi.org/10.5281/zenodo.18475390)
- **Author**: glsalierno (GitHub: [glsalierno](https://github.com/glsalierno))
- **Date**: September 2025
- **License**: MIT (see [LICENSE](LICENSE))

## Files
- `CAS2SMILES2HSPiP.m`: MATLAB script to process CAS lists, call Python for PubChem data, and integrate with HSPiP.
- `get_smiles_InChI_IUPAC_props.py`: Python script to fetch SMILES, InChI, IUPAC names, and properties from PubChem.
- `HSPiP_CLI_v7.py`: Pure Python script for batch-processing SMILES with HSPiP CLI, including retries and RDKit validation.
- `README.md`: This file.
- `LICENSE`: MIT license.

## Requirements
### Must-Haves
- **Python 3.x** … (unchanged)
- **MATLAB** … (unchanged)
- **HSPiP Software**: Must be installed with CLI mode enabled.  
  Obtain from the official HSPiP website: [HSPiP | Hansen Solubility Parameters](https://www.hansen-solubility.com/HSPiP).  
  CLI requires a specific license—contact HSPiP support if unsure.  
  For details on the Command Line Interface, see [HSPiP CLI Documentation](https://www.hansen-solubility.com/HSPiP/CLI.php).

### Optional
- A MATLAB `.mat` file with CAS numbers (see example below).

**Note**: No internet access is needed for HSPiP computation, but PubChem queries require it. Ensure compliance with [PubChem API terms](https://pubchem.ncbi.nlm.nih.gov/docs/programmatic-access).

---

> **⚠️ IMPORTANT: Update Paths Before Running**
>
> You **must** edit `HSPiP_CLI_v7.py` (and optionally `CAS2SMILES2HSPiP.m`):
>
> - **Line 151** → `hspip_path = 'PATH_TO_HSPIP_INSTALLATION'`  
>   Example: `'C:\\Program Files\\Hansen-Solubility\\HSPiP'`
>
> - **Line 251** → `fold = 'PATH_TO_WORKING_DIRECTORY'` (optional)  
>   Only needed if you want the script to return to a specific folder after processing.  
>   You can comment out or delete these two lines if everything runs in the same folder.
>
> Replace the placeholder strings with your actual paths (use double backslashes `\\` on Windows).

**Note**: No internet access is needed for HSPiP computation, but PubChem queries require it. Ensure compliance with [PubChem API terms](https://pubchem.ncbi.nlm.nih.gov/docs/programmatic-access) (e.g., rate limits: ~5 requests/second).

## Quick Start Guide
### Option 1: Full MATLAB-Python Workflow (Recommended for CAS Lists)
1. Clone this repository: `git clone https://github.com/glsalierno/cas-to-HSPiP_data.git`.
2. Place all files in the same directory.
3. Prepare input: Create a `.mat` file with a structure like this (example in MATLAB):

                        cas_list = {'50-00-0', '67-56-1'};  % Formaldehyde and Methanol CAS numbers
                        save('cas_input.mat', 'cas_list');

4. Update paths in `CAS2SMILES2HSPiP.m` and other scripts (e.g., HSPiP install dir, Python executable). **as shown in the ⚠️ IMPORTANT box above.**
5. Run in MATLAB: Open and execute `CAS2SMILES2HSPiP.m`. It will:
- Load CAS from `.mat`.
- Call `get_smiles_InChI_IUPAC_props.py` for PubChem data.
- Call `HSPiP_CLI_v7.py` for HSP computation.
6. Output: Results saved as CSV/Excel files with SMILES, properties, and HSP values.

### Option 2: Python-Only for SMILES to HSP
1. Obtain SMILES first (e.g., manually or via `get_smiles_InChI_IUPAC_props.py`).
- Example input file `smiles_input.csv`:

                      SMILES
                      C=O
                      CO

2. Update paths in `HSPiP_CLI_v7.py`.
3. Run: `python HSPiP_CLI_v7.py smiles_input.csv` (or `.txt`).
4. Output: HSP values in a new CSV, with logs for retries/failures.

## Detailed Usage Notes
- **Customizing Scripts**: Edit variables like batch size or retry counts in the Python scripts.
- **Example Output**: For CAS '50-00-0' (Formaldehyde), expect SMILES 'C=O', HSP values like δD=15.5, δP=11.0, δH=7.0 (approximate; depends on HSPiP version).
- **Integration Details**: MATLAB uses `system()` to call Python scripts—ensure Python is in your system PATH.

## Troubleshooting
- **PubChem Errors**: Check API status; add delays if hitting rate limits.
- **HSPiP Fails**: Verify CLI license and paths. Test HSPiP CLI manually first.
- **RDKit Issues**: Ensure conda installation; common for dependency conflicts.
- **MATLAB-Python Bridge Fails**: Confirm Python executable path in MATLAB (use `pyenv`).
- If stuck, open an issue with error logs.

## Notes
- This is an anonymized version; no sensitive data included.
- For large batches, monitor API usage to avoid blocks.
- Contributions welcome—fork and PR!

If you have questions or need help adapting this, feel free to ask!

## References
- Official HSPiP Website: [HSPiP | Hansen Solubility Parameters](https://www.hansen-solubility.com/HSPiP)
- HSPiP CLI Documentation: [Command Line Interface (CLI)](https://www.hansen-solubility.com/HSPiP/CLI.php)
- HSPiP CLI Guide: [HSPiP Command Line Interface.docx](https://www.hansen-solubility.com/contents/HSPiP%20Command%20Line%20Interface.docx)
- PubChem API: [Programmatic Access Documentation](https://pubchem.ncbi.nlm.nih.gov/docs/programmatic-access)
- RDKit: [Official Documentation](https://www.rdkit.org/docs/)
