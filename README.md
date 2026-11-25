[![DOI](https://zenodo.org/badge/1048602949.svg)](https://doi.org/10.5281/zenodo.18475398)
[![MIT License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Python 3.x](https://img.shields.io/badge/python-3.x-blue.svg)](https://www.python.org/downloads/)
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

## Install and configure

For PubChem retrieval, install Python 3.10+ and run
`python -m pip install -r requirements.txt` in this checkout. For the Windows
HSPiP CLI also run `python -m pip install -r requirements-hspip.txt`.
HSPiP must be licensed for CLI use; see the [official CLI documentation](https://www.hansen-solubility.com/HSPiP/CLI.php).
MATLAB is only needed for the legacy MATLAB workflow.

For the Python HSPiP workflow, set these variables in the PowerShell session:

```powershell
$env:HSPIP_PATH = 'C:\Program Files\Hansen-Solubility\HSPiP'
$env:HSPIP_WORKDIR = (Get-Location).Path
```

`HSPIP_PATH` is the directory containing `HSPiP.exe`, or the full executable path.
`HSPIP_WORKDIR` is an existing writable directory used for logs and final output;
it defaults to the current directory. The input path is resolved before switching
to the output directory. HSPiP still writes its intermediate `Out.dat` to its
installation directory, so use a writable installation and one batch at a time.

## Quick Start Guide
### Option 1: Full MATLAB-Python Workflow (Recommended for CAS Lists)
1. Clone this repository: `git clone https://github.com/glsalierno/cas-to-HSPiP_data.git`.
2. Place all files in the same directory.
3. Prepare input: Create a `.mat` file with a structure like this (example in MATLAB):

                        cas_list = {'50-00-0', '67-56-1'};  % Formaldehyde and Methanol CAS numbers
                        save('cas_input.mat', 'cas_list');

4. For this legacy MATLAB path, set the installation directory in `CAS2SMILES2HSPiP.m` and make sure the Python interpreter is on PATH.
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

2. Set `HSPIP_PATH` and optionally `HSPIP_WORKDIR` as shown above. No Python file needs editing.
3. Run: `python HSPiP_CLI_v7.py smiles_input.csv` (or `.txt`).
4. Output: HSP values in a new CSV, with logs for retries/failures.

## Detailed Usage Notes
- **Customizing Scripts**: Edit variables like batch size or retry counts in the Python scripts.
- **Example Output**: For CAS '50-00-0' (Formaldehyde), expect SMILES 'C=O', HSP values like δD=15.5, δP=11.0, δH=7.0 (approximate; depends on HSPiP version).
- **Integration Details**: MATLAB uses `system()` to call Python scripts—ensure Python is in your system PATH.

## Troubleshooting
- **PubChem Errors**: Check API status; add delays if hitting rate limits.
- **HSPiP Fails**: Verify CLI license and paths. Test HSPiP CLI manually first.
- **RDKit issues**: Install the optional requirements in a fresh supported Python environment.
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

## Offline tests

```text
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
```

The tests mock PubChem requests; CI does not call PubChem or exercise MATLAB.
Network calls are paced in each process and retry HTTP 429 or 503 at most twice.
The tool honors Retry-After. If PubChem requests a wait over 30 seconds, the
tool does not retry early or issue a new request during that cooldown. Keep
batch sizes modest and follow [PubChem's usage policy](https://pubchem.ncbi.nlm.nih.gov/docs/programmatic-access).

CI validates the Python configuration errors, but cannot execute the licensed HSPiP program.
