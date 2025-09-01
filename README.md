# cas-to-HSPiP_data
Scripts to retrieve SMILES, InChI, IUPAC names, and properties from PubChem via CAS numbers, then compute Hansen Solubility Parameters (HSP) using HSPiP. MATLAB-Python integration for chem data processing.

# CAS to SMILES and HSP Retrieval

This repository contains scripts to retrieve SMILES and other chemical properties from PubChem using a CAS number, then compute Hansen Solubility Parameters (HSP) using HSPiP software.

## Files
- `cas2smiles2HSP.m`: MATLAB script to process CAS lists and integrate with Python and HSPiP.
- `get_smiles_InChI_IUPAC_props.py`: Python script to fetch data from PubChem.

## Requirements
- **Python 3.x**: With `requests` library (`pip install requests`).
- **MATLAB**: Base installation; no additional toolboxes needed.
- **HSPiP Software**: Installed with CLI license enabled (available from official HSPiP providers; CLI mode requires specific licensing https://www.hansen-solubility.com/HSPiP/CLI.php).
- A .mat file with CAS numbers (example structure provided in code comments).

## Usage
1. Place all files in the same directory.
2. Update paths in MATLAB script as needed.
3. Run the MATLAB script.

## Notes
- Ensure compliance with PubChem API terms and HSPiP licensing.
- For questions, open an issue.

Author: glsalierno  
Date: September 2025  
GitHub: [glsalierno](https://github.com/glsalierno)
