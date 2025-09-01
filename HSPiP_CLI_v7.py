# HSPiP_CLI_v7.py
#
# This script processes a list of SMILES from an input file (CSV or TXT) using HSPiP CLI.
# It handles batch processing, retries for failures, canonicalization with RDKit,
# and saves results to CSV/Excel. Designed for integration with PubChem-retrieved SMILES.
#
# Requirements:
# - Python 3.x
# - Libraries: pandas, numpy, rdkit (install via: conda install -c conda-forge rdkit), pyperclip, argparse
# - HSPiP software with CLI license enabled (official sources; CLI requires licensing)
# - Input: CSV/TXT file with SMILES (column 'all_smiles' for CSV or plain lines for TXT)
#
# Usage:
# python HSPiP_CLI_v7.py <smiles_file.csv or .txt>
#
# Notes:
# - Replace 'PATH_TO_HSPIP_INSTALLATION' with your local HSPiP directory (e.g., 'C:\Path\To\HSPiP').
# - Anonymized paths for sharing; update as needed.
# - Handles invalid SMILES, retries, and preserves original order with NaN for failures.
# - Logs to 'hspip_processing.log' and console.
# - Batch size, retries, and delays configurable for system load.
#
# Author: glsalierno
# Date: September 2025

import os
import subprocess
import pandas as pd
import numpy as np
from pathlib import Path
from rdkit import Chem
import logging
from datetime import datetime
import time
import shutil
import pyperclip
import argparse
import sys

# Set up command-line argument parsing
parser = argparse.ArgumentParser(description="Process SMILES files with HSPiP CLI.")
parser.add_argument(
    "smiles_file",
    type=str,
    help="Path to the SMILES input file (CSV or TXT)"
)
args = parser.parse_args()

# Validate the provided file
smiles_file = args.smiles_file
file_extension = Path(smiles_file).suffix.lower()
if file_extension not in ['.csv', '.txt']:
    print(f"Error: Unsupported file extension {file_extension}. Use .csv or .txt")
    sys.exit(1)
if not os.path.exists(smiles_file):
    print(f"Error: SMILES file {smiles_file} not found")
    sys.exit(1)

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('hspip_processing.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Function to canonicalize SMILES
def canonicalize_smiles(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return Chem.MolToSmiles(mol, canonical=True)
    except Exception:
        return None

# Create a batch file to execute HSPiP commands one at a time
def create_batch_file(smi, batch_file_path, hspip_path, hspip_cmd):
    with open(batch_file_path, 'w') as f:
        f.write(f'cd "{hspip_path}"\n')
        f.write(f'{hspip_cmd} "{smi}"\n')
        f.write('exit\n')

# Save intermediate results
def save_intermediate_results(all_hsp_full, original_smiles, timestamp, batch_count):
    ordered_results = [None] * len(original_smiles)
    for df, original_index in all_hsp_full:
        ordered_results[original_index] = df

    final_dfs = []
    smiles_column = []
    for i, smi in enumerate(original_smiles):
        smiles_column.append(smi)
        if ordered_results[i] is not None:
            final_dfs.append(ordered_results[i])
        else:
            empty_df = pd.DataFrame({
                col: [np.nan] for col in [
                    "HSPiP_SMILES", "Formula", "D", "P", "H", "HDon", "HAcc", "MWt", "Density", "MVol",
                    "Area", "Ovality", "BPt", "MPt", "Tc", "Pc", "Vc", "Zc", "AntA", "AntB",
                    "AntC", "Ant1T", "LogKow", "LogS", "Henry", "LogOHR", "RI", "Hfus", "HvBPt",
                    "Trouton", "RER", "Abra", "Abrb", "EdmiW", "Parachor", "RD", "Cp", "log",
                    "Cond", "SurfTen", "HeavyAtom", "C", "H1", "Br", "Cl", "F", "I", "N", "O",
                    "P1", "S", "Si", "B", "MaxPc", "MinMc", "Sym", "MCI", "Hcomb", "Hform",
                    "Gform", "FGList"
                ]
            })
            final_dfs.append(empty_df)

    amHSP = pd.concat([pd.DataFrame({'SMILES': smiles_column}), pd.concat(final_dfs, ignore_index=True)], axis=1)

    csv_path = f'intermediate_HSP_{timestamp}.csv'
    excel_path = f'intermediate_HSP_{timestamp}.xlsx'
    amHSP.to_csv(csv_path, index=False)
    logger.info(f"Saved intermediate results to {csv_path} after {batch_count} batches")
    amHSP.to_excel(excel_path, index=False)
    logger.info(f"Saved intermediate results to {excel_path} after {batch_count} batches")

# Load all_smiles from file
try:
    if file_extension == '.csv':
        all_smiles_df = pd.read_csv(smiles_file)
        if 'all_smiles' in all_smiles_df.columns:
            all_smiles = all_smiles_df['all_smiles'].tolist()
        else:
            all_smiles_df = pd.read_csv(smiles_file, header=None)
            all_smiles = all_smiles_df[0].tolist()
    elif file_extension == '.txt':
        with open(smiles_file, 'r') as f:
            all_smiles = [line.strip() for line in f if line.strip()]
    logger.info(f"Loaded {len(all_smiles)} SMILES from {smiles_file}")
except FileNotFoundError:
    logger.error(f"SMILES file {smiles_file} not found")
    raise
except Exception as e:
    logger.error(f"Error loading SMILES file {smiles_file}: {e}")
    raise

# Filter out "CAS number not found" entries but keep track of indices
original_smiles = all_smiles.copy()  # Preserve original for final output
valid_smiles_indices = [i for i, smi in enumerate(all_smiles) if smi != "CAS number not found"]
invalid_cas_indices = [i for i, smi in enumerate(all_smiles) if smi == "CAS number not found"]
all_smiles = [smi for smi in all_smiles if smi != "CAS number not found"]
logger.info(f"Filtered out {len(invalid_cas_indices)} 'CAS number not found' entries at indices: {invalid_cas_indices}")
logger.info(f"Remaining {len(all_smiles)} SMILES to process")

# Validate SMILES and track invalid ones
invalid_smiles_indices = []
valid_smiles = []
for i, smi in enumerate(all_smiles):
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            logger.warning(f"Invalid SMILES at index {valid_smiles_indices[i]} (original index): {smi}")
            invalid_smiles_indices.append(valid_smiles_indices[i])
        else:
            logger.info(f"SMILES at index {valid_smiles_indices[i]} (original index): {smi} is valid")
            valid_smiles.append((smi, valid_smiles_indices[i]))
    except Exception as e:
        logger.warning(f"Error validating SMILES at index {valid_smiles_indices[i]} (original index): {smi}: {e}")
        invalid_smiles_indices.append(valid_smiles_indices[i])

logger.info(f"Found {len(invalid_smiles_indices)} invalid SMILES at indices: {sorted(invalid_smiles_indices)}")

# HSPiP configuration - Replace with your local HSPiP installation path
hspip_path = 'PATH_TO_HSPIP_INSTALLATION'  # e.g., 'C:\\Program Files\\Hansen-Solubility\\HSPiP'
hspip_cmd = 'HSPiP.exe Y-MBSX'

# Timestamp for output files
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

# Function to process a batch of SMILES (initial or retry round)
def process_smiles_batch(smiles_list, batch_size, max_retries, max_output_retries, retry_delay, clipboard_delay, round_desc):
    all_hsp_full = []
    failed_output_smiles = []
    batch_count = 0

    for batch_start in range(0, len(smiles_list), batch_size):
        batch_end = min(batch_start + batch_size, len(smiles_list))
        batch = smiles_list[batch_start:batch_end]
        batch_count += 1
        logger.info(f"Processing {round_desc} batch {batch_count} ({batch_start + 1} to {batch_end})")

        for smi, original_index in batch:
            start_time = pd.Timestamp.now()
            success = False
            for retry in range(max_retries):
                try:
                    # Canonicalize SMILES
                    canon_smi = canonicalize_smiles(smi)
                    if canon_smi is None:
                        logger.warning(f"Failed to canonicalize SMILES at original index {original_index}: {smi}")
                        invalid_smiles_indices.append(original_index)
                        break

                    # Create temporary batch file
                    batch_file_path = os.path.join(os.getcwd(), f'temp_hspip_{original_index}.bat')

                    create_batch_file(canon_smi, batch_file_path, hspip_path, hspip_cmd)

                    # Execute batch file
                    subprocess.run(batch_file_path, shell=True, check=True)

                    # Read output with retries for empty/invalid files
                    output_path = os.path.join(hspip_path, 'Out.dat')
                    for output_retry in range(max_output_retries):
                        if os.path.exists(output_path) and os.path.getsize(output_path) > 0:
                            try:
                                opts = pd.read_csv(output_path, delimiter='\t', skiprows=1, names=[
                                    "HSPiP_SMILES", "Formula", "D", "P", "H", "HDon", "HAcc", "MWt", "Density", "MVol",
                                    "Area", "Ovality", "BPt", "MPt", "Tc", "Pc", "Vc", "Zc", "AntA", "AntB",
                                    "AntC", "Ant1T", "LogKow", "LogS", "Henry", "LogOHR", "RI", "Hfus", "HvBPt",
                                    "Trouton", "RER", "Abra", "Abrb", "EdmiW", "Parachor", "RD", "Cp", "log",
                                    "Cond", "SurfTen", "HeavyAtom", "C", "H1", "Br", "Cl", "F", "I", "N", "O",
                                    "P1", "S", "Si", "B", "MaxPc", "MinMc", "Sym", "MCI", "Hcomb", "Hform",
                                    "Gform", "FGList"
                                ])
                                hsp_df = pd.DataFrame(opts.iloc[0]).T
                                all_hsp_full.append((hsp_df, original_index))
                                success = True
                                logger.info(f"Successfully processed SMILES at original index {original_index} ({round_desc} round): {smi}")
                                break
                            except Exception as e:
                                logger.warning(f"Error reading output (retry {output_retry + 1}/{max_output_retries}) for SMILES at original index {original_index} ({round_desc} round): {smi}: {e}")
                                time.sleep(retry_delay)
                        else:
                            logger.warning(f"Output file missing or empty (retry {output_retry + 1}/{max_output_retries}) for SMILES at original index {original_index} ({round_desc} round): {smi}")
                            time.sleep(retry_delay)

                    if success:
                        break
                    else:
                        logger.warning(f"Max output retries reached for SMILES at original index {original_index} ({round_desc} round): {smi}")
                        if round_desc == "initial":
                            failed_output_smiles.append((smi, original_index))
                        else:
                            invalid_smiles_indices.append(original_index)
                        break
                except Exception as e:
                    logger.error(f"Error processing SMILES at original index {original_index} ({round_desc} round, retry {retry + 1}/{max_retries}): {smi}: {e}")
                    time.sleep(retry_delay)
                finally:
                    # Clean up batch file
                    if os.path.exists(batch_file_path):
                        os.remove(batch_file_path)

            if not success:
                logger.error(f"Failed to process SMILES at original index {original_index} ({round_desc} round) after {max_retries} retries: {smi}")
                if round_desc == "initial":
                    failed_output_smiles.append((smi, original_index))
                else:
                    invalid_smiles_indices.append(original_index)

            # Clear clipboard to avoid contention
            pyperclip.copy('')
            time.sleep(clipboard_delay)

            elapsed = (pd.Timestamp.now() - start_time).total_seconds()
            logger.info(f"Time elapsed for SMILES at original index {original_index} ({round_desc} round): {elapsed:.2f} seconds")

        # Save intermediate results every 10 batches
        if batch_count % 10 == 0:
            logger.info(f"Saving intermediate results after {batch_count} {round_desc} batches")
            save_intermediate_results(all_hsp_full, original_smiles, timestamp, batch_count)

        # Delay between batches to further reduce system load
        if batch_start + batch_size < len(smiles_list):
            logger.info(f"Pausing between {round_desc} batches for 10 seconds...")
            time.sleep(10)

    # Save intermediate results at the end of the round if not already saved
    if batch_count % 10 != 0:
        logger.info(f"Saving intermediate results at end of {round_desc} round after {batch_count} batches")
        save_intermediate_results(all_hsp_full, original_smiles, timestamp, batch_count)

# Initial processing of valid SMILES
batch_size = 100
max_retries = 3
max_output_retries = 3
retry_delay = 2  # seconds
clipboard_delay = 2  # Increased delay to avoid clipboard contention
process_smiles_batch(valid_smiles, batch_size, max_retries, max_output_retries, retry_delay, clipboard_delay, "initial")

# Second round for failed SMILES
if failed_output_smiles:
    logger.info(f"Starting second round processing for {len(failed_output_smiles)} SMILES that failed due to empty/invalid output")
    process_smiles_batch(failed_output_smiles, batch_size, max_retries, max_output_retries, retry_delay, clipboard_delay, "retry")
else:
    logger.info("No SMILES required a second round of processing")

# Return to original directory - Replace with your working directory if needed
fold = 'PATH_TO_WORKING_DIRECTORY'  # e.g., 'C:\\Users\\YourUser\\Projects\\AI\\Retrieve'
try:
    os.chdir(fold)
    logger.info(f"Returned to original directory: {fold}")
except OSError as e:
    logger.error(f"Failed to change directory to {fold}: {e}")
    raise

# Combine results while preserving original order
if all_hsp_full:
    ordered_results = [None] * len(original_smiles)
    for df, original_index in all_hsp_full:
        ordered_results[original_index] = df

    final_dfs = []
    smiles_column = []
    for i, smi in enumerate(original_smiles):
        smiles_column.append(smi)
        if ordered_results[i] is not None:
            final_dfs.append(ordered_results[i])
        else:
            empty_df = pd.DataFrame({
                col: [np.nan] for col in [
                    "HSPiP_SMILES", "Formula", "D", "P", "H", "HDon", "HAcc", "MWt", "Density", "MVol",
                    "Area", "Ovality", "BPt", "MPt", "Tc", "Pc", "Vc", "Zc", "AntA", "AntB",
                    "AntC", "Ant1T", "LogKow", "LogS", "Henry", "LogOHR", "RI", "Hfus", "HvBPt",
                    "Trouton", "RER", "Abra", "Abrb", "EdmiW", "Parachor", "RD", "Cp", "log",
                    "Cond", "SurfTen", "HeavyAtom", "C", "H1", "Br", "Cl", "F", "I", "N", "O",
                    "P1", "S", "Si", "B", "MaxPc", "MinMc", "Sym", "MCI", "Hcomb", "Hform",
                    "Gform", "FGList"
                ]
            })
            final_dfs.append(empty_df)

    amHSP = pd.concat([pd.DataFrame({'SMILES': smiles_column}), pd.concat(final_dfs, ignore_index=True)], axis=1)

    csv_path = f'all_HSP_{timestamp}.csv'
    excel_path = f'all_HSP_{timestamp}.xlsx'
    amHSP.to_csv(csv_path, index=False)
    logger.info(f"Saved final results to {csv_path}")
    amHSP.to_excel(excel_path, index=False)
    logger.info(f"Saved final results to {excel_path}")
else:
    logger.warning("No valid HSPiP results to save")
    amHSP = pd.DataFrame({
        'SMILES': original_smiles,
        **{col: [np.nan] * len(original_smiles) for col in [
            "HSPiP_SMILES", "Formula", "D", "P", "H", "HDon", "HAcc", "MWt", "Density", "MVol",
            "Area", "Ovality", "BPt", "MPt", "Tc", "Pc", "Vc", "Zc", "AntA", "AntB",
            "AntC", "Ant1T", "LogKow", "LogS", "Henry", "LogOHR", "RI", "Hfus", "HvBPt",
            "Trouton", "RER", "Abra", "Abrb", "EdmiW", "Parachor", "RD", "Cp", "log",
            "Cond", "SurfTen", "HeavyAtom", "C", "H1", "Br", "Cl", "F", "I", "N", "O",
            "P1", "S", "Si", "B", "MaxPc", "MinMc", "Sym", "MCI", "Hcomb", "Hform",
            "Gform", "FGList"
        ]}
    })
    csv_path = f'all_HSP_{timestamp}_empty.csv'
    excel_path = f'all_HSP_{timestamp}_empty.xlsx'
    amHSP.to_csv(csv_path, index=False)
    logger.info(f"Saved empty final results to {csv_path}")
    amHSP.to_excel(excel_path, index=False)
    logger.info(f"Saved empty final results to {excel_path}")

# Log final summary
logger.info(f"Processing complete. Processed {len(all_hsp_full)} valid SMILES. Invalid or skipped indices: {sorted(set(invalid_smiles_indices + invalid_cas_indices))}")
logger.info(f"Retried {len(failed_output_smiles)} SMILES in second round")
