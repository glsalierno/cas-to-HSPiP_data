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
# - Set HSPIP_PATH to the HSPiP directory or executable. HSPIP_WORKDIR is optional.
# - Handles invalid SMILES, retries, and preserves original order with NaN for failures.
# - Logs to 'hspip_processing.log' and console.
# - Batch size, retries, and delays configurable for system load.
#
# Author: glsalierno
# Date: September 2025

import os
import subprocess
from pathlib import Path
import logging
from datetime import datetime
import time
import argparse

# Set up command-line argument parsing
parser = argparse.ArgumentParser(description="Process SMILES files with HSPiP CLI.")
parser.add_argument(
    "smiles_file",
    type=str,
    help="Path to the SMILES input file (CSV or TXT)"
)
args = parser.parse_args()

# Resolve input first so a separate output directory does not change its meaning.
smiles_file = Path(args.smiles_file).expanduser().resolve()
file_extension = smiles_file.suffix.lower()
if file_extension not in {'.csv', '.txt'}:
    parser.error(f"unsupported input extension {file_extension!r}; use .csv or .txt")
if not smiles_file.is_file():
    parser.error(f"SMILES input file does not exist: {smiles_file}")

setting = os.environ.get("HSPIP_PATH")
if not setting:
    parser.error("set HSPIP_PATH to the directory containing HSPiP.exe or to the executable itself")
hspip_exe = Path(setting).expanduser().resolve()
if hspip_exe.is_dir():
    hspip_exe = hspip_exe / "HSPiP.exe"
if not hspip_exe.is_file() or hspip_exe.name.lower() != "hspip.exe":
    parser.error(f"HSPIP_PATH does not point to HSPiP.exe or its installation directory: {hspip_exe}")
if not os.access(hspip_exe, os.X_OK):
    parser.error(f"HSPiP.exe is not executable: {hspip_exe}")
workdir = Path(os.environ.get("HSPIP_WORKDIR", ".")).expanduser().resolve()
if not workdir.is_dir() or not os.access(workdir, os.W_OK | os.X_OK):
    parser.error(f"HSPIP_WORKDIR is not a writable directory: {workdir}")
os.chdir(workdir)

import numpy as np
import pandas as pd
import pyperclip
from rdkit import Chem

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
        name = next((col for col in all_smiles_df.columns if col.casefold() in {"smiles", "all_smiles"}), None)
        if name is not None:
            all_smiles = all_smiles_df[name].tolist()
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

                    # HSPiP's CLI writes Out.dat in its installation directory.
                    output_path = hspip_exe.parent / "Out.dat"
                    previous_output = output_path.stat().st_mtime_ns if output_path.is_file() else None
                    subprocess.run(
                        [str(hspip_exe), "Y-MBSX", canon_smi],
                        cwd=hspip_exe.parent,
                        check=True,
                    )

                    for output_retry in range(max_output_retries):
                        if output_path.is_file() and output_path.stat().st_size > 0 and output_path.stat().st_mtime_ns != previous_output:
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
                        break
                except Exception as e:
                    logger.error(f"Error processing SMILES at original index {original_index} ({round_desc} round, retry {retry + 1}/{max_retries}): {smi}: {e}")
                    time.sleep(retry_delay)

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

    return all_hsp_full, failed_output_smiles

# Initial processing and one second round for failed outputs.
batch_size = 100
max_retries = 3
max_output_retries = 3
retry_delay = 2
clipboard_delay = 2
all_hsp_full, failed_output_smiles = process_smiles_batch(
    valid_smiles, batch_size, max_retries, max_output_retries, retry_delay, clipboard_delay, "initial"
)
if failed_output_smiles:
    logger.info("Starting second round for %d SMILES", len(failed_output_smiles))
    retried_results, _ = process_smiles_batch(
        failed_output_smiles, batch_size, max_retries, max_output_retries,
        retry_delay, clipboard_delay, "retry"
    )
    all_hsp_full.extend(retried_results)
else:
    logger.info("No SMILES required a second round of processing")

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
