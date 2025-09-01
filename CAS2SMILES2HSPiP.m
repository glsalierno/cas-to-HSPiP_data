% cas2smiles2HSPiP.m
% 
% This script processes a list of CAS numbers to retrieve SMILES strings from PubChem
% using a Python helper script, then uses HSPiP (Hansen Solubility Parameters in Practice)
% software to compute Hansen Solubility Parameters (HSP) and other properties for each compound.
% 
% Requirements:
% - MATLAB (tested with recent versions; base installation should suffice)
% - Python 3.x with 'requests' library installed (pip install requests)
% - HSPiP software installed with CLI (Command Line Interface) license enabled
%   (Download from official sources; CLI mode requires specific licensing)
% - A .mat file containing CAS numbers (e.g., 'CAS.mat' with variable 'am' where am(:,1) are CAS strings)
% - The Python script 'get_smiles_InChI_IUPAC_props.py' in the same directory
% 
% Usage:
% 1. Ensure HSPiP.exe is in your system's PATH or update the 'hsip_path' variable below.
% 2. Run this script in the directory containing your .mat file and the Python script.
% 
% Notes:
% - This script assumes CAS numbers are stored in a cell array 'all_cas'.
% - Error handling is basic; failed retrievals are marked as 'Error' or 'N/A'.
% - Output is saved in a .mat file and optionally to Excel (commented out).
% - HSPiP output is read from 'Out.dat' generated in the HSPiP directory.
% 
% Author: glsalierno
% Date: April 2025

close all
clear all
clc

% Load CAS numbers from .mat file
% Replace 'CASwFP.mat' with your actual file name if different
% load CAS.mat
all_cas = am(:,1);  % Assuming 'am' is a table or matrix with CAS in first column

% Initialize cell arrays to store results
all_smiles = cell(length(all_cas), 1);
all_hsp_full = cell(length(all_cas), 1);

% Store current working directory to return later
fold = pwd;

% Loop through each CAS to retrieve SMILES and properties via Python script
for i = 1:numel(all_cas)
    tic  % Start timer for this iteration
    
    cas = all_cas{i};
    if ~isempty(cas) && ~strcmp(cas, 'NaN')
        % Construct command to call Python script
        % The Python script fetches SMILES, InChI, IUPAC, and other props from PubChem
        command = sprintf('python get_smiles_InChI_IUPAC_props.py "%s"', cas);
        
        disp(['Executing command: ', command]);  % Debug: Show command
        [status, result] = system(command);
        disp(['Status: ', num2str(status)]);  % Debug: Show execution status
        disp(['Result: ', result]);  % Debug: Show raw output
        
        if status == 0
            all_smiles{i} = strtrim(result);  % Store trimmed result (properties as string)
        else
            all_smiles{i} = 'Error calling Python script';
        end
    else
        all_smiles{i} = 'N/A';
    end
    
    toc  % End timer and display elapsed time
end

% Path to HSPiP installation (update to your actual path)
% Example: hsip_path = 'C:\Path\To\HSPiP\';
hsip_path = 'PATH_TO_HSPIP_INSTALLATION';  % Replace with your HSPiP directory
cd(hsip_path);  % Change to HSPiP directory for CLI execution

% HSPiP CLI command prefix (Y-MBSX mode for SMILES processing)
HSPiPop = 'HSPiP.exe Y-MBSX ';

% Set up import options for reading HSPiP output file 'Out.dat'
opts = delimitedTextImportOptions("NumVariables", 61);
opts.DataLines = [2, Inf];
opts.Delimiter = "\t";
opts.VariableNames = ["SMILES", "Formula", "D", "P", "H", "HDon", "HAcc", "MWt", "Density", "MVol", "Area", "Ovality", "BPt", "MPt", "Tc", "Pc", "Vc", "Zc", "AntA", "AntB", "AntC", "Ant1T", "LogKow", "LogS", "Henry", "LogOHR", "RI", "Hfus", "HvBPt", "Trouton", "RER", "Abra", "Abrb", "EdmiW", "Parachor", "RD", "Cp", "log", "Cond", "SurfTen", "HeavyAtom", "C", "H1", "Br", "Cl", "F", "I", "N", "O", "P1", "S", "Si", "B", "MaxPc", "MinMc", "Sym", "MCI", "Hcomb", "Hform", "Gform", "FGList"];
opts.VariableTypes = ["char", "char", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "char"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts = setvaropts(opts, ["SMILES", "Formula", "FGList"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["SMILES", "Formula", "FGList"], "EmptyFieldRule", "auto");

% Process each SMILES via HSPiP CLI
for i = 1:numel(all_cas)
    tic  % Start timer
    
    smi = all_smiles{i};
    if ~strcmp(smi, 'N/A') && ~strcmp(smi, 'Error calling Python script')
        % Execute HSPiP CLI with SMILES input
        system([HSPiPop, smi]);
        
        % Read generated output file
        all_hsp_full{i} = readtable(fullfile(pwd, 'Out.dat'), opts);
    else
        all_hsp_full{i} = 'N/A';
    end
    
    toc  % End timer
end

% Clean up import options
clear opts

% Return to original directory
cd(fold)

% Save results to .mat file
save('CAS_SMILES_HSP.mat', 'all_cas', 'all_smiles', 'all_hsp_full');

% Optional: Export SMILES to Excel (uncomment if needed)
% xlswrite('SMILESout.xlsx', all_smiles', 'SMILES');
% disp('SMILES data has been added to SMILESout.xlsx');