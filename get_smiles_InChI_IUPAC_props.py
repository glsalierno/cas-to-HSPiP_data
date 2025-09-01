# get_smiles_InChI_IUPAC_props.py
#
# This script fetches chemical properties from PubChem for a given CAS number,
# including SMILES, InChI, IUPAC name, and various physical/chemical properties.
# It uses the PubChem PUG REST API.
#
# Requirements:
# - Python 3.x
# - requests library (install via: pip install requests)
# - json library (built-in)
#
# Usage:
# python get_smiles_InChI_IUPAC_props.py <CAS_number>
#
# Example:
# python get_smiles_InChI_IUPAC_props.py 50-00-0
#
# Output:
# Prints properties in "Name: Value" format, one per line.
# If an error occurs, prints "Error: <message>".
#
# Notes:
# - Handles network errors and JSON parsing issues.
# - Fetches a wide range of properties; see PubChem API docs for details.
# - Intended to be called from external scripts (e.g., MATLAB) for automation.
#
# Author: glsalierno
# Date: April 2025

import sys
import requests
import json

def get_compound_info(cas_number):
    """Fetches various properties from PubChem for a CAS number using PUG REST API."""
    base_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug'
    
    try:
        # Step 1: Get CID from CAS number
        cid_url = f'{base_url}/compound/name/{cas_number}/cids/JSON'
        response = requests.get(cid_url)
        data = response.json()
        
        if 'IdentifierList' not in data or 'CID' not in data['IdentifierList']:
            return {'Error': f'CAS number {cas_number} not found in PubChem'}
        
        cid = data['IdentifierList']['CID'][0]
        
        # Step 2: Get properties
        properties_url = f'{base_url}/compound/cid/{cid}/property/MolecularFormula,MolecularWeight,CanonicalSMILES,IsomericSMILES,InChI,InChIKey,IUPACName,XLogP,ExactMass,MonoisotopicMass,TPSA,Complexity,Charge,HBondDonorCount,HBondAcceptorCount,RotatableBondCount,HeavyAtomCount,IsotopeAtomCount,AtomStereoCount,DefinedAtomStereoCount,UndefinedAtomStereoCount,BondStereoCount,DefinedBondStereoCount,UndefinedBondStereoCount,CovalentUnitCount,Volume3D,XStericQuadrupole3D,YStericQuadrupole3D,ZStericQuadrupole3D,FeatureCount3D,FeatureAcceptorCount3D,FeatureDonorCount3D,FeatureAnionCount3D,FeatureCationCount3D,FeatureRingCount3D,FeatureHydrophobeCount3D,ConformerModelRMSD3D,EffectiveRotorCount3D,ConformerCount3D/JSON'
        response = requests.get(properties_url)
        data = response.json()
        
        if 'PropertyTable' not in data or 'Properties' not in data['PropertyTable']:
            return {'Error': f'Failed to retrieve properties for CID {cid}'}
        
        properties = data['PropertyTable']['Properties'][0]
        
        # Step 3: Get additional data (like boiling point, melting point, etc.)
        sections = ['Physical Properties', 'Chemical and Physical Properties']
        for section in sections:
            section_url = f'{base_url}/compound/cid/{cid}/JSON?heading={section}'
            response = requests.get(section_url)
            section_data = response.json()
            
            if 'Record' in section_data and 'Section' in section_data['Record']:
                for s in section_data['Record']['Section']:
                    if s['TOCHeading'] == section:
                        for subsection in s.get('Section', []):
                            for info in subsection.get('Information', []):
                                name = info.get('Name', '')
                                value = info.get('Value', {}).get('StringWithMarkup', [{}])[0].get('String', '')
                                if name and value:
                                    properties[name] = value
        
        return properties
    except requests.RequestException as e:
        return {'Error': f'Network error: {str(e)}'}
    except json.JSONDecodeError:
        return {'Error': 'Failed to parse JSON response'}
    except Exception as e:
        return {'Error': f'Unexpected error: {str(e)}'}

def format_property(prop_name, prop_value):
    """Formats the property value with its name."""
    return f"{prop_name}: {prop_value}"

if __name__ == "__main__":
    if len(sys.argv) > 1:
        properties = get_compound_info(sys.argv[1])
        if 'Error' not in properties:
            for name, value in properties.items():
                print(format_property(name, value))
        else:
            print(f"Error: {properties['Error']}")
    else:
        print("Please provide a CAS number as an argument.")
