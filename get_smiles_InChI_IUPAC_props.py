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
import logging
from urllib.parse import quote
from pubchem_http import pubchem_get

LOG = logging.getLogger(__name__)
TIMEOUT = 30

def get_compound_info(cas_number):
    """Fetches various properties from PubChem for a CAS number using PUG REST API."""
    base_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug'
    
    try:
        # Step 1: Get CID from CAS number
        cid_url = f"{base_url}/compound/name/{quote(cas_number, safe='')}/cids/JSON"
        response = pubchem_get(cid_url, timeout=TIMEOUT, context=f"CAS {cas_number} CID lookup")
        data = response.json()
        
        if not data.get('IdentifierList', {}).get('CID'):
            return {'Error': f'CAS number {cas_number} not found in PubChem'}
        
        cid = data['IdentifierList']['CID'][0]
        
        # Step 2: Get properties
        properties_url = f'{base_url}/compound/cid/{cid}/property/MolecularFormula,MolecularWeight,ConnectivitySMILES,SMILES,InChI,InChIKey,IUPACName,XLogP,ExactMass,MonoisotopicMass,TPSA,Complexity,Charge,HBondDonorCount,HBondAcceptorCount,RotatableBondCount,HeavyAtomCount,IsotopeAtomCount,AtomStereoCount,DefinedAtomStereoCount,UndefinedAtomStereoCount,BondStereoCount,DefinedBondStereoCount,UndefinedBondStereoCount,CovalentUnitCount,Volume3D,XStericQuadrupole3D,YStericQuadrupole3D,ZStericQuadrupole3D,FeatureCount3D,FeatureAcceptorCount3D,FeatureDonorCount3D,FeatureAnionCount3D,FeatureCationCount3D,FeatureRingCount3D,FeatureHydrophobeCount3D,ConformerModelRMSD3D,EffectiveRotorCount3D,ConformerCount3D/JSON'
        response = pubchem_get(properties_url, timeout=TIMEOUT, context=f"CAS {cas_number} CID {cid} properties")
        data = response.json()
        
        if not data.get('PropertyTable', {}).get('Properties'):
            return {'Error': f'Failed to retrieve properties for CID {cid}'}
        
        properties = data['PropertyTable']['Properties'][0]
        # Preserve the legacy property names used by the MATLAB bridge.
        if 'ConnectivitySMILES' in properties:
            properties.setdefault('CanonicalSMILES', properties['ConnectivitySMILES'])
        if 'SMILES' in properties:
            properties.setdefault('IsomericSMILES', properties['SMILES'])
        
        # Step 3: Get additional data (like boiling point, melting point, etc.)
        sections = ['Physical Properties', 'Chemical and Physical Properties']
        for section in sections:
            section_url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON?heading={quote(section)}'
            try:
                response = pubchem_get(section_url, timeout=TIMEOUT, context=f"CID {cid}: {section}")
                section_data = response.json()
            except (requests.RequestException, ValueError) as e:
                LOG.warning("CAS %s CID %s: optional section %s unavailable: %s", cas_number, cid, section, e)
                continue
            
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
    except requests.exceptions.JSONDecodeError as e:
        LOG.error("CAS %s: failed to parse PubChem JSON: %s", cas_number, e)
        return {'Error': 'Failed to parse JSON response'}
    except requests.HTTPError as e:
        LOG.warning("CAS %s: PubChem HTTP error: %s", cas_number, e)
        if e.response is not None and e.response.status_code == 404:
            return {'Error': f'CAS number {cas_number} or a related record was not found in PubChem'}
        return {'Error': f'PubChem HTTP error: {e}'}
    except requests.RequestException as e:
        LOG.error("CAS %s: PubChem request failed: %s", cas_number, e)
        return {'Error': f'Network error: {e}'}
    except (ValueError, KeyError, IndexError, TypeError, AttributeError) as e:
        LOG.error("CAS %s: invalid PubChem response: %s", cas_number, e)
        return {'Error': 'Invalid PubChem response'}

def format_property(prop_name, prop_value):
    """Formats the property value with its name."""
    return f"{prop_name}: {prop_value}"

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    if len(sys.argv) > 1:
        properties = get_compound_info(sys.argv[1])
        if 'Error' not in properties:
            for name, value in properties.items():
                print(format_property(name, value))
        else:
            print(f"Error: {properties['Error']}")
    else:
        print("Please provide a CAS number as an argument.")
