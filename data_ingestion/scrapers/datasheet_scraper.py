import os
import re
from dotenv import load_dotenv
import PyPDF2

load_dotenv()
PDF_PATH = os.getenv('PDF_PATH')

def scrape_datasheets():
    """Extract polymer datasheets from PDF using PyPDF2 with enhanced parsing."""
    records = []
    if PDF_PATH and os.path.exists(PDF_PATH):
        with open(PDF_PATH, 'rb') as f:
            reader = PyPDF2.PdfReader(f)
            text = "\n".join(page.extract_text() for page in reader.pages if page.extract_text())
            
            # Enhanced parsing for polymer data
            records = parse_polymer_data(text)
            
    return records

def parse_polymer_data(text):
    """Parse polymer data from extracted PDF text."""
    records = []
    
    # Split text into lines and clean up
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    
    # Common polymer name patterns
    polymer_patterns = [
        r'poly\w+',  # polyethylene, polystyrene, etc.
        r'PE\b',     # Polyethylene abbreviation
        r'PP\b',     # Polypropylene abbreviation
        r'PS\b',     # Polystyrene abbreviation
        r'PVC\b',    # Polyvinyl chloride
        r'PMMA\b',   # Polymethyl methacrylate
        r'PET\b',    # Polyethylene terephthalate
        r'PA\b',     # Polyamide (Nylon)
        r'PC\b',     # Polycarbonate
        r'ABS\b',    # Acrylonitrile butadiene styrene
    ]
    
    # Property patterns (with units)
    property_patterns = {
        'Tg': r'(?:Tg|glass\s+transition|T_g)\s*[:\-=]?\s*(-?\d+(?:\.\d+)?)\s*°?C?',
        'Tm': r'(?:Tm|melting\s+point|T_m)\s*[:\-=]?\s*(-?\d+(?:\.\d+)?)\s*°?C?',
        'density': r'(?:density|ρ)\s*[:\-=]?\s*(\d+(?:\.\d+)?)\s*g/cm³?',
        'molecular_weight': r'(?:MW|molecular\s+weight|M_w)\s*[:\-=]?\s*(\d+(?:\.\d+)?)\s*(?:g/mol|kDa)?',
        'tensile_strength': r'(?:tensile\s+strength|σ_t)\s*[:\-=]?\s*(\d+(?:\.\d+)?)\s*MPa?',
        'elastic_modulus': r'(?:elastic\s+modulus|Young.s\s+modulus|E)\s*[:\-=]?\s*(\d+(?:\.\d+)?)\s*(?:GPa|MPa)?'
    }
    
    current_polymer = None
    
    for i, line in enumerate(lines):
        # Look for polymer names
        for pattern in polymer_patterns:
            matches = re.finditer(pattern, line, re.IGNORECASE)
            for match in matches:
                polymer_name = match.group(0)
                
                # Get more context around the match
                context_start = max(0, i-2)
                context_end = min(len(lines), i+10)
                context_lines = lines[context_start:context_end]
                
                # Extract properties from context
                properties = extract_properties(' '.join(context_lines), property_patterns)
                
                # Generate polymer record
                polymer_record = {
                    'id': generate_polymer_id(polymer_name),
                    'name': format_polymer_name(polymer_name),
                    'psmiles': get_polymer_smiles(polymer_name),
                    'source': 'pdf_datasheet',
                    'description': f'Polymer extracted from datasheet: {polymer_name}',
                    'properties': properties,
                    'literature': [],
                    'raw_text': ' '.join(context_lines[:5])  # Store some context
                }
                
                # Only add if we found some properties
                if properties:
                    records.append(polymer_record)
    
    # Remove duplicates based on polymer ID
    unique_records = {}
    for record in records:
        polymer_id = record['id']
        if polymer_id not in unique_records:
            unique_records[polymer_id] = record
        else:
            # Merge properties if duplicate found
            existing = unique_records[polymer_id]
            for prop, value in record['properties'].items():
                if prop not in existing['properties']:
                    existing['properties'][prop] = value
    
    return list(unique_records.values())

def extract_properties(text, property_patterns):
    """Extract property values from text using regex patterns."""
    properties = {}
    
    for prop_name, pattern in property_patterns.items():
        matches = re.finditer(pattern, text, re.IGNORECASE)
        for match in matches:
            try:
                value = float(match.group(1))
                properties[prop_name] = value
                break  # Take first match
            except (ValueError, IndexError):
                continue
    
    return properties

def generate_polymer_id(polymer_name):
    """Generate a standardized ID for the polymer."""
    # Clean and standardize the name
    clean_name = re.sub(r'[^a-zA-Z0-9]', '_', polymer_name.lower())
    return clean_name

def format_polymer_name(polymer_name):
    """Format polymer name for display."""
    name_mapping = {
        'pe': 'Polyethylene',
        'pp': 'Polypropylene', 
        'ps': 'Polystyrene',
        'pvc': 'Polyvinyl Chloride',
        'pmma': 'Polymethyl Methacrylate',
        'pet': 'Polyethylene Terephthalate',
        'pa': 'Polyamide',
        'pc': 'Polycarbonate',
        'abs': 'Acrylonitrile Butadiene Styrene'
    }
    
    clean_name = polymer_name.lower().strip()
    return name_mapping.get(clean_name, polymer_name.title())

def get_polymer_smiles(polymer_name):
    """Get approximate SMILES for common polymers."""
    smiles_mapping = {
        'polyethylene': 'CC',
        'pe': 'CC',
        'polypropylene': 'CC(C)',
        'pp': 'CC(C)',
        'polystyrene': 'CCc1ccccc1',
        'ps': 'CCc1ccccc1',
        'pvc': 'CCCl',
        'pmma': 'CC(C)(C(=O)OC)',
        'pet': 'CCOc1ccc(C(=O)O)cc1',
        'pa': 'CCCCCCN',
        'pc': 'CC(C)(c1ccc(O)cc1)',
        'abs': 'CCc1ccccc1'  # Simplified
    }
    
    clean_name = polymer_name.lower().strip()
    return smiles_mapping.get(clean_name, 'CC')  # Default to simple alkane
