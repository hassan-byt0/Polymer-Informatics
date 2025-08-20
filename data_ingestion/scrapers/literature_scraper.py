import os
from dotenv import load_dotenv
import PyPDF2

load_dotenv()
PDF_PATH = os.getenv('PDF_PATH')

# Processes research papers

def scrape_literature():
    """Extract polymer property data from PDF using PyPDF2."""
    records = []
    if PDF_PATH and os.path.exists(PDF_PATH):
        with open(PDF_PATH, 'rb') as f:
            reader = PyPDF2.PdfReader(f)
            text = "\n".join(page.extract_text() for page in reader.pages if page.extract_text())
            # Simple parsing: look for property tables or lines
            for line in text.splitlines():
                if 'Property:' in line:
                    parts = line.split('Property:')[1].split(',')
                    name = parts[0].strip() if len(parts) > 0 else ''
                    value = parts[1].strip() if len(parts) > 1 else ''
                    units = parts[2].strip() if len(parts) > 2 else ''
                    records.append({
                        'id': name.replace(' ', '_').lower(),
                        'name': name,
                        'psmiles': '',
                        'source': 'literature',
                        'properties': [{
                            'name': name,
                            'value': value,
                            'units': units,
                            'conditions': ''
                        }],
                        'literature': []
                    })
                # Add more robust parsing for literature references, etc.
    return records
