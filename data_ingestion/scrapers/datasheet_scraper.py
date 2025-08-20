import os
from dotenv import load_dotenv
import PyPDF2

load_dotenv()
PDF_PATH = os.getenv('PDF_PATH')

# Extracts from online datasheets

def scrape_datasheets():
    """Extract polymer datasheets from PDF using PyPDF2."""
    records = []
    if PDF_PATH and os.path.exists(PDF_PATH):
        with open(PDF_PATH, 'rb') as f:
            reader = PyPDF2.PdfReader(f)
            text = "\n".join(page.extract_text() for page in reader.pages if page.extract_text())
            # Simple parsing: look for polymer names and properties
            for line in text.splitlines():
                if 'Polymer:' in line:
                    name = line.split('Polymer:')[1].strip()
                    records.append({
                        'id': name.replace(' ', '_').lower(),
                        'name': name,
                        'psmiles': '',
                        'source': 'datasheet',
                        'properties': [],
                        'literature': []
                    })
                # Add more robust parsing for properties, etc.
    return records
