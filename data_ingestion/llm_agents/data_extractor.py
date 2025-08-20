import os
from dotenv import load_dotenv
import PyPDF2
import re

load_dotenv()
PDF_PATH = os.getenv('PDF_PATH')

# LLM-powered data extraction

def extract_entities(records):
    """Extract named entities (polymers, properties) from PDF using regex."""
    entities = []
    if PDF_PATH and os.path.exists(PDF_PATH):
        with open(PDF_PATH, 'rb') as f:
            reader = PyPDF2.PdfReader(f)
            text = "\n".join(page.extract_text() for page in reader.pages if page.extract_text())
            # Example: extract polymer names
            for match in re.finditer(r'Polymer:\s*(\w+)', text):
                entities.append({'type': 'polymer', 'name': match.group(1)})
            # Example: extract property names and values
            for match in re.finditer(r'Property:\s*(\w+),\s*([\d\.]+),\s*(\w+)', text):
                entities.append({'type': 'property', 'name': match.group(1), 'value': match.group(2), 'units': match.group(3)})
    return entities
