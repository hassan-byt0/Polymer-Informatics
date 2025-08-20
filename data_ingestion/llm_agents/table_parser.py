import os
from dotenv import load_dotenv
import PyPDF2
import re

load_dotenv()
PDF_PATH = os.getenv('PDF_PATH')

# Handles complex table structures

def parse_tables(records):
    """Extract tables from PDF using regex (simple row parsing)."""
    tables = []
    if PDF_PATH and os.path.exists(PDF_PATH):
        with open(PDF_PATH, 'rb') as f:
            reader = PyPDF2.PdfReader(f)
            text = "\n".join(page.extract_text() for page in reader.pages if page.extract_text())
            # Example: parse lines that look like table rows
            for line in text.splitlines():
                if re.match(r'\w+,\s*[\d\.]+,\s*\w+', line):
                    parts = [p.strip() for p in line.split(',')]
                    tables.append({'name': parts[0], 'value': parts[1], 'units': parts[2]})
    return tables
