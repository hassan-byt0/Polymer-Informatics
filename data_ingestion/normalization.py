import re

# Data cleaning & deduplication

def normalize_data(records):
    """
    Normalize and deduplicate polymer records.
    - Standardizes property names and units
    - Removes duplicates by polymer id
    - Cleans up missing/invalid fields
    """
    seen = set()
    normalized = []
    for rec in records:
        rec['id'] = rec.get('id', '').replace(' ', '_').lower()
        if rec['id'] in seen:
            continue
        seen.add(rec['id'])
        # Standardize property names/units
        for prop in rec.get('properties', []):
            prop['name'] = prop['name'].strip().lower().replace(' ', '_')
            prop['units'] = prop.get('units', '').strip()
            # Convert temperature units to Celsius if needed
            if prop['units'].lower() in ['k', 'kelvin']:
                try:
                    prop['value'] = float(prop['value']) - 273.15
                    prop['units'] = 'C'
                except Exception:
                    pass
        normalized.append(rec)
    return normalized
