#!/usr/bin/env python3
import csv
import sys
import os
sys.path.append('/app')
from data_ingestion.scrapers.datasheet_scraper import scrape_datasheets

# Scrape data
records = scrape_datasheets()
print(f"Found {len(records)} polymer records")

# Export to CSV
with open('/app/scraped_polymer_data.csv', 'w', newline='', encoding='utf-8') as f:
    if records:
        # Get all unique properties
        all_props = set()
        for r in records:
            all_props.update(r.get('properties', {}).keys())
        
        fieldnames = ['id', 'name', 'psmiles', 'source', 'description'] + sorted(all_props)
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        
        for record in records:
            row = {
                'id': record['id'],
                'name': record['name'],
                'psmiles': record['psmiles'], 
                'source': record['source'],
                'description': record['description']
            }
            # Add properties
            for prop, val in record.get('properties', {}).items():
                row[prop] = val
            writer.writerow(row)

print("CSV exported successfully!")
