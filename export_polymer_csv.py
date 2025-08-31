#!/usr/bin/env python3
"""
Export scraped polymer data to CSV for review and reference
"""
import os
import sys
import csv
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from data_ingestion.scrapers.datasheet_scraper import scrape_datasheets

def export_to_csv():
    """Export scraped polymer data to CSV file"""
    
    print("Scraping polymer data from PDF...")
    records = scrape_datasheets()
    
    if not records:
        print("No polymer data found in PDF")
        return
    
    csv_filename = "scraped_polymer_data.csv"
    
    # Prepare CSV data
    csv_data = []
    
    for record in records:
        row = {
            'id': record['id'],
            'name': record['name'], 
            'psmiles': record['psmiles'],
            'smiles': record.get('smiles', ''),
            'source': record['source'],
            'description': record['description']
        }
        
        # Add all properties as separate columns
        properties = record.get('properties', {})
        for prop_name, prop_value in properties.items():
            row[f'{prop_name}'] = prop_value
            
        # Add raw text for reference
        row['raw_text_sample'] = record.get('raw_text', '')[:200] + '...' if record.get('raw_text') else ''
        
        csv_data.append(row)
    
    # Get all unique property names for column headers
    all_properties = set()
    for record in records:
        all_properties.update(record.get('properties', {}).keys())
    
    # Define CSV column order
    base_columns = ['id', 'name', 'psmiles', 'smiles', 'source', 'description']
    property_columns = sorted(list(all_properties))
    all_columns = base_columns + property_columns + ['raw_text_sample']
    
    # Write CSV file
    with open(csv_filename, 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=all_columns)
        writer.writeheader()
        
        for row in csv_data:
            # Ensure all columns exist in row
            for col in all_columns:
                if col not in row:
                    row[col] = ''
            writer.writerow(row)
    
    print(f"✅ Exported {len(records)} polymer records to {csv_filename}")
    
    # Print summary statistics
    print(f"\n📊 Summary Statistics:")
    print(f"Total polymers: {len(records)}")
    
    # Count properties
    property_counts = {}
    for record in records:
        for prop_name in record.get('properties', {}).keys():
            property_counts[prop_name] = property_counts.get(prop_name, 0) + 1
    
    print(f"\nProperty distribution:")
    for prop, count in sorted(property_counts.items(), key=lambda x: x[1], reverse=True):
        print(f"  {prop}: {count} polymers")
    
    # Show sample records
    print(f"\n🔍 Sample Records (first 3):")
    for i, record in enumerate(records[:3]):
        print(f"\n{i+1}. {record['name']} ({record['id']})")
        print(f"   pSMILES: {record['psmiles']}")
        print(f"   Properties: {list(record.get('properties', {}).keys())}")
        print(f"   Sample text: {record.get('raw_text', '')[:100]}...")

if __name__ == "__main__":
    export_to_csv()
