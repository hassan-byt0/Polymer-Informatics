#!/usr/bin/env python3
"""
Enhanced data loader that scrapes PDF data and loads it into Neo4j
"""
import os
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from graph_db.neo4j_connector import Neo4jConnector
from data_ingestion.scrapers.datasheet_scraper import scrape_datasheets
import json

def load_pdf_data():
    """Load polymer data from PDF into Neo4j"""
    print("🔄 Starting PDF data extraction and loading...")
    
    # Initialize Neo4j connection
    db = Neo4jConnector()
    
    try:
        # Extract data from PDF
        print("📄 Extracting data from PDF...")
        pdf_records = scrape_datasheets()
        print(f"✅ Extracted {len(pdf_records)} polymer records from PDF")
        
        if not pdf_records:
            print("⚠️ No polymer data found in PDF. Loading sample data instead...")
            load_sample_data(db)
            return
        
        # Clear existing data
        print("🗑️ Clearing existing polymer data...")
        with db.driver.session() as session:
            session.run("MATCH (n) DETACH DELETE n")
        
        # Load PDF data into Neo4j
        print("💾 Loading PDF data into Neo4j...")
        for record in pdf_records:
            print(f"  Loading: {record['name']} with {len(record['properties'])} properties")
            load_polymer_record(db, record)
        
        # Create similarity relationships
        print("🔗 Creating similarity relationships...")
        create_similarity_relationships(db, pdf_records)
        
        # Verify loaded data
        print("\n✅ PDF data loaded successfully!")
        verify_data(db)
        
    except Exception as e:
        print(f"❌ Error loading PDF data: {e}")
        print("🔄 Falling back to sample data...")
        load_sample_data(db)
    finally:
        db.close()

def load_polymer_record(db, record):
    """Load a single polymer record into Neo4j"""
    with db.driver.session() as session:
        # Create polymer node
        create_polymer_query = """
        CREATE (p:Polymer {
            id: $id,
            name: $name,
            psmiles: $psmiles,
            smiles: $smiles,
            source: $source,
            description: $description,
            molecular_weight: $molecular_weight,
            polar_surface_area: 0,
            logp: 1.0
        })
        """
        session.run(create_polymer_query, {
            "id": record["id"],
            "name": record["name"],
            "psmiles": record.get("psmiles", ""),
            "smiles": record.get("psmiles", "CC"),  # Use psmiles as fallback
            "source": record["source"],
            "description": record["description"],
            "molecular_weight": record["properties"].get("molecular_weight", 100.0)
        })
        
        # Create property nodes and relationships
        for prop_name, prop_value in record["properties"].items():
            create_property_query = """
            MATCH (p:Polymer {id: $polymer_id})
            CREATE (prop:Property {
                name: $prop_name,
                value: $prop_value,
                unit: $unit
            })
            CREATE (p)-[:HAS_PROPERTY]->(prop)
            """
            unit = get_property_unit(prop_name)
            session.run(create_property_query, {
                "polymer_id": record["id"],
                "prop_name": prop_name,
                "prop_value": prop_value,
                "unit": unit
            })

def get_property_unit(prop_name):
    """Get the unit for a property"""
    unit_mapping = {
        "Tg": "°C",
        "Tm": "°C", 
        "density": "g/cm³",
        "molecular_weight": "g/mol",
        "tensile_strength": "MPa",
        "elastic_modulus": "GPa"
    }
    return unit_mapping.get(prop_name, "")

def create_similarity_relationships(db, records):
    """Create similarity relationships between polymers"""
    with db.driver.session() as session:
        # Simple similarity based on polymer type
        similarity_groups = {
            'polyolefin': ['polyethylene', 'pe', 'polypropylene', 'pp'],
            'aromatic': ['polystyrene', 'ps', 'pet', 'pc'],
            'vinyl': ['pvc', 'pmma'],
            'engineering': ['pa', 'pc', 'abs']
        }
        
        for group_name, polymers in similarity_groups.items():
            # Connect polymers within the same group
            for i, p1 in enumerate(polymers):
                for p2 in polymers[i+1:]:
                    create_similarity_query = """
                    MATCH (p1:Polymer), (p2:Polymer)
                    WHERE LOWER(p1.name) CONTAINS $p1 OR LOWER(p1.id) CONTAINS $p1
                    AND LOWER(p2.name) CONTAINS $p2 OR LOWER(p2.id) CONTAINS $p2
                    AND p1 <> p2
                    CREATE (p1)-[:SIMILAR_TO {similarity: 0.8, group: $group}]->(p2)
                    CREATE (p2)-[:SIMILAR_TO {similarity: 0.8, group: $group}]->(p1)
                    """
                    session.run(create_similarity_query, {
                        "p1": p1,
                        "p2": p2,
                        "group": group_name
                    })

def load_sample_data(db):
    """Fallback to load sample data if PDF extraction fails"""
    sample_data = [
        {
            "id": "polyethylene",
            "name": "Polyethylene",
            "psmiles": "[CH2]-[CH2]",
            "source": "sample_data",
            "description": "Linear polymer of ethylene",
            "properties": {
                "Tg": -125.0,
                "molecular_weight": 28.05,
                "density": 0.92
            }
        },
        {
            "id": "polypropylene", 
            "name": "Polypropylene",
            "psmiles": "[CH2]-[CH(CH3)]",
            "source": "sample_data",
            "description": "Linear polymer of propylene",
            "properties": {
                "Tg": -10.0,
                "molecular_weight": 42.08,
                "density": 0.90
            }
        },
        {
            "id": "polystyrene",
            "name": "Polystyrene",
            "psmiles": "[CH2]-[CH(C6H5)]",
            "source": "sample_data", 
            "description": "Linear polymer of styrene",
            "properties": {
                "Tg": 100.0,
                "molecular_weight": 104.15,
                "density": 1.05
            }
        }
    ]
    
    print("📝 Loading sample polymer data...")
    for record in sample_data:
        load_polymer_record(db, record)
    
    # Create basic similarity relationships
    with db.driver.session() as session:
        similarity_pairs = [
            ("polyethylene", "polypropylene"),
            ("polypropylene", "polystyrene")
        ]
        
        for p1, p2 in similarity_pairs:
            create_similarity_query = """
            MATCH (p1:Polymer {id: $id1}), (p2:Polymer {id: $id2})
            CREATE (p1)-[:SIMILAR_TO {similarity: 0.7}]->(p2)
            CREATE (p2)-[:SIMILAR_TO {similarity: 0.7}]->(p1)
            """
            session.run(create_similarity_query, {"id1": p1, "id2": p2})

def verify_data(db):
    """Verify the loaded data"""
    polymers = db.get_all_polymers()
    print(f"\n📊 Database Summary:")
    print(f"   Total polymers loaded: {len(polymers)}")
    
    for polymer in polymers[:10]:  # Show first 10
        polymer_name = polymer.get('name', 'Unknown')
        polymer_id = polymer.get('psmiles', polymer.get('id', 'No ID'))
        print(f"   - {polymer_name}: {polymer_id}")
    
    if len(polymers) > 10:
        print(f"   ... and {len(polymers) - 10} more")

if __name__ == "__main__":
    load_pdf_data()
