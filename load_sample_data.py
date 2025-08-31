#!/usr/bin/env python3
"""
Load sample polymer data into Neo4j for testing predictions
"""
import os
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from graph_db.neo4j_connector import Neo4jConnector
import json

def load_sample_polymers():
    """Load sample polymer data with properties"""
    db = Neo4jConnector()
    
    # Sample polymer data with glass transition temperatures
    sample_data = [
        {
            "id": "[CH2]-[CH2]", # Polyethylene pSMILES
            "name": "Polyethylene",
            "smiles": "CC",
            "source": "sample_data",
            "description": "Linear polymer of ethylene",
            "properties": {
                "Tg": -125.0,  # Glass transition temperature in Celsius
                "molecular_weight": 28.05,
                "density": 0.92
            }
        },
        {
            "id": "[CH2]-[CH(CH3)]", # Polypropylene pSMILES  
            "name": "Polypropylene",
            "smiles": "CC(C)",
            "source": "sample_data", 
            "description": "Linear polymer of propylene",
            "properties": {
                "Tg": -10.0,
                "molecular_weight": 42.08,
                "density": 0.90
            }
        },
        {
            "id": "[CH2]-[CH(C6H5)]", # Polystyrene pSMILES
            "name": "Polystyrene", 
            "smiles": "CCc1ccccc1",
            "source": "sample_data",
            "description": "Linear polymer of styrene",
            "properties": {
                "Tg": 100.0,
                "molecular_weight": 104.15,
                "density": 1.05
            }
        },
        {
            "id": "[CH3]-[CH3]", # Ethane (close to input)
            "name": "Ethane-based polymer",
            "smiles": "CC", 
            "source": "sample_data",
            "description": "Simple ethane-based structure",
            "properties": {
                "Tg": -50.0,  # Estimated value
                "molecular_weight": 30.07,
                "density": 0.85
            }
        }
    ]
    
    try:
        # Clear existing data
        print("Clearing existing polymer data...")
        with db.driver.session() as session:
            session.run("MATCH (n) DETACH DELETE n")
        
        # Create polymers and properties
        print("Loading sample polymer data...")
        for polymer in sample_data:
            # Create polymer node
            with db.driver.session() as session:
                create_polymer_query = """
                CREATE (p:Polymer {
                    id: $id,
                    name: $name, 
                    smiles: $smiles,
                    source: $source,
                    description: $description,
                    molecular_weight: $molecular_weight,
                    polar_surface_area: 0,
                    logp: 1.0
                })
                """
                session.run(create_polymer_query, {
                    "id": polymer["id"],
                    "name": polymer["name"],
                    "smiles": polymer["smiles"], 
                    "source": polymer["source"],
                    "description": polymer["description"],
                    "molecular_weight": polymer["properties"]["molecular_weight"]
                })
                
                # Create property nodes and relationships
                for prop_name, prop_value in polymer["properties"].items():
                    create_property_query = """
                    MATCH (p:Polymer {id: $polymer_id})
                    CREATE (prop:Property {
                        name: $prop_name,
                        value: $prop_value,
                        unit: $unit
                    })
                    CREATE (p)-[:HAS_PROPERTY]->(prop)
                    """
                    unit = "°C" if prop_name == "Tg" else ("g/mol" if prop_name == "molecular_weight" else "g/cm³")
                    session.run(create_property_query, {
                        "polymer_id": polymer["id"],
                        "prop_name": prop_name,
                        "prop_value": prop_value,
                        "unit": unit
                    })
        
        # Create similarity relationships (based on structural similarity)
        print("Creating similarity relationships...")
        with db.driver.session() as session:
            # Connect similar polymers
            similarity_pairs = [
                ("[CH3]-[CH3]", "[CH2]-[CH2]"),  # Ethane similar to polyethylene
                ("[CH2]-[CH2]", "[CH2]-[CH(CH3)]"),  # PE similar to PP
                ("[CH2]-[CH(CH3)]", "[CH2]-[CH(C6H5)]")  # PP similar to PS
            ]
            
            for p1, p2 in similarity_pairs:
                create_similarity_query = """
                MATCH (p1:Polymer {id: $id1}), (p2:Polymer {id: $id2})
                CREATE (p1)-[:SIMILAR_TO {similarity: 0.8}]->(p2)
                CREATE (p2)-[:SIMILAR_TO {similarity: 0.8}]->(p1)
                """
                session.run(create_similarity_query, {"id1": p1, "id2": p2})
        
        print("✅ Sample data loaded successfully!")
        
        # Verify data
        print("\nVerifying loaded data:")
        polymers = db.get_all_polymers()
        print(f"Total polymers loaded: {len(polymers)}")
        for polymer in polymers:
            print(f"- {polymer.get('name', 'Unknown')}: {polymer.get('psmiles', 'No ID')}")
            
    except Exception as e:
        print(f"❌ Error loading sample data: {e}")
    finally:
        db.close()

if __name__ == "__main__":
    load_sample_polymers()
