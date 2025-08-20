# Imports data to Neo4j

def load_data(records, db):
    """
    Load normalized records into Neo4j graph DB.
    - Creates Polymer, Property, and Literature nodes
    - Creates relationships between them
    """
    for rec in records:
        db.run("""
        MERGE (p:Polymer {id: $id, name: $name, psmiles: $psmiles, source: $source})
        """, rec)
        for prop in rec.get('properties', []):
            db.run("""
            MERGE (prop:Property {name: $name, value: $value, units: $units, conditions: $conditions})
            MERGE (p:Polymer {id: $polymer_id})
            MERGE (p)-[:HAS_PROPERTY]->(prop)
            """, {
                "name": prop["name"], "value": prop["value"], "units": prop["units"], "conditions": prop["conditions"], "polymer_id": rec["id"]
            })
        for lit in rec.get('literature', []):
            db.run("""
            MERGE (lit:Literature {doi: $doi, title: $title, authors: $authors})
            MERGE (p:Polymer {id: $polymer_id})
            MERGE (p)-[:CITED_IN]->(lit)
            """, {
                "doi": lit["doi"], "title": lit["title"], "authors": lit["authors"], "polymer_id": rec["id"]
            })
