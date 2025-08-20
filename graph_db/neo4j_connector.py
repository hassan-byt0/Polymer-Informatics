# Database interface
from neo4j import GraphDatabase
import os
from dotenv import load_dotenv

load_dotenv()

class Neo4jConnector:
    def __init__(self):
        uri = os.getenv('NEO4J_URI', 'bolt://localhost:7687')
        user = os.getenv('NEO4J_USER', 'neo4j')
        password = os.getenv('NEO4J_PASSWORD', 'test')
        self.driver = GraphDatabase.driver(uri, auth=(user, password))
    def close(self):
        self.driver.close()
    def run(self, query, parameters=None):
        with self.driver.session() as session:
            return session.run(query, parameters or {})
    def get_all_polymers(self):
        query = "MATCH (p:Polymer) RETURN p.id AS psmiles, p.name AS name, p.source AS source, p.description AS description"
        result = self.run(query)
        return [r.data() for r in result]
    def get_polymer(self, polymer_id):
        query = "MATCH (p:Polymer {id: $id}) RETURN p"
        result = self.run(query, {"id": polymer_id})
        record = result.single()
        return record["p"] if record else None
