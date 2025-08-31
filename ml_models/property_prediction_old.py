import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.gaussian_process import GaussianProcessRegressor

class PropertyPredictor:
    def __init__(self, graph_db):
        self.db = graph_db  # Neo4j connection
        self.models = {
            'random_forest': RandomForestRegressor(n_estimators=100),
            'gaussian_process': GaussianProcessRegressor()
        }
    
    def retrieve_similar(self, polymer_id, radius=2):
        """Retrieve chemically similar polymers from graph DB"""
        query = f"""
        MATCH (p:Polymer {{id: '{polymer_id}'}})-[:HAS_PROPERTY]->(prop)
        WITH p, prop
        MATCH (p)-[:SIMILAR_TO*0..{radius}]-(neighbor)-[:HAS_PROPERTY]->(nprop)
        RETURN neighbor, nprop
        """
        return self.db.run(query)
    
    def predict(self, polymer_id, property_name):
        """Predict property using similar polymers"""
        # First try to find exact match
        with self.db.driver.session() as session:
            exact_match_query = """
            MATCH (p:Polymer {id: $id})-[:HAS_PROPERTY]->(prop:Property {name: $prop_name})
            RETURN prop.value as value
            """
            result = session.run(exact_match_query, {"id": polymer_id, "prop_name": property_name})
            record = result.single()
            if record:
                return float(record["value"])
        
        # If no exact match, try similar polymers
        training_data = self.retrieve_similar(polymer_id)
        
        if not training_data:
            # Fallback: get any polymers with this property
            with self.db.driver.session() as session:
                fallback_query = """
                MATCH (p:Polymer)-[:HAS_PROPERTY]->(prop:Property {name: $prop_name})
                RETURN p, prop.value as value
                """
                result = session.run(fallback_query, {"prop_name": property_name})
                fallback_data = [{"neighbor": record["p"], "value": record["value"]} for record in result]
                if fallback_data:
                    # Simple average for fallback
                    values = [float(d["value"]) for d in fallback_data]
                    return sum(values) / len(values)
            return None
        
        # Use ML model for prediction
        if len(training_data) > 50:
            model = self.models['random_forest']
        else:
            model = self.models['gaussian_process']
        
        # Extract features and target values
        X, y = [], []
        for d in training_data:
            if 'neighbor' in d:
                neighbor = d['neighbor']
                # Look for property value in the neighbor
                with self.db.driver.session() as session:
                    prop_query = """
                    MATCH (p:Polymer {id: $id})-[:HAS_PROPERTY]->(prop:Property {name: $prop_name})
                    RETURN prop.value as value
                    """
                    result = session.run(prop_query, {"id": neighbor.get("id"), "prop_name": property_name})
                    record = result.single()
                    if record:
                        X.append(self._extract_features(neighbor))
                        y.append(float(record["value"]))
        
        if not X or not y:
            return None
        
        model.fit(X, y)
        target_polymer = self.db.get_polymer(polymer_id)
        
        if not target_polymer:
            return None
        
        target_features = self._extract_features(target_polymer)
        return model.predict([target_features])[0]
    
    def _extract_features(self, polymer):
        """Convert polymer to feature vector"""
        return [
            polymer.get('molecular_weight', 0),
            polymer.get('polar_surface_area', 0),
            polymer.get('logp', 0)
            # Add more features
        ]
