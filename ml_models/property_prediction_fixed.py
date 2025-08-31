import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.gaussian_process import GaussianProcessRegressor

class PropertyPredictor:
    def __init__(self, graph_db):
        self.db = graph_db  # Neo4j connection
        self.models = {
            'random_forest': RandomForestRegressor(n_estimators=100, random_state=42),
            'gaussian_process': GaussianProcessRegressor(random_state=42)
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
        try:
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
            
            # If no exact match, get similar polymers with the property
            with self.db.driver.session() as session:
                similar_query = """
                MATCH (target:Polymer {id: $id})
                MATCH (target)-[:SIMILAR_TO*0..2]-(similar:Polymer)-[:HAS_PROPERTY]->(prop:Property {name: $prop_name})
                RETURN similar, prop.value as value, similar.molecular_weight as mw, 
                       similar.polar_surface_area as psa, similar.logp as logp
                """
                result = session.run(similar_query, {"id": polymer_id, "prop_name": property_name})
                training_data = [record.data() for record in result]
            
            if not training_data:
                # Fallback: get any polymers with this property
                with self.db.driver.session() as session:
                    fallback_query = """
                    MATCH (p:Polymer)-[:HAS_PROPERTY]->(prop:Property {name: $prop_name})
                    RETURN p, prop.value as value, p.molecular_weight as mw,
                           p.polar_surface_area as psa, p.logp as logp
                    """
                    result = session.run(fallback_query, {"prop_name": property_name})
                    training_data = [record.data() for record in result]
                    
                    if training_data:
                        # Simple average for fallback
                        values = [float(d["value"]) for d in training_data]
                        return sum(values) / len(values)
                return None
            
            # Extract features and targets
            X = []
            y = []
            for d in training_data:
                try:
                    features = [
                        float(d.get('mw', 0)),
                        float(d.get('psa', 0)), 
                        float(d.get('logp', 0))
                    ]
                    X.append(features)
                    y.append(float(d['value']))
                except (ValueError, TypeError):
                    continue
            
            if not X or not y:
                return None
            
            # Use simple average if too few samples
            if len(X) < 3:
                return sum(y) / len(y)
            
            # Use ML model for prediction
            if len(training_data) > 10:
                model = self.models['random_forest']
            else:
                model = self.models['gaussian_process']
            
            model.fit(X, y)
            
            # Get target polymer features
            target_polymer = self.db.get_polymer(polymer_id)
            if not target_polymer:
                return None
            
            target_features = [
                float(target_polymer.get('molecular_weight', 0)),
                float(target_polymer.get('polar_surface_area', 0)),
                float(target_polymer.get('logp', 0))
            ]
            
            return float(model.predict([target_features])[0])
            
        except Exception as e:
            print(f"Error in property prediction: {e}")
            return None
    
    def _extract_features(self, polymer):
        """Convert polymer to feature vector"""
        return [
            polymer.get('molecular_weight', 0),
            polymer.get('polar_surface_area', 0),
            polymer.get('logp', 0)
        ]
