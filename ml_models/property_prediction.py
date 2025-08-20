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
        return self.db.run(query).data()
    
    def predict(self, polymer_id, property_name):
        """Predict property using similar polymers"""
        training_data = self.retrieve_similar(polymer_id)
        
        if not training_data:
            return None
        
        if len(training_data) > 50:
            model = self.models['random_forest']
        else:
            model = self.models['gaussian_process']
        
        X = [self._extract_features(d['neighbor']) for d in training_data]
        y = [d['nprop'][property_name] for d in training_data if property_name in d['nprop']]
        
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
