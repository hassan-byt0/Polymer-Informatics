from flask import Flask, request, jsonify
from polymer_representation import to_psmiles
from search import HybridSearch

app = Flask(__name__)
search_engine = HybridSearch()

@app.route('/predict', methods=['POST'])
def predict_property():
    data = request.json
    polymer_smiles = data['smiles']
    property_name = data['property']
    
    # Convert to p-SMILES
    mol = Chem.MolFromSmiles(polymer_smiles)
    psmiles = to_psmiles(mol)
    
    # Retrieve similar polymers
    similar = search_engine.find_similar(
        psmiles, 
        f"polymers with {property_name}"
    )
    
    # Predict property
    prediction = property_predictor.predict(psmiles, property_name)
    
    return jsonify({
        'prediction': prediction,
        'similar_polymers': similar[:5],
        'confidence': 0.92  # Example value
    })

@app.route('/visualize', methods=['POST'])
def visualize_structure():
    smiles = request.json['smiles']
    mol = Chem.MolFromSmiles(smiles)
    img = Chem.Draw.MolToImage(mol, size=(400, 300))
    return img
