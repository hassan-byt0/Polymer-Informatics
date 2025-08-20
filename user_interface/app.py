from flask import Flask, request, jsonify, send_file
from io import BytesIO
from rdkit import Chem
from rdkit.Chem import Draw
from polymer_representation.psmiles import to_psmiles
from polymer_representation.graph_encoder import to_graph
from polymer_representation.bigsmiles import to_bigsmiles
from search.hybrid_search import HybridSearch
from ml_models.property_prediction import PropertyPredictor
from ml_models.uncertainty import estimate_uncertainty
from ml_models.active_learning import suggest_next_experiment
from graph_db.neo4j_connector import Neo4jConnector
from graph_db.data_loader import load_data
from data_ingestion.scrapers.datasheet_scraper import scrape_datasheets
from data_ingestion.scrapers.literature_scraper import scrape_literature
from data_ingestion.llm_agents.data_extractor import extract_entities
from data_ingestion.llm_agents.table_parser import parse_tables
from data_ingestion.normalization import normalize_data
from utils.cheminformatics import compute_descriptors
from utils.visualization import render_molecule

app = Flask(__name__)
db = Neo4jConnector()
search_engine = HybridSearch()
property_predictor = PropertyPredictor(db)

@app.route('/ingest', methods=['POST'])
def ingest_data():
    """Ingest data from datasheets and literature, normalize, and load to graph DB."""
    sources = request.json.get('sources', [])
    raw_data = []
    if 'datasheets' in sources:
        raw_data += scrape_datasheets()
    if 'literature' in sources:
        raw_data += scrape_literature()
    entities = extract_entities(raw_data)
    tables = parse_tables(raw_data)
    normalized = normalize_data(entities + tables)
    load_data(normalized, db)
    return jsonify({'status': 'success', 'records_loaded': len(normalized)})

@app.route('/predict', methods=['POST'])
def predict_property():
    data = request.json
    polymer_smiles = data['smiles']
    property_name = data['property']
    mol = Chem.MolFromSmiles(polymer_smiles)
    psmiles = to_psmiles(mol)
    polymer_id = psmiles  # Use pSMILES as ID for now
    prediction = property_predictor.predict(polymer_id, property_name)
    uncertainty = estimate_uncertainty(polymer_id, property_name)
    similar = search_engine.hybrid_search(
        query_psmiles=psmiles,
        query_text=f"polymers with {property_name}",
        targets=db.get_all_polymers(),
        alpha=0.7
    )
    return jsonify({
        'prediction': prediction,
        'uncertainty': uncertainty,
        'similar_polymers': [s[0] for s in similar[:5]],
        'confidence': 1 - uncertainty
    })

@app.route('/visualize', methods=['POST'])
def visualize_structure():
    smiles = request.json['smiles']
    mol = Chem.MolFromSmiles(smiles)
    img = render_molecule(mol)
    buf = BytesIO()
    img.save(buf, format='PNG')
    buf.seek(0)
    return send_file(buf, mimetype='image/png')

@app.route('/represent', methods=['POST'])
def represent_polymer():
    smiles = request.json['smiles']
    mol = Chem.MolFromSmiles(smiles)
    psmiles = to_psmiles(mol)
    bigsmiles = to_bigsmiles(mol)
    graph = to_graph(mol)
    return jsonify({'psmiles': psmiles, 'bigsmiles': bigsmiles, 'graph': graph})

@app.route('/descriptors', methods=['POST'])
def get_descriptors():
    smiles = request.json['smiles']
    mol = Chem.MolFromSmiles(smiles)
    desc = compute_descriptors(mol)
    return jsonify(desc)

@app.route('/active_learning', methods=['POST'])
def active_learning():
    property_name = request.json['property']
    suggestion = suggest_next_experiment(property_name, db)
    return jsonify({'next_experiment': suggestion})

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000)
