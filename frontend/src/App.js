import React, { useState } from 'react';
import { ingestData, predictProperty, visualizeStructure, getDescriptors, representPolymer, activeLearning } from './api';

function App() {
  const [smiles, setSmiles] = useState('');
  const [property, setProperty] = useState('Tg');
  const [prediction, setPrediction] = useState(null);
  const [imageUrl, setImageUrl] = useState(null);
  const [descriptors, setDescriptors] = useState(null);
  const [representation, setRepresentation] = useState(null);
  const [activeSuggestion, setActiveSuggestion] = useState(null);
  const [loading, setLoading] = useState(false);

  const handlePredict = async () => {
    setLoading(true);
    try {
      const res = await predictProperty(smiles, property);
      setPrediction(res.data);
    } catch (error) {
      console.error('Prediction error:', error);
    }
    setLoading(false);
  };

  const handleVisualize = async () => {
    setLoading(true);
    try {
      const res = await visualizeStructure(smiles);
      const url = URL.createObjectURL(res.data);
      setImageUrl(url);
    } catch (error) {
      console.error('Visualization error:', error);
    }
    setLoading(false);
  };

  const handleDescriptors = async () => {
    setLoading(true);
    try {
      const res = await getDescriptors(smiles);
      setDescriptors(res.data);
    } catch (error) {
      console.error('Descriptors error:', error);
    }
    setLoading(false);
  };

  const handleRepresent = async () => {
    setLoading(true);
    try {
      const res = await representPolymer(smiles);
      setRepresentation(res.data);
    } catch (error) {
      console.error('Representation error:', error);
    }
    setLoading(false);
  };

  const handleActiveLearning = async () => {
    setLoading(true);
    try {
      const res = await activeLearning(property);
      setActiveSuggestion(res.data.next_experiment);
    } catch (error) {
      console.error('Active learning error:', error);
    }
    setLoading(false);
  };

  return (
    <div style={{ 
      maxWidth: 800, 
      margin: '40px auto', 
      background: '#fff', 
      borderRadius: 16, 
      boxShadow: '0 8px 32px rgba(0,0,0,0.1)', 
      padding: 32,
      fontFamily: 'Arial, sans-serif'
    }}>
      <h1 style={{ textAlign: 'center', color: '#764ba2', marginBottom: 30 }}>
        🧪 Polymer Informatics Portal
      </h1>
      
      <div style={{ marginBottom: 24 }}>
        <div style={{ marginBottom: 16 }}>
          <label style={{ display: 'block', marginBottom: 8, fontWeight: 'bold' }}>
            SMILES String:
          </label>
          <input
            type="text"
            value={smiles}
            onChange={e => setSmiles(e.target.value)}
            placeholder="Enter SMILES (e.g., CCO)"
            style={{ 
              width: '100%', 
              padding: 12, 
              fontSize: 16, 
              border: '2px solid #ddd',
              borderRadius: 8,
              boxSizing: 'border-box'
            }}
          />
        </div>
        
        <div style={{ marginBottom: 16 }}>
          <label style={{ display: 'block', marginBottom: 8, fontWeight: 'bold' }}>
            Property:
          </label>
          <input
            type="text"
            value={property}
            onChange={e => setProperty(e.target.value)}
            placeholder="Property (e.g., Tg, molecular_weight)"
            style={{ 
              width: '100%', 
              padding: 12, 
              fontSize: 16, 
              border: '2px solid #ddd',
              borderRadius: 8,
              boxSizing: 'border-box'
            }}
          />
        </div>
      </div>
      
      <div style={{ display: 'flex', gap: 12, marginBottom: 24, flexWrap: 'wrap' }}>
        <button onClick={handlePredict} disabled={loading} style={buttonStyle}>
          Predict Property
        </button>
        <button onClick={handleVisualize} disabled={loading} style={buttonStyle}>
          Visualize Structure
        </button>
        <button onClick={handleDescriptors} disabled={loading} style={buttonStyle}>
          Get Descriptors
        </button>
        <button onClick={handleRepresent} disabled={loading} style={buttonStyle}>
          Represent Polymer
        </button>
        <button onClick={handleActiveLearning} disabled={loading} style={buttonStyle}>
          Active Learning
        </button>
      </div>
      
      {loading && (
        <div style={{ textAlign: 'center', margin: '20px 0', color: '#666' }}>
          Loading...
        </div>
      )}
      
      {imageUrl && (
        <div style={{ marginTop: 24, textAlign: 'center' }}>
          <h3>Molecular Structure</h3>
          <img src={imageUrl} alt="Molecule" style={{ maxWidth: 400, border: '1px solid #ddd', borderRadius: 8 }} />
        </div>
      )}
      
      {prediction && (
        <div style={{ marginTop: 24 }}>
          <h3>Prediction Results</h3>
          <pre style={{ background: '#f5f5f5', padding: 16, borderRadius: 8, overflow: 'auto' }}>
            {JSON.stringify(prediction, null, 2)}
          </pre>
        </div>
      )}
      
      {descriptors && (
        <div style={{ marginTop: 24 }}>
          <h3>Molecular Descriptors</h3>
          <pre style={{ background: '#f5f5f5', padding: 16, borderRadius: 8, overflow: 'auto' }}>
            {JSON.stringify(descriptors, null, 2)}
          </pre>
        </div>
      )}
      
      {representation && (
        <div style={{ marginTop: 24 }}>
          <h3>Polymer Representations</h3>
          <pre style={{ background: '#f5f5f5', padding: 16, borderRadius: 8, overflow: 'auto' }}>
            {JSON.stringify(representation, null, 2)}
          </pre>
        </div>
      )}
      
      {activeSuggestion && (
        <div style={{ marginTop: 24 }}>
          <h3>Next Experiment Suggestion</h3>
          <div style={{ background: '#e8f4f8', padding: 16, borderRadius: 8 }}>
            <strong>{activeSuggestion}</strong>
          </div>
        </div>
      )}
    </div>
  );
}

const buttonStyle = {
  padding: '12px 20px',
  fontSize: 14,
  border: 'none',
  borderRadius: 8,
  background: 'linear-gradient(135deg, #667eea 0%, #764ba2 100%)',
  color: 'white',
  cursor: 'pointer',
  fontWeight: 'bold',
  transition: 'transform 0.2s ease',
  flex: '1 1 auto',
  minWidth: '120px'
};

export default App;
