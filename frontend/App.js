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
    const res = await predictProperty(smiles, property);
    setPrediction(res.data);
    setLoading(false);
  };

  const handleVisualize = async () => {
    setLoading(true);
    const res = await visualizeStructure(smiles);
    const url = URL.createObjectURL(res.data);
    setImageUrl(url);
    setLoading(false);
  };

  const handleDescriptors = async () => {
    setLoading(true);
    const res = await getDescriptors(smiles);
    setDescriptors(res.data);
    setLoading(false);
  };

  const handleRepresent = async () => {
    setLoading(true);
    const res = await representPolymer(smiles);
    setRepresentation(res.data);
    setLoading(false);
  };

  const handleActiveLearning = async () => {
    setLoading(true);
    const res = await activeLearning(property);
    setActiveSuggestion(res.data.next_experiment);
    setLoading(false);
  };

  return (
    <div style={{ maxWidth: 700, margin: '40px auto', background: '#fff', borderRadius: 16, boxShadow: '0 8px 32px #0002', padding: 32 }}>
      <h1 style={{ textAlign: 'center', color: '#764ba2' }}>Polymer Informatics Portal</h1>
      <div style={{ marginBottom: 24 }}>
        <input
          type="text"
          value={smiles}
          onChange={e => setSmiles(e.target.value)}
          placeholder="Enter SMILES"
          style={{ width: '60%', padding: 8, fontSize: 16, marginRight: 8 }}
        />
        <input
          type="text"
          value={property}
          onChange={e => setProperty(e.target.value)}
          placeholder="Property (e.g. Tg)"
          style={{ width: '30%', padding: 8, fontSize: 16 }}
        />
      </div>
      <div style={{ display: 'flex', gap: 12, marginBottom: 24 }}>
        <button onClick={handlePredict} disabled={loading}>Predict</button>
        <button onClick={handleVisualize} disabled={loading}>Visualize</button>
        <button onClick={handleDescriptors} disabled={loading}>Descriptors</button>
        <button onClick={handleRepresent} disabled={loading}>Represent</button>
        <button onClick={handleActiveLearning} disabled={loading}>Active Learning</button>
      </div>
      {imageUrl && <img src={imageUrl} alt="Molecule" style={{ maxWidth: 400, margin: '16px auto', display: 'block' }} />}
      {prediction && (
        <div style={{ marginTop: 24 }}>
          <h3>Prediction</h3>
          <pre>{JSON.stringify(prediction, null, 2)}</pre>
        </div>
      )}
      {descriptors && (
        <div style={{ marginTop: 24 }}>
          <h3>Descriptors</h3>
          <pre>{JSON.stringify(descriptors, null, 2)}</pre>
        </div>
      )}
      {representation && (
        <div style={{ marginTop: 24 }}>
          <h3>Polymer Representation</h3>
          <pre>{JSON.stringify(representation, null, 2)}</pre>
        </div>
      )}
      {activeSuggestion && (
        <div style={{ marginTop: 24 }}>
          <h3>Next Experiment Suggestion</h3>
          <pre>{activeSuggestion}</pre>
        </div>
      )}
    </div>
  );
}

export default App;
