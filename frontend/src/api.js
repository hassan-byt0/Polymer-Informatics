import axios from 'axios';

const API_BASE = 'http://localhost:5050';

export const ingestData = async (sources) => {
  return axios.post(`${API_BASE}/ingest`, { sources });
};

export const predictProperty = async (smiles, property) => {
  return axios.post(`${API_BASE}/predict`, { smiles, property });
};

export const visualizeStructure = async (smiles) => {
  return axios.post(`${API_BASE}/visualize`, { smiles }, { responseType: 'blob' });
};

export const getDescriptors = async (smiles) => {
  return axios.post(`${API_BASE}/descriptors`, { smiles });
};

export const representPolymer = async (smiles) => {
  return axios.post(`${API_BASE}/represent`, { smiles });
};

export const activeLearning = async (property) => {
  return axios.post(`${API_BASE}/active_learning`, { property });
};
