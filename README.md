# Polymer Informatics Pipeline

This repository contains a full-stack, Dockerized pipeline for advanced polymer informatics, including data extraction from scientific PDFs, property prediction, molecular visualization, and a modern web interface.

## Features
- **Data Ingestion:** Extracts polymer data and properties from scientific PDFs using NLP and scrapers.
- **Normalization:** Cleans and standardizes polymer records.
- **Graph Database:** Stores polymers, properties, and literature in Neo4j.
- **Polymer Representation:** Converts structures to pSMILES, BigSMILES, and graph encodings.
- **Hybrid Search:** Combines structural and semantic similarity using RDKit and BERT.
- **ML Models:** Predicts properties, estimates uncertainty, and suggests experiments.
- **Web Frontend:** React-based UI for input, visualization, and results.
- **Dockerized:** All components run in containers for easy setup and reproducibility.

## Quick Start

### 1. Prerequisites
- [Docker Desktop](https://www.docker.com/products/docker-desktop) (Mac/Windows/Linux)
- [Git](https://git-scm.com/)

### 2. Clone the Repository
```zsh
git clone git@github.com:hassan-byt0/PolyInfo_v2.git
cd PolyInfo_v2/polymer-informatics
```

### 3. Configure Environment
- Place your scientific PDF (e.g. `vdoc.pub_polymers-a-property-database.pdf`) in the project root.
- Edit `.env` with your PDF filename and Neo4j credentials:
  ```
  PDF_PATH=./vdoc.pub_polymers-a-property-database.pdf
  NEO4J_URI=bolt://neo4j-polymer:7687
  NEO4J_USER=neo4j
  NEO4J_PASSWORD=polymer123
  ```

### 4. Build and Run the Pipeline
```zsh
docker-compose up --build
```

### 5. Access the Services
- **Frontend:** [http://localhost:3000](http://localhost:3000)
- **Backend API:** [http://localhost:5050](http://localhost:5050)
- **Neo4j Browser:** [http://localhost:7474](http://localhost:7474) (user: neo4j, password: polymer123)

## Usage
- Use the web UI to input SMILES, select properties, and visualize results.
- The backend will extract, normalize, and store data from your PDF, and serve predictions and visualizations.

## Development
- All source code is in subfolders (`data_ingestion`, `ml_models`, `polymer_representation`, etc.).
- Frontend code is in `frontend/src/`.
- Backend API is in `user_interface/app.py`.

## Troubleshooting
- If you see a white screen, check browser dev tools for errors.
- If containers fail, check logs with:
  ```zsh
  docker logs polymer-backend
  docker logs polymer-frontend
  docker logs neo4j-polymer
  ```
- Make sure your PDF filename matches `.env` and Dockerfile.

## License
MIT

## Contact
For questions or collaboration, open an issue or contact [hassan-byt0](https://github.com/hassan-byt0).
