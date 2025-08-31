from rdkit import Chem
from rdkit.Chem import AllChem
import torch
from transformers import AutoModel, AutoTokenizer
from rdkit import DataStructs

class HybridSearch:
    def __init__(self, llm_model="bert-base-uncased"):
        self.tokenizer = AutoTokenizer.from_pretrained(llm_model)
        self.model = AutoModel.from_pretrained(llm_model)
    
    def structural_similarity(self, query_psmiles, targets):
        """Calculate structural similarity using fingerprints"""
        query_mol = Chem.MolFromSmiles(query_psmiles)
        if query_mol is None:
            # If query molecule is invalid, return zero similarity for all
            return [(target, 0.0) for target in targets]
        
        query_fp = AllChem.GetMorganFingerprint(query_mol, radius=3)
        results = []
        for target in targets:
            try:
                target_smiles = target.get('psmiles', target.get('smiles', ''))
                target_mol = Chem.MolFromSmiles(target_smiles)
                if target_mol is None:
                    # Skip invalid molecules
                    results.append((target, 0.0))
                    continue
                
                target_fp = AllChem.GetMorganFingerprint(target_mol, radius=3)
                similarity = DataStructs.TanimotoSimilarity(query_fp, target_fp)
                results.append((target, similarity))
            except Exception as e:
                # Handle any RDKit errors gracefully
                results.append((target, 0.0))
        return sorted(results, key=lambda x: x[1], reverse=True)
    
    def semantic_similarity(self, query_text, targets):
        """Calculate semantic similarity using LLM embeddings"""
        try:
            inputs = self.tokenizer(query_text, return_tensors="pt", truncation=True, padding=True)
            query_embedding = self.model(**inputs).last_hidden_state.mean(dim=1)
            
            results = []
            for target in targets:
                try:
                    target_text = target.get('description', target.get('name', 'polymer'))
                    target_inputs = self.tokenizer(target_text, return_tensors="pt", truncation=True, padding=True)
                    target_embedding = self.model(**target_inputs).last_hidden_state.mean(dim=1)
                    similarity = torch.cosine_similarity(query_embedding, target_embedding)
                    results.append((target, similarity.item()))
                except Exception as e:
                    # Handle any tokenization/model errors gracefully
                    results.append((target, 0.0))
            return sorted(results, key=lambda x: x[1], reverse=True)
        except Exception as e:
            # If semantic similarity fails completely, return zero similarity for all
            return [(target, 0.0) for target in targets]
    
    def hybrid_search(self, query_psmiles, query_text, targets, alpha=0.7):
        """Combine structural and semantic similarity"""
        struct_results = self.structural_similarity(query_psmiles, targets)
        sem_results = self.semantic_similarity(query_text, targets)
        
        # Combine scores
        combined = []
        for i, (struct, sem) in enumerate(zip(struct_results, sem_results)):
            combined_score = alpha * struct[1] + (1-alpha) * sem[1]
            combined.append((struct[0], combined_score))
        return sorted(combined, key=lambda x: x[1], reverse=True)
