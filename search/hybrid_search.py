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
        query_fp = AllChem.GetMorganFingerprint(
            Chem.MolFromSmiles(query_psmiles), 
            radius=3
        )
        results = []
        for target in targets:
            target_fp = AllChem.GetMorganFingerprint(
                Chem.MolFromSmiles(target['psmiles']), 
                radius=3
            )
            similarity = DataStructs.TanimotoSimilarity(query_fp, target_fp)
            results.append((target, similarity))
        return sorted(results, key=lambda x: x[1], reverse=True)
    
    def semantic_similarity(self, query_text, targets):
        """Calculate semantic similarity using LLM embeddings"""
        inputs = self.tokenizer(query_text, return_tensors="pt")
        query_embedding = self.model(**inputs).last_hidden_state.mean(dim=1)
        
        results = []
        for target in targets:
            target_embedding = self.model(
                **self.tokenizer(target['description'], return_tensors="pt")
            ).last_hidden_state.mean(dim=1)
            similarity = torch.cosine_similarity(query_embedding, target_embedding)
            results.append((target, similarity.item()))
        return sorted(results, key=lambda x: x[1], reverse=True)
    
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
