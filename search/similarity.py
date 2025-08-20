from rdkit import DataStructs
from rdkit.Chem import AllChem, Chem

def tanimoto_similarity(smiles1, smiles2):
    """
    Compute Tanimoto similarity between two SMILES using RDKit fingerprints.
    """
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)
    fp1 = AllChem.GetMorganFingerprint(mol1, radius=2)
    fp2 = AllChem.GetMorganFingerprint(mol2, radius=2)
    return DataStructs.TanimotoSimilarity(fp1, fp2)
# Similarity metrics
