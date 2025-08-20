# RDKit helpers

from rdkit.Chem import Descriptors

def compute_descriptors(mol):
    """
    Compute molecular descriptors using RDKit.
    Returns a dictionary of key descriptors.
    """
    return {
        'MolWt': Descriptors.MolWt(mol),
        'LogP': Descriptors.MolLogP(mol),
        'NumAtoms': mol.GetNumAtoms(),
        'NumBonds': mol.GetNumBonds(),
        'TPSA': Descriptors.TPSA(mol),
        'NumRotatableBonds': Descriptors.NumRotatableBonds(mol)
    }
