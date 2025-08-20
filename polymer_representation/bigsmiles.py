# BigSMILES handling

from rdkit import Chem

def to_bigsmiles(mol):
    """
    Convert RDKit molecule to BigSMILES string.
    - Uses canonical SMILES as base
    - Adds curly brackets for repeat units (simple heuristic)
    """
    smiles = Chem.MolToSmiles(mol)
    # Heuristic: wrap the main chain in curly brackets
    atoms = [atom.GetSymbol() for atom in mol.GetAtoms()]
    if len(atoms) > 4:
        return '{' + smiles + '}'
    return smiles
