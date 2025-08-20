from rdkit import Chem

def to_psmiles(mol, stereo=True, branches=True):
    """
    Convert RDKit molecule to p-SMILES with stereo/branch support.
    - Assigns stereochemistry
    - Returns SMILES with explicit bonds and hydrogens
    """
    Chem.AssignStereochemistry(mol, cleanIt=True)
    return Chem.MolToSmiles(
        mol,
        isomericSmiles=stereo,
        allBondsExplicit=True,
        allHsExplicit=True
    )

# Example usage
polymer = Chem.MolFromSmiles("CCOC(=O)C(C)C")
psmiles = to_psmiles(polymer)
print(f"Enhanced p-SMILES: {psmiles}")
