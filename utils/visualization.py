# Molecule rendering

from rdkit.Chem import Draw

def render_molecule(mol):
    """
    Render molecule as PNG image using RDKit.
    Returns a PIL Image object.
    """
    return Draw.MolToImage(mol, size=(400, 300))
