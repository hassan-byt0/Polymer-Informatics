from rdkit import Chem

def to_graph(mol):
    """
    Convert RDKit molecule to graph representation.
    Returns adjacency list and atom features.
    """
    atoms = [{'idx': atom.GetIdx(), 'symbol': atom.GetSymbol(), 'atomic_num': atom.GetAtomicNum()} for atom in mol.GetAtoms()]
    bonds = [{'begin': bond.GetBeginAtomIdx(), 'end': bond.GetEndAtomIdx(), 'type': str(bond.GetBondType())} for bond in mol.GetBonds()]
    adjacency = {atom['idx']: [] for atom in atoms}
    for bond in bonds:
        adjacency[bond['begin']].append(bond['end'])
        adjacency[bond['end']].append(bond['begin'])
    return {'atoms': atoms, 'bonds': bonds, 'adjacency': adjacency}
