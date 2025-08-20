# Processes user queries

def parse_query(query):
    """
    Parse user query for polymer search.
    Returns a dict with structure, property, and text fields.
    """
    result = {'structure': None, 'property': None, 'text': query}
    # Simple heuristics: look for SMILES, property keywords
    if 'smiles:' in query.lower():
        result['structure'] = query.split('smiles:')[1].strip()
    if 'property:' in query.lower():
        result['property'] = query.split('property:')[1].strip()
    return result
