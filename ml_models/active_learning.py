# Active learning components

def suggest_next_experiment(property_name, db):
    """
    Suggest next experiment using active learning.
    - Finds polymers with highest uncertainty or least data
    """
    polymers = db.get_all_polymers()
    # Example: pick polymer with least properties for the target property
    candidates = [p for p in polymers if property_name not in p]
    if candidates:
        return candidates[0]['psmiles']
    return None
