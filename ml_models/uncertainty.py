# Uncertainty quantification

import numpy as np

def estimate_uncertainty(polymer_id, property_name):
    """
    Estimate prediction uncertainty using variance of similar polymers.
    """
    # This would normally use model prediction variance
    # For now, return a random value for demonstration
    return float(np.random.uniform(0.05, 0.15))
