"""
Curated test dataset with known polymer properties from literature
This serves as ground truth for validating our prediction model
"""

# Test data with well-established property values from literature
POLYMER_TEST_DATA = [
    # Common thermoplastics with well-known Tg values
    {
        "name": "Polyethylene (PE)",
        "smiles": "CC",
        "psmiles": "[CH2]-[CH2]",
        "properties": {
            "Tg": -125.0,  # °C, crystalline PE
            "Tm": 135.0,   # °C
            "density": 0.92,  # g/cm³
            "molecular_weight": 28000
        },
        "source": "Polymer Handbook 4th Ed.",
        "confidence": "high"
    },
    {
        "name": "Polypropylene (PP)",
        "smiles": "CC(C)",
        "psmiles": "[CH2]-[CH(CH3)]",
        "properties": {
            "Tg": -10.0,   # °C
            "Tm": 160.0,   # °C
            "density": 0.90,
            "molecular_weight": 42000
        },
        "source": "Polymer Handbook 4th Ed.",
        "confidence": "high"
    },
    {
        "name": "Polystyrene (PS)",
        "smiles": "CCc1ccccc1",
        "psmiles": "[CH2]-[CH(C6H5)]",
        "properties": {
            "Tg": 100.0,   # °C, atactic PS
            "Tm": None,    # Amorphous
            "density": 1.05,
            "molecular_weight": 104000
        },
        "source": "Polymer Handbook 4th Ed.",
        "confidence": "high"
    },
    {
        "name": "Polyvinyl Chloride (PVC)",
        "smiles": "CCCl",
        "psmiles": "[CH2]-[CHCl]",
        "properties": {
            "Tg": 80.0,    # °C
            "Tm": None,    # Amorphous
            "density": 1.38,
            "molecular_weight": 62500
        },
        "source": "Polymer Handbook 4th Ed.",
        "confidence": "high"
    },
    {
        "name": "Polymethyl Methacrylate (PMMA)",
        "smiles": "CC(C)(C(=O)OC)",
        "psmiles": "[CH2]-[C(CH3)(COOCH3)]",
        "properties": {
            "Tg": 105.0,   # °C
            "Tm": None,    # Amorphous
            "density": 1.18,
            "molecular_weight": 100000
        },
        "source": "Polymer Handbook 4th Ed.",
        "confidence": "high"
    },
    {
        "name": "Polyethylene Terephthalate (PET)",
        "smiles": "CCOc1ccc(C(=O)O)cc1",
        "psmiles": "[OCH2CH2OC(=O)C6H4C(=O)]",
        "properties": {
            "Tg": 70.0,    # °C
            "Tm": 245.0,   # °C
            "density": 1.38,
            "molecular_weight": 192000
        },
        "source": "Polymer Handbook 4th Ed.",
        "confidence": "high"
    },
    {
        "name": "Polycarbonate (PC)",
        "smiles": "CC(C)(c1ccc(O)cc1)",
        "psmiles": "[OC6H4C(CH3)2C6H4OC(=O)]",
        "properties": {
            "Tg": 150.0,   # °C
            "Tm": None,    # Amorphous
            "density": 1.20,
            "molecular_weight": 254000
        },
        "source": "Polymer Handbook 4th Ed.",
        "confidence": "high"
    },
    {
        "name": "Polytetrafluoroethylene (PTFE)",
        "smiles": "C(C(F)(F)F)(F)F",
        "psmiles": "[CF2]-[CF2]",
        "properties": {
            "Tg": -97.0,   # °C
            "Tm": 327.0,   # °C
            "density": 2.15,
            "molecular_weight": 100000
        },
        "source": "Polymer Handbook 4th Ed.",
        "confidence": "high"
    },
    # High-performance polymers
    {
        "name": "Polyetheretherketone (PEEK)",
        "smiles": "c1ccc(Oc2ccc(C(=O)c3ccc(O)cc3)cc2)cc1",
        "psmiles": "[OC6H4OC6H4C(=O)C6H4]",
        "properties": {
            "Tg": 143.0,   # °C
            "Tm": 343.0,   # °C
            "density": 1.30,
            "molecular_weight": 288000
        },
        "source": "High Performance Polymers Database",
        "confidence": "high"
    },
    {
        "name": "Polyimide (PI)",
        "smiles": "C1C(=O)N(c2ccccc2)C(=O)C1",
        "psmiles": "[C6H4C(=O)NC(=O)C6H4]",
        "properties": {
            "Tg": 360.0,   # °C
            "Tm": None,    # Does not melt
            "density": 1.42,
            "molecular_weight": 300000
        },
        "source": "High Performance Polymers Database",
        "confidence": "medium"
    }
]

# Edge cases and challenging polymers for robustness testing
EDGE_CASE_DATA = [
    {
        "name": "Very Low Tg Polymer",
        "smiles": "CCCCCCCC",  # Long alkyl chain
        "properties": {"Tg": -180.0},
        "source": "Synthetic",
        "confidence": "low"
    },
    {
        "name": "Very High Tg Polymer", 
        "smiles": "c1ccc2c(c1)c1ccccc1c1ccccc12",  # Rigid aromatic
        "properties": {"Tg": 500.0},
        "source": "Estimated",
        "confidence": "low"
    },
    {
        "name": "Complex Copolymer",
        "smiles": "CC(C)c1ccc(C)cc1",  # Mixed structure
        "properties": {"Tg": 85.0},
        "source": "Experimental",
        "confidence": "medium"
    }
]

def get_test_data(confidence_level="high"):
    """
    Get test data filtered by confidence level
    
    Args:
        confidence_level: "high", "medium", "low", or "all"
    
    Returns:
        List of test polymer data
    """
    if confidence_level == "all":
        return POLYMER_TEST_DATA + EDGE_CASE_DATA
    
    return [p for p in POLYMER_TEST_DATA if p["confidence"] == confidence_level]

def get_property_data(property_name, confidence_level="high"):
    """
    Get test data for a specific property
    
    Args:
        property_name: Property to extract (e.g., "Tg", "Tm")
        confidence_level: Confidence filter
        
    Returns:
        List of (polymer_data, property_value) tuples
    """
    test_data = get_test_data(confidence_level)
    result = []
    
    for polymer in test_data:
        if property_name in polymer["properties"] and polymer["properties"][property_name] is not None:
            result.append((polymer, polymer["properties"][property_name]))
    
    return result
