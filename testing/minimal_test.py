#!/usr/bin/env python3
"""
Minimal component test - just checks if basic imports work
"""
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

print("🔧 MINIMAL COMPONENT TEST")
print("=" * 30)

# Test 1: Basic imports
try:
    from testing.test_data.known_polymers import POLYMER_TEST_DATA
    print("✅ Test data import: OK")
    print(f"   Found {len(POLYMER_TEST_DATA)} test polymers")
except Exception as e:
    print(f"❌ Test data import failed: {e}")

# Test 2: Check if model files exist
try:
    from ml_models.property_prediction import PropertyPredictor
    print("✅ PropertyPredictor import: OK")
except Exception as e:
    print(f"❌ PropertyPredictor import failed: {e}")

# Test 3: Check polymer representation
try:
    from polymer_representation.psmiles import PSMILESParser
    print("✅ PSMILESParser import: OK")
except Exception as e:
    print(f"❌ PSMILESParser import failed: {e}")

# Test 4: Check utils
try:
    from utils.cheminformatics import ChemUtils
    print("✅ ChemUtils import: OK")
except Exception as e:
    print(f"❌ ChemUtils import failed: {e}")

# Test 5: Quick data validation
print("\n🧪 TEST DATA SAMPLE:")
print("-" * 20)
for i, polymer in enumerate(POLYMER_TEST_DATA[:3]):
    print(f"{i+1}. {polymer['name']}")
    print(f"   SMILES: {polymer['smiles']}")
    print(f"   Tg: {polymer['properties'].get('Tg', 'N/A')}°C")
    print()

print("📊 QUICK ASSESSMENT:")
print("-" * 20)

# Check for known issues
tg_values = [p['properties'].get('Tg') for p in POLYMER_TEST_DATA if p['properties'].get('Tg') is not None]
print(f"Available Tg values: {len(tg_values)}")
print(f"Tg range: {min(tg_values):.1f}°C to {max(tg_values):.1f}°C")

# Identify problematic SMILES
print("\n🔍 POTENTIAL ISSUES:")
print("-" * 20)

try:
    from rdkit import Chem
    invalid_smiles = []
    for polymer in POLYMER_TEST_DATA[:5]:
        mol = Chem.MolFromSmiles(polymer['smiles'])
        if mol is None:
            invalid_smiles.append(polymer['name'])
    
    if invalid_smiles:
        print(f"❌ Invalid SMILES found: {invalid_smiles}")
    else:
        print("✅ All test SMILES are valid")
        
except ImportError:
    print("⚠️  RDKit not available for SMILES validation")

print("\n🎯 NEXT STEPS:")
print("1. If imports work, the slow issue is likely Neo4j connection")
print("2. Consider testing with mock data instead of database")
print("3. Check if extracted PDF data has quality issues")
