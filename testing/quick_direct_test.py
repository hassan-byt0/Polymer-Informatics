#!/usr/bin/env python3
"""
Quick direct model test - bypasses API for speed
Tests the ML models directly with known polymer data
"""
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from datetime import datetime

# Import our test data
from testing.test_data.known_polymers import POLYMER_TEST_DATA

# Import model components directly
try:
    from ml_models.property_prediction import PropertyPredictor
    from polymer_representation.psmiles import PSMILESParser
    from utils.cheminformatics import ChemUtils
except ImportError as e:
    print(f"❌ Import error: {e}")
    sys.exit(1)

def quick_accuracy_test():
    """Run a fast accuracy test on key polymers"""
    print("🚀 QUICK POLYMER MODEL ACCURACY TEST")
    print("=" * 50)
    print(f"Time: {datetime.now().strftime('%H:%M:%S')}")
    print()
    
    # Initialize components
    try:
        predictor = PropertyPredictor()
        print("✅ PropertyPredictor initialized")
    except Exception as e:
        print(f"❌ Failed to initialize PropertyPredictor: {e}")
        return
    
    # Test on first 5 polymers for speed
    test_polymers = POLYMER_TEST_DATA[:5]
    
    results = {
        "Tg": {"predictions": [], "actual": [], "errors": []},
        "Tm": {"predictions": [], "actual": [], "errors": []},
    }
    
    print("\n🧪 Testing predictions...")
    print("-" * 30)
    
    for i, polymer in enumerate(test_polymers):
        print(f"\n{i+1}. {polymer['name']}")
        smiles = polymer['smiles']
        
        # Test Tg prediction
        if 'Tg' in polymer['properties'] and polymer['properties']['Tg'] is not None:
            try:
                pred_tg = predictor.predict_property(smiles, 'Tg')
                actual_tg = polymer['properties']['Tg']
                error = abs(pred_tg - actual_tg)
                
                results["Tg"]["predictions"].append(pred_tg)
                results["Tg"]["actual"].append(actual_tg)
                results["Tg"]["errors"].append(error)
                
                print(f"   Tg: Predicted {pred_tg:.1f}°C, Actual {actual_tg:.1f}°C, Error {error:.1f}°C")
                
                if error > 100:
                    print(f"   ⚠️  Large error detected!")
                elif error < 20:
                    print(f"   ✅ Good prediction!")
                    
            except Exception as e:
                print(f"   ❌ Tg prediction failed: {e}")
        
        # Test Tm prediction  
        if 'Tm' in polymer['properties'] and polymer['properties']['Tm'] is not None:
            try:
                pred_tm = predictor.predict_property(smiles, 'Tm')
                actual_tm = polymer['properties']['Tm']
                error = abs(pred_tm - actual_tm)
                
                results["Tm"]["predictions"].append(pred_tm)
                results["Tm"]["actual"].append(actual_tm)
                results["Tm"]["errors"].append(error)
                
                print(f"   Tm: Predicted {pred_tm:.1f}°C, Actual {actual_tm:.1f}°C, Error {error:.1f}°C")
                
            except Exception as e:
                print(f"   ❌ Tm prediction failed: {e}")
    
    # Calculate metrics
    print("\n📊 ACCURACY SUMMARY")
    print("=" * 50)
    
    for prop, data in results.items():
        if len(data["predictions"]) > 0:
            mae = np.mean(data["errors"])
            rmse = np.sqrt(np.mean([e**2 for e in data["errors"]]))
            
            # R² calculation
            actual = np.array(data["actual"])
            predicted = np.array(data["predictions"])
            ss_res = np.sum((actual - predicted) ** 2)
            ss_tot = np.sum((actual - np.mean(actual)) ** 2)
            r2 = 1 - (ss_res / ss_tot) if ss_tot != 0 else 0
            
            print(f"\n{prop} Property:")
            print(f"  Samples tested: {len(data['predictions'])}")
            print(f"  Mean Absolute Error: {mae:.1f}°C")
            print(f"  Root Mean Square Error: {rmse:.1f}°C")
            print(f"  R² Score: {r2:.3f}")
            
            # Accuracy within tolerance
            within_20 = sum(1 for e in data["errors"] if e <= 20)
            within_50 = sum(1 for e in data["errors"] if e <= 50)
            
            print(f"  Accuracy within ±20°C: {within_20}/{len(data['errors'])} ({100*within_20/len(data['errors']):.1f}%)")
            print(f"  Accuracy within ±50°C: {within_50}/{len(data['errors'])} ({100*within_50/len(data['errors']):.1f}%)")
            
            # Quality assessment
            if mae < 30:
                print(f"  🎯 {prop} predictions: GOOD quality")
            elif mae < 100:
                print(f"  ⚠️  {prop} predictions: MODERATE quality")
            else:
                print(f"  ❌ {prop} predictions: POOR quality")
        else:
            print(f"\n{prop} Property: No successful predictions")
    
    print(f"\n⏱️  Test completed in < 10 seconds")
    print("📝 For detailed analysis, run the full test suite when time permits")

if __name__ == "__main__":
    quick_accuracy_test()
