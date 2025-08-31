#!/usr/bin/env python3
"""
Quick accuracy test script to run a fast validation
"""
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from validation.test_accuracy import PredictionTester

def quick_test():
    """Run a quick accuracy test with a few polymers"""
    print("🚀 Quick Accuracy Test")
    print("=" * 40)
    
    tester = PredictionTester()
    
    # Test with just a few well-known polymers
    quick_tests = [
        {"smiles": "CC", "name": "Ethane", "expected_tg": -125, "property": "Tg"},
        {"smiles": "CCc1ccccc1", "name": "Styrene", "expected_tg": 100, "property": "Tg"},
        {"smiles": "CC(C)", "name": "Propylene", "expected_tg": -10, "property": "Tg"}
    ]
    
    results = []
    print(f"Testing {len(quick_tests)} polymers...")
    print()
    
    for test in quick_tests:
        try:
            response = tester.predict_property(test["smiles"], test["property"])
            
            if response and "prediction" in response:
                predicted = response["prediction"]
                actual = test["expected_tg"]
                error = abs(predicted - actual)
                
                print(f"✓ {test['name']:<15} | "
                      f"Expected: {actual:6.1f}°C | "
                      f"Predicted: {predicted:6.1f}°C | "
                      f"Error: {error:6.1f}°C")
                
                results.append({
                    "polymer": test["name"],
                    "actual": actual,
                    "predicted": predicted,
                    "error": error
                })
            else:
                print(f"✗ {test['name']:<15} | Prediction failed")
                
        except Exception as e:
            print(f"✗ {test['name']:<15} | Error: {e}")
    
    if results:
        avg_error = sum(r["error"] for r in results) / len(results)
        max_error = max(r["error"] for r in results)
        print()
        print(f"Average Error: {avg_error:.1f}°C")
        print(f"Maximum Error: {max_error:.1f}°C")
        
        if avg_error < 20:
            print("🟢 Quick test PASSED - Model appears to be working")
        else:
            print("🔴 Quick test FAILED - Model needs attention")
    else:
        print("🔴 No successful predictions made")

if __name__ == "__main__":
    quick_test()
