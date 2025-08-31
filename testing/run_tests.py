#!/usr/bin/env python3
"""
Master test runner - executes all validation tests
"""
import os
import sys
import subprocess
from datetime import datetime

def run_test_suite():
    """Run all available tests in sequence"""
    print("🧪 POLYMER INFORMATICS MODEL VALIDATION SUITE")
    print("=" * 60)
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()
    
    # List of test scripts to run
    test_scripts = [
        ("Data Quality Check", "validation/data_quality_check.py"),
        ("Quick Accuracy Test", "validation/quick_test.py"),
        ("Full Accuracy Test", "validation/test_accuracy.py"),
        ("Comprehensive Report", "validation/run_all_tests.py")
    ]
    
    results = {}
    
    for test_name, script_path in test_scripts:
        print(f"🔬 Running: {test_name}")
        print("-" * 40)
        
        try:
            # Run the test script
            result = subprocess.run([sys.executable, script_path], 
                                  capture_output=True, text=True, timeout=300)
            
            if result.returncode == 0:
                print(result.stdout)
                results[test_name] = "✅ PASSED"
            else:
                print(f"❌ FAILED: {result.stderr}")
                results[test_name] = "❌ FAILED"
                
        except subprocess.TimeoutExpired:
            print("⏰ TIMEOUT: Test took too long")
            results[test_name] = "⏰ TIMEOUT"
        except Exception as e:
            print(f"💥 ERROR: {e}")
            results[test_name] = "💥 ERROR"
        
        print()
    
    # Summary
    print("📊 TEST SUITE SUMMARY")
    print("=" * 60)
    for test_name, status in results.items():
        print(f"{test_name:<30} {status}")
    
    # Overall result
    failed_tests = [name for name, status in results.items() if "FAILED" in status or "ERROR" in status]
    
    if not failed_tests:
        print("\n🎉 ALL TESTS PASSED! Model validation successful.")
    else:
        print(f"\n⚠️  {len(failed_tests)} test(s) failed: {', '.join(failed_tests)}")
    
    print(f"\nCompleted: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

if __name__ == "__main__":
    # Change to testing directory
    testing_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(testing_dir)
    
    run_test_suite()
