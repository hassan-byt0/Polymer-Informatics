#!/usr/bin/env python3
"""
Accuracy Testing Framework for Polymer Property Predictions
Tests the model against known polymer data and calculates accuracy metrics
"""
import os
import sys
import json
import requests
import numpy as np
from datetime import datetime

# Add parent directory to path for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from test_data.known_polymers import get_property_data, POLYMER_TEST_DATA

class PredictionTester:
    def __init__(self, api_base_url="http://localhost:5050"):
        self.api_url = api_base_url
        self.results = []
        
    def test_property_predictions(self, property_name="Tg", confidence_level="high"):
        """
        Test predictions against known polymer data
        
        Args:
            property_name: Property to test (e.g., "Tg", "Tm")
            confidence_level: Quality filter for test data
            
        Returns:
            Dictionary with test results and metrics
        """
        print(f"🧪 Testing {property_name} predictions...")
        
        test_data = get_property_data(property_name, confidence_level)
        if not test_data:
            print(f"❌ No test data found for {property_name}")
            return None
            
        print(f"📊 Testing against {len(test_data)} known polymers")
        
        predictions = []
        actual_values = []
        errors = []
        
        for i, (polymer, actual_value) in enumerate(test_data):
            try:
                # Make prediction request
                response = self.predict_property(polymer["smiles"], property_name)
                
                if response and "prediction" in response and response["prediction"] is not None:
                    predicted_value = response["prediction"]
                    error = abs(predicted_value - actual_value)
                    
                    result = {
                        "polymer_name": polymer["name"],
                        "smiles": polymer["smiles"],
                        "actual_value": actual_value,
                        "predicted_value": predicted_value,
                        "error": error,
                        "confidence": response.get("confidence", 0),
                        "uncertainty": response.get("uncertainty", 0),
                        "source": polymer["source"]
                    }
                    
                    predictions.append(predicted_value)
                    actual_values.append(actual_value)
                    errors.append(error)
                    self.results.append(result)
                    
                    print(f"  {i+1:2d}. {polymer['name']:<25} | "
                          f"Actual: {actual_value:7.1f}°C | "
                          f"Predicted: {predicted_value:7.1f}°C | "
                          f"Error: {error:6.1f}°C")
                else:
                    print(f"  {i+1:2d}. {polymer['name']:<25} | ❌ Prediction failed")
                    
            except Exception as e:
                print(f"  {i+1:2d}. {polymer['name']:<25} | ❌ Error: {e}")
                continue
        
        if not predictions:
            print("❌ No successful predictions made")
            return None
            
        # Calculate metrics
        metrics = self.calculate_metrics(actual_values, predictions, errors)
        
        # Print summary
        self.print_summary(property_name, metrics)
        
        return {
            "property": property_name,
            "metrics": metrics,
            "results": self.results[-len(predictions):],  # Last batch of results
            "test_count": len(predictions),
            "timestamp": datetime.now().isoformat()
        }
    
    def predict_property(self, smiles, property_name):
        """Make a prediction request to the API"""
        try:
            response = requests.post(
                f"{self.api_url}/predict",
                json={"smiles": smiles, "property": property_name},
                timeout=30
            )
            if response.status_code == 200:
                return response.json()
            else:
                print(f"API Error: {response.status_code}")
                return None
        except Exception as e:
            print(f"Request Error: {e}")
            return None
    
    def calculate_metrics(self, actual, predicted, errors):
        """Calculate prediction accuracy metrics"""
        actual = np.array(actual)
        predicted = np.array(predicted)
        errors = np.array(errors)
        
        # Basic metrics
        mae = np.mean(errors)  # Mean Absolute Error
        rmse = np.sqrt(np.mean((actual - predicted) ** 2))  # Root Mean Square Error
        
        # R-squared
        ss_res = np.sum((actual - predicted) ** 2)
        ss_tot = np.sum((actual - np.mean(actual)) ** 2)
        r2 = 1 - (ss_res / ss_tot) if ss_tot != 0 else 0
        
        # Accuracy within tolerance
        tolerance_5 = np.mean(errors <= 5.0) * 100  # Within ±5°C
        tolerance_10 = np.mean(errors <= 10.0) * 100  # Within ±10°C
        tolerance_20 = np.mean(errors <= 20.0) * 100  # Within ±20°C
        
        # Error statistics
        max_error = np.max(errors)
        min_error = np.min(errors)
        median_error = np.median(errors)
        std_error = np.std(errors)
        
        return {
            "mae": float(mae),
            "rmse": float(rmse),
            "r2": float(r2),
            "max_error": float(max_error),
            "min_error": float(min_error),
            "median_error": float(median_error),
            "std_error": float(std_error),
            "accuracy_5c": float(tolerance_5),
            "accuracy_10c": float(tolerance_10),
            "accuracy_20c": float(tolerance_20)
        }
    
    def print_summary(self, property_name, metrics):
        """Print formatted test summary"""
        print(f"\n📈 {property_name} Prediction Accuracy Summary")
        print("=" * 50)
        print(f"Mean Absolute Error (MAE):    {metrics['mae']:8.2f}°C")
        print(f"Root Mean Square Error:       {metrics['rmse']:8.2f}°C")
        print(f"R² Score:                     {metrics['r2']:8.3f}")
        print(f"Max Error:                    {metrics['max_error']:8.2f}°C")
        print(f"Median Error:                 {metrics['median_error']:8.2f}°C")
        print(f"Error Std Deviation:          {metrics['std_error']:8.2f}°C")
        print()
        print("Accuracy within tolerance:")
        print(f"  ±5°C:   {metrics['accuracy_5c']:6.1f}% of predictions")
        print(f"  ±10°C:  {metrics['accuracy_10c']:6.1f}% of predictions")
        print(f"  ±20°C:  {metrics['accuracy_20c']:6.1f}% of predictions")
        print()
        
        # Quality assessment
        if metrics['mae'] < 10:
            quality = "🟢 Excellent"
        elif metrics['mae'] < 25:
            quality = "🟡 Good"
        elif metrics['mae'] < 50:
            quality = "🟠 Fair"
        else:
            quality = "🔴 Poor"
            
        print(f"Overall Model Quality: {quality}")
        print("=" * 50)

def main():
    """Run the prediction accuracy tests"""
    print("🧪 Polymer Property Prediction Accuracy Test")
    print("=" * 60)
    
    tester = PredictionTester()
    
    # Test Tg predictions
    results = tester.test_property_predictions("Tg", confidence_level="high")
    
    if results:
        # Save results
        os.makedirs("reports", exist_ok=True)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        report_file = f"reports/accuracy_test_{timestamp}.json"
        
        with open(report_file, 'w') as f:
            json.dump(results, f, indent=2)
            
        print(f"\n💾 Results saved to: {report_file}")
    
    return results

if __name__ == "__main__":
    main()
