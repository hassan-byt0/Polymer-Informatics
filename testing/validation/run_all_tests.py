#!/usr/bin/env python3
"""
Comprehensive testing suite for polymer informatics system
Runs all validation tests and generates a complete accuracy report
"""
import os
import sys
import json
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

# Add parent directory to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from validation.test_accuracy import PredictionTester
from test_data.known_polymers import POLYMER_TEST_DATA

class TestSuite:
    def __init__(self):
        self.tester = PredictionTester()
        self.all_results = {}
        
    def run_all_tests(self):
        """Run comprehensive test suite"""
        print("🧪 Starting Comprehensive Polymer Informatics Test Suite")
        print("=" * 70)
        
        # Test different properties
        properties_to_test = ["Tg", "Tm", "density"]
        
        for prop in properties_to_test:
            print(f"\n🔬 Testing {prop} predictions...")
            results = self.tester.test_property_predictions(prop, confidence_level="high")
            if results:
                self.all_results[prop] = results
        
        # Generate comprehensive report
        self.generate_report()
        
    def generate_report(self):
        """Generate comprehensive accuracy report"""
        if not self.all_results:
            print("❌ No test results to report")
            return
            
        print("\n📊 Generating Comprehensive Accuracy Report...")
        
        # Create reports directory
        os.makedirs("reports", exist_ok=True)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # Generate text report
        report_text = self.create_text_report()
        text_file = f"reports/comprehensive_report_{timestamp}.txt"
        with open(text_file, 'w') as f:
            f.write(report_text)
        
        # Generate JSON report
        json_file = f"reports/comprehensive_report_{timestamp}.json"
        with open(json_file, 'w') as f:
            json.dump(self.all_results, f, indent=2)
        
        # Generate visualizations
        self.create_visualizations(timestamp)
        
        print(f"📄 Text report saved: {text_file}")
        print(f"📊 JSON data saved: {json_file}")
        print(f"📈 Plots saved in reports/")
        
    def create_text_report(self):
        """Create formatted text report"""
        report = []
        report.append("POLYMER INFORMATICS MODEL ACCURACY REPORT")
        report.append("=" * 60)
        report.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        report.append(f"Test Data: {len(POLYMER_TEST_DATA)} known polymers")
        report.append("")
        
        # Summary table
        report.append("PROPERTY PREDICTION ACCURACY SUMMARY")
        report.append("-" * 60)
        report.append(f"{'Property':<12} {'Tests':<6} {'MAE':<8} {'RMSE':<8} {'R²':<8} {'±10°C':<8}")
        report.append("-" * 60)
        
        for prop, results in self.all_results.items():
            metrics = results['metrics']
            report.append(
                f"{prop:<12} "
                f"{results['test_count']:<6} "
                f"{metrics['mae']:<8.1f} "
                f"{metrics['rmse']:<8.1f} "
                f"{metrics['r2']:<8.3f} "
                f"{metrics['accuracy_10c']:<8.1f}%"
            )
        
        report.append("")
        
        # Detailed results for each property
        for prop, results in self.all_results.items():
            report.append(f"\nDETAILED RESULTS: {prop}")
            report.append("-" * 40)
            
            metrics = results['metrics']
            report.append(f"Mean Absolute Error:      {metrics['mae']:.2f}°C")
            report.append(f"Root Mean Square Error:   {metrics['rmse']:.2f}°C")
            report.append(f"R² Correlation:           {metrics['r2']:.3f}")
            report.append(f"Maximum Error:            {metrics['max_error']:.2f}°C")
            report.append(f"Median Error:             {metrics['median_error']:.2f}°C")
            report.append("")
            report.append("Accuracy within tolerance:")
            report.append(f"  ±5°C:   {metrics['accuracy_5c']:.1f}%")
            report.append(f"  ±10°C:  {metrics['accuracy_10c']:.1f}%")
            report.append(f"  ±20°C:  {metrics['accuracy_20c']:.1f}%")
            report.append("")
            
            # Individual predictions
            report.append("Individual Predictions:")
            report.append(f"{'Polymer':<25} {'Actual':<8} {'Predicted':<10} {'Error':<8}")
            report.append("-" * 55)
            
            for result in results['results'][:10]:  # Show first 10
                report.append(
                    f"{result['polymer_name']:<25} "
                    f"{result['actual_value']:<8.1f} "
                    f"{result['predicted_value']:<10.1f} "
                    f"{result['error']:<8.1f}"
                )
            
            if len(results['results']) > 10:
                report.append(f"... and {len(results['results']) - 10} more")
            
            report.append("")
        
        # Overall assessment
        report.append("OVERALL MODEL ASSESSMENT")
        report.append("-" * 40)
        
        avg_mae = np.mean([r['metrics']['mae'] for r in self.all_results.values()])
        avg_r2 = np.mean([r['metrics']['r2'] for r in self.all_results.values()])
        
        if avg_mae < 15 and avg_r2 > 0.7:
            assessment = "🟢 EXCELLENT - Model shows high accuracy across properties"
        elif avg_mae < 30 and avg_r2 > 0.5:
            assessment = "🟡 GOOD - Model performance is acceptable"
        elif avg_mae < 50 and avg_r2 > 0.3:
            assessment = "🟠 FAIR - Model needs improvement"
        else:
            assessment = "🔴 POOR - Model requires significant refinement"
        
        report.append(f"Average MAE: {avg_mae:.1f}°C")
        report.append(f"Average R²:  {avg_r2:.3f}")
        report.append(f"Assessment:  {assessment}")
        
        return "\n".join(report)
    
    def create_visualizations(self, timestamp):
        """Create accuracy visualization plots"""
        try:
            # Prediction vs Actual scatter plot
            fig, axes = plt.subplots(2, 2, figsize=(12, 10))
            fig.suptitle('Polymer Property Prediction Accuracy', fontsize=16)
            
            plot_idx = 0
            for prop, results in self.all_results.items():
                if plot_idx >= 4:
                    break
                    
                row, col = plot_idx // 2, plot_idx % 2
                ax = axes[row, col]
                
                actual = [r['actual_value'] for r in results['results']]
                predicted = [r['predicted_value'] for r in results['results']]
                
                # Scatter plot
                ax.scatter(actual, predicted, alpha=0.7)
                
                # Perfect prediction line
                min_val = min(min(actual), min(predicted))
                max_val = max(max(actual), max(predicted))
                ax.plot([min_val, max_val], [min_val, max_val], 'r--', alpha=0.8)
                
                ax.set_xlabel(f'Actual {prop} (°C)')
                ax.set_ylabel(f'Predicted {prop} (°C)')
                ax.set_title(f'{prop} Predictions (R² = {results["metrics"]["r2"]:.3f})')
                ax.grid(True, alpha=0.3)
                
                plot_idx += 1
            
            # Hide unused subplots
            for i in range(plot_idx, 4):
                row, col = i // 2, i % 2
                axes[row, col].set_visible(False)
            
            plt.tight_layout()
            plt.savefig(f'reports/accuracy_plots_{timestamp}.png', dpi=300, bbox_inches='tight')
            plt.close()
            
            # Error distribution histogram
            plt.figure(figsize=(10, 6))
            for prop, results in self.all_results.items():
                errors = [r['error'] for r in results['results']]
                plt.hist(errors, bins=15, alpha=0.7, label=f'{prop} (MAE: {results["metrics"]["mae"]:.1f}°C)')
            
            plt.xlabel('Prediction Error (°C)')
            plt.ylabel('Frequency')
            plt.title('Distribution of Prediction Errors')
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.savefig(f'reports/error_distribution_{timestamp}.png', dpi=300, bbox_inches='tight')
            plt.close()
            
        except ImportError:
            print("⚠️  Matplotlib not available, skipping plots")
        except Exception as e:
            print(f"⚠️  Error creating plots: {e}")

def main():
    """Run the complete test suite"""
    suite = TestSuite()
    suite.run_all_tests()

if __name__ == "__main__":
    main()
