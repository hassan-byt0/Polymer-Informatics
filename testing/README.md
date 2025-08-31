# Testing Framework for Polymer Property Predictions

This directory contains tools and datasets for testing the accuracy of our polymer informatics model predictions.

## Directory Structure

```
testing/
├── test_data/          # Curated test datasets with known properties
├── validation/         # Scripts for running validation tests
├── metrics/           # Accuracy calculation and statistical analysis
└── reports/           # Generated test reports and visualizations
```

## Test Data Sources

1. **Curated Polymer Database**: Known Tg, Tm values from literature
2. **Experimental Validation**: Real lab measurements for comparison
3. **Cross-validation**: Hold-out sets from our training data
4. **Edge Cases**: Unusual polymers to test model robustness

## Metrics Evaluated

- **Mean Absolute Error (MAE)**: Average prediction error
- **Root Mean Square Error (RMSE)**: Penalizes large errors
- **R² Score**: Correlation between predicted and actual values
- **Accuracy within ±X°C**: Percentage of predictions within tolerance
- **Confidence Calibration**: How well uncertainty estimates match actual errors

## Usage

```bash
# Run full validation suite
python validation/run_all_tests.py

# Test specific property
python validation/test_property.py --property Tg

# Generate accuracy report
python validation/generate_report.py
```

## Test Categories

1. **Unit Tests**: Individual component validation
2. **Integration Tests**: End-to-end pipeline testing  
3. **Accuracy Tests**: Prediction quality assessment
4. **Performance Tests**: Speed and scalability
5. **Robustness Tests**: Edge cases and error handling
