#!/usr/bin/env python3
"""
Data quality analysis for PDF-extracted polymer data
Identifies potential issues with the training data
"""
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import requests
import numpy as np
import json

class DataQualityAnalyzer:
    def __init__(self, csv_file="scraped_polymer_data.csv"):
        self.csv_file = csv_file
        self.issues = []
        
    def analyze_csv_data(self):
        """Analyze the scraped CSV data for quality issues"""
        print("🔍 Analyzing Data Quality in scraped_polymer_data.csv")
        print("=" * 60)
        
        try:
            import pandas as pd
            df = pd.read_csv(self.csv_file)
        except ImportError:
            print("❌ pandas not available, analyzing manually...")
            return self.analyze_csv_manual()
        except FileNotFoundError:
            print(f"❌ CSV file not found: {self.csv_file}")
            return
        
        # Basic statistics
        print(f"📊 Dataset Overview:")
        print(f"   Total records: {len(df)}")
        print(f"   Columns: {list(df.columns)}")
        print()
        
        # Analyze each property
        property_columns = ['Tg', 'Tm', 'elastic_modulus', 'molecular_weight', 'density']
        
        for prop in property_columns:
            if prop in df.columns:
                self.analyze_property(df, prop)
        
        # Identify outliers and anomalies
        self.identify_outliers(df)
        
        # Generate recommendations
        self.generate_recommendations()
        
    def analyze_property(self, df, property_name):
        """Analyze a specific property for quality issues"""
        prop_data = df[property_name].dropna()
        
        if len(prop_data) == 0:
            print(f"⚠️  {property_name}: No data available")
            return
        
        print(f"📈 {property_name} Analysis:")
        print(f"   Records with data: {len(prop_data)}/{len(df)} ({len(prop_data)/len(df)*100:.1f}%)")
        print(f"   Min value: {prop_data.min():.2f}")
        print(f"   Max value: {prop_data.max():.2f}")
        print(f"   Mean: {prop_data.mean():.2f}")
        print(f"   Median: {prop_data.median():.2f}")
        print(f"   Std Dev: {prop_data.std():.2f}")
        
        # Check for unrealistic values
        if property_name == "Tg":
            suspicious = prop_data[(prop_data < -200) | (prop_data > 500)]
            if len(suspicious) > 0:
                print(f"   ⚠️  {len(suspicious)} suspicious Tg values (outside -200°C to 500°C)")
                self.issues.append(f"Suspicious {property_name} values: {suspicious.tolist()[:5]}")
        
        elif property_name == "Tm":
            suspicious = prop_data[(prop_data < -100) | (prop_data > 600)]
            if len(suspicious) > 0:
                print(f"   ⚠️  {len(suspicious)} suspicious Tm values (outside -100°C to 600°C)")
                self.issues.append(f"Suspicious {property_name} values: {suspicious.tolist()[:5]}")
        
        elif property_name == "molecular_weight":
            suspicious = prop_data[(prop_data < 100) | (prop_data > 1000000)]
            if len(suspicious) > 0:
                print(f"   ⚠️  {len(suspicious)} suspicious MW values (outside 100 to 1M g/mol)")
                self.issues.append(f"Suspicious {property_name} values: {suspicious.tolist()[:5]}")
        
        print()
    
    def identify_outliers(self, df):
        """Identify statistical outliers in the data"""
        print("🎯 Outlier Detection:")
        
        numeric_columns = df.select_dtypes(include=[np.number]).columns
        
        for col in numeric_columns:
            data = df[col].dropna()
            if len(data) > 5:  # Need enough data points
                Q1 = data.quantile(0.25)
                Q3 = data.quantile(0.75)
                IQR = Q3 - Q1
                
                # Define outliers as values outside 1.5*IQR
                lower_bound = Q1 - 1.5 * IQR
                upper_bound = Q3 + 1.5 * IQR
                
                outliers = data[(data < lower_bound) | (data > upper_bound)]
                
                if len(outliers) > 0:
                    print(f"   {col}: {len(outliers)} outliers detected")
                    print(f"     Range: {outliers.min():.2f} to {outliers.max():.2f}")
                    self.issues.append(f"{col} has {len(outliers)} statistical outliers")
        
        print()
    
    def generate_recommendations(self):
        """Generate data quality improvement recommendations"""
        print("💡 Data Quality Recommendations:")
        print("-" * 40)
        
        if not self.issues:
            print("✅ No major data quality issues detected!")
            return
        
        recommendations = []
        
        # Based on identified issues
        for issue in self.issues:
            if "suspicious" in issue.lower():
                recommendations.append("🔧 Review and validate extreme property values")
            if "outliers" in issue.lower():
                recommendations.append("📊 Consider removing or investigating statistical outliers")
        
        # General recommendations
        recommendations.extend([
            "🧹 Implement data validation rules during PDF extraction",
            "📖 Cross-reference values with literature databases",
            "🔍 Add confidence scoring for extracted values",
            "📝 Include source context for better extraction accuracy",
            "⚖️ Implement unit detection and conversion",
            "🎯 Focus on high-confidence data for model training"
        ])
        
        for i, rec in enumerate(set(recommendations), 1):
            print(f"{i}. {rec}")
        
        print()
        print("🎯 Priority Actions:")
        print("1. Filter training data to remove obvious errors")
        print("2. Implement robust data validation in the scraper")
        print("3. Add literature-based ground truth validation")
    
    def analyze_csv_manual(self):
        """Manual CSV analysis when pandas is not available"""
        try:
            with open(self.csv_file, 'r') as f:
                lines = f.readlines()
            
            header = lines[0].strip().split(',')
            data_lines = lines[1:]
            
            print(f"📊 Manual Analysis:")
            print(f"   Total records: {len(data_lines)}")
            print(f"   Columns: {header}")
            
            # Look for property columns
            tg_values = []
            for line in data_lines:
                fields = line.strip().split(',')
                if len(fields) > 5 and fields[5]:  # Tg column
                    try:
                        tg_values.append(float(fields[5]))
                    except ValueError:
                        continue
            
            if tg_values:
                print(f"   Tg values found: {len(tg_values)}")
                print(f"   Tg range: {min(tg_values):.1f} to {max(tg_values):.1f}°C")
                
                # Check for suspicious values
                suspicious = [v for v in tg_values if v < -200 or v > 500]
                if suspicious:
                    print(f"   ⚠️  {len(suspicious)} suspicious Tg values")
            
        except Exception as e:
            print(f"❌ Error analyzing CSV: {e}")

def main():
    """Run data quality analysis"""
    analyzer = DataQualityAnalyzer()
    analyzer.analyze_csv_data()

if __name__ == "__main__":
    main()
