#!/usr/bin/env python3
"""
Script to create Price data CSV from BA Price data for all case studies
"""

import os
import csv

def read_csv_file(file_path):
    """Read CSV file and return data as list of dictionaries"""
    data = []
    try:
        with open(file_path, 'r', newline='', encoding='utf-8') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                data.append(row)
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
    return data

def main():
    # Base directory
    base_dir = '/Users/manish/Documents/U-smart_Projects/usmart-westmap/viz/amin_data/manish_amin_modified_data'
    
    # Case study directories
    case_studies = ['Case study_0', 'Case study_1', 'Case study_2', 'Case_study_3']
    
    # Initialize the final data
    all_data = []
    
    print("Processing price data for case studies...")
    
    for case_study in case_studies:
        print(f"Processing {case_study}...")
        
        price_file = os.path.join(base_dir, case_study, 'BA Price', 'lmp_by_ba_hour.csv')
        
        if not os.path.exists(price_file):
            print(f"Price file not found: {price_file}")
            continue
            
        # Read the price data
        price_data = read_csv_file(price_file)
        
        if not price_data:
            print(f"No price data found in {price_file}")
            continue
            
        # Process the data
        for row in price_data:
            try:
                processed_row = {
                    'Case_Study': case_study.replace('Case study_', 'Case_').replace('Case_study_', 'Case_'),
                    'Hour': int(row.get('hour', 0)),
                    'BA': row.get('ba', '').strip(),
                    'LMP': float(row.get('lmp', 0))
                }
                all_data.append(processed_row)
            except (ValueError, TypeError) as e:
                print(f"Error processing row in {case_study}: {e}")
                continue
    
    # Write to CSV
    output_file = os.path.join(base_dir, 'Price_Data.csv')
    
    headers = ['Case_Study', 'Hour', 'BA', 'LMP']
    
    with open(output_file, 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=headers)
        writer.writeheader()
        
        # Sort data by Case_Study, Hour, and BA
        all_data.sort(key=lambda x: (x['Case_Study'], x['Hour'], x['BA']))
        
        for row in all_data:
            writer.writerow(row)
    
    print(f"Created {output_file}")
    print(f"Total rows: {len(all_data)}")
    
    # Display sample data
    print("\nSample price data (first 5 rows):")
    for i, row in enumerate(all_data[:5]):
        print(f"Row {i+1}: Case={row['Case_Study']}, Hour={row['Hour']}, BA={row['BA']}, LMP=${row['LMP']:.2f}")

if __name__ == "__main__":
    main()
