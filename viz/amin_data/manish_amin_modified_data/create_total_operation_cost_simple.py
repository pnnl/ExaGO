#!/usr/bin/env python3
"""
Script to create Total_Operation_Cost.csv from individual BA hourly operation costs
This will aggregate all WECC region operational costs for all case studies
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
    
    # WECC regions list (based on the files we've seen)
    wecc_regions = [
        'AESO', 'AVA', 'AZPS', 'BANC', 'BCHA', 'BPAT', 'CENACE', 'CHPD', 'CISO', 
        'DOPD', 'EPE', 'GCPD', 'IID', 'IPCO', 'LDWP', 'NEVP', 'NWMT', 'PACE', 
        'PACW', 'PGE', 'PNM', 'PSCO', 'PSEI', 'SCL', 'SRP', 'TEPC', 'TIDC', 
        'TPWR', 'WACM', 'WALC', 'WAUW'
    ]
    
    # Initialize the final data
    all_data = []
    
    print("Processing case studies...")
    
    for case_study in case_studies:
        print(f"Processing {case_study}...")
        
        case_dir = os.path.join(base_dir, case_study, 'Balancing Authority Hourly Operation Costs')
        
        if not os.path.exists(case_dir):
            print(f"Directory not found: {case_dir}")
            continue
            
        # Process each hour (1-24)
        for hour in range(1, 25):
            row_data = {
                'Case_Study': case_study.replace('Case study_', 'Case_').replace('Case_study_', 'Case_'),
                'Hour': hour
            }
            
            # Process each WECC region
            for region in wecc_regions:
                cost_file = os.path.join(case_dir, f'{region}_hourly_operation_costs.csv')
                
                if os.path.exists(cost_file):
                    try:
                        # Read the CSV file
                        data = read_csv_file(cost_file)
                        
                        # Find the row for this hour
                        total_cost = 0
                        for row in data:
                            if int(row.get('Hour', 0)) == hour:
                                total_cost = float(row.get('BA_Total_Hourly_Cost_$', 0))
                                break
                        
                        row_data[region] = total_cost
                        
                    except Exception as e:
                        print(f"Error processing {cost_file}: {e}")
                        row_data[region] = 0
                else:
                    # File doesn't exist, set to 0
                    row_data[region] = 0
            
            all_data.append(row_data)
    
    # Write to CSV
    output_file = os.path.join(base_dir, 'Total_Operation_Cost.csv')
    
    # Prepare headers
    headers = ['Case_Study', 'Hour'] + wecc_regions
    
    with open(output_file, 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=headers)
        writer.writeheader()
        
        # Sort data by Case_Study and Hour
        all_data.sort(key=lambda x: (x['Case_Study'], x['Hour']))
        
        for row in all_data:
            writer.writerow(row)
    
    print(f"Created {output_file}")
    print(f"Total rows: {len(all_data)}")
    print(f"Columns: {len(headers)}")
    
    # Display sample data
    print("\nSample data (first 3 rows):")
    for i, row in enumerate(all_data[:3]):
        print(f"Row {i+1}: Case={row['Case_Study']}, Hour={row['Hour']}, AESO=${row.get('AESO', 0)}, CISO=${row.get('CISO', 0)}")

if __name__ == "__main__":
    main()
