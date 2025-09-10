#!/usr/bin/env python3
"""
Script to create Power Exchange data CSV from generation data for all case studies
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
    
    # WECC regions list
    wecc_regions = [
        'AESO', 'AVA', 'AZPS', 'BANC', 'BCHA', 'BPAT', 'CENACE', 'CHPD', 'CISO', 
        'DOPD', 'EPE', 'GCPD', 'IID', 'IPCO', 'LDWP', 'NEVP', 'NWMT', 'PACE', 
        'PACW', 'PGE', 'PNM', 'PSCO', 'PSEI', 'SCL', 'SRP', 'TEPC', 'TIDC', 
        'TPWR', 'WACM', 'WALC', 'WAUW'
    ]
    
    # Initialize the final data
    all_data = []
    
    print("Processing power exchange data for case studies...")
    
    for case_study in case_studies:
        print(f"Processing {case_study}...")
        
        case_dir = os.path.join(base_dir, case_study, 'Balancing Authority Power Generation')
        
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
                gen_file = os.path.join(case_dir, f'{region}_generation_by_fuel.csv')
                
                if os.path.exists(gen_file):
                    try:
                        # Read the generation file
                        data = read_csv_file(gen_file)
                        
                        # Find the row for this hour
                        import_export = 0
                        for row in data:
                            if int(row.get('Hour', 0)) == hour:
                                import_export = float(row.get('IMPORT/EXPORT_MW', 0))
                                break
                        
                        row_data[region] = import_export
                        
                    except Exception as e:
                        print(f"Error processing {gen_file}: {e}")
                        row_data[region] = 0
                else:
                    # File doesn't exist, set to 0
                    row_data[region] = 0
            
            all_data.append(row_data)
    
    # Write to CSV
    output_file = os.path.join(base_dir, 'Power_Exchange_Data.csv')
    
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
    print("\nSample power exchange data (first 3 rows):")
    for i, row in enumerate(all_data[:3]):
        print(f"Row {i+1}: Case={row['Case_Study']}, Hour={row['Hour']}, AESO={row.get('AESO', 0):.2f} MW, CISO={row.get('CISO', 0):.2f} MW")

if __name__ == "__main__":
    main()
