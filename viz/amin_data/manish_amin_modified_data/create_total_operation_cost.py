#!/usr/bin/env python3
"""
Script to create Total_Operation_Cost.csv from individual BA hourly operation costs
This will aggregate all WECC region operational costs for all case studies
"""

import os
import pandas as pd
import glob

def main():
    # Base directory
    base_dir = '/Users/manish/Documents/U-smart_Projects/usmart-westmap/viz/amin_data/manish_amin_modified_data'
    
    # Case study directories
    case_studies = ['Case study_0', 'Case study_1', 'Case study_2', 'Case_study_3']
    
    # Initialize the final dataframe
    all_data = []
    
    # WECC regions list (based on the files we've seen)
    wecc_regions = [
        'AESO', 'AVA', 'AZPS', 'BANC', 'BCHA', 'BPAT', 'CENACE', 'CHPD', 'CISO', 
        'DOPD', 'EPE', 'GCPD', 'IID', 'IPCO', 'LDWP', 'NEVP', 'NWMT', 'PACE', 
        'PACW', 'PGE', 'PNM', 'PSCO', 'PSEI', 'SCL', 'SRP', 'TEPC', 'TIDC', 
        'TPWR', 'WACM', 'WALC', 'WAUW'
    ]
    
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
                        df = pd.read_csv(cost_file)
                        
                        # Find the row for this hour
                        hour_row = df[df['Hour'] == hour]
                        
                        if not hour_row.empty:
                            # Get the total cost for this hour
                            total_cost = hour_row['BA_Total_Hourly_Cost_$'].iloc[0]
                            row_data[region] = total_cost
                        else:
                            row_data[region] = 0
                            
                    except Exception as e:
                        print(f"Error reading {cost_file}: {e}")
                        row_data[region] = 0
                else:
                    # File doesn't exist, set to 0
                    row_data[region] = 0
            
            all_data.append(row_data)
    
    # Create DataFrame
    df_final = pd.DataFrame(all_data)
    
    # Sort by Case_Study and Hour
    df_final = df_final.sort_values(['Case_Study', 'Hour'])
    
    # Save to CSV
    output_file = os.path.join(base_dir, 'Total_Operation_Cost.csv')
    df_final.to_csv(output_file, index=False)
    
    print(f"Created {output_file}")
    print(f"Total rows: {len(df_final)}")
    print(f"Columns: {list(df_final.columns)}")
    
    # Display sample data
    print("\nSample data:")
    print(df_final.head())

if __name__ == "__main__":
    main()
