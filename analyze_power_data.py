#!/usr/bin/env python3
"""
Comprehensive Power Plant Data Analysis
Calculates accurate totals from Western Power Plants CSV files
"""

import pandas as pd
import numpy as np
import json
from pathlib import Path

def load_and_analyze_usa_data():
    """Load and analyze USA power plants data"""
    print("📊 Analyzing USA Power Plants Data...")
    
    usa_file = "viz/data/Western_Power_plants_Locations-USA (2).csv"
    df = pd.read_csv(usa_file)
    
    print(f"Total USA Plants: {len(df)}")
    print(f"Total USA Capacity: {df['Total Capacity (MW)'].sum():.1f} MW")
    
    # Technology breakdown with individual technology columns
    tech_columns = ['Battery Storage', 'Geothermal', 'Hydro', 'Natural Gas', 
                   'Nuclear', 'Oil', 'Other', 'Solar', 'Wind']
    
    usa_totals = {}
    for tech in tech_columns:
        if tech in df.columns:
            total = df[tech].sum()
            count = (df[tech] > 0).sum()
            usa_totals[tech] = {'capacity': total, 'count': count}
            print(f"  {tech}: {total:.1f} MW ({count} plants)")
    
    # Primary type analysis
    primary_totals = df.groupby('Primary Type').agg({
        'Total Capacity (MW)': ['sum', 'count']
    }).round(1)
    
    print("\n📈 USA Primary Type Analysis:")
    for tech in primary_totals.index:
        capacity = primary_totals.loc[tech, ('Total Capacity (MW)', 'sum')]
        count = primary_totals.loc[tech, ('Total Capacity (MW)', 'count')]
        print(f"  {tech}: {capacity:.1f} MW ({count} plants)")
    
    # State breakdown
    state_totals = df.groupby('State').agg({
        'Total Capacity (MW)': ['sum', 'count']
    }).round(1)
    
    print("\n🗺️ USA State Breakdown:")
    for state in state_totals.index:
        capacity = state_totals.loc[state, ('Total Capacity (MW)', 'sum')]
        count = state_totals.loc[state, ('Total Capacity (MW)', 'count')]
        print(f"  {state}: {capacity:.1f} MW ({count} plants)")
    
    return df, usa_totals, primary_totals

def load_and_analyze_canada_mexico_data():
    """Load and analyze Canada-Mexico power plants data"""
    print("\n📊 Analyzing Canada-Mexico Power Plants Data...")
    
    intl_file = "viz/data/Western_Power_plants_Locations-CANADA-MEX_Update (1).csv"
    df = pd.read_csv(intl_file)
    
    print(f"Total International Plants: {len(df)}")
    print(f"Total International Capacity: {df['Total Capacity (MW)'].sum():.1f} MW")
    
    # Primary type analysis
    primary_totals = df.groupby('Primary Type').agg({
        'Total Capacity (MW)': ['sum', 'count']
    }).round(1)
    
    print("\n📈 International Primary Type Analysis:")
    for tech in primary_totals.index:
        capacity = primary_totals.loc[tech, ('Total Capacity (MW)', 'sum')]
        count = primary_totals.loc[tech, ('Total Capacity (MW)', 'count')]
        print(f"  {tech}: {capacity:.1f} MW ({count} plants)")
    
    # Country breakdown
    country_totals = df.groupby('Country').agg({
        'Total Capacity (MW)': ['sum', 'count']
    }).round(1)
    
    print("\n🌍 Country Breakdown:")
    for country in country_totals.index:
        capacity = country_totals.loc[country, ('Total Capacity (MW)', 'sum')]
        count = country_totals.loc[country, ('Total Capacity (MW)', 'count')]
        print(f"  {country}: {capacity:.1f} MW ({count} plants)")
    
    return df, primary_totals

def load_and_analyze_wecc_data():
    """Load and analyze WECC balancing authority data"""
    print("\n📊 Analyzing WECC Generation Data...")
    
    wecc_file = "viz/amin_data/new_data/WECC_31_BAs_2028_All_Clean (1).csv"
    df = pd.read_csv(wecc_file)
    
    print(f"Total WECC Balancing Authorities: {len(df)}")
    
    # Technology columns in WECC data
    tech_columns = ['Hydro', 'Nuclear', 'Coal', 'Natural Gas', 'Geothermal', 
                   'Biomass', 'Wind', 'PV', 'Pumped_Storage_MW', 'Battery_Storage_MW']
    
    wecc_totals = {}
    total_wecc_capacity = 0
    
    print("\n📈 WECC Technology Totals:")
    for tech in tech_columns:
        if tech in df.columns:
            total = df[tech].sum()
            count = (df[tech] > 0).sum()
            wecc_totals[tech] = {'capacity': total, 'count': count}
            total_wecc_capacity += total
            print(f"  {tech}: {total:.1f} MW ({count} BAs)")
    
    print(f"\nTotal WECC Capacity: {total_wecc_capacity:.1f} MW")
    
    # BA breakdown
    print("\n🏢 WECC Balancing Authority Breakdown:")
    for _, row in df.iterrows():
        ba = row['BA']
        total_ba = row['Total_MW']
        print(f"  {ba}: {total_ba:.1f} MW")
    
    return df, wecc_totals

def calculate_combined_totals(usa_df, intl_df, usa_totals):
    """Calculate combined totals from all power plant data"""
    print("\n🔄 Calculating Combined Western Power Plants Totals...")
    
    # Combine both datasets
    combined_capacity = usa_df['Total Capacity (MW)'].sum() + intl_df['Total Capacity (MW)'].sum()
    combined_count = len(usa_df) + len(intl_df)
    
    print(f"Total Combined Plants: {combined_count}")
    print(f"Total Combined Capacity: {combined_capacity:.1f} MW")
    
    # Technology totals from USA detailed breakdown
    combined_tech_totals = {}
    
    # Get USA technology totals from individual columns
    for tech, data in usa_totals.items():
        combined_tech_totals[tech] = data.copy()
    
    # Add international plants by primary type
    intl_primary = intl_df.groupby('Primary Type').agg({
        'Total Capacity (MW)': ['sum', 'count']
    })
    
    print("\n📈 Combined Technology Analysis:")
    # Map international primary types to USA technology columns
    tech_mapping = {
        'Solar': 'Solar',
        'Wind': 'Wind', 
        'Hydro': 'Hydro',
        'Natural Gas': 'Natural Gas',
        'Nuclear': 'Nuclear',
        'Geothermal': 'Geothermal',
        'Biomass': 'Other',  # Map biomass to other for now
        'Oil': 'Oil',
        'Battery': 'Battery Storage',
        'Coal': 'Other'  # Map coal to other
    }
    
    for intl_tech in intl_primary.index:
        capacity = intl_primary.loc[intl_tech, ('Total Capacity (MW)', 'sum')]
        count = intl_primary.loc[intl_tech, ('Total Capacity (MW)', 'count')]
        
        # Map to USA technology column
        usa_tech = tech_mapping.get(intl_tech, 'Other')
        
        if usa_tech in combined_tech_totals:
            combined_tech_totals[usa_tech]['capacity'] += capacity
            combined_tech_totals[usa_tech]['count'] += count
        else:
            combined_tech_totals[usa_tech] = {'capacity': capacity, 'count': count}
    
    # Print combined totals
    total_check = 0
    for tech, data in combined_tech_totals.items():
        capacity = data['capacity']
        count = data['count']
        total_check += capacity
        percentage = (capacity / combined_capacity) * 100
        print(f"  {tech}: {capacity:.1f} MW ({count} plants) - {percentage:.1f}%")
    
    print(f"\nVerification - Total from technologies: {total_check:.1f} MW")
    print(f"Difference: {abs(combined_capacity - total_check):.1f} MW")
    
    return combined_tech_totals, combined_capacity, combined_count

def generate_chart_data(tech_totals, title="Western Power Plants"):
    """Generate chart data for the frontend"""
    
    # Sort technologies by capacity
    sorted_techs = sorted(tech_totals.items(), key=lambda x: x[1]['capacity'], reverse=True)
    
    labels = []
    data = []
    colors = []
    
    # Color mapping for technologies
    color_map = {
        'Solar': 'rgba(255, 193, 7, 0.8)',
        'Wind': 'rgba(76, 175, 80, 0.8)', 
        'Hydro': 'rgba(33, 150, 243, 0.8)',
        'Natural Gas': 'rgba(255, 87, 34, 0.8)',
        'Nuclear': 'rgba(156, 39, 176, 0.8)',
        'Geothermal': 'rgba(139, 69, 19, 0.8)',
        'Battery Storage': 'rgba(255, 235, 59, 0.8)',
        'Oil': 'rgba(96, 125, 139, 0.8)',
        'Other': 'rgba(158, 158, 158, 0.8)'
    }
    
    for tech, data_dict in sorted_techs:
        if data_dict['capacity'] > 0:  # Only include technologies with capacity
            labels.append(tech)
            data.append(round(data_dict['capacity'], 1))
            colors.append(color_map.get(tech, 'rgba(158, 158, 158, 0.8)'))
    
    chart_data = {
        'labels': labels,
        'datasets': [{
            'label': f'{title} Capacity (MW)',
            'data': data,
            'backgroundColor': colors,
            'borderColor': [c.replace('0.8', '1') for c in colors],
            'borderWidth': 2
        }]
    }
    
    return chart_data

def main():
    """Main analysis function"""
    print("🚀 Starting Comprehensive Power Plant Data Analysis")
    print("=" * 60)
    
    # Analyze USA data
    usa_df, usa_totals, usa_primary = load_and_analyze_usa_data()
    
    # Analyze Canada-Mexico data  
    intl_df, intl_primary = load_and_analyze_canada_mexico_data()
    
    # Analyze WECC data
    wecc_df, wecc_totals = load_and_analyze_wecc_data()
    
    # Calculate combined totals
    combined_totals, total_capacity, total_count = calculate_combined_totals(usa_df, intl_df, usa_totals)
    
    # Generate chart data
    western_chart_data = generate_chart_data(combined_totals, "Western Power Plants")
    wecc_chart_data = generate_chart_data(wecc_totals, "WECC Generation")
    
    # Save results
    results = {
        'western_power_plants': {
            'total_capacity_mw': total_capacity,
            'total_plants': total_count,
            'average_capacity_mw': total_capacity / total_count,
            'technology_breakdown': combined_totals,
            'chart_data': western_chart_data
        },
        'wecc_generation': {
            'total_capacity_mw': sum([data['capacity'] for data in wecc_totals.values()]),
            'total_bas': len(wecc_df),
            'technology_breakdown': wecc_totals,
            'chart_data': wecc_chart_data
        }
    }
    
    # Save to JSON file
    with open('power_plant_analysis.json', 'w') as f:
        json.dump(results, f, indent=2, default=str)
    
    print("\n" + "=" * 60)
    print("📄 Analysis Complete! Results saved to power_plant_analysis.json")
    print("\n🎯 KEY FINDINGS:")
    print(f"• Western Power Plants Total: {total_capacity:.1f} MW ({total_count} plants)")
    print(f"• WECC Generation Total: {sum([data['capacity'] for data in wecc_totals.values()]):.1f} MW")
    print(f"• Average Plant Size: {total_capacity / total_count:.1f} MW")
    
    # Top 5 technologies
    top_5 = sorted(combined_totals.items(), key=lambda x: x[1]['capacity'], reverse=True)[:5]
    print("\n🏆 Top 5 Technologies by Capacity:")
    for i, (tech, data) in enumerate(top_5, 1):
        percentage = (data['capacity'] / total_capacity) * 100
        print(f"  {i}. {tech}: {data['capacity']:.1f} MW ({percentage:.1f}%)")

if __name__ == "__main__":
    main()
