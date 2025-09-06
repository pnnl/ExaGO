#!/usr/bin/env python3
"""
Convert powerlines_WUS_CAN_sgca shapefile to CSV format with WGS84 coordinates.
This script reads the shapefile, reprojects to WGS84, and converts it to a CSV with WKT geometry format.
"""

import os
import csv
import sys

def convert_with_fiona():
    """Convert shapefile to CSV using Fiona with coordinate transformation"""
    try:
        import fiona
        from fiona.crs import from_epsg
        from shapely.geometry import shape
        from shapely.ops import transform
        import pyproj
        
        csv_data = []
        
        with fiona.open("powerlines_WUS_CAN_sgca.shp", "r") as shapefile:
            print(f"Schema: {shapefile.schema}")
            print(f"Original CRS: {shapefile.crs}")
            
            # Set up coordinate transformation
            # From NAD 1927 Albers to WGS84
            transformer = pyproj.Transformer.from_crs(
                shapefile.crs,
                "EPSG:4326",  # WGS84
                always_xy=True
            )
            
            feature_count = 0
            for record in shapefile:
                feature_count += 1
                # Process ALL features for complete WECC coverage including Canada/Mexico
                # if feature_count > 5000:  # Removed limit to process all features
                #     break
                    
                # Get geometry and transform to WGS84
                geom = shape(record['geometry'])
                geom_wgs84 = transform(transformer.transform, geom)
                wkt = geom_wgs84.wkt
                
                # Get properties
                props = record['properties']
                length = props.get('LENGTH', 0)
                ca_power = props.get('CA_POWER', 'UNKNOWN')
                source_file = props.get('SOURCEFILE', 'UNKNOWN')
                
                # Create meaningful line name
                line_name = f"{ca_power}_{feature_count}"
                if source_file and source_file != 'UNKNOWN':
                    line_name = f"{source_file}_{feature_count}"
                
                # Determine voltage level based on source
                voltage = 230  # Default
                if 'transmission' in str(source_file).lower():
                    voltage = 500
                elif 'distribution' in str(source_file).lower():
                    voltage = 138
                
                # Create CSV row
                csv_row = {
                    'wkt': wkt,
                    'flow capacity': max(length / 1000, 100) if length else 100,  # Convert length to rough capacity estimate
                    'pf': 0,
                    'qf': 0,
                    'pt': 0,
                    'qt': 0,
                    'kilovolt': voltage,
                    'line_name': line_name,
                    'source': str(ca_power)[:30] if ca_power else 'UNKNOWN',
                    'target': str(source_file)[:30] if source_file else 'UNKNOWN',
                    'actual flow': 0,
                    'data_source': 'powerlines_WUS_CAN_sgca'
                }
                csv_data.append(csv_row)
                
                if feature_count % 1000 == 0:
                    print(f"Processed {feature_count} features...")
        
        # Write CSV file
        if csv_data:
            with open('powerlines_WUS_CAN_sgca_wgs84.csv', 'w', newline='', encoding='utf-8') as csvfile:
                fieldnames = ['wkt', 'flow capacity', 'pf', 'qf', 'pt', 'qt', 'kilovolt', 'line_name', 'source', 'target', 'actual flow', 'data_source']
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(csv_data)
            
            print(f"Successfully converted {len(csv_data)} features to powerlines_WUS_CAN_sgca_wgs84.csv")
            return True
        else:
            print("No features found in shapefile")
            return False
            
    except ImportError as e:
        print(f"Required packages not available: {e}")
        return False
    except Exception as e:
        print(f"Error with conversion: {e}")
        return False

def main():
    print("Converting powerlines_WUS_CAN_sgca shapefile to WGS84 CSV format...")
    
    # Check if shapefile exists
    if not os.path.exists("powerlines_WUS_CAN_sgca.shp"):
        print("Error: powerlines_WUS_CAN_sgca.shp not found")
        return
    
    success = convert_with_fiona()
    
    if not success:
        print("\nTo properly convert the shapefile, install the required packages:")
        print("pip install fiona shapely pyproj")

if __name__ == "__main__":
    main()
