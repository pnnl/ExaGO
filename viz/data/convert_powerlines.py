#!/usr/bin/env python3
"""
Convert powerlines_WUS_CAN_sgca shapefile to CSV format compatible with the app.js visualization.
This script reads the shapefile and converts it to a CSV with WKT geometry format.
"""

import os
import csv
import sys

def convert_with_ogr():
    """Convert shapefile to CSV using GDAL/OGR if available"""
    try:
        from osgeo import ogr, osr
        
        # Open the shapefile
        driver = ogr.GetDriverByName("ESRI Shapefile")
        dataset = driver.Open("powerlines_WUS_CAN_sgca.shp", 0)
        
        if dataset is None:
            print("Error: Could not open shapefile")
            return False
        
        layer = dataset.GetLayer()
        
        # Get the spatial reference system
        srs = layer.GetSpatialRef()
        if srs:
            print(f"Spatial Reference: {srs.ExportToWkt()}")
        
        # Prepare CSV output
        csv_data = []
        
        # Get layer definition
        layer_defn = layer.GetLayerDefn()
        field_names = []
        for i in range(layer_defn.GetFieldCount()):
            field_names.append(layer_defn.GetFieldDefn(i).GetName())
        
        print(f"Field names: {field_names}")
        
        # Read features
        feature_count = 0
        for feature in layer:
            feature_count += 1
            if feature_count > 10:  # Limit for testing
                break
                
            # Get geometry as WKT
            geom = feature.GetGeometryRef()
            if geom:
                wkt = geom.ExportToWkt()
                
                # Get attributes
                properties = {}
                for field_name in field_names:
                    properties[field_name] = feature.GetField(field_name)
                
                # Create CSV row similar to existing format
                csv_row = {
                    'wkt': wkt,
                    'flow capacity': properties.get('LENGTH', 0),  # Using LENGTH as a proxy
                    'pf': 0,  # Default values since not available in source
                    'qf': 0,
                    'pt': 0,
                    'qt': 0,
                    'kilovolt': 230,  # Default voltage level
                    'line_name': f"POWERLINE_{feature_count}",
                    'source': properties.get('CA_POWER', 'UNKNOWN'),
                    'target': properties.get('SOURCEFILE', 'UNKNOWN'),
                    'actual flow': 0,
                    'data_source': 'powerlines_WUS_CAN_sgca'
                }
                csv_data.append(csv_row)
        
        # Write CSV file
        if csv_data:
            with open('powerlines_WUS_CAN_sgca_converted.csv', 'w', newline='', encoding='utf-8') as csvfile:
                fieldnames = ['wkt', 'flow capacity', 'pf', 'qf', 'pt', 'qt', 'kilovolt', 'line_name', 'source', 'target', 'actual flow', 'data_source']
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(csv_data)
            
            print(f"Successfully converted {len(csv_data)} features to powerlines_WUS_CAN_sgca_converted.csv")
            return True
        else:
            print("No features found in shapefile")
            return False
            
    except ImportError:
        print("GDAL/OGR not available, trying alternative method...")
        return False
    except Exception as e:
        print(f"Error with GDAL/OGR: {e}")
        return False

def convert_with_fiona():
    """Convert shapefile to CSV using Fiona if available"""
    try:
        import fiona
        from shapely.geometry import shape
        
        csv_data = []
        
        with fiona.open("powerlines_WUS_CAN_sgca.shp", "r") as shapefile:
            print(f"Schema: {shapefile.schema}")
            print(f"CRS: {shapefile.crs}")
            
            feature_count = 0
            for record in shapefile:
                feature_count += 1
                if feature_count > 1000:  # Limit for performance
                    break
                    
                # Get geometry as WKT
                geom = shape(record['geometry'])
                wkt = geom.wkt
                
                # Get properties
                props = record['properties']
                
                # Create CSV row
                csv_row = {
                    'wkt': wkt,
                    'flow capacity': props.get('LENGTH', 0),
                    'pf': 0,
                    'qf': 0,
                    'pt': 0,
                    'qt': 0,
                    'kilovolt': 230,  # Default voltage
                    'line_name': f"POWERLINE_{feature_count}",
                    'source': str(props.get('CA_POWER', 'UNKNOWN'))[:20],  # Truncate long names
                    'target': str(props.get('SOURCEFILE', 'UNKNOWN'))[:20],
                    'actual flow': 0,
                    'data_source': 'powerlines_WUS_CAN_sgca'
                }
                csv_data.append(csv_row)
        
        # Write CSV file
        if csv_data:
            with open('powerlines_WUS_CAN_sgca_converted.csv', 'w', newline='', encoding='utf-8') as csvfile:
                fieldnames = ['wkt', 'flow capacity', 'pf', 'qf', 'pt', 'qt', 'kilovolt', 'line_name', 'source', 'target', 'actual flow', 'data_source']
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(csv_data)
            
            print(f"Successfully converted {len(csv_data)} features to powerlines_WUS_CAN_sgca_converted.csv")
            return True
        else:
            print("No features found in shapefile")
            return False
            
    except ImportError:
        print("Fiona not available, trying alternative method...")
        return False
    except Exception as e:
        print(f"Error with Fiona: {e}")
        return False

def simple_dbf_reader():
    """Simple DBF reader to extract basic information"""
    try:
        # Try to read the DBF file directly
        import struct
        
        with open('powerlines_WUS_CAN_sgca.dbf', 'rb') as f:
            # Read DBF header
            header = f.read(32)
            if len(header) < 32:
                print("Invalid DBF file")
                return False
            
            # Parse header
            version = header[0]
            num_records = struct.unpack('<I', header[4:8])[0]
            header_length = struct.unpack('<H', header[8:10])[0]
            record_length = struct.unpack('<H', header[10:12])[0]
            
            print(f"DBF Info: {num_records} records, header length: {header_length}, record length: {record_length}")
            
            # This is a basic implementation - for production use a proper DBF library
            print("Basic DBF structure detected. For full conversion, install GDAL or Fiona.")
            return False
            
    except Exception as e:
        print(f"Error reading DBF: {e}")
        return False

def main():
    print("Converting powerlines_WUS_CAN_sgca shapefile to CSV format...")
    
    # Check if shapefile exists
    if not os.path.exists("powerlines_WUS_CAN_sgca.shp"):
        print("Error: powerlines_WUS_CAN_sgca.shp not found")
        return
    
    # Try different conversion methods
    success = False
    
    print("Trying GDAL/OGR...")
    success = convert_with_ogr()
    
    if not success:
        print("Trying Fiona...")
        success = convert_with_fiona()
    
    if not success:
        print("Trying simple DBF reader...")
        simple_dbf_reader()
    
    if not success:
        print("\nTo properly convert the shapefile, install one of the following:")
        print("1. GDAL: pip install gdal")
        print("2. Fiona: pip install fiona")
        print("3. GeoPandas: pip install geopandas")

if __name__ == "__main__":
    main()
