#!/usr/bin/env python3
"""
Convert powerlines shapefile to CSV format with WKT geometry
"""

import struct
import os

def read_shp_header(shp_file):
    """Read SHP file header to get basic info"""
    with open(shp_file, 'rb') as f:
        # Read file code (should be 9994)
        file_code = struct.unpack('>I', f.read(4))[0]
        
        # Skip unused fields
        f.seek(24)
        
        # Read file length
        file_length = struct.unpack('>I', f.read(4))[0]
        
        # Read version
        version = struct.unpack('<I', f.read(4))[0]
        
        # Read shape type
        shape_type = struct.unpack('<I', f.read(4))[0]
        
        # Read bounding box
        xmin = struct.unpack('<d', f.read(8))[0]
        ymin = struct.unpack('<d', f.read(8))[0]
        xmax = struct.unpack('<d', f.read(8))[0]
        ymax = struct.unpack('<d', f.read(8))[0]
        
        return {
            'file_code': file_code,
            'file_length': file_length,
            'version': version,
            'shape_type': shape_type,
            'bounds': [xmin, ymin, xmax, ymax]
        }

def read_dbf_header(dbf_file):
    """Read DBF file header to get field information"""
    fields = []
    with open(dbf_file, 'rb') as f:
        # Read DBF header
        header = f.read(32)
        
        # Number of records
        num_records = struct.unpack('<I', header[4:8])[0]
        
        # Header length
        header_length = struct.unpack('<H', header[8:10])[0]
        
        # Record length
        record_length = struct.unpack('<H', header[10:12])[0]
        
        # Read field descriptors
        field_count = (header_length - 33) // 32
        
        for i in range(field_count):
            field_desc = f.read(32)
            field_name = field_desc[:11].decode('ascii').rstrip('\x00')
            field_type = chr(field_desc[11])
            field_length = field_desc[16]
            
            fields.append({
                'name': field_name,
                'type': field_type,
                'length': field_length
            })
    
    return {
        'num_records': num_records,
        'header_length': header_length,
        'record_length': record_length,
        'fields': fields
    }

def main():
    shp_file = 'data/powerlines_WUS_CAN_sgca.shp'
    dbf_file = 'data/powerlines_WUS_CAN_sgca.dbf'
    
    if not os.path.exists(shp_file) or not os.path.exists(dbf_file):
        print(f"Error: Required files not found")
        print(f"SHP exists: {os.path.exists(shp_file)}")
        print(f"DBF exists: {os.path.exists(dbf_file)}")
        return
    
    # Read headers
    print("Reading SHP header...")
    shp_info = read_shp_header(shp_file)
    print(f"Shape type: {shp_info['shape_type']}")
    print(f"Bounds: {shp_info['bounds']}")
    
    print("\nReading DBF header...")
    dbf_info = read_dbf_header(dbf_file)
    print(f"Number of records: {dbf_info['num_records']}")
    print(f"Fields: {[f['name'] for f in dbf_info['fields']]}")

if __name__ == '__main__':
    main()
