#!/usr/bin/env python3
"""
Database setup script for ChatGrid backend
This script creates the necessary tables and imports CSV data into PostgreSQL
"""

import pandas as pd
import psycopg2
from sqlalchemy import create_engine
import config
import os

def setup_database():
    """Setup the PostgreSQL database with required tables and data"""
    
    print("Setting up database...")
    
    # Database connection parameters
    db_params = {
        'host': 'westmap-database',
        'port': 5432,
        'database': config.database_name,
        'user': 'postgres',
        'password': config.sql_key
    }
    
    # Create SQLAlchemy engine
    engine = create_engine(f"postgresql+psycopg2://postgres:{config.sql_key}@westmap-database:5432/{config.database_name}")
    
    # Path to CSV files
    data_dir = "../data"
    
    try:
        # Test database connection
        print(f"Connecting to database: {config.database_name}")
        conn = psycopg2.connect(**db_params)
        cursor = conn.cursor()
        
        print("✓ Database connection successful!")
        
        # Check if tables exist and get their info
        cursor.execute("""
            SELECT table_name 
            FROM information_schema.tables 
            WHERE table_schema = 'public'
        """)
        existing_tables = [row[0] for row in cursor.fetchall()]
        print(f"Existing tables: {existing_tables}")
        
        # Import CSV files to database
        csv_files = {
            'generation': os.path.join(data_dir, 'generation.csv'),
            'bus': os.path.join(data_dir, 'bus.csv'), 
            'transmission_line': os.path.join(data_dir, 'transmission_line.csv'),
            'counties': os.path.join(data_dir, 'counties.csv'),
            'us_states': os.path.join(data_dir, 'us states.csv')
        }
        
        for table_name, csv_path in csv_files.items():
            if os.path.exists(csv_path):
                print(f"\nImporting {csv_path} to table '{table_name}'...")
                
                # Read CSV file
                df = pd.read_csv(csv_path)
                print(f"  - Rows: {len(df)}, Columns: {list(df.columns)}")
                
                # Import to database (replace existing table)
                df.to_sql(table_name, engine, if_exists='replace', index=False)
                print(f"  ✓ Successfully imported to '{table_name}' table")
            else:
                print(f"  ⚠ CSV file not found: {csv_path}")
        
        # Show final table list
        cursor.execute("""
            SELECT table_name 
            FROM information_schema.tables 
            WHERE table_schema = 'public'
        """)
        final_tables = [row[0] for row in cursor.fetchall()]
        print(f"\nFinal tables in database: {final_tables}")
        
        # Show sample data from generation table
        if 'generation' in final_tables:
            cursor.execute("SELECT * FROM generation LIMIT 3")
            sample_data = cursor.fetchall()
            cursor.execute("SELECT column_name FROM information_schema.columns WHERE table_name = 'generation'")
            columns = [row[0] for row in cursor.fetchall()]
            print(f"\nSample data from 'generation' table:")
            print(f"Columns: {columns}")
            for row in sample_data:
                print(f"  {row}")
        
        cursor.close()
        conn.close()
        
        print("\n✓ Database setup completed successfully!")
        print("\nYou can now run the backend server with: python server.py")
        
    except psycopg2.Error as e:
        print(f"❌ Database error: {e}")
        print("\nPlease check:")
        print("1. PostgreSQL is running")
        print("2. Database exists and credentials are correct")
        print("3. You have permission to create tables")
        
    except Exception as e:
        print(f"❌ Error: {e}")

if __name__ == "__main__":
    setup_database()
