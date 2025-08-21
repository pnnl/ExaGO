#!/usr/bin/env python3
"""
Test script to verify the ChatGrid backend is working
"""

import requests
import json

def test_backend():
    """Test the backend with a simple query"""
    
    url = "http://127.0.0.1:5000/data"
    
    # Test queries
    test_queries = [
        "How many generators are there?",
        "Show me all hydro generators",
        "What are the different generation types?"
    ]
    
    for query in test_queries:
        print(f"\n🔍 Testing query: '{query}'")
        
        try:
            response = requests.post(url, json={"inputText": query})
            
            if response.status_code == 200:
                result = response.json()
                print(f"✅ Success!")
                print(f"📝 Response: {result.get('text', 'No text response')}")
                print(f"📊 Data points: {len(result.get('result_list', []))}")
                
                # Show first few results if available
                data = result.get('result_list', [])
                if data:
                    print(f"📋 Sample results:")
                    for i, item in enumerate(data[:3]):
                        print(f"   {i+1}: {item}")
                    if len(data) > 3:
                        print(f"   ... and {len(data)-3} more results")
            else:
                print(f"❌ Error: Status code {response.status_code}")
                print(f"Response: {response.text}")
                
        except requests.exceptions.ConnectionError:
            print("❌ Cannot connect to backend server. Make sure it's running on http://127.0.0.1:5000")
        except Exception as e:
            print(f"❌ Error: {e}")

if __name__ == "__main__":
    print("🧪 Testing ChatGrid Backend...")
    test_backend()
    print("\n✅ Test completed!")
