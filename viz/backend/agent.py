from pathlib import Path
from typing import List, Optional, Union, Dict, Any
from pydantic import BaseModel, Field
from agno.agent import Agent
from agno.tools.csv_toolkit import CsvTools
from agno.models.openai import OpenAIChat
from config import openai_key

# load all csv file paths in data directory
csv_list = list(Path("data").glob("*.csv"))

# Pydantic models for structured responses matching frontend visualization requirements
class GenerationDataItem(BaseModel):
    generation_name: str = Field(..., description="Name of the power generation facility", alias="generation name")

class TransmissionDataItem(BaseModel):
    line_name: str = Field(..., description="Name of the transmission line", alias="line name")

class BusDataItem(BaseModel):
    bus_name: str = Field(..., description="Name of the electrical bus/substation", alias="bus name")

class PowerGridResponse(BaseModel):
    result_list: List[Union[GenerationDataItem, TransmissionDataItem, BusDataItem]] = Field(
        ..., 
        description="List of power grid entities (generation facilities, transmission lines, or buses/substations)"
    )
    text: str = Field(..., description="Natural language summary of the results. but in markdown")
    query_type: str = Field(
        ..., 
        description="Type of query: 'generation', 'transmission', 'bus', or 'general'",
        enum=["generation", "transmission", "bus", "general"]
    )

instruction = [
    "You are a data analysis assistant for power grid visualization. You help users find specific power system data.",
    "INTRODUCE YOURSELF LIKE: 'Hello, I am a power grid data analysis assistant. I can help you with queries about power generation, transmission lines, and electrical substations.'",
    "First always get the list of files",
    "Then check the columns in the file",
    "Then run the query to answer the question",
    "Always wrap column names with double quotes if they contain spaces or special characters",
    "Remember to escape the quotes in the JSON string (use \")",
    "Use single quotes for string values",
    """
    IMPORTANT: Your responses must be structured for visualization mapping:
    
    For GENERATION queries (about power plants, generators, wind farms, solar, etc.):
    - Return a list with 'generation name' field (exact match from generation.csv)
    - Set query_type to 'generation'
    - Include power capacity, generation type, and location info in text summary
    
    For TRANSMISSION LINE queries (about power lines, transmission, connections):
    - Return a list with 'line name' field (exact match from transmission_line.csv) 
    - Set query_type to 'transmission'
    - Include flow capacity, voltage levels, and endpoints in text summary
    
    For BUS/SUBSTATION queries (about electrical buses, substations, voltage nodes):
    - Return a list with 'bus name' field (exact match from bus.csv)
    - Set query_type to 'bus' 
    - Include voltage levels, location, and connections in text summary
    
    For GENERAL queries (statistics, summaries, counts):
    - Set query_type to 'general'
    - Provide comprehensive text summary with key statistics
    
    Data Structure Knowledge:
    - bus.csv: Electrical buses/substations with voltage levels (kilovolt), coordinates (WKT), and power flow data
    - generation.csv: Power generation facilities with capacity, output, generation type (hydro/wind/solar/coal/gas/nuclear), and voltage levels
    - transmission_line.csv: Power transmission lines with flow capacity, actual power flow (pf/qf/pt/qt), voltage ratings, and endpoints
    - counties.csv: Geographic boundaries with WKT polygons, FIPS codes, and census areas
    - us states.csv: State boundaries with WKT polygons and geographic identifiers
    - WECC_BA_*.csv: Western Electricity Coordinating Council Balancing Authority data with hourly generation by source type
    
    For geographic queries, use WKT (Well-Known Text) format for spatial data and coordinates
    For power analysis, consider voltage levels, power flows (MW), generation types, and balancing authority areas
    When analyzing power data, distinguish between generation capacity vs actual output
    For temporal analysis, note that WECC data includes hourly time series with generation mix and demand patterns
    
    Always ensure the 'name' fields in your response exactly match the names from the CSV files for proper visualization mapping.
    """
]

agent = Agent(
    model=OpenAIChat(id="gpt-4.1-mini", temperature=0, api_key=openai_key),
    tools=[CsvTools(csvs=csv_list)],
    response_model=PowerGridResponse,
    markdown=True,
    show_tool_calls=True,
    instructions=instruction,
    add_history_to_messages=True,
    num_history_responses=10,
    add_datetime_to_instructions=True
)

def run_agent(query: str) -> PowerGridResponse:
    """
    Run the agent with a query and return structured response for visualization
    
    Args:
        query: Natural language query about power grid data
        
    Returns:
        PowerGridResponse: Structured response with result_list, text, and query_type
    """
    response = agent.run(query)
    # Extract the actual PowerGridResponse from the RunResponse wrapper
    return response.content

def run_agent_json(query: str) -> dict:
    """
    Run the agent and return JSON-serializable dictionary for API responses
    
    Args:
        query: Natural language query about power grid data
        
    Returns:
        dict: JSON-serializable response matching frontend expectations
    """
    response = run_agent(query)
    return response.model_dump(by_alias=True)

if __name__ == "__main__":
    print("Power Grid Data Analysis Assistant")
    print("Available query types:")
    print("- Generation: 'show wind farms', 'find solar plants'")
    print("- Transmission: 'show transmission lines', 'find high voltage lines'")
    print("- Bus/Substations: 'show substations', 'find electrical buses'")
    print("- General: 'statistics on power generation', 'count of facilities'\n")
    
    while True:
        user_input = input("Enter your query (or 'exit' to quit): ")
        if user_input.lower() == 'exit':
            break
        try:
            response = run_agent(user_input)
            print(f"\nQuery Type: {response.query_type}")
            print(f"Results Count: {len(response.result_list)}")
            print(f"Summary: {response.text}")
            
            if response.result_list:
                print("\nDetailed Results:")
                for i, item in enumerate(response.result_list[:5], 1):  # Show first 5 results
                    if hasattr(item, 'generation_name'):
                        print(f"  {i}. Generation: {item.generation_name}")
                    elif hasattr(item, 'line_name'):
                        print(f"  {i}. Transmission Line: {item.line_name}")
                    elif hasattr(item, 'bus_name'):
                        print(f"  {i}. Bus/Substation: {item.bus_name}")
                if len(response.result_list) > 5:
                    print(f"  ... and {len(response.result_list) - 5} more")
            
            # Show JSON format for API integration
            print(f"\nJSON Response for API:")
            import json
            json_response = run_agent_json(user_input)
            print(json.dumps(json_response, indent=2))
            print("-" * 70)
        except Exception as e:
            print(f"Error: {e}")
            print("-" * 70)