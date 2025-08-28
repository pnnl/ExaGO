# ChatGrid Visualization Flow Documentation

## Overview
This document explains how the Westmap application handles chat input through the ChatGrid widget and converts backend responses into dynamic map visualizations with animations.

## Architecture Flow

```
User Input → ChatGrid Widget → Frontend Processing → Backend API → SQL/AI Processing → Response → Visualization Update → Map Animation
```

## Key Components and Files

### 1. Frontend Chat Interface
**File**: `app.js` (lines 770-850)

#### Chat Widget Setup
```javascript
<Widget
  handleNewUserMessage={handleUserInput}
  title="Westmap Chatbot"
  subtitle="What do you want to know about this power grid network?"
/>
```

#### Key Functions:
- **`handleUserInput(inputText)`**: Main handler for chat messages
- **`addResponseMessage()`**: Displays response in chat
- **`toggleMsgLoader()`**: Shows/hides loading indicator

### 2. Backend API Integration
**Files**: 
- `backend/server.py` - Flask server with CORS enabled
- `backend/sqlchain.py` - LangChain SQL agent processing

#### API Endpoint
```bash
POST /api/data
Content-Type: application/json
Body: {"inputText": "user query"}
```

#### Response Format
```json
{
    "result_list": [
        {"generation_name": "ROOSEVELT 1"},
        {"generation_name": "CLE ELUM 2"},
        // ... more results
    ],
    "text": "ROOSEVELT 1, CLE ELUM 2, GOLDENDALE 2, ELLENSBURG 3, DAYTON 2"
}
```

### 3. Data Processing and Visualization
**Files**:
- `src/dataprocess.js` - Core data processing functions
- `module_casedata.js` - Data loading and structure
- `src/color.js` - Color mapping functions

## Detailed Flow Analysis

### Step 1: User Input Processing
```javascript
const handleUserInput = (inputText) => {
    console.log(`New message incoming! ${inputText}`);
    toggleMsgLoader(); // Show loading
    
    const postData = {"inputText": inputText};
    const apiPath = '/api/data';
    
    fetch(apiPath, {
        "method": "POST",
        "body": JSON.stringify(postData),
    })
    .then(res => res.json())
    .then(chatOutput => {
        // Process response...
    });
};
```

### Step 2: Backend Processing
The backend uses LangChain to:
1. Parse natural language queries
2. Convert to SQL queries
3. Execute against PostgreSQL database
4. Return structured results

**Key Backend Function** (`sqlchain.py`):
```python
def sqlchain(input_text):
    llm = ChatOpenAI(model_name="gpt-3.5-turbo", temperature=0)
    db = SQLDatabase.from_uri(postgresql_uri)
    
    text_chain = SQLDatabaseChain.from_llm(llm, db, verbose=True)
    output = text_chain.invoke({"query": input_text})
    
    return {
        "text": output['result'],
        "result_list": query_dict
    }
```

### Step 3: Response Classification and Visualization
The frontend classifies responses based on data structure and triggers appropriate visualizations:

#### Generation Data Visualization
```javascript
if ("generation name" in chatList[0]) {
    const genNameList = chatList.map(d => d["generation name"]);
    setGenLayerActive(true);           // Enable generation layer
    setNameSelectItems(genNameList);   // Filter to specific generators
    setGenFilterValue([gendata.minPg, gendata.maxPg]);
    
    // Animate camera
    setInitialViewState(viewState => ({
        ...viewState,
        pitch: 40,
        transitionInterpolator: transitionFlyToInterpolator,
        transitionDuration: 2000,
    }));
}
```

#### Transmission Line Visualization
```javascript
if ("line name" in chatList[0]) {
    const lineNameList = chatList.map(d => d["line name"]);
    setNetLayerActive(true);
    setFlowLayerActive(true);
    setLineNameSelectItems(lineNameList);
}
```

#### Bus/Substation Visualization
```javascript
if ('bus name' in chatList[0]) {
    const busNameList = chatList.map(d => d["bus name"]);
    setNetLayerActive(true);
    setFlowLayerActive(true);
    setBusNameSelectItems(busNameList);
}
```

### Step 4: Map Layer Updates and Animations

#### Layer Configuration
The application uses Deck.GL layers for visualization:

1. **GeoJsonLayer** - Network topology (buses, lines)
2. **ColumnLayer** - Generation facilities (3D columns)
3. **FlowmapLayer** - Power flow visualization
4. **PolygonLayer** - Geographic regions (counties, areas, zones)

#### Animation System
```javascript
const transitionFlyToInterpolator = new FlyToInterpolator(['zoom']);

const zoomToGen = useCallback((lat, long) => {
    setInitialViewState(viewState => ({
        ...viewState,
        latitude: lat,
        longitude: long,
        pitch: 50,
        transitionInterpolator: transitionFlyToInterpolator,
        transitionDuration: 2000,
        zoom: 7.5,
        onTransitionEnd: activatePopup
    }));
});
```

### Step 5: Real-time Data Filtering
The application uses React's `useEffect` to update visualizations based on chat responses:

```javascript
useEffect(() => {
    const locations = [];
    const flows = [];
    
    if (lineNameSelectItems.length > 0) {
        // Process transmission line data
        data.features.forEach(feature => {
            if (feature.geometry.type === "LineString" && 
                lineNameSelectItems.includes(feature.properties.NAME)) {
                // Create flow visualization data
            }
        });
    }
    
    setFlowData({locations, flows, maxloading: 120});
}, [data, netfiltervalue, flowfiltervalue, lineNameSelectItems]);
```

## Supported Query Types and Responses

### 1. Generation Queries
**Example Input**: "tell me wind generation of 5 locations"

**Response Structure**:
```json
{
    "result_list": [
        {"generation name": "ROOSEVELT 1"},
        {"generation name": "CLE ELUM 2"}
    ],
    "text": "ROOSEVELT 1, CLE ELUM 2, ..."
}
```

**Visualization**:
- Activates generation layer (3D columns)
- Filters to specific generators
- Animates camera to show generation facilities
- Color-codes by fuel type (wind=green, solar=yellow, etc.)

### 2. Transmission Line Queries
**Example Input**: "show transmission lines between Seattle and Portland"

**Response Structure**:
```json
{
    "result_list": [
        {"line name": "SEATTLE -- PORTLAND"},
        {"line name": "SEATTLE 2 -- PORTLAND 3"}
    ],
    "text": "Lines connecting Seattle and Portland"
}
```

**Visualization**:
- Activates network layer
- Highlights specific transmission lines
- Shows power flow direction and magnitude
- Animates camera to focus on selected lines

### 3. Substation/Bus Queries
**Example Input**: "find all substations in California"

**Response Structure**:
```json
{
    "result_list": [
        {"bus name": "LOS ANGELES"},
        {"bus name": "SAN FRANCISCO"}
    ],
    "text": "California substations"
}
```

**Visualization**:
- Activates network layer
- Highlights specific substations
- Shows voltage levels and connections
- Focuses camera on selected region

### 4. Capacity Queries
**Example Input**: "show generation capacity for wind farms"

**Response Detection**: Checks if response contains "capacity" keyword

**Visualization**:
- Activates both generation power and capacity layers
- Shows dual 3D columns (power vs capacity)
- Enables capacity-specific filtering

## Animation Parameters

### Camera Transitions
- **Duration**: 2000ms (2 seconds)
- **Pitch**: 40-50 degrees for 3D visualization
- **Zoom**: 7.5 for focused view
- **Interpolator**: `FlyToInterpolator` for smooth transitions

### Layer Animations
- **Column Height**: Scaled by generation power (factor of 5)
- **Color Transitions**: Smooth interpolation between states
- **Flow Animation**: Real-time power flow visualization

### Filter Animations
- **Range Updates**: Smooth slider transitions
- **Selection Highlighting**: Immediate visual feedback
- **Data Filtering**: Real-time updates using DataFilterExtension

## Data Flow Dependencies

### Generation Data Processing
**Function**: `getGeneration(data)` in `src/dataprocess.js`
- Processes power plant data from GeoJSON
- Calculates fuel mix statistics
- Creates 3D column visualization data
- Color-codes by fuel type

### Network Data Processing
**Function**: `ExtractFlowData(data)` in `src/dataprocess.js`
- Extracts transmission line flows
- Calculates loading percentages
- Creates flow visualization data
- Handles bidirectional flows

### Geographic Data Processing
**Function**: `getCountyNodes(data)` in `src/dataprocess.js`
- Maps power system data to counties
- Aggregates load and generation by region
- Creates polygon visualization data
- Enables spatial filtering

## Configuration and Customization

### Map Styles
Multiple background map options available:
- OpenStreetMap (default)
- Satellite view
- Terrain view
- Positron light/dark themes

### Color Schemes
Fuel type color mapping:
```javascript
const colorMap = {
    'green': 'Wind',
    'yellow': 'Solar', 
    'gray': 'Coal',
    'red': 'Nuclear',
    'blue': 'Hydro',
    'orange': 'Natural Gas',
    'black': 'Other'
};
```

### Filter Ranges
- **Voltage**: 0-800 kV
- **Generation**: Dynamic based on data
- **Power Flow**: 0-120% loading
- **Load**: Dynamic based on county data

## Error Handling

### Backend Errors
```python
except Exception as error:
    if "maximum length" in str(error):
        text_result = "Check the visualization for updated results"
    elif "Rate limit" in str(error):
        text_result = "Check the visualization for updated results"  
    else:
        text_result = "Sorry, I can't find the answer to your question"
```

### Frontend Errors
```javascript
catch (error) {
    setOutputMes("Sorry I didn't find the answer to your question. Please try to rephrase it or provide more details.");
}
```

## Performance Optimizations

### Data Filtering
- Uses Deck.GL's `DataFilterExtension` for efficient filtering
- Updates only necessary layers based on query type
- Maintains filter state across interactions

### Memory Management
- Lazy loading of map data
- Efficient GeoJSON processing
- Minimal re-renders through React optimization

### Network Optimization
- Single API endpoint for all queries
- Compressed JSON responses
- Progressive data loading

## Future Enhancements

### Potential Improvements
1. **WebSocket Integration**: Real-time data updates
2. **Advanced Analytics**: Time-series visualization
3. **Multi-layer Queries**: Complex spatial-temporal queries
4. **Export Functionality**: Save visualizations and data
5. **Voice Interface**: Speech-to-text query input

### Scalability Considerations
1. **Data Pagination**: Handle large result sets
2. **Cluster Visualization**: Aggregate similar features
3. **Performance Monitoring**: Track query response times
4. **Caching Strategy**: Cache frequent queries

This documentation provides a comprehensive overview of how the ChatGrid widget integrates with the visualization system to create dynamic, animated map responses based on natural language queries about power grid data.
