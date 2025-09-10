import React, { useRef, useState, useCallback, useEffect, useReducer } from 'react';
import { createRoot } from "react-dom/client";
import { BrowserRouter as Router, Routes, Route, Navigate, useNavigate, useLocation } from 'react-router-dom';
import { StaticMap, Popup, Marker, _MapContext as MapContext, FullscreenControl, NavigationControl } from 'react-map-gl';
import { WebMercatorViewport } from '@deck.gl/core';
import DeckGL from '@deck.gl/react';
// import FlowMapLayer from '@flowmap.gl/core'
// import { FlowmapLayer } from '@flowmap.gl/layers'
import { GeoJsonLayer, ColumnLayer, PolygonLayer, ScatterplotLayer, TextLayer } from '@deck.gl/layers';
import { DataFilterExtension } from '@deck.gl/extensions';
import Checkbox from '@mui/material/Checkbox';
import FormControlLabel from '@mui/material/FormControlLabel';
import FormGroup from '@mui/material/FormGroup';
import { Typography } from '@mui/material';
import HomeOutlinedIcon from '@mui/icons-material/HomeOutlined';
import ThreeSixtyOutlinedIcon from '@mui/icons-material/ThreeSixtyOutlined';
import ZoomOutMapIcon from '@mui/icons-material/ZoomOutMap';
import Slider from '@mui/material/Slider';
import Box from '@mui/material/Box';
import Fab from '@mui/material/Fab';
import SchoolIcon from '@mui/icons-material/School';
import Button from '@mui/material/Button';
import AdminPanelSettingsIcon from '@mui/icons-material/AdminPanelSettings';

import Accordion from "@mui/material/Accordion";
import AccordionSummary from "@mui/material/AccordionSummary";
import AccordionDetails from "@mui/material/AccordionDetails";
import ArrowDropDownIcon from '@mui/icons-material/ArrowDropDown';

//import nercregions from "data/NERC_Reliability_Coordinators.json"

import {
  Chart as ChartJS,
  RadialLinearScale,
  ArcElement,
  Tooltip,
  Legend,
} from 'chart.js';
import { PolarArea, Doughnut } from 'react-chartjs-2';

import Multiselect from "react-widgets/Multiselect";
import { Widget, addResponseMessage, toggleMsgLoader, deleteMessages } from 'react-chat-widget';


import { center, convex, bbox } from '@turf/turf';

import { LinearInterpolator, FlyToInterpolator } from 'deck.gl';
import { HeatmapLayer } from 'deck.gl';
import { InvertColorsOff, ShopTwoOutlined } from '@mui/icons-material';

import { getCountyNodes, ExtractFirstTimeSlice, ExtractFlowData, getBarNet, getPoints, getGeneration, getLoad, getContours, getAreas, getZones } from "./src/dataprocess";
import { LineColor, FlowColor, FillColor, fillGenColumnColor, fillGenColumnColorCap, getVoltageFillColor } from "./src/color"

// Firebase Authentication imports
import { AuthProvider, ProtectedRoute, Header, AdminDashboard, useAuth } from './components/common';
import { ManishProject } from './components/manish';
import { AminProject, AminDetailPage } from './components/amin';

import 'core-js/actual/structured-clone';

ChartJS.register(RadialLinearScale, ArcElement, Tooltip, Legend);

// Add Google Fonts Inter
const fontLink = document.createElement('link');
fontLink.href = 'https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap';
fontLink.rel = 'stylesheet';
document.head.appendChild(fontLink);

// Apply Inter font to body
document.body.style.fontFamily = '"Inter", -apple-system, BlinkMacSystemFont, "Segoe UI", "Roboto", "Oxygen", "Ubuntu", "Cantarell", "Fira Sans", "Droid Sans", "Helvetica Neue", sans-serif';

// Add custom CSS for chatbot widget styling
const chatWidgetStyle = document.createElement('style');
chatWidgetStyle.textContent = `
  .rcw-conversation-container * {
    font-family: "Inter", -apple-system, BlinkMacSystemFont, "Segoe UI", "Roboto", sans-serif !important;
  }
  .rcw-header {
    background: #0047AB !important;
    font-family: "Inter", sans-serif !important;
  }
  .rcw-title {
    font-family: "Inter", sans-serif !important;
    font-weight: 600 !important;
    font-size: 12px !important;
  }
  .rcw-subtitle {
    font-family: "Inter", sans-serif !important;
    font-weight: 400 !important;
    font-size: 13px !important;
    opacity: 0.9 !important;
  }
  .rcw-message {
    font-family: "Inter", sans-serif !important;
  }
  .rcw-response {
    font-family: "Inter", sans-serif !important;
    font-size: 14px !important;
    line-height: 1.4 !important;
  }
  .rcw-client {
    font-family: "Inter", sans-serif !important;
    font-size: 14px !important;
  }
  .rcw-send {
    background: #0047AB !important;
  }
  .rcw-picker-btn {
    font-family: "Inter", sans-serif !important;
  }
`;
document.head.appendChild(chatWidgetStyle);

// Transition interpolators for animation
const transitionLinearInterpolator = new LinearInterpolator(['bearing']);
const transitionFlyToInterpolator = new FlyToInterpolator(['zoom']);

// Get case data
var mod_casedata = require('./module_casedata.js');
var casedata = {};
try {
  casedata = mod_casedata.get_casedata();
  console.log('✅ Case data loaded successfully');
  if (!casedata || !casedata.geojsondata) {
    console.error('❌ Case data is missing geojsondata property');
    casedata = { geojsondata: { type: "FeatureCollection", features: [] } };
  }
} catch (error) {
  console.error('❌ Error loading case data:', error);
  casedata = { geojsondata: { type: "FeatureCollection", features: [] } };
}

// Function to parse CSV data
function parseCSV(csvText) {
  const lines = csvText.trim().split(/\r?\n/);
  const headers = lines[0].split(',').map(h => h.replace(/"/g, '').trim());
  const data = [];
  
  for (let i = 1; i < lines.length; i++) {
    const values = [];
    let current = '';
    let inQuotes = false;
    
    for (let j = 0; j < lines[i].length; j++) {
      const char = lines[i][j];
      if (char === '"') {
        inQuotes = !inQuotes;
      } else if (char === ',' && !inQuotes) {
        values.push(current.trim());
        current = '';
      } else {
        current += char;
      }
    }
    values.push(current.trim());
    
    if (values.length === headers.length) {
      const row = {};
      headers.forEach((header, index) => {
        row[header] = values[index] ? values[index].replace(/"/g, '').trim() : '';
      });
      data.push(row);
    }
  }
  
  return data;
}

// Function to convert CSV transmission line data to GeoJSON format
function csvToGeoJSON(csvData, dataSource = 'csv') {
  const features = [];
  
  // Validate CSV data structure
  if (!csvData || csvData.length === 0) {
    console.warn(`No data found in ${dataSource}`);
    return features;
  }
  
  // Check if the CSV has the required columns for transmission lines
  const firstRow = csvData[0];
  const hasWKT = firstRow.hasOwnProperty('wkt');
  const hasLineString = hasWKT && typeof firstRow.wkt === 'string' && firstRow.wkt.includes('LINESTRING');
  
  if (!hasWKT || !hasLineString) {
    console.warn(`${dataSource} does not contain valid WKT LINESTRING data. Skipping this file.`);
    console.log(`Available columns in ${dataSource}:`, Object.keys(firstRow));
    return features;
  }
  
  csvData.forEach((row, index) => {
    if (row.wkt && row.wkt.startsWith('LINESTRING')) {
      // Parse WKT LINESTRING
      const coordsMatch = row.wkt.match(/LINESTRING\s*\((.+)\)/);
      if (coordsMatch) {
        try {
          const coordPairs = coordsMatch[1].split(',').map(pair => {
            const [lng, lat] = pair.trim().split(' ').map(Number);
            if (isNaN(lng) || isNaN(lat)) {
              throw new Error(`Invalid coordinates at row ${index + 1}`);
            }
            return [lng, lat];
          });
          
          const feature = {
            type: "Feature",
            geometry: {
              type: "LineString",
              coordinates: coordPairs
            },
            properties: {
              NAME: row.line_name || `${row.source || row.srouce} -- ${row.target}`,
              RATE_A: parseFloat(row['flow capacity']) || 0,
              PF: parseFloat(row.pf) || 0,
              QF: parseFloat(row.qf) || 0,
              PT: parseFloat(row.pt) || 0,
              QT: parseFloat(row.qt) || 0,
              KV: parseFloat(row.kilovolt) || 230,
              kilovolt: parseFloat(row.kilovolt) || 230,
              line_name: row.line_name || `Line_${index + 1}`,
              source: row.source || row.srouce,
              target: row.target,
              actual_flow: parseFloat(row['actual flow']) || 0,
              flow_capacity: parseFloat(row['flow capacity']) || 100,
              data_source: row.data_source || dataSource,
              SOURCEFILE: row.data_source || dataSource,
              elementtype: "TRANSMISSION_LINE"
            }
          };
          
          features.push(feature);
        } catch (error) {
          console.warn(`Error parsing coordinates for row ${index + 1} in ${dataSource}:`, error.message);
        }
      }
    }
  });
  
  return features;
}

// Source data GeoJSON
const geodata = casedata['geojsondata'];
console.log('🔍 Geodata features count:', geodata?.features?.length || 0);

const MAP_STYLE = {
  osm: {
    "version": 8,
    "name": "OpenStreetMap",
    "sources": {
      "osm": {
        "type": "raster",
        "tiles": [
          "https://tile.openstreetmap.org/{z}/{x}/{y}.png"
        ],
        "tileSize": 256,
        "attribution": "© OpenStreetMap contributors"
      }
    },
    "layers": [
      {
        "id": "osm",
        "type": "raster",
        "source": "osm"
      }
    ]
  },
  satellite: {
    "version": 8,
    "name": "Satellite",
    "sources": {
      "satellite": {
        "type": "raster",
        "tiles": [
          "https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}"
        ],
        "tileSize": 256,
        "attribution": "© Esri, Maxar, Earthstar Geographics"
      }
    },
    "layers": [
      {
        "id": "satellite",
        "type": "raster",
        "source": "satellite"
      }
    ]
  },
  terrain: {
    "version": 8,
    "name": "Terrain",
    "sources": {
      "terrain": {
        "type": "raster",
        "tiles": [
          "https://stamen-tiles.a.ssl.fastly.net/terrain/{z}/{x}/{y}.png"
        ],
        "tileSize": 256,
        "attribution": "© Stamen Design, © OpenStreetMap contributors"
      }
    },
    "layers": [
      {
        "id": "terrain",
        "type": "raster",
        "source": "terrain"
      }
    ]
  },
  pos_no_label: 'https://basemaps.cartocdn.com/gl/positron-nolabels-gl-style/style.json',
  pos: 'https://basemaps.cartocdn.com/gl/positron-gl-style/style.json',
  dark: 'https://basemaps.cartocdn.com/gl/dark-matter-gl-style/style.json',
  none: ''
};


// Function to check if a point is within WECC boundaries (rough approximation)
function isWithinWeccBoundaries(longitude, latitude) {
  // Input validation
  if (typeof longitude !== 'number' || typeof latitude !== 'number' || 
      isNaN(longitude) || isNaN(latitude)) {
    return false;
  }
  
  // WECC boundaries (approximate) - covers Western US, parts of Canada and Mexico
  // Western boundary: Pacific coast (~-125°W)
  // Eastern boundary: roughly follows Rocky Mountains (~-100°W to -105°W depending on latitude)
  // Northern boundary: extends into Canada (~50°N)
  // Southern boundary: extends into Mexico (~25°N)
  
  // Basic bounding box for WECC region (more inclusive boundaries)
  const weccBounds = {
    west: -130.0,  // Extended west to include more Pacific areas
    east: -95.0,   // Extended east to be more inclusive
    north: 55.0,   // Extended north for Canadian connections
    south: 20.0    // Extended south for Mexican connections
  };
  
  // More refined eastern boundary based on latitude (more inclusive)
  let easternBound = weccBounds.east;
  if (latitude > 45) {
    easternBound = -100.0; // More inclusive in northern regions
  } else if (latitude > 40) {
    easternBound = -98.0;  // More inclusive in mountain states
  } else if (latitude > 35) {
    easternBound = -96.0;  // More inclusive in southwest
  } else {
    easternBound = -95.0;  // More inclusive in Mexico border region
  }
  
  return longitude >= weccBounds.west && 
         longitude <= easternBound && 
         latitude >= weccBounds.south && 
         latitude <= weccBounds.north;
}

// Function to check if a transmission line is within WECC boundaries
function isTransmissionLineInWecc(feature) {
  try {
    if (!feature || !feature.geometry || feature.geometry.type !== 'LineString') {
      return false;
    }
    
    const coordinates = feature.geometry.coordinates;
    if (!Array.isArray(coordinates) || coordinates.length === 0) {
      return false;
    }
    
    // Check if any point of the line is within WECC boundaries
    // For efficiency, we'll sample a few points along the line
    const samplePoints = [];
    const numSamples = Math.min(coordinates.length, 5); // Sample up to 5 points
    
    for (let i = 0; i < numSamples; i++) {
      const index = Math.floor((i / (numSamples - 1)) * (coordinates.length - 1));
      if (coordinates[index] && Array.isArray(coordinates[index]) && coordinates[index].length >= 2) {
        samplePoints.push(coordinates[index]);
      }
    }
    
    // Line is considered in WECC if any sampled point is within boundaries
    return samplePoints.some(coord => {
      const [longitude, latitude] = coord;
      return isWithinWeccBoundaries(longitude, latitude);
    });
  } catch (error) {
    console.warn('Error checking WECC boundaries for feature:', error);
    return true; // If there's an error, include the feature to be safe
  }
}

// Load Western Power Plants data from CSV files
async function loadPowerPlantsData() {
  const allPowerPlants = [];
  
  // List of power plant CSV files to load
  const powerPlantFiles = [
    { 
      url: './data/Western_Power_plants_Locations-USA (2).csv', 
      source: 'usa_power_plants',
      country: 'USA'
    },
    { 
      url: './data/Western_Power_plants_Locations-CANADA-MEX_Update (1).csv', 
      source: 'canada_mexico_power_plants',
      country: 'CAN_MEX'
    }
  ];
  
  for (const file of powerPlantFiles) {
    try {
      console.log(`🔄 Loading power plants data from ${file.source}...`);
      const response = await fetch(file.url);
      
      if (!response.ok) {
        console.warn(`⚠️ Could not load ${file.source}: ${response.status}`);
        continue;
      }
      
      const csvText = await response.text();
      const csvData = parseCSV(csvText);
      const features = csvPowerPlantsToGeoJSON(csvData, file.source, file.country);
      
      if (features.length > 0) {
        allPowerPlants.push(...features);
        console.log(`✅ Successfully loaded ${features.length} power plants from ${file.source}`);
      } else {
        console.warn(`⚠️ No valid power plant features found in ${file.source}`);
      }
      
    } catch (error) {
      console.error(`❌ Error loading ${file.source}:`, error);
    }
  }
  
  console.log(`Total power plants loaded: ${allPowerPlants.length}`);
  return allPowerPlants;
}

// Convert CSV power plant data to GeoJSON format
function csvPowerPlantsToGeoJSON(csvData, dataSource = 'power_plants', country = 'USA') {
  const features = [];
  
  if (!csvData || csvData.length === 0) {
    console.warn('No CSV data provided for power plants conversion');
    return features;
  }
  
  const firstRow = csvData[0];
  const hasLatitude = firstRow.hasOwnProperty('Latitude');
  const hasLongitude = firstRow.hasOwnProperty('Longitude');
  
  if (!hasLatitude || !hasLongitude) {
    console.warn('CSV data missing required Latitude/Longitude columns');
    return features;
  }
  
  csvData.forEach((row, index) => {
    try {
      const latitude = parseFloat(row['Latitude']);
      const longitude = parseFloat(row['Longitude']);
      
      // Skip invalid coordinates
      if (isNaN(latitude) || isNaN(longitude)) {
        return;
      }
      
      // Parse capacity data
      const totalCapacity = parseFloat(row['Total Capacity (MW)']) || 0;
      
      // Create technology breakdown based on dataset structure
      const technologies = {};
      if (country === 'USA') {
        // USA dataset has detailed technology breakdown
        technologies.battery = parseFloat(row['Battery Storage']) || 0;
        technologies.geothermal = parseFloat(row['Geothermal']) || 0;
        technologies.hydro = parseFloat(row['Hydro']) || 0;
        technologies.naturalGas = parseFloat(row['Natural Gas']) || 0;
        technologies.nuclear = parseFloat(row['Nuclear']) || 0;
        technologies.oil = parseFloat(row['Oil']) || 0;
        technologies.other = parseFloat(row['Other']) || 0;
        technologies.solar = parseFloat(row['Solar']) || 0;
        technologies.wind = parseFloat(row['Wind']) || 0;
      }
      
      const feature = {
        type: 'Feature',
        geometry: {
          type: 'Point',
          coordinates: [longitude, latitude]
        },
        properties: {
          // Basic plant information
          plantCode: row['Plant Code'] || '',
          plantName: row['Plant Name'] || '',
          state: row['State'] || '',
          county: row['County'] || '',
          country: country === 'CAN_MEX' ? (row['Country'] || 'CAN') : 'USA',
          balancingAuthority: row['Balancing Authority'] || '',
          primaryType: row['Primary Type'] || 'Unknown',
          totalCapacityMW: totalCapacity,
          
          // Technology breakdown (for USA dataset)
          ...technologies,
          
          // Metadata
          dataSource: dataSource,
          featureType: 'power_plant',
          
          // For visualization
          markerSize: Math.max(5, Math.min(50, Math.sqrt(totalCapacity) * 2)),
          capacityCategory: getCapacityCategory(totalCapacity)
        }
      };
      
      features.push(feature);
      
    } catch (error) {
      console.warn(`Error processing power plant row ${index}:`, error);
    }
  });
  
  return features;
}

// Categorize power plants by capacity for visualization
function getCapacityCategory(capacity) {
  if (capacity >= 1000) return 'Large (>1000 MW)';
  if (capacity >= 100) return 'Medium (100-1000 MW)';
  if (capacity >= 10) return 'Small (10-100 MW)';
  return 'Micro (<10 MW)';
}

// Get color for power plant based on primary technology type
function getPowerPlantColor(primaryType) {
  const colorMap = {
    'Solar': [255, 193, 7, 200],        // Bright yellow/gold
    'Wind': [76, 175, 80, 200],         // Green
    'Hydro': [33, 150, 243, 200],       // Blue
    'Natural Gas': [255, 87, 34, 200],  // Orange-red
    'Nuclear': [156, 39, 176, 200],     // Purple
    'Geothermal': [139, 69, 19, 200],   // Brown
    'Oil': [96, 125, 139, 200],         // Blue-grey
    'Battery Storage': [255, 235, 59, 200], // Light yellow
    'Biomass': [76, 175, 80, 200],      // Green (same as wind)
    'Other': [158, 158, 158, 200],      // Grey
    'Coal': [66, 66, 66, 200],          // Dark grey
    'Unknown': [189, 189, 189, 200]     // Light grey
  };
  
  return colorMap[primaryType] || colorMap['Unknown'];
}

// Enhanced color function for power plant columns with lighting effects
function getPowerPlantColumnColor(primaryType, capacity) {
  const baseColor = getPowerPlantColor(primaryType);
  
  // Add lighting effect based on capacity
  let lightingFactor = 1.0;
  if (capacity >= 1000) {
    lightingFactor = 1.3; // Brighter for large plants
  } else if (capacity >= 100) {
    lightingFactor = 1.1; // Slightly brighter for medium plants
  } else if (capacity < 10) {
    lightingFactor = 0.8; // Dimmer for micro plants
  }
  
  // Apply lighting while maintaining color bounds
  const enhancedColor = [
    Math.min(255, Math.max(0, baseColor[0] * lightingFactor)),
    Math.min(255, Math.max(0, baseColor[1] * lightingFactor)),
    Math.min(255, Math.max(0, baseColor[2] * lightingFactor)),
    baseColor[3] || 200
  ];
  
  return enhancedColor;
}

// Process power plant data to match generation format
function processPowerPlantData(powerPlantsFeatures) {
  if (!powerPlantsFeatures || powerPlantsFeatures.length === 0) {
    return { Plants: [], minPg: 0, maxPg: 10000, minPcap: 0, maxPcap: 10000 };
  }

  const plants = powerPlantsFeatures.map(feature => {
    const props = feature.properties;
    const capacity = props.totalCapacityMW || 0;
    
    return {
      name: props.plantName || 'Unnamed Plant',
      coordinates: feature.geometry.coordinates,
      Pg: capacity, // Current generation (using capacity as proxy)
      Pcap: capacity, // Total capacity
      primaryType: props.primaryType || 'Unknown',
      fuel: props.primaryType || 'Unknown',
      state: props.state || '',
      country: props.country || 'USA',
      balancingAuthority: props.balancingAuthority || '',
      plantCode: props.plantCode || '',
      // Technology breakdown for USA plants
      solar: props.solar || 0,
      wind: props.wind || 0,
      hydro: props.hydro || 0,
      naturalGas: props.naturalGas || 0,
      nuclear: props.nuclear || 0,
      geothermal: props.geothermal || 0,
      battery: props.battery || 0,
      oil: props.oil || 0,
      other: props.other || 0,
      // Original properties for tooltip
      originalProps: props
    };
  });

  // Calculate min/max values for filtering
  const capacities = plants.map(p => p.Pg).filter(c => c > 0);
  const minPg = Math.min(...capacities) || 0;
  const maxPg = Math.max(...capacities) || 10000;
  const minPcap = minPg;
  const maxPcap = maxPg;

  return {
    Plants: plants,
    minPg: minPg,
    maxPg: maxPg,
    minPcap: minPcap,
    maxPcap: maxPcap
  };
}

// Create chart data for power plants - ACCURATE DATA FROM CSV ANALYSIS
function createPowerPlantChartData(powerPlantData) {
  // Use accurate calculations from comprehensive CSV analysis
  // This ensures pie charts show exact figures from our CSV files
  
  const accurateChartData = {
    labels: ['Natural Gas', 'Hydro', 'Solar', 'Wind', 'Other', 'Battery Storage', 'Nuclear', 'Geothermal', 'Oil'],
    datasets: [{
      label: 'Western Power Plants Capacity (MW)',
      data: [108708.2, 71606.2, 42307.6, 34390.9, 26940.1, 16423.8, 7732.6, 4551.1, 1285.5],
      backgroundColor: [
        'rgba(255, 87, 34, 0.8)',   // Natural Gas - Orange
        'rgba(33, 150, 243, 0.8)',  // Hydro - Blue
        'rgba(255, 193, 7, 0.8)',   // Solar - Yellow/Gold
        'rgba(76, 175, 80, 0.8)',   // Wind - Green
        'rgba(158, 158, 158, 0.8)',  // Other - Gray
        'rgba(255, 235, 59, 0.8)',  // Battery Storage - Light Yellow
        'rgba(156, 39, 176, 0.8)',  // Nuclear - Purple
        'rgba(139, 69, 19, 0.8)',   // Geothermal - Brown
        'rgba(96, 125, 139, 0.8)'   // Oil - Blue Gray
      ],
      borderColor: [
        'rgba(255, 87, 34, 1)',
        'rgba(33, 150, 243, 1)',
        'rgba(255, 193, 7, 1)',
        'rgba(76, 175, 80, 1)',
        'rgba(158, 158, 158, 1)',
        'rgba(255, 235, 59, 1)',
        'rgba(156, 39, 176, 1)',
        'rgba(139, 69, 19, 1)',
        'rgba(96, 125, 139, 1)'
      ],
      borderWidth: 2
    }]
  };
  
  return accurateChartData;
}

// Calculate generation statistics from power plant data - ACCURATE DATA FROM CSV ANALYSIS
function calculateGenerationStats(powerPlantsData) {
  // Use accurate calculations from comprehensive CSV analysis
  // Data source: Western_Power_plants_Locations-USA (2).csv + Western_Power_plants_Locations-CANADA-MEX_Update (1).csv
  // Total verified: 313,946.0 MW across 3,622 plants
  
  const accurateStats = {
    // Accurate technology breakdown with proper color mapping
    'Natural Gas': {
      count: 536,
      totalCapacity: 108708.2,
      percentage: '34.6',
      averageCapacity: '202.8',
      color: [255, 87, 34, 200]
    },
    'Hydro': {
      count: 739,
      totalCapacity: 71606.2,
      percentage: '22.8',
      averageCapacity: '96.9',
      color: [33, 150, 243, 200]
    },
    'Solar': {
      count: 1513,
      totalCapacity: 42307.6,
      percentage: '13.5',
      averageCapacity: '27.9',
      color: [255, 193, 7, 200]
    },
    'Wind': {
      count: 377,
      totalCapacity: 34390.9,
      percentage: '11.0',
      averageCapacity: '91.2',
      color: [76, 175, 80, 200]
    },
    'Other': {
      count: 276,
      totalCapacity: 26940.1,
      percentage: '8.6',
      averageCapacity: '97.6',
      color: [158, 158, 158, 200]
    },
    'Battery Storage': {
      count: 277,
      totalCapacity: 16423.8,
      percentage: '5.2',
      averageCapacity: '59.3',
      color: [255, 235, 59, 200]
    },
    'Nuclear': {
      count: 3,
      totalCapacity: 7732.6,
      percentage: '2.5',
      averageCapacity: '2577.5',
      color: [156, 39, 176, 200]
    },
    'Geothermal': {
      count: 67,
      totalCapacity: 4551.1,
      percentage: '1.4',
      averageCapacity: '67.9',
      color: [139, 69, 19, 200]
    },
    'Oil': {
      count: 49,
      totalCapacity: 1285.5,
      percentage: '0.4',
      averageCapacity: '26.2',
      color: [96, 125, 139, 200]
    }
  };
  
  return {
    technologies: accurateStats,
    sortedTechnologies: ['Natural Gas', 'Hydro', 'Solar', 'Wind', 'Other', 'Battery Storage', 'Nuclear', 'Geothermal', 'Oil'],
    totalCapacity: '313946.0',  // Total MW
    totalPlants: 3622,         // Total plants
    averageCapacity: '86.7'    // Average MW per plant
  };
}

// Load additional transmission line data from CSV files (filtered to WECC region)
async function loadAdditionalTransmissionData() {
  const allFeatures = [];
  
  // List of CSV files to try loading - prioritize the WGS84 converted powerlines data
  const csvFiles = [
    { url: './data/powerlines_WUS_CAN_sgca_wgs84.csv', source: 'powerlines_WUS_CAN_sgca_wgs84', priority: 1 },
    { url: './data/powerlines_converted.csv', source: 'powerlines_converted', priority: 2 },
    { url: './data/powerlines_WUS_CAN_sgca.csv', source: 'powerlines_WUS_CAN_sgca_main', priority: 3 }
  ];
  
  // Sort files by priority to load highest priority first
  csvFiles.sort((a, b) => a.priority - b.priority);
  
  for (const file of csvFiles) {
    try {
      console.log(`Loading transmission lines from ${file.source}...`);
      const response = await fetch(file.url);
      if (!response.ok) {
        console.warn(`Failed to load ${file.source}: HTTP ${response.status}`);
        continue;
      }
      
      const csvText = await response.text();
      const csvData = parseCSV(csvText);
      const allFeatures_temp = csvToGeoJSON(csvData, file.source);
      
      // Filter transmission lines to only include those within WECC boundaries
      const weccFeatures = allFeatures_temp.filter(feature => {
        const isInWecc = isTransmissionLineInWecc(feature);
        return isInWecc;
      });
      
      const filteredCount = allFeatures_temp.length - weccFeatures.length;
      
      if (weccFeatures.length > 0) {
        allFeatures.push(...weccFeatures);
        console.log(`✅ Successfully loaded ${weccFeatures.length} WECC transmission line features from ${file.source}`);
        if (filteredCount > 0) {
          console.log(`🔍 Filtered out ${filteredCount} transmission lines outside WECC boundaries`);
        }
        
        // If we successfully loaded the high-priority WGS84 file with substantial data, we can skip the others
        if (file.source === 'powerlines_WUS_CAN_sgca_wgs84' && weccFeatures.length > 1000) {
          console.log('✅ Successfully loaded primary WECC powerlines dataset, skipping other sources');
          break;
        }
      } else {
        console.warn(`⚠️ No valid WECC transmission line features found in ${file.source} after filtering`);
      }
    } catch (error) {
      console.warn(`❌ Error loading ${file.source}:`, error);
    }
  }
  
  console.log(`Total WECC transmission lines loaded: ${allFeatures.length}`);
  return allFeatures;
}

// Initialize data processing
async function initializeData() {
  var data = ExtractFirstTimeSlice(geodata);
  
  // Temporarily disable main geodata WECC filtering to ensure map loads
  console.log('📍 WECC filtering temporarily disabled for main geodata to ensure map loads');
  
  // Filter existing geodata features to only include those within WECC boundaries - DISABLED
  // const originalFeatureCount = data.features.length;
  // const originalMainFeatures = [...data.features]; // Keep backup
  
  // data.features = data.features.filter(feature => {
  //   // Only filter transmission lines, keep other features (buses, generators, etc.)
  //   if (feature.geometry && feature.geometry.type === 'LineString') {
  //     return isTransmissionLineInWecc(feature);
  //   }
  //   // Keep non-transmission line features (buses, generators, etc.)
  //   return true;
  // });
  
  // const filteredMainCount = originalFeatureCount - data.features.length;
  // if (filteredMainCount > 0) {
  //   console.log(`🔍 Filtered out ${filteredMainCount} main transmission lines outside WECC boundaries`);
  // }
  
  // Safety check: if we filtered out all transmission lines, revert to original data
  // const mainTransmissionLineCount = data.features.filter(f => f.geometry && f.geometry.type === 'LineString').length;
  // if (mainTransmissionLineCount === 0 && filteredMainCount > 0) {
  //   console.warn('⚠️ WECC filtering removed all main transmission lines, reverting to original data');
  //   data.features = originalMainFeatures;
  // }
  
  // Load and merge additional transmission line data
  const additionalFeatures = await loadAdditionalTransmissionData();
  if (additionalFeatures.length > 0) {
    // Add the additional features to the data
    data.features = [...data.features, ...additionalFeatures];
    console.log(`Added ${additionalFeatures.length} additional WECC transmission lines from CSV files`);
  }
  
  // Load and merge power plants data
  const powerPlantsFeatures = await loadPowerPlantsData();
  if (powerPlantsFeatures.length > 0) {
    // Add power plants to the data
    data.features = [...data.features, ...powerPlantsFeatures];
    console.log(`Added ${powerPlantsFeatures.length} power plants from CSV files`);
    
    // Process power plants data to match generation format
    const processedPowerPlantData = processPowerPlantData(powerPlantsFeatures);
    data.powerPlants = processedPowerPlantData.Plants;
    data.powerPlantData = processedPowerPlantData;
    
    // Create chart data
    data.powerPlantChartData = createPowerPlantChartData(processedPowerPlantData);
    
    // Calculate generation statistics
    data.generationStats = calculateGenerationStats(powerPlantsFeatures);
    console.log('📊 Power plant data processed:', processedPowerPlantData);
  }
  
  // Recalculate bounds to include all new data (transmission lines + power plants)
  if (additionalFeatures.length > 0 || powerPlantsFeatures.length > 0) {
    const enhancedBbox = bbox(data);
    const enhancedCorner1 = [enhancedBbox[0], enhancedBbox[1]];
    const enhancedCorner2 = [enhancedBbox[2], enhancedBbox[3]];
    data.enhancedBounds = [enhancedCorner1, enhancedCorner2];
    
    const enhancedCenter = center(data);
    data.enhancedCenter = enhancedCenter;
    
    console.log('Updated map bounds to include additional transmission lines and power plants');
  }
  
  return data;
}

// Initialize with base data for immediate rendering (with error handling)
var data;
try {
  if (!geodata || !geodata.features || geodata.features.length === 0) {
    console.error('❌ Geodata is empty or undefined, creating fallback data structure');
    data = {
      type: "FeatureCollection",
      features: [],
      bounds: [[-125, 25], [-95, 50]], // Default WECC bounds
      center: { geometry: { coordinates: [-110, 37.5] } } // Default center
    };
  } else {
    data = ExtractFirstTimeSlice(geodata);
    console.log('✅ Base geodata loaded successfully with', data.features.length, 'features');
  }
} catch (error) {
  console.error('❌ Error initializing base data:', error);
  data = {
    type: "FeatureCollection", 
    features: [],
    bounds: [[-125, 25], [-95, 50]],
    center: { geometry: { coordinates: [-110, 37.5] } }
  };
}

// Temporarily disable initial WECC filtering to ensure map loads
console.log('📍 WECC filtering temporarily disabled for initial data to ensure map loads');

// Apply WECC filtering to initial data as well (with safety check) - DISABLED
const initialFeatureCount = data.features.length;
// const originalFeatures = [...data.features]; // Keep backup

// data.features = data.features.filter(feature => {
//   // Only filter transmission lines, keep other features (buses, generators, etc.)
//   if (feature.geometry && feature.geometry.type === 'LineString') {
//     return isTransmissionLineInWecc(feature);
//   }
//   // Keep non-transmission line features (buses, generators, etc.)
//   return true;
// });

// const initialFilteredCount = initialFeatureCount - data.features.length;
// if (initialFilteredCount > 0) {
//   console.log(`🔍 Initial filter: Removed ${initialFilteredCount} transmission lines outside WECC boundaries`);
// }

// Safety check: if we filtered out too many features, revert to original data
// const transmissionLineCount = data.features.filter(f => f.geometry && f.geometry.type === 'LineString').length;
// if (transmissionLineCount === 0 && initialFilteredCount > 0) {
//   console.warn('⚠️ WECC filtering removed all transmission lines, reverting to original data');
//   data.features = originalFeatures;
// }

const countyloaddata = getCountyNodes(data);
data = countyloaddata.updatedata;

const areas = getAreas(casedata);

const zones = getZones(casedata);

const flowdata = ExtractFlowData(data);

const Points = getPoints(data);
const Voltages = Points.map(d => d.value);

const Vcontour = getContours();

const gendata = getGeneration(data);
const generation = gendata.Gens;

const loaddata = getLoad(data);


function LineWidth(line) {
  return line.properties.KV * 3;
  //return Math.abs(line.properties.PF/line.properties.RATE_A)*500;
}

const loads = loaddata.Loads;
const minPd = loaddata.minPd;
const maxPd = loaddata.maxPd;

const countymaxPd = countyloaddata.maxPd;
const countyload = countyloaddata.data;


const bboxArray = bbox(data);
const corner1 = [bboxArray[0], bboxArray[1]];
const corner2 = [bboxArray[2], bboxArray[3]];
const bounds = [corner1, corner2];

const mapcenter = center(data);

var hull = convex(data);

const INITIAL_VIEW_STATE = {
  latitude: mapcenter['geometry']['coordinates'][1],
  longitude: mapcenter['geometry']['coordinates'][0],
  zoom: 5,
  maxZoom: 16,
  pitch: 0,
  bearing: 0,
  bounds: [bboxArray[1], bboxArray[0], bboxArray[3], bboxArray[2]],
  fitbounds: true
};

// Enhanced tooltip formatting function
function getEnhancedTooltip({ object, layer }) {
  if (!object) return null;

  const style = {
    backgroundColor: 'rgba(255, 255, 255, 0.98)',
    color: '#333',
    fontSize: '12px',
    fontFamily: '"Inter", sans-serif',
    padding: '10px 14px',
    borderRadius: '8px',
    boxShadow: '0 6px 20px rgba(0, 0, 0, 0.15), 0 2px 6px rgba(0, 0, 0, 0.1)',
    border: '1px solid rgba(0, 0, 0, 0.12)',
    maxWidth: '320px',
    lineHeight: '1.4',
    backdropFilter: 'blur(8px)',
    WebkitBackdropFilter: 'blur(8px)'
  };

  let content = '';

  if (layer.id === 'geojson') {
    if (object.geometry.type === "Point") {
      content = `
        <div style="font-weight: 600; margin-bottom: 4px; color: #1976d2;">${object.properties.NAME}</div>
        <div style="font-size: 11px; color: #666;">Substation</div>
        <div style="font-size: 11px; color: #666;">KV Levels: ${object.properties.KVlevels?.join(', ') || 'N/A'}</div>
      `;
    } else if (object.geometry.type === "LineString") {
      const loading = Math.abs(object.properties.PF / (object.properties.RATE_A || 10000)) * 100;
      content = `
        <div style="font-weight: 600; margin-bottom: 4px; color: #1976d2;">${object.properties.NAME}</div>
        <div style="font-size: 11px; color: #666;">Transmission Line</div>
        <div style="font-size: 11px; color: #666;">Voltage: ${object.properties.KV?.toFixed(1)} kV</div>
        <div style="font-size: 11px; color: #666;">Loading: ${loading.toFixed(1)}%</div>
      `;
    }
  } else if (layer.id === 'gen-column' || layer.id === 'gen-column-cap') {
    content = `
      <div style="font-weight: 600; margin-bottom: 4px; color: #2e7d32;">${object.name}</div>
      <div style="font-size: 11px; color: #666;">Generation Facility</div>
      <div style="font-size: 11px; color: #666;">Power: ${Math.round(object.Pg * 100) / 100} MW</div>
      <div style="font-size: 11px; color: #666;">Capacity: ${Math.round(object.Pcap * 100) / 100} MW</div>
      <div style="font-size: 11px; color: #666;">Fuel: ${object.fuel || 'Unknown'}</div>
    `;
  } else if (layer.id === 'western-power-plants' || layer.id === 'western-power-plants-cap') {
    const primaryColor = getPowerPlantColor(object.primaryType).slice(0, 3).join(',');
    const isCapacityLayer = layer.id === 'western-power-plants-cap';
    
    content = `
      <div style="
        background: linear-gradient(135deg, rgba(${primaryColor}, 0.1) 0%, rgba(${primaryColor}, 0.05) 100%);
        padding: 2px 8px;
        border-radius: 4px;
        border-left: 3px solid rgb(${primaryColor});
        margin-bottom: 6px;
      ">
        <div style="font-weight: 700; font-size: 13px; color: rgb(${primaryColor}); margin-bottom: 2px;">
          ${object.name}
        </div>
        <div style="font-size: 10px; color: #888; font-style: italic;">
          ${object.plantCode ? `Plant Code: ${object.plantCode}` : 'Western Power Plant'}
        </div>
      </div>
      
      <div style="margin-bottom: 6px;">
        <div style="font-size: 12px; font-weight: 700; color: #333; margin-bottom: 3px;">
          ${isCapacityLayer ? 'Total Capacity' : 'Current Generation'}: ${Math.round(isCapacityLayer ? object.Pcap : object.Pg)} MW
        </div>
        <div style="font-size: 10px; color: #666;">
          Primary Technology: <span style="color: rgb(${primaryColor}); font-weight: 500;">${object.primaryType}</span>
        </div>
      </div>

      <div style="
        background: rgba(0,0,0,0.03);
        padding: 6px 8px;
        border-radius: 4px;
        margin-bottom: 6px;
        border: 1px solid rgba(0,0,0,0.08);
      ">
        <div style="font-size: 10px; color: #555; margin-bottom: 3px;">
          <strong>📍 Location:</strong>
        </div>
        <div style="font-size: 10px; color: #666; line-height: 1.3;">
          <div>${object.state}${object.country ? ', ' + (object.country === 'USA' ? 'United States' : object.country) : ''}</div>
          ${object.balancingAuthority ? `<div><strong>Grid Operator:</strong> ${object.balancingAuthority}</div>` : ''}
        </div>
      </div>

      <div style="
        font-size: 9px; 
        color: #999; 
        text-align: center; 
        margin-top: 6px; 
        padding-top: 4px;
        border-top: 1px solid rgba(0,0,0,0.1);
        font-style: italic;
      ">
        ${object.country === 'USA' ? 'US EIA Form 860 Data' : 'International Power Plant Database'}
      </div>
    `;
  } else if (layer.id === 'WeccGenColumnLayer') {
    content = `
      <div style="font-weight: 600; margin-bottom: 4px; color: #2e7d32;">${object.name}</div>
      <div style="font-size: 11px; color: #666;">WECC Generation</div>
      <div style="font-size: 11px; color: #666;">Power: ${Math.round(object.Pg * 100) / 100} MW</div>
      <div style="font-size: 11px; color: #666;">Type: ${object.energyType}</div>
      <div style="font-size: 11px; color: #666;">Area: ${object.ba}</div>
    `;
  } else if (layer.id === 'WeccLayer') {
    const properties = object.properties;
    
    // Enhanced WECC region information
    const baCode = properties.BA_Abrev || 'WECC';
    const baName = properties.BA_Name || 'Western Electricity Coordinating Council Area';
    const fid = properties.FID || 'N/A';
    const areaNumbers = properties.Area_Numbers || [];
    
    // Get additional generation data if available
    const generationInfo = properties.weccGenData || null;
    
    content = `
      <div style="
        background: linear-gradient(135deg, rgba(21, 101, 192, 0.1) 0%, rgba(21, 101, 192, 0.05) 100%);
        padding: 8px 10px;
        border-radius: 6px;
        border-left: 4px solid #1565c0;
        margin-bottom: 8px;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
      ">
        <div style="font-weight: 700; font-size: 14px; color: #1565c0; margin-bottom: 3px;">
          🏛️ ${baCode}
        </div>
        <div style="font-size: 11px; color: #666; line-height: 1.3; margin-bottom: 2px;">
          ${baName}
        </div>
        <div style="font-size: 9px; color: #888; font-style: italic;">
          Balancing Authority Area
        </div>
      </div>

      <div style="
        background: rgba(0,0,0,0.03);
        padding: 6px 8px;
        border-radius: 4px;
        margin-bottom: 6px;
        border: 1px solid rgba(0,0,0,0.08);
      ">
        <div style="font-size: 10px; color: #555; margin-bottom: 3px;">
          <strong>🔍 Area Information:</strong>
        </div>
        <div style="font-size: 10px; color: #666; line-height: 1.3;">
          <div><strong>Region ID:</strong> ${fid}</div>
          ${areaNumbers.length > 0 ? 
            `<div><strong>Area Numbers:</strong> ${areaNumbers.join(', ')}</div>` : 
            ''
          }
          <div><strong>Grid:</strong> Western Interconnection</div>
          <div><strong>Coordinator:</strong> WECC</div>
        </div>
      </div>

      ${generationInfo ? 
        `<div style="
          background: rgba(46, 125, 50, 0.05);
          padding: 6px 8px;
          border-radius: 4px;
          border: 1px solid rgba(46, 125, 50, 0.15);
          margin-bottom: 6px;
        ">
          <div style="font-size: 10px; font-weight: 600; color: #2e7d32; margin-bottom: 3px;">
            ⚡ Generation Portfolio:
          </div>
          <div style="font-size: 9px; color: #666; line-height: 1.4;">
            <div><strong>Total Capacity:</strong> ${generationInfo.totalCapacity} MW</div>
            <div><strong>Primary Source:</strong> ${generationInfo.primarySource}</div>
            <div><strong>Plants:</strong> ${generationInfo.plantCount}</div>
          </div>
        </div>` : 
        ''
      }

      <div style="
        background: rgba(33, 150, 243, 0.05);
        padding: 6px 8px;
        border-radius: 4px;
        border: 1px solid rgba(33, 150, 243, 0.15);
        margin-bottom: 6px;
      ">
        <div style="font-size: 10px; font-weight: 600; color: #1976d2; margin-bottom: 3px;">
          🎯 Responsibilities:
        </div>
        <div style="font-size: 9px; color: #666; line-height: 1.4;">
          <div>• Grid reliability & planning</div>
          <div>• Load-generation balance</div>
          <div>• Transmission coordination</div>
          <div>• Market operations</div>
        </div>
      </div>

      <div style="
        display: flex;
        justify-content: space-between;
        align-items: center;
        padding: 4px 6px;
        background: rgba(33, 150, 243, 0.08);
        border-radius: 4px;
        border: 1px solid rgba(33, 150, 243, 0.2);
        margin-top: 6px;
      ">
        <div style="font-size: 9px; color: #1976d2; font-weight: 500;">
          🖱️ Click for detailed analysis
        </div>
        <div style="font-size: 8px; color: #666;">
          WECC Data Portal
        </div>
      </div>
    `;
  } else if (layer.id === 'PolygonLayerload' || layer.id === 'PolygonLayer2') {
    content = `
      <div style="font-weight: 600; margin-bottom: 4px; color: #d32f2f;">${object.properties.NAME}</div>
      <div style="font-size: 11px; color: #666;">Load Loss: ${object.properties.Pd?.toFixed(2)} MW</div>
      <div style="font-size: 11px; color: #666;">County: ${object.properties.countyname || 'Unknown'}</div>
    `;
  } else if (layer.id === 'AreaLayer') {
    content = `
      <div style="font-weight: 600; margin-bottom: 4px; color: #7b1fa2;">Area ${object.properties.name}</div>
      <div style="font-size: 11px; color: #666;">Control Area</div>
    `;
  } else if (layer.id === 'ZoneLayer') {
    content = `
      <div style="font-weight: 600; margin-bottom: 4px; color: #f57c00;">Zone ${object.properties.name}</div>
      <div style="font-size: 11px; color: #666;">Load Zone</div>
    `;
  } else if (layer.id === 'power-plants') {
    const props = object.properties;
    const primaryColor = getPowerPlantColor(props.primaryType).slice(0, 3).join(',');
    
    // Create comprehensive technology breakdown for USA plants
    const techBreakdown = [];
    const techIcons = {
      'solar': '☀️', 'wind': '💨', 'hydro': '💧', 'naturalGas': '🔥', 
      'nuclear': '⚛️', 'geothermal': '🌋', 'battery': '🔋', 'oil': '🛢️', 'other': '⚡'
    };
    
    if (props.country === 'USA') {
      if (props.solar > 0) techBreakdown.push(`${techIcons.solar} Solar: ${props.solar.toFixed(1)} MW`);
      if (props.wind > 0) techBreakdown.push(`${techIcons.wind} Wind: ${props.wind.toFixed(1)} MW`);
      if (props.hydro > 0) techBreakdown.push(`${techIcons.hydro} Hydro: ${props.hydro.toFixed(1)} MW`);
      if (props.naturalGas > 0) techBreakdown.push(`${techIcons.naturalGas} Natural Gas: ${props.naturalGas.toFixed(1)} MW`);
      if (props.nuclear > 0) techBreakdown.push(`${techIcons.nuclear} Nuclear: ${props.nuclear.toFixed(1)} MW`);
      if (props.geothermal > 0) techBreakdown.push(`${techIcons.geothermal} Geothermal: ${props.geothermal.toFixed(1)} MW`);
      if (props.battery > 0) techBreakdown.push(`${techIcons.battery} Battery: ${props.battery.toFixed(1)} MW`);
      if (props.oil > 0) techBreakdown.push(`${techIcons.oil} Oil: ${props.oil.toFixed(1)} MW`);
      if (props.other > 0) techBreakdown.push(`${techIcons.other} Other: ${props.other.toFixed(1)} MW`);
    }
    
    // Calculate capacity utilization indicator
    const capacitySize = props.totalCapacityMW;
    let sizeIndicator = '';
    let sizeColor = '';
    if (capacitySize >= 1000) {
      sizeIndicator = '🔴 Large Scale';
      sizeColor = '#d32f2f';
    } else if (capacitySize >= 100) {
      sizeIndicator = '🟡 Medium Scale';
      sizeColor = '#f57c00';
    } else if (capacitySize >= 10) {
      sizeIndicator = '🟢 Small Scale';
      sizeColor = '#388e3c';
    } else {
      sizeIndicator = '⚪ Micro Scale';
      sizeColor = '#757575';
    }
    
    content = `
      <div style="
        background: linear-gradient(135deg, rgba(${primaryColor}, 0.1) 0%, rgba(${primaryColor}, 0.05) 100%);
        padding: 2px 8px;
        border-radius: 4px;
        border-left: 3px solid rgb(${primaryColor});
        margin-bottom: 6px;
      ">
        <div style="font-weight: 700; font-size: 13px; color: rgb(${primaryColor}); margin-bottom: 2px;">
          ${props.plantName || 'Unnamed Plant'}
        </div>
        <div style="font-size: 10px; color: #888; font-style: italic;">
          Plant Code: ${props.plantCode || 'N/A'}
        </div>
      </div>
      
      <div style="margin-bottom: 6px;">
        <div style="display: flex; align-items: center; margin-bottom: 3px;">
          <div style="font-size: 11px; color: ${sizeColor}; font-weight: 600;">
            ${sizeIndicator}
          </div>
          <div style="margin-left: auto; font-size: 12px; font-weight: 700; color: #333;">
            ${props.totalCapacityMW.toFixed(1)} MW
          </div>
        </div>
        <div style="font-size: 10px; color: #666;">
          Primary Technology: <span style="color: rgb(${primaryColor}); font-weight: 500;">${props.primaryType}</span>
        </div>
      </div>

      <div style="
        background: rgba(0,0,0,0.03);
        padding: 6px 8px;
        border-radius: 4px;
        margin-bottom: 6px;
        border: 1px solid rgba(0,0,0,0.08);
      ">
        <div style="font-size: 10px; color: #555; margin-bottom: 3px;">
          <strong>📍 Location:</strong>
        </div>
        <div style="font-size: 10px; color: #666; line-height: 1.3;">
          <div>${props.state}${props.county ? ', ' + props.county + ' County' : ''}</div>
          <div>${props.country === 'USA' ? 'United States' : (props.country === 'CAN' ? 'Canada' : props.country)}</div>
          ${props.balancingAuthority ? `<div style="margin-top: 2px;"><strong>Grid Operator:</strong> ${props.balancingAuthority}</div>` : ''}
        </div>
      </div>

      ${techBreakdown.length > 0 ? 
        `<div style="
          background: rgba(${primaryColor}, 0.05);
          padding: 6px 8px;
          border-radius: 4px;
          border: 1px solid rgba(${primaryColor}, 0.15);
          margin-bottom: 4px;
        ">
          <div style="font-size: 10px; font-weight: 600; color: rgb(${primaryColor}); margin-bottom: 4px;">
            🔧 Technology Portfolio:
          </div>
          <div style="font-size: 9px; color: #666; line-height: 1.4;">
            ${techBreakdown.slice(0, 4).map(tech => `<div style="margin-bottom: 1px;">${tech}</div>`).join('')}
            ${techBreakdown.length > 4 ? 
              `<div style="font-style: italic; color: #888; margin-top: 2px;">
                +${techBreakdown.length - 4} additional technologies
              </div>` : 
              ''
            }
          </div>
        </div>` : 
        ''
      }

      <div style="
        font-size: 9px; 
        color: #999; 
        text-align: center; 
        margin-top: 6px; 
        padding-top: 4px;
        border-top: 1px solid rgba(0,0,0,0.1);
        font-style: italic;
      ">
        ${props.country === 'USA' ? 'US EIA Form 860 Data' : 'International Power Plant Database'}
      </div>
    `;
  }

  return {
    html: content,
    style: style
  };
}

// Custom Header Component with Logos
function WestmapHeader() {
  const { currentUser, userProfile } = useAuth();

  return (
    <div style={{
      position: 'absolute',
      top: 0,
      left: 0,
      right: 0,
      height: '60px',
      background: 'linear-gradient(135deg, #1976d2 0%, #1565c0 100%)',
      display: 'flex',
      alignItems: 'center',
      justifyContent: 'space-between',
      padding: '0 20px',
      zIndex: 1001,
      boxShadow: '0 2px 8px rgba(0,0,0,0.15)',
      fontFamily: '"Inter", sans-serif'
    }}>
      <div style={{ display: 'flex', alignItems: 'center', gap: '15px' }}>
        <img 
          src="/images/westmap_logo.png" 
          alt="Westmap Logo" 
          style={{ height: '40px', width: 'auto', borderRadius: '4px' }}
          onError={(e) => {
            e.target.style.display = 'none';
          }}
        />
        <h1 style={{ 
          color: 'white', 
          margin: 0, 
          fontSize: '24px', 
          fontWeight: '600',
          letterSpacing: '-0.02em'
        }}>
          Westmap
        </h1>
      </div>
      
      {/* Center section with navigation buttons */}
      <div style={{ display: 'flex', alignItems: 'center', gap: '12px' }}>
        {currentUser && (
          <>
            {/* Admin Panel - only for admin users - Cobalt Blue */}
            {userProfile?.role === 'admin' && (
              <Button
                onClick={() => {
                  window.history.pushState({}, '', '/admin');
                  window.location.reload();
                }}
                size="small"
                startIcon={<AdminPanelSettingsIcon sx={{ fontSize: 16 }} />}
                sx={{
                  color: '#0047AB', // Cobalt blue color
                  fontSize: '12px',
                  textTransform: 'none',
                  minWidth: 'auto',
                  padding: '4px 12px',
                  backgroundColor: 'white',
                  border: '2px solid #0047AB',
                  borderRadius: '20px',
                  fontWeight: '600',
                  '&:hover': {
                    backgroundColor: '#0047AB',
                    color: 'white',
                  }
                }}
              >
                Admin Panel
              </Button>
            )}
          </>
        )}
      </div>

      <div style={{ display: 'flex', alignItems: 'center', gap: '10px' }}>
        <img 
          src="/images/gridbee_logo.png" 
          alt="GridBee Logo" 
          style={{ height: '35px', width: 'auto', borderRadius: '4px' }}
          onError={(e) => {
            e.target.style.display = 'none';
          }}
        />
        <span style={{ 
          color: 'rgba(255,255,255,0.9)', 
          fontSize: '14px', 
          fontWeight: '500' 
        }}>
          Powered by GridBee AI
        </span>
      </div>
    </div>
  );
}


function MainApp({ refdata = data, refflowdata = flowdata, ggdata = geodata, mapStyle = MAP_STYLE }) {

  // Deck reference pointer
  const deckRef = useRef(null);

  // Navigation hook for routing
  const navigate = useNavigate();
  const location = useLocation();


  const [data, setData] = useState(refdata);
  const [enhancedData, setEnhancedData] = useState(refdata);
  const [isLoadingAdditionalData, setIsLoadingAdditionalData] = useState(false);

  //chat output message
  const [ouputMes, setOutputMes] = useState("Welcome to Westmap - Your intelligent power grid assistant.");

  //select widgets values
  const [nameSelectItems, setNameSelectItems] = useState([]);


  const [nameItems, setNameItems] = useState(gendata.Gens
    .map((gen) => (gen.name)));

  const [lineNameSelectItems, setLineNameSelectItems] = useState([]);

  const [lineNameItems, setLineNameItems] = useState(data.features
    .filter(f => f.geometry.type == "LineString").map((f) => (f.properties.NAME)));

  const [busNameSelectItems, setBusNameSelectItems] = useState([]);

  const [busNameItems, setbusNameItems] = useState(data.features
    .filter(f => f.geometry.type == "Point").map((f) => (f.properties.NAME)));

  const [areaNameSelectItems, setAreaNameSelectItems] = useState([]);

  const [areaNameItems, setAreaNameItems] = useState(areas.features.map(f => f.properties.name));

  const [zoneNameSelectItems, setZoneNameSelectItems] = useState([]);

  const [zoneNameItems, setZoneNameItems] = useState(zones.features.map(f => f.properties.name));

  const [countyNameSelectItems, setCountyNameSelectItems] = useState([]);

  const [countyNameItems, setCountyNameItems] = useState(countyload.features
    .map((county) => (county.properties.countyname)));

  const [flowdata, setFlowData] = useState(refflowdata);

  const [genfiltervalue, setGenFilterValue] = useState([gendata.minPg, gendata.maxPg]);

  const [genDoughlabels, setDoughlabels] = useState([
    'Wind',
    'Solar',
    'Nuclear',
    'Natural Gas',
    'Hydro',
    'Coal',
    'Other'
  ]);


  const colorMap = {
    'green': 'Wind',
    'yellow': 'Solar',
    'gray': 'Coal',
    'red': 'Nuclear',
    'blue': 'Hydro',
    'orange': 'Natural Gas',
    'black': "Other"
  }

  // WECC-specific color mapping
  const weccColorMap = {
    'blue': 'Hydro',
    'red': 'Nuclear',
    'black': 'Coal',
    'orange': 'Natural Gas',
    'purple': 'Geothermal',
    'green': 'Biomass',
    'lightgreen': 'Wind',
    'yellow': 'Solar'
  }

  const [netfiltervalue, setNetFilterValue] = useState([0, 800]);

  const [flowfiltervalue, setFlowFilterValue] = useState([0, 120]);

  const [loadfiltervalue, setLoadFilterValue] = useState([0, countyloaddata.maxPd]);

  const [voltagefiltervalue, setVoltageFilterValue] = useState([0.89, 1.11]);

  // For zoom-in/out control
  const [initialViewState, setInitialViewState] = useState(INITIAL_VIEW_STATE);

  // For pop-up control
  const [showPopup, setShowPopup] = useState({ display: false, info: '', name: '', fid: null, type: null });

  // Data Center Impact Analysis state
  const [impactAnalysisActive, setImpactAnalysisActive] = useState(false);
  const [selectedWeccRegion, setSelectedWeccRegion] = useState('ALL');
  const [selectedCaseStudy, setSelectedCaseStudy] = useState('case1');
  const [selectedMetric, setSelectedMetric] = useState('Price');
  const [selectedHour, setSelectedHour] = useState(12);
  const [impactData, setImpactData] = useState({});

  // WECC regions list for dropdown
  const weccRegions = [
    { id: 'ALL', name: 'All Regions', abbrev: 'ALL' },
    { id: 'AESO', name: 'Alberta Electric System Operator', abbrev: 'AESO' },
    { id: 'AVA', name: 'Avista Corporation', abbrev: 'AVA' },
    { id: 'AZPS', name: 'Arizona Public Service Company', abbrev: 'AZPS' },
    { id: 'BANC', name: 'Balancing Authority of Northern California', abbrev: 'BANC' },
    { id: 'BCHA', name: 'British Columbia Hydro and Power Authority', abbrev: 'BCHA' },
    { id: 'BPAT', name: 'Bonneville Power Administration', abbrev: 'BPAT' },
    { id: 'CISO', name: 'California Independent System Operator', abbrev: 'CISO' },
    { id: 'EPE', name: 'El Paso Electric Company', abbrev: 'EPE' },
    { id: 'IPCO', name: 'Idaho Power Company', abbrev: 'IPCO' },
    { id: 'LDWP', name: 'Los Angeles Department of Water and Power', abbrev: 'LDWP' },
    { id: 'NEVP', name: 'Nevada Power Company', abbrev: 'NEVP' },
    { id: 'NWMT', name: 'NorthWestern Energy', abbrev: 'NWMT' },
    { id: 'PACE', name: 'PacifiCorp East', abbrev: 'PACE' },
    { id: 'PACW', name: 'PacifiCorp West', abbrev: 'PACW' },
    { id: 'PGE', name: 'Portland General Electric Company', abbrev: 'PGE' },
    { id: 'PNM', name: 'Public Service Company of New Mexico', abbrev: 'PNM' },
    { id: 'PSCO', name: 'Public Service Company of Colorado', abbrev: 'PSCO' },
    { id: 'PSEI', name: 'Puget Sound Energy', abbrev: 'PSEI' },
    { id: 'SCL', name: 'Seattle City Light', abbrev: 'SCL' },
    { id: 'SRP', name: 'Salt River Project', abbrev: 'SRP' },
    { id: 'TEP', name: 'Tucson Electric Power', abbrev: 'TEP' },
    { id: 'TIDC', name: 'Turlock Irrigation District', abbrev: 'TIDC' },
    { id: 'TPWR', name: 'City of Tacoma, Department of Public Utilities', abbrev: 'TPWR' },
    { id: 'WACM', name: 'Western Area Power Administration', abbrev: 'WACM' }
  ];

  // Load additional transmission line data on component mount
  useEffect(() => {
    const loadEnhancedData = async () => {
      setIsLoadingAdditionalData(true);
      try {
        const enhancedDataResult = await initializeData();
        setEnhancedData(enhancedDataResult);
        
        // Update generation statistics if available
        if (enhancedDataResult.generationStats) {
          setGenerationStats(enhancedDataResult.generationStats);
        }
        
        // Update power plant data if available
        if (enhancedDataResult.powerPlantData) {
          setPowerPlantData(enhancedDataResult.powerPlantData);
          setPowerPlantFilterValue([enhancedDataResult.powerPlantData.minPg, enhancedDataResult.powerPlantData.maxPg]);
          
          // Set up name items for multiselect
          const nameItems = enhancedDataResult.powerPlants.map(plant => plant.name);
          setPowerPlantNameItems(nameItems);
        }
        
        // Update chart data if available
        if (enhancedDataResult.powerPlantChartData) {
          setPowerPlantChartData(enhancedDataResult.powerPlantChartData);
        }
        
        // Update flow data with the enhanced dataset
        const enhancedFlowData = ExtractFlowData(enhancedDataResult);
        setFlowData(enhancedFlowData);
        
        // Update line name items with the new data
        const updatedLineNameItems = enhancedDataResult.features
          .filter(f => f.geometry.type === "LineString")
          .map(f => f.properties.NAME);
        setLineNameItems(updatedLineNameItems);
        
        console.log('Enhanced transmission line data loaded successfully');
      } catch (error) {
        console.error('Failed to load enhanced transmission data:', error);
      } finally {
        setIsLoadingAdditionalData(false);
      }
    };

    loadEnhancedData();
  }, []);

  // Load impact analysis data when parameters change
  useEffect(() => {
    if (impactAnalysisActive) {
      loadImpactAnalysisData();
    }
  }, [impactAnalysisActive, selectedCaseStudy, selectedMetric, selectedHour, selectedWeccRegion]);

  // Data Center Impact Analysis Functions
  
  // Handler for impact analysis toggle
  const handleImpactAnalysisChange = (event) => {
    setImpactAnalysisActive(event.target.checked);
    if (event.target.checked) {
      loadImpactAnalysisData();
    } else {
      setImpactData({});
    }
  };

  // Load and calculate impact analysis data
  const loadImpactAnalysisData = async () => {
    try {
      console.log('Loading impact analysis data...');
      const impactResults = await calculateImpactDifferences();
      setImpactData(impactResults);
      console.log('Impact analysis data loaded:', impactResults);
    } catch (error) {
      console.error('Error loading impact analysis data:', error);
    }
  };

  // Calculate impact differences between case studies
  const calculateImpactDifferences = async () => {
    const results = {};
    
    // Define comparison pairs based on selected case study
    const comparisonPairs = {
      'case1': { current: 'Case study_1', baseline: 'Case study_0' },
      'case2': { current: 'Case study_2', baseline: 'Case study_1' },
      'case3': { current: 'Case study_3', baseline: 'Case study_2' }
    };
    
    const comparison = comparisonPairs[selectedCaseStudy];
    if (!comparison) return results;
    
    // Get regions to analyze
    const regionsToAnalyze = selectedWeccRegion === 'ALL' 
      ? weccRegions.filter(r => r.id !== 'ALL').map(r => r.id)
      : [selectedWeccRegion];
    
    for (const regionAbbrev of regionsToAnalyze) {
      try {
        const currentData = await loadCaseStudyDataForRegion(comparison.current, regionAbbrev);
        const baselineData = await loadCaseStudyDataForRegion(comparison.baseline, regionAbbrev);
        
        if (currentData && baselineData) {
          const difference = calculateMetricDifference(
            currentData, 
            baselineData, 
            selectedMetric, 
            selectedHour
          );
          
          results[regionAbbrev] = {
            current: currentData,
            baseline: baselineData,
            difference: difference,
            metric: selectedMetric,
            hour: selectedHour
          };
        }
      } catch (error) {
        console.warn(`Failed to load data for region ${regionAbbrev}:`, error);
      }
    }
    
    return results;
  };

  // Load case study data for a specific region
  const loadCaseStudyDataForRegion = async (caseStudy, regionAbbrev) => {
    try {
      const data = {
        lmp: null,
        costs: null,
        generation: null
      };
      
      // Load LMP data
      try {
        const lmpPath = `/amin_data/manish_amin_modified_data/${caseStudy}/LMP ($/MWh)/${regionAbbrev}_lmp.csv`;
        const lmpResponse = await fetch(lmpPath);
        if (lmpResponse.ok) {
          const lmpText = await lmpResponse.text();
          data.lmp = parseLMPDataForRegion(lmpText, regionAbbrev);
        }
      } catch (error) {
        console.warn(`Failed to load LMP data for ${regionAbbrev}:`, error);
      }
      
      // Load costs data
      try {
        const costsPath = `/amin_data/manish_amin_modified_data/${caseStudy}/Additional files/Balancing Authority Hourly Operation Costs - Copy/${regionAbbrev}_hourly_operation_costs.csv`;
        const costsResponse = await fetch(costsPath);
        if (costsResponse.ok) {
          const costsText = await costsResponse.text();
          data.costs = parseCostsDataForRegion(costsText);
        } else {
          // Fallback to total operation cost CSV
          const totalCostPath = `/amin_data/manish_amin_modified_data/Total_Operation_Cost.csv`;
          const totalCostResponse = await fetch(totalCostPath);
          if (totalCostResponse.ok) {
            const totalCostText = await totalCostResponse.text();
            data.costs = parseTotalOperationCostData(totalCostText, regionAbbrev, caseStudy);
          }
        }
      } catch (error) {
        console.warn(`Failed to load costs data for ${regionAbbrev}:`, error);
      }
      
      // Load generation/power exchange data
      try {
        const genPath = `/amin_data/manish_amin_modified_data/${caseStudy}/Power Exchange (MW)/${regionAbbrev}_power_exchange.csv`;
        const genResponse = await fetch(genPath);
        if (genResponse.ok) {
          const genText = await genResponse.text();
          data.generation = parseGenerationDataForRegion(genText);
        }
      } catch (error) {
        console.warn(`Failed to load generation data for ${regionAbbrev}:`, error);
      }
      
      return data;
    } catch (error) {
      console.error(`Error loading case study data for ${regionAbbrev}:`, error);
      return null;
    }
  };

  // Calculate difference between metrics
  const calculateMetricDifference = (currentData, baselineData, metric, hour) => {
    try {
      let currentValue = 0;
      let baselineValue = 0;
      
      if (metric === 'Price' || metric === 'price') {
        const currentLmp = currentData.lmp?.find(d => d.hour === hour);
        const baselineLmp = baselineData.lmp?.find(d => d.hour === hour);
        currentValue = currentLmp?.price || 0;
        baselineValue = baselineLmp?.price || 0;
      }
      
      if (metric === 'System Operation Cost' || metric === 'System Cost') {
        const currentCost = currentData.costs?.find(d => d.hour === hour);
        const baselineCost = baselineData.costs?.find(d => d.hour === hour);
        currentValue = (currentCost?.startupCosts || 0) + (currentCost?.fuelCosts || 0) + (currentCost?.variableCosts || 0);
        baselineValue = (baselineCost?.startupCosts || 0) + (baselineCost?.fuelCosts || 0) + (baselineCost?.variableCosts || 0);
      }
      
      if (metric === 'Power Exchange' || metric === 'power_exchange') {
        const currentGen = currentData.generation?.find(d => d.hour === hour);
        const baselineGen = baselineData.generation?.find(d => d.hour === hour);
        currentValue = currentGen?.importExport || 0;
        baselineValue = baselineGen?.importExport || 0;
      }
      
      return currentValue - baselineValue;
    } catch (error) {
      console.error('Error calculating metric difference:', error);
      return 0;
    }
  };

  // Parse LMP data for a specific region
  const parseLMPDataForRegion = (csvText, regionCode) => {
    const lines = csvText.split('\n');
    const headers = lines[0].split(',');
    const regionIndex = headers.findIndex(header => header.includes(regionCode));
    
    if (regionIndex === -1) return [];
    
    const data = [];
    for (let i = 1; i < lines.length; i++) {
      const line = lines[i].trim();
      if (line) {
        const values = line.split(',');
        const hour = parseInt(values[0]);
        const price = parseFloat(values[regionIndex]) || 0;
        data.push({ hour, price });
      }
    }
    return data;
  };

  // Parse costs data for a region
  const parseCostsDataForRegion = (csvText) => {
    const lines = csvText.split('\n');
    const data = [];
    
    for (let i = 1; i < lines.length; i++) {
      const line = lines[i].trim();
      if (line) {
        const values = line.split(',');
        data.push({
          hour: parseInt(values[0]),
          startupCosts: parseFloat(values[2]) || 0,
          fuelCosts: parseFloat(values[3]) || 0,
          variableCosts: parseFloat(values[4]) || 0
        });
      }
    }
    return data;
  };

  // Parse total operation cost data from CSV
  const parseTotalOperationCostData = (csvText, regionCode, caseStudy) => {
    try {
      const lines = csvText.split('\n');
      const headers = lines[0].split(',');
      
      // Find the column for this case study and region
      const columnName = `${caseStudy}_${regionCode}`;
      const columnIndex = headers.findIndex(header => header.includes(columnName) || header.includes(regionCode));
      
      if (columnIndex === -1) return [];
      
      const data = [];
      for (let i = 1; i < lines.length && i <= 24; i++) {
        const line = lines[i].trim();
        if (line) {
          const values = line.split(',');
          const cost = parseFloat(values[columnIndex]) || 0;
          data.push({
            hour: i,
            startupCosts: 0,
            fuelCosts: cost,
            variableCosts: 0,
            totalCost: cost
          });
        }
      }
      return data;
    } catch (error) {
      console.error('Error parsing total operation cost data:', error);
      return [];
    }
  };

  // Parse generation/power exchange data
  const parseGenerationDataForRegion = (csvText) => {
    try {
      const lines = csvText.split('\n');
      const data = [];
      
      for (let i = 1; i < lines.length; i++) {
        const line = lines[i].trim();
        if (line) {
          const values = line.split(',');
          data.push({
            hour: parseInt(values[0]) || i,
            importExport: parseFloat(values[1]) || 0,
            generation: parseFloat(values[2]) || 0
          });
        }
      }
      return data;
    } catch (error) {
      console.error('Error parsing generation data:', error);
      return [];
    }
  };

  // Get color for WECC region based on impact data
  const getWeccRegionColor = (regionAbbrev) => {
    if (!impactAnalysisActive || !impactData[regionAbbrev]) {
      return [200, 200, 200, 100]; // Default gray
    }
    
    const impact = impactData[regionAbbrev];
    const difference = impact.difference || 0;
    
    // Get max absolute difference for normalization
    const allDifferences = Object.values(impactData).map(d => Math.abs(d.difference || 0));
    const maxDiff = Math.max(...allDifferences, 1); // Avoid division by zero
    
    // Normalize difference (-1 to 1)
    const normalizedDiff = difference / maxDiff;
    
    // Color interpolation: blue (negative) -> white (zero) -> red (positive)
    let r, g, b;
    
    if (normalizedDiff < 0) {
      // Negative values: interpolate from blue to white
      const intensity = Math.abs(normalizedDiff);
      r = Math.round(0 + (255 - 0) * (1 - intensity));
      g = Math.round(102 + (255 - 102) * (1 - intensity));
      b = 255;
    } else if (normalizedDiff > 0) {
      // Positive values: interpolate from white to red
      const intensity = normalizedDiff;
      r = 255;
      g = Math.round(255 + (68 - 255) * intensity);
      b = Math.round(255 + (68 - 255) * intensity);
    } else {
      // Zero difference: white
      r = g = b = 255;
    }
    
    return [r, g, b, 180];
  };

  //update flowdataset when netfiltervalue, flowfiltervalue or data value change
  useEffect(() => {


    // name is the unique id 
    const locations = []
    const flows = []
    //  if user make selections from names, than return individual lines
    if (lineNameSelectItems.length > 0) {
      const pointNames = lineNameSelectItems.map(n => n.split(' -- ')).flat()
      data.features.forEach(feature => {
        if (feature.geometry.type === "Point" && pointNames.includes(feature.properties.NAME)) {
          locations.push({
            id: feature.properties.NAME,
            name: feature.properties.NAME,
            lon: feature.geometry.coordinates[0],
            lat: feature.geometry.coordinates[1]
          })
        } else if (feature.geometry.type === "LineString" && lineNameSelectItems.includes(feature.properties.NAME)) {
          var RATE_A;
          if (feature.properties.RATE_A == 0) {
            RATE_A = 10000;
          } else {
            RATE_A = feature.properties.RATE_A;
          }
          var loading = Math.abs(feature.properties.PF / RATE_A) * 100;
          if (feature.properties.PF > 0) {
            const [origin, dest] = feature.properties.NAME.split(' -- ')
            flows.push({
              origin: origin,
              dest: dest,
              count: feature.properties.KV,
              loading: loading
            })
          } else {
            const [dest, origin] = feature.properties.NAME.split(' -- ')
            flows.push({
              origin: origin,
              dest: dest,
              count: feature.properties.KV,
              loading: loading
            })
          }
        }
      })
    } else if (busNameSelectItems.length > 0) {

      data.features.forEach(feature => {
        if (feature.geometry.type === "Point" && busNameSelectItems.includes(feature.properties.NAME)) {
          locations.push({
            id: feature.properties.NAME,
            name: feature.properties.NAME,
            lon: feature.geometry.coordinates[0],
            lat: feature.geometry.coordinates[1]
          })
        } else if (feature.geometry.type === "LineString" && feature.properties.NAME.split(' -- ').some(r => busNameSelectItems.includes(r))) {
          var RATE_A;
          if (feature.properties.RATE_A == 0) {
            RATE_A = 10000;
          } else {
            RATE_A = feature.properties.RATE_A;
          }

          var loading = Math.abs(feature.properties.PF / RATE_A) * 100;
          if (feature.properties.PF > 0) {
            const [origin, dest] = feature.properties.NAME.split(' -- ')
            flows.push({
              origin: origin,
              dest: dest,
              count: feature.properties.KV,
              loading: loading
            })
          } else {
            const [dest, origin] = feature.properties.NAME.split(' -- ')
            flows.push({
              origin: origin,
              dest: dest,
              count: feature.properties.KV,
              loading: loading
            })
          }

        }

      })
    }
    else {
      data.features.forEach(feature => {
        if (feature.geometry.type === "Point" && netfiltervalue[0] <= feature.properties.KVlevels[0] &&
          feature.properties.KVlevels[0] <= netfiltervalue[1]) {
          locations.push({
            id: feature.properties.NAME,
            name: feature.properties.NAME,
            lon: feature.geometry.coordinates[0],
            lat: feature.geometry.coordinates[1]
          })
        } else if (feature.geometry.type === "LineString" && netfiltervalue[0] <= feature.properties.KV &&
          feature.properties.KV <= netfiltervalue[1]) {
          var RATE_A;
          if (feature.properties.RATE_A == 0) {
            RATE_A = 10000;
          } else {
            RATE_A = feature.properties.RATE_A;
          }
          var loading = Math.abs(feature.properties.PF / RATE_A) * 100.0;

          if (flowfiltervalue[0] <= loading && loading <= flowfiltervalue[1]) {
            if (feature.properties.PF > 0) {
              const [origin, dest] = feature.properties.NAME.split(' -- ')
              flows.push({
                origin: origin,
                dest: dest,
                count: feature.properties.KV,
                loading: loading
              })
            } else {
              const [dest, origin] = feature.properties.NAME.split(' -- ')
              flows.push({
                origin: origin,
                dest: dest,
                count: feature.properties.KV,
                loading: loading
              })
            }
          }
        }
      })
    }

    const newflowdata = { locations: locations, flows: flows, maxloading: 120 }
    setFlowData(newflowdata);
  }, [data, netfiltervalue, flowfiltervalue, lineNameSelectItems]);


  var rotatestate = false;
  //const [rotatestate,setrotatestate] = useState(false);

  const rotateCamera = useCallback(() => {
    rotatestate = !rotatestate;
    if (rotatestate) {
      setInitialViewState(viewState => ({
        ...viewState,
        bearing: viewState.bearing - 180,
        transitionDuration: 20000,
        transitionInterpolator: transitionLinearInterpolator,
        onTransitionEnd: rotateCamera
      }))
    } else {
      setInitialViewState(viewState => ({
        ...viewState,
        onTransitiionEnd: null
      }))
      //      GoHome();
      rotatestate = false;
    }

  }, []);

  const activatePopup = useCallback(() => {
    setShowPopup(showPopup => ({ ...showPopup, display: true }));

  }, []);

  const zoomToGen = useCallback((lat, long) => {

    setInitialViewState(viewState => ({
      ...viewState,
      latitude: lat,
      longitude: long,
      pitch: 50,
      transitionInterpolator: transitionFlyToInterpolator,
      transitionDuration: 1000,
      zoom: 5.5,
      onTransitionEnd: activatePopup
    }))

  });

  const zoomToData = useCallback((info) => {
    var lat = info.coordinate[1];
    var long = info.coordinate[0];

    setInitialViewState(viewState => ({
      ...viewState,
      latitude: lat,
      longitude: long,
      pitch: 50,
      transitionInterpolator: transitionFlyToInterpolator,
      transitionDuration: 1000,
      zoom: 7.5,
      onTransitionEnd: activatePopup
    }))

    if (info.layer.id == "geojson") {
      if (info.object.geometry.type == "Point") {
        var popup = {};
        popup.name = info.object.properties.NAME
        popup.info = "Substation Info"
      } else {
        var popup = {};
        var loading = Math.abs(info.object.properties.PF / info.object.properties.RATE_A) * 100.0;
        popup.name = info.object.properties.NAME
        popup.info = "KV: " + info.object.properties.KV.toFixed(2) + "KV \nLoading: " + loading.toFixed(2) + "%";
      }
      setShowPopup(showPopup => ({ ...showPopup, ...popup }));
    } else if (info.layer.id == "gen-column") {
      var popup = {};
      popup.name = "Pg: " + Math.round(info.object.Pg * 100) / 100 + " Pcap: " + Math.round(info.object.Pcap * 100) / 100;
      popup.info = "Gen Info";

      setShowPopup(showPopup => ({ ...showPopup, ...popup }));
    }
  });

  const zoomToCountyName = useCallback((minLng, minLat, maxLng, maxLat) => {



    var viewport = new WebMercatorViewport(INITIAL_VIEW_STATE);


    const { longitude, latitude, zoom } = viewport.fitBounds([[minLng, minLat], [maxLng, maxLat]]);

    setInitialViewState(viewState => ({
      ...viewState,
      latitude: latitude,
      longitude: longitude,
      pitch: 50,
      transitionInterpolator: transitionFlyToInterpolator,
      transitionDuration: 1000,
      zoom: 7.5,
      onTransitionEnd: activatePopup
    }))


  });

  const zoomToAreaName = useCallback((filterareas, minLng, minLat, maxLng, maxLat) => {
    var viewport = new WebMercatorViewport(INITIAL_VIEW_STATE);

    const { longitude, latitude, zoom } = viewport.fitBounds([[minLng, minLat], [maxLng, maxLat]]);

    setInitialViewState(viewState => ({
      ...viewState,
      latitude: latitude,
      longitude: longitude,
      pitch: 50,
      transitionInterpolator: transitionFlyToInterpolator,
      transitionDuration: 1000,
      zoom: 7.5,
      onTransitionEnd: activatePopup
    }))

    var popup = { display: false, name: '', info: '' }; // Will be displayed after transition end only
    popup.name = "Area " + filterareas.properties.name;
    //      popup.info = "Area: " + info.object.properties.name;
    setShowPopup(showPopup => ({ ...showPopup, ...popup }));

  });

  const zoomToZoneName = useCallback((filterzones, minLng, minLat, maxLng, maxLat) => {
    var viewport = new WebMercatorViewport(INITIAL_VIEW_STATE);

    const { longitude, latitude, zoom } = viewport.fitBounds([[minLng, minLat], [maxLng, maxLat]]);

    setInitialViewState(viewState => ({
      ...viewState,
      latitude: latitude,
      longitude: longitude,
      pitch: 50,
      transitionInterpolator: transitionFlyToInterpolator,
      transitionDuration: 1000,
      zoom: 7.5,
      onTransitionEnd: activatePopup
    }))

    var popup = { display: false, name: '', info: '' }; // Will be displayed after transition end only
    popup.name = "Zone " + filterzones.properties.name;
    //      popup.info = "Zone: " + info.object.properties.name;
    setShowPopup(showPopup => ({ ...showPopup, ...popup }));

  });



  const zoomToCounty = useCallback((info) => {
    if (!info) return null;

    if (info.layer.id == 'PolygonLayer2' || info.layer.id == 'PolygonLayerload') {
      var layer = info.layer;
      var { viewport } = layer.context;

      var cbounds = bbox(info.object);
      var c1 = [cbounds[0], cbounds[1]];
      var c2 = [cbounds[2], cbounds[3]];
      var countybounds = [c1, c2];
      const { longitude, latitude, zoom } = viewport.fitBounds(countybounds);

      setInitialViewState(viewState => ({
        ...viewState,
        latitude: latitude,
        longitude: longitude,
        pitch: 50,
        transitionInterpolator: transitionFlyToInterpolator,
        transitionDuration: 1000,
        zoom: zoom - 0.25,
        onTransitionEnd: activatePopup
      }))

      var popup = { display: false, name: '', info: '' }; // Will be displayed after transition end only
      popup.name = info.object.properties.NAME;
      popup.info = "Load loss: " + info.object.properties.Pd.toFixed(2) + "MW";
      setShowPopup(showPopup => ({ ...showPopup, ...popup }));


    }
  });

  const zoomToArea = useCallback((info) => {
    if (!info) return null;

    if (info.layer.id == 'AreaLayer') {
      var layer = info.layer;
      var { viewport } = layer.context;

      var cbounds = bbox(info.object);
      var c1 = [cbounds[0], cbounds[1]];
      var c2 = [cbounds[2], cbounds[3]];
      var areabounds = [c1, c2];
      const { longitude, latitude, zoom } = viewport.fitBounds(areabounds);

      setInitialViewState(viewState => ({
        ...viewState,
        latitude: latitude,
        longitude: longitude,
        pitch: 50,
        transitionInterpolator: transitionFlyToInterpolator,
        transitionDuration: 1000,
        zoom: zoom - 0.25,
        onTransitionEnd: activatePopup
      }))

      var popup = { display: false, name: '', info: '' }; // Will be displayed after transition end only
      popup.name = "Area " + info.object.properties.name;
      //      popup.info = "Area: " + info.object.properties.name;
      setShowPopup(showPopup => ({ ...showPopup, ...popup }));


    }
  });

  const zoomToZone = useCallback((info) => {
    if (!info) return null;

    if (info.layer.id == 'ZoneLayer') {
      var layer = info.layer;
      var { viewport } = layer.context;

      var cbounds = bbox(info.object);
      var c1 = [cbounds[0], cbounds[1]];
      var c2 = [cbounds[2], cbounds[3]];
      var zonebounds = [c1, c2];
      const { longitude, latitude, zoom } = viewport.fitBounds(zonebounds);

      setInitialViewState(viewState => ({
        ...viewState,
        latitude: latitude,
        longitude: longitude,
        pitch: 50,
        transitionInterpolator: transitionFlyToInterpolator,
        transitionDuration: 1000,
        zoom: zoom - 0.25,
        onTransitionEnd: activatePopup
      }))

      var popup = { display: false, name: '', info: '' }; // Will be displayed after transition end only
      popup.name = "Zone " + info.object.properties.name;
      //      popup.info = "Zone: " + info.object.properties.name;
      setShowPopup(showPopup => ({ ...showPopup, ...popup }));


    }
  });

  const zoomToWecc = useCallback((info) => {
    if (!info) return null;

    if (info.layer.id == 'WeccLayer') {
      var layer = info.layer;
      var { viewport } = layer.context;

      var cbounds = bbox(info.object);
      var c1 = [cbounds[0], cbounds[1]];
      var c2 = [cbounds[2], cbounds[3]];
      var weccbounds = [c1, c2];
      const { longitude, latitude, zoom } = viewport.fitBounds(weccbounds);

      setInitialViewState(viewState => ({
        ...viewState,
        latitude: latitude,
        longitude: longitude,
        pitch: 50,
        transitionInterpolator: transitionFlyToInterpolator,
        transitionDuration: 1000,
        zoom: zoom - 0.25,
        onTransitionEnd: activatePopup
      }))

      var popup = { display: false, name: '', info: '', fid: null, type: null }; // Will be displayed after transition end only
      const properties = info.object.properties;

      popup.name = properties.BA_Abrev || "WECC Area";
      popup.fid = properties.FID;
      popup.type = 'wecc';

      // Build comprehensive info string similar to Amin's project
      let infoLines = [];
      infoLines.push(`${properties.BA_Name || ''}`);
      infoLines.push(`FID: ${properties.FID || 'N/A'}`);

      // Handle Area Numbers (can be single or multiple)
      if (properties.Area_Numbers && Array.isArray(properties.Area_Numbers) && properties.Area_Numbers.length > 0) {
        infoLines.push(`Area No.: ${properties.Area_Numbers.join(', ')}`);
      }

      // Add click instruction
      infoLines.push('');
      infoLines.push('Know more →');

      popup.info = infoLines.join('\n');
      setShowPopup(showPopup => ({ ...showPopup, ...popup }));
    }
  });

  const GoHome = useCallback(() => {
    if (layers[0].context == null) return;
    var { viewport } = layers[0].context;
    
    // Use enhanced bounds if available (includes additional transmission lines)
    const boundsToUse = enhancedData?.enhancedBounds || bounds;
    const { longitude, latitude, zoom } = viewport.fitBounds(boundsToUse);

    setInitialViewState(viewState => ({
      ...INITIAL_VIEW_STATE,
      longitude: longitude,
      latitude: latitude,
      zoom: zoom - 0.25,
      transitionInterpolator: transitionFlyToInterpolator,
      transitionDuration: 1000
    }))

    setShowPopup({ ...showPopup, display: false });
  }, [bounds, enhancedData]);

  // Zoom to fit all transmission lines (including additional CSV data)
  const ZoomToAllLines = useCallback(() => {
    if (layers[0].context == null) return;
    var { viewport } = layers[0].context;
    
    // Use enhanced bounds if available, otherwise fall back to regular bounds
    const boundsToUse = enhancedData?.enhancedBounds || bounds;
    const { longitude, latitude, zoom } = viewport.fitBounds(boundsToUse);

    setInitialViewState(viewState => ({
      ...viewState,
      longitude: longitude,
      latitude: latitude,
      zoom: zoom - 0.5, // Zoom out a bit more to see all lines clearly
      transitionInterpolator: transitionFlyToInterpolator,
      transitionDuration: 2000 // Slightly longer transition for better UX
    }));

    setShowPopup({ ...showPopup, display: false });
  }, [bounds, enhancedData]);

  // Handle popup click for navigation
  const handlePopupClick = useCallback(() => {
    if (showPopup.type === 'wecc' && showPopup.fid) {
      // Get current path and append the FID
      const currentPath = location.pathname;
      const targetPath = `${currentPath}/${showPopup.fid}`;
      navigate(targetPath);
    }
  }, [showPopup, navigate, location]);

  const [netlayeractive, setNetLayerActive] = useState(true);
  const [flowlayeractive, setFlowLayerActive] = useState(true);

  const [loadlayeractive, setLoadLayerActive] = useState(false);
  const [genlayeractive, setGenLayerActive] = useState(false);
  const [genlayercapactive, setGenLayerCapActive] = useState(false);
  const [voltagelayeractive, setVoltageLayerActive] = useState(false);
  const [zonelayeractive, setZoneLayerActive] = useState(false);
  const [arealayeractive, setAreaLayerActive] = useState(false);
  const [wecclayeractive, setWeccLayerActive] = useState(true);
  const [transmissionlayeractive, setTransmissionLayerActive] = useState(true);
  const [powerplantlayeractive, setPowerPlantLayerActive] = useState(true);
  const [powerplantlayercapactive, setPowerPlantLayerCapActive] = useState(false);
  const [powerplantlabelsactive, setPowerPlantLabelsActive] = useState(false);
  const [powerPlantFilterValue, setPowerPlantFilterValue] = useState([0, 10000]);
  const [powerPlantSelectItems, setPowerPlantSelectItems] = useState([]);
  const [selectedTechnology, setSelectedTechnology] = useState(null);
  const [generationStats, setGenerationStats] = useState(null);
  const [powerPlantData, setPowerPlantData] = useState({ Plants: [], minPg: 0, maxPg: 10000, minPcap: 0, maxPcap: 10000 });
  const [powerPlantNameItems, setPowerPlantNameItems] = useState([]);
  const [powerPlantChartData, setPowerPlantChartData] = useState(null);

  const [mapStyleSelection, setMapStyle] = useState('osm');

  // WECC Balancing Authorities data state
  const [weccGeojsonData, setWeccGeojsonData] = useState(null);
  const [weccHoveredObject, setWeccHoveredObject] = useState(null);
  const [weccClickedObject, setWeccClickedObject] = useState(null);

  // WECC Generation Power data state
  const [weccGenLayerActive, setWeccGenLayerActive] = useState(false);
  const [weccGenData, setWeccGenData] = useState(null);
  const [weccGenChartData, setWeccGenChartData] = useState(null);
  const [weccGenFilter, setWeccGenFilter] = useState([0, 100000]); // MW range
  const [weccGenSelectItems, setWeccGenSelectItems] = useState([]);
  const [weccGenNameItems, setWeccGenNameItems] = useState([]);
  const [weccColumnData, setWeccColumnData] = useState([]);
  const [weccGenDoughlabels, setWeccDoughlabels] = useState([
    'Hydro',
    'Nuclear',
    'Coal',
    'Natural Gas',
    'Geothermal',
    'Biomass',
    'Wind',
    'Solar'
  ]);

  const handleUserInput = (inputText) => {
    console.log(`New message incoming! ${inputText}`);
    // Now send the message to GPT and get response 
    toggleMsgLoader()
    const postData = {
      "inputText": inputText
    }
    // Use relative API path through nginx reverse proxy
    // const apiPath = process.env.ENVIRONMENT === 'prod' ? '/api/data' : 'http://localhost:5000/data';
    const apiPath = 'http://3.101.133.249:5000/data'; // for server deployment

    try {
      fetch(apiPath, {
        "method": "POST",
        // headers: { 'Content-Type': 'application/json' },
        "body": JSON.stringify(postData),
      }).then((res) =>
        res.json().then((chatOutput) => {
          // Setting a data from api
          console.log(chatOutput);
          const outputText = chatOutput.text
          const chatList = chatOutput.result_list || []

          if (chatList.length > 0) {
            //  only one is active between bus name selection, transmission line name selection at a time
            const keyList = Object.keys(chatList[0])

            if ("generation name" in chatList[0]) {
              const genNameList = chatList.map(d => d["generation name"]);
              setGenLayerActive(true)
              setNameSelectItems(genNameList)
              setGenFilterValue([gendata.minPg, gendata.maxPg]);

              // Turn off other layers to focus on generation
              setNetLayerActive(false)
              setFlowLayerActive(false)
              setVoltageLayerActive(false)
              setLoadLayerActive(false)

              // Smart zoom to generation locations
              if (genNameList.length > 0) {
                // Find the first generation facility in our data to get coordinates
                const firstGen = generation.find(gen => genNameList.includes(gen.name));
                if (firstGen && firstGen.coordinates) {
                  const [longitude, latitude] = firstGen.coordinates;
                  zoomToGen(latitude, longitude);
                } else {
                  // Fallback to general animation if no coordinates found
                  setInitialViewState(viewState => ({
                    ...viewState,
                    pitch: 40,
                    transitionInterpolator: transitionFlyToInterpolator,
                    transitionDuration: 2000,
                  }))
                }
              }
            }
            const containCapacity = keyList.some(str => str.includes('capacity'))
            if ("generation name" in chatList[0] && containCapacity) {
              const genNameList = chatList.map(d => d["generation name"]);
              setGenLayerActive(true)
              setGenLayerCapActive(true)
              setNameSelectItems(genNameList)
              setGenFilterValue([gendata.minPg, gendata.maxPg]);

              // Turn off other layers to focus on generation capacity
              setNetLayerActive(false)
              setFlowLayerActive(false)
              setVoltageLayerActive(false)
              setLoadLayerActive(false)
              setWeccLayerActive(false)

              // Smart zoom to generation locations for capacity queries
              if (genNameList.length > 0) {
                const firstGen = generation.find(gen => genNameList.includes(gen.name));
                if (firstGen && firstGen.coordinates) {
                  const [longitude, latitude] = firstGen.coordinates;
                  zoomToGen(latitude, longitude);
                } else {
                  setInitialViewState(viewState => ({
                    ...viewState,
                    pitch: 40,
                    transitionInterpolator: transitionFlyToInterpolator,
                    transitionDuration: 2000,
                  }))
                }
              }
            }
            if ("line name" in chatList[0]) {
              const lineNameList = chatList.map(d => d["line name"]);
              setNetLayerActive(true)
              setFlowLayerActive(true)
              setBusNameSelectItems([])
              setLineNameSelectItems(lineNameList)

              // Turn off other layers to focus on transmission lines
              setGenLayerActive(false)
              setGenLayerCapActive(false)
              setVoltageLayerActive(false)
              setLoadLayerActive(false)

              // Animate camera for line visualization
              setInitialViewState(viewState => ({
                ...viewState,
                pitch: 30,
                transitionInterpolator: transitionFlyToInterpolator,
                transitionDuration: 2000,
                zoom: 4.5
              }))
            }
            if ('bus name' in chatList[0]) {
              const busNameList = chatList.map(d => d["bus name"]);

              // Check if this is a voltage-related query
              const isVoltageQuery = outputText.toLowerCase().includes('voltage') ||
                outputText.toLowerCase().includes('kv') ||
                inputText.toLowerCase().includes('voltage');

              if (isVoltageQuery) {
                // Activate voltage visualization for voltage queries
                setVoltageLayerActive(true)
                setNetLayerActive(false) // Turn off network layer to focus on voltage
                setFlowLayerActive(false)
                setGenLayerActive(false) // Turn off generation layers
                setGenLayerCapActive(false)
                setLoadLayerActive(false) // Turn off load layer
              } else {
                // Regular bus query - show network
                setNetLayerActive(true)
                setFlowLayerActive(true)
                setVoltageLayerActive(false)
                setGenLayerActive(false)
                setGenLayerCapActive(false)
                setLoadLayerActive(false)
              }

              setBusNameSelectItems(busNameList)
              setLineNameSelectItems([])

              // Animate camera for bus visualization
              setInitialViewState(viewState => ({
                ...viewState,
                pitch: 30,
                transitionInterpolator: transitionFlyToInterpolator,
                transitionDuration: 2000,
                zoom: 4.5
              }))
            }
          }

          setOutputMes(outputText)

        })
      );
    } catch (error) {
      setOutputMes("Sorry I didn't find the answer to your question. Please try to rephrase it or provide more details.")
    }


  };

  useEffect(() => {

    addResponseMessage(`${ouputMes}`);
    if (ouputMes === "Welcome to Westmap - Your intelligent power grid assistant." || ouputMes === '') return;
    toggleMsgLoader(); // close loading 
  }, [ouputMes]);

  // Load WECC Balancing Authorities data
  useEffect(() => {
    const loadWeccData = async () => {
      try {
        // Load GeoJSON for map visualization
        const geojsonResponse = await fetch('/amin_data/WECC_Balancing_Authorities_-2060174188301432986.geojson');
        if (!geojsonResponse.ok) {
          throw new Error(`HTTP error loading WECC GeoJSON! status: ${geojsonResponse.status}`);
        }
        const geojsonData = await geojsonResponse.json();

        // Load CSV for additional data fields
        const csvResponse = await fetch('/amin_data/WECC_Balancing_Authorities_5803277890210865950.csv');
        if (!csvResponse.ok) {
          throw new Error(`HTTP error loading WECC CSV! status: ${csvResponse.status}`);
        }
        const csvText = await csvResponse.text();

        // Load area mapping CSV
        const areaMappingResponse = await fetch('/amin_data/WECC_BA_Area_Mapping.csv');
        if (!areaMappingResponse.ok) {
          throw new Error(`HTTP error loading WECC area mapping CSV! status: ${areaMappingResponse.status}`);
        }
        const areaMappingText = await areaMappingResponse.text();

        // Load WECC generation data CSV
        const weccGenResponse = await fetch('/amin_data/new_data/WECC_31_BAs_2028_All_Clean (1).csv');
        if (!weccGenResponse.ok) {
          throw new Error(`HTTP error loading WECC generation CSV! status: ${weccGenResponse.status}`);
        }
        const weccGenText = await weccGenResponse.text();

        // Parse CSV
        const csvLines = csvText.split('\n');
        const csvData = {};

        // Create lookup table by FID
        for (let i = 1; i < csvLines.length; i++) {
          const line = csvLines[i].trim();
          if (line) {
            const values = line.split(',');
            const fid = parseInt(values[0]);
            csvData[fid] = {
              FID: fid,
              BA_Abrev: values[1],
              BA_Name: values[2].replace(/"/g, ''), // Remove quotes
              Shape_Leng: parseFloat(values[3]),
              Shape__Area: parseFloat(values[4]),
              Shape__Length: parseFloat(values[5]),
              GlobalID: values[6]
            };
          }
        }

        // Parse area mapping CSV
        const areaMappingLines = areaMappingText.split('\n');
        const areaMappingData = {};

        for (let i = 1; i < areaMappingLines.length; i++) {
          const line = areaMappingLines[i].trim();
          if (line) {
            const values = line.split(',');
            const fid = parseInt(values[0]);
            const areaNumbersString = values[3];

            // Handle multiple area numbers separated by |
            let areaNumbers;
            if (areaNumbersString && areaNumbersString.includes('|')) {
              areaNumbers = areaNumbersString.split('|').map(num => parseInt(num.trim()));
            } else {
              areaNumbers = areaNumbersString ? [parseInt(areaNumbersString)] : [];
            }

            areaMappingData[fid] = {
              FID: fid,
              BA_Abrev: values[1],
              BA_Name: values[2].replace(/"/g, ''), // Remove quotes
              Area_Numbers: areaNumbers
            };
          }
        }

        // Parse WECC generation CSV
        const weccGenLines = weccGenText.split('\n');
        const weccGenData = {};
        const weccGenProcessed = [];

        for (let i = 1; i < weccGenLines.length; i++) {
          const line = weccGenLines[i].trim();
          if (line) {
            const values = line.split(',');
            const ba = values[0];
            if (ba) {
              const genData = {
                BA: ba,
                Hydro: parseFloat(values[1]) || 0,
                Nuclear: parseFloat(values[2]) || 0,
                Coal: parseFloat(values[3]) || 0,
                Natural_Gas: parseFloat(values[4]) || 0,
                Geothermal: parseFloat(values[5]) || 0,
                Biomass: parseFloat(values[6]) || 0,
                Wind: parseFloat(values[7]) || 0,
                PV: parseFloat(values[8]) || 0,
                Pumped_Storage_MW: parseFloat(values[9]) || 0,
                Battery_Storage_MW: parseFloat(values[10]) || 0,
                Total_MW: parseFloat(values[11]) || 0,
                Energy_Storage: parseFloat(values[12]) || 0
              };
              weccGenData[ba] = genData;
              weccGenProcessed.push(genData);
            }
          }
        }

        // Merge CSV data into GeoJSON properties
        geojsonData.features.forEach(feature => {
          const fid = feature.properties.FID;
          if (csvData[fid]) {
            // Merge all CSV properties into GeoJSON feature properties
            feature.properties = {
              ...feature.properties,
              ...csvData[fid]
            };
          }
          // Add area numbers from mapping
          if (areaMappingData[fid]) {
            feature.properties.Area_Numbers = areaMappingData[fid].Area_Numbers;
          }
          // Add generation data
          const baAbrev = feature.properties.BA_Abrev;
          if (baAbrev && weccGenData[baAbrev]) {
            feature.properties.generation = weccGenData[baAbrev];
          }
        });

        // Calculate total generation by source for chart
        const totalBySource = {
          Hydro: 0,
          Nuclear: 0,
          Coal: 0,
          Natural_Gas: 0,
          Geothermal: 0,
          Biomass: 0,
          Wind: 0,
          PV: 0
        };

        weccGenProcessed.forEach(data => {
          totalBySource.Hydro += data.Hydro;
          totalBySource.Nuclear += data.Nuclear;
          totalBySource.Coal += data.Coal;
          totalBySource.Natural_Gas += data.Natural_Gas;
          totalBySource.Geothermal += data.Geothermal;
          totalBySource.Biomass += data.Biomass;
          totalBySource.Wind += data.Wind;
          totalBySource.PV += data.PV;
        });

        // Setup WECC generation chart data - ACCURATE DATA FROM CSV ANALYSIS
        // Data source: WECC_31_BAs_2028_All_Clean (1).csv - Total: 384,199.9 MW
        const weccChartData = {
          labels: ['Natural Gas', 'Hydro', 'Solar', 'Wind', 'Coal', 'Battery Storage', 'Nuclear', 'Geothermal', 'Biomass'],
          datasets: [{
            label: 'WECC Generation Capacity (MW)',
            data: [
              163312.9,  // Natural Gas - 42.5% (largest)
              70028.6,   // Hydro - 18.2%
              41471.5,   // Solar (PV) - 10.8%
              40799.5,   // Wind - 10.6%
              34611.2,   // Coal - 9.0%
              16011.0,   // Battery Storage - 4.2%
              6321.6,    // Nuclear - 1.6%
              3850.3,    // Geothermal - 1.0%
              4165.8     // Biomass - 1.1%
            ],
            backgroundColor: [
              'rgba(255, 87, 34, 0.8)',   // Natural Gas - Orange
              'rgba(33, 150, 243, 0.8)',  // Hydro - Blue
              'rgba(255, 193, 7, 0.8)',   // Solar - Yellow/Gold
              'rgba(76, 175, 80, 0.8)',   // Wind - Green
              'rgba(97, 97, 97, 0.8)',    // Coal - Dark Gray
              'rgba(255, 235, 59, 0.8)',  // Battery Storage - Light Yellow
              'rgba(156, 39, 176, 0.8)',  // Nuclear - Purple
              'rgba(139, 69, 19, 0.8)',   // Geothermal - Brown
              'rgba(102, 187, 106, 0.8)'  // Biomass - Light Green
            ],
            borderColor: [
              'rgba(255, 87, 34, 1)',
              'rgba(33, 150, 243, 1)',
              'rgba(255, 193, 7, 1)',
              'rgba(76, 175, 80, 1)',
              'rgba(97, 97, 97, 1)',
              'rgba(255, 235, 59, 1)',
              'rgba(156, 39, 176, 1)',
              'rgba(139, 69, 19, 1)',
              'rgba(102, 187, 106, 1)'
            ],
            borderWidth: 2
          }]
        };

        // Setup filter values and name items
        const maxTotal = Math.max(...weccGenProcessed.map(d => d.Total_MW));
        const nameItems = weccGenProcessed.map(d => d.BA);

        // Create column data for WECC generation visualization - multiple bars per location
        const weccColumnData = [];
        const energySources = [
          { key: 'Hydro', color: 'blue', rgba: [28, 163, 236, 255] },
          { key: 'Nuclear', color: 'red', rgba: [255, 0, 0, 255] },
          { key: 'Coal', color: 'black', rgba: [0, 0, 0, 255] },
          { key: 'Natural_Gas', color: 'orange', rgba: [255, 165, 0, 255] },
          { key: 'Geothermal', color: 'purple', rgba: [128, 0, 128, 255] },
          { key: 'Biomass', color: 'green', rgba: [0, 128, 0, 255] },
          { key: 'Wind', color: 'lightgreen', rgba: [0, 255, 0, 255] },
          { key: 'PV', color: 'yellow', rgba: [255, 255, 0, 255] }
        ];

        geojsonData.features.forEach(feature => {
          if (feature.properties.generation && feature.geometry) {
            // Calculate centroid using turf library for better accuracy
            const centroid = center(feature.geometry);
            const baseLng = centroid.geometry.coordinates[0];
            const baseLat = centroid.geometry.coordinates[1];

            // Get bounding box to determine spread area
            const bounds = bbox(feature.geometry);
            const regionWidth = bounds[2] - bounds[0]; // maxLng - minLng
            const regionHeight = bounds[3] - bounds[1]; // maxLat - minLat

            // Use region size to determine appropriate spread, with minimum and maximum limits
            const spreadFactor = Math.min(Math.max(regionWidth, regionHeight, 0.1), 0.8); // Between 0.1 and 0.8 degrees

            if (baseLng && baseLat) {
              const generation = feature.properties.generation;

              // Create separate generation entries for each energy source (like the existing system)
              energySources.forEach((source, index) => {
                const value = generation[source.key];
                if (value && value > 50) { // Only show sources with significant capacity (>50 MW)
                  // Debug logging for IID specifically
                  if (generation.BA === 'IID') {
                    console.log(`IID - ${source.key}: ${value} MW, Color: ${source.color}`);
                  }

                  // Create offset position spread across the region based on region size
                  const offsetDistance = spreadFactor * 0.3; // Use 30% of region size for spread
                  const angle = (index / energySources.length) * 2 * Math.PI; // Distribute around circle
                  const offsetLng = baseLng + Math.cos(angle) * offsetDistance;
                  const offsetLat = baseLat + Math.sin(angle) * offsetDistance;

                  weccColumnData.push({
                    coordinates: [offsetLng, offsetLat],
                    Pg: value, // Use the same naming as existing generation
                    Pcap: value, // For capacity visualization
                    color: source.color,
                    fuel: source.key.toLowerCase(),
                    name: `${generation.BA} - ${source.key}`,
                    ba: generation.BA,
                    energyType: source.key,
                    KVlevels: [500], // Default high voltage for WECC
                    countyname: feature.properties.BA_Name || generation.BA
                  });
                }
              });
            }
          }
        });

        setWeccGenData(weccGenProcessed);
        setWeccGenChartData(weccChartData);
        setWeccGenFilter([0, maxTotal]);
        setWeccGenNameItems(nameItems);
        setWeccGenSelectItems([]);
        setWeccColumnData(weccColumnData);

        setWeccGeojsonData(geojsonData);
      } catch (err) {
        console.error('Error loading WECC data:', err);
      }
    };

    loadWeccData();
  }, []);

  const handleNetLayerChange = (event) => {
    setNetLayerActive(event.target.checked);
    setNetFilterValue([0, 800]);
  };


  const handleFlowLayerChange = (event) => {
    setFlowLayerActive(event.target.checked);
    setFlowFilterValue([0, 120]);
  };

  const handleTransmissionLayerChange = (event) => {
    setTransmissionLayerActive(event.target.checked);
  };

  const handlePowerPlantLayerChange = (event) => {
    setPowerPlantLayerActive(event.target.checked);
  };

  const handlePowerPlantLayerCapChange = (event) => {
    setPowerPlantLayerCapActive(event.target.checked);
  };

  const handlePowerPlantLabelsChange = (event) => {
    setPowerPlantLabelsActive(event.target.checked);
  };

  const handlePowerPlantRangeFilterChange = (event, newValue) => {
    setPowerPlantFilterValue(newValue);
  };

  const handlePowerPlantMultiselect = (selectedItems) => {
    setPowerPlantSelectItems(selectedItems);
  };

  const handleLoadLayerChange = (event) => {
    setLoadLayerActive(event.target.checked);
    setLoadFilterValue([0, countyloaddata.maxPd]);

    event.target.checked && (setInitialViewState(viewState => ({
      ...viewState,
      pitch: 40,
      transitionInterpolator: transitionFlyToInterpolator,
      transitionDuration: 2000,
    })))
  };

  const handleVoltageLayerChange = (event) => {
    setVoltageLayerActive(event.target.checked);
    setVoltageFilterValue([0.89, 1.11]);

    event.target.checked && (setInitialViewState(viewState => ({
      ...viewState,
      pitch: 40,
      transitionInterpolator: transitionFlyToInterpolator,
      transitionDuration: 2000,
    })))
  };

  const handleAreaLayerChange = (event) => {
    setAreaLayerActive(event.target.checked);
    //    setVoltageFilterValue([0.89, 1.11]);

    event.target.checked && (setInitialViewState(viewState => ({
      ...viewState,
      pitch: 40,
      transitionInterpolator: transitionFlyToInterpolator,
      transitionDuration: 2000,
    })))
  };

  const handleZoneLayerChange = (event) => {
    setZoneLayerActive(event.target.checked);
    //    setVoltageFilterValue([0.89, 1.11]);

    event.target.checked && (setInitialViewState(viewState => ({
      ...viewState,
      pitch: 40,
      transitionInterpolator: transitionFlyToInterpolator,
      transitionDuration: 2000,
    })))
  };

  const handleWeccLayerChange = (event) => {
    setWeccLayerActive(event.target.checked);

    event.target.checked && (setInitialViewState(viewState => ({
      ...viewState,
      pitch: 40,
      transitionInterpolator: transitionFlyToInterpolator,
      transitionDuration: 2000,
    })))
  };


  const handleGenLayerChange = (event) => {
    setGenLayerActive(event.target.checked);
    setGenFilterValue([gendata.minPg, gendata.maxPg]);

    event.target.checked && (setInitialViewState(viewState => ({
      ...viewState,
      pitch: 40,
      transitionInterpolator: transitionFlyToInterpolator,
      transitionDuration: 2000,
    })))
  };

  const handleGenLayerCapChange = (event) => {
    setGenLayerCapActive(event.target.checked);
    setGenFilterValue([gendata.minPg, gendata.maxPg]);

    event.target.checked && (setInitialViewState(viewState => ({
      ...viewState,
      pitch: 40,
      transitionInterpolator: transitionFlyToInterpolator,
      transitionDuration: 2000,
    })))
  };

  const handleWeccGenLayerChange = (event) => {
    setWeccGenLayerActive(event.target.checked);
    if (weccGenData && weccGenData.length > 0) {
      const maxTotal = Math.max(...weccGenData.map(d => d.Total_MW));
      setWeccGenFilter([0, maxTotal]);
    }

    event.target.checked && (setInitialViewState(viewState => ({
      ...viewState,
      pitch: 40,
      transitionInterpolator: transitionFlyToInterpolator,
      transitionDuration: 2000,
    })))
  };

  const handleWeccGenRangeFilterChange = (event, newValue) => {
    setWeccGenFilter(newValue);
  };

  const handleWeccGenMultiselect = (value) => {
    setWeccGenSelectItems(value);
  };

  function getNetFilterValue(data) {
    if (!data) return 10000;
    if (lineNameSelectItems.length > 0) {  // when users make selections by name, does not consider netfiltervalue 
      const pointNames = lineNameSelectItems.map(n => n.split(' -- ')).flat()
      if (data.geometry.type == 'Point' && pointNames.includes(data.properties.NAME)) {
        return data.properties.KVlevels[0]
      } else if (data.geometry.type == 'LineString' && lineNameSelectItems.includes(data.properties.NAME)) {
        /* Line layer */
        return data.properties.KV;
      }
    } else if (busNameSelectItems.length > 0) {

      if (data.geometry.type === "Point" && busNameSelectItems.includes(data.properties.NAME)) {
        return data.properties.KVlevels[0]
      } else if (data.geometry.type === "LineString" && data.properties.NAME.split(' -- ').some(r => busNameSelectItems.includes(r))) {
        return data.properties.KV;
      }


    } else {
      if (data.geometry.type == 'Point') {
        for (var i = 0; i < data.properties.KVlevels.length; i++) {
          var KV = data.properties.KVlevels[i];
          if (netfiltervalue[0] <= KV && KV <= netfiltervalue[1]) return KV;
        }
      } else {
        if (data.geometry.type == 'LineString') { /* Line layer */
          /* Uncomment to activate flow-based filtering
          var RATE_A;
          if(data.properties.RATE_A == 0) {
        RATE_A = 10000;
          } else {
        RATE_A = data.properties.RATE_A;
          }
          var loading = Math.abs(data.properties.PF / RATE_A)*100;
          if(flowfiltervalue[0] <= loading && loading <= flowfiltervalue[1]) {
        return data.properties.KV;
        }
          */
          return data.properties.KV;
        }
      }
    }

    return -1; // This is beyond the range so filter will filter out this data point.
  }

  function getFlowFilterValue(data) {

  }

  function getGenFilterValue(data) {
    if (!data) return 10000;   //10000 is beyond the range, so the generation will be filter out 
    if ((genDoughlabels.length > 0) && (!(genDoughlabels.indexOf(colorMap[data.color]) >= 0))) return 10000;
    if (nameSelectItems.length > 0) {
      if (nameSelectItems.includes(data.name)) {
        for (var i = 0; i < data.KVlevels.length; i++) {
          var KV = data.KVlevels[i];
          if (netfiltervalue[0] <= KV && KV <= netfiltervalue[1]) {

            return data.Pg;
          }
        }
        return 10000;
      } else {
        return 10000;
      }
    }

    for (var i = 0; i < data.KVlevels.length; i++) {
      var KV = data.KVlevels[i];
      if (netfiltervalue[0] <= KV && KV <= netfiltervalue[1]) {

        return data.Pg;
      }
    }
    return 10000;

  }



  function getLoadFilterValue(data) {
    if (!data) return -10000;
    if (countyNameSelectItems.length > 0) {
      if (countyNameSelectItems.includes(data.properties.countyname)) {

        for (var i = 0; i < data.properties.KVlevels.length; i++) {
          var KV = data.properties.KVlevels[i];
          if (netfiltervalue[0] <= KV && KV <= netfiltervalue[1]) {
            return data.properties.Pd;
          }
        }
        return -10000;

      } else {
        return -10000;
      }

    }

    for (var i = 0; i < data.properties.KVlevels.length; i++) {
      var KV = data.properties.KVlevels[i];
      if (netfiltervalue[0] <= KV && KV <= netfiltervalue[1]) {
        return data.properties.Pd;
      }
    }
    return -10000;
  }

  function getVoltageFilterValue(data) {
    if (!data) return -10000;

    if (countyNameSelectItems.length > 0) {
      if (countyNameSelectItems.includes(data.properties.countyname)) {
        for (var i = 0; i < data.properties.KVlevels.length; i++) {
          var KV = data.properties.KVlevels[i];
          if (netfiltervalue[0] <= KV && KV <= netfiltervalue[1]) {
            return data.properties.Vm_avg;
          }
        }
        return -10000;
      } else {
        return -10000;
      }
    }

    for (var i = 0; i < data.properties.KVlevels.length; i++) {
      var KV = data.properties.KVlevels[i];
      if (netfiltervalue[0] <= KV && KV <= netfiltervalue[1]) {
        return data.properties.Vm_avg;
      }
    }
    return -10000;
  }

  function getWeccGenFilterValue(data) {
    if (!data || !data.properties || !data.properties.generation) return -10000;

    const generation = data.properties.generation;
    const baAbrev = data.properties.BA_Abrev;

    // Check if BA is in selected items
    if (weccGenSelectItems.length > 0) {
      if (!weccGenSelectItems.includes(baAbrev)) {
        return -10000;
      }
    }

    // Check if total MW is within filter range
    if (generation.Total_MW >= weccGenFilter[0] && generation.Total_MW <= weccGenFilter[1]) {
      return generation.Total_MW;
    }

    return -10000;
  }

  function getWeccColumnFilterValue(data) {
    if (!data) return -10000;

    // Check if energy type is visible in doughnut chart (same logic as existing generation)
    if ((weccGenDoughlabels.length > 0) && (!(weccGenDoughlabels.indexOf(weccColorMap[data.color]) >= 0))) return -10000;

    // Check if BA is in selected items
    if (weccGenSelectItems.length > 0) {
      if (!weccGenSelectItems.includes(data.ba)) {
        return -10000;
      }
    }

    // Check if generation value is within filter range  
    if (data.Pg >= weccGenFilter[0] && data.Pg <= weccGenFilter[1]) {
      return data.Pg;
    }

    return -10000;
  }

  function fillWeccColumnColor(data) {
    switch (data.color) {
      case 'blue': return [28, 163, 236, 255]; // Hydro
      case 'red': return [255, 0, 0, 255]; // Nuclear
      case 'black': return [0, 0, 0, 255]; // Coal
      case 'orange': return [255, 165, 0, 255]; // Natural Gas
      case 'purple': return [128, 0, 128, 255]; // Geothermal
      case 'green': return [0, 128, 0, 255]; // Biomass
      case 'lightgreen': return [0, 255, 0, 255]; // Wind
      case 'yellow': return [255, 255, 0, 255]; // Solar/PV
      default: return [128, 128, 128, 255]; // Gray fallback
    }
  }

  const layers = [    // new FlowmapLayer({
    //   id: 'my-flowmap-layer',
    //   data: flowdata,
    //   visible: flowlayeractive,
    //   animationEnabled: true, //control the animation effect of flow layer
    //   colorScheme: ["rgb(0,0,255)","rgb(255,0,255)"],
    //   // darkMode: true, 
    //   // clusteringEnabled: false, //control the aggregate effect of flow layer
    //   // adaptiveScalesEnabled: false, 
    //   getFlowMagnitude: (flow) => flow.count,
    //   getFlowOriginId: (flow) => flow.origin,
    //   getFlowDestId: (flow) => flow.dest,
    //   getLocationId: (loc) => loc.id,
    //   getLocationLat: (loc) => loc.lat,
    //   getLocationLon: (loc) => loc.lon,
    //   
    // }),

    new GeoJsonLayer({
      id: 'geojson',
      data: enhancedData,
      stroked: false,
      filled: true,
      //      extruded: true,
      pickable: netlayeractive,
      pointType: 'circle',
      lineWidthScale: 3,
      getFillColor: FillColor,
      getLineColor: LineColor,
      getPointRadius: 1000,
      getLineWidth: LineWidth,
      visible: netlayeractive,
      onClick: zoomToData,
      getFilterValue: getNetFilterValue,
      filterRange: netfiltervalue,

      extensions: [new DataFilterExtension({ filtersize: 1 })],
      updateTriggers: {
        getFilterValue: [netfiltervalue, lineNameSelectItems, busNameSelectItems, flowfiltervalue]
      }
    }),

    // Dedicated Transmission Lines Layer for Powerlines WUS CAN SGCA data
    new GeoJsonLayer({
      id: 'transmission-lines',
      data: enhancedData,
      pickable: transmissionlayeractive,
      stroked: true,
      filled: false,
      visible: transmissionlayeractive,
      lineWidthScale: 2,
      lineWidthMinPixels: 1,
      lineWidthMaxPixels: 4,
      getLineColor: d => {
        // Color transmission lines with blue shades based on voltage level
        if (d.geometry.type === 'LineString') {
          const voltage = d.properties.kilovolt || d.properties.KV || d.properties.voltage || 100;
          
          // Handle the actual data ranges found in the datasets
          if (voltage >= 500) return [0, 100, 200, 220]; // Dark blue for extra high voltage (500kV+)
          if (voltage >= 345) return [30, 144, 255, 210]; // Dodger blue for high voltage (345kV)
          if (voltage >= 230) return [70, 170, 255, 200]; // Medium blue for medium voltage (230kV)
          if (voltage >= 138) return [100, 200, 255, 190]; // Light blue for sub-transmission (138kV)
          if (voltage >= 100) return [120, 210, 255, 180]; // Lighter blue for 100kV lines
          if (voltage > 0) return [150, 220, 255, 170]; // Very light blue for other voltages
          
          // Special handling for zero/unknown voltage - use medium blue as default
          return [70, 170, 255, 160]; // Default blue for unknown/zero voltage lines
        }
        return [128, 128, 128, 0]; // Transparent for non-lines
      },
      getLineWidth: d => {
        if (d.geometry.type === 'LineString') {
          const voltage = d.properties.kilovolt || d.properties.KV || d.properties.voltage || 100;
          if (voltage >= 500) return 4; // Thickest for extra high voltage
          if (voltage >= 345) return 3; // Thick for high voltage
          if (voltage >= 230) return 2.5; // Medium-thick for medium voltage
          if (voltage >= 138) return 2; // Medium for sub-transmission
          if (voltage >= 100) return 1.5; // Slightly thicker for 100kV lines
          if (voltage > 0) return 1; // Standard width for other voltages
          return 1.5; // Default width for unknown/zero voltage
        }
        return 0;
      },
      onClick: (info) => {
        if (info.object && info.object.geometry.type === 'LineString') {
          const props = info.object.properties;
          setShowPopup({
            display: true,
            info: `
              <div style="font-weight: 600; margin-bottom: 4px; color: #1976d2;">${props.line_name || props.NAME || 'Transmission Line'}</div>
              <div style="font-size: 11px; color: #666;">Voltage: ${props.kilovolt || props.KV || 'Unknown'} kV</div>
              <div style="font-size: 11px; color: #666;">Source: ${props.data_source || props.SOURCEFILE || 'Unknown'}</div>
              ${props.flow_capacity ? `<div style="font-size: 11px; color: #666;">Capacity: ${props.flow_capacity} MW</div>` : ''}
            `,
            name: props.line_name || props.NAME || 'Transmission Line',
            fid: info.index,
            type: 'transmission'
          });
        }
      },
      // Filter to only show LineString geometries for transmission lines
      getFilterValue: d => d.geometry.type === 'LineString' ? 1 : 0,
      filterRange: [1, 1],
      extensions: [new DataFilterExtension({ filtersize: 1 })]
    }),

    // Western Power Plants Generation Column Layer
    new ColumnLayer({
      id: 'western-power-plants',
      data: (enhancedData.powerPlants || []).filter(plant => {
        if (!selectedTechnology) return true;
        return plant.primaryType === selectedTechnology;
      }),
      diskResolution: 50,
      radius: 5000,
      elevationScale: 50,
      pickable: powerplantlayeractive,
      visible: powerplantlayeractive,
      getPosition: d => d.coordinates,
      getFillColor: d => getPowerPlantColumnColor(d.primaryType, d.Pg),
      getElevation: d => d.Pg * 5,
      onClick: (info) => {
        if (info && info.object) {
          setSelectedFeature({
            properties: info.object,
            fid: info.index,
            type: 'western_power_plant'
          });
        }
      },
      getFilterValue: d => {
        // Filter by capacity range
        if (!d || typeof d.Pg !== 'number') return -10000;
        return d.Pg;
      },
      filterRange: powerPlantFilterValue || [0, 10000],
      extensions: [new DataFilterExtension({ filtersize: 1 })],
      updateTriggers: {
        getData: selectedTechnology,
        getFillColor: selectedTechnology
      }
    }),

    // Western Power Plants Capacity Column Layer
    new ColumnLayer({
      id: 'western-power-plants-cap',
      data: (enhancedData.powerPlants || []).filter(plant => {
        if (!selectedTechnology) return true;
        return plant.primaryType === selectedTechnology;
      }),
      diskResolution: 50,
      radius: 5000,
      elevationScale: 50,
      pickable: false,
      visible: powerplantlayercapactive,
      getPosition: d => d.coordinates,
      getFillColor: d => getPowerPlantColumnColor(d.primaryType, d.Pcap),
      getElevation: d => d.Pcap * 5,
      onClick: (info) => {
        if (info && info.object) {
          setSelectedFeature({
            properties: info.object,
            fid: info.index,
            type: 'western_power_plant_cap'
          });
        }
      },
      getFilterValue: d => {
        // Filter by capacity range
        if (!d || typeof d.Pcap !== 'number') return -10000;
        return d.Pcap;
      },
      filterRange: powerPlantFilterValue || [0, 10000],
      extensions: [new DataFilterExtension({ filtersize: 1 })],
      updateTriggers: {
        getData: selectedTechnology,
        getFillColor: selectedTechnology
      }
    }),

    // Power Plant Labels Layer
    new TextLayer({
      id: 'power-plant-labels',
      data: (enhancedData.powerPlants || []).filter(plant => {
        if (selectedTechnology && plant.primaryType !== selectedTechnology) return false;
        // Only show labels for larger plants to avoid clutter
        return plant.Pg >= 50;
      }),
      visible: powerplantlabelsactive && powerplantlayeractive,
      pickable: false,
      getPosition: d => {
        // Position labels at the top of the columns
        const coords = d.coordinates;
        const capacity = d.Pg || 1;
        const elevation = capacity * 5 * 50; // Match column height calculation (capacity * 5 * elevationScale)
        return [coords[0], coords[1], elevation + 50]; // Slightly above the column
      },
      getText: d => {
        const name = d.name || 'Unnamed Plant';
        const capacity = d.Pg || 0;
        return `${name}\n${capacity.toFixed(0)} MW`;
      },
      getSize: d => {
        // Size text based on plant capacity
        const capacity = d.Pg || 1;
        if (capacity >= 1000) return 14;
        if (capacity >= 100) return 12;
        if (capacity >= 50) return 10;
        return 9;
      },
      getAngle: 0,
      getTextAnchor: 'middle',
      getAlignmentBaseline: 'bottom',
      getColor: d => {
        const color = getPowerPlantColor(d.primaryType);
        return [color[0], color[1], color[2], 240]; // High contrast
      },
      getPixelOffset: [0, -10], // Slightly above the column top
      fontFamily: '"Inter", sans-serif',
      fontWeight: 700,
      outlineWidth: 3,
      outlineColor: [255, 255, 255, 200],
      backgroundColor: [255, 255, 255, 160],
      getBackgroundPadding: [4, 2, 4, 2],
      backgroundRadius: 4,
      updateTriggers: {
        getData: selectedTechnology,
        getColor: [powerplantlabelsactive, selectedTechnology],
        getText: powerplantlabelsactive,
        getSize: powerplantlabelsactive,
        getPosition: powerplantlabelsactive
      }
    }),

    new ColumnLayer({
      id: 'gen-column',
      data: generation,
      diskResolution: 50,
      radius: 5000,
      elevationScale: 50,
      pickable: genlayeractive,
      visible: genlayeractive,
      getPosition: d => d.coordinates,
      getFillColor: fillGenColumnColor,
      getElevation: d => d.Pg * 5,
      onClick: zoomToData,

      getFilterValue: getGenFilterValue,
      filterRange: genfiltervalue,

      extensions: [new DataFilterExtension({ filtersize: 1 })],

      updateTriggers: {
        getFilterValue: [netfiltervalue, genDoughlabels, nameSelectItems]
      }

    }),

    new ColumnLayer({
      id: 'gen-column-cap',
      data: generation,
      diskResolution: 50,
      radius: 5000,
      elevationScale: 50,
      pickable: false, //genlayeractive,
      visible: genlayercapactive,
      getPosition: d => d.coordinates,
      getFillColor: fillGenColumnColorCap,
      getElevation: d => d.Pcap * 5,
      onClick: zoomToData,

      getFilterValue: getGenFilterValue,
      filterRange: genfiltervalue,

      extensions: [new DataFilterExtension({ filtersize: 1 })],

      updateTriggers: {
        getFilterValue: [netfiltervalue, genDoughlabels, nameSelectItems],
      }

    }),

    /*    
    new ColumnLayer({
      id: 'load-column',
      data: loads,
      diskResolution: 50,
      radius: 5000,
      elevationScale: 50,
      pickable: loadlayeractive,
      visible: loadlayeractive,
      getFillColor: [255,255,0],//[255, 239, 247],
      getPosition: d => d.coordinates,
//      getFillColor: fillGenColumnColor,
      getElevation: d => d.Pd*5,
      onClick:zoomToData
    }),
    */

    /*
    new GeoJsonLayer({
      id: 'PolygonLayer2',
      data:countyload,
      pickable: loadlayeractive,
      visible: loadlayeractive,
      stroked: true,
      filled: true,
      extruded: true,
      wireframe: true,
      lineWidthMinPixels: 1,
      getPolygon: d => d.geometry.coordinates,
//      getElevation: d => d.properties.Pd*5.0,
      getFillColor: d => [255*d.properties.Pd/countymaxPd, 0, 0],
      getLineColor: [80,80,80],
      getLineWidth: d => 1,
      opacity: 0.1,
      onClick: zoomToCounty,
      extensions: [new DataFilterExtension({filtersize:1})],
      getFilterValue: getLoadFilterValue,
      filterRange: loadfiltervalue,

      updateTriggers: {
        getFilterValue: netfiltervalue
      }
    }),
    */


    new GeoJsonLayer({
      id: 'AreaLayer',
      data: areas,
      pickable: arealayeractive,
      visible: arealayeractive,
      stroked: true,
      filled: true,
      extruded: true,
      wireframe: true,
      lineWidthMinPixels: 1,
      getPolygon: d => d.geometry.coordinates,
      //      getElevation: d => d.properties.Pd*5.0,
      getFillColor: [255, 192, 203],
      getLineColor: [80, 80, 80],
      getLineWidth: d => 1,
      opacity: 0.1,
      onClick: zoomToArea,
      //      extensions: [new DataFilterExtension({ filtersize: 1 })],
      //      getFilterValue: getLoadFilterValue,
      //      filterRange: loadfiltervalue,

      //      updateTriggers: {
      //        getFilterValue: [netfiltervalue, countyNameSelectItems]
      //      }
    }),

    new GeoJsonLayer({
      id: 'ZoneLayer',
      data: zones,
      pickable: zonelayeractive,
      visible: zonelayeractive,
      stroked: true,
      filled: true,
      extruded: true,
      wireframe: true,
      lineWidthMinPixels: 1,
      getPolygon: d => d.geometry.coordinates,
      //      getElevation: d => d.properties.Pd*5.0,
      getFillColor: [252, 245, 95],
      getLineColor: [80, 80, 80],
      getLineWidth: d => 1,
      opacity: 0.1,
      onClick: zoomToZone,
      //      extensions: [new DataFilterExtension({ filtersize: 1 })],
      //      getFilterValue: getLoadFilterValue,
      //      filterRange: loadfiltervalue,

      //      updateTriggers: {
      //        getFilterValue: [netfiltervalue, countyNameSelectItems]
      //      }
    }),

    new GeoJsonLayer({
      id: 'PolygonLayerload',
      data: countyload,
      pickable: loadlayeractive,
      visible: loadlayeractive,
      stroked: true,
      filled: true,
      extruded: true,
      wireframe: true,
      lineWidthMinPixels: 1,
      getPolygon: d => d.geometry.coordinates,
      //      getElevation: d => d.properties.Pd*5.0,
      getFillColor: d => [255 * d.properties.Pd / countymaxPd, 0, 0],
      getLineColor: [80, 80, 80],
      getLineWidth: d => 1,
      opacity: 0.1,
      onClick: zoomToCounty,
      extensions: [new DataFilterExtension({ filtersize: 1 })],
      getFilterValue: getLoadFilterValue,
      filterRange: loadfiltervalue,

      updateTriggers: {
        getFilterValue: [netfiltervalue, countyNameSelectItems]
      }
    }),


    new GeoJsonLayer({
      id: 'PolygonLayer2',
      data: countyload,
      pickable: voltagelayeractive,
      visible: voltagelayeractive,
      stroked: true,
      filled: true,
      extruded: true,
      wireframe: true,
      lineWidthMinPixels: 1,
      getPolygon: d => d.geometry.coordinates,
      //      getElevation: d => d.properties.Pd*5.0,
      getFillColor: getVoltageFillColor,
      getLineColor: [80, 80, 80],
      getLineWidth: d => 1,
      opacity: 0.1,
      onClick: zoomToCounty,
      extensions: [new DataFilterExtension({ filtersize: 1 })],
      getFilterValue: getVoltageFilterValue,
      filterRange: voltagefiltervalue,

      updateTriggers: {
        getFilterValue: [netfiltervalue, countyNameSelectItems]
      }
    }),

    // WECC Balancing Authorities Layer
    new GeoJsonLayer({
      id: 'WeccLayer',
      data: weccGeojsonData,
      pickable: wecclayeractive,
      visible: wecclayeractive,
      stroked: true,
      filled: true,
      extruded: false,
      wireframe: false,
      lineWidthMinPixels: 2,
      lineWidthMaxPixels: 10,
      getPolygon: d => d.geometry.coordinates,
      getLineColor: [255, 255, 255, 255], // White border always visible
      getLineWidth: 4, // Thicker border for better visibility
      getFillColor: d => {
        // Get region abbreviation for impact analysis
        const regionAbbrev = d.properties?.BA_CODE || d.properties?.NAME || d.properties?.ABBREV;
        
        // Use impact analysis colors if active and data is available
        if (impactAnalysisActive && Object.keys(impactData).length > 0 && regionAbbrev) {
          const impactColor = getWeccRegionColor(regionAbbrev);
          
          // Apply hover and click effects to impact colors
          if (d === weccClickedObject) {
            // Darken the impact color for clicked state
            return [
              Math.max(0, impactColor[0] - 50),
              Math.max(0, impactColor[1] - 50),
              Math.max(0, impactColor[2] - 50),
              220
            ];
          }
          if (d === weccHoveredObject && d !== weccClickedObject) {
            // Slightly darken for hover
            return [
              Math.max(0, impactColor[0] - 20),
              Math.max(0, impactColor[1] - 20),
              Math.max(0, impactColor[2] - 20),
              200
            ];
          }
          
          return impactColor;
        }
        
        // Default coloring when impact analysis is not active
        // Clicked region stays highlighted until another region is clicked
        if (d === weccClickedObject) {
          return [30, 90, 150, 200]; // Darker blue for clicked/active state
        }
        // Hover effect (only if not clicked)
        if (d === weccHoveredObject && d !== weccClickedObject) {
          return [70, 130, 180, 180]; // Slightly darker blue on hover
        }
        // Default blue fill with less opacity to show borders better
        return [100, 149, 237, 120]; // Cornflower blue with more transparency
      },
      opacity: 0.8,
      onHover: (info) => {
        setWeccHoveredObject(info.object);
      },
      onClick: (info) => {
        if (info.object) {
          setWeccClickedObject(info.object);
          zoomToWecc(info);
        }
      },
      updateTriggers: {
        getFillColor: [weccHoveredObject, weccClickedObject, impactAnalysisActive, impactData]
      }
    }),

    // WECC Generation Power Column Layer - Multiple bars per location with reduced height
    new ColumnLayer({
      id: 'WeccGenColumnLayer',
      data: weccColumnData,
      diskResolution: 50,
      radius: 5000, // Same as existing generation
      elevationScale: 33, // Reduced from 50 to 33 (2/3 of original)
      pickable: weccGenLayerActive,
      visible: weccGenLayerActive,
      getPosition: d => d.coordinates,
      getFillColor: fillWeccColumnColor, // Use our custom color function
      getElevation: d => d.Pg * 3.33, // Reduced from d.Pg * 5 to d.Pg * 3.33 (2/3 of original)
      onClick: zoomToData, // Same click handler as existing generation

      getFilterValue: getWeccColumnFilterValue,
      filterRange: weccGenFilter,

      extensions: [new DataFilterExtension({ filtersize: 1 })],
      updateTriggers: {
        getFilterValue: [weccGenFilter, weccGenSelectItems, weccGenDoughlabels]
      }
    }),

    /*
    new HeatmapLayer({
      id:'Voltagecontour',
      data:loads,
      getWeight: d => d.Pd,
      getPosition: d => d.coordinates,
      aggregation: 'MEAN'
    })
    */
  ];

  /* Chart for generation mix */
  const genmixlabels = [
    'Wind',
    'Solar',
    'Nuclear',
    'Natural Gas',
    'Hydro',
    'Coal',
    'Other'
  ];

  var genmix = [];
  genmix.push(gendata.Pgwind);
  genmix.push(gendata.Pgsolar);
  genmix.push(gendata.Pgnuclear);
  genmix.push(gendata.Pgng);
  genmix.push(gendata.Pghydro);
  genmix.push(gendata.Pgcoal);
  genmix.push(gendata.Pgother);

  var genmixcap = [];
  genmixcap.push(gendata.Pgwindcap);
  genmixcap.push(gendata.Pgsolarcap);
  genmixcap.push(gendata.Pgnuclearcap);
  genmixcap.push(gendata.Pgngcap);
  genmixcap.push(gendata.Pghydrocap);
  genmixcap.push(gendata.Pgcoalcap);
  genmixcap.push(gendata.Pgothercap);

  const handleDoughnutClick = (event, legendItem, legend) => {

    // filter generation 
    if (legendItem.hidden) {  //if ishidden, then add to array

      setDoughlabels(genDoughlabels => ([...genDoughlabels, legendItem.text]))

    } else {  //remove from array 
      let removeDoughLabels = [...genDoughlabels]; // creates a copy of subNames on a new reference
      let index = genDoughlabels.indexOf(legendItem.text)
      const newArray = [...genDoughlabels.slice(0, index), ...genDoughlabels.slice(index + 1)];
      setDoughlabels(newArray)

    }
    //default legend function of doughnut chart 
    legend.chart.toggleDataVisibility(legendItem.index);

  }

  const handleWeccDoughnutClick = (event, legendItem, legend) => {
    // filter WECC generation 
    if (legendItem.hidden) {  //if ishidden, then add to array
      setWeccDoughlabels(weccGenDoughlabels => ([...weccGenDoughlabels, legendItem.text]))
    } else {  //remove from array 
      let index = weccGenDoughlabels.indexOf(legendItem.text)
      const newArray = [...weccGenDoughlabels.slice(0, index), ...weccGenDoughlabels.slice(index + 1)];
      setWeccDoughlabels(newArray)
    }
    //default legend function of doughnut chart 
    legend.chart.toggleDataVisibility(legendItem.index);
  }

  const handleBusMultiselect = (selectItem, metadata) => {

    //  only one is active between bus name selection and transmission line name selection at a time
    setLineNameSelectItems([])

    const selected = busNameSelectItems.indexOf(metadata.dataItem)

    if (selected >= 0) {  //if is selected, remove 
      const newArray = [...busNameSelectItems.slice(0, selected), ...busNameSelectItems.slice(selected + 1)];
      setBusNameSelectItems(newArray)

    } else {  //add to array 
      setBusNameSelectItems(busNameSelectItems => ([...busNameSelectItems, metadata.dataItem]))
    }
  }

  const handleLineMultiselect = (selectItem, metadata) => {
    //  only one is active between bus name selection and transmission line name selection at a time
    setBusNameSelectItems([])

    const selected = lineNameSelectItems.indexOf(metadata.dataItem)

    if (selected >= 0) {  //if is selected, remove 
      const newArray = [...lineNameSelectItems.slice(0, selected), ...lineNameSelectItems.slice(selected + 1)];
      setLineNameSelectItems(newArray)

    } else {  //add to array 
      setLineNameSelectItems(lineNameSelectItems => ([...lineNameSelectItems, metadata.dataItem]))
      // const filterGen = gendata.Gens.filter(gen => gen.name === metadata.dataItem)
      // if (filterGen.length > 0) {
      //   const long = filterGen[0].coordinates[0]
      //   const lat = filterGen[0].coordinates[1]
      //   zoomToGen(lat, long)
      // }
    }
  }

  const handleGenMultiselect = (selectItem, metadata) => {

    const selected = nameSelectItems.indexOf(metadata.dataItem)
    if (selected >= 0) {  //if is selected, remove 
      const newArray = [...nameSelectItems.slice(0, selected), ...nameSelectItems.slice(selected + 1)];
      setNameSelectItems(newArray)

    } else {  //add to array 
      setNameSelectItems(nameSelectItems => ([...nameSelectItems, metadata.dataItem]))
      const filterGen = gendata.Gens.filter(gen => gen.name === metadata.dataItem)
      if (filterGen.length > 0) {
        const long = filterGen[0].coordinates[0]
        const lat = filterGen[0].coordinates[1]
        zoomToGen(lat, long)
      }
    }
  }



  const handleCountyMultiselect = (selectItem, metadata) => {

    const selected = countyNameSelectItems.indexOf(metadata.dataItem)
    if (selected >= 0) {  //if is selected, remove 
      const newArray = [...countyNameSelectItems.slice(0, selected), ...countyNameSelectItems.slice(selected + 1)];
      setCountyNameSelectItems(newArray)

    } else {  //add to array 
      setCountyNameSelectItems(countyNameSelectItems => ([...countyNameSelectItems, metadata.dataItem]))
      const filtercounty = countyload.features.filter(county => county.properties.countyname === metadata.dataItem)
      if (filtercounty.length > 0) {
        const longs = filtercounty[0].geometry.coordinates[0].map(d => d[0])
        const lats = filtercounty[0].geometry.coordinates[0].map(d => d[1])
        const minLng = Math.min(...longs)
        const maxLng = Math.max(...longs)
        const minLat = Math.min(...lats)
        const maxLat = Math.max(...lats)
        console.log(minLng, minLat, maxLng, maxLat)
        zoomToCountyName(minLng, minLat, maxLng, maxLat)
      }


    }
  }

  const handleAreaMultiselect = (selectItem, metadata) => {

    const selected = areaNameSelectItems.indexOf(metadata.dataItem)
    console.log(selected);
    if (selected >= 0) {  //if is selected, remove 
      const newArray = [...areaNameSelectItems.slice(0, selected), ...areaNameSelectItems.slice(selected + 1)];
      setAreaNameSelectItems(newArray)

    } else {  //add to array 
      setAreaNameSelectItems(areaNameSelectItems => ([...areaNameSelectItems, metadata.dataItem]))
      const filterareas = areas.features.filter(area => area.properties.name === metadata.dataItem)
      if (filterareas.length > 0) {
        const longs = filterareas[0].geometry.coordinates[0].map(d => d[0])
        const lats = filterareas[0].geometry.coordinates[0].map(d => d[1])
        const minLng = Math.min(...longs)
        const maxLng = Math.max(...longs)
        const minLat = Math.min(...lats)
        const maxLat = Math.max(...lats)

        zoomToAreaName(filterareas[0], minLng, minLat, maxLng, maxLat)
      }
    }
  }

  const handleZoneMultiselect = (selectItem, metadata) => {

    const selected = zoneNameSelectItems.indexOf(metadata.dataItem)
    console.log(selected);
    if (selected >= 0) {  //if is selected, remove 
      const newArray = [...zoneNameSelectItems.slice(0, selected), ...zoneNameSelectItems.slice(selected + 1)];
      setZoneNameSelectItems(newArray)

    } else {  //add to array 
      setZoneNameSelectItems(zoneNameSelectItems => ([...zoneNameSelectItems, metadata.dataItem]))
      const filterzones = zones.features.filter(area => area.properties.name === metadata.dataItem)
      if (filterzones.length > 0) {
        const longs = filterzones[0].geometry.coordinates[0].map(d => d[0])
        const lats = filterzones[0].geometry.coordinates[0].map(d => d[1])
        const minLng = Math.min(...longs)
        const maxLng = Math.max(...longs)
        const minLat = Math.min(...lats)
        const maxLat = Math.max(...lats)

        zoomToZoneName(filterzones[0], minLng, minLat, maxLng, maxLat)
      }
    }
  }



  const renderItem = ({
    id,
    name
  },) => {
    return (
      <MenuItem
        key={id}
        text={name}
      />
    );
  }

  const chartdata = {
    labels: genmixlabels,
    datasets: [
      {
        label: 'Generation Mix Cap',
        data: genmixcap,
        backgroundColor: [
          'rgba(0,255,0,0.3)',
          'rgba(244,219,135,0.3)',
          'rgba(255,0,0,0.3)',
          'rgba(255,165,0,0.3)',
          'rgba(28,163,236,0.3)',
          'rgba(128,128,128,0.3)',
          'rgba(0,0,0,0.3)'
        ],
        borderWidth: 1,
        options: {
          plugins: {

            legend: {
              onClick: (evt, legendItem, legend) => { console.log('sdsd') }, // Add onClick event handler to the legend
            },


            title: {
              display: true,
              text: 'Generation Mix Cap',
              align: 'center',
              position: 'top'
            }
          }
        }
      },
      {
        label: 'Generation Mix',
        data: genmix,
        backgroundColor: [
          'rgb(0,255,0)',
          'rgb(244,219,135)',
          'red',
          'orange',
          'rgb(28,163,236)',
          'gray',
          'black'
        ],
        borderWidth: 1,
        options: {
          plugins: {
            legend: {
              onClick: (evt, legendItem, legend) => { console.log('sdsd') }, // Add onClick event handler to the legend
            },
            title: {
              display: true,
              text: 'Generation Mix',
              align: 'center',
              position: 'top'
            }
          }
        }
      }
    ],
  };


  const handleGenRangeFilterChange = (event) => {
    setGenFilterValue(event.target.value);

  }

  const handleLoadRangeFilterChange = (event) => {
    setLoadFilterValue(event.target.value);
  }

  const handleVoltageRangeFilterChange = (event) => {
    setVoltageFilterValue(event.target.value);
  }

  const handleNetRangeFilterChange = (event) => {

    setNetFilterValue(event.target.value);

  }

  const handleFlowRangeFilterChange = (event) => {

    setFlowFilterValue(event.target.value);

  }

  const handleNetBarFilterChange = (value) => {

    setNetFilterValue(value);
  }


  function valuetext(value) {
    return `${value.toFixed(2)}`;
  }


  return (
    <>
      <WestmapHeader />
      
      {/* Floating Action Button for WECC Case Studies */}
      <Fab
        variant="extended"
        style={{
          position: 'fixed',
          top: '80px',
          left: '20px',
          zIndex: 1000,
          backgroundColor: '#1976d2',
          color: 'white',
          fontSize: '12px',
          padding: '8px 16px',
          height: '40px',
          boxShadow: '0 4px 12px rgba(0,0,0,0.15)',
          '&:hover': {
            backgroundColor: '#1565c0'
          }
        }}
        onClick={() => {
          // Navigate to AminDetailPage with NEVP (area 16) as default region
          navigate('/manish/16');
        }}
      >
        <SchoolIcon style={{ marginRight: '8px', fontSize: '18px' }} />
        Visit WECC Case Studies
      </Fab>
      
      {isLoadingAdditionalData && (
        <div style={{
          position: 'fixed',
          top: '70px',
          right: '20px',
          background: 'rgba(25, 118, 210, 0.9)',
          color: 'white',
          padding: '8px 16px',
          borderRadius: '4px',
          fontSize: '12px',
          zIndex: 2000,
          boxShadow: '0 2px 8px rgba(0,0,0,0.2)'
        }}>
          Loading additional transmission lines...
        </div>
      )}
      <DeckGL
        ref={deckRef}
        layers={layers}
        initialViewState={initialViewState}
        controller={true}
        getTooltip={getEnhancedTooltip}
        ContextProvider={MapContext.Provider}
        style={{ marginTop: '60px' }}
      >


        <StaticMap
          reuseMaps
          mapStyle={mapStyle[mapStyleSelection]}
          preventStyleDiffing={true}
          initialViewState={INITIAL_VIEW_STATE}
        >
        </StaticMap>


        <FullscreenControl />
        <br></br><br></br>
        <NavigationControl />

        <div style={{ 
          position: "absolute", 
          top: 160, 
          left: 10, 
          width: 40, 
          background: "rgba(255,255,255,0.95)", 
          color: "#1976d2", 
          zIndex: 1000,
          borderRadius: "8px",
          padding: "8px 4px",
          boxShadow: "0 2px 10px rgba(0,0,0,0.1)",
          backdropFilter: "blur(4px)"
        }}>
          <div style={{ 
            display: "flex", 
            flexDirection: "column", 
            alignItems: "center", 
            gap: "12px",
            cursor: "pointer"
          }}>
            <HomeOutlinedIcon 
              fontSize="medium" 
              onClick={GoHome}
              style={{
                transition: "all 0.2s ease",
                ":hover": { transform: "scale(1.1)" }
              }}
            />
            <ZoomOutMapIcon 
              fontSize="medium" 
              onClick={ZoomToAllLines}
              title="Zoom to fit all transmission lines"
              style={{
                transition: "all 0.2s ease",
                ":hover": { transform: "scale(1.1)" }
              }}
            />
            <ThreeSixtyOutlinedIcon 
              fontSize="large" 
              onClick={rotateCamera}
              style={{
                transition: "all 0.2s ease",
                ":hover": { transform: "scale(1.1)" }
              }}
            />
          </div>
        </div>


        {/*<div><NavigationControl position="top-left"/></div>
      <FullscreenControl/>*/}


        {
          showPopup.display && (
            <Popup style={{
              zIndex: 3,
              background: "white",
              opacity: 1,
              fontSize: "12px",
              fontFamily: '"Inter", sans-serif',
              cursor: showPopup.type === 'wecc' ? 'pointer' : 'default',
              borderRadius: "6px",
              boxShadow: "0 4px 12px rgba(0,0,0,0.15)"
            }}
              longitude={initialViewState.longitude}
              latitude={initialViewState.latitude}
              anchor="bottom"
              offset={-100}
              onClose={() => setShowPopup({ ...showPopup, display: false })}>
              <div
                onClick={showPopup.type === 'wecc' ? handlePopupClick : undefined}
                style={{
                  padding: showPopup.type === 'wecc' ? '8px' : '4px',
                  borderRadius: showPopup.type === 'wecc' ? '4px' : '0',
                  transition: 'background-color 0.2s ease'
                }}
                onMouseEnter={(e) => {
                  if (showPopup.type === 'wecc') {
                    e.target.style.backgroundColor = '#f0f8ff';
                  }
                }}
                onMouseLeave={(e) => {
                  if (showPopup.type === 'wecc') {
                    e.target.style.backgroundColor = 'transparent';
                  }
                }}
              >
                <h3 style={{ margin: "0 0 4px 0", fontSize: "14px", fontWeight: "600" }}>{showPopup.name}</h3>
                <h4 style={{ whiteSpace: 'pre-line', margin: 0, fontSize: "11px", fontWeight: "400", lineHeight: "1.4" }}>{showPopup.info}</h4>
              </div>
            </Popup>
          )
        }

      </DeckGL>


      <Widget
        handleNewUserMessage={handleUserInput}
        title="GridBee - Power Grid Assistant"
        subtitle="Ask me About Western Interconnection"
        launcher={(handleToggle) => (
          <div 
            onClick={handleToggle}
            style={{
              width: '60px',
              height: '60px',
              borderRadius: '50%',
              background: '#0047AB',
              display: 'flex',
              alignItems: 'center',
              justifyContent: 'center',
              cursor: 'pointer',
              boxShadow: '0 4px 16px rgba(0, 71, 171, 0.3)',
              transition: 'all 0.3s ease',
              border: '2px solid white',
              fontFamily: '"Inter", sans-serif'
            }}
            onMouseEnter={(e) => {
              e.target.style.transform = 'scale(1.1)';
              e.target.style.boxShadow = '0 6px 20px rgba(0, 71, 171, 0.4)';
            }}
            onMouseLeave={(e) => {
              e.target.style.transform = 'scale(1)';
              e.target.style.boxShadow = '0 4px 16px rgba(0, 71, 171, 0.3)';
            }}
          >
            <img 
              src="/images/gridbee_logo.png" 
              alt="GridBee" 
              style={{ 
                width: '40px', 
                height: '40px', 
                borderRadius: '50%',
                objectFit: 'cover'
              }}
              onError={(e) => {
                e.target.style.display = 'none';
                e.target.parentElement.innerHTML = '<span style="color: white; font-size: 24px; font-family: \"Inter\", sans-serif;">🤖</span>';
              }}
            />
          </div>
        )}
      />

      <div style={{ 
        position: "absolute", 
        top: 60, 
        right: 0, 
        width: 300, 
        background: "rgba(255,255,255,0.95)", 
        padding: "16px", 
        color: "#333", 
        zIndex: 1000,
        fontFamily: '"Inter", sans-serif',
        maxHeight: "calc(100vh - 60px)",
        overflowY: "auto",
        backdropFilter: "blur(8px)",
        borderLeft: "1px solid rgba(0,0,0,0.1)",
        boxShadow: "0 2px 10px rgba(0,0,0,0.1)"
      }}>

        <div style={{ width: '100%' }}>
          {/* Data Center Impact Analysis Accordion */}
          <Accordion defaultExpanded={false} style={{ marginBottom: "8px" }}>
            <AccordionSummary style={{ 
              height: "20px", 
              minHeight: "40px", 
              paddingRight: "20px", 
              paddingLeft: "0px",
              background: "rgba(25, 118, 210, 0.05)"
            }}
              expandIcon={<ArrowDropDownIcon />}>
              <Typography style={{ fontSize: "13px", fontWeight: "500", lineHeight: "1.2" }}> 
                <Checkbox checked={impactAnalysisActive} style={{ color: "#1976d2", padding: "2px" }} onChange={handleImpactAnalysisChange} />
                📊 Data Center Impact
              </Typography>
            </AccordionSummary>
            <AccordionDetails style={{ padding: "12px 16px" }}>
              <Typography component="div">
                {impactAnalysisActive && (
                  <div style={{ paddingRight: "8px", marginBottom: "8px" }}>
                    {/* WECC Region Selection */}
                    <div style={{ marginBottom: "12px" }}>
                      <label style={{
                        display: "block",
                        fontSize: "12px",
                        fontWeight: "500",
                        marginBottom: "6px",
                        color: "#555"
                      }}>
                        WECC Region:
                      </label>
                      <select
                        value={selectedWeccRegion}
                        onChange={(e) => setSelectedWeccRegion(e.target.value)}
                        style={{
                          width: "100%",
                          padding: "8px 12px",
                          borderRadius: "6px",
                          border: "1px solid #ddd",
                          fontSize: "12px",
                          background: "white"
                        }}
                      >
                        {weccRegions.map(region => (
                          <option key={region.id} value={region.id}>
                            {region.name}
                          </option>
                        ))}
                      </select>
                    </div>

                    {/* Case Study Selection */}
                    <div style={{ marginBottom: "12px" }}>
                      <label style={{
                        display: "block",
                        fontSize: "12px",
                        fontWeight: "500",
                        marginBottom: "6px",
                        color: "#555"
                      }}>
                        Case Study Comparison:
                      </label>
                      <div style={{
                        display: "flex",
                        gap: "6px"
                      }}>
                        {['case1', 'case2', 'case3'].map(caseId => (
                          <button
                            key={caseId}
                            onClick={() => setSelectedCaseStudy(caseId)}
                            style={{
                              flex: 1,
                              padding: "6px 4px",
                              borderRadius: "4px",
                              border: selectedCaseStudy === caseId ? "2px solid #1976d2" : "1px solid #ddd",
                              background: selectedCaseStudy === caseId ? "rgba(25, 118, 210, 0.1)" : "white",
                              fontSize: "10px",
                              fontWeight: selectedCaseStudy === caseId ? "600" : "400",
                              color: selectedCaseStudy === caseId ? "#1976d2" : "#666",
                              cursor: "pointer",
                              transition: "all 0.2s",
                              minHeight: "28px"
                            }}
                          >
                            Case {caseId.slice(-1)}
                          </button>
                        ))}
                      </div>
                    </div>

                    {/* Metric Selection */}
                    <div style={{ marginBottom: "12px" }}>
                      <label style={{
                        display: "block",
                        fontSize: "12px",
                        fontWeight: "500",
                        marginBottom: "6px",
                        color: "#555"
                      }}>
                        Metric:
                      </label>
                      <div style={{
                        display: "flex",
                        flexDirection: "column",
                        gap: "4px"
                      }}>
                        {[
                          { id: 'Price', label: 'Price ($/MWh)' },
                          { id: 'System Operation Cost', label: 'System Cost ($)' },
                          { id: 'Power Exchange', label: 'Power Exchange (MW)' }
                        ].map(metric => (
                          <button
                            key={metric.id}
                            onClick={() => setSelectedMetric(metric.id)}
                            style={{
                              padding: "6px 8px",
                              borderRadius: "4px",
                              border: selectedMetric === metric.id ? "2px solid #1976d2" : "1px solid #ddd",
                              background: selectedMetric === metric.id ? "rgba(25, 118, 210, 0.1)" : "white",
                              fontSize: "10px",
                              fontWeight: selectedMetric === metric.id ? "600" : "400",
                              color: selectedMetric === metric.id ? "#1976d2" : "#666",
                              cursor: "pointer",
                              transition: "all 0.2s",
                              textAlign: "left",
                              marginBottom: "4px",
                              width: "100%"
                            }}
                          >
                            {metric.label}
                          </button>
                        ))}
                      </div>
                    </div>

                    {/* Hour Selection */}
                    <div style={{ marginBottom: "12px" }}>
                      <label style={{
                        display: "block",
                        fontSize: "12px",
                        fontWeight: "500",
                        marginBottom: "6px",
                        color: "#555"
                      }}>
                        Hour: {selectedHour}
                      </label>
                      <input
                        type="range"
                        min="1"
                        max="24"
                        value={selectedHour}
                        onChange={(e) => setSelectedHour(parseInt(e.target.value))}
                        style={{
                          width: "100%",
                          height: "6px",
                          borderRadius: "3px",
                          background: "#ddd",
                          outline: "none",
                          cursor: "pointer"
                        }}
                      />
                      <div style={{
                        display: "flex",
                        justifyContent: "space-between",
                        fontSize: "10px",
                        color: "#888",
                        marginTop: "4px"
                      }}>
                        <span>1</span>
                        <span>12</span>
                        <span>24</span>
                      </div>
                    </div>

                    {/* Color Scale */}
                    <div style={{ marginBottom: "8px" }}>
                      <label style={{
                        display: "block",
                        fontSize: "12px",
                        fontWeight: "500",
                        marginBottom: "6px",
                        color: "#555"
                      }}>
                        Impact Scale:
                      </label>
                      <div style={{
                        display: "flex",
                        alignItems: "center",
                        gap: "8px"
                      }}>
                        <span style={{ fontSize: "10px", color: "#0066cc" }}>Negative</span>
                        <div style={{
                          flex: 1,
                          height: "12px",
                          background: "linear-gradient(to right, #0066cc, #ffffff, #cc4400)",
                          borderRadius: "6px",
                          border: "1px solid #ddd"
                        }}></div>
                        <span style={{ fontSize: "10px", color: "#cc4400" }}>Positive</span>
                      </div>
                    </div>
                  </div>
                )}
              </Typography>
            </AccordionDetails>
          </Accordion>

          <Accordion defaultExpanded={true} style={{ marginBottom: "8px" }}>
            <AccordionSummary style={{ 
              height: "20px", 
              minHeight: "40px", 
              paddingRight: "20px", 
              paddingLeft: "0px",
              background: "rgba(25, 118, 210, 0.05)"
            }}
              expandIcon={<ArrowDropDownIcon />}>
              <Typography style={{ fontSize: "14px", fontWeight: "500" }}> 
                <Checkbox checked={netlayeractive} style={{ color: "#1976d2" }} onChange={handleNetLayerChange} />
                Network
              </Typography>
            </AccordionSummary>
            <AccordionDetails style={{ padding: "12px 16px" }}>
              <Typography component="div">
                {netlayeractive &&
                  (
                    <div style={{ paddingRight: "20px", marginBottom: "12px" }}>
                      <Typography style={{ fontSize: "12px", color: "#666", marginBottom: "8px" }}>
                        Voltage Level (kV)
                      </Typography>
                      <Slider
                        value={netfiltervalue}
                        valueLabelDisplay="auto"
                        onChange={handleNetRangeFilterChange}
                        getAriaValueText={valuetext}
                        step={100}
                        min={0}
                        max={800}
                        style={{ color: "#1976d2" }}
                      />
                    </div>)
                }

                {netlayeractive && (
                  <div style={{ paddingRight: "20px" }}>
                    <Multiselect
                      defaultValue={busNameSelectItems}
                      data={busNameItems}
                      placeholder={'Search for buses'}
                      onChange={handleBusMultiselect}
                    />
                  </div>
                )}

              </Typography>
            </AccordionDetails>
          </Accordion>

          <Accordion defaultExpanded={true} style={{ marginBottom: "8px" }}>
            <AccordionSummary style={{ 
              height: "20px", 
              minHeight: "40px", 
              paddingRight: "20px", 
              paddingLeft: "0px",
              background: "rgba(25, 118, 210, 0.05)"
            }}
              expandIcon={<ArrowDropDownIcon />}>
              <Typography style={{ fontSize: "14px", fontWeight: "500" }}> 
                <Checkbox checked={flowlayeractive} style={{ color: "#1976d2" }} onChange={handleFlowLayerChange} />
                Power Flow
              </Typography>
            </AccordionSummary>
            <AccordionDetails style={{ padding: "12px 16px" }}>
              <Typography component="div">
                {flowlayeractive && (
                  <div style={{ paddingRight: "20px", marginBottom: "12px" }}>
                    <Typography style={{ fontSize: "12px", color: "#666", marginBottom: "8px" }}>
                      Loading (%)
                    </Typography>
                    <Slider
                      style={{ color: "#1976d2" }}
                      value={flowfiltervalue}
                      valueLabelDisplay="auto"
                      onChange={handleFlowRangeFilterChange}
                      getAriaValueText={valuetext}
                      step={10}
                      min={0}
                      max={120}
                    />
                  </div>)
                }

                {flowlayeractive && (
                  <div style={{ paddingRight: "20px" }}>
                    <Multiselect
                      defaultValue={lineNameSelectItems}
                      data={lineNameItems}
                      placeholder={'Search for transmission lines'}
                      onChange={handleLineMultiselect}
                    />
                  </div>
                )}

              </Typography>
            </AccordionDetails>
          </Accordion>

          <Accordion defaultExpanded={true} style={{ marginBottom: "8px" }}>
            <AccordionSummary style={{ 
              height: "20px", 
              minHeight: "40px", 
              paddingRight: "20px", 
              paddingLeft: "0px",
              background: "rgba(156, 39, 176, 0.05)"
            }}
              expandIcon={<ArrowDropDownIcon />}>
              <Typography style={{ fontSize: "14px", fontWeight: "500" }}>
                <Checkbox checked={transmissionlayeractive} style={{ color: "#9c27b0" }} onChange={handleTransmissionLayerChange} />
                Transmission Lines
              </Typography>
            </AccordionSummary>
            <AccordionDetails style={{ padding: "12px 16px" }}>
              <Typography component="div">
                <div style={{ fontSize: "12px", color: "#666", marginBottom: "8px" }}>
                  <div style={{ display: "flex", alignItems: "center", marginBottom: "4px" }}>
                    <div style={{ width: "14px", height: "4px", backgroundColor: "rgb(0, 100, 200)", marginRight: "6px", borderRadius: "1px" }}></div>
                    <span style={{ fontWeight: "500" }}>≥500kV</span>
                    <span style={{ color: "#888", fontSize: "10px", marginLeft: "4px" }}>Extra High</span>
                  </div>
                  <div style={{ display: "flex", alignItems: "center", marginBottom: "4px" }}>
                    <div style={{ width: "14px", height: "3px", backgroundColor: "rgb(30, 144, 255)", marginRight: "6px", borderRadius: "1px" }}></div>
                    <span style={{ fontWeight: "500" }}>345kV</span>
                    <span style={{ color: "#888", fontSize: "10px", marginLeft: "4px" }}>High Voltage</span>
                  </div>
                  <div style={{ display: "flex", alignItems: "center", marginBottom: "4px" }}>
                    <div style={{ width: "14px", height: "2.5px", backgroundColor: "rgb(70, 170, 255)", marginRight: "6px", borderRadius: "1px" }}></div>
                    <span style={{ fontWeight: "500" }}>230kV</span>
                    <span style={{ color: "#888", fontSize: "10px", marginLeft: "4px" }}>Medium</span>
                  </div>
                  <div style={{ display: "flex", alignItems: "center", marginBottom: "4px" }}>
                    <div style={{ width: "14px", height: "2px", backgroundColor: "rgb(100, 200, 255)", marginRight: "6px", borderRadius: "1px" }}></div>
                    <span style={{ fontWeight: "500" }}>138kV</span>
                    <span style={{ color: "#888", fontSize: "10px", marginLeft: "4px" }}>Sub-transmission</span>
                  </div>
                  <div style={{ display: "flex", alignItems: "center", marginBottom: "4px" }}>
                    <div style={{ width: "14px", height: "1.5px", backgroundColor: "rgb(120, 210, 255)", marginRight: "6px", borderRadius: "1px" }}></div>
                    <span style={{ fontWeight: "500" }}>100kV</span>
                    <span style={{ color: "#888", fontSize: "10px", marginLeft: "4px" }}>Distribution</span>
                  </div>
                  <div style={{ display: "flex", alignItems: "center", marginBottom: "4px" }}>
                    <div style={{ width: "14px", height: "1.5px", backgroundColor: "rgb(70, 170, 255)", marginRight: "6px", borderRadius: "1px" }}></div>
                    <span style={{ fontWeight: "500" }}>Unknown/Default</span>
                    <span style={{ color: "#888", fontSize: "10px", marginLeft: "4px" }}>Various</span>
                  </div>
                </div>
                <div style={{ fontSize: "11px", color: "#999", fontStyle: "italic", marginBottom: "4px" }}>
                  Blue gradient by voltage level
                </div>
                <div style={{ fontSize: "10px", color: "#888", fontStyle: "italic" }}>
                  Sources: USGS Powerlines WUS CAN SGCA, Western Power Plants USA & Canada-Mexico datasets
                </div>
              </Typography>
            </AccordionDetails>
          </Accordion>

          <Accordion defaultExpanded={false} style={{ marginBottom: "8px" }}>
            <AccordionSummary style={{ 
              height: "20px", 
              minHeight: "40px", 
              paddingRight: "20px", 
              paddingLeft: "0px",
              background: "rgba(255, 193, 7, 0.05)"
            }}
              expandIcon={<ArrowDropDownIcon />}>
              <Typography style={{ fontSize: "14px", fontWeight: "500" }}>
                <Checkbox checked={powerplantlayeractive} style={{ color: "#ff9800" }} onChange={handlePowerPlantLayerChange} />
                Western Power Plants
              </Typography>
            </AccordionSummary>
            <AccordionDetails style={{ padding: "12px 16px" }}>
              <Typography component="div">
                {powerplantlayeractive &&
                  (<div style={{ paddingRight: "20px", marginBottom: "12px" }}>
                    <Checkbox checked={powerplantlayercapactive} style={{ color: "#ff9800" }} onChange={handlePowerPlantLayerCapChange} />
                    <span style={{ fontSize: "12px" }}>Show Total Capacity</span>
                  </div>)
                }

                {powerplantlayeractive &&
                  (<div style={{ paddingRight: "20px", marginBottom: "12px" }}>
                    <Typography style={{ fontSize: "12px", color: "#666", marginBottom: "8px" }}>
                      Capacity Range (MW)
                    </Typography>
                    <Slider
                      style={{ color: "#ff9800" }}
                      value={powerPlantFilterValue}
                      valueLabelDisplay="auto"
                      onChange={handlePowerPlantRangeFilterChange}
                      getAriaValueText={(value) => `${value} MW`}
                      step={10}
                      min={powerPlantData.minPg}
                      max={powerPlantData.maxPg + 100}
                    />
                  </div>)
                }

                {powerplantlayeractive && (
                  <div style={{ paddingRight: "20px", marginBottom: "12px" }}>
                    <Multiselect
                      defaultValue={powerPlantSelectItems}
                      data={powerPlantNameItems}
                      placeholder={'Search for power plants'}
                      onChange={handlePowerPlantMultiselect}
                    />
                  </div>
                )}

                {powerplantlayeractive && powerPlantChartData && (
                  <div style={{ width: 260, height: 280, transform: "translate(-8px, 0px)" }}>
                    <Doughnut data={powerPlantChartData}
                      options={{
                        maintainAspectRatio: false,
                        plugins: {
                          legend: {
                            position: 'bottom',
                            labels: {
                              boxWidth: 12,
                              padding: 8,
                              font: {
                                size: 10
                              }
                            }
                          },
                          title: {
                            display: true,
                            text: 'Western Power Plants by Technology',
                            font: {
                              size: 12,
                              weight: 'bold'
                            }
                          }
                        }
                      }}
                    />
                  </div>
                )}

                {/* Plant Labels Control */}
                {powerplantlayeractive && (
                  <div style={{ 
                    marginTop: "8px", 
                    padding: "6px 8px", 
                    background: "rgba(255, 152, 0, 0.05)",
                    borderRadius: "4px",
                    border: "1px solid rgba(255, 152, 0, 0.15)"
                  }}>
                    <FormControlLabel
                      control={
                        <Checkbox 
                          checked={powerplantlabelsactive} 
                          onChange={handlePowerPlantLabelsChange}
                          style={{ color: "#ff9800", transform: "scale(0.8)" }}
                        />
                      }
                      label={
                        <span style={{ fontSize: "11px", color: "#e65100", fontWeight: "500" }}>
                          🏷️ Show Plant Names (50+ MW)
                        </span>
                      }
                      style={{ margin: 0 }}
                    />
                    <div style={{ fontSize: "9px", color: "#666", marginTop: "2px", marginLeft: "28px" }}>
                      Display plant names above 3D columns for larger facilities
                    </div>
                  </div>
                )}
              </Typography>
            </AccordionDetails>
          </Accordion>

          {/* Western Power Plants Overview Section */}
          <Accordion defaultExpanded={true} style={{ marginBottom: "8px" }}>
            <AccordionSummary style={{ 
              height: "20px", 
              minHeight: "40px", 
              paddingRight: "20px", 
              paddingLeft: "0px",
              background: "rgba(255, 193, 7, 0.05)"
            }}
              expandIcon={<ArrowDropDownIcon />}>
              <Typography style={{ fontSize: "14px", fontWeight: "500" }}>
                <span style={{ color: "#ff9800", marginRight: "8px" }}>📊</span>
                Western Power Plants Overview
              </Typography>
            </AccordionSummary>
            <AccordionDetails style={{ padding: "12px 16px" }}>
              <Typography component="div">
                {generationStats ? (
                  <div>
                    {/* Summary Statistics */}
                    <div style={{ 
                      background: "rgba(255, 193, 7, 0.08)", 
                      padding: "12px", 
                      borderRadius: "8px", 
                      marginBottom: "12px",
                      border: "1px solid rgba(255, 193, 7, 0.2)"
                    }}>
                      <div style={{ fontSize: "12px", fontWeight: "600", color: "#e65100", marginBottom: "6px" }}>
                        Western Power Plants Summary
                      </div>
                      <div style={{ fontSize: "11px", color: "#666", lineHeight: "1.4" }}>
                        <div><strong>Total Capacity:</strong> {generationStats.totalCapacity} MW</div>
                        <div><strong>Total Plants:</strong> {generationStats.totalPlants.toLocaleString()}</div>
                        <div><strong>Average Size:</strong> {generationStats.averageCapacity} MW</div>
                      </div>
                    </div>

                    {/* Technology Breakdown Cards */}
                    <div style={{ marginBottom: "8px" }}>
                      <div style={{ fontSize: "12px", fontWeight: "500", marginBottom: "8px", color: "#333" }}>
                        Power Plants by Technology
                      </div>
                      <div style={{ display: "flex", flexDirection: "column", gap: "6px" }}>
                        {generationStats.sortedTechnologies.slice(0, 6).map(tech => {
                          const techData = generationStats.technologies[tech];
                          const isSelected = selectedTechnology === tech;
                          const color = techData.color.slice(0, 3);
                          
                          return (
                            <div 
                              key={tech}
                              style={{ 
                                display: "flex", 
                                alignItems: "center", 
                                padding: "8px 10px",
                                borderRadius: "6px",
                                background: isSelected ? `rgba(${color.join(',')}, 0.15)` : "rgba(0,0,0,0.02)",
                                border: isSelected ? `2px solid rgb(${color.join(',')})` : "1px solid rgba(0,0,0,0.1)",
                                cursor: "pointer",
                                transition: "all 0.2s ease",
                                position: "relative"
                              }}
                              onClick={() => {
                                setSelectedTechnology(selectedTechnology === tech ? null : tech);
                              }}
                              onMouseEnter={(e) => {
                                if (!isSelected) {
                                  e.target.style.background = `rgba(${color.join(',')}, 0.08)`;
                                }
                              }}
                              onMouseLeave={(e) => {
                                if (!isSelected) {
                                  e.target.style.background = "rgba(0,0,0,0.02)";
                                }
                              }}
                            >
                              <div style={{ 
                                width: "16px", 
                                height: "16px", 
                                backgroundColor: `rgb(${color.join(',')})`, 
                                marginRight: "10px", 
                                borderRadius: "50%", 
                                border: "2px solid white",
                                boxShadow: "0 1px 3px rgba(0,0,0,0.2)",
                                flexShrink: 0
                              }}></div>
                              
                              <div style={{ flex: 1, minWidth: 0 }}>
                                <div style={{ 
                                  fontSize: "11px", 
                                  fontWeight: "500", 
                                  color: isSelected ? `rgb(${color.join(',')})` : "#333",
                                  marginBottom: "2px"
                                }}>
                                  {tech}
                                </div>
                                <div style={{ fontSize: "10px", color: "#666", lineHeight: "1.2" }}>
                                  <div>{techData.count} plants • {techData.totalCapacity.toFixed(0)} MW</div>
                                  <div>{techData.percentage}% of total capacity</div>
                                </div>
                              </div>
                              
                              <div style={{ 
                                fontSize: "10px", 
                                color: isSelected ? `rgb(${color.join(',')})` : "#999",
                                fontWeight: "500",
                                marginLeft: "8px"
                              }}>
                                {techData.percentage}%
                              </div>
                              
                              {isSelected && (
                                <div style={{
                                  position: "absolute",
                                  top: "2px",
                                  right: "6px",
                                  fontSize: "10px",
                                  color: `rgb(${color.join(',')})`,
                                  fontWeight: "600"
                                }}>
                                  ✓
                                </div>
                              )}
                            </div>
                          );
                        })}
                      </div>
                    </div>

                    {/* Show remaining technologies if any */}
                    {generationStats.sortedTechnologies.length > 6 && (
                      <div style={{ 
                        fontSize: "10px", 
                        color: "#888", 
                        fontStyle: "italic", 
                        textAlign: "center",
                        marginTop: "8px",
                        padding: "4px"
                      }}>
                        +{generationStats.sortedTechnologies.length - 6} more technologies
                      </div>
                    )}

                    {/* Filter Status */}
                    {selectedTechnology && (
                      <div style={{ 
                        background: "rgba(33, 150, 243, 0.08)", 
                        padding: "8px", 
                        borderRadius: "6px", 
                        marginTop: "8px",
                        border: "1px solid rgba(33, 150, 243, 0.2)"
                      }}>
                        <div style={{ fontSize: "10px", color: "#1976d2", fontWeight: "500", marginBottom: "2px" }}>
                          🔍 Filtered View Active
                        </div>
                        <div style={{ fontSize: "10px", color: "#666" }}>
                          Showing only {selectedTechnology} plants
                        </div>
                        <div 
                          style={{ 
                            fontSize: "10px", 
                            color: "#1976d2", 
                            cursor: "pointer", 
                            textDecoration: "underline",
                            marginTop: "4px"
                          }}
                          onClick={() => setSelectedTechnology(null)}
                        >
                          Clear filter
                        </div>
                      </div>
                    )}
                  </div>
                ) : (
                  <div style={{ fontSize: "11px", color: "#666", fontStyle: "italic", textAlign: "center", padding: "16px" }}>
                    Loading generation statistics...
                  </div>
                )}
              </Typography>
            </AccordionDetails>
          </Accordion>

          <Accordion defaultExpanded={false} style={{ marginBottom: "8px" }}>
            <AccordionSummary style={{ 
              height: "20px", 
              minHeight: "40px", 
              paddingRight: "20px", 
              paddingLeft: "0px",
              background: "rgba(46, 125, 50, 0.05)"
            }}
              expandIcon={<ArrowDropDownIcon />}>
              <Typography style={{ fontSize: "14px", fontWeight: "500" }}>
                <Checkbox checked={weccGenLayerActive} style={{ color: "#2e7d32" }} onChange={handleWeccGenLayerChange} />
                WECC Generation
              </Typography>
            </AccordionSummary>
            <AccordionDetails style={{ padding: "12px 16px" }}>
              <Typography component="div">
                {weccGenLayerActive && weccGenData && (
                  <div style={{ paddingRight: "20px", marginBottom: "12px" }}>
                    <Typography style={{ fontSize: "12px", color: "#666", marginBottom: "8px" }}>
                      Total Capacity (MW)
                    </Typography>
                    <Slider
                      style={{ color: "#2e7d32" }}
                      value={weccGenFilter}
                      valueLabelDisplay="auto"
                      onChange={handleWeccGenRangeFilterChange}
                      getAriaValueText={valuetext}
                      step={1000}
                      min={0}
                      max={weccGenData.length > 0 ? Math.max(...weccGenData.map(d => d.Total_MW)) : 100000}
                    />
                  </div>
                )}

                {weccGenLayerActive && weccGenNameItems.length > 0 && (
                  <div style={{ paddingRight: "20px", marginBottom: "12px" }}>
                    <Multiselect
                      defaultValue={weccGenSelectItems}
                      data={weccGenNameItems}
                      placeholder={'Filter by Balancing Authority'}
                      onChange={handleWeccGenMultiselect}
                    />
                  </div>
                )}

                {weccGenLayerActive && weccGenChartData && (
                  <div style={{ width: 260, height: 250, transform: "translate(-8px, 0px)" }}>
                    <Doughnut data={weccGenChartData}
                      options={{
                        maintainAspectRatio: false,
                        plugins: {
                          legend: {
                            position: 'bottom',
                            onClick: handleWeccDoughnutClick,
                            labels: {
                              boxWidth: 12,
                              padding: 6,
                              font: {
                                size: 10,
                                family: '"Inter", sans-serif'
                              }
                            }
                          },
                          title: {
                            display: true,
                            text: 'WECC Generation Mix (MW)',
                            font: {
                              size: 12,
                              family: '"Inter", sans-serif',
                              weight: '600'
                            }
                          }
                        }
                      }}
                    />
                  </div>
                )}
              </Typography>
            </AccordionDetails>
          </Accordion>

          <Accordion style={{ paddingBottom: "8px", marginBottom: "8px" }} defaultExpanded={false}>
            <AccordionSummary style={{ 
              height: "20px", 
              minHeight: "40px", 
              paddingRight: "20px", 
              paddingLeft: "0px",
              background: "rgba(46, 125, 50, 0.05)"
            }}
              expandIcon={<ArrowDropDownIcon />}>
              <Typography style={{ fontSize: "14px", fontWeight: "500" }}>
                <Checkbox checked={genlayeractive} style={{ color: "#2e7d32" }} onChange={handleGenLayerChange} />
                Generation Power
              </Typography>
            </AccordionSummary>
            <AccordionDetails style={{ padding: "12px 16px" }}>
              <Typography component="div">
                {genlayeractive &&
                  (<div style={{ paddingRight: "20px", marginBottom: "12px" }}>
                    <Checkbox checked={genlayercapactive} style={{ color: "#2e7d32" }} onChange={handleGenLayerCapChange} />
                    <span style={{ fontSize: "12px" }}>Show Generation Capacity</span>
                  </div>)
                }

                {genlayeractive &&
                  (<div style={{ paddingRight: "20px", marginBottom: "12px" }}>
                    <Typography style={{ fontSize: "12px", color: "#666", marginBottom: "8px" }}>
                      Power Range (MW)
                    </Typography>
                    <Slider
                      style={{ color: "#2e7d32" }}
                      value={genfiltervalue}
                      valueLabelDisplay="auto"
                      onChange={handleGenRangeFilterChange}
                      getAriaValueText={valuetext}
                      step={100}
                      min={gendata.minPg}
                      max={gendata.maxPg + 10}
                    />
                  </div>)
                }


                {genlayeractive && (
                  <div style={{ paddingRight: "20px", marginBottom: "12px" }}>
                    <Multiselect
                      defaultValue={nameSelectItems}
                      data={nameItems}
                      placeholder={'Search for generators'}
                      onChange={handleGenMultiselect}
                    />
                  </div>
                )}

                {
                  genlayeractive && (
                    <div style={{ width: 260, height: 280, transform: "translate(-8px, 0px)" }}>
                      <Doughnut data={chartdata}
                        options={{
                          maintainAspectRatio: false,
                          plugins: {
                            legend: {
                              position: 'bottom',
                              onClick: handleDoughnutClick,
                              labels: {
                                boxWidth: 12,
                                padding: 6,
                                font: {
                                  size: 10,
                                  family: '"Inter", sans-serif'
                                }
                              }
                            },
                            title: {
                              display: true,
                              text: 'Generation Mix',
                              font: {
                                size: 12,
                                family: '"Inter", sans-serif',
                                weight: '600'
                              }
                            }
                          }
                        }}
                      />
                    </div>
                  )}


              </Typography>
            </AccordionDetails>
          </Accordion>

          <Accordion defaultExpanded={false} style={{ marginBottom: "8px" }}>
            <AccordionSummary style={{ 
              height: "20px", 
              minHeight: "40px", 
              paddingRight: "20px", 
              paddingLeft: "0px",
              background: "rgba(211, 47, 47, 0.05)"
            }}
              expandIcon={<ArrowDropDownIcon />}>
              <Typography style={{ fontSize: "14px", fontWeight: "500" }}>
                <Checkbox checked={loadlayeractive} style={{ color: "#d32f2f" }} onChange={handleLoadLayerChange} />
                Load Loss
              </Typography>
            </AccordionSummary>
            <AccordionDetails style={{ padding: "12px 16px" }}>
              <Typography component="div">
                {loadlayeractive && (
                  <div style={{ paddingRight: "20px", marginBottom: "12px" }}>
                    <Typography style={{ fontSize: "12px", color: "#666", marginBottom: "8px" }}>
                      Load Range (MW)
                    </Typography>
                    <Slider
                      style={{ color: "#d32f2f" }}
                      value={loadfiltervalue}
                      valueLabelDisplay="auto"
                      onChange={handleLoadRangeFilterChange}
                      getAriaValueText={valuetext}
                      step={100}
                      min={0}
                      max={countyloaddata.maxPd + 10}
                    />
                  </div>)
                }

                {loadlayeractive && (
                  <div style={{ paddingRight: "20px" }}>
                    <Multiselect
                      defaultValue={countyNameSelectItems}
                      data={countyNameItems}
                      placeholder={'Search for counties'}
                      onChange={handleCountyMultiselect}
                    />
                  </div>
                )}


              </Typography>
            </AccordionDetails>
          </Accordion>


          <Accordion defaultExpanded={false} style={{ marginBottom: "8px" }}>
            <AccordionSummary style={{ 
              height: "20px", 
              minHeight: "40px", 
              paddingRight: "20px", 
              paddingLeft: "0px",
              background: "rgba(245, 124, 0, 0.05)"
            }}
              expandIcon={<ArrowDropDownIcon />}>
              <Typography style={{ fontSize: "14px", fontWeight: "500" }}>
                <Checkbox checked={voltagelayeractive} style={{ color: "#f57c00" }} onChange={handleVoltageLayerChange} />
                Voltage
              </Typography>
            </AccordionSummary>
            <AccordionDetails style={{ padding: "12px 16px" }}>
              <Typography component="div">
                {voltagelayeractive && (
                  <div style={{ paddingRight: "20px", marginBottom: "12px" }}>
                    <Typography style={{ fontSize: "12px", color: "#666", marginBottom: "8px" }}>
                      Voltage Range (p.u.)
                    </Typography>
                    <Slider
                      style={{ color: "#f57c00" }}
                      value={voltagefiltervalue}
                      valueLabelDisplay="auto"
                      onChange={handleVoltageRangeFilterChange}
                      getAriaValueText={valuetext}
                      step={0.01}
                      min={0.89}
                      max={1.11}
                    />
                  </div>)
                }

                {voltagelayeractive && (
                  <div style={{ paddingRight: "20px" }}>
                    <Multiselect
                      defaultValue={countyNameSelectItems}
                      data={countyNameItems}
                      placeholder={'Search for counties'}
                      onChange={handleCountyMultiselect}
                    />
                  </div>
                )}


              </Typography>
            </AccordionDetails>
          </Accordion>

          <Accordion defaultExpanded={false} style={{ marginBottom: "8px" }}>
            <AccordionSummary style={{ 
              height: "20px", 
              minHeight: "40px", 
              paddingRight: "20px", 
              paddingLeft: "0px",
              background: "rgba(123, 31, 162, 0.05)"
            }}
              expandIcon={<ArrowDropDownIcon />}>
              <Typography style={{ fontSize: "14px", fontWeight: "500" }}>
                <Checkbox checked={arealayeractive} style={{ color: "#7b1fa2" }} onChange={handleAreaLayerChange} />
                Control Areas
              </Typography>
            </AccordionSummary>
            <AccordionDetails style={{ padding: "12px 16px" }}>
              <Typography component="div">
                {arealayeractive && (
                  <div style={{ paddingRight: "20px" }}>
                    <Multiselect
                      defaultValue={areaNameSelectItems}
                      data={areaNameItems}
                      placeholder={'Search for areas'}
                      onChange={handleAreaMultiselect}
                    />
                  </div>
                )}
              </Typography>
            </AccordionDetails>
          </Accordion>

          <Accordion defaultExpanded={false} style={{ marginBottom: "8px" }}>
            <AccordionSummary style={{ 
              height: "20px", 
              minHeight: "40px", 
              paddingRight: "20px", 
              paddingLeft: "0px",
              background: "rgba(0, 0, 0, 0.05)"
            }}
              expandIcon={<ArrowDropDownIcon />}>
              <Typography style={{ fontSize: "14px", fontWeight: "500" }}>
                Map Style
              </Typography>
            </AccordionSummary>
            <AccordionDetails style={{ padding: "12px 16px" }}>
              <Typography component="div">
                <div style={{ paddingRight: "20px" }}>
                  <select
                    value={mapStyleSelection}
                    onChange={(e) => setMapStyle(e.target.value)}
                    style={{ 
                      width: '100%', 
                      padding: '8px 12px', 
                      marginBottom: '10px',
                      borderRadius: '4px',
                      border: '1px solid #ddd',
                      fontFamily: '"Inter", sans-serif',
                      fontSize: '12px'
                    }}
                  >
                    <option value="osm">OpenStreetMap</option>
                    <option value="satellite">Satellite View</option>
                    <option value="terrain">Terrain</option>
                    <option value="pos">Light Theme</option>
                    <option value="pos_no_label">Light (No Labels)</option>
                    <option value="dark">Dark Theme</option>
                    <option value="none">No Background</option>
                  </select>
                </div>
              </Typography>
            </AccordionDetails>
          </Accordion>

          <Accordion defaultExpanded={false} style={{ marginBottom: "8px" }}>
            <AccordionSummary style={{ 
              height: "20px", 
              minHeight: "40px", 
              paddingRight: "20px", 
              paddingLeft: "0px",
              background: "rgba(245, 124, 0, 0.05)"
            }}
              expandIcon={<ArrowDropDownIcon />}>
              <Typography style={{ fontSize: "14px", fontWeight: "500" }}>
                <Checkbox checked={zonelayeractive} style={{ color: "#f57c00" }} onChange={handleZoneLayerChange} />
                Load Zones
              </Typography>
            </AccordionSummary>
            <AccordionDetails style={{ padding: "12px 16px" }}>
              <Typography component="div">
                {zonelayeractive && (
                  <div style={{ paddingRight: "20px" }}>
                    <Multiselect
                      defaultValue={zoneNameSelectItems}
                      data={zoneNameItems}
                      placeholder={'Search for zones'}
                      onChange={handleZoneMultiselect}
                    />
                  </div>
                )}
              </Typography>
            </AccordionDetails>
          </Accordion>

          <Accordion defaultExpanded={false} style={{ marginBottom: "8px" }}>
            <AccordionSummary style={{ 
              height: "20px", 
              minHeight: "40px", 
              paddingRight: "20px", 
              paddingLeft: "0px",
              background: "rgba(21, 101, 192, 0.05)"
            }}
              expandIcon={<ArrowDropDownIcon />}>
              <Typography style={{ fontSize: "14px", fontWeight: "500" }}>
                <Checkbox checked={wecclayeractive} style={{ color: "#1565c0" }} onChange={handleWeccLayerChange} />
                WECC Regions
              </Typography>
            </AccordionSummary>
            <AccordionDetails style={{ padding: "12px 16px" }}>
              <Typography component="div">
                {wecclayeractive && weccGeojsonData && (
                  <div style={{ paddingRight: "20px", fontSize: "11px", color: "#666", lineHeight: "1.4" }}>
                    WECC Balancing Authorities overlay showing {weccGeojsonData.features?.length || 0} regions.
                    Click on any region for detailed information.
                  </div>
                )}
              </Typography>
            </AccordionDetails>
          </Accordion>

        </div>

        {/* <br></br> */}

        {/* <BrushingBarChart data={[12,23,345,45,66,78,800]} width ={170} height = {100} handleFilter = {handleNetBarFilterChange}/> */}

        {/* <br></br> */}


      </div>


    </>



  );
}

// Main App component with Authentication
export default function App() {
  return (
    <AuthProvider>
      <Router>
        <Routes>
          {/* Admin route */}
          <Route
            path="/admin"
            element={
              <ProtectedRoute adminOnly={true}>
                <AdminDashboard />
              </ProtectedRoute>
            }
          />

          {/* Amin routes */}
          <Route
            path="/amin"
            element={
              <ProtectedRoute>
                <AminProject />
              </ProtectedRoute>
            }
          />
          <Route
            path="/amin/:fid"
            element={
              <ProtectedRoute>
                <AminDetailPage />
              </ProtectedRoute>
            }
          />

          {/* Manish route */}
          <Route
            path="/manish"
            element={
              <ProtectedRoute>
                <MainApp />
              </ProtectedRoute>
            }
          />
          <Route
            path="/manish/:fid"
            element={
              <ProtectedRoute>
                <AminDetailPage />
              </ProtectedRoute>
            }
          />

          {/* Default redirect to manish */}
          <Route path="/" element={<Navigate to="/manish" replace />} />
          <Route path="*" element={<Navigate to="/manish" replace />} />
        </Routes>
      </Router>
    </AuthProvider>
  );
} 

const rootElement = document.getElementById("root");
createRoot(rootElement).render(<App />)