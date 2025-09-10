import React, { useState, useEffect } from 'react';
import { useNavigate } from 'react-router-dom';
import { DeckGL } from '@deck.gl/react';
import { GeoJsonLayer } from '@deck.gl/layers';
import { MapView } from '@deck.gl/core';
import { StaticMap } from 'react-map-gl';
import ArrowBackIcon from '@mui/icons-material/ArrowBack';
import Accordion from "@mui/material/Accordion";
import AccordionSummary from "@mui/material/AccordionSummary";
import AccordionDetails from "@mui/material/AccordionDetails";
import ArrowDropDownIcon from '@mui/icons-material/ArrowDropDown';
import Typography from '@mui/material/Typography';
import Checkbox from '@mui/material/Checkbox';

// Mapbox token
const MAPBOX_ACCESS_TOKEN = 'pk.eyJ1IjoidXNtYXJ0LXdlc3RtYXAiLCJhIjoiY2tvazV6MzU2MDE4YjJ0bXd5ZDcwdm16ciJ9.q2BIGvGPAJjw1X9CdvyKSA';

// Initial view state for WECC regions
const INITIAL_VIEW_STATE = {
  longitude: -115.0,
  latitude: 44.0,
  zoom: 4.5,
  pitch: 0,
  bearing: 0
};

// WECC regions list for dropdown
const WECC_REGIONS = [
  { id: 'ALL', name: 'All Regions', abbrev: 'ALL' },
  { id: 'AESO', name: 'Alberta Electric System Operator', abbrev: 'AESO' },
  { id: 'AVA', name: 'Avista Corporation', abbrev: 'AVA' },
  { id: 'AZPS', name: 'Arizona Public Service Company', abbrev: 'AZPS' },
  { id: 'BANC', name: 'Balancing Authority of Northern California', abbrev: 'BANC' },
  { id: 'BCHA', name: 'British Columbia Hydro and Power Authority', abbrev: 'BCHA' },
  { id: 'BPAT', name: 'Bonneville Power Administration', abbrev: 'BPAT' },
  { id: 'CENACE', name: 'Centro Nacional de Control de Energía', abbrev: 'CENACE' },
  { id: 'CHPD', name: 'PUD No. 1 of Chelan County', abbrev: 'CHPD' },
  { id: 'CISO', name: 'California Independent System Operator', abbrev: 'CISO' },
  { id: 'DOPD', name: 'PUD No. 1 of Douglas County', abbrev: 'DOPD' },
  { id: 'EPE', name: 'El Paso Electric Company', abbrev: 'EPE' },
  { id: 'GCPD', name: 'Grant County PUD No. 2', abbrev: 'GCPD' },
  { id: 'IID', name: 'Imperial Irrigation District', abbrev: 'IID' },
  { id: 'IPCO', name: 'Idaho Power Company', abbrev: 'IPCO' },
  { id: 'LDWP', name: 'Los Angeles Department of Water and Power', abbrev: 'LDWP' },
  { id: 'NEVP', name: 'Nevada Power Company', abbrev: 'NEVP' },
  { id: 'NWMT', name: 'NorthWestern Corporation', abbrev: 'NWMT' },
  { id: 'PACE', name: 'PacifiCorp East', abbrev: 'PACE' },
  { id: 'PACW', name: 'PacifiCorp West', abbrev: 'PACW' },
  { id: 'PGE', name: 'Portland General Electric Company', abbrev: 'PGE' },
  { id: 'PNM', name: 'Public Service Company of New Mexico', abbrev: 'PNM' },
  { id: 'PSCO', name: 'Public Service Company of Colorado', abbrev: 'PSCO' },
  { id: 'PSEI', name: 'Puget Sound Energy', abbrev: 'PSEI' },
  { id: 'SCL', name: 'Seattle City Light', abbrev: 'SCL' },
  { id: 'SRP', name: 'Salt River Project', abbrev: 'SRP' },
  { id: 'TEPC', name: 'Tucson Electric Power Company', abbrev: 'TEPC' },
  { id: 'TIDC', name: 'Turlock Irrigation District', abbrev: 'TIDC' },
  { id: 'TPWR', name: 'City of Tacoma, Department of Public Utilities', abbrev: 'TPWR' },
  { id: 'WACM', name: 'Western Area Power Administration, Colorado-Missouri Region', abbrev: 'WACM' },
  { id: 'WALC', name: 'Western Area Power Administration, Lower Colorado Region', abbrev: 'WALC' },
  { id: 'WAUW', name: 'Western Area Power Administration, Upper Great Plains West', abbrev: 'WAUW' }
];

const DataCenterImpactPage = () => {
  const navigate = useNavigate();
  
  // State management
  const [viewState, setViewState] = useState(INITIAL_VIEW_STATE);
  const [weccGeojsonData, setWeccGeojsonData] = useState(null);
  const [hoveredObject, setHoveredObject] = useState(null);
  const [clickedObject, setClickedObject] = useState(null);
  
  // Data Center Impact Analysis state
  const [selectedWeccRegion, setSelectedWeccRegion] = useState('ALL');
  const [selectedCaseStudy, setSelectedCaseStudy] = useState('case1');
  const [selectedMetric, setSelectedMetric] = useState('System Operation Cost');
  const [selectedHour, setSelectedHour] = useState(12);
  const [impactData, setImpactData] = useState({});
  const [isLoading, setIsLoading] = useState(false);

  // Load WECC GeoJSON data
  useEffect(() => {
    const loadWeccData = async () => {
      try {
        const response = await fetch('/amin_data/WECC_Balancing_Authorities_-2060174188301432986.geojson');
        const geojsonData = await response.json();
        setWeccGeojsonData(geojsonData);
        console.log('WECC GeoJSON data loaded for Data Center Impact page');
      } catch (error) {
        console.error('Error loading WECC data:', error);
      }
    };
    
    loadWeccData();
  }, []);

  // Load impact data when parameters change
  useEffect(() => {
    if (weccGeojsonData) {
      loadImpactAnalysisData();
    }
  }, [selectedCaseStudy, selectedMetric, selectedHour, selectedWeccRegion, weccGeojsonData]);

  // Calculate impact differences between case studies
  const loadImpactAnalysisData = async () => {
    try {
      setIsLoading(true);
      console.log(`Loading impact data - Case: ${selectedCaseStudy}, Metric: ${selectedMetric}, Hour: ${selectedHour}`);
      
      // Determine which data file to load based on selected metric
      let dataFileName = '';
      if (selectedMetric === 'System Operation Cost') {
        dataFileName = 'Total_Operation_Cost.csv';
      } else if (selectedMetric === 'Price') {
        dataFileName = 'Price_Data.csv';
      } else if (selectedMetric === 'Power Exchange') {
        dataFileName = 'Power_Exchange_Data.csv';
      } else {
        dataFileName = 'Total_Operation_Cost.csv';
      }
      
      // Load the appropriate data file
      const response = await fetch(`/amin_data/manish_amin_modified_data/${dataFileName}`);
      if (!response.ok) {
        console.error(`Failed to load ${dataFileName}`);
        setIsLoading(false);
        return;
      }
      
      const csvText = await response.text();
      const lines = csvText.trim().split('\n');
      const headers = lines[0].split(',');
      const data = {};
      
      // Organize data by case study
      for (let i = 1; i < lines.length; i++) {
        const row = lines[i].split(',');
        const caseStudy = row[0];
        const hour = parseInt(row[1]);
        
        if (!data[caseStudy]) {
          data[caseStudy] = {};
        }
        if (!data[caseStudy][hour]) {
          data[caseStudy][hour] = {};
        }
        
        // Store values for each WECC region
        for (let j = 2; j < headers.length; j++) {
          const regionCode = headers[j].trim();
          const value = parseFloat(row[j]) || 0;
          data[caseStudy][hour][regionCode] = value;
        }
      }
      
      // Define comparison pairs for case studies
      const comparisonPairs = {
        'case1': { current: 'Case_1', baseline: 'Case_0' },
        'case2': { current: 'Case_2', baseline: 'Case_1' },
        'case3': { current: 'Case_3', baseline: 'Case_2' }
      };
      
      const comparison = comparisonPairs[selectedCaseStudy];
      if (!comparison || !data[comparison.current] || !data[comparison.baseline]) {
        console.warn('Case study data not found:', comparison);
        setIsLoading(false);
        return;
      }
      
      // Get regions to analyze
      const regionsToAnalyze = selectedWeccRegion === 'ALL' 
        ? WECC_REGIONS.filter(r => r.id !== 'ALL').map(r => r.id)
        : [selectedWeccRegion];
      
      const results = {};
      
      // Calculate differences for each region
      for (const regionAbbrev of regionsToAnalyze) {
        try {
          const currentValue = data[comparison.current][selectedHour]?.[regionAbbrev] || 0;
          const baselineValue = data[comparison.baseline][selectedHour]?.[regionAbbrev] || 0;
          
          const difference = currentValue - baselineValue;
          const percentChange = baselineValue !== 0 ? ((currentValue - baselineValue) / baselineValue) * 100 : 0;
          
          results[regionAbbrev] = {
            current: currentValue,
            baseline: baselineValue,
            difference: difference,
            percentChange: percentChange,
            region: regionAbbrev,
            metric: selectedMetric,
            hour: selectedHour,
            unit: selectedMetric === 'System Operation Cost' ? '$' : 
                  selectedMetric === 'Price' ? '$/MWh' : 'MW'
          };
        } catch (error) {
          console.warn(`Failed to calculate difference for region ${regionAbbrev}:`, error);
        }
      }
      
      setImpactData(results);
      updateWeccRegionColors(results);
      setIsLoading(false);
      
      console.log(`Impact analysis completed - ${Object.keys(results).length} regions processed`);
    } catch (error) {
      console.error('Error loading impact analysis data:', error);
      setIsLoading(false);
    }
  };

  // Update WECC region colors based on impact analysis results
  const updateWeccRegionColors = (impactResults) => {
    try {
      if (!weccGeojsonData || !impactResults) return;
      
      const differences = Object.values(impactResults).map(result => result.difference);
      const minDiff = Math.min(...differences);
      const maxDiff = Math.max(...differences);
      
      console.log(`Impact Range: ${minDiff.toFixed(2)} to ${maxDiff.toFixed(2)}`);
      
      // Create 20-step linear color scale
      const getLinearImpactColor = (difference) => {
        if (minDiff === maxDiff) return [200, 200, 200, 180];
        
        const normalizedValue = (difference - minDiff) / (maxDiff - minDiff);
        
        const colorSteps = [
          [8, 48, 107], [8, 81, 156], [33, 113, 181], [66, 146, 198], [107, 174, 214],
          [158, 202, 225], [198, 219, 239], [222, 235, 247], [247, 251, 255], [255, 255, 255],
          [255, 245, 240], [254, 224, 210], [252, 187, 161], [252, 146, 114], [251, 106, 74],
          [239, 59, 44], [203, 24, 29], [165, 15, 21], [103, 0, 13], [77, 0, 9]
        ];
        
        const stepIndex = Math.min(Math.floor(normalizedValue * colorSteps.length), colorSteps.length - 1);
        const color = colorSteps[stepIndex];
        
        return [color[0], color[1], color[2], 200];
      };
      
      // Update features with impact colors
      const updatedFeatures = weccGeojsonData.features.map(feature => {
        const regionAbbrev = feature.properties?.BA_Abrev || feature.properties?.BA_CODE;
        const impactResult = impactResults[regionAbbrev];
        
        if (impactResult) {
          const impactColor = getLinearImpactColor(impactResult.difference);
          const colorStep = Math.min(Math.floor(((impactResult.difference - minDiff) / (maxDiff - minDiff)) * 20), 19);
          
          return {
            ...feature,
            properties: {
              ...feature.properties,
              impactColor: impactColor,
              impactValue: impactResult.difference,
              impactData: impactResult,
              colorStep: colorStep
            }
          };
        }
        
        return {
          ...feature,
          properties: {
            ...feature.properties,
            impactColor: [200, 200, 200, 100],
            impactValue: null,
            impactData: null,
            colorStep: null
          }
        };
      });
      
      setWeccGeojsonData({
        ...weccGeojsonData,
        features: updatedFeatures
      });
    } catch (error) {
      console.error('Error updating WECC colors:', error);
    }
  };

  // Tooltip content
  const getTooltip = ({ object }) => {
    if (!object) return null;
    
    const properties = object.properties;
    const baCode = properties.BA_Abrev || 'WECC';
    const baName = properties.BA_Name || 'Western Electricity Coordinating Council Area';
    const impactResult = properties.impactData;
    
    let content = `
      <div style="
        background: linear-gradient(135deg, rgba(21, 101, 192, 0.1) 0%, rgba(21, 101, 192, 0.05) 100%);
        padding: 10px 12px;
        border-radius: 8px;
        border-left: 4px solid #1565c0;
        margin-bottom: 8px;
        box-shadow: 0 4px 12px rgba(0,0,0,0.15);
        max-width: 300px;
      ">
        <div style="font-weight: 700; font-size: 16px; color: #1565c0; margin-bottom: 4px;">
          🏛️ ${baCode}
        </div>
        <div style="font-size: 12px; color: #666; line-height: 1.3; margin-bottom: 8px;">
          ${baName}
        </div>
    `;

    if (impactResult) {
      content += `
        <div style="
          background: linear-gradient(135deg, rgba(25, 118, 210, 0.1) 0%, rgba(25, 118, 210, 0.05) 100%);
          padding: 10px;
          border-radius: 6px;
          border-left: 4px solid #1976d2;
          margin-bottom: 8px;
        ">
          <div style="font-weight: 700; font-size: 14px; color: #1976d2; margin-bottom: 6px;">
            📊 Data Center Impact Analysis
          </div>
          <div style="font-size: 11px; color: #666; line-height: 1.4;">
            <div style="margin-bottom: 3px;"><strong>Case Study:</strong> ${selectedCaseStudy.toUpperCase()}</div>
            <div style="margin-bottom: 3px;"><strong>Metric:</strong> ${impactResult.metric}</div>
            <div style="margin-bottom: 3px;"><strong>Hour:</strong> ${impactResult.hour}</div>
            <div style="margin-bottom: 6px;"><strong>Difference:</strong> 
              <span style="color: ${impactResult.difference > 0 ? '#d32f2f' : '#1976d2'}; font-weight: 700; font-size: 12px;">
                ${impactResult.difference > 0 ? '+' : ''}${impactResult.difference.toFixed(2)} ${impactResult.unit}
              </span>
            </div>
            <div style="font-size: 10px; color: #888; margin-bottom: 3px;">
              ${impactResult.difference > 0 ? '🔴 Higher impact' : '🔵 Lower impact'} vs baseline
            </div>
            <div style="font-size: 9px; color: #666; font-style: italic;">
              Color Step: ${properties.colorStep !== null ? (properties.colorStep + 1) : 'N/A'} of 20
            </div>
          </div>
        </div>
      `;
    }

    content += `
      </div>
    `;

    return {
      html: content,
      style: {
        backgroundColor: 'rgba(255, 255, 255, 0.95)',
        fontSize: '12px',
        fontFamily: '"Inter", sans-serif',
        borderRadius: '8px',
        boxShadow: '0 4px 20px rgba(0,0,0,0.15)',
        border: '1px solid rgba(0,0,0,0.1)',
        backdropFilter: 'blur(8px)'
      }
    };
  };

  // Create map layers
  const layers = [
    new GeoJsonLayer({
      id: 'DataCenterImpactLayer',
      data: weccGeojsonData,
      pickable: true,
      stroked: true,
      filled: true,
      extruded: false,
      lineWidthMinPixels: 2,
      lineWidthMaxPixels: 8,
      getLineColor: [255, 255, 255, 255],
      getLineWidth: 3,
      getFillColor: d => {
        if (d.properties?.impactColor) {
          // Apply hover and click effects
          const baseColor = d.properties.impactColor;
          
          if (d === clickedObject) {
            return [
              Math.max(0, baseColor[0] - 30),
              Math.max(0, baseColor[1] - 30),
              Math.max(0, baseColor[2] - 30),
              255
            ];
          }
          
          if (d === hoveredObject && d !== clickedObject) {
            return [
              Math.min(255, baseColor[0] + 15),
              Math.min(255, baseColor[1] + 15),
              Math.min(255, baseColor[2] + 15),
              255
            ];
          }
          
          return baseColor;
        }
        
        // Default color for regions without impact data
        return [200, 200, 200, 120];
      },
      onHover: (info) => {
        setHoveredObject(info.object);
      },
      onClick: (info) => {
        setClickedObject(info.object);
      },
      updateTriggers: {
        getFillColor: [hoveredObject, clickedObject, impactData]
      }
    })
  ];

  return (
    <div style={{ width: '100vw', height: '100vh', position: 'relative' }}>
      {/* Header */}
      <div style={{
        position: 'absolute',
        top: 0,
        left: 0,
        right: 0,
        height: '60px',
        background: 'linear-gradient(135deg, #1e3a8a 0%, #3b82f6 100%)',
        display: 'flex',
        alignItems: 'center',
        padding: '0 20px',
        zIndex: 1000,
        boxShadow: '0 2px 10px rgba(0,0,0,0.1)'
      }}>
        <button
          onClick={() => navigate('/')}
          style={{
            background: 'rgba(255,255,255,0.1)',
            border: '1px solid rgba(255,255,255,0.2)',
            borderRadius: '8px',
            padding: '8px 16px',
            color: 'white',
            fontSize: '14px',
            fontWeight: '500',
            cursor: 'pointer',
            display: 'flex',
            alignItems: 'center',
            gap: '8px',
            transition: 'all 0.2s'
          }}
          onMouseEnter={(e) => {
            e.target.style.background = 'rgba(255,255,255,0.2)';
          }}
          onMouseLeave={(e) => {
            e.target.style.background = 'rgba(255,255,255,0.1)';
          }}
        >
          <ArrowBackIcon style={{ fontSize: '18px' }} />
          Back to Westmap
        </button>
        
        <div style={{ flex: 1, textAlign: 'center' }}>
          <h1 style={{
            margin: 0,
            color: 'white',
            fontSize: '24px',
            fontWeight: '700',
            textShadow: '0 2px 4px rgba(0,0,0,0.3)'
          }}>
            📊 Data Center Impact Analysis
          </h1>
          <p style={{
            margin: '4px 0 0 0',
            color: 'rgba(255,255,255,0.9)',
            fontSize: '14px',
            fontWeight: '400'
          }}>
            Visualizing WECC Region Impact Differences Across Case Studies
          </p>
        </div>
        
        {isLoading && (
          <div style={{
            background: 'rgba(255,255,255,0.1)',
            padding: '6px 12px',
            borderRadius: '6px',
            color: 'white',
            fontSize: '12px'
          }}>
            Loading data...
          </div>
        )}
      </div>

      {/* Map */}
      <DeckGL
        initialViewState={INITIAL_VIEW_STATE}
        controller={true}
        layers={layers}
        getTooltip={getTooltip}
        onViewStateChange={({ viewState }) => setViewState(viewState)}
      >
        <StaticMap
          mapboxApiAccessToken={MAPBOX_ACCESS_TOKEN}
          mapStyle="mapbox://styles/mapbox/light-v10"
        />
      </DeckGL>

      {/* Controls Panel */}
      <div style={{
        position: 'absolute',
        top: '80px',
        right: '20px',
        width: '320px',
        background: 'rgba(255,255,255,0.95)',
        padding: '20px',
        borderRadius: '12px',
        boxShadow: '0 4px 20px rgba(0,0,0,0.15)',
        backdropFilter: 'blur(8px)',
        border: '1px solid rgba(0,0,0,0.1)',
        zIndex: 1000,
        maxHeight: 'calc(100vh - 100px)',
        overflowY: 'auto'
      }}>
        <div style={{ marginBottom: '20px' }}>
          <h2 style={{
            margin: '0 0 8px 0',
            fontSize: '18px',
            fontWeight: '600',
            color: '#1976d2'
          }}>
            Impact Analysis Controls
          </h2>
          <p style={{
            margin: 0,
            fontSize: '12px',
            color: '#666',
            lineHeight: '1.4'
          }}>
            Adjust parameters to see real-time color changes on WECC regions
          </p>
        </div>

        {/* WECC Region Selection */}
        <div style={{ marginBottom: '16px' }}>
          <label style={{
            display: 'block',
            fontSize: '13px',
            fontWeight: '500',
            marginBottom: '8px',
            color: '#333'
          }}>
            WECC Region:
          </label>
          <select
            value={selectedWeccRegion}
            onChange={(e) => setSelectedWeccRegion(e.target.value)}
            style={{
              width: '100%',
              padding: '10px 12px',
              borderRadius: '6px',
              border: '1px solid #ddd',
              fontSize: '13px',
              background: 'white'
            }}
          >
            {WECC_REGIONS.map(region => (
              <option key={region.id} value={region.id}>
                {region.name}
              </option>
            ))}
          </select>
        </div>

        {/* Case Study Selection */}
        <div style={{ marginBottom: '16px' }}>
          <label style={{
            display: 'block',
            fontSize: '13px',
            fontWeight: '500',
            marginBottom: '8px',
            color: '#333'
          }}>
            Case Study Comparison:
          </label>
          <div style={{ display: 'flex', gap: '8px' }}>
            {['case1', 'case2', 'case3'].map(caseId => (
              <button
                key={caseId}
                onClick={() => setSelectedCaseStudy(caseId)}
                style={{
                  flex: 1,
                  padding: '10px 8px',
                  borderRadius: '6px',
                  border: selectedCaseStudy === caseId ? '2px solid #1976d2' : '1px solid #ddd',
                  background: selectedCaseStudy === caseId ? 'rgba(25, 118, 210, 0.1)' : 'white',
                  fontSize: '12px',
                  fontWeight: selectedCaseStudy === caseId ? '600' : '400',
                  color: selectedCaseStudy === caseId ? '#1976d2' : '#666',
                  cursor: 'pointer',
                  transition: 'all 0.2s'
                }}
              >
                Case {caseId.slice(-1)}
              </button>
            ))}
          </div>
        </div>

        {/* Metric Selection */}
        <div style={{ marginBottom: '16px' }}>
          <label style={{
            display: 'block',
            fontSize: '13px',
            fontWeight: '500',
            marginBottom: '8px',
            color: '#333'
          }}>
            Metric:
          </label>
          <div style={{ display: 'flex', flexDirection: 'column', gap: '6px' }}>
            {[
              { id: 'System Operation Cost', label: 'System Operation Cost ($)' },
              { id: 'Price', label: 'Price ($/MWh)' },
              { id: 'Power Exchange', label: 'Power Exchange (MW)' }
            ].map(metric => (
              <button
                key={metric.id}
                onClick={() => setSelectedMetric(metric.id)}
                style={{
                  padding: '10px 12px',
                  borderRadius: '6px',
                  border: selectedMetric === metric.id ? '2px solid #1976d2' : '1px solid #ddd',
                  background: selectedMetric === metric.id ? 'rgba(25, 118, 210, 0.1)' : 'white',
                  fontSize: '12px',
                  fontWeight: selectedMetric === metric.id ? '600' : '400',
                  color: selectedMetric === metric.id ? '#1976d2' : '#666',
                  cursor: 'pointer',
                  transition: 'all 0.2s',
                  textAlign: 'left'
                }}
              >
                {metric.label}
              </button>
            ))}
          </div>
        </div>

        {/* Hour Selection */}
        <div style={{ marginBottom: '16px' }}>
          <label style={{
            display: 'block',
            fontSize: '13px',
            fontWeight: '500',
            marginBottom: '8px',
            color: '#333'
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
              width: '100%',
              height: '8px',
              borderRadius: '4px',
              background: '#ddd',
              outline: 'none',
              cursor: 'pointer'
            }}
          />
          <div style={{
            display: 'flex',
            justifyContent: 'space-between',
            fontSize: '11px',
            color: '#888',
            marginTop: '6px'
          }}>
            <span>1</span>
            <span>12</span>
            <span>24</span>
          </div>
        </div>

        {/* Color Scale Legend */}
        <div style={{ marginBottom: '16px' }}>
          <label style={{
            display: 'block',
            fontSize: '13px',
            fontWeight: '500',
            marginBottom: '8px',
            color: '#333'
          }}>
            20-Step Linear Color Scale:
          </label>
          <div style={{
            display: 'flex',
            alignItems: 'center',
            gap: '8px',
            marginBottom: '6px'
          }}>
            <span style={{ fontSize: '10px', color: '#08306b', fontWeight: '600' }}>Lowest</span>
            <div style={{
              flex: 1,
              height: '16px',
              background: 'linear-gradient(to right, #08306b, #08519c, #2171b5, #4292c6, #6baed6, #9ecae1, #c6dbef, #deebf7, #f7fbff, #ffffff, #fff5f0, #fee0d2, #fcbba1, #fc9272, #fb6a4a, #ef3b2c, #cb181d, #a50f15, #67000d, #4d000a)',
              borderRadius: '8px',
              border: '1px solid #ddd',
              boxShadow: 'inset 0 1px 3px rgba(0,0,0,0.1)'
            }}></div>
            <span style={{ fontSize: '10px', color: '#4d000a', fontWeight: '600' }}>Highest</span>
          </div>
          <div style={{
            fontSize: '10px',
            color: '#888',
            fontStyle: 'italic',
            textAlign: 'center'
          }}>
            Each region gets one of 20 distinct colors based on relative impact
          </div>
        </div>

        {/* Case Study Explanations */}
        <div style={{
          background: 'rgba(25, 118, 210, 0.05)',
          padding: '12px',
          borderRadius: '8px',
          border: '1px solid rgba(25, 118, 210, 0.1)'
        }}>
          <div style={{
            fontSize: '12px',
            fontWeight: '600',
            color: '#1976d2',
            marginBottom: '8px'
          }}>
            📋 Case Study Comparisons:
          </div>
          <div style={{
            fontSize: '11px',
            color: '#666',
            lineHeight: '1.4'
          }}>
            <div style={{ marginBottom: '2px' }}>• <strong>Case 1:</strong> vs Base Case (Case 0)</div>
            <div style={{ marginBottom: '2px' }}>• <strong>Case 2:</strong> vs Case 1</div>
            <div style={{ marginBottom: '2px' }}>• <strong>Case 3:</strong> vs Case 2</div>
          </div>
        </div>
      </div>

      {/* Statistics Panel */}
      {Object.keys(impactData).length > 0 && (
        <div style={{
          position: 'absolute',
          bottom: '20px',
          left: '20px',
          background: 'rgba(255,255,255,0.95)',
          padding: '16px',
          borderRadius: '12px',
          boxShadow: '0 4px 20px rgba(0,0,0,0.15)',
          backdropFilter: 'blur(8px)',
          border: '1px solid rgba(0,0,0,0.1)',
          minWidth: '300px'
        }}>
          <h3 style={{
            margin: '0 0 12px 0',
            fontSize: '16px',
            fontWeight: '600',
            color: '#1976d2'
          }}>
            📈 Impact Summary
          </h3>
          <div style={{
            fontSize: '12px',
            color: '#666',
            lineHeight: '1.4'
          }}>
            <div><strong>Regions Analyzed:</strong> {Object.keys(impactData).length}</div>
            <div><strong>Positive Impact:</strong> {Object.values(impactData).filter(r => r.difference > 0).length} regions</div>
            <div><strong>Negative Impact:</strong> {Object.values(impactData).filter(r => r.difference < 0).length} regions</div>
            <div><strong>Max Increase:</strong> {Math.max(...Object.values(impactData).map(r => r.difference)).toFixed(2)}</div>
            <div><strong>Max Decrease:</strong> {Math.min(...Object.values(impactData).map(r => r.difference)).toFixed(2)}</div>
          </div>
        </div>
      )}
    </div>
  );
};

export default DataCenterImpactPage;
