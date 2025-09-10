import React, { useState, useEffect } from 'react';
import { useParams, useNavigate, useLocation } from 'react-router-dom';
import { DeckGL } from '@deck.gl/react';
import { GeoJsonLayer } from '@deck.gl/layers';
import { MapView } from '@deck.gl/core';
import { StaticMap } from 'react-map-gl';
import ArrowBackIcon from '@mui/icons-material/ArrowBack';
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer, AreaChart, Area, BarChart, Bar, ComposedChart, PieChart, Pie, Cell } from 'recharts';

// Mapbox token
const MAPBOX_ACCESS_TOKEN = 'pk.eyJ1IjoidXNtYXJ0LXdlc3RtYXAiLCJhIjoiY2tvazV6MzU2MDE4YjJ0bXd5ZDcwdm16ciJ9.q2BIGvGPAJjw1X9CdvyKSA';

// OpenStreetMap style
const OSM_MAP_STYLE = {
    "version": 8,
    "name": "OpenStreetMap",
    "sources": {
        "osm": {
            "type": "raster",
            "tiles": [
                "https://tile.openstreetmap.org/{z}/{x}/{y}.png"
            ],
            "tileSize": 256,
            "attribution": "© OpenStreetMap contributors, University of Utah 2025"
        }
    },
    "layers": [
        {
            "id": "osm",
            "type": "raster",
            "source": "osm"
        }
    ]
};

// Case study configurations - Updated based on requirements
const CASE_STUDIES = [
    {
        id: 'case0',
        name: 'Case Study 0',
        subtitle: 'Base Case',
        description: 'No data center flexibility - Current operational baseline',
        color: '#64748B',
        icon: '📊',
        dataPath: 'Case study_0'
    },
    {
        id: 'case1',
        name: 'Case Study 1',
        subtitle: 'Data Center Temporal Flexibility',
        description: 'Time-shifted data center operations with flexible scheduling',
        color: '#3B82F6',
        icon: '⏰',
        dataPath: 'Case study_1'
    },
    {
        id: 'case2',
        name: 'Case Study 2',
        subtitle: 'Data Center Spatio-Temporal Flexibility',
        description: 'Geographic and temporal load distribution optimization',
        color: '#10B981',
        icon: '🌐',
        dataPath: 'Case study_2'
    },
    {
        id: 'case3',
        name: 'Case Study 3',
        subtitle: 'Data Center Spatio-Temporal Flexibility with On-Site Energy Resources',
        description: 'Advanced flexibility with integrated renewable energy and storage',
        color: '#F59E0B',
        icon: '🔋',
        dataPath: 'Case_study_3'
    },
    {
        id: 'comparison',
        name: 'Comparison',
        subtitle: 'Comparative Analysis',
        description: 'Side-by-side analysis and comparison of all scenarios',
        color: '#8B5CF6',
        icon: '📈',
        dataPath: 'comparison'
    }
];

// Balancing Authorities list with full names for dropdown
const BALANCING_AUTHORITIES = [
    { code: 'AESO', name: 'Alberta Electric System Operator' },
    { code: 'AVA', name: 'Avista Corporation' },
    { code: 'AZPS', name: 'Arizona Public Service Company' },
    { code: 'BANC', name: 'Balancing Authority of Northern California' },
    { code: 'BCHA', name: 'British Columbia Hydro and Power Authority' },
    { code: 'BPAT', name: 'Bonneville Power Administration' },
    { code: 'CENACE', name: 'Centro Nacional de Control de Energía' },
    { code: 'CHPD', name: 'PUD No. 1 of Chelan County' },
    { code: 'CISO', name: 'California Independent System Operator' },
    { code: 'DOPD', name: 'PUD No. 1 of Douglas County' },
    { code: 'EPE', name: 'El Paso Electric Company' },
    { code: 'GCPD', name: 'Grant County PUD No. 2' },
    { code: 'IID', name: 'Imperial Irrigation District' },
    { code: 'IPCO', name: 'Idaho Power Company' },
    { code: 'LDWP', name: 'Los Angeles Department of Water and Power' },
    { code: 'NEVP', name: 'Nevada Power Company' },
    { code: 'NWMT', name: 'NorthWestern Corporation' },
    { code: 'PACE', name: 'PacifiCorp East' },
    { code: 'PACW', name: 'PacifiCorp West' },
    { code: 'PGE', name: 'Portland General Electric Company' },
    { code: 'PNM', name: 'Public Service Company of New Mexico' },
    { code: 'PSCO', name: 'Public Service Company of Colorado' },
    { code: 'PSEI', name: 'Puget Sound Energy' },
    { code: 'SCL', name: 'Seattle City Light' },
    { code: 'SRP', name: 'Salt River Project' },
    { code: 'TEPC', name: 'Tucson Electric Power Company' },
    { code: 'TIDC', name: 'Turlock Irrigation District' },
    { code: 'TPWR', name: 'City of Tacoma, Department of Public Utilities' },
    { code: 'WACM', name: 'Western Area Power Administration, Colorado-Missouri Region' },
    { code: 'WALC', name: 'Western Area Power Administration, Lower Colorado Region' },
    { code: 'WAUW', name: 'Western Area Power Administration, Upper Great Plains West' }
];

// Helper function to format large numbers
const formatValue = (value, unit = '') => {
    if (value >= 1000000) {
        return `${(value / 1000000).toFixed(1)}M${unit}`;
    } else if (value >= 1000) {
        return `${(value / 1000).toFixed(1)}K${unit}`;
    }
    return `${value.toFixed(1)}${unit}`;
};

const AminDetailPage = () => {
    const { fid } = useParams();
    const navigate = useNavigate();
    const location = useLocation();
    
    // State management
    // Get WECC region from URL params or default to NEVP
    const urlParams = new URLSearchParams(location.search);
    const initialWeccRegion = urlParams.get('region') || 'NEVP';
    
    const [selectedWeccRegion, setSelectedWeccRegion] = useState(initialWeccRegion);
    const [balancingAuthority, setBalancingAuthority] = useState(null);
    const [selectedCaseStudy, setSelectedCaseStudy] = useState('case0');
    const [caseStudyData, setCaseStudyData] = useState({});
    const [capacityData, setCapacityData] = useState({});
    const [plantLocations, setPlantLocations] = useState([]);
    const [selectedHour, setSelectedHour] = useState(12);
    const [loading, setLoading] = useState(true);
    const [error, setError] = useState(null);
    
    // Map view state
    const [viewState, setViewState] = useState({
        longitude: -116.5,
        latitude: 37.5,
        zoom: 6,
        minZoom: 3,
        maxZoom: 12,
        pitch: 0,
        bearing: 0
    });

    // Add Inter font
    useEffect(() => {
        const link = document.createElement('link');
        link.href = 'https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700;800&display=swap';
        link.rel = 'stylesheet';
        document.head.appendChild(link);

        return () => {
            if (document.head.contains(link)) {
                document.head.removeChild(link);
            }
        };
    }, []);

    // Set font family
    useEffect(() => {
        const originalFontFamily = document.body.style.fontFamily;
        document.body.style.fontFamily = '"Inter", system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif';
        
        return () => {
            document.body.style.fontFamily = originalFontFamily;
        };
    }, []);

    // Load data on component mount
    useEffect(() => {
        const loadData = async () => {
            try {
                // Load GeoJSON data
                const geojsonResponse = await fetch('/amin_data/WECC_Balancing_Authorities_-2060174188301432986.geojson');
                if (!geojsonResponse.ok) {
                    throw new Error(`HTTP error loading GeoJSON! status: ${geojsonResponse.status}`);
                }
                const geojsonData = await geojsonResponse.json();

                // Load main CSV for shape data
                const csvResponse = await fetch('/amin_data/WECC_Balancing_Authorities_5803277890210865950.csv');
                if (!csvResponse.ok) {
                    throw new Error(`HTTP error loading main CSV! status: ${csvResponse.status}`);
                }
                const csvText = await csvResponse.text();

                // Parse main CSV for shape data
                const csvLines = csvText.split('\n');
                const csvData = {};

                for (let i = 1; i < csvLines.length; i++) {
                    const line = csvLines[i].trim();
                    if (line) {
                        const values = line.split(',');
                        const fidValue = parseInt(values[0]);
                        csvData[fidValue] = {
                            FID: fidValue,
                            BA_Abrev: values[1],
                            BA_Name: values[2].replace(/"/g, ''),
                            Shape_Leng: parseFloat(values[3]),
                            Shape__Area: parseFloat(values[4]),
                            Shape__Length: parseFloat(values[5]),
                            GlobalID: values[6]
                        };
                    }
                }

                if (fid) {
                    // Find the specific balancing authority by FID
                    const targetFid = parseInt(fid);
                    const feature = geojsonData.features.find(f => f.properties.FID === targetFid);

                    if (!feature) {
                        throw new Error(`Balancing authority with FID ${fid} not found`);
                    }

                    // Merge data from all sources
                    const baData = {
                        ...feature.properties,
                        ...csvData[targetFid]
                    };

                    // Update selected WECC region based on the loaded data
                    if (baData.BA_Abrev && baData.BA_Abrev !== selectedWeccRegion) {
                        setSelectedWeccRegion(baData.BA_Abrev);
                    }

                    // Calculate bounds for the specific feature to center the map
                    let coordinates;
                    if (feature.geometry.type === 'Polygon') {
                        coordinates = feature.geometry.coordinates[0];
                    } else if (feature.geometry.type === 'MultiPolygon') {
                        coordinates = feature.geometry.coordinates[0][0];
                    } else {
                        throw new Error(`Unsupported geometry type: ${feature.geometry.type}`);
                    }

                    let minLng = Infinity, maxLng = -Infinity, minLat = Infinity, maxLat = -Infinity;
                    let validCoordinates = false;

                    coordinates.forEach(coord => {
                        const [lng, lat] = coord;
                        if (typeof lng === 'number' && typeof lat === 'number' &&
                            !isNaN(lng) && !isNaN(lat) &&
                            lng >= -180 && lng <= 180 &&
                            lat >= -90 && lat <= 90) {
                            minLng = Math.min(minLng, lng);
                            maxLng = Math.max(maxLng, lng);
                            minLat = Math.min(minLat, lat);
                            maxLat = Math.max(maxLat, lat);
                            validCoordinates = true;
                        }
                    });

                    if (validCoordinates) {
                        const centerLng = (minLng + maxLng) / 2;
                        const centerLat = (minLat + maxLat) / 2;
                        const lngDiff = maxLng - minLng;
                        const latDiff = maxLat - minLat;
                        const maxDiff = Math.max(lngDiff, latDiff);

                        let zoom = 6;
                        if (maxDiff < 1) zoom = 8;
                        else if (maxDiff < 2) zoom = 7;
                        else if (maxDiff < 4) zoom = 6;
                        else zoom = 5;

                        setViewState(prev => ({
                            ...prev,
                            longitude: centerLng,
                            latitude: centerLat,
                            zoom: zoom
                        }));
                    }

                    setBalancingAuthority({
                        ...baData,
                        feature: feature,
                        bounds: validCoordinates ? { minLng, maxLng, minLat, maxLat } : null
                    });
                } else if (selectedWeccRegion) {
                    // Load data directly for WECC region without fid
                    const selectedBA = BALANCING_AUTHORITIES.find(ba => ba.code === selectedWeccRegion);
                    if (selectedBA) {
                        setBalancingAuthority({
                            BA_Abrev: selectedBA.code,
                            BA_Name: selectedBA.name,
                            FID: null,
                            feature: null,
                            bounds: null
                        });
                    }
                }

                // Load all case study data
                loadAllCaseStudyData();
                setLoading(false);
            } catch (err) {
                console.error('Error loading data:', err);
                setError(err.message);
                setLoading(false);
            }
        };

        if (fid || selectedWeccRegion) {
            loadData();
        }
    }, [fid, selectedWeccRegion]);

    // Comprehensive data loading function
    const loadAllCaseStudyData = async () => {
        try {
            setLoading(true);
            const allData = {};

            // Load data for all case studies
            for (const caseStudy of CASE_STUDIES) {
                if (caseStudy.id === 'comparison') continue; // Handle comparison separately
                
                console.log(`Loading data for ${caseStudy.name}...`);
                allData[caseStudy.id] = await loadCaseStudyData(caseStudy);
            }

            // Load comparison data (combines all case studies)
            allData.comparison = await loadComparisonData();

            setCaseStudyData(allData);
            setLoading(false);
            console.log('All case study data loaded:', allData);
        } catch (error) {
            console.error('Error loading case study data:', error);
            setError(error.message);
            setLoading(false);
        }
    };

    // Load data for a specific case study
    const loadCaseStudyData = async (caseStudy) => {
        const data = {
            prices: [],
            demand: [],
            dataCenterDemand: [],
            generation: [],
            costs: [],
            interchange: [],
            flexibility: []
        };

        try {
            // 1. Load BA Prices
            const priceResponse = await fetch(`/amin_data/manish_amin_modified_data/${caseStudy.dataPath}/BA Price/lmp_by_ba_hour.csv`);
            if (priceResponse.ok) {
                const priceText = await priceResponse.text();
                data.prices = parsePriceData(priceText, selectedWeccRegion);
            }

            // 2. Load BA Demand
            const demandResponse = await fetch(`/amin_data/manish_amin_modified_data/${caseStudy.dataPath}/Balancing Authority Demand (MW).csv`);
            if (demandResponse.ok) {
                const demandText = await demandResponse.text();
                data.demand = parseDemandData(demandText, selectedWeccRegion);
            }

            // 3. Load Data Center Demand
            const dcDemandResponse = await fetch(`/amin_data/manish_amin_modified_data/${caseStudy.dataPath}/Data Center Demand (MW).csv`);
            if (dcDemandResponse.ok) {
                const dcDemandText = await dcDemandResponse.text();
                data.dataCenterDemand = parseDataCenterDemandData(dcDemandText, selectedWeccRegion);
            }

            // 4. Load Generation Data
            const genResponse = await fetch(`/amin_data/manish_amin_modified_data/${caseStudy.dataPath}/Balancing Authority Power Generation/${selectedWeccRegion}_generation_by_fuel.csv`);
            if (genResponse.ok) {
                const genText = await genResponse.text();
                data.generation = parseGenerationData(genText);
                data.interchange = parseInterchangeData(genText); // Extract import/export from generation data
            }

            // 5. Load Operation Costs
            const costsResponse = await fetch(`/amin_data/manish_amin_modified_data/${caseStudy.dataPath}/Balancing Authority Hourly Operation Costs/${selectedWeccRegion}_hourly_operation_costs.csv`);
            if (costsResponse.ok) {
                const costsText = await costsResponse.text();
                data.costs = parseCostsData(costsText);
            }

            // 6. Load Data Center Flexibility (for case studies 1, 2, 3)
            if (caseStudy.id !== 'case0') {
                const flexPath = caseStudy.id === 'case1' 
                    ? `Data Center Flexibility/dc_flexibility_hourly_by_ba.csv`
                    : `Data Center Flexibility/dc_flexibility_hourly_by_ba.csv`;
                
                const flexResponse = await fetch(`/amin_data/manish_amin_modified_data/${caseStudy.dataPath}/${flexPath}`);
                if (flexResponse.ok) {
                    const flexText = await flexResponse.text();
                    data.flexibility = parseFlexibilityData(flexText, selectedWeccRegion);
                }
            }

            return data;
        } catch (error) {
            console.error(`Error loading ${caseStudy.name} data:`, error);
            return data;
        }
    };

    // Load comparison data (combines all case studies) - Fixed to work independently
    const loadComparisonData = async () => {
        const comparisonData = {
            prices: [],
            demand: [],
            dataCenterDemand: [],
            generation: [],
            costs: [],
            interchange: [],
            flexibility: [],
            totalSystemCosts: []
        };

        try {
            // Load total system operation costs
            const totalCostsResponse = await fetch('/amin_data/manish_amin_modified_data/Total_Operation_Cost.csv');
            if (totalCostsResponse.ok) {
                const totalCostsText = await totalCostsResponse.text();
                comparisonData.totalSystemCosts = parseTotalSystemCosts(totalCostsText, selectedWeccRegion);
            }

            // Load data directly for each case study for comparison
            const caseStudyPaths = ['Case study_0', 'Case study_1', 'Case study_2', 'Case_study_3'];
            const caseNames = ['Case Study 0', 'Case Study 1', 'Case Study 2', 'Case Study 3'];
            
            for (let hour = 1; hour <= 24; hour++) {
                const priceData = { hour };
                const demandData = { hour };
                const dcDemandData = { hour };
                const genData = { hour };
                const costData = { hour };
                const interchangeData = { hour };
                const flexData = { hour };

                // Load data for each case study for this comparison
                for (let i = 0; i < caseStudyPaths.length; i++) {
                    const casePath = caseStudyPaths[i];
                    const caseName = caseNames[i];
                    
                    try {
                        // Load price data
                        const priceResponse = await fetch(`/amin_data/manish_amin_modified_data/${casePath}/BA Price/lmp_by_ba_hour.csv`);
                        if (priceResponse.ok) {
                            const priceText = await priceResponse.text();
                            const priceDataParsed = parsePriceData(priceText, selectedWeccRegion);
                            const hourPrice = priceDataParsed.find(d => d.hour === hour);
                            priceData[caseName] = hourPrice?.price || 0;
                        }

                        // Load demand data
                        const demandResponse = await fetch(`/amin_data/manish_amin_modified_data/${casePath}/Balancing Authority Demand (MW).csv`);
                        if (demandResponse.ok) {
                            const demandText = await demandResponse.text();
                            const demandDataParsed = parseDemandData(demandText, selectedWeccRegion);
                            const hourDemand = demandDataParsed.find(d => d.hour === hour);
                            demandData[caseName] = hourDemand?.demand || 0;
                        }

                        // Load generation data
                        const genResponse = await fetch(`/amin_data/manish_amin_modified_data/${casePath}/Balancing Authority Power Generation/${selectedWeccRegion}_generation_by_fuel.csv`);
                        if (genResponse.ok) {
                            const genText = await genResponse.text();
                            const genDataParsed = parseGenerationData(genText);
                            const hourGen = genDataParsed.find(d => d.hour === hour);
                            genData[caseName] = hourGen?.total || 0;
                            
                            const interchangeDataParsed = parseInterchangeData(genText);
                            const hourInterchange = interchangeDataParsed.find(d => d.hour === hour);
                            interchangeData[caseName] = hourInterchange?.netInterchange || 0;
                        }

                        // Load cost data
                        const costsResponse = await fetch(`/amin_data/manish_amin_modified_data/${casePath}/Balancing Authority Hourly Operation Costs/${selectedWeccRegion}_hourly_operation_costs.csv`);
                        if (costsResponse.ok) {
                            const costsText = await costsResponse.text();
                            const costsDataParsed = parseCostsData(costsText);
                            const hourCost = costsDataParsed.find(d => d.hour === hour);
                            costData[caseName] = hourCost?.totalCosts || 0;
                        }

                        // Load flexibility data (for cases 1, 2, 3)
                        if (i > 0) { // Skip case 0
                            const flexPath = i === 1 ? 'Data Center Flexibility' : 'Data Center Flexibility';
                            const flexResponse = await fetch(`/amin_data/manish_amin_modified_data/${casePath}/${flexPath}/dc_flexibility_hourly_by_ba.csv`);
                            if (flexResponse.ok) {
                                const flexText = await flexResponse.text();
                                const flexDataParsed = parseFlexibilityData(flexText, selectedWeccRegion);
                                const hourFlex = flexDataParsed.find(d => d.hour === hour);
                                flexData[caseName] = hourFlex?.flexibility || 0;
                            }
                        }
                    } catch (error) {
                        console.warn(`Error loading data for ${caseName} hour ${hour}:`, error);
                    }
                }

                comparisonData.prices.push(priceData);
                comparisonData.demand.push(demandData);
                comparisonData.dataCenterDemand.push(dcDemandData);
                comparisonData.generation.push(genData);
                comparisonData.costs.push(costData);
                comparisonData.interchange.push(interchangeData);
                comparisonData.flexibility.push(flexData);
            }

            return comparisonData;
        } catch (error) {
            console.error('Error loading comparison data:', error);
            return comparisonData;
        }
    };

    // Data parsing functions
    const parsePriceData = (csvText, region) => {
        const lines = csvText.trim().split('\n');
        const headers = lines[0].split(',');
        const data = [];
        
        const regionColumn = headers.findIndex(h => h.includes(`${region}_LMP_$/MWh`));
        if (regionColumn === -1) return data;

        for (let i = 1; i < lines.length; i++) {
            const values = lines[i].split(',');
            const hour = parseInt(values[0]);
            const price = parseFloat(values[regionColumn]) || 0;
            data.push({ hour, price });
        }
        return data;
    };

    const parseDemandData = (csvText, region) => {
        const lines = csvText.trim().split('\n');
        const headers = lines[0].split(',');
        const data = [];
        
        const regionColumn = headers.findIndex(h => h === region);
        if (regionColumn === -1) return data;

        for (let i = 1; i < lines.length; i++) {
            const values = lines[i].split(',');
            const hour = parseInt(values[3]); // Period column is the hour
            const demand = parseFloat(values[regionColumn]) || 0;
            data.push({ hour, demand });
        }
        return data;
    };

    const parseDataCenterDemandData = (csvText, region) => {
        const lines = csvText.trim().split('\n');
        const data = [];
        
        for (let i = 1; i < lines.length; i++) {
            const values = lines[i].split(',');
            if (values[1] === region) {
                const hour = parseInt(values[0]);
                const demand = parseFloat(values[7]) || 0; // Base_Total_Load_MW column
                data.push({ hour, demand });
            }
        }
        
        // Sort by hour and ensure we have all 24 hours
        data.sort((a, b) => a.hour - b.hour);
        return data;
    };

    const parseGenerationData = (csvText) => {
        const lines = csvText.trim().split('\n');
        const data = [];
        
        for (let i = 1; i < lines.length; i++) {
            const values = lines[i].split(',');
            const hour = parseInt(values[0]);
            const genData = {
                hour,
                naturalGas: parseFloat(values[1]) || 0,
                geothermal: parseFloat(values[2]) || 0,
                biomass: parseFloat(values[3]) || 0,
                nuclear: parseFloat(values[4]) || 0,
                coal: parseFloat(values[5]) || 0,
                wind: parseFloat(values[6]) || 0,
                solar: parseFloat(values[7]) || 0,
                hydro: parseFloat(values[8]) || 0,
                battery: parseFloat(values[9]) || 0,
                total: 0
            };
            
            // Calculate total generation
            genData.total = genData.naturalGas + genData.geothermal + genData.biomass + 
                           genData.nuclear + genData.coal + genData.wind + 
                           genData.solar + genData.hydro + genData.battery;
            
            data.push(genData);
        }
        return data;
    };

    const parseInterchangeData = (csvText) => {
        const lines = csvText.trim().split('\n');
        const data = [];
        
        for (let i = 1; i < lines.length; i++) {
            const values = lines[i].split(',');
            const hour = parseInt(values[0]);
            const netInterchange = parseFloat(values[10]) || 0; // IMPORT/EXPORT_MW column
            data.push({ hour, netInterchange });
        }
        return data;
    };

    const parseCostsData = (csvText) => {
        const lines = csvText.trim().split('\n');
        const data = [];
        
        for (let i = 1; i < lines.length; i++) {
            const values = lines[i].split(',');
            const hour = parseInt(values[0]);
            const totalCosts = parseFloat(values[6]) || 0; // BA_Total_Hourly_Cost_$ column
            data.push({ hour, totalCosts: totalCosts / 1000000 }); // Convert to millions
        }
        return data;
    };

    const parseFlexibilityData = (csvText, region) => {
        const lines = csvText.trim().split('\n');
        const headers = lines[0].split(',');
        const data = [];
        
        const regionColumn = headers.findIndex(h => h === region);
        if (regionColumn === -1) return data;

        for (let i = 1; i < lines.length; i++) {
            const values = lines[i].split(',');
            const hour = parseInt(values[0]);
            const flexibility = parseFloat(values[regionColumn]) || 0;
            data.push({ hour, flexibility });
        }
        return data;
    };

    const parseTotalSystemCosts = (csvText, region) => {
        const lines = csvText.trim().split('\n');
        const headers = lines[0].split(',');
        const data = [];
        
        const regionColumn = headers.findIndex(h => h === region);
        if (regionColumn === -1) return data;

        for (let i = 1; i < lines.length; i++) {
            const values = lines[i].split(',');
            const caseStudy = values[0];
            const hour = parseInt(values[1]);
            const cost = parseFloat(values[regionColumn]) || 0;
            
            data.push({ 
                hour, 
                caseStudy, 
                cost: cost / 1000000 // Convert to millions
            });
        }
        return data;
    };

    // Map rendering function
    const renderMap = () => {
        if (!balancingAuthority && !selectedWeccRegion) return null;
        
        // If we don't have balancingAuthority but have selectedWeccRegion, show a placeholder
        if (!balancingAuthority) {
            return (
                <div style={{
                    width: '100%',
                    height: '100%',
                    display: 'flex',
                    alignItems: 'center',
                    justifyContent: 'center',
                    background: 'linear-gradient(135deg, #f8fafc 0%, #e2e8f0 100%)',
                    borderRadius: '16px',
                    border: '2px dashed #cbd5e1',
                    color: '#64748b',
                    fontSize: '16px',
                    fontWeight: 500
                }}>
                    <div style={{ textAlign: 'center' }}>
                        <div style={{ fontSize: '48px', marginBottom: '12px' }}>🗺️</div>
                        <div>Map view for {selectedWeccRegion}</div>
                        <div style={{ fontSize: '14px', marginTop: '8px', opacity: 0.7 }}>
                            Geographic boundary data loading...
                        </div>
                    </div>
                </div>
            );
        }

        const layers = [
            new GeoJsonLayer({
                id: 'selected-balancing-authority',
                data: {
                    type: 'FeatureCollection',
                    features: [balancingAuthority.feature]
                },
                pickable: false,
                stroked: true,
                filled: true,
                extruded: false,
                wireframe: false,
                getLineColor: [30, 90, 150, 255],
                getLineWidth: 3,
                lineWidthMinPixels: 2,
                getFillColor: [100, 149, 237, 180]
            })
        ];

        // Removed power plant and data center bubbles for clean visualization - showing only region boundaries

        return (
            <div style={{ width: '100%', height: '100%', position: 'relative', borderRadius: '16px', overflow: 'hidden' }}>
                <DeckGL
                    viewState={viewState}
                    controller={true}
                    layers={layers}
                    views={new MapView({ id: 'map' })}
                    width="100%"
                    height="100%"
                    onViewStateChange={({ viewState }) => setViewState(viewState)}
                >
                    <StaticMap
                        mapboxApiAccessToken={MAPBOX_ACCESS_TOKEN}
                        mapStyle={OSM_MAP_STYLE}
                    />
                </DeckGL>
                <div style={{
                    position: 'absolute',
                    bottom: '10px',
                    right: '10px',
                    background: 'rgba(255, 255, 255, 0.9)',
                    padding: '8px 12px',
                    borderRadius: '6px',
                    fontSize: '12px',
                    fontWeight: 500
                }}>
                    © University of Utah 2025
                </div>
            </div>
        );
    };

    // Chart rendering functions
    const renderPriceChart = (caseId) => {
        const data = caseId === 'comparison' 
            ? caseStudyData.comparison?.prices || []
            : caseStudyData[caseId]?.prices || [];
        
        if (!data.length) return <div>No price data available</div>;

        return (
            <ResponsiveContainer width="100%" height={400}>
                <LineChart data={data} margin={{ top: 30, right: 40, left: 60, bottom: 100 }}>
                    <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                    <XAxis
                        dataKey="hour"
                        label={{ 
                            value: 'Hour of Day', 
                            position: 'insideBottom', 
                            offset: -20, 
                            style: { fontSize: '14px', fontWeight: 600, textAnchor: 'middle', fill: '#1e293b' } 
                        }}
                        tick={{ fontSize: 12, fill: '#64748b' }}
                        stroke="#94a3b8"
                    />
                    <YAxis
                        label={{ 
                            value: 'Price ($/MWh)', 
                            angle: -90, 
                            position: 'insideLeft', 
                            style: { fontSize: '14px', fontWeight: 600, textAnchor: 'middle', fill: '#1e293b' } 
                        }}
                        tick={{ fontSize: 12, fill: '#64748b' }}
                        stroke="#94a3b8"
                    />
                    <Tooltip contentStyle={{ backgroundColor: 'rgba(255, 255, 255, 0.98)', border: '1px solid #e2e8f0', borderRadius: '12px' }} />
                    <Legend />
                    {caseId === 'comparison' ? (
                        <>
                            <Line type="monotone" dataKey="Case Study 0" stroke="#64748B" strokeWidth={3} dot={{ fill: '#64748B', r: 4 }} />
                            <Line type="monotone" dataKey="Case Study 1" stroke="#3B82F6" strokeWidth={3} dot={{ fill: '#3B82F6', r: 4 }} />
                            <Line type="monotone" dataKey="Case Study 2" stroke="#10B981" strokeWidth={3} dot={{ fill: '#10B981', r: 4 }} />
                            <Line type="monotone" dataKey="Case Study 3" stroke="#F59E0B" strokeWidth={3} dot={{ fill: '#F59E0B', r: 4 }} />
                        </>
                    ) : (
                        <Line type="monotone" dataKey="price" stroke="#3b82f6" strokeWidth={3} dot={{ fill: '#3b82f6', r: 4 }} name="Price" />
                    )}
                </LineChart>
            </ResponsiveContainer>
        );
    };

    const renderDemandChart = (caseId) => {
        const data = caseId === 'comparison' 
            ? caseStudyData.comparison?.demand || []
            : caseStudyData[caseId]?.demand || [];
        
        if (!data.length) return <div>No demand data available</div>;

        return (
            <ResponsiveContainer width="100%" height={400}>
                <AreaChart data={data} margin={{ top: 30, right: 40, left: 60, bottom: 100 }}>
                    <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                    <XAxis
                        dataKey="hour"
                        label={{ value: 'Hour of Day', position: 'insideBottom', offset: -20, style: { fontSize: '14px', fontWeight: 600, textAnchor: 'middle', fill: '#1e293b' } }}
                        tick={{ fontSize: 12, fill: '#64748b' }}
                        stroke="#94a3b8"
                    />
                    <YAxis
                        label={{ value: 'Demand (MW)', angle: -90, position: 'insideLeft', style: { fontSize: '14px', fontWeight: 600, textAnchor: 'middle', fill: '#1e293b' } }}
                        tick={{ fontSize: 12, fill: '#64748b' }}
                        stroke="#94a3b8"
                    />
                    <Tooltip contentStyle={{ backgroundColor: 'rgba(255, 255, 255, 0.98)', border: '1px solid #e2e8f0', borderRadius: '12px' }} />
                    <Legend />
                    {caseId === 'comparison' ? (
                        <>
                            <Area type="monotone" dataKey="Case Study 0" stackId="1" stroke="#64748B" fill="rgba(100, 116, 139, 0.6)" />
                            <Area type="monotone" dataKey="Case Study 1" stackId="2" stroke="#3B82F6" fill="rgba(59, 130, 246, 0.6)" />
                            <Area type="monotone" dataKey="Case Study 2" stackId="3" stroke="#10B981" fill="rgba(16, 185, 129, 0.6)" />
                            <Area type="monotone" dataKey="Case Study 3" stackId="4" stroke="#F59E0B" fill="rgba(245, 158, 11, 0.6)" />
                        </>
                    ) : (
                        <Area type="monotone" dataKey="demand" stroke="#10b981" fill="rgba(16, 185, 129, 0.6)" name="Demand" />
                    )}
                </AreaChart>
            </ResponsiveContainer>
        );
    };

    const renderDataCenterDemandChart = (caseId) => {
        const data = caseId === 'comparison' 
            ? caseStudyData.comparison?.dataCenterDemand || []
            : caseStudyData[caseId]?.dataCenterDemand || [];
        
        if (!data.length) return <div>No data center demand data available</div>;

        return (
            <ResponsiveContainer width="100%" height={400}>
                <AreaChart data={data} margin={{ top: 30, right: 40, left: 60, bottom: 100 }}>
                    <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                    <XAxis
                        dataKey="hour"
                        label={{ value: 'Hour of Day', position: 'insideBottom', offset: -20, style: { fontSize: '14px', fontWeight: 600, textAnchor: 'middle', fill: '#1e293b' } }}
                        tick={{ fontSize: 12, fill: '#64748b' }}
                        stroke="#94a3b8"
                    />
                    <YAxis
                        label={{ value: 'Data Center Demand (MW)', angle: -90, position: 'insideLeft', style: { fontSize: '14px', fontWeight: 600, textAnchor: 'middle', fill: '#1e293b' } }}
                        tick={{ fontSize: 12, fill: '#64748b' }}
                        stroke="#94a3b8"
                    />
                    <Tooltip contentStyle={{ backgroundColor: 'rgba(255, 255, 255, 0.98)', border: '1px solid #e2e8f0', borderRadius: '12px' }} />
                    <Legend />
                    {caseId === 'comparison' ? (
                        <>
                            <Area type="monotone" dataKey="Case Study 0" stroke="#64748B" fill="rgba(100, 116, 139, 0.6)" />
                            <Area type="monotone" dataKey="Case Study 1" stroke="#3B82F6" fill="rgba(59, 130, 246, 0.6)" />
                            <Area type="monotone" dataKey="Case Study 2" stroke="#10B981" fill="rgba(16, 185, 129, 0.6)" />
                            <Area type="monotone" dataKey="Case Study 3" stroke="#F59E0B" fill="rgba(245, 158, 11, 0.6)" />
                        </>
                    ) : (
                        <Area type="monotone" dataKey="demand" stroke="#f59e0b" fill="rgba(245, 158, 11, 0.6)" name="Data Center Demand" />
                    )}
                </AreaChart>
            </ResponsiveContainer>
        );
    };

    const renderGenerationChart = (caseId) => {
        const data = caseId === 'comparison' 
            ? caseStudyData.comparison?.generation || []
            : caseStudyData[caseId]?.generation || [];
        
        if (!data.length) return <div>No generation data available</div>;

        return (
            <ResponsiveContainer width="100%" height={450}>
                <AreaChart data={data} margin={{ top: 30, right: 40, left: 60, bottom: 100 }}>
                    <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                    <XAxis
                        dataKey="hour"
                        label={{ value: 'Hour of Day', position: 'insideBottom', offset: -20, style: { fontSize: '14px', fontWeight: 600, textAnchor: 'middle', fill: '#1e293b' } }}
                        tick={{ fontSize: 12, fill: '#64748b' }}
                        stroke="#94a3b8"
                    />
                    <YAxis
                        label={{ value: 'Generation (MW)', angle: -90, position: 'insideLeft', style: { fontSize: '14px', fontWeight: 600, textAnchor: 'middle', fill: '#1e293b' } }}
                        tick={{ fontSize: 12, fill: '#64748b' }}
                        stroke="#94a3b8"
                    />
                    <Tooltip contentStyle={{ backgroundColor: 'rgba(255, 255, 255, 0.98)', border: '1px solid #e2e8f0', borderRadius: '12px' }} />
                    <Legend />
                    
                    {caseId === 'comparison' ? (
                        <>
                            <Area type="monotone" dataKey="Case Study 0" stroke="#64748B" fill="rgba(100, 116, 139, 0.6)" />
                            <Area type="monotone" dataKey="Case Study 1" stroke="#3B82F6" fill="rgba(59, 130, 246, 0.6)" />
                            <Area type="monotone" dataKey="Case Study 2" stroke="#10B981" fill="rgba(16, 185, 129, 0.6)" />
                            <Area type="monotone" dataKey="Case Study 3" stroke="#F59E0B" fill="rgba(245, 158, 11, 0.6)" />
                        </>
                    ) : (
                        <>
                            <Area type="monotone" dataKey="wind" stackId="1" stroke="#10b981" fill="rgba(16, 185, 129, 0.8)" name="Wind" />
                            <Area type="monotone" dataKey="solar" stackId="1" stroke="#f59e0b" fill="rgba(245, 158, 11, 0.8)" name="Solar" />
                            <Area type="monotone" dataKey="hydro" stackId="1" stroke="#3b82f6" fill="rgba(59, 130, 246, 0.8)" name="Hydro" />
                            <Area type="monotone" dataKey="nuclear" stackId="1" stroke="#ef4444" fill="rgba(239, 68, 68, 0.8)" name="Nuclear" />
                            <Area type="monotone" dataKey="naturalGas" stackId="1" stroke="#8b5cf6" fill="rgba(139, 92, 246, 0.8)" name="Natural Gas" />
                            <Area type="monotone" dataKey="coal" stackId="1" stroke="#6b7280" fill="rgba(107, 114, 128, 0.8)" name="Coal" />
                            <Area type="monotone" dataKey="geothermal" stackId="1" stroke="#dc2626" fill="rgba(220, 38, 38, 0.8)" name="Geothermal" />
                            <Area type="monotone" dataKey="biomass" stackId="1" stroke="#059669" fill="rgba(5, 150, 105, 0.8)" name="Biomass" />
                            <Area type="monotone" dataKey="battery" stackId="1" stroke="#ec4899" fill="rgba(236, 72, 153, 0.9)" name="🔋 Battery" strokeWidth={3} />
                        </>
                    )}
                </AreaChart>
            </ResponsiveContainer>
        );
    };

    const renderCostsChart = (caseId) => {
        const data = caseId === 'comparison' 
            ? caseStudyData.comparison?.costs || []
            : caseStudyData[caseId]?.costs || [];
        
        if (!data.length) return <div>No costs data available</div>;

        return (
            <ResponsiveContainer width="100%" height={400}>
                <AreaChart data={data} margin={{ top: 30, right: 40, left: 60, bottom: 100 }}>
                    <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                    <XAxis
                        dataKey="hour"
                        label={{ value: 'Hour of Day', position: 'insideBottom', offset: -20, style: { fontSize: '14px', fontWeight: 600, textAnchor: 'middle', fill: '#1e293b' } }}
                        tick={{ fontSize: 12, fill: '#64748b' }}
                        stroke="#94a3b8"
                    />
                    <YAxis
                        label={{ value: 'Operation Cost (M$)', angle: -90, position: 'insideLeft', style: { fontSize: '14px', fontWeight: 600, textAnchor: 'middle', fill: '#1e293b' } }}
                        tick={{ fontSize: 12, fill: '#64748b' }}
                        stroke="#94a3b8"
                    />
                    <Tooltip 
                        contentStyle={{ backgroundColor: 'rgba(255, 255, 255, 0.98)', border: '1px solid #e2e8f0', borderRadius: '12px' }}
                        formatter={(value) => [`$${value.toFixed(2)}M`, 'Cost']}
                    />
                    <Legend />
                    {caseId === 'comparison' ? (
                        <>
                            <Area type="monotone" dataKey="Case Study 0" stroke="#64748B" fill="rgba(100, 116, 139, 0.6)" />
                            <Area type="monotone" dataKey="Case Study 1" stroke="#3B82F6" fill="rgba(59, 130, 246, 0.6)" />
                            <Area type="monotone" dataKey="Case Study 2" stroke="#10B981" fill="rgba(16, 185, 129, 0.6)" />
                            <Area type="monotone" dataKey="Case Study 3" stroke="#F59E0B" fill="rgba(245, 158, 11, 0.6)" />
                        </>
                    ) : (
                        <Area type="monotone" dataKey="totalCosts" stroke="#1d4ed8" fill="rgba(29, 78, 216, 0.6)" name="Total Cost" />
                    )}
                </AreaChart>
            </ResponsiveContainer>
        );
    };

    const renderInterchangeChart = (caseId) => {
        const data = caseId === 'comparison' 
            ? caseStudyData.comparison?.interchange || []
            : caseStudyData[caseId]?.interchange || [];
        
        if (!data.length) return <div>No interchange data available</div>;

        return (
            <ResponsiveContainer width="100%" height={400}>
                <AreaChart data={data} margin={{ top: 30, right: 40, left: 60, bottom: 100 }}>
                    <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                    <XAxis
                        dataKey="hour"
                        label={{ value: 'Hour of Day', position: 'insideBottom', offset: -20, style: { fontSize: '14px', fontWeight: 600, textAnchor: 'middle', fill: '#1e293b' } }}
                        tick={{ fontSize: 12, fill: '#64748b' }}
                        stroke="#94a3b8"
                    />
                    <YAxis
                        label={{ value: 'Net Interchange (MW)', angle: -90, position: 'insideLeft', style: { fontSize: '14px', fontWeight: 600, textAnchor: 'middle', fill: '#1e293b' } }}
                        tick={{ fontSize: 12, fill: '#64748b' }}
                        stroke="#94a3b8"
                    />
                    <Tooltip contentStyle={{ backgroundColor: 'rgba(255, 255, 255, 0.98)', border: '1px solid #e2e8f0', borderRadius: '12px' }} />
                    <Legend />
                    {caseId === 'comparison' ? (
                        <>
                            <Area type="monotone" dataKey="Case Study 0" stroke="#64748B" fill="rgba(100, 116, 139, 0.6)" />
                            <Area type="monotone" dataKey="Case Study 1" stroke="#3B82F6" fill="rgba(59, 130, 246, 0.6)" />
                            <Area type="monotone" dataKey="Case Study 2" stroke="#10B981" fill="rgba(16, 185, 129, 0.6)" />
                            <Area type="monotone" dataKey="Case Study 3" stroke="#F59E0B" fill="rgba(245, 158, 11, 0.6)" />
                        </>
                    ) : (
                        <Area type="monotone" dataKey="netInterchange" stroke="#8b5cf6" fill="rgba(139, 92, 246, 0.6)" name="Net Interchange" />
                    )}
                </AreaChart>
            </ResponsiveContainer>
        );
    };

    const renderFlexibilityChart = (caseId) => {
        if (caseId === 'case0') return <div>No flexibility data for Base Case</div>;
        
        const data = caseId === 'comparison' 
            ? caseStudyData.comparison?.flexibility || []
            : caseStudyData[caseId]?.flexibility || [];
        
        if (!data.length) return <div>No flexibility data available</div>;

        return (
            <ResponsiveContainer width="100%" height={400}>
                <BarChart data={data} margin={{ top: 30, right: 40, left: 60, bottom: 100 }}>
                    <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                    <XAxis
                        dataKey="hour"
                        label={{ value: 'Hour of Day', position: 'insideBottom', offset: -20, style: { fontSize: '14px', fontWeight: 600, textAnchor: 'middle', fill: '#1e293b' } }}
                        tick={{ fontSize: 12, fill: '#64748b' }}
                        stroke="#94a3b8"
                    />
                    <YAxis
                        label={{ value: 'Flexibility (MW)', angle: -90, position: 'insideLeft', style: { fontSize: '14px', fontWeight: 600, textAnchor: 'middle', fill: '#1e293b' } }}
                        tick={{ fontSize: 12, fill: '#64748b' }}
                        stroke="#94a3b8"
                    />
                    <Tooltip contentStyle={{ backgroundColor: 'rgba(255, 255, 255, 0.98)', border: '1px solid #e2e8f0', borderRadius: '12px' }} />
                    <Legend />
                    {caseId === 'comparison' ? (
                        <>
                            <Bar dataKey="Case Study 1" fill="rgba(59, 130, 246, 0.8)" />
                            <Bar dataKey="Case Study 2" fill="rgba(16, 185, 129, 0.8)" />
                            <Bar dataKey="Case Study 3" fill="rgba(245, 158, 11, 0.8)" />
                        </>
                    ) : (
                        <Bar dataKey="flexibility" fill="rgba(236, 72, 153, 0.8)" name="Data Center Flexibility" />
                    )}
                </BarChart>
            </ResponsiveContainer>
        );
    };

    const renderTotalSystemCostsChart = () => {
        const data = caseStudyData.comparison?.totalSystemCosts || [];
        if (!data.length) return <div>No total system costs data available</div>;

        // Group data by case study
        const groupedData = [];
        for (let hour = 1; hour <= 24; hour++) {
            const hourData = { hour };
            data.filter(d => d.hour === hour).forEach(d => {
                hourData[d.caseStudy] = d.cost;
            });
            groupedData.push(hourData);
        }

        return (
            <ResponsiveContainer width="100%" height={400}>
                <AreaChart data={groupedData} margin={{ top: 30, right: 40, left: 60, bottom: 100 }}>
                    <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                    <XAxis
                        dataKey="hour"
                        label={{ value: 'Hour of Day', position: 'insideBottom', offset: -20, style: { fontSize: '14px', fontWeight: 600, textAnchor: 'middle', fill: '#1e293b' } }}
                        tick={{ fontSize: 12, fill: '#64748b' }}
                        stroke="#94a3b8"
                    />
                    <YAxis
                        label={{ value: 'Total System Cost (M$)', angle: -90, position: 'insideLeft', style: { fontSize: '14px', fontWeight: 600, textAnchor: 'middle', fill: '#1e293b' } }}
                        tick={{ fontSize: 12, fill: '#64748b' }}
                        stroke="#94a3b8"
                    />
                    <Tooltip 
                        contentStyle={{ backgroundColor: 'rgba(255, 255, 255, 0.98)', border: '1px solid #e2e8f0', borderRadius: '12px' }}
                        formatter={(value) => [`$${value.toFixed(2)}M`, 'Cost']}
                    />
                    <Legend />
                    <Area type="monotone" dataKey="Case_0" stroke="#64748B" fill="rgba(100, 116, 139, 0.6)" name="Case 0" />
                    <Area type="monotone" dataKey="Case_1" stroke="#3B82F6" fill="rgba(59, 130, 246, 0.6)" name="Case 1" />
                    <Area type="monotone" dataKey="Case_2" stroke="#10B981" fill="rgba(16, 185, 129, 0.6)" name="Case 2" />
                    <Area type="monotone" dataKey="Case_3" stroke="#F59E0B" fill="rgba(245, 158, 11, 0.6)" name="Case 3" />
                </AreaChart>
            </ResponsiveContainer>
        );
    };

    // Pie chart for generation mix
    const renderGenerationPieChart = (caseId) => {
        const data = caseStudyData[caseId]?.generation || [];
        if (!data.length) return <div>No generation data available</div>;

        // Calculate total generation by fuel type
        const fuelTotals = data.reduce((totals, hourData) => {
            totals.naturalGas = (totals.naturalGas || 0) + (hourData.naturalGas || 0);
            totals.hydro = (totals.hydro || 0) + (hourData.hydro || 0);
            totals.wind = (totals.wind || 0) + (hourData.wind || 0);
            totals.solar = (totals.solar || 0) + (hourData.solar || 0);
            totals.nuclear = (totals.nuclear || 0) + (hourData.nuclear || 0);
            totals.coal = (totals.coal || 0) + (hourData.coal || 0);
            totals.geothermal = (totals.geothermal || 0) + (hourData.geothermal || 0);
            totals.biomass = (totals.biomass || 0) + (hourData.biomass || 0);
            totals.battery = (totals.battery || 0) + (hourData.battery || 0);
            return totals;
        }, {});

        const pieData = Object.entries(fuelTotals)
            .filter(([_, value]) => value > 0)
            .map(([fuel, value]) => ({
                name: fuel.charAt(0).toUpperCase() + fuel.slice(1).replace(/([A-Z])/g, ' $1'),
                value: Math.round(value),
                percentage: Math.round((value / Object.values(fuelTotals).reduce((a, b) => a + b, 0)) * 100)
            }))
            .sort((a, b) => b.value - a.value);

        const COLORS = {
            'Natural Gas': '#8b5cf6', 'Hydro': '#3b82f6', 'Wind': '#10b981',
            'Solar': '#f59e0b', 'Nuclear': '#ef4444', 'Coal': '#6b7280',
            'Geothermal': '#dc2626', 'Biomass': '#059669', 'Battery': '#ec4899'
        };

        return (
            <ResponsiveContainer width="100%" height={350}>
                <PieChart>
                    <Pie
                        data={pieData}
                        cx="50%"
                        cy="50%"
                        outerRadius={100}
                        innerRadius={50}
                        paddingAngle={1}
                        dataKey="value"
                    >
                        {pieData.map((entry, index) => (
                            <Cell key={`cell-${index}`} fill={COLORS[entry.name] || '#8884d8'} />
                        ))}
                    </Pie>
                    <Tooltip 
                        contentStyle={{ backgroundColor: 'rgba(255, 255, 255, 0.98)', border: '1px solid #e2e8f0', borderRadius: '12px' }}
                        formatter={(value) => [`${value.toLocaleString()} MW`, 'Generation']}
                    />
                    <Legend verticalAlign="bottom" height={60} />
                </PieChart>
            </ResponsiveContainer>
        );
    };

    // Chart container component
    const ChartContainer = ({ title, icon, children }) => (
        <div style={{
            background: 'linear-gradient(135deg, #ffffff 0%, #f8fafc 100%)',
            borderRadius: '20px',
            padding: '32px',
            boxShadow: '0 8px 32px rgba(0, 0, 0, 0.08)',
            border: '1px solid rgba(226, 232, 240, 0.8)',
            backdropFilter: 'blur(10px)',
            marginBottom: '32px'
        }}>
            <div style={{ display: 'flex', alignItems: 'center', gap: '12px', marginBottom: '24px' }}>
                <span style={{ fontSize: '24px' }}>{icon}</span>
                <h3 style={{ 
                    margin: 0, 
                    fontSize: '1.5rem', 
                    fontWeight: 700,
                    color: '#1e293b'
                }}>
                    {title}
                </h3>
            </div>
            {children}
        </div>
    );

    // Main render function
    if (loading) {
        return (
            <div style={{ 
                display: 'flex', 
                justifyContent: 'center', 
                alignItems: 'center', 
                height: '100vh',
                background: 'linear-gradient(135deg, #f8fafc 0%, #e2e8f0 100%)'
            }}>
                <div style={{
                    background: 'rgba(255, 255, 255, 0.9)',
                    padding: '40px',
                    borderRadius: '20px',
                    boxShadow: '0 8px 32px rgba(0, 0, 0, 0.1)',
                    textAlign: 'center'
                }}>
                    <div style={{ fontSize: '48px', marginBottom: '16px' }}>⚡</div>
                    <h2 style={{ margin: '0 0 8px 0', color: '#1e293b' }}>Loading WECC Analytics</h2>
                    <p style={{ margin: 0, color: '#64748b' }}>Preparing comprehensive data analysis...</p>
                </div>
            </div>
        );
    }

    if (error) {
        return (
            <div style={{ 
                display: 'flex', 
                justifyContent: 'center', 
                alignItems: 'center', 
                height: '100vh',
                background: 'linear-gradient(135deg, #f8fafc 0%, #e2e8f0 100%)'
            }}>
                <div style={{
                    background: 'rgba(255, 255, 255, 0.9)',
                    padding: '40px',
                    borderRadius: '20px',
                    boxShadow: '0 8px 32px rgba(0, 0, 0, 0.1)',
                    textAlign: 'center'
                }}>
                    <div style={{ fontSize: '48px', marginBottom: '16px' }}>❌</div>
                    <h2 style={{ margin: '0 0 8px 0', color: '#dc2626' }}>Error Loading Data</h2>
                    <p style={{ margin: 0, color: '#64748b' }}>{error}</p>
                </div>
            </div>
        );
    }

    return (
        <div style={{
            background: 'linear-gradient(135deg, #f8fafc 0%, #e2e8f0 100%)',
            minHeight: '100vh',
            overflow: 'auto',
            fontFamily: 'Inter, system-ui, sans-serif'
        }}>
            {/* Enhanced Header */}
            <div style={{
                background: 'linear-gradient(135deg, #1e3a8a 0%, #3b82f6 100%)',
                padding: '24px 40px',
                boxShadow: '0 4px 20px rgba(0, 0, 0, 0.1)'
            }}>
                <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between' }}>
                    <button
                        onClick={() => navigate('/')}
                        style={{
                            padding: '12px 20px',
                            background: 'rgba(255,255,255,0.1)',
                            color: 'white',
                            border: '1px solid rgba(255,255,255,0.2)',
                            borderRadius: '12px',
                            cursor: 'pointer',
                            fontSize: '14px',
                            fontWeight: 500,
                            transition: 'all 0.2s ease',
                            display: 'flex',
                            alignItems: 'center',
                            gap: '8px'
                        }}
                    >
                        <ArrowBackIcon style={{ fontSize: 16 }} />
                        Back to Westmap
                    </button>
                    
                    <div style={{ textAlign: 'center', flex: 1 }}>
                        <h1 style={{ 
                            margin: 0, 
                            fontSize: '2.2rem', 
                            fontWeight: 800,
                            color: 'white',
                            textShadow: '0 2px 4px rgba(0,0,0,0.1)'
                        }}>
                            🏛️ WECC Analytics Portal
                        </h1>
                        <p style={{ 
                            margin: '8px 0 0 0', 
                            color: 'rgba(255, 255, 255, 0.9)', 
                            fontSize: '1.1rem',
                            fontWeight: 500
                        }}>
                            Comprehensive Data Center Impact Analysis
                        </p>
                    </div>
                    
                    {/* WECC Region Selector */}
                    <select
                        value={selectedWeccRegion}
                        onChange={(e) => setSelectedWeccRegion(e.target.value)}
                        style={{
                            padding: '12px 16px',
                            borderRadius: '12px',
                            border: '1px solid rgba(255,255,255,0.2)',
                            background: 'rgba(255,255,255,0.1)',
                            color: 'white',
                            fontSize: '14px',
                            fontWeight: 500,
                            cursor: 'pointer',
                            minWidth: '200px'
                        }}
                    >
                        {BALANCING_AUTHORITIES.map(ba => (
                            <option key={ba.code} value={ba.code} style={{ color: '#1e293b' }}>
                                {ba.code} - {ba.name}
                            </option>
                        ))}
                    </select>
                </div>
            </div>

            {/* Main Content */}
            <div style={{ padding: '40px' }}>
                {/* Case Study Selector */}
                <div style={{
                    display: 'grid',
                    gridTemplateColumns: 'repeat(auto-fit, minmax(300px, 1fr))',
                    gap: '20px',
                    marginBottom: '40px'
                }}>
                    {CASE_STUDIES.map(caseStudy => (
                        <div
                            key={caseStudy.id}
                            onClick={() => setSelectedCaseStudy(caseStudy.id)}
                            style={{
                                background: selectedCaseStudy === caseStudy.id 
                                    ? `linear-gradient(135deg, ${caseStudy.color}20 0%, ${caseStudy.color}10 100%)`
                                    : 'linear-gradient(135deg, #ffffff 0%, #f8fafc 100%)',
                                border: selectedCaseStudy === caseStudy.id 
                                    ? `2px solid ${caseStudy.color}`
                                    : '2px solid transparent',
                                borderRadius: '20px',
                                padding: '24px',
                                cursor: 'pointer',
                                transition: 'all 0.3s ease',
                                boxShadow: selectedCaseStudy === caseStudy.id 
                                    ? '0 12px 32px rgba(0, 0, 0, 0.15)'
                                    : '0 8px 32px rgba(0, 0, 0, 0.06)',
                                transform: selectedCaseStudy === caseStudy.id ? 'translateY(-4px)' : 'translateY(0)'
                            }}
                        >
                            <div style={{ display: 'flex', alignItems: 'center', gap: '12px', marginBottom: '12px' }}>
                                <span style={{ fontSize: '24px' }}>{caseStudy.icon}</span>
                                <div>
                                    <h3 style={{ 
                                        margin: 0, 
                                        fontSize: '1.25rem', 
                                        fontWeight: 700,
                                        color: selectedCaseStudy === caseStudy.id ? caseStudy.color : '#1e293b'
                                    }}>
                                        {caseStudy.name}
                                    </h3>
                                    <p style={{ 
                                        margin: '2px 0 0 0', 
                                        fontSize: '0.9rem', 
                                        fontWeight: 600,
                                        color: selectedCaseStudy === caseStudy.id ? caseStudy.color : '#64748b'
                                    }}>
                                        {caseStudy.subtitle}
                                    </p>
                                </div>
                            </div>
                            <p style={{ 
                                margin: 0, 
                                color: '#64748b', 
                                fontSize: '0.95rem',
                                lineHeight: '1.5'
                            }}>
                                {caseStudy.description}
                            </p>
                        </div>
                    ))}
                </div>

                {/* Hero Section: Map and Info */}
                <div style={{
                    display: 'grid',
                    gridTemplateColumns: '2fr 1fr',
                    gap: '32px',
                    marginBottom: '60px'
                }}>
                    {/* Enhanced Map */}
                    <div style={{
                        background: 'linear-gradient(135deg, #ffffff 0%, #f8fafc 100%)',
                        borderRadius: '20px',
                        padding: '24px',
                        boxShadow: '0 8px 32px rgba(0, 0, 0, 0.08)',
                        border: '1px solid rgba(226, 232, 240, 0.8)',
                        backdropFilter: 'blur(10px)',
                        height: '500px'
                    }}>
                        <h3 style={{ 
                            margin: '0 0 20px 0', 
                            fontSize: '1.5rem', 
                            fontWeight: 700,
                            color: '#1e293b'
                        }}>
                            {balancingAuthority?.BA_Abrev || selectedWeccRegion} Territory
                        </h3>
                        <div style={{ height: '440px' }}>
                            {renderMap()}
                        </div>
                    </div>

                    {/* WECC Generation Pie Chart - Aligned with Map */}
                    <div style={{
                        background: 'linear-gradient(135deg, rgba(16, 185, 129, 0.08) 0%, rgba(52, 211, 153, 0.04) 100%)',
                        borderRadius: '20px',
                        padding: '24px',
                        boxShadow: '0 8px 32px rgba(0, 0, 0, 0.08)',
                        border: '1px solid rgba(16, 185, 129, 0.1)',
                        backdropFilter: 'blur(10px)',
                        height: '500px'
                    }}>
                        <div style={{ display: 'flex', alignItems: 'center', gap: '12px', marginBottom: '20px' }}>
                            <div style={{
                                background: 'linear-gradient(135deg, #10b981 0%, #34d399 100%)',
                                borderRadius: '12px',
                                padding: '12px',
                                boxShadow: '0 4px 12px rgba(16, 185, 129, 0.3)'
                            }}>
                                <span style={{ fontSize: '24px', filter: 'drop-shadow(0 2px 4px rgba(0,0,0,0.1))' }}>⚡</span>
                            </div>
                            <div>
                                <h3 style={{ 
                                    margin: 0, 
                                    fontSize: '1.5rem', 
                                    fontWeight: 700,
                                    color: '#1e293b'
                                }}>
                                    Generation Portfolio
                                </h3>
                                <p style={{
                                    margin: '4px 0 0 0',
                                    fontSize: '1rem',
                                    color: '#64748b',
                                    fontWeight: 500
                                }}>
                                    {balancingAuthority?.BA_Abrev || selectedWeccRegion} - {selectedCaseStudy === 'comparison' ? 'Comparison' : CASE_STUDIES.find(cs => cs.id === selectedCaseStudy)?.subtitle}
                                </p>
                            </div>
                        </div>
                        <div style={{ 
                            height: '420px',
                            background: 'rgba(255, 255, 255, 0.4)',
                            borderRadius: '12px',
                            padding: '16px'
                        }}>
                            {renderGenerationPieChart(selectedCaseStudy)}
                        </div>
                    </div>
                </div>

                {/* Charts Section */}
                <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '32px' }}>
                    {/* Left Column - Main Charts */}
                    <div>
                        <ChartContainer title="Balancing Authority Price" icon="💰">
                            {renderPriceChart(selectedCaseStudy)}
                        </ChartContainer>

                        <ChartContainer title="Balancing Authority Demand" icon="📊">
                            {renderDemandChart(selectedCaseStudy)}
                        </ChartContainer>

                        <ChartContainer title="Data Center Demand" icon="🏢">
                            {renderDataCenterDemandChart(selectedCaseStudy)}
                        </ChartContainer>

                        <ChartContainer title="Power Generation by Technology" icon="⚡">
                            {renderGenerationChart(selectedCaseStudy)}
                        </ChartContainer>
                    </div>

                    {/* Right Column - Additional Charts */}
                    <div>
                        <ChartContainer title="Generation Portfolio Mix" icon="🥧">
                            {renderGenerationPieChart(selectedCaseStudy)}
                        </ChartContainer>

                        <ChartContainer title="Total Hourly Operation Costs" icon="💸">
                            {renderCostsChart(selectedCaseStudy)}
                        </ChartContainer>

                        <ChartContainer title="Net Electricity Interchange" icon="🔄">
                            {renderInterchangeChart(selectedCaseStudy)}
                        </ChartContainer>

                        {selectedCaseStudy !== 'case0' && (
                            <ChartContainer title="Data Center Energy Flexibility" icon="🔧">
                                {renderFlexibilityChart(selectedCaseStudy)}
                            </ChartContainer>
                        )}

                        {selectedCaseStudy === 'comparison' && (
                            <ChartContainer title="Total System Operation Cost" icon="🏦">
                                {renderTotalSystemCostsChart()}
                            </ChartContainer>
                        )}
                    </div>
                </div>
            </div>
        </div>
    );
};

export default AminDetailPage;
