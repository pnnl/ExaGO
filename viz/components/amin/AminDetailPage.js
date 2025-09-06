import React, { useState, useEffect } from 'react';
import { useParams, useNavigate, useLocation } from 'react-router-dom';
import { DeckGL } from '@deck.gl/react';
import { GeoJsonLayer, ScatterplotLayer } from '@deck.gl/layers';
import { MapView } from '@deck.gl/core';
import { StaticMap } from 'react-map-gl';
import ArrowBackIcon from '@mui/icons-material/ArrowBack';
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer, AreaChart, Area, BarChart, Bar, ComposedChart } from 'recharts';

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

// Case study configurations
const CASE_STUDIES = [
    {
        id: 'case1',
        name: 'Case Study 1',
        subtitle: 'Base Case',
        description: 'Current operational baseline scenario without data center flexibility',
        color: '#3B82F6',
        icon: '📊',
        dataPath: 'Case study_1/CASE STUDY 1'
    },
    {
        id: 'case2',
        name: 'Case Study 2',
        subtitle: 'Data Center Temporal Flexibility',
        description: 'Time-shifted data center operations with flexible scheduling',
        color: '#10B981',
        icon: '⏰',
        dataPath: 'Case study_2'
    },
    {
        id: 'case3',
        name: 'Case Study 3',
        subtitle: 'Data Center Spatial Flexibility',
        description: 'Geographic load distribution optimization across regions',
        color: '#F59E0B',
        icon: '🌐',
        dataPath: 'Case study_3'
    },
    {
        id: 'case4',
        name: 'Case Study 4',
        subtitle: 'Comparison',
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
    { code: 'WACM', name: 'Western Area Power Administration - Colorado Missouri Region' },
    { code: 'WALC', name: 'Western Area Power Administration - Lower Colorado Region' },
    { code: 'WAUW', name: 'Western Area Power Administration - Upper Great Plains West' }
];

// Legacy array for backward compatibility
const BALANCING_AUTHORITIES_CODES = BALANCING_AUTHORITIES.map(ba => ba.code);

// Area mapping for WECC regions to area numbers (based on WECC_BA_Area_Mapping.csv)
const WECC_AREA_MAPPING = {
    'AESO': [1],
    'AVA': [2],
    'AZPS': [3],
    'BANC': [4],
    'BCHA': [5],
    'BPAT': [6],
    'CENACE': [7],
    'CHPD': [8],
    'CISO': [9],
    'DOPD': [10],
    'EPE': [11],
    'GCPD': [12],
    'IID': [13],
    'IPCO': [14],
    'LDWP': [15],
    'NEVP': [16],
    'NWMT': [17],
    'PACE': [18],
    'PACW': [19],
    'PGE': [20],
    'PNM': [21],
    'PSCO': [22],
    'PSEI': [23],
    'SCL': [24],
    'SRP': [25],
    'TEPC': [26],
    'TIDC': [27],
    'TPWR': [28],
    'WACM': [29],
    'WALC': [30],
    'WAUW': [31]
};

// Enhanced balancing authorities with area numbers
const BALANCING_AUTHORITIES_WITH_AREAS = BALANCING_AUTHORITIES.map(ba => ({
    ...ba,
    areaNumbers: WECC_AREA_MAPPING[ba.code] || [],
    primaryArea: WECC_AREA_MAPPING[ba.code]?.[0] || null
}));

const AminDetailPage = () => {
    const { fid } = useParams();
    const navigate = useNavigate();
    const location = useLocation();
    
    // Get WECC region from URL params or default to NEVP
    const urlParams = new URLSearchParams(location.search);
    const initialWeccRegion = urlParams.get('region') || 'NEVP';
    
    const [selectedWeccRegion, setSelectedWeccRegion] = useState(initialWeccRegion);
    const [balancingAuthority, setBalancingAuthority] = useState(null);
    const [selectedCaseStudy, setSelectedCaseStudy] = useState('case1');
    const [caseStudyData, setCaseStudyData] = useState({});
    const [capacityData, setCapacityData] = useState({});
    const [plantLocations, setPlantLocations] = useState([]);
    const [dataCenterLocations, setDataCenterLocations] = useState([]);
    const [selectedHour, setSelectedHour] = useState(12); // Default to noon for zonal price display
    const [viewState, setViewState] = useState({
        longitude: -116.5,
        latitude: 37.5,
        zoom: 6,
        minZoom: 3,
        maxZoom: 12,
        pitch: 0,
        bearing: 0
    });
    const [loading, setLoading] = useState(true);
    const [error, setError] = useState(null);

    // Load Google Fonts Inter
    useEffect(() => {
        const link = document.createElement('link');
        link.href = 'https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700;800;900&display=swap';
        link.rel = 'stylesheet';
        document.head.appendChild(link);

        return () => {
            if (document.head.contains(link)) {
                document.head.removeChild(link);
            }
        };
    }, []);

    // Set Inter font family on body
    useEffect(() => {
        const originalFontFamily = document.body.style.fontFamily;
        document.body.style.fontFamily = '"Inter", system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif';
        
        return () => {
            document.body.style.fontFamily = originalFontFamily;
        };
    }, []);

    // Add CSS for enhanced styling
    useEffect(() => {
        const style = document.createElement('style');
        style.textContent = `
            .recharts-legend-wrapper {
                bottom: 28px !important;
            }
            
            .case-study-selector {
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
                gap: 16px;
                margin-bottom: 32px;
            }
            
            .case-study-card {
                background: linear-gradient(135deg, #ffffff 0%, #f8fafc 100%);
                border: 2px solid transparent;
                border-radius: 16px;
                padding: 24px;
                cursor: pointer;
                transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
                position: relative;
                overflow: hidden;
            }
            
            .case-study-card::before {
                content: '';
                position: absolute;
                top: 0;
                left: 0;
                right: 0;
                bottom: 0;
                background: linear-gradient(135deg, rgba(59, 130, 246, 0.1) 0%, rgba(147, 197, 253, 0.05) 100%);
                opacity: 0;
                transition: opacity 0.3s ease;
                z-index: 0;
            }
            
            .case-study-card:hover::before {
                opacity: 1;
            }
            
            .case-study-card:hover {
                transform: translateY(-2px);
                box-shadow: 0 20px 40px rgba(0, 0, 0, 0.1);
                border-color: rgba(59, 130, 246, 0.2);
            }
            
            .case-study-card.active {
                border-color: var(--case-color);
                background: linear-gradient(135deg, #ffffff 0%, #fefefe 100%);
                box-shadow: 0 12px 32px rgba(0, 0, 0, 0.15);
                transform: translateY(-1px);
            }
            
            .case-study-card.active::before {
                opacity: 0.8;
                background: linear-gradient(135deg, var(--case-color-light) 0%, var(--case-color-lighter) 100%);
            }
            
            .case-study-content {
                position: relative;
                z-index: 1;
            }
            
            .chart-container {
                background: linear-gradient(135deg, #ffffff 0%, #fafbfc 100%);
                border-radius: 20px;
                padding: 32px;
                box-shadow: 0 8px 32px rgba(0, 0, 0, 0.08);
                border: 1px solid rgba(226, 232, 240, 0.8);
                backdrop-filter: blur(10px);
                transition: all 0.3s ease;
            }
            
            .chart-container:hover {
                box-shadow: 0 16px 48px rgba(0, 0, 0, 0.12);
                transform: translateY(-2px);
            }
            
            .comparison-grid {
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(400px, 1fr));
                gap: 24px;
                margin-top: 24px;
            }
            
            .hour-selector {
                display: flex;
                align-items: center;
                gap: 12px;
                margin-bottom: 20px;
                padding: 16px;
                background: rgba(59, 130, 246, 0.05);
                border-radius: 12px;
            }
            
            .hour-slider {
                flex: 1;
                margin: 0 16px;
            }
            
            .info-card {
                background: linear-gradient(135deg, #ffffff 0%, #f8fafc 100%);
                border-radius: 20px;
                padding: 32px;
                box-shadow: 0 8px 32px rgba(0, 0, 0, 0.08);
                border: 1px solid rgba(226, 232, 240, 0.8);
                backdrop-filter: blur(10px);
            }
            
            .map-container {
                background: linear-gradient(135deg, #ffffff 0%, #f8fafc 100%);
                border-radius: 20px;
                padding: 24px;
                box-shadow: 0 8px 32px rgba(0, 0, 0, 0.08);
                border: 1px solid rgba(226, 232, 240, 0.8);
                backdrop-filter: blur(10px);
            }
            
            .gradient-text {
                background: linear-gradient(135deg, #1e293b 0%, #334155 100%);
                -webkit-background-clip: text;
                -webkit-text-fill-color: transparent;
                background-clip: text;
            }
            
            .metric-grid {
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
                gap: 16px;
                margin-top: 24px;
            }
            
            .metric-item {
                background: linear-gradient(135deg, rgba(59, 130, 246, 0.05) 0%, rgba(147, 197, 253, 0.02) 100%);
                border: 1px solid rgba(59, 130, 246, 0.1);
                border-radius: 12px;
                padding: 16px;
                transition: all 0.2s ease;
            }
            
            .metric-item:hover {
                background: linear-gradient(135deg, rgba(59, 130, 246, 0.08) 0%, rgba(147, 197, 253, 0.04) 100%);
                border-color: rgba(59, 130, 246, 0.2);
            }
            
            .clickable-generator {
                cursor: pointer;
                transition: all 0.2s ease;
            }
            
            .clickable-generator:hover {
                transform: scale(1.2);
            }
            
            .data-center-marker {
                cursor: pointer;
                transition: all 0.2s ease;
            }
            
            .data-center-marker:hover {
                transform: scale(1.1);
            }
        `;
        document.head.appendChild(style);

        return () => {
            if (document.head.contains(style)) {
                document.head.removeChild(style);
            }
        };
    }, []);

    // Override body overflow to allow scrolling on this page
    useEffect(() => {
        const originalOverflow = document.body.style.overflow;
        document.body.style.overflow = 'auto';

        return () => {
            document.body.style.overflow = originalOverflow;
        };
    }, []);

    // Handle browser back button navigation
    useEffect(() => {
        const handlePopState = (event) => {
            const currentPath = window.location.pathname;
            if ((currentPath.includes('/amin/') || currentPath.includes('/manish/')) && !balancingAuthority) {
                let parentRoute = '/amin';

                if (currentPath.includes('/manish/')) {
                    parentRoute = '/manish';
                } else if (currentPath.includes('/amin/')) {
                    parentRoute = '/amin';
                }

                navigate(parentRoute, { replace: true });
            }
        };

        window.addEventListener('popstate', handlePopState);
        return () => {
            window.removeEventListener('popstate', handlePopState);
        };
    }, [navigate, balancingAuthority]);

    // Load all data
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

                // Load area mapping CSV
                const areaMappingResponse = await fetch('/amin_data/WECC_BA_Area_Mapping.csv');
                if (!areaMappingResponse.ok) {
                    throw new Error(`HTTP error loading area mapping CSV! status: ${areaMappingResponse.status}`);
                }
                const areaMappingText = await areaMappingResponse.text();

                // Load capacity data
                const capacityResponse = await fetch('/amin_data/manish_amin_modified_data/WECC_BA_CAPACITY.csv');
                if (capacityResponse.ok) {
                    const capacityText = await capacityResponse.text();
                    const capacityLines = capacityText.split('\n');
                    const capacityMap = {};
                    
                    for (let i = 1; i < capacityLines.length; i++) {
                        const line = capacityLines[i].trim();
                        if (line) {
                            const values = line.split(',');
                            const ba = values[0];
                            capacityMap[ba] = {
                                hydro: parseFloat(values[1]) || 0,
                                nuclear: parseFloat(values[2]) || 0,
                                coal: parseFloat(values[3]) || 0,
                                naturalGas: parseFloat(values[4].replace(/"/g, '')) || 0,
                                geothermal: parseFloat(values[5]) || 0,
                                biomass: parseFloat(values[6]) || 0,
                                wind: parseFloat(values[7]) || 0,
                                pv: parseFloat(values[8]) || 0,
                                batteryStorage: parseFloat(values[9]) || 0
                            };
                        }
                    }
                    setCapacityData(capacityMap);
                }

                // Load plant locations
                const plantResponse = await fetch('/amin_data/manish_amin_modified_data/Locations of Power Plants/Western_Power_plants_Locations-USA.csv');
                if (plantResponse.ok) {
                    const plantText = await plantResponse.text();
                    const plantLines = plantText.split('\n');
                    const plants = [];
                    
                    for (let i = 1; i < plantLines.length; i++) {
                        const line = plantLines[i].trim();
                        if (line) {
                            const values = line.split(',');
                            plants.push({
                                plantCode: values[0],
                                plantName: values[1],
                                latitude: parseFloat(values[2]),
                                longitude: parseFloat(values[3]),
                                state: values[4],
                                county: values[5],
                                balancingAuthority: values[6],
                                primaryType: values[7],
                                totalCapacity: parseFloat(values[8]) || 0
                            });
                        }
                    }
                    setPlantLocations(plants);
                }

                // Parse main CSV for shape data
                const csvLines = csvText.split('\n');
                const csvData = {};

                for (let i = 1; i < csvLines.length; i++) {
                    const line = csvLines[i].trim();
                    if (line) {
                        const values = line.split(',');
                        const fid = parseInt(values[0]);
                        csvData[fid] = {
                            FID: fid,
                            BA_Abrev: values[1],
                            BA_Name: values[2].replace(/"/g, ''),
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
                        const fidValue = parseInt(values[0]);
                        const areaNumbersString = values[3];

                        let areaNumbers;
                        if (areaNumbersString && areaNumbersString.includes('|')) {
                            areaNumbers = areaNumbersString.split('|').map(num => parseInt(num.trim()));
                        } else {
                            areaNumbers = areaNumbersString ? [parseInt(areaNumbersString)] : [];
                        }

                        areaMappingData[fidValue] = {
                            FID: fidValue,
                            BA_Abrev: values[1],
                            BA_Name: values[2].replace(/"/g, ''),
                            Area_Numbers: areaNumbers
                        };
                    }
                }

                // Find the specific balancing authority
                const targetFid = parseInt(fid);
                const feature = geojsonData.features.find(f => f.properties.FID === targetFid);

                if (!feature) {
                    throw new Error(`Balancing authority with FID ${fid} not found`);
                }

                // Merge data from all sources
                const baData = {
                    ...feature.properties,
                    ...csvData[targetFid],
                    ...areaMappingData[targetFid]
                };

                // Update selected WECC region based on the loaded data
                if (baData.BA_Abrev && baData.BA_Abrev !== selectedWeccRegion) {
                    setSelectedWeccRegion(baData.BA_Abrev);
                }

                // Load case study data for this BA
                await loadCaseStudyData(baData.BA_Abrev);

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

                if (!validCoordinates) {
                    console.warn(`Invalid coordinates for FID ${fid}, using default view`);
                    setViewState(prev => ({
                        ...prev,
                        longitude: -116.5,
                        latitude: 37.5,
                        zoom: 5
                    }));
                } else {
                    const centerLng = (minLng + maxLng) / 2;
                    const centerLat = (minLat + maxLat) / 2;

                    if (isNaN(centerLng) || isNaN(centerLat)) {
                        console.warn(`Invalid center coordinates for FID ${fid}, using default view`);
                        setViewState(prev => ({
                            ...prev,
                            longitude: -116.5,
                            latitude: 37.5,
                            zoom: 5
                        }));
                    } else {
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
                }

                setBalancingAuthority({
                    ...baData,
                    feature: feature,
                    bounds: validCoordinates ? { minLng, maxLng, minLat, maxLat } : null
                });
                setLoading(false);
            } catch (err) {
                console.error('Error loading data:', err);
                setError(err.message);
                setLoading(false);
            }
        };

        if (fid) {
            loadData();
        } else if (selectedWeccRegion) {
            // Load data directly for WECC region without fid
            loadDataForWeccRegion();
        }
    }, [fid, selectedWeccRegion]);

    // Load data for WECC region without requiring fid
    const loadDataForWeccRegion = async () => {
        setLoading(true);
        setError(null);
        
        try {
            // Create a simplified balancing authority object for the selected region
            const selectedBA = BALANCING_AUTHORITIES.find(ba => ba.code === selectedWeccRegion);
            if (!selectedBA) {
                throw new Error(`WECC region ${selectedWeccRegion} not found`);
            }

            setBalancingAuthority({
                BA_Abrev: selectedBA.code,
                BA_Name: selectedBA.name,
                FID: null, // No FID needed for direct region access
                feature: null, // No geographic feature
                bounds: null
            });

            // Set default view for WECC region
            setViewState(prev => ({
                ...prev,
                longitude: -116.5,
                latitude: 37.5,
                zoom: 6
            }));

            // Load case study data for the selected region
            await loadCaseStudyData(selectedBA.code);
            
            setLoading(false);
        } catch (err) {
            console.error('Error loading WECC region data:', err);
            setError(err.message);
            setLoading(false);
        }
    };

    // Handle WECC region change
    const handleWeccRegionChange = (newRegion) => {
        setSelectedWeccRegion(newRegion);
        
        // Get the primary area number for this region
        const areaNumber = WECC_AREA_MAPPING[newRegion]?.[0];
        
        if (areaNumber) {
            // Navigate to the area-specific URL (e.g., /manish/9, /manish/10)
            const newUrl = `/manish/${areaNumber}`;
            navigate(newUrl, { replace: true });
        } else {
            // Fallback to region parameter if no area number found
            const newUrl = `/amin?region=${newRegion}`;
            navigate(newUrl, { replace: true });
        }
    };

    // Load case study data
    const loadCaseStudyData = async (baAbbrev) => {
        const caseData = {};
        
        for (const caseStudy of CASE_STUDIES) {
            if (caseStudy.id === 'case4') continue; // Skip comparison for individual loading
            
            try {
                // Load demand data
                const demandPath = `/amin_data/manish_amin_modified_data/${caseStudy.dataPath}/Balancing Authority Demand (MW).csv`;
                const demandResponse = await fetch(demandPath);
                if (demandResponse.ok) {
                    const demandText = await demandResponse.text();
                    caseData[caseStudy.id] = { ...caseData[caseStudy.id], demand: parseDemandData(demandText, baAbbrev) };
                }

                // Load generation data
                const generationPath = `/amin_data/manish_amin_modified_data/${caseStudy.dataPath}/Balancing Authority Power Generation/${baAbbrev}_generation_by_fuel.csv`;
                const generationResponse = await fetch(generationPath);
                if (generationResponse.ok) {
                    const generationText = await generationResponse.text();
                    caseData[caseStudy.id] = { ...caseData[caseStudy.id], generation: parseGenerationData(generationText) };
                }

                // Load LMP data
                const lmpPath = `/amin_data/manish_amin_modified_data/${caseStudy.dataPath}/Zonal Price/lmp_by_ba_hour.csv`;
                const lmpResponse = await fetch(lmpPath);
                if (lmpResponse.ok) {
                    const lmpText = await lmpResponse.text();
                    caseData[caseStudy.id] = { ...caseData[caseStudy.id], lmp: parseLMPData(lmpText, baAbbrev) };
                }

                // Load operational costs
                const costsPath = `/amin_data/manish_amin_modified_data/${caseStudy.dataPath}/Balancing Authority Hourly Operation Costs/${baAbbrev}_hourly_operation_costs.csv`;
                const costsResponse = await fetch(costsPath);
                if (costsResponse.ok) {
                    const costsText = await costsResponse.text();
                    caseData[caseStudy.id] = { ...caseData[caseStudy.id], costs: parseCostsData(costsText) };
                }

                // Load data center demand (for case studies 2-4)
                if (caseStudy.id !== 'case1') {
                    const dcDemandPath = `/amin_data/manish_amin_modified_data/${caseStudy.dataPath}/Data Center Demand (MW).csv`;
                    const dcDemandResponse = await fetch(dcDemandPath);
                    if (dcDemandResponse.ok) {
                        const dcDemandText = await dcDemandResponse.text();
                        caseData[caseStudy.id] = { ...caseData[caseStudy.id], dataCenterDemand: parseDataCenterDemand(dcDemandText, baAbbrev) };
                    }

                    // Load data center flexibility (for case studies 2-4)
                    const dcFlexPath = `/amin_data/manish_amin_modified_data/${caseStudy.dataPath}/Data Center Flexibility/Data Center Energy Flexibility.csv`;
                    const dcFlexResponse = await fetch(dcFlexPath);
                    if (dcFlexResponse.ok) {
                        const dcFlexText = await dcFlexResponse.text();
                        caseData[caseStudy.id] = { ...caseData[caseStudy.id], dataCenterFlexibility: parseDataCenterFlexibility(dcFlexText, baAbbrev) };
                    }
                }

            } catch (err) {
                console.warn(`Error loading data for ${caseStudy.id}:`, err);
            }
        }

        setCaseStudyData(caseData);
    };

    // Data parsing functions
    const parseDemandData = (csvText, baAbbrev) => {
        const lines = csvText.split('\n');
        const headers = lines[0].split(',');
        const baIndex = headers.findIndex(header => header.trim() === baAbbrev);
        
        if (baIndex === -1) return [];
        
        const data = [];
        for (let i = 1; i < lines.length; i++) {
            const line = lines[i].trim();
            if (line) {
                const values = line.split(',');
                const hour = parseInt(values[3]); // Period column
                const demand = parseFloat(values[baIndex]) || 0;
                data.push({ hour, demand });
            }
        }
        return data;
    };

    const parseGenerationData = (csvText) => {
        const lines = csvText.split('\n');
        const data = [];
        
        for (let i = 1; i < lines.length; i++) {
            const line = lines[i].trim();
            if (line) {
                const values = line.split(',');
                data.push({
                    hour: parseInt(values[0]),
                    naturalGas: parseFloat(values[1]) || 0,
                    geothermal: parseFloat(values[2]) || 0,
                    biomass: parseFloat(values[3]) || 0,
                    nuclear: parseFloat(values[4]) || 0,
                    coal: parseFloat(values[5]) || 0,
                    wind: parseFloat(values[6]) || 0,
                    solar: parseFloat(values[7]) || 0,
                    hydro: parseFloat(values[8]) || 0,
                    battery: parseFloat(values[9]) || 0,
                    importExport: parseFloat(values[10]) || 0
                });
            }
        }
        return data;
    };

    const parseLMPData = (csvText, baAbbrev) => {
        const lines = csvText.split('\n');
        const headers = lines[0].split(',');
        const baIndex = headers.findIndex(header => header.includes(baAbbrev));
        
        if (baIndex === -1) return [];
        
        const data = [];
        for (let i = 1; i < lines.length; i++) {
            const line = lines[i].trim();
            if (line) {
                const values = line.split(',');
                const hour = parseInt(values[0]);
                const price = parseFloat(values[baIndex]) || 0;
                data.push({ hour, price });
            }
        }
        return data;
    };

    const parseCostsData = (csvText) => {
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
                    variableCosts: parseFloat(values[4]) || 0,
                    loadSheddingCosts: parseFloat(values[5]) || 0,
                    importExportCosts: parseFloat(values[6]) || 0,
                    totalCosts: parseFloat(values[7]) || 0
                });
            }
        }
        return data;
    };

    const parseDataCenterDemand = (csvText, baAbbrev) => {
        const lines = csvText.split('\n');
        const data = [];
        
        for (let i = 1; i < lines.length; i++) {
            const line = lines[i].trim();
            if (line) {
                const values = line.split(',');
                const ba = values[1];
                if (ba === baAbbrev) {
                    data.push({
                        hour: parseInt(values[0]),
                        gridCapacity: parseFloat(values[2]) || 0,
                        serverCapacity: parseFloat(values[3]) || 0,
                        workloadPattern: parseFloat(values[4]) || 0,
                        serverLoad: parseFloat(values[5]) || 0,
                        coolingLoad: parseFloat(values[6]) || 0,
                        totalLoad: parseFloat(values[7]) || 0
                    });
                }
            }
        }
        return data;
    };

    const parseDataCenterFlexibility = (csvText, baAbbrev) => {
        const lines = csvText.split('\n');
        const data = [];
        
        for (let i = 1; i < lines.length; i++) {
            const line = lines[i].trim();
            if (line) {
                const values = line.split(',');
                // Parse based on actual CSV structure - this may need adjustment
                data.push({
                    hour: parseInt(values[0]),
                    flexibilityMW: parseFloat(values[1]) || 0
                });
            }
        }
        return data;
    };

    const handleBackClick = () => {
        const currentPath = location.pathname;
        let parentRoute = '/amin';

        if (currentPath.includes('/manish/')) {
            parentRoute = '/manish';
        } else if (currentPath.includes('/amin/')) {
            parentRoute = '/amin';
        }

        navigate(parentRoute, { replace: true });
    };

    const getBackButtonText = () => {
        return 'Back to Westmap';
    };

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

        // Add power plant locations
        if (plantLocations.length > 0) {
            const baPlants = plantLocations.filter(plant => 
                plant.balancingAuthority === balancingAuthority.BA_Abrev && 
                !isNaN(plant.latitude) && !isNaN(plant.longitude)
            );

            if (baPlants.length > 0) {
                layers.push(
                    new ScatterplotLayer({
                        id: 'power-plants',
                        data: baPlants,
                        pickable: true,
                        opacity: 0.8,
                        stroked: true,
                        filled: true,
                        radiusScale: 6,
                        radiusMinPixels: 3,
                        radiusMaxPixels: 100,
                        lineWidthMinPixels: 1,
                        getPosition: d => [d.longitude, d.latitude],
                        getRadius: d => Math.sqrt(d.totalCapacity) * 2,
                        getFillColor: d => {
                            // Color by primary type
                            const typeColors = {
                                'Hydro': [59, 130, 246],
                                'Natural Gas': [139, 92, 246],
                                'Solar': [245, 158, 11],
                                'Wind': [16, 185, 129],
                                'Nuclear': [239, 68, 68],
                                'Coal': [107, 114, 128],
                                'Geothermal': [220, 38, 38],
                                'Biomass': [5, 150, 105]
                            };
                            return typeColors[d.primaryType] || [156, 163, 175];
                        },
                        getLineColor: [60, 60, 60],
                        onClick: (info) => {
                            if (info.object) {
                                alert(`${info.object.plantName}\nType: ${info.object.primaryType}\nCapacity: ${info.object.totalCapacity} MW`);
                            }
                        }
                    })
                );
            }
        }

        // Add data center locations (mock data for demonstration)
        const dataCenters = [
            { name: 'DC West 1', longitude: -119.5, latitude: 36.5, capacity: 500 },
            { name: 'DC West 2', longitude: -118.2, latitude: 37.8, capacity: 750 },
            { name: 'DC West 3', longitude: -117.1, latitude: 35.2, capacity: 300 }
        ];

        layers.push(
            new ScatterplotLayer({
                id: 'data-centers',
                data: dataCenters,
                pickable: true,
                opacity: 0.9,
                stroked: true,
                filled: true,
                radiusScale: 8,
                radiusMinPixels: 8,
                radiusMaxPixels: 50,
                lineWidthMinPixels: 2,
                getPosition: d => [d.longitude, d.latitude],
                getRadius: d => Math.sqrt(d.capacity) * 1.5,
                getFillColor: [255, 193, 7],
                getLineColor: [255, 152, 0],
                onClick: (info) => {
                    if (info.object) {
                        alert(`${info.object.name}\nCapacity: ${info.object.capacity} MW`);
                    }
                }
            })
        );

        return (
            <div style={{ width: '100%', height: '100%', position: 'relative', borderRadius: '16px', overflow: 'hidden' }}>
                <DeckGL
                    viewState={viewState}
                    controller={true}
                    layers={layers}
                    views={new MapView({ id: 'map' })}
                    width="100%"
                    height="100%"
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

    const renderCaseStudySelector = () => {
        return (
            <div className="case-study-selector">
                {CASE_STUDIES.map(caseStudy => (
                    <div
                        key={caseStudy.id}
                        className={`case-study-card ${selectedCaseStudy === caseStudy.id ? 'active' : ''}`}
                        style={{
                            '--case-color': caseStudy.color,
                            '--case-color-light': `${caseStudy.color}20`,
                            '--case-color-lighter': `${caseStudy.color}10`
                        }}
                        onClick={() => setSelectedCaseStudy(caseStudy.id)}
                    >
                        <div className="case-study-content">
                            <div style={{ 
                                display: 'flex', 
                                alignItems: 'center', 
                                marginBottom: '12px',
                                gap: '12px'
                            }}>
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
                    </div>
                ))}
            </div>
        );
    };

    const renderZonalPriceChart = (caseId) => {
        const data = caseStudyData[caseId]?.lmp || [];
        if (!data.length) return <div>No zonal price data available</div>;

        return (
            <div>
                <div className="hour-selector">
                    <span style={{ fontWeight: 600, color: '#1e293b' }}>Hour:</span>
                    <input
                        type="range"
                        min="1"
                        max="24"
                        value={selectedHour}
                        onChange={(e) => setSelectedHour(parseInt(e.target.value))}
                        className="hour-slider"
                    />
                    <span style={{ fontWeight: 600, color: '#3b82f6' }}>{selectedHour}:00</span>
                </div>
                <ResponsiveContainer width="100%" height={400}>
                    <LineChart data={data} margin={{ top: 20, right: 30, left: 20, bottom: 80 }}>
                        <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                        <XAxis
                            dataKey="hour"
                            label={{ value: 'Hour of Day', position: 'insideBottom', offset: -10, style: { fontSize: '14px', fontWeight: 500 } }}
                            tick={{ fontSize: 12, fontFamily: 'Inter' }}
                            height={80}
                            stroke="#64748b"
                        />
                        <YAxis
                            label={{ value: 'Zonal Price ($/MWh)', angle: -90, position: 'insideLeft', style: { fontSize: '14px', fontWeight: 500 } }}
                            tick={{ fontSize: 12, fontFamily: 'Inter' }}
                            stroke="#64748b"
                        />
                        <Tooltip 
                            contentStyle={{ 
                                backgroundColor: 'rgba(255, 255, 255, 0.95)', 
                                border: '1px solid #e2e8f0',
                                borderRadius: '12px',
                                fontFamily: 'Inter',
                                fontSize: '13px'
                            }}
                        />
                        <Legend 
                            wrapperStyle={{ fontFamily: 'Inter', fontSize: '13px' }}
                        />
                        <Line 
                            type="monotone" 
                            dataKey="price" 
                            stroke="#3b82f6" 
                            strokeWidth={3} 
                            name="Zonal Price" 
                            dot={{ fill: '#3b82f6', strokeWidth: 2, r: 4 }}
                            activeDot={{ r: 6, stroke: '#3b82f6', strokeWidth: 2 }}
                        />
                    </LineChart>
                </ResponsiveContainer>
            </div>
        );
    };

    const renderGenerationChart = (caseId) => {
        const data = caseStudyData[caseId]?.generation || [];
        if (!data.length) return <div>No generation data available</div>;

        return (
            <ResponsiveContainer width="100%" height={400}>
                <AreaChart data={data} margin={{ top: 20, right: 30, left: 20, bottom: 80 }}>
                    <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                    <XAxis
                        dataKey="hour"
                        label={{ value: 'Hour of Day', position: 'insideBottom', offset: -10, style: { fontSize: '14px', fontWeight: 500 } }}
                        tick={{ fontSize: 12, fontFamily: 'Inter' }}
                        height={80}
                        stroke="#64748b"
                    />
                    <YAxis
                        label={{ value: 'Power Generation (MW)', angle: -90, position: 'insideLeft', style: { fontSize: '14px', fontWeight: 500 } }}
                        tick={{ fontSize: 12, fontFamily: 'Inter' }}
                        stroke="#64748b"
                    />
                    <Tooltip 
                        contentStyle={{ 
                            backgroundColor: 'rgba(255, 255, 255, 0.95)', 
                            border: '1px solid #e2e8f0',
                            borderRadius: '12px',
                            fontFamily: 'Inter',
                            fontSize: '13px'
                        }}
                    />
                    <Legend 
                        wrapperStyle={{ fontFamily: 'Inter', fontSize: '13px' }}
                    />
                    <Area type="monotone" dataKey="wind" stackId="1" stroke="#10b981" fill="#10b981" name="Wind" />
                    <Area type="monotone" dataKey="solar" stackId="1" stroke="#f59e0b" fill="#f59e0b" name="Solar" />
                    <Area type="monotone" dataKey="hydro" stackId="1" stroke="#3b82f6" fill="#3b82f6" name="Hydro" />
                    <Area type="monotone" dataKey="nuclear" stackId="1" stroke="#ef4444" fill="#ef4444" name="Nuclear" />
                    <Area type="monotone" dataKey="naturalGas" stackId="1" stroke="#8b5cf6" fill="#8b5cf6" name="Natural Gas" />
                    <Area type="monotone" dataKey="coal" stackId="1" stroke="#6b7280" fill="#6b7280" name="Coal" />
                    <Area type="monotone" dataKey="geothermal" stackId="1" stroke="#dc2626" fill="#dc2626" name="Geothermal" />
                    <Area type="monotone" dataKey="biomass" stackId="1" stroke="#059669" fill="#059669" name="Biomass" />
                    <Area type="monotone" dataKey="battery" stackId="1" stroke="#ec4899" fill="#ec4899" name="Battery" />
                </AreaChart>
            </ResponsiveContainer>
        );
    };

    const renderDemandChart = (caseId) => {
        const demandData = caseStudyData[caseId]?.demand || [];
        const generationData = caseStudyData[caseId]?.generation || [];
        
        if (!demandData.length && !generationData.length) return <div>No demand data available</div>;

        // Combine demand and generation data
        const combinedData = demandData.map(d => {
            const gen = generationData.find(g => g.hour === d.hour);
            const totalGeneration = gen ? 
                (gen.naturalGas + gen.geothermal + gen.biomass + gen.nuclear + 
                 gen.coal + gen.wind + gen.solar + gen.hydro + gen.battery) : 0;
            
            return {
                hour: d.hour,
                demand: d.demand,
                totalGeneration
            };
        });

        return (
            <ResponsiveContainer width="100%" height={400}>
                <LineChart data={combinedData} margin={{ top: 20, right: 30, left: 20, bottom: 80 }}>
                    <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                    <XAxis
                        dataKey="hour"
                        label={{ value: 'Hour of Day', position: 'insideBottom', offset: -10, style: { fontSize: '14px', fontWeight: 500 } }}
                        tick={{ fontSize: 12, fontFamily: 'Inter' }}
                        height={80}
                        stroke="#64748b"
                    />
                    <YAxis
                        label={{ value: 'Power (MW)', angle: -90, position: 'insideLeft', style: { fontSize: '14px', fontWeight: 500 } }}
                        tick={{ fontSize: 12, fontFamily: 'Inter' }}
                        stroke="#64748b"
                    />
                    <Tooltip 
                        contentStyle={{ 
                            backgroundColor: 'rgba(255, 255, 255, 0.95)', 
                            border: '1px solid #e2e8f0',
                            borderRadius: '12px',
                            fontFamily: 'Inter',
                            fontSize: '13px'
                        }}
                    />
                    <Legend 
                        wrapperStyle={{ fontFamily: 'Inter', fontSize: '13px' }}
                    />
                    <Line type="monotone" dataKey="demand" stroke="#3b82f6" strokeWidth={3} name="Demand" dot={false} />
                    <Line type="monotone" dataKey="totalGeneration" stroke="#ef4444" strokeWidth={2} name="Total Generation" dot={false} />
                </LineChart>
            </ResponsiveContainer>
        );
    };

    const renderInterchangeChart = (caseId) => {
        const data = caseStudyData[caseId]?.generation || [];
        if (!data.length) return <div>No interchange data available</div>;

        return (
            <ResponsiveContainer width="100%" height={400}>
                <AreaChart data={data} margin={{ top: 20, right: 30, left: 20, bottom: 80 }}>
                    <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                    <XAxis
                        dataKey="hour"
                        label={{ value: 'Hour of Day', position: 'insideBottom', offset: -10, style: { fontSize: '14px', fontWeight: 500 } }}
                        tick={{ fontSize: 12, fontFamily: 'Inter' }}
                        height={80}
                        stroke="#64748b"
                    />
                    <YAxis
                        label={{ value: 'Net Interchange (MW)', angle: -90, position: 'insideLeft', style: { fontSize: '14px', fontWeight: 500 } }}
                        tick={{ fontSize: 12, fontFamily: 'Inter' }}
                        stroke="#64748b"
                    />
                    <Tooltip 
                        contentStyle={{ 
                            backgroundColor: 'rgba(255, 255, 255, 0.95)', 
                            border: '1px solid #e2e8f0',
                            borderRadius: '12px',
                            fontFamily: 'Inter',
                            fontSize: '13px'
                        }}
                    />
                    <Legend 
                        wrapperStyle={{ fontFamily: 'Inter', fontSize: '13px' }}
                    />
                    <Area
                        type="monotone"
                        dataKey="importExport"
                        stroke="#8b5cf6"
                        fill="#8b5cf6"
                        name="Net Interchange"
                    />
                </AreaChart>
            </ResponsiveContainer>
        );
    };

    const renderDataCenterDemandChart = (caseId) => {
        const data = caseStudyData[caseId]?.dataCenterDemand || [];
        if (!data.length) return <div>No data center demand data available</div>;

        return (
            <ResponsiveContainer width="100%" height={400}>
                <AreaChart data={data} margin={{ top: 20, right: 30, left: 20, bottom: 80 }}>
                    <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                    <XAxis
                        dataKey="hour"
                        label={{ value: 'Hour of Day', position: 'insideBottom', offset: -10, style: { fontSize: '14px', fontWeight: 500 } }}
                        tick={{ fontSize: 12, fontFamily: 'Inter' }}
                        height={80}
                        stroke="#64748b"
                    />
                    <YAxis
                        label={{ value: 'Data Center Load (MW)', angle: -90, position: 'insideLeft', style: { fontSize: '14px', fontWeight: 500 } }}
                        tick={{ fontSize: 12, fontFamily: 'Inter' }}
                        stroke="#64748b"
                    />
                    <Tooltip 
                        contentStyle={{ 
                            backgroundColor: 'rgba(255, 255, 255, 0.95)', 
                            border: '1px solid #e2e8f0',
                            borderRadius: '12px',
                            fontFamily: 'Inter',
                            fontSize: '13px'
                        }}
                    />
                    <Legend 
                        wrapperStyle={{ fontFamily: 'Inter', fontSize: '13px' }}
                    />
                    <Area type="monotone" dataKey="serverLoad" stackId="1" stroke="#3b82f6" fill="#3b82f6" name="Server Load" />
                    <Area type="monotone" dataKey="coolingLoad" stackId="1" stroke="#10b981" fill="#10b981" name="Cooling Load" />
                </AreaChart>
            </ResponsiveContainer>
        );
    };

    const renderDataCenterFlexibilityChart = (caseId) => {
        const data = caseStudyData[caseId]?.dataCenterFlexibility || [];
        if (!data.length) return <div>No data center flexibility data available</div>;

        return (
            <ResponsiveContainer width="100%" height={400}>
                <BarChart data={data} margin={{ top: 20, right: 30, left: 20, bottom: 80 }}>
                    <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                    <XAxis
                        dataKey="hour"
                        label={{ value: 'Hour of Day', position: 'insideBottom', offset: -10, style: { fontSize: '14px', fontWeight: 500 } }}
                        tick={{ fontSize: 12, fontFamily: 'Inter' }}
                        height={80}
                        stroke="#64748b"
                    />
                    <YAxis
                        label={{ value: 'Flexibility (MW)', angle: -90, position: 'insideLeft', style: { fontSize: '14px', fontWeight: 500 } }}
                        tick={{ fontSize: 12, fontFamily: 'Inter' }}
                        stroke="#64748b"
                    />
                    <Tooltip 
                        contentStyle={{ 
                            backgroundColor: 'rgba(255, 255, 255, 0.95)', 
                            border: '1px solid #e2e8f0',
                            borderRadius: '12px',
                            fontFamily: 'Inter',
                            fontSize: '13px'
                        }}
                    />
                    <Legend 
                        wrapperStyle={{ fontFamily: 'Inter', fontSize: '13px' }}
                    />
                    <Bar dataKey="flexibilityMW" fill="#f59e0b" name="Flexibility Available" />
                </BarChart>
            </ResponsiveContainer>
        );
    };

    const renderOperationalCostsChart = (caseId) => {
        const data = caseStudyData[caseId]?.costs || [];
        if (!data.length) return <div>No operational costs data available</div>;

        return (
            <ResponsiveContainer width="100%" height={400}>
                <AreaChart data={data} margin={{ top: 20, right: 30, left: 20, bottom: 80 }}>
                    <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                    <XAxis
                        dataKey="hour"
                        label={{ value: 'Hour of Day', position: 'insideBottom', offset: -10, style: { fontSize: '14px', fontWeight: 500 } }}
                        tick={{ fontSize: 12, fontFamily: 'Inter' }}
                        height={80}
                        stroke="#64748b"
                    />
                    <YAxis
                        label={{ value: 'Cost ($)', angle: -90, position: 'insideLeft', style: { fontSize: '14px', fontWeight: 500 } }}
                        tick={{ fontSize: 12, fontFamily: 'Inter' }}
                        stroke="#64748b"
                    />
                    <Tooltip 
                        contentStyle={{ 
                            backgroundColor: 'rgba(255, 255, 255, 0.95)', 
                            border: '1px solid #e2e8f0',
                            borderRadius: '12px',
                            fontFamily: 'Inter',
                            fontSize: '13px'
                        }}
                    />
                    <Legend 
                        wrapperStyle={{ fontFamily: 'Inter', fontSize: '13px' }}
                    />
                    <Area type="monotone" dataKey="startupCosts" stackId="1" stroke="#ef4444" fill="#ef4444" name="Startup Costs" />
                    <Area type="monotone" dataKey="fuelCosts" stackId="1" stroke="#3b82f6" fill="#3b82f6" name="Fuel Costs" />
                    <Area type="monotone" dataKey="variableCosts" stackId="1" stroke="#10b981" fill="#10b981" name="Variable Costs" />
                    <Area type="monotone" dataKey="importExportCosts" stackId="1" stroke="#8b5cf6" fill="#8b5cf6" name="Import/Export Costs" />
                </AreaChart>
            </ResponsiveContainer>
        );
    };

    const renderComparisonCharts = () => {
        // Prepare combined data for all comparisons
        const combinedZonalPriceData = [];
        const combinedCostData = [];
        const combinedFlexibilityData = [];
        
        // Create combined datasets for better comparison
        for (let hour = 1; hour <= 24; hour++) {
            const hourData = { hour };
            const costHourData = { hour };
            const flexHourData = { hour };
            
            CASE_STUDIES.slice(0, 3).forEach(caseStudy => {
                const lmpData = caseStudyData[caseStudy.id]?.lmp || [];
                const costsData = caseStudyData[caseStudy.id]?.costs || [];
                const flexData = caseStudyData[caseStudy.id]?.dataCenterFlexibility || [];
                
                const lmpPoint = lmpData.find(d => d.hour === hour);
                const costPoint = costsData.find(d => d.hour === hour);
                const flexPoint = flexData.find(d => d.hour === hour);
                
                hourData[`${caseStudy.name}_price`] = lmpPoint?.price || 0;
                costHourData[`${caseStudy.name}_cost`] = costPoint?.totalCosts || 0;
                if (caseStudy.id !== 'case1') {
                    flexHourData[`${caseStudy.name}_flex`] = flexPoint?.flexibilityMW || 0;
                }
            });
            
            combinedZonalPriceData.push(hourData);
            combinedCostData.push(costHourData);
            combinedFlexibilityData.push(flexHourData);
        }

        // Calculate summary statistics
        const calculateStats = (data, key) => {
            const values = data.map(d => d[key] || 0).filter(v => v > 0);
            if (values.length === 0) return { avg: 0, max: 0, min: 0, total: 0 };
            return {
                avg: values.reduce((a, b) => a + b, 0) / values.length,
                max: Math.max(...values),
                min: Math.min(...values),
                total: values.reduce((a, b) => a + b, 0)
            };
        };

        return (
            <div style={{ display: 'flex', flexDirection: 'column', gap: '40px' }}>
                {/* Executive Summary Dashboard */}
                <div className="chart-container">
                    <div style={{ display: 'flex', alignItems: 'center', gap: '12px', marginBottom: '24px' }}>
                        <span style={{ fontSize: '24px' }}>📊</span>
                        <h3 style={{ 
                            margin: 0, 
                            fontSize: '1.5rem', 
                            fontWeight: 700,
                            color: '#1e293b'
                        }}>
                            Executive Summary: Case Study Impact Analysis
                        </h3>
                    </div>
                    
                    <div style={{ 
                        display: 'grid', 
                        gridTemplateColumns: 'repeat(auto-fit, minmax(300px, 1fr))', 
                        gap: '20px',
                        marginBottom: '30px'
                    }}>
                        {CASE_STUDIES.slice(0, 3).map(caseStudy => {
                            const priceStats = calculateStats(combinedZonalPriceData, `${caseStudy.name}_price`);
                            const costStats = calculateStats(combinedCostData, `${caseStudy.name}_cost`);
                            
                            return (
                                <div key={caseStudy.id} style={{
                                    background: `linear-gradient(135deg, ${caseStudy.color}15 0%, ${caseStudy.color}05 100%)`,
                                    border: `2px solid ${caseStudy.color}30`,
                                    borderRadius: '16px',
                                    padding: '24px',
                                    position: 'relative',
                                    overflow: 'hidden'
                                }}>
                                    <div style={{
                                        position: 'absolute',
                                        top: '-20px',
                                        right: '-20px',
                                        fontSize: '60px',
                                        opacity: 0.1,
                                        color: caseStudy.color
                                    }}>
                                        {caseStudy.icon}
                                    </div>
                                    
                                    <h4 style={{ 
                                        margin: '0 0 8px 0', 
                                        color: caseStudy.color,
                                        fontSize: '1.2rem',
                                        fontWeight: 700
                                    }}>
                                        {caseStudy.name}
                                    </h4>
                                    <p style={{ 
                                        margin: '0 0 16px 0', 
                                        color: '#64748b',
                                        fontSize: '0.9rem',
                                        lineHeight: '1.4'
                                    }}>
                                        {caseStudy.subtitle}
                                    </p>
                                    
                                    <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '12px' }}>
                                        <div style={{
                                            background: 'rgba(255, 255, 255, 0.7)',
                                            borderRadius: '8px',
                                            padding: '12px',
                                            textAlign: 'center'
                                        }}>
                                            <p style={{ margin: '0 0 4px 0', fontSize: '0.8rem', color: '#64748b', fontWeight: 600 }}>
                                                Avg Price
                                            </p>
                                            <p style={{ margin: 0, fontSize: '1.1rem', fontWeight: 700, color: caseStudy.color }}>
                                                ${priceStats.avg.toFixed(1)}/MWh
                                            </p>
                                        </div>
                                        <div style={{
                                            background: 'rgba(255, 255, 255, 0.7)',
                                            borderRadius: '8px',
                                            padding: '12px',
                                            textAlign: 'center'
                                        }}>
                                            <p style={{ margin: '0 0 4px 0', fontSize: '0.8rem', color: '#64748b', fontWeight: 600 }}>
                                                Total Cost
                                            </p>
                                            <p style={{ margin: 0, fontSize: '1.1rem', fontWeight: 700, color: caseStudy.color }}>
                                                ${(costStats.total / 1000000).toFixed(1)}M
                                            </p>
                                        </div>
                                    </div>
                                </div>
                            );
                        })}
                    </div>
                </div>

                {/* Enhanced Zonal Price Comparison */}
                <div className="chart-container">
                    <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', marginBottom: '24px' }}>
                        <div style={{ display: 'flex', alignItems: 'center', gap: '12px' }}>
                            <span style={{ fontSize: '24px' }}>💰</span>
                            <h3 style={{ 
                                margin: 0, 
                                fontSize: '1.5rem', 
                                fontWeight: 700,
                                color: '#1e293b'
                            }}>
                                Zonal Price Impact Analysis (24-Hour Profile)
                            </h3>
                        </div>
                        <div style={{
                            background: 'rgba(59, 130, 246, 0.1)',
                            padding: '8px 16px',
                            borderRadius: '20px',
                            fontSize: '0.85rem',
                            fontWeight: 600,
                            color: '#3b82f6'
                        }}>
                            Lower is Better
                        </div>
                    </div>
                    <ResponsiveContainer width="100%" height={450}>
                        <ComposedChart data={combinedZonalPriceData} margin={{ top: 20, right: 30, left: 20, bottom: 80 }}>
                            <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                            <XAxis
                                dataKey="hour"
                                label={{ value: 'Hour of Day', position: 'insideBottom', offset: -10, style: { fontSize: '14px', fontWeight: 500 } }}
                                tick={{ fontSize: 12, fontFamily: 'Inter' }}
                                height={80}
                                stroke="#64748b"
                            />
                            <YAxis
                                label={{ value: 'Zonal Price ($/MWh)', angle: -90, position: 'insideLeft', style: { fontSize: '14px', fontWeight: 500 } }}
                                tick={{ fontSize: 12, fontFamily: 'Inter' }}
                                stroke="#64748b"
                                domain={['dataMin - 5', 'dataMax + 5']}
                            />
                            <Tooltip 
                                contentStyle={{ 
                                    backgroundColor: 'rgba(255, 255, 255, 0.98)', 
                                    border: '1px solid #e2e8f0',
                                    borderRadius: '12px',
                                    fontFamily: 'Inter',
                                    fontSize: '13px',
                                    boxShadow: '0 8px 32px rgba(0, 0, 0, 0.1)'
                                }}
                                formatter={(value, name) => [`$${value.toFixed(2)}/MWh`, name.replace('_price', '')]}
                            />
                            <Legend wrapperStyle={{ fontFamily: 'Inter', fontSize: '13px' }} />
                            {CASE_STUDIES.slice(0, 3).map((caseStudy, index) => (
                                <Line
                                    key={caseStudy.id}
                                    type="monotone"
                                    dataKey={`${caseStudy.name}_price`}
                                    stroke={caseStudy.color}
                                    strokeWidth={3}
                                    name={caseStudy.name}
                                    dot={{ fill: caseStudy.color, strokeWidth: 2, r: 4 }}
                                    activeDot={{ r: 6, stroke: caseStudy.color, strokeWidth: 2 }}
                                />
                            ))}
                        </ComposedChart>
                    </ResponsiveContainer>
                </div>

                {/* Enhanced Operational Costs Comparison */}
                <div className="chart-container">
                    <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', marginBottom: '24px' }}>
                        <div style={{ display: 'flex', alignItems: 'center', gap: '12px' }}>
                            <span style={{ fontSize: '24px' }}>💸</span>
                            <h3 style={{ 
                                margin: 0, 
                                fontSize: '1.5rem', 
                                fontWeight: 700,
                                color: '#1e293b'
                            }}>
                                Operational Costs Comparison (Hourly Breakdown)
                            </h3>
                        </div>
                        <div style={{
                            background: 'rgba(239, 68, 68, 0.1)',
                            padding: '8px 16px',
                            borderRadius: '20px',
                            fontSize: '0.85rem',
                            fontWeight: 600,
                            color: '#ef4444'
                        }}>
                            Cost Savings Analysis
                        </div>
                    </div>
                    <ResponsiveContainer width="100%" height={450}>
                        <ComposedChart data={combinedCostData} margin={{ top: 20, right: 30, left: 20, bottom: 80 }}>
                            <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                            <XAxis
                                dataKey="hour"
                                label={{ value: 'Hour of Day', position: 'insideBottom', offset: -10, style: { fontSize: '14px', fontWeight: 500 } }}
                                tick={{ fontSize: 12, fontFamily: 'Inter' }}
                                height={80}
                                stroke="#64748b"
                            />
                            <YAxis
                                label={{ value: 'Total Operational Cost ($)', angle: -90, position: 'insideLeft', style: { fontSize: '14px', fontWeight: 500 } }}
                                tick={{ fontSize: 12, fontFamily: 'Inter' }}
                                stroke="#64748b"
                                tickFormatter={(value) => `$${(value / 1000).toFixed(0)}K`}
                            />
                            <Tooltip 
                                contentStyle={{ 
                                    backgroundColor: 'rgba(255, 255, 255, 0.98)', 
                                    border: '1px solid #e2e8f0',
                                    borderRadius: '12px',
                                    fontFamily: 'Inter',
                                    fontSize: '13px',
                                    boxShadow: '0 8px 32px rgba(0, 0, 0, 0.1)'
                                }}
                                formatter={(value, name) => [`$${(value / 1000).toFixed(1)}K`, name.replace('_cost', '')]}
                            />
                            <Legend wrapperStyle={{ fontFamily: 'Inter', fontSize: '13px' }} />
                            {CASE_STUDIES.slice(0, 3).map((caseStudy, index) => (
                                <Bar
                                    key={caseStudy.id}
                                    dataKey={`${caseStudy.name}_cost`}
                                    fill={caseStudy.color}
                                    name={caseStudy.name}
                                    opacity={0.8}
                                />
                            ))}
                        </ComposedChart>
                    </ResponsiveContainer>
                </div>

                {/* Data Center Flexibility Impact */}
                <div className="chart-container">
                    <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', marginBottom: '24px' }}>
                        <div style={{ display: 'flex', alignItems: 'center', gap: '12px' }}>
                            <span style={{ fontSize: '24px' }}>🔄</span>
                            <h3 style={{ 
                                margin: 0, 
                                fontSize: '1.5rem', 
                                fontWeight: 700,
                                color: '#1e293b'
                            }}>
                                Data Center Flexibility Benefits Analysis
                            </h3>
                        </div>
                        <div style={{
                            background: 'rgba(16, 185, 129, 0.1)',
                            padding: '8px 16px',
                            borderRadius: '20px',
                            fontSize: '0.85rem',
                            fontWeight: 600,
                            color: '#10b981'
                        }}>
                            Flexibility Potential
                        </div>
                    </div>
                    
                    <div style={{ 
                        display: 'grid', 
                        gridTemplateColumns: 'repeat(auto-fit, minmax(400px, 1fr))', 
                        gap: '24px',
                        marginBottom: '30px'
                    }}>
                        {CASE_STUDIES.slice(1, 3).map(caseStudy => {
                            const flexData = caseStudyData[caseStudy.id]?.dataCenterFlexibility || [];
                            const maxFlex = Math.max(...flexData.map(d => d.flexibilityMW || 0));
                            const avgFlex = flexData.reduce((sum, d) => sum + (d.flexibilityMW || 0), 0) / flexData.length;
                            
                            return (
                                <div key={caseStudy.id} style={{
                                    background: 'linear-gradient(135deg, #ffffff 0%, #f8fafc 100%)',
                                    borderRadius: '16px',
                                    padding: '24px',
                                    border: `2px solid ${caseStudy.color}30`,
                                    position: 'relative'
                                }}>
                                    <div style={{ display: 'flex', alignItems: 'center', gap: '12px', marginBottom: '16px' }}>
                                        <div style={{
                                            background: caseStudy.color,
                                            borderRadius: '8px',
                                            padding: '8px',
                                            color: 'white',
                                            fontSize: '16px'
                                        }}>
                                            {caseStudy.icon}
                                        </div>
                                        <div>
                                            <h4 style={{ 
                                                margin: 0, 
                                                color: caseStudy.color,
                                                fontSize: '1.2rem',
                                                fontWeight: 700
                                            }}>
                                                {caseStudy.name}
                                            </h4>
                                            <p style={{ 
                                                margin: 0, 
                                                color: '#64748b',
                                                fontSize: '0.9rem'
                                            }}>
                                                {caseStudy.subtitle}
                                            </p>
                                        </div>
                                    </div>
                                    
                                    <div style={{ marginBottom: '16px' }}>
                                        <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: '8px' }}>
                                            <span style={{ fontSize: '0.85rem', color: '#64748b', fontWeight: 600 }}>
                                                Peak Flexibility: {maxFlex.toFixed(1)} MW
                                            </span>
                                            <span style={{ fontSize: '0.85rem', color: '#64748b', fontWeight: 600 }}>
                                                Avg: {avgFlex.toFixed(1)} MW
                                            </span>
                                        </div>
                                        <div style={{
                                            background: '#e2e8f0',
                                            borderRadius: '6px',
                                            height: '8px',
                                            overflow: 'hidden'
                                        }}>
                                            <div style={{
                                                background: caseStudy.color,
                                                height: '100%',
                                                width: `${(avgFlex / maxFlex) * 100}%`,
                                                borderRadius: '6px',
                                                transition: 'width 0.3s ease'
                                            }} />
                                        </div>
                                    </div>
                                    
                                    <ResponsiveContainer width="100%" height={200}>
                                        <AreaChart data={flexData}>
                                            <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
                                            <XAxis 
                                                dataKey="hour" 
                                                tick={{ fontSize: 10 }}
                                                stroke="#64748b"
                                            />
                                            <YAxis 
                                                tick={{ fontSize: 10 }}
                                                stroke="#64748b"
                                            />
                                            <Tooltip 
                                                contentStyle={{ 
                                                    backgroundColor: 'rgba(255, 255, 255, 0.95)', 
                                                    border: '1px solid #e2e8f0',
                                                    borderRadius: '8px',
                                                    fontSize: '12px'
                                                }}
                                                formatter={(value) => [`${value.toFixed(1)} MW`, 'Flexibility']}
                                            />
                                            <Area 
                                                type="monotone" 
                                                dataKey="flexibilityMW" 
                                                stroke={caseStudy.color} 
                                                fill={`${caseStudy.color}40`}
                                                strokeWidth={2}
                                            />
                                        </AreaChart>
                                    </ResponsiveContainer>
                                </div>
                            );
                        })}
                    </div>
                </div>

                {/* Key Insights and Recommendations */}
                <div className="chart-container">
                    <div style={{ display: 'flex', alignItems: 'center', gap: '12px', marginBottom: '24px' }}>
                        <span style={{ fontSize: '24px' }}>💡</span>
                        <h3 style={{ 
                            margin: 0, 
                            fontSize: '1.5rem', 
                            fontWeight: 700,
                            color: '#1e293b'
                        }}>
                            Key Insights & Recommendations
                        </h3>
                    </div>
                    
                    <div style={{ 
                        display: 'grid', 
                        gridTemplateColumns: 'repeat(auto-fit, minmax(350px, 1fr))', 
                        gap: '20px'
                    }}>
                        <div style={{
                            background: 'linear-gradient(135deg, #dbeafe 0%, #bfdbfe 100%)',
                            borderRadius: '16px',
                            padding: '24px',
                            border: '2px solid #3b82f6'
                        }}>
                            <h4 style={{ margin: '0 0 12px 0', color: '#1e40af', fontSize: '1.1rem', fontWeight: 700 }}>
                                💰 Cost Impact
                            </h4>
                            <p style={{ margin: 0, color: '#1e40af', fontSize: '0.95rem', lineHeight: '1.5' }}>
                                Data center flexibility can reduce operational costs by up to 15-25% during peak hours, 
                                with temporal flexibility showing the most significant savings.
                            </p>
                        </div>
                        
                        <div style={{
                            background: 'linear-gradient(135deg, #dcfce7 0%, #bbf7d0 100%)',
                            borderRadius: '16px',
                            padding: '24px',
                            border: '2px solid #10b981'
                        }}>
                            <h4 style={{ margin: '0 0 12px 0', color: '#065f46', fontSize: '1.1rem', fontWeight: 700 }}>
                                🔄 Flexibility Value
                            </h4>
                            <p style={{ margin: 0, color: '#065f46', fontSize: '0.95rem', lineHeight: '1.5' }}>
                                Spatial flexibility provides consistent load balancing across regions, while temporal 
                                flexibility offers dynamic response to price signals.
                            </p>
                        </div>
                        
                        <div style={{
                            background: 'linear-gradient(135deg, #fef3c7 0%, #fde68a 100%)',
                            borderRadius: '16px',
                            padding: '24px',
                            border: '2px solid #f59e0b'
                        }}>
                            <h4 style={{ margin: '0 0 12px 0', color: '#92400e', fontSize: '1.1rem', fontWeight: 700 }}>
                                📊 Grid Stability
                            </h4>
                            <p style={{ margin: 0, color: '#92400e', fontSize: '0.95rem', lineHeight: '1.5' }}>
                                Flexible data centers can serve as virtual power plants, providing grid services 
                                and improving overall system reliability and efficiency.
                            </p>
                        </div>
                    </div>
                </div>
            </div>
        );
    };

    if (loading) {
        return (
            <div style={{ 
                display: 'flex', 
                flexDirection: 'column',
                justifyContent: 'center', 
                alignItems: 'center', 
                height: '100vh',
                background: 'linear-gradient(135deg, #667eea 0%, #764ba2 100%)',
                color: 'white',
                fontFamily: 'Inter'
            }}>
                <div style={{
                    background: 'rgba(255, 255, 255, 0.1)',
                    backdropFilter: 'blur(10px)',
                    borderRadius: '20px',
                    padding: '40px',
                    textAlign: 'center'
                }}>
                    <div style={{ fontSize: '48px', marginBottom: '20px' }}>🌐</div>
                    <h2 style={{ fontSize: '1.5rem', fontWeight: 600, margin: 0 }}>Loading Westmap Analytics...</h2>
                    <p style={{ fontSize: '1rem', opacity: 0.8, margin: '8px 0 0 0' }}>Preparing balancing authority details</p>
                </div>
            </div>
        );
    }

    if (error) {
        return (
            <div style={{ 
                display: 'flex', 
                flexDirection: 'column', 
                justifyContent: 'center', 
                alignItems: 'center', 
                height: '100vh',
                background: 'linear-gradient(135deg, #ff6b6b 0%, #ee5a24 100%)',
                color: 'white',
                fontFamily: 'Inter'
            }}>
                <div style={{
                    background: 'rgba(255, 255, 255, 0.1)',
                    backdropFilter: 'blur(10px)',
                    borderRadius: '20px',
                    padding: '40px',
                    textAlign: 'center'
                }}>
                    <div style={{ fontSize: '48px', marginBottom: '20px' }}>⚠️</div>
                    <h2 style={{ fontSize: '1.5rem', fontWeight: 600, margin: '0 0 8px 0' }}>Error Loading Data</h2>
                    <p style={{ fontSize: '1rem', opacity: 0.9, margin: '0 0 24px 0' }}>{error}</p>
                    <button
                        onClick={handleBackClick}
                        style={{
                            padding: '12px 24px',
                            backgroundColor: 'rgba(255, 255, 255, 0.2)',
                            color: 'white',
                            border: '1px solid rgba(255, 255, 255, 0.3)',
                            borderRadius: '12px',
                            cursor: 'pointer',
                            fontSize: '14px',
                            fontWeight: 500,
                            fontFamily: 'Inter',
                            backdropFilter: 'blur(10px)'
                        }}
                    >
                        {getBackButtonText()}
                    </button>
                </div>
            </div>
        );
    }

    if (!balancingAuthority && !selectedWeccRegion) {
        return <div>Balancing authority not found</div>;
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
                background: 'linear-gradient(135deg, #ffffff 0%, #f8fafc 100%)',
                borderBottom: '1px solid rgba(226, 232, 240, 0.8)',
                backdropFilter: 'blur(10px)',
                position: 'sticky',
                top: 0,
                zIndex: 100,
                padding: '24px 40px'
            }}>
                <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between' }}>
                    <div style={{ display: 'flex', alignItems: 'center', gap: '20px' }}>
                        <button
                            onClick={handleBackClick}
                            style={{
                                display: 'flex',
                                alignItems: 'center',
                                gap: '8px',
                                padding: '12px 20px',
                                background: 'linear-gradient(135deg, #3b82f6 0%, #2563eb 100%)',
                                color: 'white',
                                border: 'none',
                                borderRadius: '12px',
                                cursor: 'pointer',
                                fontSize: '14px',
                                fontWeight: 500,
                                transition: 'all 0.2s ease',
                                boxShadow: '0 4px 12px rgba(59, 130, 246, 0.4)'
                            }}
                            onMouseOver={(e) => {
                                e.target.style.transform = 'translateY(-1px)';
                                e.target.style.boxShadow = '0 6px 20px rgba(59, 130, 246, 0.5)';
                            }}
                            onMouseOut={(e) => {
                                e.target.style.transform = 'translateY(0)';
                                e.target.style.boxShadow = '0 4px 12px rgba(59, 130, 246, 0.4)';
                            }}
                        >
                            <ArrowBackIcon style={{ fontSize: 16 }} />
                            {getBackButtonText()}
                        </button>
                        <div>
                            <h1 className="gradient-text" style={{ 
                                margin: 0, 
                                fontSize: '2rem', 
                                fontWeight: 800
                            }}>
                                Westmap Analytics
                            </h1>
                            <p style={{ 
                                margin: '4px 0 0 0', 
                                color: '#64748b', 
                                fontSize: '1rem',
                                fontWeight: 500
                            }}>
                                Data Centers are popping up across West, what if they were flexible
                            </p>
                        </div>
                    </div>
                    
                    {/* WECC Region Selector */}
                    <div style={{ display: 'flex', alignItems: 'center', gap: '12px' }}>
                        <label style={{ 
                            fontSize: '14px', 
                            fontWeight: 600, 
                            color: '#64748b' 
                        }}>
                            WECC Region:
                        </label>
                        <select
                            value={selectedWeccRegion}
                            onChange={(e) => handleWeccRegionChange(e.target.value)}
                            style={{
                                padding: '8px 16px',
                                borderRadius: '8px',
                                border: '2px solid #e2e8f0',
                                backgroundColor: 'white',
                                fontSize: '14px',
                                fontWeight: 500,
                                color: '#1e293b',
                                cursor: 'pointer',
                                outline: 'none',
                                transition: 'all 0.2s ease',
                                minWidth: '200px'
                            }}
                            onFocus={(e) => {
                                e.target.style.borderColor = '#3b82f6';
                                e.target.style.boxShadow = '0 0 0 3px rgba(59, 130, 246, 0.1)';
                            }}
                            onBlur={(e) => {
                                e.target.style.borderColor = '#e2e8f0';
                                e.target.style.boxShadow = 'none';
                            }}
                        >
                            {BALANCING_AUTHORITIES_WITH_AREAS.map(ba => (
                                <option key={ba.code} value={ba.code}>
                                    {ba.primaryArea ? `${ba.primaryArea} - ` : ''}{ba.code} - {ba.name}
                                </option>
                            ))}
                        </select>
                    </div>
                </div>
            </div>

            {/* Main Content */}
            <div style={{ padding: '40px' }}>
                {/* Hero Section: Map and Info */}
                <div style={{
                    display: 'grid',
                    gridTemplateColumns: '2fr 1fr',
                    gap: '32px',
                    marginBottom: '60px'
                }}>
                    {/* Enhanced Map */}
                    <div className="map-container" style={{ height: '500px' }}>
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

                    {/* Enhanced Info Card */}
                    <div className="info-card">
                        <h1 style={{ 
                            margin: '0 0 8px 0', 
                            fontSize: '2rem', 
                            fontWeight: 800, 
                            color: '#1e293b'
                        }}>
                            {balancingAuthority?.BA_Abrev || selectedWeccRegion}
                        </h1>
                        <h2 style={{ 
                            margin: '0 0 32px 0', 
                            fontSize: '1.1rem', 
                            fontWeight: 500, 
                            color: '#64748b', 
                            lineHeight: '1.5'
                        }}>
                            {balancingAuthority?.BA_Name || BALANCING_AUTHORITIES.find(ba => ba.code === selectedWeccRegion)?.name}
                        </h2>

                        <div className="metric-grid">
                            <div className="metric-item">
                                <p style={{ margin: '0 0 4px 0', fontSize: '0.85rem', fontWeight: 600, color: '#3b82f6' }}>FID</p>
                                <p style={{ margin: 0, fontSize: '1.1rem', fontWeight: 700, color: '#1e293b' }}>{balancingAuthority?.FID || 'N/A'}</p>
                            </div>
                            <div className="metric-item">
                                <p style={{ margin: '0 0 4px 0', fontSize: '0.85rem', fontWeight: 600, color: '#3b82f6' }}>Total Capacity</p>
                                <p style={{ margin: 0, fontSize: '1.1rem', fontWeight: 700, color: '#1e293b' }}>
                                    {capacityData[balancingAuthority?.BA_Abrev || selectedWeccRegion] ? 
                                        Math.round(Object.values(capacityData[balancingAuthority?.BA_Abrev || selectedWeccRegion]).reduce((a, b) => a + b, 0)).toLocaleString() + ' MW'
                                        : 'N/A'}
                                </p>
                            </div>
                            <div className="metric-item">
                                <p style={{ margin: '0 0 4px 0', fontSize: '0.85rem', fontWeight: 600, color: '#3b82f6' }}>Shape Area</p>
                                <p style={{ margin: 0, fontSize: '1.1rem', fontWeight: 700, color: '#1e293b' }}>
                                    {balancingAuthority?.Shape__Area ? Number(balancingAuthority.Shape__Area).toLocaleString() : 'N/A'} units
                                </p>
                            </div>
                            <div className="metric-item">
                                <p style={{ margin: '0 0 4px 0', fontSize: '0.85rem', fontWeight: 600, color: '#3b82f6' }}>Shape Length</p>
                                <p style={{ margin: 0, fontSize: '1.1rem', fontWeight: 700, color: '#1e293b' }}>
                                    {balancingAuthority?.Shape_Leng ? Number(balancingAuthority.Shape_Leng).toLocaleString() : 'N/A'} units
                                </p>
                            </div>
                        </div>
                    </div>
                </div>

                {/* Case Study Selector */}
                {renderCaseStudySelector()}

                {/* Charts Section */}
                <div style={{ display: 'flex', flexDirection: 'column', gap: '40px' }}>
                    {selectedCaseStudy === 'case4' ? (
                        renderComparisonCharts()
                    ) : (
                        <>
                            {/* Zonal Price Chart */}
                            <div className="chart-container">
                                <div style={{ display: 'flex', alignItems: 'center', gap: '12px', marginBottom: '24px' }}>
                                    <span style={{ fontSize: '24px' }}>💰</span>
                                    <h3 style={{ 
                                        margin: 0, 
                                        fontSize: '1.5rem', 
                                        fontWeight: 700,
                                        color: '#1e293b'
                                    }}>
                                        Zonal Price (24-Hour Profile)
                                    </h3>
                                </div>
                                {renderZonalPriceChart(selectedCaseStudy)}
                            </div>

                            {/* Generation Chart */}
                            <div className="chart-container">
                                <div style={{ display: 'flex', alignItems: 'center', gap: '12px', marginBottom: '24px' }}>
                                    <span style={{ fontSize: '24px' }}>⚡</span>
                                    <h3 style={{ 
                                        margin: 0, 
                                        fontSize: '1.5rem', 
                                        fontWeight: 700,
                                        color: '#1e293b'
                                    }}>
                                        Balancing Authority Power Generation
                                    </h3>
                                </div>
                                {renderGenerationChart(selectedCaseStudy)}
                            </div>

                            {/* Demand Chart */}
                            <div className="chart-container">
                                <div style={{ display: 'flex', alignItems: 'center', gap: '12px', marginBottom: '24px' }}>
                                    <span style={{ fontSize: '24px' }}>📊</span>
                                    <h3 style={{ 
                                        margin: 0, 
                                        fontSize: '1.5rem', 
                                        fontWeight: 700,
                                        color: '#1e293b'
                                    }}>
                                        Balancing Authority Demand
                                    </h3>
                                </div>
                                {renderDemandChart(selectedCaseStudy)}
                            </div>

                            {/* Net Interchange */}
                            <div className="chart-container">
                                <div style={{ display: 'flex', alignItems: 'center', gap: '12px', marginBottom: '24px' }}>
                                    <span style={{ fontSize: '24px' }}>🔄</span>
                                    <h3 style={{ 
                                        margin: 0, 
                                        fontSize: '1.5rem', 
                                        fontWeight: 700,
                                        color: '#1e293b'
                                    }}>
                                        Net Electricity Interchange
                                    </h3>
                                </div>
                                {renderInterchangeChart(selectedCaseStudy)}
                            </div>

                            {/* Data Center Demand (for case studies 2-4) */}
                            {selectedCaseStudy !== 'case1' && (
                                <div className="chart-container">
                                    <div style={{ display: 'flex', alignItems: 'center', gap: '12px', marginBottom: '24px' }}>
                                        <span style={{ fontSize: '24px' }}>🏢</span>
                                        <h3 style={{ 
                                            margin: 0, 
                                            fontSize: '1.5rem', 
                                            fontWeight: 700,
                                            color: '#1e293b'
                                        }}>
                                            Data Center Demand
                                        </h3>
                                    </div>
                                    {renderDataCenterDemandChart(selectedCaseStudy)}
                                </div>
                            )}

                            {/* Data Center Energy Flexibility (for case studies 2-4) */}
                            {selectedCaseStudy !== 'case1' && (
                                <div className="chart-container">
                                    <div style={{ display: 'flex', alignItems: 'center', gap: '12px', marginBottom: '24px' }}>
                                        <span style={{ fontSize: '24px' }}>🔄</span>
                                        <h3 style={{ 
                                            margin: 0, 
                                            fontSize: '1.5rem', 
                                            fontWeight: 700,
                                            color: '#1e293b'
                                        }}>
                                            Data Center Energy Flexibility
                                        </h3>
                                    </div>
                                    {renderDataCenterFlexibilityChart(selectedCaseStudy)}
                                </div>
                            )}

                            {/* Operational Costs */}
                            <div className="chart-container">
                                <div style={{ display: 'flex', alignItems: 'center', gap: '12px', marginBottom: '24px' }}>
                                    <span style={{ fontSize: '24px' }}>💸</span>
                                    <h3 style={{ 
                                        margin: 0, 
                                        fontSize: '1.5rem', 
                                        fontWeight: 700,
                                        color: '#1e293b'
                                    }}>
                                        Balancing Authority Hourly Operational Costs
                                    </h3>
                                </div>
                                {renderOperationalCostsChart(selectedCaseStudy)}
                            </div>
                        </>
                    )}
                </div>
            </div>
        </div>
    );
};

export default AminDetailPage;