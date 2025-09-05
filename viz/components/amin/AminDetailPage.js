import React, { useState, useEffect } from 'react';
import { useParams, useNavigate, useLocation } from 'react-router-dom';
import { DeckGL } from '@deck.gl/react';
import { GeoJsonLayer } from '@deck.gl/layers';
import { MapView } from '@deck.gl/core';
import { StaticMap } from 'react-map-gl';
import ArrowBackIcon from '@mui/icons-material/ArrowBack';
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer, AreaChart, Area, BarChart, Bar } from 'recharts';

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
};

// Case study configurations
const CASE_STUDIES = [
    {
        id: 'base',
        name: 'Base Case',
        description: 'Current operational baseline scenario',
        color: '#3B82F6',
        icon: '📊'
    },
    {
        id: 'temporal',
        name: 'Data Center Temporal Flexibility',
        description: 'Time-shifted data center operations',
        color: '#10B981',
        icon: '⏰'
    },
    {
        id: 'spatial',
        name: 'Data Center Spatial Flexibility',
        description: 'Geographic load distribution optimization',
        color: '#F59E0B',
        icon: '🌐'
    },
    {
        id: 'comparison',
        name: 'Comparison',
        description: 'Side-by-side analysis of all scenarios',
        color: '#8B5CF6',
        icon: '📈'
    }
];

const AminDetailPage = () => {
    const { fid } = useParams();
    const navigate = useNavigate();
    const location = useLocation();
    const [balancingAuthority, setBalancingAuthority] = useState(null);
    const [areaData, setAreaData] = useState({});
    const [selectedCaseStudy, setSelectedCaseStudy] = useState('base');
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

                // Load actual generation data based on BA abbreviation
                const actualGenerationData = {};

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

                // Load actual generation data based on BA abbreviation
                if (baData.BA_Abrev) {
                    try {
                        const genResponse = await fetch(`/amin_data/new_data/power_gen_data_areawise_24hr/${baData.BA_Abrev}_generation_by_fuel.csv`);
                        if (genResponse.ok) {
                            const genText = await genResponse.text();
                            const genLines = genText.split('\n');

                            if (baData.Area_Numbers && baData.Area_Numbers.length > 0) {
                                baData.Area_Numbers.forEach(areaNumber => {
                                    actualGenerationData[areaNumber] = [];

                                    for (let i = 1; i < genLines.length; i++) {
                                        const line = genLines[i].trim();
                                        if (line) {
                                            const values = line.split(',');
                                            const hour = parseInt(values[0]);

                                            const hourData = {
                                                hour: hour,
                                                naturalGas: parseFloat(values[1]) || 0,
                                                geothermal: parseFloat(values[2]) || 0,
                                                biomass: parseFloat(values[3]) || 0,
                                                nuclear: parseFloat(values[4]) || 0,
                                                coal: parseFloat(values[5]) || 0,
                                                wind: parseFloat(values[6]) || 0,
                                                solar: parseFloat(values[7]) || 0,
                                                hydro: parseFloat(values[8]) || 0,
                                                battery: parseFloat(values[9]) || 0,
                                                importExport: parseFloat(values[10]) || 0,
                                                totalGeneration: (parseFloat(values[1]) || 0) +
                                                    (parseFloat(values[2]) || 0) +
                                                    (parseFloat(values[3]) || 0) +
                                                    (parseFloat(values[4]) || 0) +
                                                    (parseFloat(values[5]) || 0) +
                                                    (parseFloat(values[6]) || 0) +
                                                    (parseFloat(values[7]) || 0) +
                                                    (parseFloat(values[8]) || 0) +
                                                    (parseFloat(values[9]) || 0),
                                                demand: (parseFloat(values[1]) || 0) +
                                                    (parseFloat(values[2]) || 0) +
                                                    (parseFloat(values[3]) || 0) +
                                                    (parseFloat(values[4]) || 0) +
                                                    (parseFloat(values[5]) || 0) +
                                                    (parseFloat(values[6]) || 0) +
                                                    (parseFloat(values[7]) || 0) +
                                                    (parseFloat(values[8]) || 0) +
                                                    (parseFloat(values[9]) || 0) +
                                                    Math.abs(parseFloat(values[10]) || 0),
                                                lmpPrice: 45 + Math.sin((hour - 1) * Math.PI / 12) * 15 + Math.random() * 10,
                                                netInterchange: parseFloat(values[10]) || 0
                                            };

                                            actualGenerationData[areaNumber].push(hourData);
                                        }
                                    }
                                });
                            }
                        } else {
                            console.warn(`Generation data not found for ${baData.BA_Abrev}`);
                        }
                    } catch (err) {
                        console.warn(`Error loading generation data for ${baData.BA_Abrev}:`, err);
                    }
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
                setAreaData(actualGenerationData);
                setLoading(false);
            } catch (err) {
                console.error('Error loading data:', err);
                setError(err.message);
                setLoading(false);
            }
        };

        loadData();
    }, [fid]);

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
        if (!balancingAuthority) return null;

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

        return (
            <div style={{ width: '100%', height: '100%', position: 'relative', borderRadius: '16px', overflow: 'hidden' }}>
                <DeckGL
                    viewState={viewState}
                    controller={false}
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
                                <h3 style={{ 
                                    margin: 0, 
                                    fontSize: '1.25rem', 
                                    fontWeight: 700,
                                    color: selectedCaseStudy === caseStudy.id ? caseStudy.color : '#1e293b'
                                }}>
                                    {caseStudy.name}
                                </h3>
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

    const renderGenerationChart = (areaNumber) => {
        const data = areaData[areaNumber];
        if (!data) return null;

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

    const renderDemandChart = (areaNumber) => {
        const data = areaData[areaNumber];
        if (!data) return null;

        return (
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

    const renderInterchangeChart = (areaNumber) => {
        const data = areaData[areaNumber];
        if (!data) return null;

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
                        dataKey="netInterchange"
                        stroke="#8b5cf6"
                        fill="#8b5cf6"
                        name="Net Interchange"
                    />
                </AreaChart>
            </ResponsiveContainer>
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

    if (!balancingAuthority) {
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
                                Advanced Grid Intelligence Platform
                            </p>
                        </div>
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
                            {balancingAuthority.BA_Abrev} Territory
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
                            {balancingAuthority.BA_Abrev}
                        </h1>
                        <h2 style={{ 
                            margin: '0 0 32px 0', 
                            fontSize: '1.1rem', 
                            fontWeight: 500, 
                            color: '#64748b', 
                            lineHeight: '1.5'
                        }}>
                            {balancingAuthority.BA_Name}
                        </h2>

                        <div className="metric-grid">
                            <div className="metric-item">
                                <p style={{ margin: '0 0 4px 0', fontSize: '0.85rem', fontWeight: 600, color: '#3b82f6' }}>FID</p>
                                <p style={{ margin: 0, fontSize: '1.1rem', fontWeight: 700, color: '#1e293b' }}>{balancingAuthority.FID}</p>
                            </div>
                            <div className="metric-item">
                                <p style={{ margin: '0 0 4px 0', fontSize: '0.85rem', fontWeight: 600, color: '#3b82f6' }}>Area Numbers</p>
                                <p style={{ margin: 0, fontSize: '1.1rem', fontWeight: 700, color: '#1e293b' }}>
                                    {balancingAuthority.Area_Numbers && balancingAuthority.Area_Numbers.length > 0
                                        ? balancingAuthority.Area_Numbers.join(', ')
                                        : 'N/A'}
                                </p>
                            </div>
                            <div className="metric-item">
                                <p style={{ margin: '0 0 4px 0', fontSize: '0.85rem', fontWeight: 600, color: '#3b82f6' }}>Shape Area</p>
                                <p style={{ margin: 0, fontSize: '1.1rem', fontWeight: 700, color: '#1e293b' }}>
                                    {balancingAuthority.Shape__Area ? Number(balancingAuthority.Shape__Area).toLocaleString() : 'N/A'} units
                                </p>
                            </div>
                            <div className="metric-item">
                                <p style={{ margin: '0 0 4px 0', fontSize: '0.85rem', fontWeight: 600, color: '#3b82f6' }}>Shape Length</p>
                                <p style={{ margin: 0, fontSize: '1.1rem', fontWeight: 700, color: '#1e293b' }}>
                                    {balancingAuthority.Shape_Leng ? Number(balancingAuthority.Shape_Leng).toLocaleString() : 'N/A'} units
                                </p>
                            </div>
                        </div>
                    </div>
                </div>

                {/* Case Studies Section */}
                {balancingAuthority.Area_Numbers && balancingAuthority.Area_Numbers.map(areaNumber => (
                    <div key={areaNumber} style={{ marginBottom: '80px' }}>

                        {/* Case Study Selector */}
                        {renderCaseStudySelector()}

                        {/* Charts Section */}
                        <div style={{ display: 'flex', flexDirection: 'column', gap: '40px' }}>
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
                                {renderGenerationChart(areaNumber)}
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
                                {renderDemandChart(areaNumber)}
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
                                {renderInterchangeChart(areaNumber)}
                            </div>
                        </div>
                    </div>
                ))}
            </div>
        </div>
    );
};

export default AminDetailPage;