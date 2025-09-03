import React, { useState, useEffect } from 'react';
import { useParams, useNavigate, useLocation } from 'react-router-dom';
import { DeckGL } from '@deck.gl/react';
import { GeoJsonLayer } from '@deck.gl/layers';
import { MapView } from '@deck.gl/core';
import { StaticMap } from 'react-map-gl';
import ArrowBackIcon from '@mui/icons-material/ArrowBack';
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer, AreaChart, Area } from 'recharts';

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

const AminDetailPage = () => {
    const { fid } = useParams();
    const navigate = useNavigate();
    const location = useLocation();
    const [balancingAuthority, setBalancingAuthority] = useState(null);
    const [areaData, setAreaData] = useState({});
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

    // Add CSS to override recharts legend positioning
    useEffect(() => {
        const style = document.createElement('style');
        style.textContent = `
            .recharts-legend-wrapper {
                bottom: 28px !important;
            }
        `;
        document.head.appendChild(style);

        return () => {
            document.head.removeChild(style);
        };
    }, []);

    // Override body overflow to allow scrolling on this page
    useEffect(() => {
        // Store original overflow value
        const originalOverflow = document.body.style.overflow;

        // Enable scrolling
        document.body.style.overflow = 'auto';

        // Cleanup: restore original overflow when component unmounts
        return () => {
            document.body.style.overflow = originalOverflow;
        };
    }, []);

    // Handle browser back button navigation
    useEffect(() => {
        const handlePopState = (event) => {
            // If user navigates back and we're on a detail page, ensure we go to the correct parent
            const currentPath = window.location.pathname;
            if ((currentPath.includes('/amin/') || currentPath.includes('/manish/')) && !balancingAuthority) {
                let parentRoute = '/amin'; // Default fallback

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

                // Load demo data CSV
                // const demoDataResponse = await fetch('/amin_data/WECC_BA_AREA_DEMODATA.csv');
                // if (!demoDataResponse.ok) {
                //     throw new Error(`HTTP error loading demo data CSV! status: ${demoDataResponse.status}`);
                // }
                // const demoDataText = await demoDataResponse.text();

                // Parse main CSV for shape data
                const csvLines = csvText.split('\n');
                const csvData = {};

                // Create lookup table by FID from main CSV
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

                // Parse demo data CSV - COMMENTED OUT TO USE ACTUAL DATA
                // const demoDataLines = demoDataText.split('\n');
                // const headers = demoDataLines[0].split(',');
                // const demoData = {};

                // for (let i = 1; i < demoDataLines.length; i++) {
                //     const line = demoDataLines[i].trim();
                //     if (line) {
                //         const values = line.split(',');
                //         const areaNumber = parseInt(values[0]);
                //         const hour = parseInt(values[1]);

                //         if (!demoData[areaNumber]) {
                //             demoData[areaNumber] = [];
                //         }

                //         demoData[areaNumber].push({
                //             hour: hour,
                //             wind: parseFloat(values[2]),
                //             solar: parseFloat(values[3]),
                //             hydro: parseFloat(values[4]),
                //             nuclear: parseFloat(values[5]),
                //             naturalGas: parseFloat(values[6]),
                //             coal: parseFloat(values[7]),
                //             other: parseFloat(values[8]),
                //             totalGeneration: parseFloat(values[9]),
                //             demand: parseFloat(values[10]),
                //             lmpPrice: parseFloat(values[11]),
                //             netInterchange: parseFloat(values[12])
                //         });
                //     }
                // }

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
                    ...csvData[targetFid], // Add shape data from main CSV
                    ...areaMappingData[targetFid] // Add area numbers from mapping CSV
                };

                // Load actual generation data based on BA abbreviation
                if (baData.BA_Abrev) {
                    try {
                        const genResponse = await fetch(`/amin_data/new_data/power_gen_data_areawise_24hr/${baData.BA_Abrev}_generation_by_fuel.csv`);
                        if (genResponse.ok) {
                            const genText = await genResponse.text();
                            const genLines = genText.split('\n');

                            // Process each area number that this BA serves
                            if (baData.Area_Numbers && baData.Area_Numbers.length > 0) {
                                baData.Area_Numbers.forEach(areaNumber => {
                                    actualGenerationData[areaNumber] = [];

                                    for (let i = 1; i < genLines.length; i++) {
                                        const line = genLines[i].trim();
                                        if (line) {
                                            const values = line.split(',');
                                            const hour = parseInt(values[0]);

                                            // Map CSV columns to chart data structure
                                            const hourData = {
                                                hour: hour,
                                                naturalGas: parseFloat(values[1]) || 0,      // NG_MW
                                                geothermal: parseFloat(values[2]) || 0,      // GEO_MW
                                                biomass: parseFloat(values[3]) || 0,         // BIO_MW (maps to "other")
                                                nuclear: parseFloat(values[4]) || 0,         // NUCLEAR_MW
                                                coal: parseFloat(values[5]) || 0,            // COAL_MW
                                                wind: parseFloat(values[6]) || 0,            // WIND_MW
                                                solar: parseFloat(values[7]) || 0,           // PV_MW (Solar)
                                                hydro: parseFloat(values[8]) || 0,           // HYDRO_MW
                                                battery: parseFloat(values[9]) || 0,         // BATTERY_MW
                                                importExport: parseFloat(values[10]) || 0,   // IMPORT/EXPORT_MW
                                                // Calculate totals
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
                                                    Math.abs(parseFloat(values[10]) || 0), // Include import/export in demand
                                                // For LMP price, we'll use a simulated curve since it's not in the CSV
                                                lmpPrice: 45 + Math.sin((hour - 1) * Math.PI / 12) * 15 + Math.random() * 10,
                                                // Net interchange is the import/export value
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
                    // For MultiPolygon, use the first polygon's outer ring
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
                    // Fallback to default coordinates if geometry is invalid
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

                    // Validate calculated center
                    if (isNaN(centerLng) || isNaN(centerLat)) {
                        console.warn(`Invalid center coordinates for FID ${fid}, using default view`);
                        setViewState(prev => ({
                            ...prev,
                            longitude: -116.5,
                            latitude: 37.5,
                            zoom: 5
                        }));
                    } else {
                        // Calculate zoom level to fit the entire region
                        const lngDiff = maxLng - minLng;
                        const latDiff = maxLat - minLat;
                        const maxDiff = Math.max(lngDiff, latDiff);

                        // Adjust zoom based on the size of the region
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
                setAreaData(actualGenerationData); // Use actual generation data instead of demo data
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
        // Determine the parent route based on current path
        const currentPath = location.pathname;
        let parentRoute = '/amin'; // Default fallback

        if (currentPath.includes('/manish/')) {
            parentRoute = '/manish';
        } else if (currentPath.includes('/amin/')) {
            parentRoute = '/amin';
        }

        // Use replace instead of navigate to avoid creating additional history entries
        navigate(parentRoute, { replace: true });
    };

    const getBackButtonText = () => {
        // const currentPath = location.pathname;
        // if (currentPath.includes('/manish/')) {
        //     return 'Back to Manish\'s Map';
        // } else if (currentPath.includes('/amin/')) {
        //     return 'Back to Amin\'s Map';
        // }
        return 'Back to Map';
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
            <div style={{ width: '100%', height: '100%', position: 'relative' }}>
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

    const renderGenerationChart = (areaNumber) => {
        const data = areaData[areaNumber];
        if (!data) return null;

        return (
            <ResponsiveContainer width="100%" height={350}>
                <AreaChart data={data} margin={{ top: 20, right: 30, left: 20, bottom: 60 }}>
                    <CartesianGrid strokeDasharray="3 3" />
                    <XAxis
                        dataKey="hour"
                        label={{ value: 'Hour', position: 'insideBottom', offset: -5 }}
                        tick={{ fontSize: 12 }}
                        height={60}
                    />
                    <YAxis
                        label={{ value: 'MW', angle: -90, position: 'insideLeft' }}
                        tick={{ fontSize: 12 }}
                    />
                    <Tooltip />
                    <Legend />
                    <Area type="monotone" dataKey="wind" stackId="1" stroke="#87CEEB" fill="#87CEEB" name="Wind" />
                    <Area type="monotone" dataKey="solar" stackId="1" stroke="#FFD700" fill="#FFD700" name="Solar" />
                    <Area type="monotone" dataKey="hydro" stackId="1" stroke="#4682B4" fill="#4682B4" name="Hydro" />
                    <Area type="monotone" dataKey="nuclear" stackId="1" stroke="#FF6347" fill="#FF6347" name="Nuclear" />
                    <Area type="monotone" dataKey="naturalGas" stackId="1" stroke="#DDA0DD" fill="#DDA0DD" name="Natural Gas" />
                    <Area type="monotone" dataKey="coal" stackId="1" stroke="#696969" fill="#696969" name="Coal" />
                    <Area type="monotone" dataKey="geothermal" stackId="1" stroke="#8B4513" fill="#8B4513" name="Geothermal" />
                    <Area type="monotone" dataKey="biomass" stackId="1" stroke="#228B22" fill="#228B22" name="Biomass" />
                    <Area type="monotone" dataKey="battery" stackId="1" stroke="#FF1493" fill="#FF1493" name="Battery" />
                </AreaChart>
            </ResponsiveContainer>
        );
    };

    const renderDemandChart = (areaNumber) => {
        const data = areaData[areaNumber];
        if (!data) return null;

        return (
            <ResponsiveContainer width="100%" height={350}>
                <LineChart data={data} margin={{ top: 20, right: 30, left: 20, bottom: 60 }}>
                    <CartesianGrid strokeDasharray="3 3" />
                    <XAxis
                        dataKey="hour"
                        label={{ value: 'Hour', position: 'insideBottom', offset: -5 }}
                        tick={{ fontSize: 12 }}
                        height={60}
                    />
                    <YAxis
                        label={{ value: 'MW', angle: -90, position: 'insideLeft' }}
                        tick={{ fontSize: 12 }}
                    />
                    <Tooltip />
                    <Legend />
                    <Line type="monotone" dataKey="demand" stroke="#2563eb" strokeWidth={3} name="Demand" />
                    <Line type="monotone" dataKey="totalGeneration" stroke="#dc2626" strokeWidth={2} name="Total Generation" />
                </LineChart>
            </ResponsiveContainer>
        );
    };

    // LMP Price chart removed as requested
    // const renderPriceChart = (areaNumber) => {
    //     const data = areaData[areaNumber];
    //     if (!data) return null;

    //     return (
    //         <ResponsiveContainer width="100%" height={350}>
    //             <LineChart data={data} margin={{ top: 20, right: 30, left: 20, bottom: 60 }}>
    //                 <CartesianGrid strokeDasharray="3 3" />
    //                 <XAxis
    //                     dataKey="hour"
    //                     label={{ value: 'Hour', position: 'insideBottom', offset: -5 }}
    //                     tick={{ fontSize: 12 }}
    //                     height={60}
    //                 />
    //                 <YAxis
    //                     label={{ value: '$/MWh', angle: -90, position: 'insideLeft' }}
    //                     tick={{ fontSize: 12 }}
    //                 />
    //                 <Tooltip />
    //                 <Legend />
    //                 <Line type="monotone" dataKey="lmpPrice" stroke="#059669" strokeWidth={3} name="LMP Price" />
    //             </LineChart>
    //         </ResponsiveContainer>
    //     );
    // };

    const renderInterchangeChart = (areaNumber) => {
        const data = areaData[areaNumber];
        if (!data) return null;

        return (
            <ResponsiveContainer width="100%" height={350}>
                <AreaChart data={data} margin={{ top: 20, right: 30, left: 20, bottom: 60 }}>
                    <CartesianGrid strokeDasharray="3 3" />
                    <XAxis
                        dataKey="hour"
                        label={{ value: 'Hour', position: 'insideBottom', offset: -5 }}
                        tick={{ fontSize: 12 }}
                        height={60}
                    />
                    <YAxis
                        label={{ value: 'MW', angle: -90, position: 'insideLeft' }}
                        tick={{ fontSize: 12 }}
                    />
                    <Tooltip />
                    <Legend />
                    <Area
                        type="monotone"
                        dataKey="netInterchange"
                        stroke="#8884d8"
                        fill="#8884d8"
                        name="Net Interchange"
                    />
                </AreaChart>
            </ResponsiveContainer>
        );
    };

    if (loading) {
        return (
            <div style={{ display: 'flex', justifyContent: 'center', alignItems: 'center', height: '100vh' }}>
                <h2>Loading balancing authority details...</h2>
            </div>
        );
    }

    if (error) {
        return (
            <div style={{ display: 'flex', flexDirection: 'column', justifyContent: 'center', alignItems: 'center', height: '100vh' }}>
                <h2 style={{ color: 'red' }}>Error: {error}</h2>
                <button
                    onClick={handleBackClick}
                    style={{
                        marginTop: '20px',
                        padding: '10px 20px',
                        backgroundColor: '#2563eb',
                        color: 'white',
                        border: 'none',
                        borderRadius: '8px',
                        cursor: 'pointer'
                    }}
                >
                    {getBackButtonText()}
                </button>
            </div>
        );
    }

    if (!balancingAuthority) {
        return <div>Balancing authority not found</div>;
    }

    return (
        <div style={{
            backgroundColor: '#f8fafc',
            minHeight: '100vh',          // Use minHeight instead of height
            height: 'auto',              // Allow content to expand
            overflow: 'auto',            // Force scrolling capability
            overflowY: 'scroll',         // Ensure vertical scrolling
            position: 'relative',        // Override any absolute positioning
            width: '100%',
            maxWidth: '100%',
            boxSizing: 'border-box'
        }}>
            {/* Header with back button */}
            <div style={{
                padding: '20px',
                backgroundColor: 'white',
                borderBottom: '1px solid #e2e8f0',
                display: 'flex',
                alignItems: 'center',
                position: 'relative',       // Ensure proper positioning
                zIndex: 1
            }}>
                <button
                    onClick={handleBackClick}
                    style={{
                        display: 'flex',
                        alignItems: 'center',
                        gap: '8px',
                        padding: '8px 16px',
                        backgroundColor: 'transparent',
                        border: '1px solid #cbd5e0',
                        borderRadius: '8px',
                        cursor: 'pointer',
                        fontSize: '14px',
                        color: '#4a5568'
                    }}
                >
                    <ArrowBackIcon style={{ fontSize: 16 }} />
                    {getBackButtonText()}
                </button>
            </div>

            {/* Section 1: Map and Info */}
            <div style={{ padding: '40px', position: 'relative' }}>
                <div style={{
                    display: 'grid',
                    gridTemplateColumns: '3fr 1fr', // 3:1 ratio
                    gap: '40px',
                    marginBottom: '60px'
                }}>
                    {/* Left: Large Map (3/4 of space) */}
                    <div style={{
                        backgroundColor: 'white',
                        borderRadius: '12px',
                        padding: '20px',
                        boxShadow: '0 4px 6px rgba(0, 0, 0, 0.1)',
                        height: '500px' // Increased height
                    }}>
                        <h3 style={{ margin: '0 0 20px 0', fontSize: '1.5rem', fontWeight: 600 }}>
                            {balancingAuthority.BA_Abrev} Territory
                        </h3>
                        <div style={{ height: '450px', borderRadius: '8px', overflow: 'hidden' }}>
                            {renderMap()}
                        </div>
                    </div>

                    {/* Right: Info (1/4 of space) */}
                    <div style={{
                        backgroundColor: 'white',
                        borderRadius: '12px',
                        padding: '30px',
                        boxShadow: '0 4px 6px rgba(0, 0, 0, 0.1)',
                        display: 'flex',
                        flexDirection: 'column'
                    }}>
                        <h1 style={{ margin: '0 0 8px 0', fontSize: '1.75rem', fontWeight: 700, color: '#1a202c' }}>
                            {balancingAuthority.BA_Abrev}
                        </h1>
                        <h2 style={{ margin: '0 0 30px 0', fontSize: '1rem', fontWeight: 500, color: '#4a5568', lineHeight: '1.4' }}>
                            {balancingAuthority.BA_Name}
                        </h2>

                        {/* All metadata in column format */}
                        <div style={{ display: 'flex', flexDirection: 'column', gap: '12px' }}>
                            <p style={{ margin: '0', color: '#4a5568', fontSize: '0.95rem' }}>
                                <strong>FID:</strong> {balancingAuthority.FID}
                            </p>
                            <p style={{ margin: '0', color: '#4a5568', fontSize: '0.95rem' }}>
                                <strong>Area Number{balancingAuthority.Area_Numbers && balancingAuthority.Area_Numbers.length > 1 ? 's' : ''}:</strong> {
                                    balancingAuthority.Area_Numbers && balancingAuthority.Area_Numbers.length > 0
                                        ? balancingAuthority.Area_Numbers.join(', ')
                                        : 'N/A'
                                }
                            </p>
                            <p style={{ margin: '0', color: '#4a5568', fontSize: '0.95rem' }}>
                                <strong>BA Abbreviation:</strong> {balancingAuthority.BA_Abrev}
                            </p>
                            <p style={{ margin: '0', color: '#4a5568', fontSize: '0.95rem' }}>
                                <strong>BA Name:</strong> {balancingAuthority.BA_Name}
                            </p>
                            <p style={{ margin: '0', color: '#4a5568', fontSize: '0.95rem' }}>
                                <strong>Shape Length:</strong> {balancingAuthority.Shape_Leng ? Number(balancingAuthority.Shape_Leng).toLocaleString() : 'N/A'} units
                            </p>
                            <p style={{ margin: '0', color: '#4a5568', fontSize: '0.95rem' }}>
                                <strong>Shape Area:</strong> {balancingAuthority.Shape__Area ? Number(balancingAuthority.Shape__Area).toLocaleString() : 'N/A'} units
                            </p>
                            <p style={{ margin: '0', color: '#4a5568', fontSize: '0.95rem' }}>
                                <strong>Shape Length (Alt):</strong> {balancingAuthority.Shape__Length ? Number(balancingAuthority.Shape__Length).toLocaleString() : 'N/A'} units
                            </p>
                            <p style={{ margin: '0', color: '#4a5568', fontSize: '0.85rem' }}>
                                <strong>Global ID:</strong> {balancingAuthority.GlobalID ? balancingAuthority.GlobalID : 'N/A'}
                            </p>
                        </div>
                    </div>
                </div>

                {/* Sections for each area number */}
                {balancingAuthority.Area_Numbers && balancingAuthority.Area_Numbers.map(areaNumber => (
                    <div key={areaNumber} style={{ marginBottom: '60px' }}>
                        <h2 style={{
                            fontSize: '1.75rem',
                            fontWeight: 700,
                            color: '#1a202c',
                            marginBottom: '30px',
                            paddingBottom: '10px',
                            borderBottom: '2px solid #e2e8f0'
                        }}>
                            Area {areaNumber} - Energy Data
                        </h2>

                        <div style={{ display: 'flex', flexDirection: 'column', gap: '30px' }}>
                            {/* Generation by Source */}
                            <div style={{
                                backgroundColor: 'white',
                                borderRadius: '12px',
                                padding: '25px',
                                boxShadow: '0 4px 6px rgba(0, 0, 0, 0.1)',
                                width: '100%'
                            }}>
                                <h3 style={{ margin: '0 0 20px 0', fontSize: '1.25rem', fontWeight: 600 }}>
                                    Electricity Generation by Energy Source
                                </h3>
                                {renderGenerationChart(areaNumber)}
                            </div>

                            {/* Demand vs Generation */}
                            <div style={{
                                backgroundColor: 'white',
                                borderRadius: '12px',
                                padding: '25px',
                                boxShadow: '0 4px 6px rgba(0, 0, 0, 0.1)',
                                width: '100%'
                            }}>
                                <h3 style={{ margin: '0 0 20px 0', fontSize: '1.25rem', fontWeight: 600 }}>
                                    Electricity Demand vs Total Generation
                                </h3>
                                {renderDemandChart(areaNumber)}
                            </div>

                            {/* Net Interchange */}
                            <div style={{
                                backgroundColor: 'white',
                                borderRadius: '12px',
                                padding: '25px',
                                boxShadow: '0 4px 6px rgba(0, 0, 0, 0.1)',
                                width: '100%'
                            }}>
                                <h3 style={{ margin: '0 0 20px 0', fontSize: '1.25rem', fontWeight: 600 }}>
                                    Net Electricity Interchange
                                </h3>
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
