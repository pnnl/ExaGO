import React, { useState, useEffect } from 'react';
import { useNavigate } from 'react-router-dom';
import { DeckGL } from '@deck.gl/react';
import { GeoJsonLayer } from '@deck.gl/layers';
import { MapView } from '@deck.gl/core';
import { StaticMap } from 'react-map-gl';
import OpenInNewIcon from '@mui/icons-material/OpenInNew';

// Mapbox token (same as main app)
const MAPBOX_ACCESS_TOKEN = 'pk.eyJ1IjoidXNtYXJ0LXdlc3RtYXAiLCJhIjoiY2tvazV6MzU2MDE4YjJ0bXd5ZDcwdm16ciJ9.q2BIGvGPAJjw1X9CdvyKSA';

// OpenStreetMap style (like in Manish's project)
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

const INITIAL_VIEW_STATE = {
    longitude: -116.5,
    latitude: 37.5,
    zoom: 5,
    minZoom: 3,
    maxZoom: 12,
    pitch: 0,
    bearing: 0
};

const AminMap = () => {
    const navigate = useNavigate();
    const [viewState, setViewState] = useState(INITIAL_VIEW_STATE);
    const [geojsonData, setGeojsonData] = useState(null);
    const [hoveredObject, setHoveredObject] = useState(null);
    const [clickedObject, setClickedObject] = useState(null);
    const [hoverInfo, setHoverInfo] = useState({});
    const [loading, setLoading] = useState(true);
    const [error, setError] = useState(null);
    const [areaMapping, setAreaMapping] = useState({});

    useEffect(() => {
        // Load both WECC GeoJSON and CSV data, then merge them
        const loadData = async () => {
            try {
                // Load GeoJSON for map visualization
                const geojsonResponse = await fetch('/amin_data/WECC_Balancing_Authorities_-2060174188301432986.geojson');
                if (!geojsonResponse.ok) {
                    throw new Error(`HTTP error loading GeoJSON! status: ${geojsonResponse.status}`);
                }
                const geojsonData = await geojsonResponse.json();

                // Load CSV for additional data fields
                const csvResponse = await fetch('/amin_data/WECC_Balancing_Authorities_5803277890210865950.csv');
                if (!csvResponse.ok) {
                    throw new Error(`HTTP error loading CSV! status: ${csvResponse.status}`);
                }
                const csvText = await csvResponse.text();

                // Load area mapping CSV
                const areaMappingResponse = await fetch('/amin_data/WECC_BA_Area_Mapping.csv');
                if (!areaMappingResponse.ok) {
                    throw new Error(`HTTP error loading area mapping CSV! status: ${areaMappingResponse.status}`);
                }
                const areaMappingText = await areaMappingResponse.text();

                // Parse CSV
                const csvLines = csvText.split('\n');
                const headers = csvLines[0].split(',');
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
                });

                setGeojsonData(geojsonData);
                setAreaMapping(areaMappingData);
                setLoading(false);
            } catch (err) {
                console.error('Error loading data:', err);
                setError(err.message);
                setLoading(false);
            }
        };

        loadData();
    }, []);

    const layers = [
        new GeoJsonLayer({
            id: 'wecc-balancing-authorities',
            data: geojsonData,
            pickable: true,
            stroked: true,
            filled: true,
            extruded: false,
            wireframe: false,
            getLineColor: [255, 255, 255, 255], // White border always visible
            getLineWidth: 4, // Thicker border for better visibility
            lineWidthMinPixels: 2, // Minimum width in pixels
            lineWidthMaxPixels: 10, // Maximum width in pixels
            getFillColor: d => {
                // Clicked region stays highlighted until another region is clicked
                if (d === clickedObject) {
                    return [30, 90, 150, 200]; // Darker blue for clicked/active state
                }
                // Hover effect (only if not clicked)
                if (d === hoveredObject && d !== clickedObject) {
                    return [70, 130, 180, 180]; // Slightly darker blue on hover
                }
                // Default blue fill with less opacity to show borders better
                return [100, 149, 237, 120]; // Cornflower blue with more transparency
            },
            onHover: (info) => {
                setHoveredObject(info.object);
                setHoverInfo(info);
            },
            onClick: (info) => {
                if (info.object) {
                    setClickedObject(info.object);
                }
            },
            updateTriggers: {
                getFillColor: [hoveredObject, clickedObject]
            }
        })
    ];

    const renderTooltip = () => {
        if (!hoveredObject || !hoverInfo.x || !hoverInfo.y) {
            return null;
        }

        const { BA_Abrev, BA_Name, Shape_Leng } = hoveredObject.properties;

        return (
            <div
                style={{
                    position: 'absolute',
                    zIndex: 1,
                    pointerEvents: 'none',
                    left: hoverInfo.x,
                    top: hoverInfo.y,
                    backgroundColor: 'rgba(0, 0, 0, 0.8)',
                    color: 'white',
                    padding: '8px',
                    borderRadius: '4px',
                    fontSize: '12px',
                    maxWidth: '300px'
                }}
            >
                <div><strong>{BA_Abrev}</strong></div>
                <div>{BA_Name}</div>
                <div>Perimeter: {Shape_Leng ? Number(Shape_Leng).toLocaleString() : 'N/A'} units</div>
            </div>
        );
    };

    const renderSidebar = () => {
        if (!clickedObject) {
            return (
                <div style={{
                    padding: '20px',
                    backgroundColor: 'rgba(255, 255, 255, 0.95)',
                    borderRadius: '12px',
                    boxShadow: '0 8px 32px rgba(0,0,0,0.1)',
                    backdropFilter: 'blur(10px)',
                    border: '1px solid rgba(255,255,255,0.2)'
                }}>
                    <h3 style={{ margin: '0 0 16px 0', fontSize: '1.25rem', fontWeight: 600, color: '#2c3e50' }}>
                        WECC Balancing Authorities
                    </h3>
                    <p style={{ margin: '0 0 16px 0', color: '#666', fontSize: '0.875rem' }}>
                        Click on a balancing authority region to see detailed information.
                    </p>
                    <p style={{ margin: '0', color: '#666', fontSize: '0.875rem' }}>
                        This map shows the Western Electricity Coordinating Council (WECC)
                        balancing authority areas. Each colored region represents a different
                        balancing authority responsible for maintaining electrical grid reliability
                        in their area.
                    </p>
                </div>
            );
        }

        const { FID, BA_Abrev, BA_Name, Shape_Leng, Shape__Area, Shape__Length, GlobalID, Area_Numbers } = clickedObject.properties;

        const handleCardClick = () => {
            // Navigate to the detail page using the FID
            navigate(`/amin/${FID}`);
        };

        return (
            <div
                style={{
                    position: 'relative',
                    padding: '20px',
                    backgroundColor: 'rgba(255, 255, 255, 0.95)',
                    borderRadius: '12px',
                    border: '1px solid rgba(255,255,255,0.2)',
                    cursor: 'pointer',
                    transition: 'all 0.2s ease'
                }}
                onClick={handleCardClick}
            >
                {/* External Link Icon in top right */}
                <div style={{
                    position: 'absolute',
                    top: '16px',
                    right: '16px',
                    opacity: 0.6,
                    transition: 'opacity 0.2s ease'
                }}>
                    <OpenInNewIcon style={{ fontSize: 18, color: '#666' }} />
                </div>

                <h3 style={{ margin: '0 0 8px 0', fontSize: '1.5rem', fontWeight: 600, color: '#2c3e50' }}>
                    {BA_Abrev}
                </h3>
                <p style={{ margin: '0 0 16px 0', fontSize: '1.1rem', fontWeight: 500, color: '#34495e' }}>
                    {BA_Name}
                </p>

                {/* All the detailed information */}
                <div style={{ borderTop: '1px solid #eee', paddingTop: '12px' }}>
                    <p style={{ margin: '0 0 8px 0', color: '#666', fontSize: '0.875rem' }}>
                        <strong>FID:</strong> {FID}
                    </p>
                    <p style={{ margin: '0 0 8px 0', color: '#666', fontSize: '0.875rem' }}>
                        <strong>Area Number{Area_Numbers && Array.isArray(Area_Numbers) && Area_Numbers.length > 1 ? 's' : ''}:</strong> {
                            Area_Numbers && Array.isArray(Area_Numbers) && Area_Numbers.length > 0
                                ? Area_Numbers.join(', ')
                                : 'N/A'
                        }
                    </p>
                    <p style={{ margin: '0 0 8px 0', color: '#666', fontSize: '0.875rem' }}>
                        <strong>BA Abbreviation:</strong> {BA_Abrev}
                    </p>
                    <p style={{ margin: '0 0 8px 0', color: '#666', fontSize: '0.875rem' }}>
                        <strong>BA Name:</strong> {BA_Name}
                    </p>
                    <p style={{ margin: '0 0 8px 0', color: '#666', fontSize: '0.875rem' }}>
                        <strong>Shape Length:</strong> {Shape_Leng ? Number(Shape_Leng).toLocaleString() : 'N/A'} units
                    </p>
                    <p style={{ margin: '0 0 8px 0', color: '#666', fontSize: '0.875rem' }}>
                        <strong>Shape Area:</strong> {Shape__Area ? Number(Shape__Area).toLocaleString() : 'N/A'} units
                    </p>
                    <p style={{ margin: '0 0 8px 0', color: '#666', fontSize: '0.875rem' }}>
                        <strong>Shape Length (Alt):</strong> {Shape__Length ? Number(Shape__Length).toLocaleString() : 'N/A'} units
                    </p>
                    <p style={{ margin: '0', color: '#666', fontSize: '0.75rem' }}>
                        <strong>Global ID:</strong> {GlobalID}
                    </p>
                </div>
            </div>
        );
    };

    if (loading) {
        return (
            <div style={{ display: 'flex', justifyContent: 'center', alignItems: 'center', height: '100vh' }}>
                <h2>Loading WECC Balancing Authorities...</h2>
            </div>
        );
    }

    if (error) {
        return (
            <div style={{ display: 'flex', justifyContent: 'center', alignItems: 'center', height: '100vh' }}>
                <h2 style={{ color: 'red' }}>
                    Error loading data: {error}
                </h2>
            </div>
        );
    }

    return (
        <div style={{ position: 'relative', height: '100vh' }}>
            {/* Full screen map */}
            <DeckGL
                viewState={viewState}
                onViewStateChange={({ viewState }) => setViewState(viewState)}
                controller={true}
                layers={layers}
                views={new MapView({ id: 'map' })}
            >
                <StaticMap
                    mapboxApiAccessToken={MAPBOX_ACCESS_TOKEN}
                    mapStyle={OSM_MAP_STYLE}
                />
            </DeckGL>
            {renderTooltip()}

            {/* Floating info block in top-right corner */}
            <div style={{
                position: 'absolute',
                top: '20px',
                right: '20px',
                width: '320px',
                maxHeight: '80vh',
                overflowY: 'auto',
                zIndex: 1000,
                pointerEvents: 'auto'
            }}>
                {renderSidebar()}
            </div>
        </div>
    );
};

export default AminMap;
