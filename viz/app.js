import React, { useRef, useState, useCallback, useEffect, useReducer } from 'react';
import { createRoot } from "react-dom/client";
import { BrowserRouter as Router, Routes, Route, Navigate, useNavigate, useLocation } from 'react-router-dom';
import { StaticMap, Popup, Marker, _MapContext as MapContext, FullscreenControl, NavigationControl } from 'react-map-gl';
import { WebMercatorViewport } from '@deck.gl/core';
import DeckGL from '@deck.gl/react';
// import FlowMapLayer from '@flowmap.gl/core'
// import { FlowmapLayer } from '@flowmap.gl/layers'
import { GeoJsonLayer, ColumnLayer, PolygonLayer } from '@deck.gl/layers';
import { DataFilterExtension } from '@deck.gl/extensions';
import Checkbox from '@mui/material/Checkbox';
import FormControlLabel from '@mui/material/FormControlLabel';
import FormGroup from '@mui/material/FormGroup';
import { Typography } from '@mui/material';
import HomeOutlinedIcon from '@mui/icons-material/HomeOutlined';
import ThreeSixtyOutlinedIcon from '@mui/icons-material/ThreeSixtyOutlined';
import Slider from '@mui/material/Slider';
import Box from '@mui/material/Box';

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
import { AuthProvider, ProtectedRoute, Header, AdminDashboard } from './components/common';
import { ManishProject } from './components/manish';
import { AminProject, AminDetailPage } from './components/amin';

import 'core-js/actual/structured-clone';

ChartJS.register(RadialLinearScale, ArcElement, Tooltip, Legend);

// Transition interpolators for animation
const transitionLinearInterpolator = new LinearInterpolator(['bearing']);
const transitionFlyToInterpolator = new FlyToInterpolator(['zoom']);

// Get case data
var mod_casedata = require('./module_casedata.js');
var casedata = {};
casedata = mod_casedata.get_casedata();

// Source data GeoJSON
const geodata = casedata['geojsondata']

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


var data = ExtractFirstTimeSlice(geodata);

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


function MainApp({ refdata = data, refflowdata = flowdata, ggdata = geodata, mapStyle = MAP_STYLE }) {

  // Deck reference pointer
  const deckRef = useRef(null);

  // Navigation hook for routing
  const navigate = useNavigate();
  const location = useLocation();


  const [data, setData] = useState(refdata);

  //chat output message
  const [ouputMes, setOutputMes] = useState("Welcome to Westmap.");

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

      // if (properties.Shape_Leng) {
      //   infoLines.push(`Shape Length: ${Number(properties.Shape_Leng).toLocaleString()} units`);
      // }
      // if (properties.Shape__Area) {
      //   infoLines.push(`Shape Area: ${Number(properties.Shape__Area).toLocaleString()} units`);
      // }
      // if (properties.GlobalID) {
      //   infoLines.push(`Global ID: ${properties.GlobalID}`);
      // }

      popup.info = infoLines.join('\n');
      setShowPopup(showPopup => ({ ...showPopup, ...popup }));
    }
  });

  const GoHome = useCallback(() => {
    if (layers[0].context == null) return;
    var { viewport } = layers[0].context;
    const { longitude, latitude, zoom } = viewport.fitBounds(bounds);

    setInitialViewState(viewState => ({
      ...INITIAL_VIEW_STATE,
      longitude: longitude,
      latitude: latitude,
      zoom: zoom - 0.25,
      transitionInterpolator: transitionFlyToInterpolator,
      transitionDuration: 1000
    }))

    setShowPopup({ ...showPopup, display: false });
  });

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
    const apiPath = 'http://localhost:5000/data'; // for development
    // const apiPath = '/api/data'; // for production

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
    if (ouputMes === "Welcome to Westmap." || ouputMes === '') return;
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

        // Setup WECC generation chart data
        const weccChartData = {
          labels: ['Hydro', 'Nuclear', 'Coal', 'Natural Gas', 'Geothermal', 'Biomass', 'Wind', 'Solar'],
          datasets: [{
            label: 'WECC Generation Capacity (MW)',
            data: [
              totalBySource.Hydro,
              totalBySource.Nuclear,
              totalBySource.Coal,
              totalBySource.Natural_Gas,
              totalBySource.Geothermal,
              totalBySource.Biomass,
              totalBySource.Wind,
              totalBySource.PV
            ],
            backgroundColor: [
              'rgba(28,163,236,0.8)', // Hydro - Blue
              'rgba(255,0,0,0.8)',    // Nuclear - Red
              'rgba(0,0,0,0.8)',      // Coal - Black
              'rgba(255,165,0,0.8)',  // Natural Gas - Orange
              'rgba(128,0,128,0.8)',  // Geothermal - Purple
              'rgba(0,128,0,0.8)',    // Biomass - Green
              'rgba(0,255,0,0.8)',    // Wind - Light Green
              'rgba(255,255,0,0.8)'   // Solar - Yellow
            ],
            borderWidth: 1
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
      data: data,
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
        getFillColor: [weccHoveredObject, weccClickedObject]
      }
    }),

    // WECC Generation Power Column Layer - Multiple bars per location
    new ColumnLayer({
      id: 'WeccGenColumnLayer',
      data: weccColumnData,
      diskResolution: 50,
      radius: 5000, // Same as existing generation
      elevationScale: 50, // Same as existing generation  
      pickable: weccGenLayerActive,
      visible: weccGenLayerActive,
      getPosition: d => d.coordinates,
      getFillColor: fillWeccColumnColor, // Use our custom color function
      getElevation: d => d.Pg * 5, // Same pattern as existing generation
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
      <Header />
      <DeckGL
        ref={deckRef}
        layers={layers}
        initialViewState={initialViewState}
        controller={true}
        getTooltip={({ object }) => object && object.NAME}
        ContextProvider={MapContext.Provider}
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

        <div style={{ position: "absolute", top: 100, left: 0, "width": 30, background: "#fff", color: " #6b6b76", zIndex: 1000 }}>
          <HomeOutlinedIcon fontSize="medium" onClick={GoHome}></HomeOutlinedIcon>
          <br></br>
          {<ThreeSixtyOutlinedIcon fontSize="large" onClick={rotateCamera}>Rotate</ThreeSixtyOutlinedIcon>}
          <br></br>
        </div>


        {/*<div><NavigationControl position="top-left"/></div>
      <FullscreenControl/>*/}


        {
          showPopup.display && (
            <Popup style={{
              zIndex: 3,
              background: "white",
              opacity: 1,
              fontSize: "11px",
              cursor: showPopup.type === 'wecc' ? 'pointer' : 'default'
            }}
              longitude={initialViewState.longitude}
              latitude={initialViewState.latitude}
              anchor="bottom"
              offset={-100}
              onClose={() => setShowPopup({ ...showPopup, display: false })}>
              <div
                onClick={showPopup.type === 'wecc' ? handlePopupClick : undefined}
                style={{
                  padding: showPopup.type === 'wecc' ? '4px' : '0',
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
                <h3>{showPopup.name}</h3>
                <h4 style={{ whiteSpace: 'pre-line' }}>{showPopup.info}</h4>
              </div>
            </Popup>
          )
        }

      </DeckGL>


      <Widget
        handleNewUserMessage={handleUserInput}
        title="Westmap Chatbot"
        subtitle="What do you want to know about this power grid network?"
      />

      <div style={{ position: "absolute", top: 0, right: 0, "width": 250, background: "#fff", padding: "12px 12px", color: " #6b6b76", zIndex: 1000 }}>

        <div style={{ width: 300 }}>
          <Accordion defaultExpanded={true}>
            <AccordionSummary style={{ height: "20px", minHeight: "30px", paddingRight: "40px", paddingLeft: "0px" }}
              expandIcon={<ArrowDropDownIcon />}>
              <Typography> <Checkbox checked={netlayeractive} style={{ color: "primary" }} onChange={handleNetLayerChange} />Net</Typography>
            </AccordionSummary>
            <AccordionDetails>
              <Typography component="div">
                {netlayeractive &&
                  (
                    <div style={{ paddingRight: "40px" }}>

                      {/* <text> Voltage Level</text>
                      <br></br> */}
                      <Slider

                        value={netfiltervalue}
                        valueLabelDisplay="auto"
                        onChange={handleNetRangeFilterChange}
                        getAriaValueText={valuetext}
                        step={100}
                        min={0}
                        max={800}
                      >
                      </Slider>
                    </div>)
                }

                {netlayeractive && (<div style={{ paddingRight: "40px" }}><Multiselect
                  defaultValue={busNameSelectItems}
                  data={busNameItems}
                  placeholder={'Search for buses'}
                  onChange={handleBusMultiselect}
                /></div>)}

              </Typography>
            </AccordionDetails>
          </Accordion>
          <Accordion defaultExpanded={true}>
            <AccordionSummary style={{ height: "20px", minHeight: "30px", paddingRight: "40px", paddingLeft: "0px" }}
              expandIcon={<ArrowDropDownIcon />}>
              <Typography> <Checkbox checked={flowlayeractive} style={{ color: "primary" }} onChange={handleFlowLayerChange} />Flow
              </Typography>
            </AccordionSummary>
            <AccordionDetails>
              <Typography component="div">
                {flowlayeractive && (
                  <div style={{ paddingRight: "40px" }}>
                    <Slider
                      style={{ padding: 2 }}
                      value={flowfiltervalue}
                      valueLabelDisplay="auto"
                      onChange={handleFlowRangeFilterChange}
                      getAriaValueText={valuetext}
                      step={10}
                      min={0}
                      max={120}
                    >
                    </Slider></div>)
                }

                {flowlayeractive && (<div style={{ paddingRight: "40px" }}><Multiselect
                  defaultValue={lineNameSelectItems}
                  data={lineNameItems}
                  placeholder={'Search for transmission lines'}
                  onChange={handleLineMultiselect}
                /></div>)}

              </Typography>
            </AccordionDetails>
          </Accordion>

          <Accordion style={{ paddingBottom: "10px" }} defaultExpanded={false}>
            <AccordionSummary style={{ height: "20px", minHeight: "30px", paddingRight: "40px", paddingLeft: "0px" }}
              expandIcon={<ArrowDropDownIcon />}>
              <Typography>
                <Checkbox checked={genlayeractive} style={{ color: "primary" }} onChange={handleGenLayerChange} />Generation Power


              </Typography>
            </AccordionSummary>
            <AccordionDetails>
              <Typography component="div">
                {genlayeractive &&
                  (<div style={{ paddingRight: "40px" }}>
                    <Checkbox checked={genlayercapactive} style={{ color: "primary" }} onChange={handleGenLayerCapChange} />Generation Capacity
                  </div>)
                }

                {genlayeractive &&
                  (<div style={{ paddingRight: "40px" }}>
                    <Slider
                      style={{ padding: 2 }}
                      value={genfiltervalue}
                      valueLabelDisplay="auto"
                      onChange={handleGenRangeFilterChange}
                      getAriaValueText={valuetext}
                      step={100}
                      min={gendata.minPg}
                      max={gendata.maxPg + 10}
                    >
                    </Slider></div>)
                }


                {genlayeractive && (<div style={{ paddingRight: "40px" }}><Multiselect
                  defaultValue={nameSelectItems}
                  data={nameItems}
                  placeholder={'Search for generations'}
                  onChange={handleGenMultiselect}
                /></div>)}

                {
                  genlayeractive && (
                    <div style={{ width: 300, height: 300, transform: "translate(-1.5vw, 10px)" }}>
                      <Doughnut data={chartdata}
                        options={{
                          plugins: {
                            legend: {
                              onClick: handleDoughnutClick
                            }
                          }
                        }}
                      />
                    </div>
                  )}


              </Typography>
            </AccordionDetails>
          </Accordion>

          <Accordion defaultExpanded={false}>
            <AccordionSummary style={{ height: "20px", minHeight: "30px", paddingRight: "40px", paddingLeft: "0px" }}
              expandIcon={<ArrowDropDownIcon />}>
              <Typography>
                <Checkbox checked={loadlayeractive} style={{ color: "primary" }} onChange={handleLoadLayerChange} />Load loss
              </Typography>
            </AccordionSummary>
            <AccordionDetails>
              <Typography component="div">
                {loadlayeractive && (<div style={{ paddingRight: "40px" }}>
                  <Slider
                    style={{ padding: 2 }}
                    value={loadfiltervalue}
                    valueLabelDisplay="auto"
                    onChange={handleLoadRangeFilterChange}
                    getAriaValueText={valuetext}
                    step={100}
                    min={0}
                    max={countyloaddata.maxPd + 10}
                  >
                  </Slider></div>)
                }

                {loadlayeractive && (<div style={{ paddingRight: "40px" }}><Multiselect
                  defaultValue={countyNameSelectItems}
                  data={countyNameItems}
                  placeholder={'Search for counties'}
                  onChange={handleCountyMultiselect}
                /></div>)}


              </Typography>
            </AccordionDetails>
          </Accordion>


          <Accordion defaultExpanded={false}>
            <AccordionSummary style={{ height: "20px", minHeight: "30px", paddingRight: "40px", paddingLeft: "0px" }}
              expandIcon={<ArrowDropDownIcon />}>
              <Typography>
                <Checkbox checked={voltagelayeractive} style={{ color: "primary" }} onChange={handleVoltageLayerChange} />Voltage
              </Typography>
            </AccordionSummary>
            <AccordionDetails>
              <Typography component="div">
                {voltagelayeractive && (<div style={{ paddingRight: "40px" }}>
                  <Slider
                    style={{ padding: 2 }}
                    value={voltagefiltervalue}
                    valueLabelDisplay="auto"
                    onChange={handleVoltageRangeFilterChange}
                    getAriaValueText={valuetext}
                    step={0.01}
                    min={0.89}
                    max={1.11}
                  >
                  </Slider></div>)
                }

                {voltagelayeractive && (<div style={{ paddingRight: "40px" }}><Multiselect
                  defaultValue={countyNameSelectItems}
                  data={countyNameItems}
                  placeholder={'Search for counties'}
                  onChange={handleCountyMultiselect}
                /></div>)}


              </Typography>
            </AccordionDetails>
          </Accordion>

          <Accordion defaultExpanded={false}>
            <AccordionSummary style={{ height: "20px", minHeight: "30px", paddingRight: "40px", paddingLeft: "0px" }}
              expandIcon={<ArrowDropDownIcon />}>
              <Typography>
                <Checkbox checked={arealayeractive} style={{ color: "primary" }} onChange={handleAreaLayerChange} />Show Areas
              </Typography>
            </AccordionSummary>
            <AccordionDetails>
              <Typography component="div">
                {arealayeractive && (<div style={{ paddingRight: "40px" }}><Multiselect
                  defaultValue={areaNameSelectItems}
                  data={areaNameItems}
                  placeholder={'Search for areas'}
                  onChange={handleAreaMultiselect}
                /></div>)}
              </Typography>
            </AccordionDetails>
          </Accordion>

          <Accordion defaultExpanded={false}>
            <AccordionSummary style={{ height: "20px", minHeight: "30px", paddingRight: "40px", paddingLeft: "0px" }}
              expandIcon={<ArrowDropDownIcon />}>
              <Typography>
                Map Background
              </Typography>
            </AccordionSummary>
            <AccordionDetails>
              <Typography component="div">
                <div style={{ paddingRight: "40px" }}>
                  <select
                    value={mapStyleSelection}
                    onChange={(e) => setMapStyle(e.target.value)}
                    style={{ width: '100%', padding: '5px', marginBottom: '10px' }}
                  >
                    <option value="osm">OpenStreetMap (Free)</option>
                    <option value="satellite">Satellite View (Free)</option>
                    <option value="terrain">Terrain (Free)</option>
                    <option value="pos">Positron Light (Free)</option>
                    <option value="pos_no_label">Positron No Labels (Free)</option>
                    <option value="dark">Dark Matter (Free)</option>
                    <option value="none">No Background</option>
                  </select>
                </div>
              </Typography>
            </AccordionDetails>
          </Accordion>

          <Accordion defaultExpanded={false}>
            <AccordionSummary style={{ height: "20px", minHeight: "30px", paddingRight: "40px", paddingLeft: "0px" }}
              expandIcon={<ArrowDropDownIcon />}>
              <Typography>
                <Checkbox checked={zonelayeractive} style={{ color: "primary" }} onChange={handleZoneLayerChange} />Show Zones
              </Typography>
            </AccordionSummary>
            <AccordionDetails>
              <Typography component="div">
                {zonelayeractive && (<div style={{ paddingRight: "40px" }}><Multiselect
                  defaultValue={zoneNameSelectItems}
                  data={zoneNameItems}
                  placeholder={'Search for zones'}
                  onChange={handleZoneMultiselect}
                /></div>)}
              </Typography>
            </AccordionDetails>
          </Accordion>

          <Accordion defaultExpanded={false}>
            <AccordionSummary style={{ height: "20px", minHeight: "30px", paddingRight: "40px", paddingLeft: "0px" }}
              expandIcon={<ArrowDropDownIcon />}>
              <Typography>
                <Checkbox checked={wecclayeractive} style={{ color: "primary" }} onChange={handleWeccLayerChange} />Show WECC Data
              </Typography>
            </AccordionSummary>
            <AccordionDetails>
              <Typography component="div">
                {wecclayeractive && weccGeojsonData && (
                  <div style={{ paddingRight: "40px", fontSize: "12px", color: "#666" }}>
                    WECC Balancing Authorities overlay showing {weccGeojsonData.features?.length || 0} regions.
                    Click on any region for details.
                  </div>
                )}
              </Typography>
            </AccordionDetails>
          </Accordion>

          <Accordion defaultExpanded={false}>
            <AccordionSummary style={{ height: "20px", minHeight: "30px", paddingRight: "40px", paddingLeft: "0px" }}
              expandIcon={<ArrowDropDownIcon />}>
              <Typography>
                <Checkbox checked={weccGenLayerActive} style={{ color: "primary" }} onChange={handleWeccGenLayerChange} />WECC Generation Power
              </Typography>
            </AccordionSummary>
            <AccordionDetails>
              <Typography component="div">
                {weccGenLayerActive && weccGenData && (
                  <div style={{ paddingRight: "40px" }}>
                    <Slider
                      style={{ padding: 2 }}
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
                  <div style={{ paddingRight: "40px" }}>
                    <Multiselect
                      defaultValue={weccGenSelectItems}
                      data={weccGenNameItems}
                      placeholder={'Search based on Area Name'}
                      onChange={handleWeccGenMultiselect}
                    />
                  </div>
                )}

                {weccGenLayerActive && weccGenChartData && (
                  <div style={{ width: 280, height: 250, transform: "translate(-1.2vw, 10px)" }}>
                    <Doughnut data={weccGenChartData}
                      options={{
                        maintainAspectRatio: false,
                        plugins: {
                          legend: {
                            position: 'bottom',
                            onClick: handleWeccDoughnutClick,
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
                            text: 'WECC Generation Mix (MW)',
                            font: {
                              size: 12
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
} const rootElement = document.getElementById("root");

createRoot(rootElement).render(<App />)
