import React, { useRef, useState, useCallback, useEffect, useReducer } from 'react';
import { createRoot } from "react-dom/client";
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

const style='pos';
const MAP_STYLE = {
  pos_no_label: 'https://basemaps.cartocdn.com/gl/positron-nolabels-gl-style/style.json',
  pos: 'https://basemaps.cartocdn.com/gl/positron-gl-style/style.json',
  dark: 'https://basemaps.cartocdn.com/gl/dark-matter-gl-style/style.json',
  none:''
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


export default function App({ refdata = data, refflowdata = flowdata, ggdata = geodata, mapStyle = MAP_STYLE }) {

  // Deck reference pointer
  const deckRef = useRef(null);


  const [data, setData] = useState(refdata);

  //chat output message
  const [ouputMes, setOutputMes] = useState("Welcome to ChatGrid.");

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

  const [netfiltervalue, setNetFilterValue] = useState([0, 800]);

  const [flowfiltervalue, setFlowFilterValue] = useState([0, 120]);

  const [loadfiltervalue, setLoadFilterValue] = useState([0, countyloaddata.maxPd]);

  const [voltagefiltervalue, setVoltageFilterValue] = useState([0.89, 1.11]);

  // For zoom-in/out control
  const [initialViewState, setInitialViewState] = useState(INITIAL_VIEW_STATE);

  // For pop-up control
  const [showPopup, setShowPopup] = useState({ display: false, info: '', name: '' });

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
	    if(feature.properties.RATE_A == 0) {
		RATE_A = 10000;
	    } else {
		RATE_A = feature.properties.RATE_A;
	    }
	  var loading = Math.abs(feature.properties.PF / RATE_A)*100;
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
	    if(feature.properties.RATE_A == 0) {
		RATE_A = 10000;
	    } else {
		RATE_A = feature.properties.RATE_A;
	    }

	  var loading = Math.abs(feature.properties.PF / RATE_A)*100;
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
	    if(feature.properties.RATE_A == 0) {
		RATE_A = 10000;
	    } else {
		RATE_A = feature.properties.RATE_A;
	    }
	    var loading = Math.abs(feature.properties.PF/RATE_A)*100.0;

	    if(flowfiltervalue[0] <= loading && loading <= flowfiltervalue[1]) {
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
      traansitionInterpolator: transitionFlyToInterpolator,
      transitionDuration: 2000,
      zoom: 7.5,
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
      traansitionInterpolator: transitionFlyToInterpolator,
      transitionDuration: 2000,
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
	  var loading = Math.abs(info.object.properties.PF / info.object.properties.RATE_A)*100.0;
          popup.name = info.object.properties.NAME
          popup.info = "KV: "+info.object.properties.KV.toFixed(2)+"KV \nLoading: "+loading.toFixed(2)+ "%"; 
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
      traansitionInterpolator: transitionFlyToInterpolator,
      transitionDuration: 5000,
      zoom: 7.5,
      onTransitionEnd: activatePopup
    }))

  
  });

    const zoomToAreaName = useCallback((filterareas,minLng, minLat, maxLng, maxLat) => {
    var viewport = new WebMercatorViewport(INITIAL_VIEW_STATE);

    const { longitude, latitude, zoom } = viewport.fitBounds([[minLng, minLat], [maxLng, maxLat]]);

    setInitialViewState(viewState => ({
      ...viewState,
      latitude: latitude,
      longitude: longitude,
      pitch: 50,
      traansitionInterpolator: transitionFlyToInterpolator,
      transitionDuration: 5000,
      zoom: 7.5,
      onTransitionEnd: activatePopup
    }))

      var popup = { display: false, name: '', info: '' }; // Will be displayed after transition end only
      popup.name = "Area " + filterareas.properties.name;
//      popup.info = "Area: " + info.object.properties.name;
      setShowPopup(showPopup => ({ ...showPopup, ...popup }));
    
    });

    const zoomToZoneName = useCallback((filterzones,minLng, minLat, maxLng, maxLat) => {
    var viewport = new WebMercatorViewport(INITIAL_VIEW_STATE);

    const { longitude, latitude, zoom } = viewport.fitBounds([[minLng, minLat], [maxLng, maxLat]]);

    setInitialViewState(viewState => ({
      ...viewState,
      latitude: latitude,
      longitude: longitude,
      pitch: 50,
      traansitionInterpolator: transitionFlyToInterpolator,
      transitionDuration: 5000,
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
        traansitionInterpolator: transitionFlyToInterpolator,
        transitionDuration: 5000,
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
        traansitionInterpolator: transitionFlyToInterpolator,
        transitionDuration: 5000,
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
        traansitionInterpolator: transitionFlyToInterpolator,
        transitionDuration: 5000,
        zoom: zoom - 0.25,
        onTransitionEnd: activatePopup
      }))

      var popup = { display: false, name: '', info: '' }; // Will be displayed after transition end only
      popup.name = "Zone " + info.object.properties.name;
//      popup.info = "Zone: " + info.object.properties.name;
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
      transitionDuration: 2000
    }))

    setShowPopup({ ...showPopup, display: false });
  });

  const [netlayeractive, setNetLayerActive] = useState(true);
  const [flowlayeractive, setFlowLayerActive] = useState(true);

  const [loadlayeractive, setLoadLayerActive] = useState(false);
  const [genlayeractive, setGenLayerActive] = useState(false);
  const [genlayercapactive, setGenLayerCapActive] = useState(false);
  const [voltagelayeractive, setVoltageLayerActive] = useState(false);
  const [zonelayeractive, setZoneLayerActive] = useState(false);
  const [arealayeractive, setAreaLayerActive] = useState(false);
    
  const handleUserInput = (inputText) => {
    console.log(`New message incoming! ${inputText}`);
    // Now send the message to GPT and get response 
    toggleMsgLoader()
    const postData = {
      "inputText": inputText
    }
    try{
      fetch(`http://localhost:5000/data`, {
        "method": "POST",
        // headers: { 'Content-Type': 'application/json' },
        "body": JSON.stringify(postData),
      }).then((res) =>
        res.json().then((chatOutput) => {
          // Setting a data from api
          console.log(chatOutput);
          const outputText = chatOutput.text
          const chatList = chatOutput.result_list
          const keyList = Object.keys(chatList[0])
          if (chatList.length > 0) {
            //  only one is active between bus name selection, transmission line name selection at a time
  
            if ("generation name" in chatList[0]) {
              const genNameList = chatList.map(d => d["generation name"]);
              setGenLayerActive(true)
              setNameSelectItems(genNameList)
              setGenFilterValue([gendata.minPg, gendata.maxPg]);
  
              setInitialViewState(viewState => ({
                ...viewState,
                pitch: 40,
                traansitionInterpolator: transitionFlyToInterpolator,
                transitionDuration: 2000,
              }))
            }
            const containCapacity =  keyList.some(str => str.includes('capacity'))
            if ("generation name" in chatList[0] && containCapacity){
              const genNameList = chatList.map(d => d["generation name"]);
              setGenLayerActive(true)
              setGenLayerCapActive(true)
              setNameSelectItems(genNameList)
              setGenFilterValue([gendata.minPg, gendata.maxPg]);
  
              setInitialViewState(viewState => ({
                ...viewState,
                pitch: 40,
                traansitionInterpolator: transitionFlyToInterpolator,
                transitionDuration: 2000,
              }))
            }
            if ("line name" in chatList[0]) {
              const lineNameList = chatList.map(d => d["line name"]);
              setNetLayerActive(true)
              setFlowLayerActive(true)
              setBusNameSelectItems([])
              setLineNameSelectItems(lineNameList)
            }
            if ('bus name' in chatList[0]) {
              const busNameList = chatList.map(d => d["bus name"]);
              setNetLayerActive(true)
              setFlowLayerActive(true)
              setBusNameSelectItems(busNameList)
              setLineNameSelectItems([])
            }
  
  
          }
  
          setOutputMes(outputText)
  
        })
      );
    } catch (error){
      setOutputMes("Sorry I didn't find the answer to your question. Please try to rephrase it or provide more details.")
    }
   

  };

  useEffect(() => {

    addResponseMessage(`${ouputMes}`);
    if (ouputMes === "Welcome to ChatGrid." || ouputMes ==='') return;
    toggleMsgLoader(); // close loading 
  }, [ouputMes]);

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
      traansitionInterpolator: transitionFlyToInterpolator,
      transitionDuration: 2000,
    })))
  };

  const handleVoltageLayerChange = (event) => {
    setVoltageLayerActive(event.target.checked);
    setVoltageFilterValue([0.89, 1.11]);

    event.target.checked && (setInitialViewState(viewState => ({
      ...viewState,
      pitch: 40,
      traansitionInterpolator: transitionFlyToInterpolator,
      transitionDuration: 2000,
    })))
  };

    const handleAreaLayerChange = (event) => {
    setAreaLayerActive(event.target.checked);
//    setVoltageFilterValue([0.89, 1.11]);

    event.target.checked && (setInitialViewState(viewState => ({
      ...viewState,
      pitch: 40,
      traansitionInterpolator: transitionFlyToInterpolator,
      transitionDuration: 2000,
    })))
    };

    const handleZoneLayerChange = (event) => {
    setZoneLayerActive(event.target.checked);
//    setVoltageFilterValue([0.89, 1.11]);

    event.target.checked && (setInitialViewState(viewState => ({
      ...viewState,
      pitch: 40,
      traansitionInterpolator: transitionFlyToInterpolator,
      transitionDuration: 2000,
    })))
  };


  const handleGenLayerChange = (event) => {
    setGenLayerActive(event.target.checked);
    setGenFilterValue([gendata.minPg, gendata.maxPg]);

    event.target.checked && (setInitialViewState(viewState => ({
      ...viewState,
      pitch: 40,
      traansitionInterpolator: transitionFlyToInterpolator,
      transitionDuration: 2000,
    })))
  };

  const handleGenLayerCapChange = (event) => {
    setGenLayerCapActive(event.target.checked);
    setGenFilterValue([gendata.minPg, gendata.maxPg]);

    event.target.checked && (setInitialViewState(viewState => ({
      ...viewState,
      pitch: 40,
      traansitionInterpolator: transitionFlyToInterpolator,
      transitionDuration: 2000,
    })))
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
	  if(data.geometry.type == 'LineString') { /* Line layer */
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

  const layers = [

    // new FlowmapLayer({
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

        zoomToAreaName(filterareas[0],minLng, minLat, maxLng, maxLat)
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

        zoomToZoneName(filterzones[0],minLng, minLat, maxLng, maxLat)
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
      <DeckGL
        ref={deckRef}
        layers={layers}
        initialViewState={initialViewState}
        controller={true}
        getTooltip={({ object }) => object && object.NAME}
        ContextProvider={MapContext.Provider}
      >


        <StaticMap reuseMaps
          mapStyle={mapStyle[style]}
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
            <Popup style={{ zIndex: 3, background: "white", opacity: 1, fontSize: "11px" }} longitude={initialViewState.longitude} latitude={initialViewState.latitude}
              anchor="bottom"
              offset={-100}
              onClose={() => setShowPopup({ ...showPopup, display: false })}>
              <h2>{showPopup.name}</h2><h3>{showPopup.info}</h3>
            </Popup>
          )
        }

      </DeckGL>


      <Widget
        handleNewUserMessage={handleUserInput}
        title="ChatGrid"
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
                { genlayeractive &&
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

        </div>

        {/* <br></br> */}

        {/* <BrushingBarChart data={[12,23,345,45,66,78,800]} width ={170} height = {100} handleFilter = {handleNetBarFilterChange}/> */}

        {/* <br></br> */}


      </div>


    </>



  );
}

const rootElement = document.getElementById("root");

createRoot(rootElement).render(<App />)
