import React, {useRef, useState, useCallback,useEffect} from 'react';
import { createRoot } from "react-dom/client";
import {StaticMap, Popup,Marker,_MapContext as MapContext,FullscreenControl,NavigationControl} from 'react-map-gl';
import DeckGL from '@deck.gl/react';
import {GeoJsonLayer, ColumnLayer, PolygonLayer} from '@deck.gl/layers';
import {DataFilterExtension} from '@deck.gl/extensions';
import Checkbox from '@mui/material/Checkbox';
import FormControlLabel from '@mui/material/FormControlLabel';
import FormGroup from '@mui/material/FormGroup';
import { Typography } from '@mui/material';
import HomeOutlinedIcon from '@mui/icons-material/HomeOutlined';
import ThreeSixtyOutlinedIcon from '@mui/icons-material/ThreeSixtyOutlined';
import Slider from '@mui/material/Slider';
import Box from '@mui/material/Box';
import {
  Chart as ChartJS,
  RadialLinearScale,
  ArcElement,
  Tooltip,
  Legend,
} from 'chart.js';
import { PolarArea, Doughnut } from 'react-chartjs-2';

					//import casedata from './data/case_ACTIVSg200.json';
//import casedata500 from './data/case_ACTIVSg500.json';
					//import casedata9 from './data/case9.json'
//import casedata200 from './data/case_ACTIVSg200.json'
					//import casedata2000 from './data/case2000_opflow.json'
//import casedata2k from './data/case_ACTIVSg2000.json';
//import casedata10k from './data/case_ACTIVSg10k.json';
//import casedata70k from './data/case_ACTIVSg70k.json';
import countydata from "./data/counties.json";

import {center, convex, bbox} from '@turf/turf';

import { LinearInterpolator, FlyToInterpolator} from 'deck.gl';
import { HeatmapLayer } from 'deck.gl';
import { InvertColorsOff, ShopTwoOutlined } from '@mui/icons-material';


//var casedata = {};
//casedata.geojsondata = {};
//casedata.geojsondata.type = "FeatureCollection";
//casedata.geojsondata.features = [...casedata10k.geojsondata.features,...casedata2k.geojsondata.features,...casedata70k.geojsondata.features];

var mod_casedata = require('./module_casedata.js');

var casedata = {};

casedata = mod_casedata.get_casedata();

ChartJS.register(RadialLinearScale, ArcElement, Tooltip, Legend);

// Transition interpolators for animation
const transitionLinearInterpolator = new LinearInterpolator(['bearing']);
const transitionFlyToInterpolator = new FlyToInterpolator(['zoom']);

// Source data GeoJSON
const geodata = casedata['geojsondata']

// Get county data with loads
function getCountyNodes(data)
{
    var ncounties = countydata.features.length;
    var countygeojson = {type:"FeatureCollection",features:[],maxPd:0.0};
    for(var j=0; j < ncounties; j++) {
    //    if(countydata.features[j].properties.STATE !== "13" & countydata.features[j].properties.STATE !== "25")
    //        continue;
        var polygon = countydata.features[j];
        var box = bbox(polygon);

        countydata.features[j].properties.Pd = 0.0;
        countydata.features[j].properties.Vm_avg = 0.0;
        countydata.features[j].properties.KVlevels = [];
        var countyhassubst = false;
        var nbuscounty = 0
        for(var i=0; i < data.features.length; i++) {
          if(data.features[i].geometry.type == 'Point') {
            var lng = data.features[i].geometry.coordinates[0];
            var lat = data.features[i].geometry.coordinates[1];

            if((box[0] <= lng) & (box[2] >= lng) & (box[1] <= lat) & (box[3] >= lat)) {
              countyhassubst = true;
              var subst = data.features[i].properties;
              var nbus = subst.nbus;

              countydata.features[j].properties.KVlevels = [...countydata.features[j].properties.KVlevels,...data.features[i].properties.KVlevels];
              for(var k=0; k < nbus; k++) {
                var bus = subst.bus[k];
                countydata.features[j].properties.Pd += bus.PD;
                countydata.features[j].properties.Vm_avg += bus.VM;
                nbuscounty++;

              }
            }
          }
        }
        if(countyhassubst) {
          countydata.features[j].properties.Vm_avg /= nbuscounty;
          if(countydata.features[j].properties.Pd > countygeojson.maxPd) countygeojson.maxPd = countydata.features[j].properties.Pd;
          countygeojson.features.push(countydata.features[j]);
        }
    }
    return {maxPd:countygeojson.maxPd,data:countygeojson}
}


// Extract data for first time-slice
function ExtractFirstTimeSlice(data)
{
    var features = data.features.filter(function(feature) {
	if(!("start" in feature.properties) || (feature.properties.start == 0)) return feature;
    })
    var gdata = {"type": 'FeatureCollection',"features":features};
    return gdata;
}

const MAP_STYLE = {pos_no_label: 'https://basemaps.cartocdn.com/gl/positron-nolabels-gl-style/style.json',
pos:'https://basemaps.cartocdn.com/gl/positron-gl-style/style.json',
dark:'https://basemaps.cartocdn.com/gl/dark-matter-gl-style/style.json'};


function LineColor(line)
{
  var loading = Math.abs(line.properties.PF/line.properties.RATE_A);
  var r = Math.min(255,255*loading);
  var g = 0;
  var b = Math.max(0,255*(1-loading));

  return [r,g,b];
}

function FillColor(subst)
{
  var Vm = subst.properties.bus[0].VM;
  var r = Math.min(255,255*(1.1-Vm)/0.2);
  var g = b;
  var b = Math.max(0,255*(Vm-0.9)/0.2);

  return [r,g,b];
}

// Get substation voltages
function getPoints(data)
{
  var Pointsi;
  var Points = [];
  var elev;
  var i,k=0;
  for(i=0; i < data.features.length; i++) {
    if(data.features[i].geometry.type == 'Point') {
      Pointsi = {coordinates: data.features[i].geometry.coordinates, value: data.features[i].properties.bus[0].VM};
      Points.push(Pointsi);
    }
  }
  return Points;
}

// Get generation
function getGeneration(data)
{
  var Geni;
  var Gens = [];
  var minPg=1000.0;
  var maxPg=0.0;
  var elev;
  var i,j,k;
  var Pgcoal=0.0,Pghydro = 0.0,Pgnuclear = 0.0,Pgng = 0.0,Pgsolar = 0.0,Pgwind = 0.0,Pgother = 0.0;
  var Pgcoalcap = 0.0,Pghydrocap = 0.0,Pgnuclearcap = 0.0,Pgngcap = 0.0,Pgsolarcap = 0.0,Pgwindcap = 0.0,Pgothercap = 0.0;
  for(i=0; i < data.features.length; i++) {
    if(data.features[i].geometry.type == 'Point') {
      var subst = data.features[i].properties;
      var nbus = subst.nbus;
      var Pg = 0.0;
      var Pcap = 0.0;
      var gen_fuel;
      var ngen = 0;
      var KV = []
      for(j=0; j < nbus; j++) {
        var bus = subst.bus[j];
        KV.push(bus.BASE_KV);
        for(k=0; k < bus.ngen; k++) {
          var gen = bus.gen[k];
          Pg += gen.GEN_STATUS*gen.PG;
          Pcap += gen.GEN_STATUS*gen.PMAX;
          gen_fuel = gen.GEN_FUEL.toLowerCase();

          if(gen_fuel == 'wind') {
            Pgwind += gen.GEN_STATUS*gen.PG;
            Pgwindcap += gen.GEN_STATUS*gen.PMAX;
          }
          else if(gen_fuel == 'solar') {
            Pgsolar += gen.GEN_STATUS*gen.PG;
            Pgsolarcap += gen.GEN_STATUS*gen.PMAX;
          }
          else if(gen_fuel == 'coal') {
            Pgcoal += gen.GEN_STATUS*gen.PG;
            Pgcoalcap += gen.GEN_STATUS*gen.PMAX;
          } else if(gen_fuel == 'nuclear') {
            Pgnuclear += gen.GEN_STATUS*gen.PG;
            Pgnuclearcap += gen.GEN_STATUS*gen.PMAX;
          } else if(gen_fuel == 'hydro') {
            Pghydro += gen.GEN_STATUS*gen.PG;
            Pghydrocap += gen.GEN_STATUS*gen.PMAX;
          } else if(gen_fuel == 'ng') {
            Pgng += gen.GEN_STATUS*gen.PG;
            Pgngcap += gen.GEN_STATUS*gen.PMAX;
          } else {
            Pgother += gen.GEN_STATUS*gen.PG;
            Pgothercap += gen.GEN_STATUS*gen.PMAX;
          }
          ngen++;
        }
      }
      if(ngen) {
        var color;
        if(gen_fuel == 'wind') color = 'green';else if(gen_fuel == 'solar') color = 'yellow';
        else if(gen_fuel == 'coal') color = 'gray';
        else if(gen_fuel == 'nuclear') color = 'red';
        else if(gen_fuel == 'hydro') color = 'blue';
        else if(gen_fuel == 'ng') color = 'orange'
        else color = 'black';
        if(Pg <= minPg) minPg = Pg;
        if(Pg >= maxPg) maxPg = Pg;
        Geni = {coordinates: data.features[i].geometry.coordinates, Pg: Pg,Pcap:Pcap,KVlevels: KV, color:color};
        Gens.push(Geni);
      }
    }
  }

  return {minPg:minPg,maxPg:maxPg,Gens:Gens,Pgwind: Pgwind, Pgsolar: Pgsolar, Pgnuclear: Pgnuclear,Pghydro:Pghydro,Pgng:Pgng,Pgcoal:Pgcoal,Pgother:Pgother,Pgwindcap: Pgwindcap, Pgsolarcap: Pgsolarcap, Pgnuclearcap: Pgnuclearcap,Pghydrocap:Pghydrocap,Pgngcap:Pgngcap,Pgcoalcap:Pgcoalcap,Pgothercap:Pgothercap};
}

// Get Load
function getLoad(data)
{
  var Loadi;
  var Loads = [];
  var minPd = 1000.0;
  var maxPd = 0.0;
  var elev;
  var i,j;
  for(i=0; i < data.features.length; i++) {
    if(data.features[i].geometry.type == 'Point') {
      var subst = data.features[i].properties;
      var nbus = subst.nbus;
      var Pd = 0.0;
      var Qd = 0.0;
      for(j=0; j < nbus; j++) {
        var bus = subst.bus[j];
        Pd += bus.PD;
        Qd += bus.QD;
      }
      if(Pd > 0) {
        if(Pd <= minPd) minPd = Pd;
        if(Pd >= maxPd) maxPd = Pd;
        Loadi = {coordinates: data.features[i].geometry.coordinates, Pd: Pd};
        Loads.push(Loadi);
      }
    }
  }
  return {minPd:minPd,maxPd:maxPd,Loads:Loads};
}


function fillGenColumnColor(data)
{
  if(data.color == 'red') return [255,0,0];
  else if(data.color == 'green') return [0,255,0];
  else if(data.color == 'yellow') return [244,219,135];
  else if(data.color == 'gray') return [128,128,128];
  else if(data.color == 'blue') return [28,163,236];
  else if(data.color == 'orange') return [255,165,0];
  else if(data.color == 'black') return [0,0,0];

}

function fillGenColumnColorCap(data)
{
  var color;

  color = fillGenColumnColor(data);

  color = [...color,255*0.3];

  return color;
}


function getVoltageFillColor(data)
{
  var Vm = data.properties.Vm_avg;
  var r = Math.min(255,255*(1.1-Vm)/0.2);
  var b = 0;
  var g = Math.max(0,255*(Vm-0.9)/0.2);

  return [r,g,b];
}

function getContours()
{
  var Vmin = 0.9;
  var Vmax = 1.1;

  var contours = [
    {threshold: [0.9,0.98], color:[255,0,0]},
    {threshold: [0.98,1.02], color:[0,255,0]},
    {threshold: [1.02,1.1], color:[0,0,255]}
  ];
  return contours;
}

var data = ExtractFirstTimeSlice(geodata);

const Points = getPoints(data);
const Voltages = Points.map(d => d.value);

const Vcontour = getContours();

const gendata = getGeneration(data);
const generation = gendata.Gens;

const loaddata = getLoad(data);

const countyloaddata = getCountyNodes(data);

function LineWidth(line) {
  return line.properties.KV*3;
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
    bounds: [bboxArray[1], bboxArray[0],bboxArray[3],bboxArray[2]],
    fitbounds: true
};


export default function App({ggdata=eodatga,mapStyle = MAP_STYLE}) {

  // Deck reference pointer
  const deckRef = useRef(null);

  const [genfiltervalue,setGenFilterValue] = useState([gendata.minPg,gendata.maxPg]);

  const [netfiltervalue,setNetFilterValue] = useState([0,800]);

  const [loadfiltervalue,setLoadFilterValue] = useState([0,countyloaddata.maxPd]);

  const [voltagefiltervalue,setVoltageFilterValue] = useState([0.89,1.11]);

  // For zoom-in/out control
  const [initialViewState,setInitialViewState] = useState(INITIAL_VIEW_STATE);

  // For pop-up control
  const [showPopup, setShowPopup] = useState({display:false,info:'',name:''});

  var rotatestate = false;
  //const [rotatestate,setrotatestate] = useState(false);

  const rotateCamera = useCallback(() => {
    rotatestate = !rotatestate;
    if(rotatestate) {
      setInitialViewState(viewState => ({
        ...viewState,
        bearing: viewState.bearing - 180,
        transitionDuration: 20000,
        transitionInterpolator:transitionLinearInterpolator,
        onTransitionEnd: rotateCamera
      }))
    } else {
      setInitialViewState(viewState => ({
        ...viewState,
        onTransitiionEnd:null
      }))
//      GoHome();
      rotatestate = false;
    }

  }, []);

  const activatePopup = useCallback(() => {
    setShowPopup(showPopup => ({...showPopup,display:true}));

  },[]);

  const zoomToData = useCallback((info) => {
      var lat  = info.coordinate[1];
      var long = info.coordinate[0];

      setInitialViewState(viewState =>({
        ...viewState,
        latitude: lat,
        longitude: long,
        pitch: 50,
        traansitionInterpolator:transitionFlyToInterpolator,
        transitionDuration: 2000,
        zoom: 7.5,
        onTransitionEnd: activatePopup
      }))

      if(info.layer.id == "geojson") {
        if(info.object.geometry.type == "Point") {
          var popup = {};
          popup.name = info.object.properties.NAME
          popup.info = "Substation Info"
        } else {
          var popup = {};
          popup.name = info.object.properties.NAME
          popup.info = "Line Info"
        }
        setShowPopup(showPopup => ({...showPopup,...popup}));
      } else if(info.layer.id == "gen-column") {
        var popup = {};
        popup.name = "";
        popup.info = "Gen Info";

        setShowPopup(showPopup => ({...showPopup,...popup}));
      }
  });

  const zoomToCounty = useCallback((info) => {
  if(!info) return null;

  if(info.layer.id == 'PolygonLayer2') {
     var layer = info.layer;
     var {viewport} = layer.context;

     var cbounds = bbox(info.object);
     var c1 = [cbounds[0], cbounds[1]];
     var c2 = [cbounds[2], cbounds[3]];
     var countybounds = [c1, c2];
     const {longitude, latitude, zoom} = viewport.fitBounds(countybounds);

     setInitialViewState(viewState =>({
      ...viewState,
      latitude: latitude,
      longitude: longitude,
      pitch: 50,
      traansitionInterpolator:transitionFlyToInterpolator,
      transitionDuration: 5000,
      zoom: zoom-0.25,
      onTransitionEnd: activatePopup
    }))

    var popup = {display:false,name:'',info:''}; // Will be displayed after transition end only
    popup.name = info.object.properties.NAME;
    popup.info = "Load: "+info.object.properties.Pd.toFixed(2)+"MW";
    setShowPopup(showPopup => ({...showPopup,...popup}));


   }
});

  const GoHome = useCallback(() => {
    if(layers[0].context == null) return;
    var {viewport} = layers[0].context;
      const {longitude, latitude, zoom} = viewport.fitBounds(bounds);

    setInitialViewState(viewState =>({
      ...INITIAL_VIEW_STATE,
     longitude: longitude,
      latitude: latitude,
      zoom: zoom-0.25,
      transitionInterpolator:transitionFlyToInterpolator,
      transitionDuration: 2000
    }))

    setShowPopup({...showPopup,display:false});
  });

  const [netlayeractive, setNetLayerActive] = useState(true);

  const [loadlayeractive,setLoadLayerActive] = useState(false);
  const [genlayeractive,setGenLayerActive] = useState(false);
  const [voltagelayeractive,setVoltageLayerActive] = useState(false);

  const handleNetLayerChange = (event) => {
    setNetLayerActive(event.target.checked);
    setNetFilterValue([0,800]);
  };

  const handleLoadLayerChange = (event) => {
    setLoadLayerActive(event.target.checked);
    setLoadFilterValue([0,countyloaddata.maxPd]);

    event.target.checked && (setInitialViewState(viewState =>({
      ...viewState,
      pitch: 40,
      traansitionInterpolator:transitionFlyToInterpolator,
      transitionDuration: 2000,
    })))
  };

  const handleVoltageLayerChange = (event) => {
    setVoltageLayerActive(event.target.checked);
    setVoltageFilterValue([0.89,1.11]);

    event.target.checked && (setInitialViewState(viewState =>({
      ...viewState,
      pitch: 40,
      traansitionInterpolator:transitionFlyToInterpolator,
      transitionDuration: 2000,
    })))
  };

  const handleGenLayerChange = (event) => {
    setGenLayerActive(event.target.checked);
    setGenFilterValue([gendata.minPg,gendata.maxPg]);

    event.target.checked && (setInitialViewState(viewState =>({
      ...viewState,
      pitch: 40,
      traansitionInterpolator:transitionFlyToInterpolator,
      transitionDuration: 2000,
    })))
  };

  function getNetFilterValue(data) {
    if(!data) return 10000;
    if(data.geometry.type == 'Point') {
      for(var i=0; i < data.properties.KVlevels.length; i++) {
        var KV = data.properties.KVlevels[i];
        if(netfiltervalue[0] <= KV && KV <= netfiltervalue[1]) return KV;
      }
    } else {
      /* Line layer */
      return data.properties.KV;
    }
    return -1; // This is beyond the range so filter will filter out this data point.
  }

  function getGenFilterValue(data) {
    if(!data) return 10000;
    for(var i=0; i < data.KVlevels.length; i++) {
        var KV = data.KVlevels[i];
        if(netfiltervalue[0] <= KV && KV <= netfiltervalue[1]) {
          return data.Pg;
        }
    }
    return 10000;
  }

  function getLoadFilterValue(data) {
    if(!data) return -10000;
    for(var i=0; i < data.properties.KVlevels.length; i++) {
        var KV = data.properties.KVlevels[i];
        if(netfiltervalue[0] <= KV && KV <= netfiltervalue[1]) {
          return data.properties.Pd;
        }
    }
    return -10000;
  }

  function getVoltageFilterValue(data) {
    if(!data) return -10000;
    for(var i=0; i < data.properties.KVlevels.length; i++) {
        var KV = data.properties.KVlevels[i];
        if(netfiltervalue[0] <= KV && KV <= netfiltervalue[1]) {
          return data.properties.Vm_avg;
        }
    }
    return -10000;
  }

  const layers = [

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
      onClick:zoomToData,
      getFilterValue: getNetFilterValue,
      filterRange: netfiltervalue,

      extensions: [new DataFilterExtension({filtersize:1})]

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
      getElevation: d => d.Pg*5,
      onClick:zoomToData,

      getFilterValue: getGenFilterValue,
      filterRange: genfiltervalue,

      extensions: [new DataFilterExtension({filtersize:1})],

      updateTriggers: {
        getFilterValue: netfiltervalue
      }

    }),

    new ColumnLayer({
      id: 'gen-column-cap',
      data: generation,
      diskResolution: 50,
      radius: 5000,
      elevationScale: 50,
      pickable: false, //genlayeractive,
      visible: genlayeractive,
      getPosition: d => d.coordinates,
      getFillColor: fillGenColumnColorCap,
      getElevation: d => d.Pcap*5,
      onClick:zoomToData,

      getFilterValue: getGenFilterValue,
      filterRange: genfiltervalue,

      extensions: [new DataFilterExtension({filtersize:1})],

      updateTriggers: {
        getFilterValue: netfiltervalue
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
      getFillColor: getVoltageFillColor,
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
      id: 'PolygonLayer2',
      data:countyload,
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
      getLineColor: [80,80,80],
      getLineWidth: d => 1,
      opacity: 0.1,
      onClick: zoomToCounty,
      extensions: [new DataFilterExtension({filtersize:1})],
      getFilterValue: getVoltageFilterValue,
      filterRange: voltagefiltervalue,

      updateTriggers: {
        getFilterValue: netfiltervalue
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
          'green',
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


  function valuetext(value) {
    return `${value.toFixed(2)}`;
  }


  return (
    <>
      <DeckGL
        ref = {deckRef}
        layers={layers}
        initialViewState={initialViewState}
        controller={true}
        ContextProvider={MapContext.Provider}
      >


      <StaticMap reuseMaps
        mapStyle={mapStyle['pos']}
        preventStyleDiffing={true}
        initialViewState={INITIAL_VIEW_STATE}
      >
      </StaticMap>


      <FullscreenControl/>
      <br></br><br></br>
      <NavigationControl/>

      <div style={{position:"absolute", top:100, left:0, "width":30,background:"#fff",color:" #6b6b76",zIndex:1000}}>
      <HomeOutlinedIcon fontSize="medium" onClick={GoHome}></HomeOutlinedIcon>
      <br></br>
      {<ThreeSixtyOutlinedIcon fontSize="large" onClick={rotateCamera}>Rotate</ThreeSixtyOutlinedIcon>}
      <br></br>
      </div>


      {/*<div><NavigationControl position="top-left"/></div>
      <FullscreenControl/>*/}


      {
        showPopup.display && (
        <Popup style={{zIndex:3,background:"white",opacity:1}} longitude={initialViewState.longitude} latitude={initialViewState.latitude}
        anchor="bottom"
        offset={-100}
        onClose={() => setShowPopup({...showPopup,display:false})}>
        <h2>{showPopup.name}</h2><h3>{showPopup.info}</h3>
        </Popup>
        )
      }

      </DeckGL>

      {
        genlayeractive && (
        <div style={{width:300, height:300, position:"absolute",right:0,bottom:0,background:"white",padding:2,zIndex:1000}}>
        <Doughnut data={chartdata} />
        </div>
      )}

      <div style={{position:"absolute", top:0, right:0, "width":150,background:"#fff",padding:"12px 12px",color:" #6b6b76",zIndex:1000}}>
        <Checkbox checked={netlayeractive} style={{color:"primary"}} onChange={handleNetLayerChange} />Net
        <br></br>

        {netlayeractive && (
          <Slider
          style={{padding:2}}
          value={netfiltervalue}
          valueLabelDisplay="auto"
          onChange={handleNetRangeFilterChange}
          getAriaValueText={valuetext}
          step={100}
          min={0}
          max={800}
          >
          </Slider>)
        }

        <Checkbox checked={genlayeractive} style={{color:"primary"}} onChange={handleGenLayerChange} />Generation


        {genlayeractive && (
        <Slider
        style={{padding:2}}
        value={genfiltervalue}
        valueLabelDisplay="auto"
        onChange={handleGenRangeFilterChange}
        getAriaValueText={valuetext}
        step={100}
        min={gendata.minPg}
        max={gendata.maxPg+10}
        >
        </Slider>)
        }

       {/*
        <br></br>
        <Checkbox checked={loadlayeractive} style={ {color:"primary"}} onChange={handleLoadLayerChange} />Load

        {loadlayeractive && (
        <Slider
        style={{padding:2}}
        value={loadfiltervalue}
        valueLabelDisplay="auto"
        onChange={handleLoadRangeFilterChange}
        getAriaValueText={valuetext}
        step={100}
        min={0}
        max={countyloaddata.maxPd+10}
        >
        </Slider>)
        }
        */}

        <br></br>
        <Checkbox checked={voltagelayeractive} style={ {color:"primary"}} onChange={handleVoltageLayerChange} />Voltage

        {voltagelayeractive && (
        <Slider
        style={{padding:2}}
        value={voltagefiltervalue}
        valueLabelDisplay="auto"
        onChange={handleVoltageRangeFilterChange}
        getAriaValueText={valuetext}
        step={0.01}
        min={0.89}
        max={1.11}
        >
        </Slider>)
        }

      </div>

    </>



  );
}

const rootElement = document.getElementById("root");

createRoot(rootElement).render(<App />)
