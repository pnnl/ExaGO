import React, {useRef, useState, useCallback} from 'react';
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
import casedata500 from './data/case_ACTIVSg500.json';
import casedata2k from './data/case_ACTIVSg2000.json';
import casedata10k from './data/case_ACTIVSg10k.json';
import casedata70k from './data/case_ACTIVSg70k.json';
import countydata from "./data/counties.json";

import {center, convex, bbox} from '@turf/turf';

import { LinearInterpolator, FlyToInterpolator} from 'deck.gl';
import { HeatmapLayer } from 'deck.gl';
import { InvertColorsOff, ShopTwoOutlined } from '@mui/icons-material';

ChartJS.register(RadialLinearScale, ArcElement, Tooltip, Legend);

var casedata = {};
casedata.geojsondata = {};
casedata.geojsondata.type = "FeatureCollection";
casedata.geojsondata.features = [...casedata2k.geojsondata.features];

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
	if(feature.properties.start == 0) return feature;
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

function LineWidth(line)
{
  return line.properties.KV*3;
  //return Math.abs(line.properties.PF/line.properties.RATE_A)*500;
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
      var gen_fuel;
      var ngen = 0;
      for(j=0; j < nbus; j++) {
        var bus = subst.bus[j];
        for(k=0; k < bus.ngen; k++) {
          var gen = bus.gen[k];
          Pg += gen.GEN_STATUS*gen.PG;
          gen_fuel = gen.GEN_FUEL;

          if(gen_fuel == 'wind') Pgwind += gen.GEN_STATUS*gen.PG;
          else if(gen_fuel == 'solar') Pgsolar += gen.GEN_STATUS*gen.PG;
          else if(gen_fuel == 'coal') Pgcoal += gen.GEN_STATUS*gen.PG;
          else if(gen_fuel == 'nuclear') Pgnuclear += gen.GEN_STATUS*gen.PG;
          else if(gen_fuel == 'hydro') Pghydro += gen.GEN_STATUS*gen.PG;
          else if(gen_fuel == 'ng') Pgng += gen.GEN_STATUS*gen.PG;
          else Pgother += gen.GEN_STATUS*gen.PG;

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
        Geni = {coordinates: data.features[i].geometry.coordinates, Pg: Pg,color:color};
        Gens.push(Geni);
      }
    }
  }
  
  return {minPg:minPg,maxPg:maxPg,Gens:Gens,Pgwind: Pgwind, Pgsolar: Pgsolar, Pgnuclear: Pgnuclear,Pghydro:Pghydro,
    Pgng:Pgng,Pgcoal:Pgcoal,Pgother:Pgother};
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

}

function fillVoltageColumnColor(data)
{
  var Vm = data.value;
  var r = Math.min(255,255*(1.1-Vm)/0.2);
  var g = 0;
  var b = Math.max(0,255*(Vm-0.9)/0.2);

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

export default function App({ggdata=geodata,mapStyle = MAP_STYLE}) {

  const deckRef = useRef(null);
  var data = ExtractFirstTimeSlice(ggdata);

  const Points = getPoints(data);

  const Voltages = Points.map(d => d.value);

  const Vcontour = getContours();

  const gendata = getGeneration(data);
  const generation = gendata.Gens;
  const [genfiltervalue,setGenFilterValue] = useState([gendata.minPg,gendata.maxPg]);

  const loaddata = getLoad(data);

  const loads = loaddata.Loads;
  const minPd = loaddata.minPd;
  const maxPd = loaddata.maxPd;

  const countyloaddata = getCountyNodes(data);
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
        bearing: viewState.bearing + 120,
        transitionDuration: 10000,
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
    console.log(showPopup);
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
      }
  });

  const zoomToCounty = useCallback((info) => {
  if(!info) return null;
   
  if(info.layer.id == 'PolygonLayer2') {
     var layer = info.layer;
     var {viewport} = layer.context;

     console.log(info);
     var cbounds = bbox(info.object);
     var c1 = [cbounds[0], cbounds[1]];
     var c2 = [cbounds[2], cbounds[3]];
     var countybounds = [c1, c2];
     console.log(countybounds);
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
    console.log(info.object)
    var popup = {display:false,name:'',info:''}; // Will be displayed after transition end only
    popup.name = info.object.properties.NAME;
    popup.info = "Load: "+info.object.properties.Pd.toFixed(2)+"MW";
    setShowPopup(showPopup => ({...showPopup,...popup}));
    console.log(showPopup);
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

  const [checked, setChecked] = useState(true);

  const [loadlayeractive,setLoadLayerActive] = useState(false);
  const [genlayeractive,setGenLayerActive] = useState(false);

  const handleChange = (event) => {
    setChecked(event.target.checked);
  };

  const handleLoadLayerChange = (event) => {
    setLoadLayerActive(event.target.checked);

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


  const layers = [
    
    new GeoJsonLayer({
      id: 'geojson', 
      data: data,
      stroked: false,
      filled: checked,
      extruded: true,
      pickable: checked,
      pointType: 'circle',
      lineWidthScale: 3,
      getFillColor: FillColor,
      getLineColor: LineColor,
      getPointRadius: 1000,
      getLineWidth: LineWidth,
      visible: checked,
      onClick:zoomToData
//      extensions: [new ClipExtension()],
 //     clipBounds: bounds,
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

      getFilterValue: d => d.Pg,
      filterRange: genfiltervalue,

      extensions: [new DataFilterExtension({filtersize:1})]

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
      getElevation: d => d.Pd*10,
      onClick:zoomToData
    }),
    */
    
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
      getElevation: d => d.properties.Pd*5.0,
      getFillColor: d => [255*d.properties.Pd/countymaxPd, 0, 0],
      getLineColor: [80,80,80],
      getLineWidth: d => 1,
      opacity: 0.1,
      onClick: zoomToCounty
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

  console.log(genmix)

  const chartdata = {
    labels: genmixlabels,
    datasets: [
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
      },
    ],
  };
  

  const handleGenRangeFilterChange = (event) => {
    setGenFilterValue(event.target.value);
  }

  function valuetext(value) {
    return `${value.toFixed(2)}`;
  }

  
  return (
    <><DeckGL 
      style={{zIndex: 1}}
      ref = {deckRef}
      layers={layers}
      initialViewState={initialViewState}
      controller={true}
//      getTooltip={Tooltip}
//      onClick={zoomToData}
//      onAfterRender={GoHome}
      ContextProvider={MapContext.Provider}
    >

    <StaticMap reuseMaps 
    style={{zIndex: 2}}
    mapStyle={mapStyle['pos']} 
    preventStyleDiffing={true} 
    initialViewState={INITIAL_VIEW_STATE}>

    {showPopup.display && (
      <Popup style={{zIndex:3,background:"white",opacity:1}} longitude={initialViewState.longitude} latitude={initialViewState.latitude}
        anchor="bottom"
        offset={-100}
        onClose={() => setShowPopup({...showPopup,display:false})}><h1>{showPopup.name}</h1><h2>{showPopup.info}</h2>
      </Popup>
      )
    }
    
    </StaticMap>

    {/*<div><NavigationControl position="top-left"/></div>
    <FullscreenControl/>*/}

    <FullscreenControl/>
    <br></br><br></br>
    <NavigationControl/>
    
    <div style={{position:"absolute", top:0, right:0, "width":120,background:"#fff",padding:"12px 24px",color:" #6b6b76"}}>

  
    <HomeOutlinedIcon color="primary" fontSize="large" onClick={GoHome}></HomeOutlinedIcon>
    <br></br>
  {<ThreeSixtyOutlinedIcon color="primary" fontSize="large" onClick={rotateCamera}></ThreeSixtyOutlinedIcon>}
<br></br>
  <Checkbox checked={checked} style={{color:"primary"}} onChange={handleChange} />Net
    <br></br>
    <Checkbox checked={genlayeractive} style={{color:"primary"}} onChange={handleGenLayerChange}/>Gen

    {/*
    {genlayeractive && (
    <Slider
      style={{padding:2}}
      value={genfiltervalue}
      valueLabelDisplay="auto"
      onChange={handleGenRangeFilterChange}
      getAriaValueText={valuetext}
      min={gendata.minPg}
      max={gendata.maxPg+10}
    >
    </Slider>)
    }
    */}
    
    <br></br>
    <Checkbox checked={loadlayeractive} style={{color:"primary"}} onChange={handleLoadLayerChange} />Load
    <br></br>
  
  </div>

{/*
  {genlayeractive && (
    <div style={{width:120, position:"absolute",right:0,bottom:"50%",background:"white",padding:2,"font-size":"18"}}>
    <div style={{background:"green",color:"white","font-size":"18"}}>Solar/Wind</div>
    <div style={{background:"gray",color:"white",padding:2,"font-size":"18"}}>Coal</div>
    <div style={{background:"red",color:"white",padding:2,"font-size":"18"}}>Nuclear</div>
    <div style={{background:"orange",color:"white",padding:2,"font-size":"18"}}>Natural Gas</div>
    <div style={{background:"blue",color:"white",padding:2,"font-size":"18"}}>Hydro</div>
    </div>
  )}
  */}

  {genlayeractive && (
    <div style={{width:300, height:300, position:"absolute",right:0,bottom:0,background:"white",padding:2}}>
      <Doughnut data={chartdata} />
    </div>
  )}

    </DeckGL>

    </>


    
  );
}

const rootElement = document.getElementById("root");

createRoot(rootElement).render(<App />)
