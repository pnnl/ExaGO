import countydata from "../data/counties.json";
import { center, convex, bbox } from '@turf/turf';

var us = require('us')
const codeDict = us.mapping('fips', 'abbr')

function getCountyNodes(data) {
  var ncounties = countydata.features.length;
  var countygeojson = { type: "FeatureCollection", features: [], maxPd: 0.0 };
  for (var j = 0; j < ncounties; j++) {
    //    if(countydata.features[j].properties.STATE !== "13" & countydata.features[j].properties.STATE !== "25")
    //        continue;
    var polygon = countydata.features[j];
    var box = bbox(polygon);

    countydata.features[j].properties.Pd = 0.0;
    countydata.features[j].properties.Vm_avg = 0.0;
    countydata.features[j].properties.KVlevels = [];
    countydata.features[j].properties.countyname = countydata.features[j].properties.NAME + ", " +codeDict[countydata.features[j].properties.STATE];
    var countyhassubst = false;
    var nbuscounty = 0
    for (var i = 0; i < data.features.length; i++) {
      if (data.features[i].geometry.type == 'Point') {
        var lng = data.features[i].geometry.coordinates[0];
        var lat = data.features[i].geometry.coordinates[1];

        if ((box[0] <= lng) & (box[2] >= lng) & (box[1] <= lat) & (box[3] >= lat)) {
          //update data with county name, similar to spatial join
          data.features[i].properties.countyname = countydata.features[j].properties.NAME + ", " +codeDict[countydata.features[j].properties.STATE]

          countyhassubst = true;
          var subst = data.features[i].properties;
          var nbus = subst.nbus;

          countydata.features[j].properties.KVlevels = [...countydata.features[j].properties.KVlevels, ...data.features[i].properties.KVlevels];
          for (var k = 0; k < nbus; k++) {
            var bus = subst.bus[k];
            countydata.features[j].properties.Pd += bus.PD;
            countydata.features[j].properties.Vm_avg += bus.VM;
            nbuscounty++;

          }
        }
      }else{
        // add source and target info to linestring 
        if(data.features[i].properties.PF > 0){
          const [src, trg] = data.features[i].properties.NAME.split(' -- ')
          data.features[i].properties.source = src
          data.features[i].properties.target = trg
        }else{
          const [trg, src] = data.features[i].properties.NAME.split(' -- ')
          data.features[i].properties.source = src
          data.features[i].properties.target = trg
        }
        
      }
    }
    if (countyhassubst) {
      countydata.features[j].properties.Vm_avg /= nbuscounty;
      if (countydata.features[j].properties.Pd > countygeojson.maxPd) countygeojson.maxPd = countydata.features[j].properties.Pd;
      countygeojson.features.push(countydata.features[j]);
    }
  }
  return { maxPd: countygeojson.maxPd, data: countygeojson, updatedata: data }
}


// Extract data for first time-slice
function ExtractFirstTimeSlice(data) {
  var features = data.features.filter(function (feature) {
    if (!("start" in feature.properties) || (feature.properties.start == 0)) return feature;
  })
  var gdata = { "type": 'FeatureCollection', "features": features };
  return gdata;
}

function ExtractFlowData(data) {


  // name is the unique id 
  const locations = []
  const flows = []
  // const locationDict = {}
  data.features.forEach(feature => {
    if (feature.geometry.type === "Point") {
      locations.push({
        id: feature.properties.NAME,
        name: feature.properties.NAME,
        lon: feature.geometry.coordinates[0],
        lat: feature.geometry.coordinates[1]
      })
    } else if (feature.geometry.type === "LineString") {
      if (feature.properties.PF > 0) {
        const [origin, dest] = feature.properties.NAME.split(' -- ')
        flows.push({
          origin: origin,
          dest: dest,
          count: feature.properties.KV 
        })
      } else {
        const [dest, origin] = feature.properties.NAME.split(' -- ')
        flows.push({
          origin: origin,
          dest: dest,
          count: feature.properties.KV 
        })
      }

    }

  })

  // remove duplicated flows 
  const uniq = new Set(flows.map(e => JSON.stringify(e)));
  const res = Array.from(uniq).map(e => JSON.parse(e));
  return ({
    locations: locations, flows: res
  })

}

//get net data for bar
function getBarNet(data) {
  const netBarValue = data.foreach((line) => {
    line.properties.KV
  })
  return netBarValue
}

// Get substation voltages
function getPoints(data) {
  var Pointsi;
  var Points = [];
  var elev;
  var i, k = 0;
  for (i = 0; i < data.features.length; i++) {
    if (data.features[i].geometry.type == 'Point') {
      Pointsi = { coordinates: data.features[i].geometry.coordinates, value: data.features[i].properties.bus[0].VM };
      Points.push(Pointsi);
    }
  }
  return Points;
}

// Get generation
function getGeneration(data) {
  var Geni;
  var Gens = [];
  var minPg = 1000.0;
  var maxPg = 0.0;
  var elev;
  var i, j, k;
  var Pgcoal = 0.0, Pghydro = 0.0, Pgnuclear = 0.0, Pgng = 0.0, Pgsolar = 0.0, Pgwind = 0.0, Pgother = 0.0;
  var Pgcoalcap = 0.0, Pghydrocap = 0.0, Pgnuclearcap = 0.0, Pgngcap = 0.0, Pgsolarcap = 0.0, Pgwindcap = 0.0, Pgothercap = 0.0;
  for (i = 0; i < data.features.length; i++) {
    if (data.features[i].geometry.type == 'Point') {
      var subst = data.features[i].properties;
      var nbus = subst.nbus;
      var Pg = 0.0;
      var Pcap = 0.0;
      var gen_fuel;
      var ngen = 0;
      var KV = [];
      var name = subst.NAME;
      var countyname = subst.countyname
      for (j = 0; j < nbus; j++) {
        var bus = subst.bus[j];
        KV.push(bus.BASE_KV);
        for (k = 0; k < bus.ngen; k++) {
          var gen = bus.gen[k];
          Pg += gen.GEN_STATUS * gen.PG;
          Pcap += gen.GEN_STATUS * gen.PMAX;
          gen_fuel = gen.GEN_FUEL.toLowerCase();

          if (gen_fuel == 'wind') {
            Pgwind += gen.GEN_STATUS * gen.PG;
            Pgwindcap += gen.GEN_STATUS * gen.PMAX;
          }
          else if (gen_fuel == 'solar') {
            Pgsolar += gen.GEN_STATUS * gen.PG;
            Pgsolarcap += gen.GEN_STATUS * gen.PMAX;
          }
          else if (gen_fuel == 'coal') {
            Pgcoal += gen.GEN_STATUS * gen.PG;
            Pgcoalcap += gen.GEN_STATUS * gen.PMAX;
          } else if (gen_fuel == 'nuclear') {
            Pgnuclear += gen.GEN_STATUS * gen.PG;
            Pgnuclearcap += gen.GEN_STATUS * gen.PMAX;
          } else if (gen_fuel == 'hydro') {
            Pghydro += gen.GEN_STATUS * gen.PG;
            Pghydrocap += gen.GEN_STATUS * gen.PMAX;
          } else if (gen_fuel == 'ng') {
            Pgng += gen.GEN_STATUS * gen.PG;
            Pgngcap += gen.GEN_STATUS * gen.PMAX;
          } else {
            Pgother += gen.GEN_STATUS * gen.PG;
            Pgothercap += gen.GEN_STATUS * gen.PMAX;
          }
          ngen++;
        }
      }
      if (ngen) {
        var color;
        if (gen_fuel == 'wind') color = 'green'; else if (gen_fuel == 'solar') color = 'yellow';
        else if (gen_fuel == 'coal') color = 'gray';
        else if (gen_fuel == 'nuclear') color = 'red';
        else if (gen_fuel == 'hydro') color = 'blue';
        else if (gen_fuel == 'ng') color = 'orange'
        else color = 'black';
        if (Pg <= minPg) minPg = Pg;
        if (Pg >= maxPg) maxPg = Pg;
        Geni = { coordinates: data.features[i].geometry.coordinates, Pg: Pg, Pcap: Pcap, KVlevels: KV, color: color, 
         fuel: gen_fuel, name: name, countyname: countyname};
        Gens.push(Geni);
      }
    }
  }

  return { minPg: minPg, maxPg: maxPg, Gens: Gens, Pgwind: Pgwind, Pgsolar: Pgsolar, Pgnuclear: Pgnuclear, Pghydro: Pghydro, Pgng: Pgng, Pgcoal: Pgcoal, Pgother: Pgother, Pgwindcap: Pgwindcap, Pgsolarcap: Pgsolarcap, Pgnuclearcap: Pgnuclearcap, Pghydrocap: Pghydrocap, Pgngcap: Pgngcap, Pgcoalcap: Pgcoalcap, Pgothercap: Pgothercap };
}

// Get Load
function getLoad(data) {
  var Loadi;
  var Loads = [];
  var minPd = 1000.0;
  var maxPd = 0.0;
  var elev;
  var i, j;
  for (i = 0; i < data.features.length; i++) {
    if (data.features[i].geometry.type == 'Point') {
      var subst = data.features[i].properties;
      var nbus = subst.nbus;
      var Pd = 0.0;
      var Qd = 0.0;
      for (j = 0; j < nbus; j++) {
        var bus = subst.bus[j];
        Pd += bus.PD;
        Qd += bus.QD;
      }
      if (Pd > 0) {
        if (Pd <= minPd) minPd = Pd;
        if (Pd >= maxPd) maxPd = Pd;
        Loadi = { coordinates: data.features[i].geometry.coordinates, Pd: Pd };
        Loads.push(Loadi);
      }
    }
  }
  return { minPd: minPd, maxPd: maxPd, Loads: Loads };
}

function getContours() {
  var Vmin = 0.9;
  var Vmax = 1.1;

  var contours = [
    { threshold: [0.9, 0.98], color: [255, 0, 0] },
    { threshold: [0.98, 1.02], color: [0, 255, 0] },
    { threshold: [1.02, 1.1], color: [0, 0, 255] }
  ];
  return contours;
}

export { getCountyNodes, ExtractFirstTimeSlice, ExtractFlowData, getBarNet, getPoints, getGeneration, getLoad, getContours };
