// ExaGo Viz Input File

import inputcasedata from "../data/opflowout-70K.json" with { type: "json" };

export default {
  get_casedata() {
    var casedata0 = {};
    
    casedata0.geojsondata = {};
    casedata0.nareas = inputcasedata.nareas;
    casedata0.nzones = inputcasedata.nzones;
    casedata0.areas = inputcasedata.areas;
    casedata0.zones = inputcasedata.zones;
    casedata0.geojsondata.type = "FeatureCollection";
    casedata0.geojsondata.features = [...inputcasedata.geojsondata.features];
    return casedata0;
  },
};
