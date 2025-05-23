// ExaGo Viz Input File

module.exports = {

	get_casedata: function () {
	    // var inputcasedata = require("./data/WECC.json");
	                        var inputcasedata = require("./data/case2000.json");
				var casedata0 = {};
				casedata0.geojsondata = {};
				casedata0.nareas = inputcasedata.nareas;
				casedata0.nzones = inputcasedata.nzones;
				casedata0.areas = inputcasedata.areas;
				casedata0.zones = inputcasedata.zones;
				casedata0.geojsondata.type = "FeatureCollection";
				casedata0.geojsondata.features = [...inputcasedata.geojsondata.features];
				return casedata0;
			}

};
