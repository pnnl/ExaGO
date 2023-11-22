// ExaGo Viz Input File

module.exports = {

	get_casedata: function () {
				var inputcasedata = require("./data/case_ACTIVSg10k.json");

				var casedata0 = {};
				casedata0.geojsondata = {};
				casedata0.geojsondata.type = "FeatureCollection";
				casedata0.geojsondata.features = [...inputcasedata.geojsondata.features];
				return casedata0;
			}

};
