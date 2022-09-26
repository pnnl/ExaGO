// ExaGo Viz Input File

module.exports = {

	get_casedata: function () {

				var casedata500 = require('./data/case_ACTIVSg500.json');
				var casedata200 = require('./data/case_ACTIVSg200.json');
				var casedata2k = require('./data/case_ACTIVSg2000.json');
				var casedata10k = require('./data/case_ACTIVSg10k.json');
				var casedata70k = require('./data/case_ACTIVSg70k.json');

				var casedata0 = {};
				casedata0.geojsondata = {};
				casedata0.geojsondata.type = "FeatureCollection";
				casedata0.geojsondata.features = [...casedata10k.geojsondata.features,...casedata2k.geojsondata.features,...casedata70k.geojsondata.features];
				return casedata0;
			}

};