import sys

filename = sys.argv[1]

with open('module_casedata.js', 'w') as f:
    f.write('// ExaGo Viz Input File\n')
    f.write('\n')
    f.write('module.exports = {\n')
    f.write('\n')
    f.write('\tget_casedata: function () {\n')
    f.write('\t\t\t\tvar inputcasedata = require("./data/' + filename + '");\n')
    f.write('\n')
    f.write('\t\t\t\tvar casedata0 = {};\n')
    f.write('\t\t\t\tcasedata0.geojsondata = {};\n')
    f.write('\t\t\t\tcasedata0.nareas = inputcasedata.nareas;\n')
    f.write('\t\t\t\tcasedata0.nzones = inputcasedata.nzones;\n')
    f.write('\t\t\t\tcasedata0.areas = inputcasedata.areas;\n')
    f.write('\t\t\t\tcasedata0.zones = inputcasedata.zones;\n')
    f.write('\t\t\t\tcasedata0.geojsondata.type = "FeatureCollection";\n')
    f.write(
        '\t\t\t\tcasedata0.geojsondata.features = [...inputcasedata.geojsondata.features];\n')
    f.write('\t\t\t\treturn casedata0;\n')
    f.write('\t\t\t}\n')
    f.write('\n')
    f.write('};\n')
