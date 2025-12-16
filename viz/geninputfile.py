import sys

filename = sys.argv[1]

with open("src/module_casedata.js", "w") as f:

    f.write("// ExaGo Viz Input File\n")
    f.write("\n")
    f.write('import inputcasedata from "../data/' + filename + '" with { type: "json" };\n')
    f.write("\n")
    f.write("export default {\n")
    f.write("  get_casedata() {\n")
    f.write("    var casedata0 = {};\n")
    f.write("\n")
    f.write("    casedata0.geojsondata = {};\n")
    f.write("    casedata0.nareas = inputcasedata.nareas;\n")
    f.write("    casedata0.nzones = inputcasedata.nzones;\n")
    f.write("    casedata0.areas = inputcasedata.areas;\n")
    f.write("    casedata0.zones = inputcasedata.zones;\n")
    f.write('    casedata0.geojsondata.type = "FeatureCollection";\n')
    f.write("    casedata0.geojsondata.features = [...inputcasedata.geojsondata.features];\n")
    f.write("    return casedata0;\n")
    f.write("  },\n")
    f.write("};\n")
