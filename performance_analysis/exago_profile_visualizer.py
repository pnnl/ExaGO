import json
import os
import pickle
import sys
import toml

#
# Print Debug levels
# 0: basic, will only print default output like execution time.
# 1: moderate, Additionally prints ExaGO output.
# 2: all, Additionally prints error and success messages.
DEBUG = 2


# Defining the debug level.
# Just a wrapper to print function, with the conditional debug level.
def printDebug(level, log_str):
    if DEBUG >= level:
        if level != 1:
            print('Auto Profiler Visualizer Log ======> ', end='')
        print(log_str)
        
def parsePickleFile(in_file):
    testsuite = toml.load(in_file)
    if 'application' not in testsuite:
        printDebug(0, "Please provide application")
    if 'presets' not in testsuite:
        printDebug(0, "Please provide presets")

    app_name = testsuite['application']
    preset_params = testsuite['presets']
    testsuite_name = "sample_testsuite"
    if 'testsuite_name' in testsuite:
        testsuite_name = testsuite['testsuite_name']
    
    t_file = testsuite_name + '.pkl'
    if os.path.exists(t_file):
        results_data = []
        with open(t_file, 'rb') as fp:
            try:
                while True:
                    results_data.append(pickle.load(fp))
            except EOFError:
                pass
        with open(testsuite_name + ".json", "w") as outfile:
            json.dump(results_data, outfile)
    else:
        printDebug(0, in_file + ' not found')
    

if __name__ == '__main__':
    in_file = "sample_testsuite.toml"
    if len(sys.argv) > 1:
        in_file = sys.argv[1]
    else:
        printDebug(0, 'No toml file provided. Using default file: ' + in_file)
    if os.path.exists(in_file):
        parsePickleFile(in_file)
    else:
        printDebug(0, in_file + ' not found')