import json
import os
import pickle
import sys
import toml
import matplotlib.pyplot as plt
import pandas as pd

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
        
def convertPklToJson(testsuite_name):
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
        printDebug(0, testsuite_name + ' not found')

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
    convertPklToJson(testsuite_name)

def visualizeProfiledData(in_file):
    profiledData = json.load(open(in_file, 'r'))
    ## BAR GRAPH ##
    xAxis = []
    cpu_values = []
    gpu_values = []
    xGroups = ['CPU', 'GPU']
    for e_data in profiledData:
        if e_data['opflow_model'] == 'POWER_BALANCE_HIOP' and \
        e_data['netfile'] == 'datafiles/case_ACTIVSg200.m':
            if e_data['compute_mode'] == 'CPU':
                cpu_values.append(e_data['petsc_solve_time'])
                xAxis.append(e_data['max_iter'])
            else:
                gpu_values.append(e_data['petsc_solve_time'])

    df = pd.DataFrame({'CPU': cpu_values,
                    'GPU': gpu_values}, index=xAxis)
    ax = df.plot.bar(rot=0, color={"CPU": "green", "GPU": "red"})

    plt.show()




if __name__ == '__main__':
    # convertPklToJson('OPFLOW_frontier_profiling')
    in_file = "profile_dumps/OPFLOW_frontier_profiling_2000.json"
    if len(sys.argv) > 1:
        in_file = sys.argv[1]
    else:
        printDebug(0, 'No JSON file provided. Using default file: ' + in_file)
    if os.path.exists(in_file):
        visualizeProfiledData(in_file)
    else:
        printDebug(0, in_file + ' not found')