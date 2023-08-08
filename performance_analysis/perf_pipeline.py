#!/usr/bin/env python3
import json
import pickle
import toml
import sys
import shutil
import os
import io
import time
import re
import math
import copy
from subprocess import call, Popen, PIPE

#
# Print Debug levels
# 0: basic, will only print default output like execution time.
# 1: moderate, Additionally prints ExaGO output.
# 2: all, Additionally prints error and success messages.
DEBUG = 1


# Defining the debug level.
# Just a wrapper to print function, with the conditional debug level.
def printDebug(level, log_str):
    if DEBUG >= level:
        if level != 1:
            print('Auto Profiler Log ======> ', end='')
        print(log_str)


# populate key and values for hiop.options or ipopt.opt files
# Create necessary input file for Exago solvers
def updateHiopOptions(options, app):
    hiops = dict()
    delList = list()
    app_prefix = ''
    op_file_name = ''
    st_ind = 1
    if options[app + '_solver'] == 'HIOP' or \
            options[app + '_solver'] == 'HIOPSPARSE':
        app_prefix = 'hiop_'
        op_file_name = 'hiop.options'
        st_ind = 5
    elif options[app + '_solver'] == 'IPOPT':
        app_prefix = 'ipopt_'
        op_file_name = 'ipopt.opt'
        st_ind = 6

    if st_ind == 1:
        hiops['max_iter'] = 1
    else:
        for key in options:
            if key.startswith(app_prefix):
                hiop_key = key[st_ind:]
                hiop_value = options[key]
                hiops[hiop_key] = hiop_value
                delList.append(key)
        with open(op_file_name, 'w') as hof:
            for key in hiops:
                hof_str = key + ' ' + str(hiops[key])
                hof.write(hof_str + '\n')
        for dl in delList:
            # hiop_mem_space is also argument to exago
            if not dl == 'hiop_mem_space':
                options.pop(dl, None)
    return hiops


# Using regex to find the OPFLOWSolve line from the PETSc log
def getSolveTimeFromPetsc(line, appName):
    pr = '^(OPFLOWSolve)\\s+(\\d+)\\s+(\\d+\\.\\d+)\\s+(\\d+\\.\\d+e[+-]\\d+)'
    if appName.casefold() == 'opflow' or appName.casefold() == 'sopflow':
        solver_line_parser = re.compile(pr)
        solver_line_match = solver_line_parser.match(line)
        if solver_line_match is not None:
            printDebug(2, 'solver line found')
            return float(solver_line_match.group(4))
    return -1


# Printing logs, iteration count, total and average time
def printTimingInfo(params,
                    hiop_options,
                    time_lists,
                    petsc_time,
                    app_name,
                    tool_name):
    avg_tm = sum(time_lists) / len(time_lists)
    var = sum(pow(x-avg_tm, 2) for x in time_lists) / len(time_lists)
    avg_tm_str = str(round(avg_tm, 5))
    avg_std = str(round(math.sqrt(var), 5))
    cm = 'CPU'

    if 'compute_mode' in hiop_options:
        cm = hiop_options['compute_mode']
    printDebug(0, "With " + tool_name + " Total Iterations: " +
               str(len(time_lists)) + ", " + cm + " Average time per testcase: " +
               avg_tm_str + " seconds, std: " + avg_std)
    if len(petsc_time) > 0:
        avg_petsc = sum(petsc_time) / len(petsc_time)
        avg_petsc_str = str(round(avg_petsc, 5))
        printDebug(0, "PETSc reported Solve Time per testcase: " +
                   avg_petsc_str + " seconds")

    if 'max_iter' in hiop_options:
        hiop_avg_time = str(round(avg_tm / int(hiop_options['max_iter']), 5))
        solver_name = params[app_name + '_solver']
        printDebug(0, "Total " + solver_name + " iterations: " +
                   str(hiop_options['max_iter']) + ", Average time per " +
                   solver_name + " iteration: " + hiop_avg_time + " seconds")
        if len(petsc_time) > 0:
            p_avg_time = str(round(avg_petsc / int(hiop_options['max_iter']), 5))
            printDebug(0, "PETSc reported Solve time per iteration: "
                       + p_avg_time)

# writing results to a python pickle file
# Used pickle format, since I needed to append to the file
# Also for a larger run, this script could terminate in the middle and generated pickle file would still be here.
def dumpResultsToFile(testsuite_name, params, hiop_options, time_lists, petsc_time, app_name):
    write_filename = testsuite_name + ".pkl"
    dict_dump = dict()
    dict_dump = copy.deepcopy(params)
    dict_dump['application'] = app_name
    dict_dump.update(hiop_options)
    dict_dump['petsc_solve_time'] = sum(petsc_time) / len(petsc_time)
    dict_dump['total_solve_time'] = sum(time_lists) / len(time_lists)
    dict_dump['profiling_time'] = time.strftime('%m/%d/%Y %H:%M:%S')
    with open(write_filename,'ab+') as fp:
        pickle.dump(dict_dump, fp)

def flattenNestedTestCases(params):
    c_tests = dict()
    keysList = list(params.keys())
    row_length = len(keysList)
    commands_list = list()
    def traversingNestingDict(row_check):
        if row_check == row_length:
            commands_list.append(copy.deepcopy(c_tests))
            return
        if isinstance(params[keysList[row_check]], str):
            t_envs = params[keysList[row_check]].split()
        else:
            t_envs = [params[keysList[row_check]]]
        for env in t_envs:
            c_tests[keysList[row_check]] = env
            traversingNestingDict(row_check+1)

    traversingNestingDict(0)
    return commands_list
 
def checkAppStartCondition(app_name, line):
    sopflow_start_string = 'SOPFLOW: Application created'
    opflow_start_string = '[ExaGO] Creating OPFlow'
    if app_name.casefold() == 'sopflow':
        if sopflow_start_string in line:
            return True
        else:
            return False
    elif app_name.casefold() == 'opflow':
        if opflow_start_string in line:
            return True
        else:
            return False
    is_app_started = True

def executeCommandAndMeasureTime(command, my_env, iterations, app_name, tool_name):
    suc_str = 'Finalizing ' + app_name + ' application.'
    suc_str_term = '-log_view'
    alt_suc_str_term = 'Maximum number of iterations reached.'
    ON_POSIX = 'posix' in sys.builtin_module_names

    time_lists = list()
    petsc_time = list()
    for i in range(iterations):
        timeDelta = 0

        input_fd, output_fd = os.pipe()
        proc = Popen(command, stdout=output_fd,
                        universal_newlines=True, bufsize=1,
                        env=my_env, close_fds=ON_POSIX)
        os.close(output_fd)
        totalPetSCTime = 0.0
        exago_success_run = False
        is_app_started = False
        with io.open(input_fd, 'r', buffering=1) as ff:
            for line in ff:
                if is_app_started is False:
                    is_app_started = checkAppStartCondition(app_name, line)
                    if is_app_started is True: # its important to keep this nested
                        timeStarted = time.time()
                else:
                    timeDelta = timeDelta + time.time() - timeStarted
                    printDebug(1, line.rstrip())
                    petSCTime = getSolveTimeFromPetsc(line, app_name)
                    if float(petSCTime) > -1:
                        totalPetSCTime = totalPetSCTime + float(petSCTime)
                    if exago_success_run is False and suc_str in line:
                        exago_success_run = True
                    if exago_success_run is False and suc_str_term in line:
                        exago_success_run = True
                    if exago_success_run is False and alt_suc_str_term in line:
                        exago_success_run = True
                    timeStarted = time.time()

        if exago_success_run:
            petsc_time.append(totalPetSCTime)
            printDebug(2, app_name + " runs successfully")
            time_lists.append(timeDelta)
            printDebug(2, "Total measured time with " +
                        tool_name + ": " + str(round(timeDelta, 5)) +
                        " seconds.")
        else:
            printDebug(0, app_name + " did NOT run with " + tool_name)
    return time_lists, petsc_time

# doing automated performance measurement
# this will read from TOML file, parse it, and populate
# the data into a dictionary. Then it will execute the application
# as a subprocess with the provided command line arguments from TOML file.
def doPerfMeasure(in_file):
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
    iterations = 1
    if 'iterations' in testsuite:
        iterations = testsuite['iterations']

    mpi_cmd = None
    mpi_start = "mpiexec"
    if "mpi_start" in testsuite:
        mpi_start = str(testsuite['mpi_start'])

    if 'mpi_args' in testsuite:
        mpi_cmd = [mpi_start]
        mpi_cmd.extend(testsuite['mpi_args'].split())
#    if "mpi_rank" in testsuite:
#        mpi_cmd = [mpi_start, "-n", str(testsuite['mpi_rank'])]
    
    profiler_cmd_list = list()
    my_env = os.environ.copy()
    if "profiler" in testsuite:
        for prof in testsuite['profiler']:
            if 'tool' not in prof:
                prof_tool = 'no tool'
                profiler_cmd_list.append([])
                continue
            prof_tool = prof['tool']
            if shutil.which(prof_tool) is None:
                printDebug(0, prof_tool + ' not installed')
            else:
                printDebug(2, prof_tool + ' found')
                pcmd = [prof_tool]
                if 'tool_args' in prof:
                    pcmd.extend(prof['tool_args'].split())
                if profiler_cmd_list is None:
                    profiler_cmd_list = list()
                profiler_cmd_list.append(pcmd)

                if "tool_envs" in prof:
                    t_envs = prof["tool_envs"].split()
                    for env in t_envs:
                        kv = env.split('=')
                        printDebug(2, 'env key: ' + kv[0] +
                                   ' env value: ' + kv[1])
                        my_env[kv[0]] = kv[1]
                    printDebug(2, 'profile dir set')
        printDebug(2, profiler_cmd_list)

    if shutil.which(app_name) is None:
        printDebug(0, app_name + ' not installed')
        return

    test_number = 1
    for tests in testsuite['testcase']:
        for profiler_cmd in profiler_cmd_list:
            tests_list = flattenNestedTestCases(tests)
            for each_test in tests_list:
                command = list()
                if mpi_cmd is not None:
                    command.extend(mpi_cmd)

                tool_name = 'no tool'
                if profiler_cmd:
                    tool_name = profiler_cmd[0]
                    command.extend(profiler_cmd)

                command.append(app_name)
                params = dict()
                params = copy.deepcopy(preset_params)
                params.update(each_test)
                hiop_options = updateHiopOptions(params, app_name)

                for key in params:
                    if key == 'argument_list':
                        cls = params[key].split()
                        for cl in cls:
                            command.append('-' + cl)
                    else:
                        command.append('-' + key)
                        command.append(str(params[key]))

                printDebug(0, "For testcase " + str(test_number))
                test_number = test_number + 1
                printDebug(2, command)
                printDebug(2, '----')
                time_lists, petsc_time = executeCommandAndMeasureTime(command, my_env, iterations, app_name, tool_name)
                if len(time_lists) > 0:
                    printTimingInfo(params, hiop_options, time_lists, petsc_time, app_name, tool_name)
                    dumpResultsToFile(testsuite_name, params, hiop_options, time_lists, petsc_time, app_name)

def parsePickleFileToJSON(in_file):
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
        os.system('rm ' + t_file)
    else:
        printDebug(0, t_file + ' not found')

# Main file. This checks if a file name is provided
# if not it will try to read the default file
# and calls the doPerfMeasure
if __name__ == '__main__':
    in_file = "sample_testsuite.toml"
    if len(sys.argv) > 1:
        in_file = sys.argv[1]
    else:
        printDebug(0, 'No toml file provided. Using default file: ' + in_file)
    if os.path.exists(in_file):
        doPerfMeasure(in_file)
        # converting it to JSON for easy readability
        parsePickleFileToJSON(in_file)
    else:
        printDebug(0, in_file + ' not found')
