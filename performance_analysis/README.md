Automated Performance Analysis
----------------------------------------------------

Documentation on how to automate the performance analysis pipeline. As conducting profiling and tracing have [different meanings](https://vampir.eu/tutorial/manual/introduction), throughout this documentation, the term **profiler** is used to indicate performance analysis tool.

This considers the **ExaGO** application and relevant performance analysis tools (like 
nvprof, HPCTookit, etc.) are already pre-installed (possibly using [Spack](https://spack.readthedocs.io/en/latest/features.html)). A sample **Spack** [environment](https://spack-tutorial.readthedocs.io/en/latest/tutorial_environments.html) file with Exago, [HPCToolkit](http://hpctoolkit.org/), and [TAU](https://hpc.llnl.gov/software/development-environment-software/tau-tuning-and-analysis-utilities) is provided [here](./spack.yaml). Note that, in the environment file, `cuda_arch=70` is a platform specific value and it should be changed to the appropriate *cuda_arch* version of the target test platform.

After doing the necessary installation, go to the [performance_analysis](./performance_analysis) directory from root directory of ExaGO source and simply run the performance analysis script:

```shell
python3 perf_pipeline.py
```

This script takes one command line argument parameter - a TOML file. If no argument is provided, it will read from a default file `sample_testsuite.toml`. 
To run the script with any other TOML file:

```shell
python3 perf_pipeline.py opflow_testsuite.toml
```

Note that, the script could be executed from any directory, considering relative paths are provided in the TOML file accordingly. Furthermore, this script also provides a setup to run ExaGO in a set of configurations without any performance analysis tools. Therefore, it is not necessary to use a performance analysis tools to use this script.

### Configuring TOML file
By default, the following keys should be provided in the TOML file:

- `testsuite_name` : Value could be any customized string for naming this testsuite.
- `application` : A string with an ExaGO application name, ex. `opflow` or `scopflow`.
- `mpi_rank` : An integer with the number of mpi ranks
- `iterations`: Number of times each testcase should be run. This is **NOT** hiop or ipopt iterations.

#### Providing ExaGO arguments
Also, **ExaGO** arguments which are supposed to be same across the testcases should be provided in a hash table titled `[presets]`. For example, defining `-netfile`, `-hiop_verbosity_level`, and  `-print_output` for all testcases can be done as follows,

```yaml
[presets]
netfile = '../datafiles/case_ACTIVSg200.m'
hiop_verbosity_level = 3
print_output = 0
argument_list = 'log_view'
```

Then each testcase should be provided with an array of tables titled `[[testcase]]`. For example, defining two testcases with different solvers can be done 
as follows,
```yaml
[[testcase]]
opflow_solver = 'HIOP'
opflow_model = 'POWER_BALANCE_HIOP'

[[testcase]]
opflow_solver = 'IPOPT'
opflow_model = 'POWER_BALANCE_POLAR'
```

Here, the keys for `presets` and `testcases` are merged into a single dictionary. Therefore, values for duplicate keys in `presets` will get overwritten with what is defined in `testcase`. Single ExaGO arguments (without any value) like `log_view` should be provided with key `argument_list` and values as string with space separated arguments. Note that, to show execution time reported by the PETSc log (for example OPFLOWSolve), the `log_view` argument will also enable collecting that value.

To set HIOP and IPOPT parameters, add `hiop_` and `ipopt_` before each key. For example,
```yaml
[[testcase]]
ipopt_print_level = 0
ipopt_max_iter = 100

[[testcase]]
hiop_compute_mode = 'CPU'
hiop_max_iter = 100
```
For the first testcase, it will create `ipopt.opt` file in the current directory and write the follows,
```
print_level = 0
max_iter = 100
```
For the second testcase, it will write to `hiop.options` file without the `hiop_` prefix.

#### Providing Performance Analysis Tool
To, execute ExaGO with a performance analysis tool, use an array of tables titled `[[profiler]]`. Here keys should be as follows,

- `tool` : Executable for the performance analysis tool. For example: `nvprof`, `hpcrun`, `tau_exec`, etc.
- `tool_args`: Arguments for the profiler. For example, to output nvprof to a log file, use `'--log-file nvprofOutput'`.
- `tool_envs`: Environment variables that are required to export for the profiler. This should be a string with space separated words, where each word is a key=value for each environment variable. For example, `'TAU_PROFILE=1 PROFILEDIR=./profiler_dumps'` set the environment variables `TAU_PROIFLE` and 
  `PROFILEDIR`. 

Use blank array table to execute without any profiler. For example, 
```yaml
[[profiler]]

[[profiler]]
tool = 'nvprof'
tool_args = '--log-file nvprofOutput'

[[profiler]]
tool = 'hpcrun'
tool_args = '--disable-auditor -t -o ./profiler_dumps'

[[profiler]]
tool = 'tau_exec'
tool_args = '-io'
tool_envs = 'TAU_PROFILE=1 PROFILEDIR=./profiler_dumps'
```

This will run ExaGO with
- no tool
- Nvidia nvprof
- HPCToolkit
- TAU

### What is being executed
This will run **ExaGO** application as a subprocess with each testcase by the defined number of `iterations` with each profiler, as follows:
```shell
mpiexec -n <mpi_rank> <profiler.tool> <profiler.tool_args> <application> <presets> <testcase>
```

If a blank array table with `[[profiler]]` is used, it will execute the following as a subprocess,

```shell
mpiexec -n <mpi_rank> <application> <presets> <testcase>
```
It allows running ExaGO in a set of configurations without any profiler tool.

### Output Log
The output log is printed directly in the standard output with the prefix `Auto Profiler Log ======> `. Output log can be controlled with the `DEBUG` 
variable which supports 3 values:
- 0: basic, will only print default output like execution time.
- 1: moderate, Additionally prints ExaGO output.
- 2: all, Additionally reports the parsed TOML keys, command line arguments, success, and error logs.

In the basic mode it prints the following,

```shell
With <profiler tool name> Total Iterations: <iteration>, <CPU/GPU> Average time per testcase: <time> seconds, std: <time standard deviation>
PETSc reported Solve Time per testcase: <time> seconds
Total <solver> iteration: <iteration>, Average time per <solver> iteration: <time> seconds
PETSc reported Solve time per iteration: <time>
```

The total wall-clock time of the subprocess is measured. Then, it is divided by the `iteration` number and reported in the first line of the output along with the standard deviation. If `log_view` is enabled, the reported *OPFLOWSolve* time is fetched from the output log of PETSc. Then, it is divided by the `iteration` number and reported in the second line of the output. If `hiop_max_iter` or `ipopt_max_iter` is defined, the measured wall-clock time is divided by `max_iter` and reported in the third line. In the fourth line, the PETSc reported *OPFLOWSolve* time is divided by the `max_iter` and reported if `log_view` is enabled.

Note that, the total measured wall-clock time should always be higher than the PETSc reported *OPFLOWSolve* time. Measured wall-clock time from an execution of ExaGO consists of Main Stage, Reading Data, Set up, and Solve stages. Whereas, *OPFLOWSolve* time reports time spent on the Solve stage.

### Analyzing Performance Tools Output
If the profiler tool is configured to print measurements in the standard output (for example, not providing `--log-file` or `--quiet` argument to `nvprof`), those output logs from the tool will directly be printed in the standard output (even in the basic mode with `DEBUG=0`). If specific output file is defined for the profiler tool, that output file can further be investigated to conduct more in-detailed analysis. For example, with the following setting,

```yaml
[[profiler]]
tool = 'nvprof'
tool_args = '--log-file nvprofOutput'
```

In this case, **nvprof** will print the measurement in file `nvprofOutput` in the current working directory. The top functions spent most time in executing can be found by running the following command

```shell
head -10 nvprofOutput
```

Furthermore, that output file could be transfered to local machine. Then a terminal based analysis tool (like [ParaProf](https://www.cs.uoregon.edu/research/tau/docs/paraprof/)) or any visualization tool (For example, [Jumpshot](https://www.anl.gov/mcs/jumpshot-performance-visualization-tool), [Vampir](https://vampir.eu/), etc.) could be utilized for further investigation.
