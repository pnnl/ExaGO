from mpi4py import MPI
import exago
import os

# Parallel communicator
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

# Input files
netfile = '/Users/abhy245/software/ExaGO/datafiles/case9/case9mod.m'
ctgcfile = '/Users/abhy245/software/ExaGO/datafiles/case9/case9.cont'
scenfile = '/Users/abhy245/software/ExaGO/datafiles/case9/10_scenarios_9bus.csv'

# Initialize ExaGO
exago.initialize("app", comm)

# Create Stochastic optimal power flow (SOPF) object
sopf = exago.SOPFLOW()

# Set network data
sopf.set_network_data(netfile)

# Set contingency data file and input format
sopf.set_contingency_data(ctgcfile, exago.ContingencyFileInputFormat.NATIVE)

# Set scenario data
sopf.set_scenario_data(scenfile, exago.ScenarioFileInputFormat.NATIVE_SINGLEPERIOD,
                       exago.ScenarioUncertaintyType.WIND)

# Set number of scenarios, we select all (-1)
sopf.set_num_scenarios(-1)

# Set number of contingencies, we select all (-1)
sopf.set_num_contingencies(-1)

# Enable multi-contingency
sopf.enable_multi_contingency(True)

# Set solver
sopf.set_solver("EMPAR")

# Solve SCOPF
sopf.solve()

# Print base-case solution
sopf.print_solution(0)

# Save all contingency solutions to directory case9scopf, use MATPOWER formatted files
sopf.save_solution_all(exago.OutputFormat.MATPOWER, 'case9sopf')

# delete object
del sopf

# Close ExaGO
exago.finalize()
