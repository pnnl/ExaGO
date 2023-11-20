from mpi4py import MPI
import exago
import os

# Parallel communicator
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

# Initialize ExaGO
exago.initialize("app", comm)

# Create security-constrained ACOPF (SCOPF) object
scopf = exago.SCOPFLOW()

netfile = '/Users/abhy245/software/ExaGO/datafiles/case9/case9mod.m'
ctgcfile = '/Users/abhy245/software/ExaGO/datafiles/case9/case9.cont'

# Set network data
scopf.set_network_data(netfile)

# Set contingency data file and input format
scopf.set_contingency_data(ctgcfile, exago.ContingencyFileInputFormat.NATIVE)

# Set number of contingencies, we select all (-1)
scopf.set_num_contingencies(-1)

# Set solver
scopf.set_solver("EMPAR")

# Solve SCOPF
scopf.solve()

# Print base-case solution
scopf.print_solution(0)

# Save all contingency solutions to directory case9scopf, use MATPOWER formatted files
scopf.save_solution_all(exago.OutputFormat.MATPOWER, 'case9scopf')

# delete scopf object
del scopf

# Close ExaGO
exago.finalize()
