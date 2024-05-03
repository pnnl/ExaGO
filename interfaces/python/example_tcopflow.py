import os
import exago

# Initialize a given application to be run
exago.initialize("tcopflow_test")
tcopf = exago.TCOPFLOW()
path = exago.prefix()
tcopf.set_network_data(
    os.path.join(path, 'share', 'exago', 'datafiles', 'case9', 'case9mod.m'))
tcopf.set_solver("IPOPT")
tcopf.solve()
tcopf.print_solution(0)

# Delete instance before finalization (try, at least)
del tcopf

exago.finalize()