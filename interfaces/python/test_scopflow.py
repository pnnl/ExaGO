import exago
import os

# Initialize a given application to be run
exago.initialize("scopflow_test")
scopf = exago.SCOPFLOW()
path = exago.prefix()
scopf.set_network_data(
    os.path.join(path, 'share', 'exago', 'datafiles', 'case9', 'case9mod.m'))
scopf.set_contingency_data(
    os.path.join(path, 'share', 'exago', 'datafiles', 'case9', 'case9.cont'),
    exago.ContingencyFileInputFormat.NATIVE)
scopf.set_num_contingencies(4)
scopf.set_solver("IPOPT")
scopf.solve()
scopf.print_solution(0)

# Delete instance before finalization (try, at least)
del scopf

exago.finalize()
