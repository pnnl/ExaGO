import exago
import os

# Initialize a given application to be run
exago.initialize("opflow")
opf = exago.OPFLOW()
path = exago.prefix()
opf.read_mat_power_data(
    os.path.join(path, 'share', 'exago', 'datafiles', 'case9', 'case9mod.m'))
opf.solve()
# opf.save_solution(exago.OutputFormat.MATPOWER, 'soln_from_test.m')
opf.print_solution()
print(f'OBJECTIE FUNCTION VALUE : {opf.get_objective()}')

# Delete instance before finalization
del opf
exago.finalize()
