import exago
import os

# Initialize a given application to be run
exago.initialize("opflow")
opf = exago.opf()
path = exago.prefix()
opf.read_mat_power_data(
    os.path.join(path, 'share', 'exago', 'datafiles', 'case9', 'case9mod.m'))
opf.solve()
opf.print_solution()
print(f'OBJECTIE FUNCTION VALUE : {opf.get_objective()}')
exago.finalize()
