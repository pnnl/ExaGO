import exago
import os

# Initialize a given application to be run
exago.initialize("pflow")
pf = exago.PFLOW()
path = exago.prefix()
pf.read_mat_power_data(
    os.path.join(path, 'share', 'exago', 'datafiles', 'case9', 'case9mod.m'))
pf.solve()

# Delete instance before finalization
del pf
exago.finalize()
