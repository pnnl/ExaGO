from exago.opflow import OPFLOW
import os

opf = OPFLOW()
opf.read_mat_power_data(
    os.path.join(exago.prefix(), 'share', 'exago', 'datafiles', 'case9',
                 'case9mod.m'))
opf.solve()
print(f'OBJECTIVE VAL = {opf.objective_function}')
