from exago.opflow import OPFLOW 
import os
from exago import config

opf = OPFLOW()
opf.read_mat_power_data(
        os.path.join(config.prefix(), 'share', 'exago', 'datafiles', 'case9', 'case9mod.m'))
opf.solve()
opf.solution_to_ps()
print(f'OBJECTIVE VAL = {opf.objective_function}')
