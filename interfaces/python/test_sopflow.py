import exago
import os

# Initialize a given application to be run
exago.initialize("sopflow_test")

sopf = exago.SOPFLOW()
path = exago.prefix()

sopf.set_tolerance(1.0E-03)
sopf.set_network_data(os.path.join(
    path, 'share', 'exago', 'datafiles', 'case9', 'case9mod_gen3_wind.m'))

sopf.set_num_scenarios(3)
sopf.set_scenario_data(os.path.join(
    path, 'share', 'exago', 'datafiles', 'case9', '10_scenarios_9bus.csv'),
    exago.ScenarioFileInputFormat.NATIVE_SINGLEPERIOD,
    exago.ScenarioUncertaintyType.WIND
)

sopf.enable_multi_contingency(True)


sopf.set_contingency_data(os.path.join(
    path, 'share', 'exago', 'datafiles', 'case9', 'case9.cont'),
    exago.ContingencyFileInputFormat.NATIVE
)
sopf.set_num_contingencies(1)

sopf.set_solver('IPOPT')

sopf.setup()
sopf.solve()

# Delete instance before finalization (try, at least)
del sopf

exago.finalize()
