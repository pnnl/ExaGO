import os
import sys
from . import config
from exago.typedefs import opflow 
from ctypes import *
import logging


class OPFLOW:
    """Python bindings for ExaGO's OPFLOW library"""
    
    # This makes sure finalize is only called once the last object is deleted
    ref_counter = 0

    def __init__(self):
        self.finalize = True
        OPFLOW.ref_counter += 1 
        logging.basicConfig(
            format=f'[ExaGO {os.path.split(__file__)[-1]} %(levelname)s]: %(message)s',
            level=logging.DEBUG)
        logging.debug('Loading dll %s' % config.libraries()['opflow'])
        self.opflowlib = cdll.LoadLibrary(config.libraries()['opflow'])

        # Get self MPI communicator
        self.comm = c_long()
        argc = c_int(0)
        logging.debug('Calling ExaGOGetSelfCommunicator')
        self.opflowlib.ExaGOGetSelfCommunicator(byref(self.comm))

        logging.debug('Calling ExaGOInitialize')
        self.opflowlib.ExaGOInitialize(self.comm, byref(
            argc), None, "opflow".encode('ascii'), None)

        # Create OPFLOW application
        self.opflow = c_longlong(0)
        logging.debug('Calling OPFLOWCreate')
        self.opflowlib.OPFLOWCreate(self.comm, byref(self.opflow))

    # Use this to not call MPI finalize when __del__ is called
    #   - This overrides in-built stack of opflow objects
    def dont_finalize(self):
        self.finalize = False
    
    def __del__(self):
        OPFLOW.ref_counter -= 1
        logging.debug('Deleting OPFLOW object')
        self.opflowlib.OPFLOWDestroy(byref(self.opflow))
        if (OPFLOW.ref_counter == 0 and self.finalize):
          logging.debug('Calling ExaGOFinalize')
          self.opflowlib.ExaGOFinalize()

    def read_mat_power_data(self, filename):
        logging.debug('Reading MatPower file %s', filename)
        self.opflowlib.OPFLOWReadMatPowerData(self.opflow, filename.encode('ascii'))

    def setup_ps(self):
        logging.debug('Calling OPFLOWSetUpPS')
        self.opflowlib.OPFLOWSetUpPS(self.opflow)

    def get_ps(self):
        ps = c_longlong(0)
        logging.debug('Calling OPFLOWGetPS')
        self.opflowlib.OPFLOWGetPS(self.opflow, byref(ps))
        return ps
    
    def ps_set_gen_power_limits(self, gbus, gid, pt, pb, qt, qb):
        ps = self.get_ps() 
        logging.debug('Setting PS object power limits')
        self.opflowlib.PSSetGenPowerLimits(ps,gbus,gid.encode('ascii') ,c_double(pt), c_double(pb),c_double(qt),c_double(qb))

    def solution_to_ps(self):
        logging.debug('Storing opflow solution to PS object')
        self.opflowlib.OPFLOWSolutionToPS(self.opflow)

    @property
    def objective_function(self):
        logging.debug('Getting objective function solution')
        obj = c_double(0.0)
        self.opflowlib.OPFLOWGetObjective(self.opflow,byref(obj))
        return obj.value

    def get_gen_dispatch(self, gbus, gid):
        logging.debug('Getting generator dispatch from opflow')
        ps = self.get_ps()
        pg = c_double()
        qg = c_double()
        self.opflowlib.PSGetGenDispatch(ps,gbus,gid.encode('ascii'),byref(pg), byref(qg));
        return pg.value, qg.value

    def set_model(self, modelname):
        logging.debug(f'Setting opflow model to {modelname}')
        if modelname not in opflow.modelnames():
          raise TypeError(f'No OPFLOW model named {modelname}. Try any of {modelnames}.')
        self.opflowlib.OPFLOWSetModel(self.opflow, modelname.encode('ascii'))

    def set_solver(self, solvername):
        logging.debug(f'Setting opflow solver to {solvername}')
        solvernames = opflow.solvernames()
        if solvername not in solvernames:
          raise TypeError(f'No OPFLOW solver named {solvername}. Try any of {solvernames}')
        self.opflowlib.OPFLOWSetSolver(self.opflow, solvername.encode('ascii'))

    def set_initialization(self, initialization_type):
        logging.debug(f'Setting opflow initialization to type {initialization_type}')
        initialization_types = opflow.initialization_types() 
        if initialization_type not in initialization_types.keys():
          raise TypeError(f'No OPFLOW initialization type matches {initialization_type}. Try any of {list(initialization_types.keys())}')
        self.opflowlib.OPFLOWSetInitializationType(self.opflow, initialization_types[initialization_type])

    def set_ignore_lineflow_constraints(self, ignore_constraints):
        logging.debug(f'Setting opflow ignore lineflow constraints to {ignore_constraints}')
        if not isinstance(ignore_constraints, (bool)):
          raise TypeError(f'Use a boolean to set opflow ignore lineflow constraints. {ignore_constraints} is invalid.')
        self.opflowlib.OPFLOWIgnoreLineflowConstraints(self.opflow, ignore_constraints)

    def set_include_loadloss(self, loadloss):
        logging.debug(f'Setting opflow include loadloss to {loadloss}')
        if not isinstance(loadloss, (bool)):
          raise TypeError(f'Use a boolean to set opflow include loadloss. {loadloss} is invalid.')
        self.opflowlib.OPFLOWHasLoadLoss(self.opflow, loadloss)

    def set_include_powerimbalance(self, powerimbalance):
        logging.debug(f'Setting opflow include powerimbalance to {powerimbalance}')
        if not isinstance(powerimbalance, (bool)):
          raise TypeError(f'Use a boolean to set opflow include powerimbalance. {powerimbalance} is invalid.')
        self.opflowlib.OPFLOWHasBusPowerImbalance(self.opflow, powerimbalance)

    def set_loadloss_penalty(self, penalty):
        logging.debug(f'Setting opflow loadloss penalty to {penalty}')
        self.opflowlib.OPFLOWSetLoadLossPenalty(self.opflow, c_double(penalty))

    def set_powerimbalance_penalty(self, penalty):
        logging.debug(f'Setting opflow power imbalance penalty to {penalty}')
        self.opflowlib.OPFLOWSetBusPowerImbalancePenalty(self.opflow, c_double(penalty))

    def set_genbusvoltage(self, control_mode):
        logging.debug(f'Setting opflow genbusvoltage to type {control_mode}')
        genbusvoltage_types = opflow.genbusvoltage_types()
        if control_mode not in genbusvoltage_types.keys():
          raise TypeError(f'No OPFLOW genbusvoltage type matches {control_mode}. Try any of {list(genbusvoltage_types.keys())}')
        self.opflowlib.OPFLOWSetGenBusVoltageType(self.opflow, genbusvoltage_types[control_mode])

    def set_has_gensetpoint(self, gensetpoint):
        logging.debug(f'Setting opflow include gensetpoint to {gensetpoint}')
        if not isinstance(gensetpoint, (bool)):
          raise TypeError(f'Use a boolean to set opflow has gensetpoint. {gensetpoint} is invalid.')
        self.opflowlib.OPFLOWHasGenSetPoint(self.opflow, gensetpoint)

    def set_objective_function(self, function):
        logging.debug(f'Setting opflow objective function to type {function}')
        objective_function_types = opflow.objective_function_types() 
        if function not in objective_function_types.keys():
          raise TypeError(f'No OPFLOW objective function type matches {function}. Try any of {list(objective_function_types.keys())}')
        self.opflowlib.OPFLOWSetObjectiveType(self.opflow, objective_function_types[function])

    def set_use_agc(self, agc):
        logging.debug(f'Setting opflow AGC for generator real power dispatch to {agc}')
        if not isinstance(agc, (bool)):
          raise TypeError(f'Use a boolean to set opflow AGC for generator real power dispatch. {agc} is invalid.')
        self.opflowlib.OPFLOWUseAGC(self.opflow, agc)

    def set_tolerance(self, tol):
        logging.debug(f'Setting tolerance for optimization solver to {tol}')
        self.opflowlib.OPFLOWSetTolerance(self.opflow, c_double(tol))

    def set_hiop_compute_mode(self, mode):
        logging.debug(f'Setting hiop compute mode to {mode}')
        modes = opflow.hiop_compute_modes()
        if mode not in modes:
          raise TypeError(f'No hiop compute mode named {mode}. Try and of {modes}')
        self.opflowlib.OPFLOWSetHIOPComputeMode(self.opflow, mode.encode('ascii'))

    def set_hiop_verbosity(self, level):
        logging.debug(f'Setting hiop verbosity level to {level}')
        if not isinstance(level, (int)) or not (0 <= level <= 10):
          raise TypeError(f'Unable to set hiop verbosity level to {level}. Please specify an integer between 0 and 12 inclusive.')
        self.opflowlib.OPFLOWSetHIOPVerbosityLevel(self.opflow, c_int(level))

    def print_solution(self):
        logging.debug(f'Printing opflow output')
        self.opflowlib.OPFLOWPrintSolution(self.opflow)

    def save_solution(self, fmt, outfile):
        logging.debug(f'Saving solution to {outfile} in {fmt} format')
        formats = opflow.output_formats() 
        if fmt not in formats.keys():
          raise TypeError(f'Unable to use format {fmt} for output. Use one of {list(formats.keys())}.')
        self.opflowlib.OPFLOWSaveSolution(self.opflow, formats[fmt], outfile.encode('ascii'))

    def solve(self):
        logging.debug('Calling OPFLOWSolve')
        self.opflowlib.OPFLOWSolve(self.opflow)

