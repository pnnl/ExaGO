from typing import Dict

_model_names = ('POWER_BALANCE_POLAR', 'POWER_BALANCE_CARTESIAN', 'POWER_BALANCE_HIOP', 'PBPOLRAJAHIOP')
_solver_names = ('HIOP', 'TAO', 'IPOPT', 'HIOPSPARSE')
_initialization_types: Dict[str, int] = {
  'MIDPOINT' : 0, 
  'FROMFILE' : 1, 
  'ACPF' : 2, 
  'FLATSTART' : 3
}
_genbusvoltage_types: Dict[str, int] = {
  'VARIABLE_WITHIN_BOUNDS' : 0, 
  'FIXED_WITHIN_QBOUNDS' : 1, 
  'FIXED_AT_SETPOINT' : 2
}
_objective_function_types: Dict[str, int] = {
  'MIN_GEN_COST' : 0,
  'MIN_GENSETPOINT_DEVIATION' : 1,
  'NO_OBJ' : 2
}
_hiop_compute_modes = ('auto', 'cpu', 'hybrid', 'gpu')
_output_formats: Dict[str, int] = {
  'CSV' : 0,
  'MATPOWER' : 1,
}

def modelnames() -> (str):
    '''Return all opflow model names for ExaGO'''
    return _model_names

def solvernames() -> (str):
    '''Return all opflow solver names for ExaGO'''
    return _solver_names

def initialization_types() -> Dict[str, int]:
    '''Return all opflow intialization types for ExaGO'''
    return _initialization_types

def genbusvoltage_types() -> Dict[str, int]:
    '''Return all opflow genbusvoltage types for ExaGO'''
    return _genbusvoltage_types

def objective_function_types() -> Dict[str, int]:
    '''Return all opflow objective function types for ExaGO'''
    return _objective_function_types 

def hiop_compute_modes() -> (str):
    '''Return all hiop compute modes for ExaGO'''
    return _hiop_compute_modes

def output_formats() -> Dict[str, int]:
    '''Return all output formats for ExaGO output'''
    return _output_formats

