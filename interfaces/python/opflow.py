import os
import sys
from . import config
from ctypes import *
import logging


class OPFLOW:
    """Python bindings for ExaGO's OPFLOW library"""
    
    # This makes sure finalize is only called once the last object is deleted
    list = [] 

    def __init__(self):
        self.finalize = True
        OPFLOW.list.append(0) 
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
        OPFLOW.list.pop()
        logging.debug('Deleting OPFLOW object')
        self.opflowlib.OPFLOWDestroy(byref(self.opflow))
        if (len(OPFLOW.list) == 0 and self.finalize):
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

    def solve(self):
        logging.debug('Calling OPFLOWSolve')
        self.opflowlib.OPFLOWSolve(self.opflow)

