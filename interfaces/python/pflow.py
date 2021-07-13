import os
import sys
from . import config
from ctypes import *
import logging


class PFLOW:
    """Python bindings for ExaGO's PFLOW library"""

    list = []
          
    def __init__(self):
        self.finalize = True
        PFLOW.list.append(0)
        logging.basicConfig(
            format=f'[ExaGO {os.path.split(__file__)[-1]} %(levelname)s]: %(message)s',
            level=logging.DEBUG)

        logging.debug('Loading dll %s' % config.libraries()['pflow'])
        self.pflowlib = cdll.LoadLibrary(config.libraries()['pflow'])

        # Get self MPI communicator
        comm = c_long()
        argc = c_int(0)
        logging.debug('Calling ExaGOGetSelfCommunicator')
        self.pflowlib.ExaGOGetSelfCommunicator(byref(comm))

        logging.debug('Calling ExaGOInitialize')
        self.pflowlib.ExaGOInitialize(comm, byref(
            argc), None, "pflow".encode('ascii'), None)

        # Create PFLOW application
        self.pflow = c_longlong(0)
        logging.debug('Calling PFLOWCreate')
        ierr = self.pflowlib.PFLOWCreate(comm, byref(self.pflow))

    def read_mat_power_data(self, filename):
        logging.debug('Calling PFLOWReadMatPowerData')
        self.pflowlib.PFLOWReadMatPowerData(self.pflow, filename.encode('ascii'))

    def solve(self):
        logging.debug('Calling PFLOWSolve')
        self.pflowlib.PFLOWSolve(self.pflow)

    # Use this to not call MPI finalize when __del__ is called
    #   - This overrides in-built stack of opflow objects
    def dont_finalize(self):
        self.finalize = False

    def __del__(self):
        PFLOW.list.pop()
        logging.debug('Deleting PFLOW object')
        self.pflowlib.PFLOWDestroy(byref(self.pflow))
        if (len(PFLOW.list) == 0 and self.finalize):
          logging.debug('Calling ExaGOFinalize')
          self.pflowlib.ExaGOFinalize()
