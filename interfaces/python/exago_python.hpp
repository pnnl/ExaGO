

#ifndef _exago_python_hpp_
#define _exago_python_hpp_

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#ifdef EXAGO_ENABLE_MPI
#include <mpi.h>
#include <mpi4py/mpi4py.h>
#endif

#include <utils.h>
#include <string.h>

#endif
