#include "exago_python.hpp"

static inline auto prefix() -> const char * { return CMAKE_INSTALL_PREFIX; }

int initialize(char *appname, MPI_Comm comm) {
  PetscErrorCode ierr;
  ierr = ExaGOInitialize(comm, nullptr, nullptr, appname, nullptr);
  ExaGOCheckError(ierr);
  return ierr;
}

int initialize_no_comm(char *appname) {
  PetscErrorCode ierr;
  MPI_Comm communicator;
  ierr = ExaGOGetSelfCommunicator(&communicator);
  ExaGOCheckError(ierr);
  ierr = ExaGOInitialize(communicator, nullptr, nullptr, appname, nullptr);
  ExaGOCheckError(ierr);
  return ierr;
}

/*! Return a MPI communicator from mpi4py communicator object. */
MPI_Comm *get_mpi_comm(pybind11::object py_comm) {
  auto comm_ptr = PyMPIComm_Get(py_comm.ptr());

  if (!comm_ptr)
    throw pybind11::error_already_set();

  return comm_ptr;
}

PYBIND11_MODULE(exago, m) {
  // Initialize mpi4py's C-API
  if (import_mpi4py() < 0) {
    // mpi4py calls the Python C API
    // we let pybind11 give us the detailed traceback
    throw pybind11::error_already_set();
  }

  m.doc() = "Python wrapper for ExaGO.";

  /* Bindings for top-level utilities, such as initialization, finalization,
   * and helpers for interacting with enums and macros */
  // Initialize with pybind11 comm
  m.def("initialize", [](char *appname, pybind11::object py_comm) {
    auto comm = get_mpi_comm(py_comm);
    return initialize(appname, *comm);
  });
  m.def("initialize", &initialize_no_comm);
  m.def("finalize", &ExaGOFinalize);
  m.def("prefix", &prefix);

  extern void init_exago_pflow(pybind11::module & m);
  init_exago_pflow(m);

#if defined(EXAGO_ENABLE_IPOPT) || defined(EXAGO_ENABLE_HIOP)
  extern void init_exago_opflow(pybind11::module & m);
  init_exago_opflow(m);

  extern void init_exago_scopflow(pybind11::module & m);
  init_exago_scopflow(m);

  extern void init_exago_sopflow(pybind11::module & m);
  init_exago_sopflow(m);
#endif
}
