#include <pybind11/pybind11.h>
#include <mpi.h>
#include <opflow.h>
#include <private/opflowimpl.h>
#include <utils.hpp>
#include <string.h>

std::string get_datafile_path() {
  return std::string(EXAGO_OPTIONS_DIR) + "/../datafiles/";
}

int initialize(char *appname) {
  PetscErrorCode ierr;
  MPI_Comm communicator;
  ierr = ExaGOGetSelfCommunicator(&communicator);
  ExaGOCheckError(ierr);
  ierr = ExaGOInitialize(communicator, NULL, NULL, appname, NULL);
  ExaGOCheckError(ierr);
  return ierr;
}

PYBIND11_MODULE(exago, m) {
  m.doc() = "Python wrapper for ExaGO.";
  m.def("initialize", &initialize);
  m.def("finalize", &ExaGOFinalize);
  m.def("get_datafile_path", &get_datafile_path);

  pybind11::class_<_p_OPFLOW>(m, "opf")
      .def(pybind11::init([]() {
        PetscErrorCode ierr;
        OPFLOW opf;
        MPI_Comm communicator;
        ierr = ExaGOGetSelfCommunicator(&communicator);
        ExaGOCheckError(ierr);
        char *appname = (char *)"opflow";
        ierr = OPFLOWCreate(communicator, &opf);
        ExaGOCheckError(ierr);
        return opf;
      }))

      .def("solve",
           [](_p_OPFLOW &opf) {
             PetscErrorCode ierr;
             ierr = OPFLOWSolve(&opf);
             ExaGOCheckError(ierr);
           })

      .def("print_solution",
           [](_p_OPFLOW &opf) {
             PetscErrorCode ierr;
             ierr = OPFLOWPrintSolution(&opf);
             ExaGOCheckError(ierr);
           })

      .def("set_loadloss_penalty",
           [](_p_OPFLOW &opf, double p) {
             PetscErrorCode ierr;
             ierr = OPFLOWSetLoadLossPenalty(&opf, p);
             ExaGOCheckError(ierr);
           })

      .def("set_powerimbalance_penalty",
           [](_p_OPFLOW &opf, double p) {
             PetscErrorCode ierr;
             ierr = OPFLOWSetBusPowerImbalancePenalty(&opf, p);
             ExaGOCheckError(ierr);
           })

      .def("set_tolerance",
           [](_p_OPFLOW &opf, double tol) {
             PetscErrorCode ierr;
             ierr = OPFLOWSetTolerance(&opf, tol);
             ExaGOCheckError(ierr);
           })

      .def("get_tolerance",
           [](_p_OPFLOW &opf) -> double {
             PetscErrorCode ierr;
             double tol;
             ierr = OPFLOWGetTolerance(&opf, &tol);
             ExaGOCheckError(ierr);
             return tol;
           })

      .def("read_mat_power_data", [](_p_OPFLOW &opf, std::string filename) {
        PetscErrorCode ierr;
        ierr = OPFLOWReadMatPowerData(&opf, filename.c_str());
        ExaGOCheckError(ierr);
      });
}
