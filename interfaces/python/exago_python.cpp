#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <mpi.h>
#include <pflow.h>
#include <private/pflowimpl.h>
#include <opflow.h>
#include <private/opflowimpl.h>
#include <utils.h>
#include <string.h>

static inline auto prefix() -> const char * { return CMAKE_INSTALL_PREFIX; }

int initialize(char *appname) {
  PetscErrorCode ierr;
  MPI_Comm communicator;
  ierr = ExaGOGetSelfCommunicator(&communicator);
  ExaGOCheckError(ierr);
  ierr = ExaGOInitialize(communicator, nullptr, nullptr, appname, nullptr);
  ExaGOCheckError(ierr);
  return ierr;
}

PYBIND11_MODULE(exago, m) {
  m.doc() = "Python wrapper for ExaGO.";

  /* Bindings for top-level utilities, such as initialization, finalization,
   * and helpers for interacting with enums and macros */
  m.def("initialize", &initialize);
  m.def("finalize", &ExaGOFinalize);
  m.def("prefix", &prefix);

  pybind11::class_<_p_PFLOW>(m, "pflow")
      .def(pybind11::init([]() {
        PetscErrorCode ierr;
        PFLOW pflow;
        MPI_Comm communicator;
        ierr = ExaGOGetSelfCommunicator(&communicator);
        ExaGOCheckError(ierr);
        ierr = PFLOWCreate(communicator, &pflow);
        ExaGOCheckError(ierr);
        return pflow;
      }))
      .def("read_mat_power_data",
           [](_p_PFLOW &pflow, std::string filename) {
             PetscErrorCode ierr;
             ierr = PFLOWReadMatPowerData(&pflow, filename.c_str());
             ExaGOCheckError(ierr);
           })
      .def("solve", [](_p_PFLOW &pflow) {
        PetscErrorCode ierr;
        ierr = PFLOWSolve(&pflow);
        ExaGOCheckError(ierr);
      });

  pybind11::class_<_p_OPFLOW>(m, "opf")
      .def(pybind11::init([]() {
        PetscErrorCode ierr;
        OPFLOW opf;
        MPI_Comm communicator;
        ierr = ExaGOGetSelfCommunicator(&communicator);
        ExaGOCheckError(ierr);
        ierr = OPFLOWCreate(communicator, &opf);
        ExaGOCheckError(ierr);
        return opf;
      }))

      /* Setters */
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

      /* Getters */
      .def("get_tolerance",
           [](_p_OPFLOW &opf) -> double {
             PetscErrorCode ierr;
             double tol;
             ierr = OPFLOWGetTolerance(&opf, &tol);
             ExaGOCheckError(ierr);
             return tol;
           })

      .def("get_objective",
           [](_p_OPFLOW opf) -> double {
             PetscErrorCode ierr;
             double obj;
             ierr = OPFLOWGetObjective(&opf, &obj);
             ExaGOCheckError(ierr);
             return obj;
           })

      /* Everything else */
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

      .def("read_mat_power_data", [](_p_OPFLOW &opf, std::string filename) {
        PetscErrorCode ierr;
        ierr = OPFLOWReadMatPowerData(&opf, filename.c_str());
        ExaGOCheckError(ierr);
      });
}
