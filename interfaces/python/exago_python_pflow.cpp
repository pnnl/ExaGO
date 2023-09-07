
#include "exago_python_pflow.hpp"

// -------------------------------------------------------------
//  class PFLOW_wrapper
// -------------------------------------------------------------
class PFLOW_wrapper {
public:
  PFLOW pflow;
  PFLOW_wrapper(void) {
    PetscErrorCode ierr;
    MPI_Comm communicator;
    ierr = ExaGOGetSelfCommunicator(&communicator);
    ExaGOCheckError(ierr);
    ierr = PFLOWCreate(communicator, &pflow);
    ExaGOCheckError(ierr);
  }
  ~PFLOW_wrapper(void) {
    PetscErrorCode ierr;
    ierr = PFLOWDestroy(&pflow);
    ExaGOCheckError(ierr);
  }
};

// -------------------------------------------------------------
// init_exago_pflow
// -------------------------------------------------------------
void init_exago_pflow(pybind11::module &m) {

  pybind11::class_<PFLOW_wrapper>(m, "PFLOW")
      .def(pybind11::init())
      .def("read_mat_power_data",
           [](PFLOW_wrapper &w, std::string filename) {
             PetscErrorCode ierr;
             ierr = PFLOWReadMatPowerData(w.pflow, filename.c_str());
             ExaGOCheckError(ierr);
           })
      .def("solve", [](PFLOW_wrapper &w) {
        PetscErrorCode ierr;
        ierr = PFLOWSolve(w.pflow);
        ExaGOCheckError(ierr);
      });
}
