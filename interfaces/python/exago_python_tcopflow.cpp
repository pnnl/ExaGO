#include "exago_python_tcopflow.hpp"

// -------------------------------------------------------------
// class TCOPFLOW_wrapper
//
// Wrap to make sure allocation and destruction is handled correctly
// -------------------------------------------------------------
class TCOPFLOW_wrapper {
public:
  TCOPFLOW tcopf;
  TCOPFLOW_wrapper(void) {
    PetscErrorCode ierr;
    MPI_Comm communicator;
    ierr = ExaGOGetSelfCommunicator(&communicator);
    ExaGOCheckError(ierr);
    ierr = TCOPFLOWCreate(communicator, &tcopf);
    ExaGOCheckError(ierr);
  }
  ~TCOPFLOW_wrapper(void) {
    PetscErrorCode ierr;
    ierr = TCOPFLOWDestroy(&tcopf);
    ExaGOCheckError(ierr);
  }
};

// -------------------------------------------------------------
// init_exago_tcopflow
// -------------------------------------------------------------
void init_exago_tcopflow(pybind11::module &m) {

  pybind11::class_<TCOPFLOW_wrapper>(m, "TCOPFLOW").def(pybind11::init())

      /* Setters */

      //   Example
      //     .def("set_model",
      //        [](SOPFLOW_wrapper &w, std::string model) {
      //          PetscErrorCode ierr;
      //          ierr = SOPFLOWSetModel(w.sopf, model.c_str());
      //          ExaGOCheckError(ierr);
      //        })
      ;
}
