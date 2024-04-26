#include "exago_python_tcopflow.hpp"

// -------------------------------------------------------------
// class SOPFLOW_wrapper
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

  pybind11::class_<TCOPFLOW_wrapper>(m, "TCOPFLOW")
          .def(pybind11::init())
          .def("set_tolerance",
           [](TCOPFLOW_wrapper &w, double tol) {
             PetscErrorCode ierr;
             ierr = TCOPFLOWSetTolerance(w.tcopf, tol);
             ExaGOCheckError(ierr);
           })
          .def("get_tolerance",
           [](TCOPFLOW_wrapper &w) -> double {
             PetscErrorCode ierr;
             PetscReal tol;
             ierr = TCOPFLOWGetTolerance(w.tcopf, &tol);
             ExaGOCheckError(ierr);
             return tol;
           });
    }



      
    
    

    

