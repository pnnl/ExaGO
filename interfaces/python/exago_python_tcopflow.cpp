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
          .def("set_network_data",
          [](TCOPFLOW_wrapper &w, std::string filename) {
          PetscErrorCode ierr;
          ierr = TCOPFLOWSetNetworkData(w.tcopf, filename.c_str());
          ExaGOCheckError(ierr);
          })
          .def("set_solver", [](TCOPFLOW_wrapper &w, std::string solver) {
          PetscErrorCode ierr;
          ierr = TCOPFLOWSetSolver(w.tcopf, solver.c_str());
       ExaGOCheckError(ierr);
    })
    .def("set_tolerance",
           [](TCOPFLOW_wrapper &w, double tol) {
             PetscErrorCode ierr;
             ierr = TCOPFLOWSetTolerance(w.tcopf, tol);
             ExaGOCheckError(ierr);
           })
    .def("setup", [](TCOPFLOW_wrapper &w) {
      PetscErrorCode ierr;
      ierr = TCOPFLOWSetUp(w.tcopf);
      ExaGOCheckError(ierr);
    })
        .def("solve",
           [](TCOPFLOW_wrapper &w) {
             PetscErrorCode ierr;
             ierr = TCOPFLOWSolve(w.tcopf);
             ExaGOCheckError(ierr);
           })

    .def("get_convergence_status",
         [](TCOPFLOW_wrapper &w) -> bool {
           PetscErrorCode ierr;
           PetscBool flag;
           ierr = TCOPFLOWGetConvergenceStatus(w.tcopf, &flag);
           ExaGOCheckError(ierr);
           return flag;
         })
        .def("get_objective",
             [](TCOPFLOW_wrapper &w) -> double {
               PetscErrorCode ierr;
               double obj;
               ierr = TCOPFLOWGetObjective(w.tcopf, &obj);
               ExaGOCheckError(ierr);
               return obj;
             })
        .def("get_num_iterations",
             [](TCOPFLOW_wrapper &w) -> int {
               PetscErrorCode ierr;
               PetscInt n;
               ierr = TCOPFLOWGetNumIterations(w.tcopf, &n);
               ExaGOCheckError(ierr);
               return n;
             });
       /*.def("save_solution",
             [](TCOPFLOW_wrapper &w, OutputFormat fmt, std::string outfile) {
               PetscErrorCode ierr;
               ierr = TCOPFLOWSaveSolution(w.tcopf, fmt, outfile.c_str());
               ExaGOCheckError(ierr);
             })
    .def("save_solution_all",
         [](TCOPFLOW_wrapper &w, OutputFormat fmt, std::string outdir) {
           PetscErrorCode ierr;
           ierr = TCCOPFLOWSaveSolutionAll(w.tcopf, fmt, outdir.c_str());
           ExaGOCheckError(ierr);
         })*/
    
    }
