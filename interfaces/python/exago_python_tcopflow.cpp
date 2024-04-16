#include "exago_python_tcopflow.hpp"

// -------------------------------------------------------------
// class SOPFLOW_wrapper
//
// Wrap to make sure allocation and destruction is handled correctly
// -------------------------------------------------------------
class TCOPFLOW_wrapper {
public:
  TCOPFLOW sopf;
  TCPFLOW_wrapper(void) {
    PetscErrorCode ierr;
    MPI_Comm communicator;
    ierr = ExaGOGetSelfCommunicator(&communicator);
    ExaGOCheckError(ierr);
    ierr = TCOPFLOWCreate(communicator, &sopf);
    ExaGOCheckError(ierr);
  }
  ~TCOPFLOW_wrapper(void) {
    PetscErrorCode ierr;
    ierr = TCOPFLOWDestroy(&sopf);
    ExaGOCheckError(ierr);
  }
};

// -------------------------------------------------------------
// init_exago_tcopflow
// -------------------------------------------------------------
void init_exago_tcopflow(pybind11::module &m) {

  pybind11::class_<TCOPFLOW_wrapper>(m, "TCOPFLOW")
    .def("set_model",
           [](SOPFLOW_wrapper &w, std::string model) {
             PetscErrorCode ierr;
             ierr = TCOPFLOWSetModel(w.tcopf, model.c_str());
             ExaGOCheckError(ierr);

    .def("set_solver",
           [](TCOPFLOW_wrapper &w, std::string solver) {
             PetscErrorCode ierr;
             ierr = TCOPFLOWSetSolver(w.tcopf, solver.c_str());
             ExaGOCheckError(ierr);
           })
    .def("set_network_data",
           [](TCOPFLOW_wrapper &w, std::string filename) {
             PetscErrorCode ierr;
             ierr = TCOPFLOWSetNetworkData(w.tcopf, filename.c_str());
             ExaGOCheckError(ierr);
    .def("setup",
           [](TCOPFLOW_wrapper &w) {
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
    .def("set_load_profiles",
           [](TCOPFLOW_wrapper &w, const std::string &pload,
              const std::string &qload) {
             PetscErrorCode ierr;
             ierr =
                 TCOPFLOWSetLoadProfiles(w.tcopf, pload.c_str(), qload.c_str());
             ExaGOCheckError(ierr);
           })
    .def("set_time_step_and_duration",
           [](TCOPFLOW_wrapper &w, double dt, double tmax) {
             PetscErrorCode ierr;
             ierr = TCOPFLOWSetTimeStepandDuration(w.tcopf, dt, tmax);
             ExaGOCheckError(ierr);
           })
     .def("set_wind_gen_profile",
           [](TCOPFLOW_wrapper &w, std::string filename) {
             PetscErrorCode ierr;
             ierr = TCOPFLOWSetWindGenProfile(w.tcopf, filename.c_str());
             ExaGOCheckError(ierr);
           })
      .def("get_objective",
           [](TCOPFLOW_wrapper &w) -> double {
             PetscErrorCode ierr;
             double obj;
             ierr = TCOPFLOWGetObjective(w.tcopf, &obj);
             ExaGOCheckError(ierr);
             return obj;
           })
      .def("get_solution",
           [](TCOPFLOW_wrapper &w) {
             PetscErrorCode ierr;
             ierr = OPFLOWGETSolution(w.tcopf);
             ExaGOCheckError(ierr);
           }) 
      .def("print_solution",
           [](TCOPFLOW_wrapper &w) {
             PetscErrorCode ierr;
             ierr = OPFLOWPrintSolution(w.tcopf);
             ExaGOCheckError(ierr);
           })
       .def("save_solution",
           [](TCOPFLOW_wrapper &w, OutputFormat fmt, std::string outfile) {
             PetscErrorCode ierr;
             ierr = OPFLOWSaveSolution(w.tcopf, fmt, outfile.c_str());
             ExaGOCheckError(ierr);
           })
         .def("save_solution_all",
           [](TCCOPFLOW_wrapper &w, OutputFormat fmt, std::string outdir) {
             PetscErrorCode ierr;
             ierr = SCOPFLOWSaveSolutionAll(w.tcopf, fmt, outdir.c_str());
             ExaGOCheckError(ierr);
           })
      
        /* Getters */
          .def("get_convergence_status",
           [](TCOPFLOW_wrapper &w) -> bool {
             PetscErrorCode ierr;
             PetscBool flag;
             ierr = TCCOPFLOWGetConvergenceStatus(w.tcopf, &flag);
             ExaGOCheckError(ierr);
             return flag;
           })
           .def("get_num_iterations",
           [](TCCOPFLOW_wrapper &w) -> int {
             PetscErrorCode ierr;
             PetscInt n;
             ierr = TCOPFLOWGetNumIterations(w.tcopf, &n);
             ExaGOCheckError(ierr);
             return n;
           })
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



      
    
    

    

