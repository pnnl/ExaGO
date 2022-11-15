#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <mpi.h>
#include <mpi4py/mpi4py.h>
#include <pflow.h>
#include <private/pflowimpl.h>
#include <opflow.h>
#include <private/opflowimpl.h>
#include <sopflow.h>
#include <private/sopflowimpl.h>
#include <scopflow.h>
#include <private/scopflowimpl.h>
#include <sopflow.h>
#include <private/sopflowimpl.h>
#include <utils.h>
#include <string.h>

/*
 * Calculate the number of types in the opflow type enums. An immediately
 * invoked lambda is used to keep the final variable const, even thought we have
 * to perform some operations to calculate that value.
 */
constexpr int num_objective_types = 0
#define OPFLOW_OBJ_DEF(x) +1
#include "private/opflow_enum.def"
#undef OPFLOW_OBJ_DEF
    ;

constexpr int num_initialization_types = 0
#define OPFLOW_INIT_DEF(x) +1
#include "private/opflow_enum.def"
#undef OPFLOW_INIT_DEF
    ;

constexpr int num_gen_bus_voltage_types = 0
#define OPFLOW_GBV(x) +1
#include "private/opflow_enum.def"
#undef OPFLOW_GBV_DEF
    ;

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

  // -------------------------------------------------------------
  // Common Types
  // -------------------------------------------------------------

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

  /* Bindings for OPFLOW enums */

#define OPFLOW_OBJ_DEF(x) .value(#x, x)
  pybind11::enum_<OPFLOWObjectiveType>(m, "OPFLOWObjectiveType")
#include "private/opflow_enum.def"
      .export_values();
#undef OPFLOW_OBJ_DEF

#define OPFLOW_INIT_DEF(x) .value(#x, x)
  pybind11::enum_<OPFLOWInitializationType>(m, "OPFLOWInitializationType")
#include "private/opflow_enum.def"
      .export_values();
#undef OPFLOW_INIT_DEF

#define OPFLOW_GBV_DEF(x) .value(#x, x)
  pybind11::enum_<OPFLOWGenBusVoltageType>(m, "OPFLOWGenBusVoltageType")
#include "private/opflow_enum.def"
      .export_values();
#undef OPFLOW_GBV_DEF

  pybind11::enum_<OutputFormat>(m, "OutputFormat")
      .value("CSV", CSV)
      .value("MATPOWER", MATPOWER)
      .export_values();

  // -------------------------------------------------------------
  //  class OPFLOW_wrapper
  // -------------------------------------------------------------
  class OPFLOW_wrapper {
  public:
    OPFLOW opf;
    OPFLOW_wrapper(void) {
      PetscErrorCode ierr;
      MPI_Comm communicator;
      ierr = ExaGOGetSelfCommunicator(&communicator);
      ExaGOCheckError(ierr);
      ierr = OPFLOWCreate(communicator, &opf);
      ExaGOCheckError(ierr);
    }
    ~OPFLOW_wrapper(void) {
      PetscErrorCode ierr;
      ierr = OPFLOWDestroy(&opf);
      ExaGOCheckError(ierr);
    }
  };

  /* OPFLOW class */
  pybind11::class_<OPFLOW_wrapper>(m, "OPFLOW")
      .def(pybind11::init())
      /* Setters */
      .def("set_model",
           [](OPFLOW_wrapper &w, std::string model) {
             PetscErrorCode ierr;
             ierr = OPFLOWSetModel(w.opf, model.c_str());
             ExaGOCheckError(ierr);
           })

      .def("set_solver",
           [](OPFLOW_wrapper &w, std::string solver) {
             PetscErrorCode ierr;
             ierr = OPFLOWSetSolver(w.opf, solver.c_str());
             ExaGOCheckError(ierr);
           })

      .def("set_objective_type",
           [](OPFLOW_wrapper &w, std::string func) {
             PetscErrorCode ierr;
             OPFLOWObjectiveType type;
             bool found = false;
// boolean for found, throw error if not
#define OPFLOW_OBJ_DEF(x)                                                      \
  if (func == #x) {                                                            \
    found = true;                                                              \
    type = x;                                                                  \
  }
#include "private/opflow_enum.def"
#undef OPFLOW_OBJ_DEF
             if (!found) {
               throw ExaGOError(
                   fmt::format(
                       "Invalid value '{}' for OPFLOWObjectiveType enum", func)
                       .c_str());
             }
             ierr = OPFLOWSetObjectiveType(w.opf, type);
             ExaGOCheckError(ierr);
           })

      .def("set_objective_type",
           [](OPFLOW_wrapper &w, int func) {
             PetscErrorCode ierr;
             OPFLOWObjectiveType type;
             bool found = false;
#define OPFLOW_OBJ_DEF(x)                                                      \
  if (func == x) {                                                             \
    found = true;                                                              \
    type = x;                                                                  \
  }
#include "private/opflow_enum.def"
#undef OPFLOW_OBJ_DEF
             if (!found) {
               throw ExaGOError(
                   fmt::format(
                       "Invalid value '{}' for OPFLOWObjectiveType enum", func)
                       .c_str());
             }
             ierr = OPFLOWSetObjectiveType(w.opf, type);
             ExaGOCheckError(ierr);
           })

      .def("set_initialization_type",
           [](OPFLOW_wrapper &w, std::string init) {
             PetscErrorCode ierr;
             OPFLOWInitializationType type;
             bool found = false;
#define OPFLOW_INIT_DEF(x)                                                     \
  if (init == #x) {                                                            \
    found = true;                                                              \
    type = x;                                                                  \
  }
#include "private/opflow_enum.def"
#undef OPFLOW_INIT_DEF
             if (!found) {
               throw ExaGOError(
                   fmt::format(
                       "Invalid value '{}' for OPFLOWInitializationType enum",
                       init)
                       .c_str());
             }
             ierr = OPFLOWSetInitializationType(w.opf, type);
             ExaGOCheckError(ierr);
           })

      .def("set_initialization_type",
           [](OPFLOW_wrapper &w, int init) {
             PetscErrorCode ierr;
             OPFLOWInitializationType type;
             bool found = false;
#define OPFLOW_INIT_DEF(x)                                                     \
  if (init == x) {                                                             \
    found = true;                                                              \
    type = x;                                                                  \
  }
#include "private/opflow_enum.def"
#undef OPFLOW_INIT_DEF
             if (!found) {
               throw ExaGOError(
                   fmt::format(
                       "Invalid value '{}' for OPFLOWInitializationType enum",
                       init)
                       .c_str());
             }
             ierr = OPFLOWSetInitializationType(w.opf, type);
             ExaGOCheckError(ierr);
           })

      .def("set_gen_bus_voltage_type",
           [](OPFLOW_wrapper &w, std::string volt) {
             PetscErrorCode ierr;
             OPFLOWGenBusVoltageType type;
             bool found = false;
#define OPFLOW_GBV_DEF(x)                                                      \
  if (volt == #x) {                                                            \
    found = true;                                                              \
    type = x;                                                                  \
  }
#include "private/opflow_enum.def"
#undef OPFLOW_GBV_DEF
             if (!found) {
               throw ExaGOError(
                   fmt::format(
                       "Invalid value '{}' for OPFLOWGenBusVoltageType enum",
                       volt)
                       .c_str());
             }
             ierr = OPFLOWSetGenBusVoltageType(w.opf, type);
             ExaGOCheckError(ierr);
           })

      .def("set_gen_bus_voltage_type",
           [](OPFLOW_wrapper &w, int volt) {
             PetscErrorCode ierr;
             OPFLOWGenBusVoltageType type;
             bool found = false;
#define OPFLOW_GBV_DEF(x)                                                      \
  if (volt == x) {                                                             \
    found = true;                                                              \
    type = x;                                                                  \
  }
#include "private/opflow_enum.def"
#undef OPFLOW_GBV_DEF
             if (!found) {
               throw ExaGOError(
                   fmt::format(
                       "Invalid value '{}' for OPFLOWGenBusVoltageType enum",
                       volt)
                       .c_str());
             }
             ierr = OPFLOWSetGenBusVoltageType(w.opf, type);
             ExaGOCheckError(ierr);
           })

      .def("set_has_gen_set_point",
           [](OPFLOW_wrapper &w, bool setpt) {
             PetscErrorCode ierr;
             ierr = OPFLOWHasGenSetPoint(w.opf, (PetscBool)setpt);
             ExaGOCheckError(ierr);
           })

      .def("set_hiop_compute_mode",
           [](OPFLOW_wrapper &w, std::string mode) {
             PetscErrorCode ierr;
             ierr = OPFLOWSetHIOPComputeMode(w.opf, mode.c_str());
             ExaGOCheckError(ierr);
           })

      .def("set_has_loadloss",
           [](OPFLOW_wrapper &w, bool loadloss) {
             PetscErrorCode ierr;
             ierr = OPFLOWHasLoadLoss(w.opf, (PetscBool)loadloss);
             ExaGOCheckError(ierr);
           })

      .def("set_ignore_lineflow_constraints",
           [](OPFLOW_wrapper &w, bool ignore) {
             PetscErrorCode ierr;
             ierr = OPFLOWIgnoreLineflowConstraints(w.opf, (PetscBool)ignore);
             ExaGOCheckError(ierr);
           })

      .def("set_has_bus_power_imbalance",
           [](OPFLOW_wrapper &w, bool powerimbalance) {
             PetscErrorCode ierr;
             ierr =
                 OPFLOWHasBusPowerImbalance(w.opf, (PetscBool)powerimbalance);
             ExaGOCheckError(ierr);
           })

      .def("set_use_agc",
           [](OPFLOW_wrapper &w, bool agc) {
             PetscErrorCode ierr;
             ierr = OPFLOWUseAGC(w.opf, (PetscBool)agc);
             ExaGOCheckError(ierr);
           })

      .def("set_hiop_verbosity_level",
           [](OPFLOW_wrapper &w, int level) {
             PetscErrorCode ierr;
             ierr = OPFLOWSetHIOPVerbosityLevel(w.opf, level);
             ExaGOCheckError(ierr);
           })

      .def("set_loadloss_penalty",
           [](OPFLOW_wrapper &w, double p) {
             PetscErrorCode ierr;
             ierr = OPFLOWSetLoadLossPenalty(w.opf, p);
             ExaGOCheckError(ierr);
           })

      .def("set_bus_power_imbalance_penalty",
           [](OPFLOW_wrapper &w, double p) {
             PetscErrorCode ierr;
             ierr = OPFLOWSetBusPowerImbalancePenalty(w.opf, p);
             ExaGOCheckError(ierr);
           })

      .def("set_tolerance",
           [](OPFLOW_wrapper &w, double tol) {
             PetscErrorCode ierr;
             ierr = OPFLOWSetTolerance(w.opf, tol);
             ExaGOCheckError(ierr);
           })

      .def("ps_set_gen_power_limits",
           [](OPFLOW_wrapper &w, int gbus, std::string gid, double pt,
              double pb, double qt, double qb) {
             PetscErrorCode ierr;
             _p_PS *ps;
             ierr = OPFLOWGetPS(w.opf, &ps);
             ExaGOCheckError(ierr);
             PSSetGenPowerLimits(ps, gbus, gid.c_str(), pt, pb, qt, qb);
           })

      /* Getters */
      .def("get_tolerance",
           [](OPFLOW_wrapper &w) -> double {
             PetscErrorCode ierr;
             double tol;
             ierr = OPFLOWGetTolerance(w.opf, &tol);
             ExaGOCheckError(ierr);
             return tol;
           })

      .def("get_hiop_compute_mode",
           [](OPFLOW_wrapper &w) -> std::string {
             PetscErrorCode ierr;
             std::string mode(PETSC_MAX_PATH_LEN, '\0');
             ierr = OPFLOWGetHIOPComputeMode(w.opf, &mode[0]);
             ExaGOCheckError(ierr);
             return mode.c_str();
           })

      .def("get_model",
           [](OPFLOW_wrapper &w) -> std::string {
             PetscErrorCode ierr;
             std::string model(PETSC_MAX_PATH_LEN, '\0');
             ierr = OPFLOWGetModel(w.opf, &model[0]);
             ExaGOCheckError(ierr);
             return model.c_str();
           })

      .def("get_solver",
           [](OPFLOW_wrapper &w) -> std::string {
             PetscErrorCode ierr;
             std::string solver(PETSC_MAX_PATH_LEN, '\0');
             ierr = OPFLOWGetSolver(w.opf, &solver[0]);
             ExaGOCheckError(ierr);
             return solver.c_str();
           })

      .def("get_convergence_status",
           [](OPFLOW_wrapper &w) -> bool {
             PetscErrorCode ierr;
             PetscBool status;
             ierr = OPFLOWGetConvergenceStatus(w.opf, &status);
             ExaGOCheckError(ierr);
             return status;
           })

      .def("get_objective_type",
           [](OPFLOW_wrapper &w) -> OPFLOWObjectiveType {
             PetscErrorCode ierr;
             OPFLOWObjectiveType obj_type;
             ierr = OPFLOWGetObjectiveType(w.opf, &obj_type);
             ExaGOCheckError(ierr);
             return obj_type;
           })

      .def("get_initialization_type",
           [](OPFLOW_wrapper &w) -> OPFLOWInitializationType {
             PetscErrorCode ierr;
             OPFLOWInitializationType init;
             ierr = OPFLOWGetInitializationType(w.opf, &init);
             ExaGOCheckError(ierr);
             return init;
           })

      .def("get_gen_bus_voltage_type",
           [](OPFLOW_wrapper &w) -> OPFLOWGenBusVoltageType {
             PetscErrorCode ierr;
             OPFLOWGenBusVoltageType volt;
             ierr = OPFLOWGetGenBusVoltageType(w.opf, &volt);
             ExaGOCheckError(ierr);
             return volt;
           })

      .def("get_has_gen_set_point",
           [](OPFLOW_wrapper &w) -> bool {
             PetscErrorCode ierr;
             PetscBool setpt;
             ierr = OPFLOWGetHasGenSetPoint(w.opf, &setpt);
             ExaGOCheckError(ierr);
             return setpt;
           })

      .def("get_loadloss_penalty",
           [](OPFLOW_wrapper &w) -> double {
             PetscErrorCode ierr;
             double p;
             ierr = OPFLOWGetLoadLossPenalty(w.opf, &p);
             ExaGOCheckError(ierr);
             return p;
           })

      .def("get_ignore_lineflow_constraints",
           [](OPFLOW_wrapper &w) -> bool {
             PetscErrorCode ierr;
             PetscBool ig;
             ierr = OPFLOWGetIgnoreLineflowConstraints(w.opf, &ig);
             ExaGOCheckError(ierr);
             return ig;
           })

      .def("get_has_loadloss",
           [](OPFLOW_wrapper &w) -> bool {
             PetscErrorCode ierr;
             PetscBool loadloss;
             ierr = OPFLOWGetHasLoadLoss(w.opf, &loadloss);
             ExaGOCheckError(ierr);
             return loadloss;
           })

      .def("get_has_bus_power_imbalance",
           [](OPFLOW_wrapper &w) -> bool {
             PetscErrorCode ierr;
             PetscBool powerimbalance;
             ierr = OPFLOWGetHasBusPowerImbalance(w.opf, &powerimbalance);
             ExaGOCheckError(ierr);
             return powerimbalance;
           })

      .def("get_use_agc",
           [](OPFLOW_wrapper &w) -> bool {
             PetscErrorCode ierr;
             PetscBool agc;
             ierr = OPFLOWGetUseAGC(w.opf, &agc);
             ExaGOCheckError(ierr);
             return agc;
           })

      .def("get_hiop_verbosity_level",
           [](OPFLOW_wrapper &w) -> int {
             PetscErrorCode ierr;
             int level;
             ierr = OPFLOWGetHIOPVerbosityLevel(w.opf, &level);
             ExaGOCheckError(ierr);
             return level;
           })

      .def("get_bus_power_imbalance_penalty",
           [](OPFLOW_wrapper &w) -> double {
             PetscErrorCode ierr;
             double p;
             ierr = OPFLOWGetBusPowerImbalancePenalty(w.opf, &p);
             ExaGOCheckError(ierr);
             return p;
           })

      .def("get_gen_dispatch",
           [](OPFLOW_wrapper &w, int gbus, std::string gid) -> pybind11::tuple {
             _p_PS *ps;
             PetscErrorCode ierr;
             ierr = OPFLOWGetPS(w.opf, &ps);
             ExaGOCheckError(ierr);
             double pg;
             double qg;
             ierr = PSGetGenDispatch(ps, gbus, gid.c_str(), &pg, &qg);
             ExaGOCheckError(ierr);
             return pybind11::make_tuple(pg, qg);
           })

      .def("get_objective_types",
           [](OPFLOW_wrapper &w) -> std::vector<OPFLOWObjectiveType> {
             std::vector<OPFLOWObjectiveType> types;
#define OPFLOW_OBJ_DEF(x) types.push_back(x);
#include "private/opflow_enum.def"
#undef OPFLOW_OBJ_DEF
             return types;
           })

      .def("get_initialization_types",
           [](OPFLOW_wrapper &w) -> std::vector<OPFLOWInitializationType> {
             std::vector<OPFLOWInitializationType> types;
#define OPFLOW_INIT_DEF(x) types.push_back(x);
#include "private/opflow_enum.def"
#undef OPFLOW_INIT_DEF
             return types;
           })

      .def("get_gen_bus_voltage_types",
           [](OPFLOW_wrapper &w) -> std::vector<OPFLOWGenBusVoltageType> {
             std::vector<OPFLOWGenBusVoltageType> types;
#define OPFLOW_GBV_DEF(x) types.push_back(x);
#include "private/opflow_enum.def"
#undef OPFLOW_GBV_DEF
             return types;
           })

      .def("get_objective",
           [](OPFLOW_wrapper &w) -> double {
             PetscErrorCode ierr;
             double obj;
             ierr = OPFLOWGetObjective(w.opf, &obj);
             ExaGOCheckError(ierr);
             return obj;
           })

      /* Everything else */
      .def("solve",
           [](OPFLOW_wrapper &w) {
             PetscErrorCode ierr;
             ierr = OPFLOWSolve(w.opf);
             ExaGOCheckError(ierr);
           })

      .def("print_solution",
           [](OPFLOW_wrapper &w) {
             PetscErrorCode ierr;
             ierr = OPFLOWPrintSolution(w.opf);
             ExaGOCheckError(ierr);
           })

      .def("solution_to_ps",
           [](OPFLOW_wrapper &w) {
             PetscErrorCode ierr;
             ierr = OPFLOWSolutionToPS(w.opf);
             ExaGOCheckError(ierr);
           })

      .def("set_up_ps",
           [](OPFLOW_wrapper &w) {
             PetscErrorCode ierr;
             ierr = OPFLOWSetUpPS(w.opf);
             ExaGOCheckError(ierr);
           })

      .def("save_solution",
           [](OPFLOW_wrapper &w, OutputFormat fmt, std::string outfile) {
             PetscErrorCode ierr;
             ierr = OPFLOWSaveSolution(w.opf, fmt, outfile.c_str());
             ExaGOCheckError(ierr);
           })

      .def("read_mat_power_data", [](OPFLOW_wrapper &w, std::string filename) {
        PetscErrorCode ierr;
        ierr = OPFLOWReadMatPowerData(w.opf, filename.c_str());
        ExaGOCheckError(ierr);
      });

  // -------------------------------------------------------------
  /* SCOPFLOW class */
  // -------------------------------------------------------------

  // -------------------------------------------------------------
  // class SCOPFLOW_wrapper
  //
  // Wrap to make sure allocation and destruction is handled correctly
  // -------------------------------------------------------------
  class SCOPFLOW_wrapper {
  public:
    SCOPFLOW scopf;
    SCOPFLOW_wrapper(void) {
      PetscErrorCode ierr;
      MPI_Comm communicator;
      ierr = ExaGOGetSelfCommunicator(&communicator);
      ExaGOCheckError(ierr);
      ierr = SCOPFLOWCreate(communicator, &scopf);
      ExaGOCheckError(ierr);
    }
    ~SCOPFLOW_wrapper(void) {
      PetscErrorCode ierr;
      ierr = SCOPFLOWDestroy(&scopf);
      ExaGOCheckError(ierr);
    }
  };

  pybind11::enum_<ContingencyFileInputFormat>(m, "ContingencyFileInputFormat")
      .value("NATIVE", NATIVE)
      .value("PSSE", PSSE)
      .export_values();

  pybind11::class_<SCOPFLOW_wrapper>(m, "SCOPFLOW")
      .def(pybind11::init())
      /* Setters */
      .def("set_model",
           [](SCOPFLOW_wrapper &w, std::string model) {
             PetscErrorCode ierr;
             ierr = SCOPFLOWSetModel(w.scopf, model.c_str());
             ExaGOCheckError(ierr);
           })
      .def("set_network_data",
           [](SCOPFLOW_wrapper &w, std::string filename) {
             PetscErrorCode ierr;
             ierr = SCOPFLOWSetNetworkData(w.scopf, filename.c_str());
             ExaGOCheckError(ierr);
           })

      .def("set_num_contingencies",
           [](SCOPFLOW_wrapper &w, int n) {
             PetscErrorCode ierr;
             ierr = SCOPFLOWSetNumContingencies(w.scopf, n);
             ExaGOCheckError(ierr);
           })
      .def("set_contingency_data",
           [](SCOPFLOW_wrapper &w, std::string filename,
              ContingencyFileInputFormat fmt) {
             PetscErrorCode ierr;
             ierr = SCOPFLOWSetContingencyData(w.scopf, fmt, filename.c_str());
             ExaGOCheckError(ierr);
           })
      .def("set_pload_data",
           [](SCOPFLOW_wrapper &w, std::string filename) {
             PetscErrorCode ierr;
             ierr = SCOPFLOWSetPLoadData(w.scopf, filename.c_str());
             ExaGOCheckError(ierr);
           })
      .def("set_qload_data",
           [](SCOPFLOW_wrapper &w, std::string filename) {
             PetscErrorCode ierr;
             ierr = SCOPFLOWSetQLoadData(w.scopf, filename.c_str());
             ExaGOCheckError(ierr);
           })
      .def("set_wind_gen_profile",
           [](SCOPFLOW_wrapper &w, std::string filename) {
             PetscErrorCode ierr;
             ierr = SCOPFLOWSetWindGenProfile(w.scopf, filename.c_str());
             ExaGOCheckError(ierr);
           })
      .def("set_time_step",
           [](SCOPFLOW_wrapper &w, double dt) {
             PetscErrorCode ierr;
             ierr = SCOPFLOWSetTimeStep(w.scopf, dt);
             ExaGOCheckError(ierr);
           })
      .def("set_duration",
           [](SCOPFLOW_wrapper &w, double d) {
             PetscErrorCode ierr;
             ierr = SCOPFLOWSetDuration(w.scopf, d);
             ExaGOCheckError(ierr);
           })
      .def("set_tolerance",
           [](SCOPFLOW_wrapper &w, double tol) {
             PetscErrorCode ierr;
             ierr = SCOPFLOWSetTolerance(w.scopf, tol);
             ExaGOCheckError(ierr);
           })
      .def("set_verbosity_level",
           [](SCOPFLOW_wrapper &w, int level) {
             PetscErrorCode ierr;
             ierr = SCOPFLOWSetVerbosityLevel(w.scopf, level);
             ExaGOCheckError(ierr);
           })
      .def("set_compute_mode",
           [](SCOPFLOW_wrapper &w, std::string mode) {
             PetscErrorCode ierr;
             ierr = SCOPFLOWSetComputeMode(w.scopf, mode.c_str());
             ExaGOCheckError(ierr);
           })
      .def("set_mode",
           [](SCOPFLOW_wrapper &w, int mode) {
             PetscErrorCode ierr;
             ierr = SCOPFLOWSetMode(w.scopf, mode);
             ExaGOCheckError(ierr);
           })
      .def("set_solver",
           [](SCOPFLOW_wrapper &w, std::string solver) {
             PetscErrorCode ierr;
             ierr = SCOPFLOWSetSolver(w.scopf, solver.c_str());
             ExaGOCheckError(ierr);
           })
      .def("set_subproblem_model",
           [](SCOPFLOW_wrapper &w, std::string model) {
             PetscErrorCode ierr;
             ierr = SCOPFLOWSetSubproblemModel(w.scopf, model.c_str());
             ExaGOCheckError(ierr);
           })
      .def("set_subproblem_solver",
           [](SCOPFLOW_wrapper &w, std::string solver) {
             PetscErrorCode ierr;
             ierr = SCOPFLOWSetSubproblemSolver(w.scopf, solver.c_str());
             ExaGOCheckError(ierr);
           })
      .def("set_initialization_type",
           [](SCOPFLOW_wrapper &w, std::string init) {
             PetscErrorCode ierr;
             OPFLOWInitializationType type;
             bool found = false;
#define OPFLOW_INIT_DEF(x)                                                     \
  if (init == #x) {                                                            \
    found = true;                                                              \
    type = x;                                                                  \
  }
#include "private/opflow_enum.def"
#undef OPFLOW_INIT_DEF
             if (!found) {
               throw ExaGOError(
                   fmt::format(
                       "Invalid value '{}' for OPFLOWInitializationType enum",
                       init)
                       .c_str());
             }
             ierr = SCOPFLOWSetInitilizationType(w.scopf, type);
             ExaGOCheckError(ierr);
           })
      .def("set_initialization_type",
           [](SCOPFLOW_wrapper &w, int init) {
             PetscErrorCode ierr;
             OPFLOWInitializationType type;
             bool found = false;
#define OPFLOW_INIT_DEF(x)                                                     \
  if (init == x) {                                                             \
    found = true;                                                              \
    type = x;                                                                  \
  }
#include "private/opflow_enum.def"
#undef OPFLOW_INIT_DEF
             if (!found) {
               throw ExaGOError(
                   fmt::format(
                       "Invalid value '{}' for SCOPFLOWInitializationType enum",
                       init)
                       .c_str());
             }
             ierr = SCOPFLOWSetInitilizationType(w.scopf, type);
             ExaGOCheckError(ierr);
           })

      .def("set_gen_bus_voltage_type",
           [](SCOPFLOW_wrapper &w, int volt) {
             PetscErrorCode ierr;
             OPFLOWGenBusVoltageType type;
             bool found = false;
#define OPFLOW_GBV_DEF(x)                                                      \
  if (volt == x) {                                                             \
    found = true;                                                              \
    type = x;                                                                  \
  }
#include "private/opflow_enum.def"
#undef OPFLOW_GBV_DEF
             if (!found) {
               throw ExaGOError(
                   fmt::format(
                       "Invalid value '{}' for SCOPFLOWGenBusVoltageType enum",
                       volt)
                       .c_str());
             }
             ierr = SCOPFLOWSetGenBusVoltageType(w.scopf, type);
             ExaGOCheckError(ierr);
           })

      .def("set_gen_bus_voltage_type",
           [](SCOPFLOW_wrapper &w, std::string volt) {
             PetscErrorCode ierr;
             OPFLOWGenBusVoltageType type;
             bool found = false;
#define OPFLOW_GBV_DEF(x)                                                      \
  if (volt == #x) {                                                            \
    found = true;                                                              \
    type = x;                                                                  \
  }
#include "private/opflow_enum.def"
#undef OPFLOW_GBV_DEF
             if (!found) {
               throw ExaGOError(
                   fmt::format(
                       "Invalid value '{}' for SCOPFLOWGenBusVoltageType enum",
                       volt)
                       .c_str());
             }
             ierr = SCOPFLOWSetGenBusVoltageType(w.scopf, type);
             ExaGOCheckError(ierr);
           })

      .def("enable_multi_period",
           [](SCOPFLOW_wrapper &w, bool flag) {
             PetscErrorCode ierr;
             ierr = SCOPFLOWEnableMultiPeriod(w.scopf, (PetscBool)flag);
             ExaGOCheckError(ierr);
           })
      .def("enable_power_imbalance_variables",
           [](SCOPFLOW_wrapper &w, bool flag) {
             PetscErrorCode ierr;
             ierr = SCOPFLOWEnablePowerImbalanceVariables(w.scopf,
                                                          (PetscBool)flag);
             ExaGOCheckError(ierr);
           })
      .def("ignore_lineflow_constraints",
           [](SCOPFLOW_wrapper &w, bool flag) {
             PetscErrorCode ierr;
             ierr = SCOPFLOWIgnoreLineflowConstraints(w.scopf, (PetscBool)flag);
             ExaGOCheckError(ierr);
           })

      /* Getters */
      .def("get_tolerance",
           [](SCOPFLOW_wrapper &w) -> double {
             PetscErrorCode ierr;
             PetscReal tol;
             ierr = SCOPFLOWGetTolerance(w.scopf, &tol);
             ExaGOCheckError(ierr);
             return tol;
           })
      .def("get_num_iterations",
           [](SCOPFLOW_wrapper &w) -> int {
             PetscErrorCode ierr;
             PetscInt n;
             ierr = SCOPFLOWGetNumIterations(w.scopf, &n);
             ExaGOCheckError(ierr);
             return n;
           })
      .def("get_convergence_status",
           [](SCOPFLOW_wrapper &w) -> bool {
             PetscErrorCode ierr;
             PetscBool flag;
             ierr = SCOPFLOWGetConvergenceStatus(w.scopf, &flag);
             ExaGOCheckError(ierr);
             return flag;
           })
      .def("get_total_objective",
           [](SCOPFLOW_wrapper &w) -> double {
             PetscErrorCode ierr;
             PetscReal obj;
             ierr = SCOPFLOWGetTotalObjective(w.scopf, &obj);
             ExaGOCheckError(ierr);
             return obj;
           })
      .def("get_base_objective",
           [](SCOPFLOW_wrapper &w) -> double {
             PetscErrorCode ierr;
             PetscReal obj;
             ierr = SCOPFLOWGetBaseObjective(w.scopf, &obj);
             ExaGOCheckError(ierr);
             return obj;
           })

      /* Solution */
      .def("set_up",
           [](SCOPFLOW_wrapper &w) {
             PetscErrorCode ierr;
             ierr = SCOPFLOWSetUp(w.scopf);
             ExaGOCheckError(ierr);
           })
      .def("solve",
           [](SCOPFLOW_wrapper &w) {
             PetscErrorCode ierr;
             ierr = SCOPFLOWSolve(w.scopf);
             ExaGOCheckError(ierr);
           })

      .def("print_solution",
           [](SCOPFLOW_wrapper &w, int n) {
             PetscErrorCode ierr;
             ierr = SCOPFLOWPrintSolution(w.scopf, n);
             ExaGOCheckError(ierr);
           })

      .def("save_solution",
           [](SCOPFLOW_wrapper &w, int n, OutputFormat fmt,
              std::string outfile) {
             PetscErrorCode ierr;
             ierr = SCOPFLOWSaveSolution(w.scopf, n, fmt, outfile.c_str());
             ExaGOCheckError(ierr);
           })

      ;

  // -------------------------------------------------------------
  /* SOPFLOW class */
  // -------------------------------------------------------------

  // -------------------------------------------------------------
  // class SOPFLOW_wrapper
  //
  // Wrap to make sure allocation and destruction is handled correctly
  // -------------------------------------------------------------
  class SOPFLOW_wrapper {
  public:
    SOPFLOW sopf;
    SOPFLOW_wrapper(void) {
      PetscErrorCode ierr;
      MPI_Comm communicator;
      ierr = ExaGOGetSelfCommunicator(&communicator);
      ExaGOCheckError(ierr);
      ierr = SOPFLOWCreate(communicator, &sopf);
      ExaGOCheckError(ierr);
    }
    ~SOPFLOW_wrapper(void) {
      PetscErrorCode ierr;
      ierr = SOPFLOWDestroy(&sopf);
      ExaGOCheckError(ierr);
    }
  };

  pybind11::enum_<ScenarioFileInputFormat>(m, "ScenarioFileInputFormat")
      .value("NATIVE_SINGLEPERIOD", SOPFLOW_NATIVE_SINGLEPERIOD)
      .value("NATIVE_MULTIPERIOD", SOPFLOW_NATIVE_MULTIPERIOD)
      .export_values();

  pybind11::enum_<ScenarioUncertaintyType>(m, "ScenarioUncertaintyType")
      .value("NONE", NONE)
      .value("WIND", WIND)
      .value("LOAD", LOAD)
      .export_values();

  pybind11::class_<SOPFLOW_wrapper>(m, "SOPFLOW")
      .def(pybind11::init())

      /* Setters */

      .def("set_model",
           [](SOPFLOW_wrapper &w, std::string model) {
             PetscErrorCode ierr;
             ierr = SOPFLOWSetModel(w.sopf, model.c_str());
             ExaGOCheckError(ierr);
           })
      .def("set_network_data",
           [](SOPFLOW_wrapper &w, std::string filename) {
             PetscErrorCode ierr;
             ierr = SOPFLOWSetNetworkData(w.sopf, filename.c_str());
             ExaGOCheckError(ierr);
           })
      .def("set_contingency_data",
           [](SOPFLOW_wrapper &w, std::string filename,
              ContingencyFileInputFormat fmt = NATIVE) {
             PetscErrorCode ierr;
             ierr = SOPFLOWSetContingencyData(w.sopf, fmt, filename.c_str());
             ExaGOCheckError(ierr);
           })
      .def("set_num_contingencies",
           [](SOPFLOW_wrapper &w, int n) {
             PetscErrorCode ierr;
             ierr = SOPFLOWSetNumContingencies(w.sopf, n);
             ExaGOCheckError(ierr);
           })
      .def("set_scenario_data",
           [](SOPFLOW_wrapper &w, std::string filename,
              ScenarioFileInputFormat sfmt, ScenarioUncertaintyType ut) {
             PetscErrorCode ierr;
             ierr = SOPFLOWSetScenarioData(w.sopf, sfmt, ut, filename.c_str());
             ExaGOCheckError(ierr);
           })
      .def("set_num_scenarios",
           [](SOPFLOW_wrapper &w, int n) {
             PetscErrorCode ierr;
             ierr = SOPFLOWSetNumScenarios(w.sopf, n);
             ExaGOCheckError(ierr);
           })
      .def("set_wind_gen_profile",
           [](SOPFLOW_wrapper &w, std::string filename) {
             PetscErrorCode ierr;
             ierr = SOPFLOWSetWindGenProfile(w.sopf, filename.c_str());
             ExaGOCheckError(ierr);
           })
      .def("set_time_step_and_duration",
           [](SOPFLOW_wrapper &w, double dt, double tmax) {
             PetscErrorCode ierr;
             ierr = SOPFLOWSetTimeStepandDuration(w.sopf, dt, tmax);
             ExaGOCheckError(ierr);
           })
      .def("set_tolerance",
           [](SOPFLOW_wrapper &w, double tol) {
             PetscErrorCode ierr;
             ierr = SOPFLOWSetTolerance(w.sopf, tol);
             ExaGOCheckError(ierr);
           })
      .def("set_subproblem_verbosity_level",
           [](SOPFLOW_wrapper &w, int level) {
             PetscErrorCode ierr;
             ierr = SOPFLOWSetSubproblemVerbosityLevel(w.sopf, level);
             ExaGOCheckError(ierr);
           })
      .def("set_subproblem_compute_mode",
           [](SOPFLOW_wrapper &w, std::string mode) {
             PetscErrorCode ierr;
             ierr = SOPFLOWSetSubproblemComputeMode(w.sopf, mode.c_str());
             ExaGOCheckError(ierr);
           })
      .def("set_subproblem_model",
           [](SOPFLOW_wrapper &w, std::string model) {
             PetscErrorCode ierr;
             ierr = SOPFLOWSetSubproblemModel(w.sopf, model.c_str());
             ExaGOCheckError(ierr);
           })
      .def("set_subproblem_solver",
           [](SOPFLOW_wrapper &w, std::string solver) {
             PetscErrorCode ierr;
             ierr = SOPFLOWSetSubproblemSolver(w.sopf, solver.c_str());
             ExaGOCheckError(ierr);
           })
      .def("set_solver",
           [](SOPFLOW_wrapper &w, std::string solver) {
             PetscErrorCode ierr;
             ierr = SOPFLOWSetSolver(w.sopf, solver.c_str());
             ExaGOCheckError(ierr);
           })

      .def("set_initialization_type",
           [](SOPFLOW_wrapper &w, std::string init) {
             PetscErrorCode ierr;
             OPFLOWInitializationType type;
             bool found = false;
#define OPFLOW_INIT_DEF(x)                                                     \
  if (init == #x) {                                                            \
    found = true;                                                              \
    type = x;                                                                  \
  }
#include "private/opflow_enum.def"
#undef OPFLOW_INIT_DEF
             if (!found) {
               throw ExaGOError(
                   fmt::format(
                       "Invalid value '{}' for OPFLOWInitializationType enum",
                       init)
                       .c_str());
             }
             ierr = SOPFLOWSetInitializationType(w.sopf, type);
             ExaGOCheckError(ierr);
           })
      .def("set_initialization_type",
           [](SOPFLOW_wrapper &w, int init) {
             PetscErrorCode ierr;
             OPFLOWInitializationType type;
             bool found = false;
#define OPFLOW_INIT_DEF(x)                                                     \
  if (init == x) {                                                             \
    found = true;                                                              \
    type = x;                                                                  \
  }
#include "private/opflow_enum.def"
#undef OPFLOW_INIT_DEF
             if (!found) {
               throw ExaGOError(
                   fmt::format(
                       "Invalid value '{}' for SOPFLOWInitializationType enum",
                       init)
                       .c_str());
             }
             ierr = SOPFLOWSetInitializationType(w.sopf, type);
             ExaGOCheckError(ierr);
           })
      .def("set_gen_bus_voltage_type",
           [](SOPFLOW_wrapper &w, int volt) {
             PetscErrorCode ierr;
             OPFLOWGenBusVoltageType type;
             bool found = false;
#define OPFLOW_GBV_DEF(x)                                                      \
  if (volt == x) {                                                             \
    found = true;                                                              \
    type = x;                                                                  \
  }
#include "private/opflow_enum.def"
#undef OPFLOW_GBV_DEF
             if (!found) {
               throw ExaGOError(
                   fmt::format(
                       "Invalid value '{}' for SOPFLOWGenBusVoltageType enum",
                       volt)
                       .c_str());
             }
             ierr = SOPFLOWSetGenBusVoltageType(w.sopf, type);
             ExaGOCheckError(ierr);
           })

      .def("set_gen_bus_voltage_type",
           [](SOPFLOW_wrapper &w, std::string volt) {
             PetscErrorCode ierr;
             OPFLOWGenBusVoltageType type;
             bool found = false;
#define OPFLOW_GBV_DEF(x)                                                      \
  if (volt == #x) {                                                            \
    found = true;                                                              \
    type = x;                                                                  \
  }
#include "private/opflow_enum.def"
#undef OPFLOW_GBV_DEF
             if (!found) {
               throw ExaGOError(
                   fmt::format(
                       "Invalid value '{}' for SOPFLOWGenBusVoltageType enum",
                       volt)
                       .c_str());
             }
             ierr = SOPFLOWSetGenBusVoltageType(w.sopf, type);
             ExaGOCheckError(ierr);
           })
      .def("enable_multi_contingency",
           [](SOPFLOW_wrapper &w, bool flag) {
             PetscErrorCode ierr;
             ierr = SOPFLOWEnableMultiContingency(w.sopf, (PetscBool)flag);
             ExaGOCheckError(ierr);
           })

      /* Getters */

      .def("get_num_scenarios",
           [](SOPFLOW_wrapper &w, const std::string &filename,
              ScenarioFileInputFormat fmt) -> int {
             PetscErrorCode ierr;
             PetscInt n;
             ierr = SOPFLOWGetNumScenarios(w.sopf, fmt, filename.c_str(), &n);
             ExaGOCheckError(ierr);
             return n;
           })
      .def("get_num_iterations",
           [](SOPFLOW_wrapper &w) -> int {
             PetscErrorCode ierr;
             PetscInt n;
             ierr = SOPFLOWGetNumIterations(w.sopf, &n);
             ExaGOCheckError(ierr);
             return n;
           })
      .def("get_convergence_status",
           [](SOPFLOW_wrapper &w) -> bool {
             PetscErrorCode ierr;
             PetscBool flag;
             ierr = SOPFLOWGetConvergenceStatus(w.sopf, &flag);
             ExaGOCheckError(ierr);
             return flag;
           })
      .def("get_total_objective",
           [](SOPFLOW_wrapper &w) -> double {
             PetscErrorCode ierr;
             double obj;
             ierr = SOPFLOWGetTotalObjective(w.sopf, &obj);
             ExaGOCheckError(ierr);
             return obj;
           })
      .def("get_base_objective",
           [](SOPFLOW_wrapper &w) -> double {
             PetscErrorCode ierr;
             double obj;
             ierr = SOPFLOWGetBaseObjective(w.sopf, &obj);
             ExaGOCheckError(ierr);
             return obj;
           })
      .def("get_converged_status",
           [](SOPFLOW_wrapper &w) -> bool {
             PetscErrorCode ierr;
             PetscBool flag;
             ierr = SOPFLOWGetConvergenceStatus(w.sopf, &flag);
             return (bool)flag;
           })
      .def("get_tolerance",
           [](SOPFLOW_wrapper &w) -> double {
             PetscErrorCode ierr;
             double tol;
             ierr = SOPFLOWGetTolerance(w.sopf, &tol);
             ExaGOCheckError(ierr);
             return tol;
           })

      /* Solution */

      .def("setup",
           [](SOPFLOW_wrapper &w) {
             PetscErrorCode ierr;
             ierr = SOPFLOWSetUp(w.sopf);
             ExaGOCheckError(ierr);
           })
      .def("solve",
           [](SOPFLOW_wrapper &w) {
             PetscErrorCode ierr;
             ierr = SOPFLOWSolve(w.sopf);
             ExaGOCheckError(ierr);
           })

      ;
}
