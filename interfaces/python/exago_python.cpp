#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <mpi.h>
#include <pflow.h>
#include <private/pflowimpl.h>
#include <opflow.h>
#include <private/opflowimpl.h>
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

  pybind11::class_<_p_PFLOW>(m, "PFLOW")
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

  /* OPFLOW class */
  pybind11::class_<_p_OPFLOW>(m, "OPFLOW")
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
      .def("set_model",
           [](_p_OPFLOW &opf, std::string model) {
             PetscErrorCode ierr;
             ierr = OPFLOWSetModel(&opf, model.c_str());
             ExaGOCheckError(ierr);
           })

      .def("set_solver",
           [](_p_OPFLOW &opf, std::string solver) {
             PetscErrorCode ierr;
             ierr = OPFLOWSetSolver(&opf, solver.c_str());
             ExaGOCheckError(ierr);
           })

      .def("set_objective_type",
           [](_p_OPFLOW &opf, std::string func) {
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
             ierr = OPFLOWSetObjectiveType(&opf, type);
             ExaGOCheckError(ierr);
           })

      .def("set_objective_type",
           [](_p_OPFLOW &opf, int func) {
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
             ierr = OPFLOWSetObjectiveType(&opf, type);
             ExaGOCheckError(ierr);
           })

      .def("set_initialization_type",
           [](_p_OPFLOW &opf, std::string init) {
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
             ierr = OPFLOWSetInitializationType(&opf, type);
             ExaGOCheckError(ierr);
           })

      .def("set_initialization_type",
           [](_p_OPFLOW &opf, int init) {
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
             ierr = OPFLOWSetInitializationType(&opf, type);
             ExaGOCheckError(ierr);
           })

      .def("set_gen_bus_voltage_type",
           [](_p_OPFLOW &opf, std::string volt) {
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
             ierr = OPFLOWSetGenBusVoltageType(&opf, type);
             ExaGOCheckError(ierr);
           })

      .def("set_gen_bus_voltage_type",
           [](_p_OPFLOW &opf, int volt) {
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
             ierr = OPFLOWSetGenBusVoltageType(&opf, type);
             ExaGOCheckError(ierr);
           })

      .def("set_has_gen_set_point",
           [](_p_OPFLOW &opf, bool setpt) {
             PetscErrorCode ierr;
             ierr = OPFLOWHasGenSetPoint(&opf, (PetscBool)setpt);
             ExaGOCheckError(ierr);
           })

      .def("set_hiop_compute_mode",
           [](_p_OPFLOW &opf, std::string mode) {
             PetscErrorCode ierr;
             ierr = OPFLOWSetHIOPComputeMode(&opf, mode.c_str());
             ExaGOCheckError(ierr);
           })

      .def("set_has_loadloss",
           [](_p_OPFLOW &opf, bool loadloss) {
             PetscErrorCode ierr;
             ierr = OPFLOWHasLoadLoss(&opf, (PetscBool)loadloss);
             ExaGOCheckError(ierr);
           })

      .def("set_ignore_lineflow_constraints",
           [](_p_OPFLOW &opf, bool ignore) {
             PetscErrorCode ierr;
             ierr = OPFLOWIgnoreLineflowConstraints(&opf, (PetscBool)ignore);
             ExaGOCheckError(ierr);
           })

      .def("set_has_bus_power_imbalance",
           [](_p_OPFLOW &opf, bool powerimbalance) {
             PetscErrorCode ierr;
             ierr = OPFLOWHasBusPowerImbalance(&opf, (PetscBool)powerimbalance);
             ExaGOCheckError(ierr);
           })

      .def("set_use_agc",
           [](_p_OPFLOW &opf, bool agc) {
             PetscErrorCode ierr;
             ierr = OPFLOWUseAGC(&opf, (PetscBool)agc);
             ExaGOCheckError(ierr);
           })

      .def("set_hiop_verbosity_level",
           [](_p_OPFLOW &opf, int level) {
             PetscErrorCode ierr;
             ierr = OPFLOWSetHIOPVerbosityLevel(&opf, level);
             ExaGOCheckError(ierr);
           })

      .def("set_loadloss_penalty",
           [](_p_OPFLOW &opf, double p) {
             PetscErrorCode ierr;
             ierr = OPFLOWSetLoadLossPenalty(&opf, p);
             ExaGOCheckError(ierr);
           })

      .def("set_bus_power_imbalance_penalty",
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

      .def("ps_set_gen_power_limits",
           [](_p_OPFLOW &opf, int gbus, std::string gid, double pt, double pb,
              double qt, double qb) {
             PetscErrorCode ierr;
             _p_PS *ps;
             ierr = OPFLOWGetPS(&opf, &ps);
             ExaGOCheckError(ierr);
             PSSetGenPowerLimits(ps, gbus, gid.c_str(), pt, pb, qt, qb);
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

      .def("get_hiop_compute_mode",
           [](_p_OPFLOW &opf) -> std::string {
             PetscErrorCode ierr;
             std::string mode(PETSC_MAX_PATH_LEN, '\0');
             ierr = OPFLOWGetHIOPComputeMode(&opf, &mode[0]);
             ExaGOCheckError(ierr);
             return mode.c_str();
           })

      .def("get_model",
           [](_p_OPFLOW &opf) -> std::string {
             PetscErrorCode ierr;
             std::string model(PETSC_MAX_PATH_LEN, '\0');
             ierr = OPFLOWGetModel(&opf, &model[0]);
             ExaGOCheckError(ierr);
             return model.c_str();
           })

      .def("get_solver",
           [](_p_OPFLOW &opf) -> std::string {
             PetscErrorCode ierr;
             std::string solver(PETSC_MAX_PATH_LEN, '\0');
             ierr = OPFLOWGetSolver(&opf, &solver[0]);
             ExaGOCheckError(ierr);
             return solver.c_str();
           })

      .def("get_convergence_status",
           [](_p_OPFLOW opf) -> bool {
             PetscErrorCode ierr;
             PetscBool status;
             ierr = OPFLOWGetConvergenceStatus(&opf, &status);
             ExaGOCheckError(ierr);
             return status;
           })

      .def("get_objective_type",
           [](_p_OPFLOW opf) -> OPFLOWObjectiveType {
             PetscErrorCode ierr;
             OPFLOWObjectiveType obj_type;
             ierr = OPFLOWGetObjectiveType(&opf, &obj_type);
             ExaGOCheckError(ierr);
             return obj_type;
           })

      .def("get_initialization_type",
           [](_p_OPFLOW opf) -> OPFLOWInitializationType {
             PetscErrorCode ierr;
             OPFLOWInitializationType init;
             ierr = OPFLOWGetInitializationType(&opf, &init);
             ExaGOCheckError(ierr);
             return init;
           })

      .def("get_gen_bus_voltage_type",
           [](_p_OPFLOW &opf) -> OPFLOWGenBusVoltageType {
             PetscErrorCode ierr;
             OPFLOWGenBusVoltageType volt;
             ierr = OPFLOWGetGenBusVoltageType(&opf, &volt);
             ExaGOCheckError(ierr);
             return volt;
           })

      .def("get_has_gen_set_point",
           [](_p_OPFLOW &opf) -> bool {
             PetscErrorCode ierr;
             PetscBool setpt;
             ierr = OPFLOWGetHasGenSetPoint(&opf, &setpt);
             ExaGOCheckError(ierr);
             return setpt;
           })

      .def("get_loadloss_penalty",
           [](_p_OPFLOW &opf) -> double {
             PetscErrorCode ierr;
             double p;
             ierr = OPFLOWGetLoadLossPenalty(&opf, &p);
             ExaGOCheckError(ierr);
             return p;
           })

      .def("get_ignore_lineflow_constraints",
           [](_p_OPFLOW &opf) -> bool {
             PetscErrorCode ierr;
             PetscBool ig;
             ierr = OPFLOWGetIgnoreLineflowConstraints(&opf, &ig);
             ExaGOCheckError(ierr);
             return ig;
           })

      .def("get_has_loadloss",
           [](_p_OPFLOW &opf) -> bool {
             PetscErrorCode ierr;
             PetscBool loadloss;
             ierr = OPFLOWGetHasLoadLoss(&opf, &loadloss);
             ExaGOCheckError(ierr);
             return loadloss;
           })

      .def("get_has_bus_power_imbalance",
           [](_p_OPFLOW &opf) -> bool {
             PetscErrorCode ierr;
             PetscBool powerimbalance;
             ierr = OPFLOWGetHasBusPowerImbalance(&opf, &powerimbalance);
             ExaGOCheckError(ierr);
             return powerimbalance;
           })

      .def("get_use_agc",
           [](_p_OPFLOW &opf) -> bool {
             PetscErrorCode ierr;
             PetscBool agc;
             ierr = OPFLOWGetUseAGC(&opf, &agc);
             ExaGOCheckError(ierr);
             return agc;
           })

      .def("get_hiop_verbosity_level",
           [](_p_OPFLOW &opf) -> int {
             PetscErrorCode ierr;
             int level;
             ierr = OPFLOWGetHIOPVerbosityLevel(&opf, &level);
             ExaGOCheckError(ierr);
             return level;
           })

      .def("get_bus_power_imbalance_penalty",
           [](_p_OPFLOW &opf) -> double {
             PetscErrorCode ierr;
             double p;
             ierr = OPFLOWGetBusPowerImbalancePenalty(&opf, &p);
             ExaGOCheckError(ierr);
           })

      .def("get_gen_dispatch",
           [](_p_OPFLOW opf, int gbus, std::string gid) -> pybind11::tuple {
             _p_PS *ps;
             PetscErrorCode ierr;
             ierr = OPFLOWGetPS(&opf, &ps);
             ExaGOCheckError(ierr);
             double pg;
             double qg;
             ierr = PSGetGenDispatch(ps, gbus, gid.c_str(), &pg, &qg);
             ExaGOCheckError(ierr);
             return pybind11::make_tuple(pg, qg);
           })

      .def("get_objective_types",
           [](_p_OPFLOW &opf) -> std::vector<OPFLOWObjectiveType> {
             std::vector<OPFLOWObjectiveType> types;
#define OPFLOW_OBJ_DEF(x) types.push_back(x);
#include "private/opflow_enum.def"
#undef OPFLOW_OBJ_DEF
             return types;
           })

      .def("get_initialization_types",
           [](_p_OPFLOW &opf) -> std::vector<OPFLOWInitializationType> {
             std::vector<OPFLOWInitializationType> types;
#define OPFLOW_INIT_DEF(x) types.push_back(x);
#include "private/opflow_enum.def"
#undef OPFLOW_INIT_DEF
             return types;
           })

      .def("get_gen_bus_voltage_types",
           [](_p_OPFLOW &opf) -> std::vector<OPFLOWGenBusVoltageType> {
             std::vector<OPFLOWGenBusVoltageType> types;
#define OPFLOW_GBV_DEF(x) types.push_back(x);
#include "private/opflow_enum.def"
#undef OPFLOW_GBV_DEF
             return types;
           })

      .def("get_objective",
           [](_p_OPFLOW &opf) -> double {
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

      .def("solution_to_ps",
           [](_p_OPFLOW &opf) {
             PetscErrorCode ierr;
             ierr = OPFLOWSolutionToPS(&opf);
             ExaGOCheckError(ierr);
           })

      .def("set_up_ps",
           [](_p_OPFLOW &opf) {
             PetscErrorCode ierr;
             ierr = OPFLOWSetUpPS(&opf);
             ExaGOCheckError(ierr);
           })

      .def("save_solution",
           [](_p_OPFLOW &opf, OutputFormat fmt, std::string outfile) {
             PetscErrorCode ierr;
             ierr = OPFLOWSaveSolution(&opf, fmt, outfile.c_str());
             ExaGOCheckError(ierr);
           })

      .def("read_mat_power_data", [](_p_OPFLOW &opf, std::string filename) {
        PetscErrorCode ierr;
        ierr = OPFLOWReadMatPowerData(&opf, filename.c_str());
        ExaGOCheckError(ierr);
      });
}
