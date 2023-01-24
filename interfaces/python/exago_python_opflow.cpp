
#include "exago_python_opflow.hpp"

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

// -------------------------------------------------------------
// init_exago_opflow
// -------------------------------------------------------------
void init_exago_opflow(pybind11::module &m) {
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
#ifdef EXAGO_ENABLE_LOGGING
               throw ExaGOError(
                   fmt::format(
                       "Invalid value '{}' for OPFLOWObjectiveType enum", func)
                       .c_str());
#else
             char sbuf[256];
             sprintf(sbuf,"Invalid value %s for OPFLOWObjectiveType enum",func.c_str());
             throw ExaGOError(sbuf);
#endif
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
#ifdef EXAGO_ENABLE_LOGGING
               throw ExaGOError(
                   fmt::format(
                       "Invalid value '{}' for OPFLOWObjectiveType enum", func)
                       .c_str());
#else
             char sbuf[256];
             sprintf(sbuf,"Invalid value %d for OPFLOWObjectiveType enum",func);
             throw ExaGOError(sbuf);
#endif
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
#ifdef EXAGO_ENABLE_LOGGING
               throw ExaGOError(
                   fmt::format(
                       "Invalid value '{}' for OPFLOWInitializationType enum",
                       init)
                       .c_str());
#else
             char sbuf[256];
             sprintf(sbuf,"Invalid value %s for OPFLOWInitializationType enum",
                     init.c_str());
             throw ExaGOError(sbuf);
#endif
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
#ifdef EXAGO_ENABLE_LOGGING
               throw ExaGOError(
                   fmt::format(
                       "Invalid value '{}' for OPFLOWInitializationType enum",
                       init)
                       .c_str());
#else
             char sbuf[256];
             sprintf(sbuf,"Invalid value %d for OPFLOWInitializationType enum",
                     init);
             throw ExaGOError(sbuf);
#endif
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
#ifdef EXAGO_ENABLE_LOGGING
               throw ExaGOError(
                   fmt::format(
                       "Invalid value '{}' for OPFLOWGenBusVoltageType enum",
                       volt)
                       .c_str());
#else
             char sbuf[256];
             sprintf(sbuf,"Invalid value %s for OPFLOWGenBusVoltageType enum",
                     volt.c_str());
             throw ExaGOError(sbuf);
#endif
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
#ifdef EXAGO_ENABLE_LOGGING
               throw ExaGOError(
                   fmt::format(
                       "Invalid value '{}' for OPFLOWGenBusVoltageType enum",
                       volt)
                       .c_str());
#else
             char sbuf[256];
             sprintf(sbuf,"Invalid value %d for OPFLOWGenBusVoltageType enum",
                     volt);
             throw ExaGOError(sbuf);
#endif
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
             ierr = OPFLOWSetHIOPComputeMode(w.opf, mode);
             ExaGOCheckError(ierr);
           })

      .def("set_hiop_mem_space",
           [](OPFLOW_wrapper &w, std::string mem_space) {
             PetscErrorCode ierr;
             ierr = OPFLOWSetHIOPMemSpace(w.opf, mem_space);
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

      .def("set_weight",
           [](OPFLOW_wrapper &w, double wt) {
             PetscErrorCode ierr;
             ierr = OPFLOWSetWeight(w.opf, wt);
             ExaGOCheckError(ierr);
           })

      .def("ps_set_gen_power_limits",
           [](OPFLOW_wrapper &w, int gbus, std::string gid, double pt,
              double pb, double qt, double qb) {
             PetscErrorCode ierr;
             _p_PS *ps;
             ierr = OPFLOWGetPS(w.opf, &ps);
             PSSetGenPowerLimits(ps, gbus, gid.c_str(), pt, pb, qt, qb);
             ExaGOCheckError(ierr);
           })

      .def("set_lines_monitored",
           [](OPFLOW_wrapper &w, const std::vector<double> &kvlevels) {
             PetscErrorCode ierr;
             std::vector<PetscScalar> tmp(kvlevels.size());
             std::copy(kvlevels.begin(), kvlevels.end(),
                       std::back_inserter(tmp));
             ierr = OPFLOWSetLinesMonitored(
                 w.opf, static_cast<PetscInt>(kvlevels.size()), &tmp[0], NULL);
             ExaGOCheckError(ierr);
           })

      //  Note that, at the time of writing, the C/C++ implementation
      //  of this exists, but always fails. Leaving this here for
      //  completeness.
      .def("set_lines_monitored",
           [](OPFLOW_wrapper &w, const int &nkvlevels,
              const std::string &monitorfile) {
             PetscErrorCode ierr;
             ierr = OPFLOWSetLinesMonitored(w.opf, (PetscInt)nkvlevels, NULL,
                                            monitorfile.c_str());
             ExaGOCheckError(ierr);
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
             ierr = OPFLOWGetHIOPComputeMode(w.opf, &mode);
             ExaGOCheckError(ierr);
             return mode.c_str();
           })

      .def("get_hiop_mem_space",
           [](OPFLOW_wrapper &w) -> std::string {
             PetscErrorCode ierr;
             std::string mem_space(PETSC_MAX_PATH_LEN, '\0');
             ierr = OPFLOWGetHIOPMemSpace(w.opf, &mem_space);
             ExaGOCheckError(ierr);
             return mem_space.c_str();
           })

      .def("get_model",
           [](OPFLOW_wrapper &w) -> std::string {
             PetscErrorCode ierr;
             std::string model(PETSC_MAX_PATH_LEN, '\0');
             ierr = OPFLOWGetModel(w.opf, &model);
             ExaGOCheckError(ierr);
             return model.c_str();
           })

      .def("get_solver",
           [](OPFLOW_wrapper &w) -> std::string {
             PetscErrorCode ierr;
             std::string solver(PETSC_MAX_PATH_LEN, '\0');
             ierr = OPFLOWGetSolver(w.opf, &solver);
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

      .def("get_num_iterations",
           [](OPFLOW_wrapper &w) -> int {
             PetscErrorCode ierr;
             int n;
             ierr = OPFLOWGetNumIterations(w.opf, &n);
             ExaGOCheckError(ierr);
             return n;
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

      .def("skip_options",
           [](OPFLOW_wrapper &w, bool skip) {
             PetscErrorCode ierr;
             ierr = OPFLOWUseAGC(w.opf, (PetscBool)skip);
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
}
