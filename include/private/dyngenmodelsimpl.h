/**
 * @file dyngenmodelsimpl.h
 * @brief Private header file that defines the base generator model object
 */
#ifndef DYNGENMODELSIMPL_H
#define DYNGENMODELSIMPL_H

#include <ps.h>
#include <dyngenmodels.h>
#include <dynexcmodels.h>
#include <dynturbgovmodels.h>

#define DYNGENMODELSMAX 10

#define DYNGENMONITORSMAX 10

extern PetscBool DYNGenModelsRegisterAllCalled;

extern PetscInt ngenmodelsregistered;

/**
 * @brief private struct to hold a list of registered dynamic generator models
 */
struct _p_DYNGenModelList{
  char     name[16]; /**< Name of the model */
  PetscInt key; /**< Key returned by DMNetworkRegisterComponent */
  PetscInt sizeofstruct; /**< size of struct obtained by sizeof() */
  PetscErrorCode (*create)(DYNGenModel); /**< Create routine */
};

extern struct _p_DYNGenModelList DYNGenModelList[DYNGENMODELSMAX];

/**
 * @brief operations available with dynamic generator model
 */
struct _p_DYNGenModelOps{
  PetscErrorCode (*getnvar)(DYNGenModel,PetscInt*); /**< get number of variables associated with this model */
  PetscErrorCode (*readdata)(DYNGenModel,char*,PetscScalar,PetscScalar); /** < read data for this model */
  PetscErrorCode (*getbusnumid)(DYNGenModel,PetscInt*,char**); /**< get the bus number and generator id */
  PetscErrorCode (*destroy)(DYNGenModel); /**< destroy the generator model */
  PetscErrorCode (*getsizeof)(DYNGenModel,PetscInt*); /**< return size of the implementation specific data struct (available through getsizeof()) */
  PetscErrorCode (*setinitialconditions)(DYNGenModel,PetscScalar,PetscScalar,PetscScalar,PetscScalar,PetscScalar*); /**< set initial conditions for dynamics simulation */
  PetscErrorCode (*setinitialconditionsp)(DYNGenModel,PetscScalar,PetscScalar,PetscScalar,PetscScalar,PetscInt,PetscInt,PetscInt,PetscInt,Mat,PetscInt); /**< set differentiation of the initial conditions for dynamcs simulations */
  PetscErrorCode (*setinitialfieldvoltagediff)(DYNGenModel,PetscScalar,PetscScalar,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*,PetscScalar*,PetscScalar*); /**< set differentiation of the initial conditions for dynamcs simulations */
  PetscErrorCode (*daerhsfunction)(DYNGenModel,PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*,PetscScalar*,PetscScalar*); /**< return the differential-algebraic equations for this model alogng with the generator currents */
  PetscErrorCode (*daerhsjacobian)(DYNGenModel,Mat,PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscInt,PetscInt[],PetscInt[]); /**< set the jacobian portion for this model */
  PetscErrorCode (*getequationtypes)(DYNGenModel,PetscInt*,PetscInt*,PetscInt*); /**< Return the number of differential and algebraic equations */
  PetscErrorCode (*getinitialfieldvoltage)(DYNGenModel,PetscScalar*); /**< get the initial field voltage */
  PetscErrorCode (*getinitialmechanicalpower)(DYNGenModel,PetscScalar*); /**< get initial mechanical power */
  PetscErrorCode (*getfieldcurrent)(DYNGenModel,PetscScalar*,PetscScalar*); /**< return field current */
  PetscErrorCode (*getfrequency)(DYNGenModel,PetscReal,const PetscScalar*,PetscScalar*); /**< return generator frequency */
  PetscErrorCode (*setfrequency)(DYNGenModel,PetscScalar*,PetscScalar); /**< set generator frequency */
  PetscErrorCode (*getspeeddeviation)(DYNGenModel,PetscReal,const PetscScalar*,PetscScalar*,PetscScalar*); /**< get speed deviation (pu) */
  PetscErrorCode (*getspeeddeviationlocation)(DYNGenModel,PetscInt*); /** < get location of speed deviation variable in the vector of variables for this model */
  PetscErrorCode (*getdfreqdstate)(DYNGenModel,PetscReal,const PetscScalar*,PetscScalar*,PetscInt*); /**< get partial derivatives of the frequency w.r.t state variables for this generator */
  PetscErrorCode (*daerhsjacobianp)(DYNGenModel,PetscReal,const PetscScalar*,Mat,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt); /**< Jacobian w.r.t parameters (used in adjoint sensitivity) */
  PetscErrorCode (*daefwdrhsjacobianp)(DYNGenModel,PetscReal,const PetscScalar*,Vec*,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt); /**< jacobian w.r.t. parameters (used in forward sensitivity */
  PetscErrorCode (*setevent)(DYNGenModel,PetscInt*,PetscInt*,PetscBool*,PetscErrorCode (**)(DYNGenModel,PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*),PetscErrorCode (**)(DYNGenModel,PetscInt,PetscInt[],PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscBool,PetscBool*),PetscErrorCode (**)(DYNGenModel,PetscInt,PetscInt[],PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*,PetscScalar*,PetscInt,PetscInt,PetscInt,PetscInt)); /**< set the events for this generator model */
  PetscErrorCode (*getcurrent)(DYNGenModel,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*,PetscScalar*);
  PetscErrorCode (*setvoltage)(DYNGenModel,PetscScalar,PetscScalar);
};

typedef struct _p_DYNGenModelOps DYNGenModelOps;

/**
 * @brief private struct for base dynamic generator model
 */
struct _p_DYNGenModel{
  char             type[16]; /**< model type */
  void*            data;     /**< model data */
  DYNGenModelOps   ops;      /**< common operations for each model */
  PetscInt         nvar;     /**< Number of variables for the model */
  PetscInt         startloc; /**< Starting location for the variables for this model */
  PetscInt         ndiff;    /**< Number of differential equations for this model */
  PetscInt         nalg;     /**< Number of algebraic equations for this model */
  PetscInt         *eqtypes; /**< Type of equation: DIFF_EQ for Differential, ALG_EQ for algebraic */
  DYNExcModel      dynexc;   /**< The exciter model associated with this generator model */
  DYNTurbgovModel  dynturbgov; /**< The turbine governor model associated with this generator model */
  PSBUS            bus;      /**< The bus on which this generator is incident */
  PSGEN            psgen;      /**< The PSGEN object that holds this dyngen model */
  PetscInt         numMonitors; /**< Number of event monitors */
  PetscInt         direction[DYNGENMONITORSMAX]; /**< zero-crossing direction for event location */
  PetscBool        terminate[DYNGENMONITORSMAX]; /** < event termination flag */
  PetscErrorCode   (*eventfcn)(DYNGenModel,PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*); /**< event handler function for this generator model */
  PetscErrorCode   (*posteventfcn)(DYNGenModel,PetscInt,PetscInt[],PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscBool,PetscBool*); /**< postevent handler function */
  PetscErrorCode   (*posteventdgamma)(DYNGenModel,PetscInt,PetscInt[],PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*,PetscScalar*,PetscInt,PetscInt,PetscInt,PetscInt); /**< Used in sensitivity analysis */

  /* Max and Min. frequency limits */
  PetscScalar   freq_max; /**< Max. allowable gen. frequency */
  PetscScalar   freq_min; /**< Min. allowable gen. frequency */

  PetscBool     iscv; /**< Flag to check if this generator is to be treated as a constant voltage generator */
};
/**
 * @brief Parses a line in dyr file and reads the model type (string identifier). If the line does not have the registered generator model then the type is set to NONE.
 * @param [in] char* line - the line to parse
 * @param [out] char[] modeltype - type of the model (string name of the model) 
 */
extern PetscErrorCode DYNGenModelReadModelType(char*,char[]);
/**
 * @brief Copies the dyngen model
 * @param [in] DYNGenModel dyngenin - the dyngen model to copy from
 * @param [out] DYNGenModel dyngenout - the dyngen model to copy to
 * Notes:
 * If a new method is created for the dyngenmodel then that should be copied too.
 * The underlying implementation (dyngen->data) data pointer is copied.
 */
extern PetscErrorCode DYNGenModelCopy(DYNGenModel,DYNGenModel);
/**
 * @brief Registers all built-in dynamic generator models
 */
extern PetscErrorCode DYNGenModelRegisterAll(void);
/**
 * @brief Set event monitors for this generator model
 * Notes:
 * Currently, only min. and max. frequency limits implemented
*/
extern PetscErrorCode DYNGenModelSetEventMonitor(DYNGenModel);

#endif
