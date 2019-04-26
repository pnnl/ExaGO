/**
 * @file dynstabmodelsimpl.h
 * @brief Private header file that defines data and structures for the base stabilizer model
 */

#ifndef DYNSTABMODELSIMPL_H
#define DYNSTABMODELSIMPL_H

#include <ps.h>
#include <dynstabmodels.h>
#include <dyngenmodels.h>
#include <dynexcmodels.h>

#define DYNSTABMODELSMAX 10

#define DYNSTABMONITORSMAX 10

extern PetscInt nstabmodelsregistered; /**< Default = 0 */
extern PetscBool DYNStabModelRegisterAllCalled; /**< Default = PETSC_FALSE */

/**
 * @brief private struct to hold a list of registered dynamic stabilizer models
 */
struct _p_DYNStabModelList{
  char     name[16]; /**<  Name of the model */
  PetscInt key; /**<  Key returned by DMNetworkRegisterComponent */
  PetscInt sizeofstruct; /**<  size of struct obtained by sizeof() */
  PetscErrorCode (*create)(DYNStabModel); /**<  Create routine */
};

extern struct _p_DYNStabModelList DYNStabModelList[DYNSTABMODELSMAX];

/**
 * @brief operations available with dynamic stabilizer model
 */
struct _p_DYNStabModelOps{
  PetscErrorCode (*getnvar)(DYNStabModel,PetscInt*); /**< get number of variables */
  PetscErrorCode (*readdata)(DYNStabModel,char*); /**< read data */
  PetscErrorCode (*getbusnumid)(DYNStabModel,PetscInt*,char**); /**< get the bus number and generator id */
  PetscErrorCode (*destroy)(DYNStabModel); /**< Destroy model */
  PetscErrorCode (*getsizeof)(DYNStabModel,PetscInt*); /**< get the size of implementation specific struct (through sizeof()) */
  PetscErrorCode (*setinitialconditions)(DYNStabModel,PetscScalar,PetscScalar,PetscScalar*); /**< set initial conditions */
  PetscErrorCode (*daerhsfunction)(DYNStabModel,PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*); /**< differential-algebraic equations for this model */
  PetscErrorCode (*daerhsjacobian)(DYNStabModel,Mat,PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscInt,PetscInt[]); /**< Jacobian portion for this model */
  PetscErrorCode (*daerhsjacobianp)(DYNStabModel,PetscReal,const PetscScalar*,Mat,PetscInt,PetscInt); /**< Jacobian w.r.t. parameters (used in adjoint sensitivity */
  PetscErrorCode (*daefwdrhsjacobianp)(DYNStabModel,PetscReal,const PetscScalar*,Vec*,PetscInt,PetscInt,PetscInt); /**< Jacobian w.r.t. parameters (used in forward sensitivity) */
  PetscErrorCode (*getvothsg)(DYNStabModel,PetscReal,PetscScalar*,PetscScalar*,PetscInt*,PetscScalar[],PetscInt[]); /**< Get stabilizer signal and its partial derivatives (optional) */
  PetscErrorCode (*getequationtypes)(DYNStabModel,PetscInt*,PetscInt*,PetscInt*); /**< Get the equation types (differential or algebraic) */
  PetscErrorCode (*setevent)(DYNStabModel,PetscInt*,PetscInt*,PetscBool*,PetscErrorCode (**)(DYNStabModel,PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*),PetscErrorCode (**)(DYNStabModel,PetscInt,PetscInt[],PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscBool,PetscBool*),PetscErrorCode (**)(DYNStabModel,PetscInt,PetscInt[],PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*,PetscScalar*,PetscInt,PetscInt,PetscInt,PetscInt)); /**< set events for this model */
};

typedef struct _p_DYNStabModelOps DYNStabModelOps;

/**
 * @brief private struct for base dynamic stabilizer model
 */
struct _p_DYNStabModel{
  char             type[16]; /**<  model type */
  void*            data; /**<  model data */
  DYNStabModelOps   ops; /**<  common operations for each model */
  PetscInt         nvar; /**<  Number of variables for the model */
  PetscInt         startloc; /**<  Starting location for the variables for this model */
  PetscInt         ndiff; /**<  Number of differential equations for this model */
  PetscInt         nalg; /**<  Number of algebraic equations for this model */
  PetscInt         *eqtypes; /**<  Type of equation: DIFF_EQ for Differential, ALG_EQ for algebraic */
  DYNGenModel      dyngen; /**<  The generator model associated with this stabilizer governor model */
  DYNExcModel      dynexc; /**<  The exciter model associated with this stabilizer model */
  PetscInt         numMonitors; /**<  Number of event monitors */
  PetscInt         direction[DYNSTABMONITORSMAX]; /**< zero-crossing direction for event location */
  PetscBool        terminate[DYNSTABMONITORSMAX]; /**< event termination flag */
  PetscErrorCode   (*eventfcn)(DYNStabModel,PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*); /**< event handler function */
  PetscErrorCode   (*posteventfcn)(DYNStabModel,PetscInt,PetscInt[],PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscBool,PetscBool*); /**< post event handler function */
  PetscErrorCode   (*posteventdgamma)(DYNStabModel,PetscInt,PetscInt[],PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*,PetscScalar*,PetscInt,PetscInt,PetscInt,PetscInt); /**< used in sensitivity analysis */
};
/**
 * @brief Parses a line in dyr file and reads the model type (string identifier). If the line does not have the registered stab model then the type is set to NONE.
 * @param [in] char* line - the line to parse
 * @param [out] char[] modeltype - type of the model (string name of the model)
 */
extern PetscErrorCode DYNStabModelReadModelType(char*,char[]);
/**
 * @brief Copies the dynstab model
 * @param [in] DYNStabModel dynstabin - the dynstab model to copy from
 * @param [out] DYNStabModel dynstabout - the dynstab model to copy to
 * Notes:
 * If a new method is created for the dynstabmodel then that should be copied too.
 * The underlying implementation (dynstab->data) data pointer is copied.
 */
extern PetscErrorCode DYNStabModelCopy(DYNStabModel,DYNStabModel);
/**
 * @brief Registers all built-in stab models
 */
extern PetscErrorCode DYNStabModelRegisterAll(void);
/**
 * @brief Sets the event info for the stab model
 * @param [in] DYNStabModel dynstab - the dynamic stab model
 */
extern PetscErrorCode DYNStabModelSetEventMonitor(DYNStabModel);

#endif
