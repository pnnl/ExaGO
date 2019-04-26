/**
 * @file dynturbgovmodelsimpl.h
 * @brief Private header files that define data and structures for base turbine-governor model
 */
#ifndef DYNTURBGOVMODELSIMPL_H
#define DYNTURBGOVMODELSIMPL_H

#include <ps.h>
#include <dynturbgovmodels.h>
#include <dyngenmodels.h>

#define DYNTURBGOVMODELSMAX 10

#define DYNTURBGOVMONITORSMAX 10

extern PetscInt nturbgovmodelsregistered;
extern PetscBool DYNTurbgovModelRegisterAllCalled;

/**
 * @brief private struct to hold a list of registered dynamic turbine-governor models
 */
struct _p_DYNTurbgovModelList{
  char     name[16]; /**< Name of the model */
  PetscInt key; /**< Key returned by DMNetworkRegisterComponent */
  PetscInt sizeofstruct; /**< size of struct obtained by sizeof() */
  PetscErrorCode (*create)(DYNTurbgovModel); /**< Create routine */
};

extern struct _p_DYNTurbgovModelList DYNTurbgovModelList[DYNTURBGOVMODELSMAX];

/**
 * @brief operations available with dynamic turbine-governor model
 */
struct _p_DYNTurbgovModelOps{
  PetscErrorCode (*getnvar)(DYNTurbgovModel,PetscInt*); /**< get number of variables (equations) */
  PetscErrorCode (*readdata)(DYNTurbgovModel,char*,PetscScalar,PetscScalar); /**< read data */
  PetscErrorCode (*getbusnumid)(DYNTurbgovModel,PetscInt*,char**); /**< get the bus number and generator id */
  PetscErrorCode (*destroy)(DYNTurbgovModel); /**< Destroy model */
  PetscErrorCode (*getsizeof)(DYNTurbgovModel,PetscInt*); /**< get size of implementation specific struct (through sizeof()) */
  PetscErrorCode (*getmechanicalpower)(DYNTurbgovModel,PetscReal,PetscScalar*,PetscScalar*,PetscScalar*,PetscInt*,PetscScalar[],PetscInt[]); /**< get the mechanical power and (optionally) the derivatives of mechanical power w.r.t. generator state variables */
  PetscErrorCode (*setinitialconditions)(DYNTurbgovModel,PetscScalar,PetscScalar,PetscScalar*); /**< set initial conditions */
  PetscErrorCode (*setinitialconditionsp)(DYNTurbgovModel,PetscScalar,PetscScalar,PetscScalar,PetscScalar,PetscInt,PetscInt,PetscInt,PetscInt,Mat,PetscInt); /**< differentiate initial conditions to parameters */
  PetscErrorCode (*daerhsfunction)(DYNTurbgovModel,PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*); /**< differential-algebraic equations for the turbine-governor model */
  PetscErrorCode (*daerhsjacobian)(DYNTurbgovModel,Mat,PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscInt,PetscInt[]); /**< Jacobian porition for turbine-governor model */
  PetscErrorCode (*daerhsjacobianp)(DYNTurbgovModel,PetscReal,const PetscScalar*,Mat,PetscInt,PetscInt); /**< Jacobian w.r.t. parameters (used in adjoint sensitivity) */
  PetscErrorCode (*daefwdrhsjacobianp)(DYNTurbgovModel,PetscReal,const PetscScalar*,Vec*,PetscInt,PetscInt,PetscInt); /**< Jacobian w.r.t parameters (used in forward sensitivity) */
  PetscErrorCode (*getequationtypes)(DYNTurbgovModel,PetscInt*,PetscInt*,PetscInt*); /**< return the type of equations (differential or algebraic) for this model */
  PetscErrorCode (*setevent)(DYNTurbgovModel,PetscInt*,PetscInt*,PetscBool*,PetscErrorCode (**)(DYNTurbgovModel,PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*),PetscErrorCode (**)(DYNTurbgovModel,PetscInt,PetscInt[],PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscBool,PetscBool*),PetscErrorCode (**)(DYNTurbgovModel,PetscInt,PetscInt[],PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*,PetscScalar*,PetscInt,PetscInt,PetscInt,PetscInt)); /**< set events for this model */
};

typedef struct _p_DYNTurbgovModelOps DYNTurbgovModelOps;

/**
 * @brief private struct for base dynamic turbine-governor model
 */
struct _p_DYNTurbgovModel{
  char             type[16]; /**< model type */
  void*            data; /**< model data */
  DYNTurbgovModelOps   ops;      /* common operations for each model */
  PetscInt         nvar; /**< Number of variables for the model */
  PetscInt         startloc; /**< Starting location for the variables for this model */
  PetscInt         ndiff; /**< Number of differential equations for this model */
  PetscInt         nalg; /**< Number of algebraic equations for this model */
  PetscInt         *eqtypes; /**< Type of equation: DIFF_EQ for Differential, ALG_EQ for algebraic */
  DYNGenModel      dyngen; /**< The generator model associated with this turbine governor model */
  PetscInt         numMonitors; /**< Number of event monitors */
  PetscInt         direction[DYNTURBGOVMONITORSMAX]; /**< zero-crossing direction for event location */
  PetscBool        terminate[DYNTURBGOVMONITORSMAX]; /**< event termination flag */
  PetscErrorCode   (*eventfcn)(DYNTurbgovModel,PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*); /** event handler function */
  PetscErrorCode   (*posteventfcn)(DYNTurbgovModel,PetscInt,PetscInt[],PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscBool,PetscBool*); /** post event handler function */
  PetscErrorCode   (*posteventdgamma)(DYNTurbgovModel,PetscInt,PetscInt[],PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*,PetscScalar*,PetscInt,PetscInt,PetscInt,PetscInt); /**< used in sensitivity analysis */
};
/**
 * @brief Parses a line in dyr file and reads the model type (string identifier). If the line does not have the registered turbgov model then the type is set to NONE.
 * @param [in] char* line - the line to parse
 * @param [out] char[] modeltype - type of the model (string name of the model)
 */
extern PetscErrorCode DYNTurbgovModelReadModelType(char*,char[]);
/**
 * @brief Copies the dynturbgov model
 * @param [in] DYNTurbgovModel dynturbgovin - the dynturbgov model to copy from
 * @param [out] DYNTurbgovModel dynturbgovout - the dynturbgov model to copy to
 * Notes:
 * If a new method is created for the dynturbgovmodel then that should be copied too.
 * The underlying implementation (dynturbgov->data) data pointer is copied.
 */
extern PetscErrorCode DYNTurbgovModelCopy(DYNTurbgovModel,DYNTurbgovModel);
/**
 * @brief Registers all built-in turbgov models
 */
extern PetscErrorCode DYNTurbgovModelRegisterAll(void);
/**
 * @brief Registers all built-in turbgov models
 * @param [in] DYNTurbgovModel dynturbgov - the dynamic turbgov model
 */
extern PetscErrorCode DYNTurbgovModelSetEventMonitor(DYNTurbgovModel);

#endif
