/**
 * @file dynloadmodelsimpl.h
 * @brief Private header file that defines the base load model object
 */
#ifndef DYNLOADMODELSIMPL_H
#define DYNLOADMODELSIMPL_H

#include <ps.h>
#include <dynloadmodels.h>

#define DYNLOADMODELSMAX 10

#define DYNLOADMONITORSMAX 10

extern PetscBool DYNLoadModelsRegisterAllCalled;

extern PetscInt nloadmodelsregistered;

/**
 * @brief private struct to hold a list of registered dynamic load models
 */
struct _p_DYNLoadModelList{
  char     name[16]; /**< Name of the model */
  PetscInt key; /**< Key returned by DMNetworkRegisterComponent */
  PetscInt sizeofstruct; /**< size of struct obtained by sizeof() */
  PetscErrorCode (*create)(DYNLoadModel); /**< Create routine */
};

extern struct _p_DYNLoadModelList DYNLoadModelList[DYNLOADMODELSMAX];

/**
 * @brief operations available with dynamic load model
 */
struct _p_DYNLoadModelOps{
  PetscErrorCode (*getnvar)(DYNLoadModel,PetscInt*); /**< get number of variables associated with this model */
  PetscErrorCode (*readdata)(DYNLoadModel,char*,PetscScalar,PetscScalar); /** < read data for this model */
  PetscErrorCode (*getbusnumid)(DYNLoadModel,PetscInt*,char**); /**< get the bus number and load id */
  PetscErrorCode (*destroy)(DYNLoadModel); /**< destroy the load model */
  PetscErrorCode (*getsizeof)(DYNLoadModel,PetscInt*); /**< return size of the implementation specific data struct (available through getsizeof()) */
  PetscErrorCode (*setinitialconditions)(DYNLoadModel,PetscScalar,PetscScalar,PetscScalar,PetscScalar,PetscScalar*); /**< set initial conditions for dynamics simulation */
  PetscErrorCode (*daerhsfunction)(DYNLoadModel,PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*,PetscScalar*,PetscScalar*); /**< return the differential-algebraic equations for this model alogng with the load currents */
  PetscErrorCode (*daerhsjacobian)(DYNLoadModel,Mat,PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscInt,PetscInt[],PetscInt[]); /**< set the jacobian portion for this model */
  PetscErrorCode (*getequationtypes)(DYNLoadModel,PetscInt*,PetscInt*,PetscInt*); /**< Return the number of differential and algebraic equations */
  PetscErrorCode (*getinitialmechanicalpower)(DYNLoadModel,PetscScalar*); /**< get initial mechanical power */
  PetscErrorCode (*getfrequency)(DYNLoadModel,PetscReal,const PetscScalar*,PetscScalar*); /**< return load frequency */
  PetscErrorCode (*setfrequency)(DYNLoadModel,PetscScalar*,PetscScalar); /**< set load frequency */
  PetscErrorCode (*getspeeddeviation)(DYNLoadModel,PetscReal,const PetscScalar*,PetscScalar*,PetscScalar*); /**< get speed deviation (pu) */
  PetscErrorCode (*getspeeddeviationlocation)(DYNLoadModel,PetscInt*); /** < get location of speed deviation variable in the vector of variables for this model */
  PetscErrorCode (*getdfreqdstate)(DYNLoadModel,PetscReal,const PetscScalar*,PetscScalar*,PetscInt*); /**< get partial derivatives of the frequency w.r.t state variables for this load */
  PetscErrorCode (*daerhsjacobianp)(DYNLoadModel,PetscReal,const PetscScalar*,Mat,PetscInt,PetscInt,PetscInt); /**< Jacobian w.r.t parameters (used in adjoint sensitivity) */
  PetscErrorCode (*daefwdrhsjacobianp)(DYNLoadModel,PetscReal,const PetscScalar*,Vec*,PetscInt,PetscInt,PetscInt); /**< jacobian w.r.t. parameters (used in forward sensitivity */
  PetscErrorCode (*setevent)(DYNLoadModel,PetscInt*,PetscInt*,PetscBool*,PetscErrorCode (**)(DYNLoadModel,PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*),PetscErrorCode (**)(DYNLoadModel,PetscInt,PetscInt[],PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscBool,PetscBool*),PetscErrorCode (**)(DYNLoadModel,PetscInt,PetscInt[],PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*,PetscScalar*,PetscInt,PetscInt,PetscInt,PetscInt)); /**< set the events for this load model */
  PetscErrorCode (*setconstantpowerload)(DYNLoadModel,PetscScalar,PetscScalar);
};

typedef struct _p_DYNLoadModelOps DYNLoadModelOps;

/**
 * @brief private struct for base dynamic load model
 */
struct _p_DYNLoadModel{
  char             type[16]; /**< model type */
  void*            data;     /**< model data */
  DYNLoadModelOps   ops;      /**< common operations for each model */
  PetscInt         nvar;     /**< Number of variables for the model */
  PetscInt         startloc; /**< Starting location for the variables for this model */
  PetscInt         ndiff;    /**< Number of differential equations for this model */
  PetscInt         nalg;     /**< Number of algebraic equations for this model */
  PetscInt         *eqtypes; /**< Type of equation: DIFF_EQ for Differential, ALG_EQ for algebraic */
  PSBUS            bus;      /**< The bus on which this load is incident */
  PSLOAD            psload;      /**< The PSLOAD object that holds this dynload model */
  PetscInt         numMonitors; /**< Number of event monitors */
  PetscInt         direction[DYNLOADMONITORSMAX]; /**< zero-crossing direction for event location */
  PetscBool        terminate[DYNLOADMONITORSMAX]; /** < event termination flag */
  PetscErrorCode   (*eventfcn)(DYNLoadModel,PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*); /**< event handler function for this load model */
  PetscErrorCode   (*posteventfcn)(DYNLoadModel,PetscInt,PetscInt[],PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscBool,PetscBool*); /**< postevent handler function */
  PetscErrorCode   (*posteventdgamma)(DYNLoadModel,PetscInt,PetscInt[],PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*,PetscScalar*,PetscInt,PetscInt,PetscInt,PetscInt); /**< Used in sensitivity analysis */

  /* Max and Min. frequency limits */
  PetscScalar   freq_max; /**< Max. allowable load. frequency */
  PetscScalar   freq_min; /**< Min. allowable load. frequency */
  
  PetscBool     iszip; /**< Flag to check if this load is a zip load model */
};
/**
 * @brief Parses a line in dyr file and reads the model type (string identifier). If the line does not have the registered load model then the type is set to NONE.
 * @param [in] char* line - the line to parse
 * @param [out] char[] modeltype - type of the model (string name of the model) 
 */
extern PetscErrorCode DYNLoadModelReadModelType(char*,char[]);
/**
 * @brief Copies the dynload model
 * @param [in] DYNLoadModel dynloadin - the dynload model to copy from
 * @param [out] DYNLoadModel dynloadout - the dynload model to copy to
 * Notes:
 * If a new method is created for the dynloadmodel then that should be copied too.
 * The underlying implementation (dynload->data) data pointer is copied.
 */
extern PetscErrorCode DYNLoadModelCopy(DYNLoadModel,DYNLoadModel);
/**
 * @brief Registers all built-in dynamic load models
 */
extern PetscErrorCode DYNLoadModelRegisterAll(void);
/**
 * @brief Set event monitors for this load model
 * Notes:
 * Currently, only min. and max. frequency limits implemented
*/
extern PetscErrorCode DYNLoadModelSetEventMonitor(DYNLoadModel);

#endif
