/**
 * @file dynexcmodelsimpl.h
 * @brief Private header files that define the base exciter model object
 */
#ifndef DYNEXCMODELSIMPL_H
#define DYNEXCMODELSIMPL_H

#include <ps.h>
#include <dynexcmodels.h>
#include <dyngenmodels.h>
#include <dynstabmodels.h>

#define DYNEXCMODELSMAX 10

#define DYNEXCMONITORSMAX 10

extern PetscInt nexcmodelsregistered;
extern PetscBool DYNExcModelRegisterAllCalled;

/** 
 * @brief Struct to keep track of registered exciter models. 
 */
struct _p_DYNExcModelList{
  char     name[16]; /**< Name of the model */
  PetscInt key; /**< Key returned by DMNetworkRegisterComponent */
  PetscInt sizeofstruct; /**< size of struct obtained by sizeof() */
  PetscErrorCode (*create)(DYNExcModel); /**< Create/Instantiation routine */
};

extern struct _p_DYNExcModelList DYNExcModelList[DYNEXCMODELSMAX];
/** 
 * @brief Operations available with exciter models 
 */
struct _p_DYNExcModelOps{
  PetscErrorCode (*getnvar)(DYNExcModel,PetscInt*); /**< Return number of variables for the exciter model */
  PetscErrorCode (*readdata)(DYNExcModel,char*); /**< Read data for this model */
  PetscErrorCode (*getbusnumid)(DYNExcModel,PetscInt*,char**); /**< Get the bus number and generator id */
  PetscErrorCode (*destroy)(DYNExcModel); /**< Destroy this exciter model */
  PetscErrorCode (*getsizeof)(DYNExcModel,PetscInt*);/**< Get size of the data structure (obtained by sizeof()) */
  PetscErrorCode (*setinitialconditions)(DYNExcModel,PetscScalar,PetscScalar,PetscScalar*); /**< Set initial conditions for the model */
  PetscErrorCode (*setinitialconditionsp)(DYNExcModel,PetscScalar,PetscScalar,PetscScalar,PetscScalar,PetscInt,PetscInt,PetscInt,PetscInt,Mat,PetscInt); /**< Set the differentiation of initial conditions for the model */
  PetscErrorCode (*daerhsfunction)(DYNExcModel,PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*); /**< Return the differential-algebraic equations for this exciter model */
  PetscErrorCode (*daerhsjacobian)(DYNExcModel,Mat,PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscInt,PetscInt[]); /**< Set the Jacobian portion for this exciter model */
  PetscErrorCode (*daerhsjacobianp)(DYNExcModel,PetscReal,const PetscScalar*,Mat,PetscInt,PetscInt,PetscInt); /**< Set the Jacobian for the adjoint sensitivity equations for this exciter model */
  PetscErrorCode (*daefwdrhsjacobianp)(DYNExcModel,PetscReal,const PetscScalar*,Vec*,PetscInt,PetscInt,PetscInt); /**< Set the Jacobian for the forward sensitivity equations for this exciter model */
  PetscErrorCode (*getequationtypes)(DYNExcModel,PetscInt*,PetscInt*,PetscInt*); /**< Set the number of differential and algebraic equations */
  PetscErrorCode (*getfieldvoltage)(DYNExcModel,PetscReal,PetscScalar*,PetscScalar*,PetscInt*,PetscScalar[],PetscInt[]); /**< Return the field voltage and its derivatives (optional) */
  PetscErrorCode (*setevent)(DYNExcModel,PetscInt*,PetscInt*,PetscBool*,PetscErrorCode (**)(DYNExcModel,PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*),PetscErrorCode (**)(DYNExcModel,PetscInt,PetscInt[],PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscBool,PetscBool*),PetscErrorCode (**)(DYNExcModel,PetscInt,PetscInt[],PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*,PetscScalar*,PetscInt,PetscInt,PetscInt,PetscInt)); /**< Set events for this exciter model */
};

typedef struct _p_DYNExcModelOps DYNExcModelOps;

/**
 * @brief private struct for base dynamic exciter model
 */
struct _p_DYNExcModel{
  char             type[16]; /**< model type */
  void*            data;     /**< model data */
  DYNExcModelOps   ops;      /**< common operations for each model */
  PetscInt         nvar;     /**< Number of variables for the model */
  PetscInt         startloc; /**< Starting location for the variables for this model */
  PetscInt         ndiff;    /**< Number of differential equations for this model */
  PetscInt         nalg;     /**< Number of algebraic equations for this model */
  PetscInt         *eqtypes; /**< Type of equation: DIFF_EQ for Differential, ALG_EQ for algebraic */
  DYNGenModel      dyngen;   /**< The generator model associated with this exciter model */
  DYNStabModel     dynstab;  /**< The stabilizer model associated with this exciter model */
  PetscInt         numMonitors; /**< Number of event monitors */
  PetscInt         direction[DYNEXCMONITORSMAX];
  PetscBool        terminate[DYNEXCMONITORSMAX];
  PetscErrorCode   (*eventfcn)(DYNExcModel,PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*);
  PetscErrorCode   (*posteventfcn)(DYNExcModel,PetscInt,PetscInt[],PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscBool,PetscBool*);
  PetscErrorCode   (*posteventdgamma)(DYNExcModel,PetscInt,PetscInt[],PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*,PetscScalar*,PetscInt,PetscInt,PetscInt,PetscInt);
};
/**
 * @brief ParsesParses a line in dyr file and reads the model type (string identifier). If the line does not have the registered exciter model then the type is set to NONE.
 * @param [in] char* line - the line to parse
 * @param [out] char[] modeltype - type of the model (string name of the model) (IMPLIMENT FILE USE DIFFERENT NAME: gentype)
 */
extern PetscErrorCode DYNExcModelReadModelType(char*,char[]);
/**
 * @brief Copies the dynexc model
 * @param [in] DYNExcModel dynexcin - the dynexc model to copy from
 * @param [out] DYNExcModel dynexcout - the dynexc model to copy to
 * Notes:
 * If a new method is created for the dynexcmodel then that should be copied too.
 * The underlying implementation (dynexc->data) data pointer is copied.
 */
extern PetscErrorCode DYNExcModelCopy(DYNExcModel,DYNExcModel);
/**
 * @brief Registers all built-in exciter models
 */
extern PetscErrorCode DYNExcModelRegisterAll(void);
/**
 * @brief Sets the event info for the exciter model
 * @param [in] DYNExcModel dynexc - the dynamic exciter model
 */
extern PetscErrorCode DYNExcModelSetEventMonitor(DYNExcModel);

#endif
