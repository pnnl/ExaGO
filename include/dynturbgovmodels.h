/**
 * @file dynturbgovmodels.h
 * @brief Public header file defining API for dynamic turbine-governor models
 */

#ifndef DYNTURBGOVMODELS_H
#define DYNTURBGOVMODELS_H

#include <petsc.h>
#include <common.h>
#include <dyngenmodels.h>

typedef const char* DYNTurbgovModelType;

#define DYNTURBGOVNONE "'NONE'"

typedef struct _p_DYNTurbgovModel *DYNTurbgovModel;

/**
 * @brief Gets the bus number and turbgov ID for this dynamic turbgov model from the data file. Calls the underlying implementation
 * @param [in] DYNTurbgovModel dynturbgov - the dynamic turbgov model
 * @param [out] PetscInt* busnum - bus number
 * @param [out] char** genid  - the turbgov ID
 */
extern PetscErrorCode DYNTurbgovModelGetBusnumID(DYNTurbgovModel,PetscInt*,char**);
/**
 * @brief NOT IMPLIMENTED
 */
extern PetscErrorCode DYNTurbgovModelGetBusnumIDType(char*,PetscInt*,char*,char[]);
/**
 * @brief Sets the type of the dynamic turbgov model
 * @param [in] DYNTurbgovModel dynturbgov - the dynamic turbgov model
 * @param [in] DYNTurbgovModelType modeltype - type of the model
 */
extern PetscErrorCode DYNTurbgovModelSetType(DYNTurbgovModel,DYNTurbgovModelType);
/**
 * @brief Parses the given line and reads the data for the dynamic turbgov model
 * @param [in] DYNTurbgovModel dynturbgov - the dynamic turbgov model
 * @param [in] char* line   - the file from the line
 * @param [unknown] PetscScalar
 * @param [unknown] PetscScalar
 * Notes: The type of turbgov model should be set prior to this call via DYNTurbgovModelSetType()
 */
extern PetscErrorCode DYNTurbgovModelReadData(DYNTurbgovModel,char*,PetscScalar,PetscScalar);
/**
 * @brief Registers a dynamic turbgov model
 * @param [in] const char[] sname     - model name (string)
 * @param [in] PetscErrorCode (*)(DYNTurbgovModel) createfunction  - the class constructor
 * @param [in] PetscInt sizeofstruct  - size of object (obtained with sizeof())
 */
extern PetscErrorCode DYNTurbgovModelRegister(const char[],PetscErrorCode (*)(DYNTurbgovModel),PetscInt);
/**
 * @brief Destroys the dynturbgov model object.
 * @param [in] DYNTurbgovModel dynturbgov - the dynamic turbgov model
 */
extern PetscErrorCode DYNTurbgovModelDestroy(DYNTurbgovModel);
/**
 * @brief Gets the size of the implementation type (obtained using sizeof())
 * @param [in] DYNTurbgovModel dynturbgov - the dynamic turbgov model
 * @param [out] PetscInt* size - the size of the derived class object
 */
extern PetscErrorCode DYNTurbgovModelGetSizeof(DYNTurbgovModel,PetscInt*);
/**
 * @brief Returns the number of variables for this dynamic turbgov model
 * @param [in] DYNTurbgovModel dynturbgov - the dynamic turbgov model
 * @param [out] PetscInt* Nvar - number of variables (dofs) for this turbgov model
 */
extern PetscErrorCode DYNTurbgovModelGetNvar(DYNTurbgovModel,PetscInt*);
/**
 * @brief Returns the mechanical power input for the machine and its partial derivative w.r.t to the machine speed deviation
 * @param [in] DYNTurbgovModel dynturbgov - the dynamic turbgov model
 * @param [in] PetscReal t - the current time
 * @param [in] PetscScalar* xdyn - array of the variables for this bus
 * @param [out] PetscScalar* Pmech - Mechanical power input for the machine
 * @param [out] PetscScalar* dPmechddw- Partial derivative of Pmech w.r.t. machine speed deviation
 * @param [out] PetscInt* dPmechdXtgov_num - Number of partial derivatives computed.
 * @param [out] PetscScalar[] dPmechdXtgov - An array of partial derivatives of the mechanical power w.r.t. the number of turbine governor states given by dPmechdXtgov_num
 * @param [out] PetscInt[] dPmechdXtgov_loc - Global locations for the partial derivatives
 * Notes: dPmechdXtgov_num, dPmechdXtgov and dPmechdXtgov_loc are optionally returned if the passed-in arrays are not NULL.
 */
extern PetscErrorCode DYNTurbgovModelGetMechanicalPower(DYNTurbgovModel,PetscReal,PetscScalar*,PetscScalar*,PetscScalar*,PetscInt*,PetscScalar[],PetscInt[]);
/**
 * @brief Sets the initial conditions for this turbgov model
 * @param [in] DYNTurbgovModel dynturbgov - the dynamic turbgov model
 * @param [in] PetscReal VD     - real component of the bus voltage
 * @param [in] PetscScalar VQ     - imaginary component of the bus voltage
 * @param [out] PetscScalar* xdyn   - the initial conditions for this turbgov model
 * Notes: 
 * The initial conditions are populated in the xdyn array
 * Other constants, parameters can be also set during this function call in
 * its implementation struct (for e.g. setting reference voltage Vref and others)
 */
extern PetscErrorCode DYNTurbgovModelSetInitialConditions(DYNTurbgovModel,PetscScalar,PetscScalar,PetscScalar*);
extern PetscErrorCode DYNTurbgovModelSetInitialConditionsP(DYNTurbgovModel,PetscScalar,PetscScalar,PetscScalar,PetscScalar,PetscInt,PetscInt,PetscInt,PetscInt,Mat,PetscInt);
/**
 * @brief Computes the RHS of the DAE function for the turbgov model [f(x,y);g(x,y)] and also returns the turbgov currents in network refernce frame
 * @param [in] DYNTurbgovModel dynturbgov - the dynamic turbgov model
 * @param [in] PetscReal t - the current time
 * @param [in] PetscScalar x - array of variables for the bus on which the generator is incident
 * @param [in] PetscScalar VD - real-part of bus voltage
 * @param [in] PetscScalar* VQ - imaginary part of bus voltage
 * @param [out] PetscScalar* f      - the residual array in which the turbgov equation residuals are inserted. The location
 *         of the turbgov variables can be obtained using DYNTurbgovModelGetFirstVariableLocation(). Using
 *         the first location, the residuals can be inserted as
 *         f[loc] = 1st turbgov equation residual
 *         f[loc+1] = 2nd ..
 *       ... and so on
 */
extern PetscErrorCode DYNTurbgovModelDAERHSFunction(DYNTurbgovModel,PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*);
/**
 * @brief Computes the RHS Jacobian of the DAE equations
 * @param [in] DYNTurbgovModel dynturbgov - the dynamic turbgov model
 * @param [out] Mat J - the Jacobian matrix
 * @param [in] PetscReal t - the current time
 * @param [in] PetscScalar VD - real-part of complex bus voltage
 * @param [in] PetscScalar VQ - imaginary-part of complex bus voltage
 * @param [in] PetscScalar* xdynturbgov - turbgov variables
 * @param [in] PetscInt dynlocglob - starting location of turbgov variables in the global vector 
 * @param [in] PetscInt[] V_loc - global location of VD and VQ V_loc[0] = VD_loc V_loc[1] = VQ_loc
 * Notes: For each bus, the network equations are ordered as [IQ_loc;ID_loc]
 */
extern PetscErrorCode DYNTurbgovModelDAERHSJacobian(DYNTurbgovModel,Mat,PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscInt,PetscInt[]);
/**
 * @brief Computes the RHS Jacobian of the DAE equations w.r.t. parameters
 * @param [in] DYNTurbgovModel dynturbgov - the dynamic turbgov model
 * @param [in] PetscReal t - the current time
 * @param [in] const PetscScalar* xdynturbgov - turbgov variables
 * @param [out] Mat jacP - matrix of partial derivatives of DAE equations w.r.t parameter PG
 * @param [in] PetscInt dynlocglob - starting location of turbgov variables in the global vector 
 * @param [in] PetscInt PGloc      - location of PG (generator real power) in the parameter vector
 */
extern PetscErrorCode DYNTurbgovModelDAERHSJacobianP(DYNTurbgovModel,PetscReal,const PetscScalar*,Mat,PetscInt,PetscInt);
/**
 * @brief Gets the number and indices of differential and algebraic equations for the turbgov model
 * @param [in] DYNTurbgovModel dynturbgov - the dynamic turbgov model
 * @param [out] PetscInt* ndiff - number of differential equations
 * @param [out] PetscInt* nalg  - number of algebraic equations
 * @param [out] PetscInt* eqtype - an array of size nvar that has [eqtype[i] = DIFF_EQ if equation i is a differential] or [eqtype[i] = ALG_EQ if equation  i is an algebraic]
 * NOTES: The user does not need to create the array vartype
 */
extern PetscErrorCode DYNTurbgovModelGetEquationTypes(DYNTurbgovModel,PetscInt*,PetscInt*,PetscInt*);
/**
 * @brief Returns the model data associated with the particular implementation type
 * @param [in] DYNTurbgovModel dynturbgov - the dynamic turbgov model
 * @param [out] void** modeldata - the model data struct
 */
extern PetscErrorCode DYNTurbgovModelGetModelData(DYNTurbgovModel,void**);
/**
 * @brief Returns the generator associated with this turbgov
 * @param [in] DYNTurbgovModel dynturbgov - the dynamic turbgov model
 * @param [out] DYNGenModel* dyngen - the dyngen object
 */
extern PetscErrorCode DYNTurbgovModelGetDynGen(DYNTurbgovModel,DYNGenModel*);
/**
 * @brief Gets the location for the first variable for the dynturbgov model from the bus variables array
 * @param [in] DYNTurbgovModel dynturbgov - the dynamic turbgov model
 * @param [out] PetscInt* loc - the location for the first variable 
 */
extern PetscErrorCode DYNTurbgovModelGetFirstVariableLocation(DYNTurbgovModel,PetscInt*);

#endif
