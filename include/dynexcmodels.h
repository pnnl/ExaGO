/**
 * @file dynexcmodels.h
 * @brief Public header file defining the public API for dynamic exciter models
 */

#ifndef DYNEXCMODELS_H
#define DYNEXCMODELS_H

#include <petsc.h>
#include <common.h>
#include <dyngenmodels.h>
#include <dynstabmodels.h>

typedef const char* DYNExcModelType;

#define DYNEXCNONE "'NONE'"

typedef struct _p_DYNExcModel *DYNExcModel;

/**
 * @brief Gets the bus number and exciter ID for this dynamic exciter model from the data file. Calls the underlying implementation
 * @param [in] DYNExcModel dynexc - the dynamic exciter model
 * @param [out] PetscInt* busnum - bus number
 * @param [out] char** genid  - the exciter ID
 */
extern PetscErrorCode DYNExcModelGetBusnumID(DYNExcModel,PetscInt*,char**);
/**
 * @brief NOT IMPLIMENTED
 */
extern PetscErrorCode DYNExcModelGetBusnumIDType(char*,PetscInt*,char*,char[]);
/**
 * @brief Sets the type of the dynamic exciter model
 * @param [in] DYNExcModel dynexc - the dynamic exciter model
 * @param [in] DYNExcModelType modeltype - type of the model
 */
extern PetscErrorCode DYNExcModelSetType(DYNExcModel,DYNExcModelType);
/**
 * @brief Parses the given line and reads the data for the dynamic exciter model
 * @param [in] DYNExcModel dynexc - the dynamic exciter model
 * @param [in] char* line   - the file from the line
 */
extern PetscErrorCode DYNExcModelReadData(DYNExcModel,char*);
/**
 * @brief Registers a dynamic exciter model
 * @param [in] const char[] sname     - model name (string)
 * @param [in] PetscErrorCode (*)(DYNExcModel) createfunction  - the class constructor
 * @param [in] PetscInt sizeofstruct  - size of object (obtained with sizeof())
 */
extern PetscErrorCode DYNExcModelRegister(const char[],PetscErrorCode (*)(DYNExcModel),PetscInt);
/**
 * @brief Destroys the dynexc model object.
 * @param [in] DYNExcModel dynexc - the dynamic exciter model
 * Notes: 
 * This routine is only called for dynamic simulation. It calls the destroy method on
 * the child class to free its data
 */
extern PetscErrorCode DYNExcModelDestroy(DYNExcModel);
/**
 * @brief Gets the size of the implementation type (obtained using sizeof())
 * @param [in] DYNExcModel dynexc - the dynamic exciter model
 * @param [out] PetscInt* size - the size of the derived class object
 */
extern PetscErrorCode DYNExcModelGetSizeof(DYNExcModel,PetscInt*);
/**
 * @brief Returns the number of variables for this dynamic exciter model
 * @param [in] DYNExcModel dynexc - the dynamic exciter model
 * @param [out] PetscInt* Nvar - number of variables (dofs) for this exciter model
 */
extern PetscErrorCode DYNExcModelGetNvar(DYNExcModel,PetscInt*);
/**
 * @brief Sets the initial conditions for this exciter model
 * @param [in] DYNExcModel dynexc - the dynamic exciter model
 * @param [in] PetscScalar VD     - real component of the bus voltage
 * @param [in] PetscScalar VQ     - imaginary component of the bus voltage
 * @param [out] PetscScalar* xdyn   - the initial conditions for this exciter model
 * Notes: 
 * The initial conditions are populated in the xdyn array
 * Other constants, parameters can be also set during this function call in
 * its implementation struct (for e.g. setting reference voltage Vref and others)
 */
extern PetscErrorCode DYNExcModelSetInitialConditions(DYNExcModel,PetscScalar,PetscScalar,PetscScalar*);
extern PetscErrorCode DYNExcModelSetInitialConditionsP(DYNExcModel,PetscScalar,PetscScalar,PetscScalar,PetscScalar,PetscInt,PetscInt,PetscInt,PetscInt,Mat,PetscInt);
/**Computes the RHS of the DAE function for the generator model [f(x,y);g(x,y)] and also returns the generator currents in network refernce frame
 * @param [in] DYNExcModel dynexc     - the dynamic exciter model
 * @param [in] PetscReal t      - the current time
 * @param [in] PetscScalar VD         - real-part of complex bus voltage
 * @param [in] PetscScalar VQ         - imaginary-part of complex bus voltage
 * @param [in] PetscScalar* x      - array of the variables for the bus on which this generator is incident
 * @param [out] PetscScalar* f      - the residual array in which the exciter equation residuals are inserted.\n
 * The location of the exciter variables can be obtained using DYNExcModelGetFirstVariableLocation(). Using the first location, the residuals can be inserted as f[loc] = 1st exciter equation residual, f[loc+1] = 2nd ... and so on
 */
extern PetscErrorCode DYNExcModelDAERHSFunction(DYNExcModel,PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*);
/**
 * @brief Computes the RHS Jacobian of the DAE equations
 * @param [in] DYNExcModel dynexc     - the dynamic exciter model
 * @param [out] Mat J          - the Jacobian matrix
 * @param [in] PetscReal t          - the current time
 * @param [in] PetscScalar VD         - real-part of complex bus voltage
 * @param [in] PetscScalar VQ         - imaginary-part of complex bus voltage
 * @param [in] PetscScalar* xdyngen    - generator variables
 * @param [in] PetscInt dynlocglob - starting location of generator variables in the global vector 
 * @param [in] PetscInt[] V_loc      - global location of VD and VQ V_loc[0] = VD_loc V_loc[1] = VQ_loc
 * Notes: 
 * For each bus, the network equations are ordered as [IQ_loc;ID_loc]
 */
extern PetscErrorCode DYNExcModelDAERHSJacobian(DYNExcModel,Mat,PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscInt,PetscInt[]);
/**
 * @brief Computes the RHS Jacobian of the DAE equations w.r.t. parameters
 * @param [in] DYNExcModel dynexc     - the dynamic exciter model
 * @param [in] PetscReal t          - the current time
 * @param [in] const PetscScalar* xdynexc    - exciter variables
 * @param [out] Mat jacP - Matrix of partial derivatives of DYN DAE equations w.r.t. parameters
 * @param [in] PetscInt dynlocglob - starting location of the variables for this dyngen model
 * @param [in] PetscInt Valoc  - location of parameter voltage angle (Va)
 * @param [in] PetscInt Vmloc  - location of parameter voltage magnitude (Vm)
 */
extern PetscErrorCode DYNExcModelDAERHSJacobianP(DYNExcModel,PetscReal,const PetscScalar*,Mat,PetscInt,PetscInt,PetscInt);
/**
 * @brief Computes the RHS Jacobian of the DAE equations w.r.t. parameters
 * @param [in] DYNExcModel dynexc     - the dynamic exciter model
 * @param [in] PetscReal t          - the current time
 * @param [in] const PetscScalar* xdynexc    - exciter variables
 * @param [out] Vec* jacP - Matrix of partial derivatives of DYN DAE equations w.r.t. parameters
 * @param [in] PetscInt dynlocglob - starting location of the variables for this dyngen model
 * @param [in] PetscInt Valoc  - location of parameter voltage angle (Va)
 * @param [in] PetscInt Vmloc  - location of parameter voltage magnitude (Vm)
 */
extern PetscErrorCode DYNExcModelDAEFWDRHSJacobianP(DYNExcModel,PetscReal,const PetscScalar*,Vec*,PetscInt,PetscInt,PetscInt);
/**
 * @brief Gets the number and indices of differential and algebraic equations for the exciter model
 * @param [in] DYNExcModel dynexc     - the dynamic exciter model
 * @param [out] PetscInt* ndiff - number of differential equations
 * @param [out] PetscInt* nalg  - number of algebraic equations
 * @param [out] PetscInt* eqtype - an array of size nvar that has [eqtype[i] = DIFF_EQ if equation i is a differential] or [eqtype[i] = ALG_EQ if equation  i is an algebraic]
 * Notes: 
 * The user does not need to create the array vartype
 */
extern PetscErrorCode DYNExcModelGetEquationTypes(DYNExcModel,PetscInt*,PetscInt*,PetscInt*);
/**
 * @brief Returns the model data associated with the particular implementation type
 * @param [in] DYNExcModel dynexc - the dynamic exciter model
 * @param [out] void** modeldata - the model data struct
 */
extern PetscErrorCode DYNExcModelGetModelData(DYNExcModel,void**);
/**
 * @brief Returns the generator associated with this exciter
 * @param [in] DYNExcModel dynexc - the dynamic exciter model
 * @param [out] DYNGenModel* dyngen - the dyngen object
 */
extern PetscErrorCode DYNExcModelGetDynGen(DYNExcModel,DYNGenModel*);
/**
 * @brief Returns the stabilizer associated with this exciter
 * @param [in] DYNExcModel dynexc - the dynamic exciter model
 * @param [out] DYNStabModel* dynstab - the dynstab object
 */
extern PetscErrorCode DYNExcModelGetDynStab(DYNExcModel,DYNStabModel*);
/**
 * @brief Gets the location for the first variable for the dynexc model from the bus variables array
 * @param [in] DYNExcModel dynexc - the dynamic exciter model
 * @param [out] PetscInt* loc - the location for the first variable 
 */
extern PetscErrorCode DYNExcModelGetFirstVariableLocation(DYNExcModel,PetscInt*);
/**
 * @brief Returns the exciter field voltage and its partial derivative w.r.t to exciter variables (optional)
 * @param [in] DYNExcModel dynexc - the dynamic exciter model
 * @param [in] PetscReal t      - the current time
 * @param [in] PetscScalar* xdyn   - array of the variables for this bus
 * @param [out] PetscScalar* Efd           - Field voltage for the exciter
 * @param [out] PetscInt* dEfddXexc_num - Number of partial derivatives of Efd w.r.t. exciter variables.
 * @param [out] PetscScalar[] dEfddXexc     - An array of partial derivatives of the field voltage w.r.t. the number of exciter states given by dEfddXexc_num
 * @param [out] PetscInt[] dEfddXexc_loc - Global locations for the partial derivatives
 * Notes: dEfddXexc_num, dEfddXexc and dEfddXexc_loc are optionally returned if the passed-in arrays are not NULL.
 */
extern PetscErrorCode DYNExcModelGetFieldVoltage(DYNExcModel,PetscReal,PetscScalar*,PetscScalar*,PetscInt*,PetscScalar[],PetscInt[]);

#endif
