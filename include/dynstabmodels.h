/**
 * @file dynstabmodels.h
 * @brief Public header file defining the public API for dynamic stabilizer models
 */

#ifndef DYNSTABMODELS_H
#define DYNSTABMODELS_H

#include <petsc.h>
#include <common.h>
#include <dynexcmodels.h>

typedef const char* DYNStabModelType;

#define DYNSTABNONE "'NONE'"

typedef struct _p_DYNStabModel *DYNStabModel;

/**
 * @brief Gets the bus number and stab ID for this dynamic stab model from the data file. Calls the underlying implementation
 * @param [in] DYNStabModel dynstab - the DYNStabModel object
 * @param [out] PetscInt* busnum - bus number
 * @param [out] char** genid  - the stab ID
 */
extern PetscErrorCode DYNStabModelGetBusnumID(DYNStabModel,PetscInt*,char**);
/**
 * @brief NOT IMPLIMENTED
 */
extern PetscErrorCode DYNStabModelGetBusnumIDType(char*,PetscInt*,char*,char[]);
/**
 * @brief Sets the type of the dynamic stab model
 * @param [in] DYNStabModel dynstab - the DYNStabModel object
 * @param [in] DYNStabModelType modeltype - type of the model
 */
extern PetscErrorCode DYNStabModelSetType(DYNStabModel,DYNStabModelType);
/**
 * @brief Parses the given line and reads the data for the dynamic stab model
 * @param [in] DYNStabModel dynstab - the DYNStabModel object
 * @param [in] char* line   - the file from the line
 * Notes: The type of stab model should be set prior to this call via DYNStabModelSetType()
 */
extern PetscErrorCode DYNStabModelReadData(DYNStabModel,char*);
/**
 * @brief Registers a dynamic stab model
 * @param [in] const char[] sname     - model name (string)
 * @param [in] PetscErrorCode (*)(DYNStabModel) createfunction  - the class constructor
 * @param [in] PetscInt sizeofstruct  - size of object (obtained with sizeof())
 */
extern PetscErrorCode DYNStabModelRegister(const char[],PetscErrorCode (*)(DYNStabModel),PetscInt);
/**
 * @brief Gets the size of the implementation type (obtained using sizeof())
 * @param [in] DYNStabModel dynstab - the dynamic stab model
 * Note: This routine is only called for dynamic simulation. It calls the destroy method on the child class to free its data
 */
extern PetscErrorCode DYNStabModelDestroy(DYNStabModel);
/**
 * @brief Gets the bus number and stab ID for this dynamic stab model from the data file. Calls the underlying implementation
 * @param [in] DYNStabModel dynstab - the DYNStabModel object
 * @param [out] PetscInt* size - the size of the derived class object
 */
extern PetscErrorCode DYNStabModelGetSizeof(DYNStabModel,PetscInt*);
/**
 * @brief Returns the number of variables for this dynamic stab model
 * @param [in] DYNStabModel dynstab - the DYNStabModel object
 * @param [out] PetscInt* Nvar - number of variables (dofs) for this stab model
 */
extern PetscErrorCode DYNStabModelGetNvar(DYNStabModel,PetscInt*);
/**
 * @brief Returns the stabilizer signal for the exciter and optionally returns partial derivative w.r.t to the stabilizer variables and speed deviation
 * @param [in] DYNStabModel dynstab - the dynamic stabilizer model
 * @param [in] PetscReal t      - the current time
 * @param [in] PetscScalar* xdyn   - array of the variables for this bus
 * @param [out] PetscScalar* VOTHSG           - Stabilizer signal output
 * @param [out] PetscInt* dVOTHSGdXdyn_num - Number of partial derivatives of VOTHSG w.r.t. the dynamic variables at the bus 
 * @param [out] PetscScalar[] dVOTHSGdXdyn     - An array of partial derivatives of the stabilizer output signal VOTHSG w.r.t. the dynamic variables at the bus
 * @param [out] PetscInt[] dVOTHSGdXdyn_loc - Location of partial derivatives of the stabilizer output signal VOTHSG w.r.t. the dynamic variables at the bus. Need to get the first variable location using DYNStabModelGetFirstVariableLocation() and then add the displacement.
 * Notes: dVOTHSGdXdyn_num, dVOTHSGdXdyn and dVOTHSGdXdyn_loc are optionally returned if the passed-in arrays are not NULL.
 */
extern PetscErrorCode DYNStabModelGetVOTHSG(DYNStabModel,PetscReal,PetscScalar*,PetscScalar*,PetscInt*,PetscScalar[],PetscInt[]);
/**
 * @brief Sets the initial conditions for this stab model
 * @param [in] DYNStabModel dynstab - the DYNStabModel object
 * @param [in] PetscScalar VD     - real component of the bus voltage
 * @param [in] PetscScalar VQ     - imaginary component of the bus voltage
 * @param [out] PetscScalar* xdyn   - the initial conditions for this stab model
 * Notes: The initial conditions are populated in the xdyn array\n
 * Other constants, parameters can be also set during this function call in its implementation struct (for e.g. setting reference voltage Vref and others)
 */
extern PetscErrorCode DYNStabModelSetInitialConditions(DYNStabModel,PetscScalar,PetscScalar,PetscScalar*);
/**
 * @brief Computes the RHS of the DAE function for the stab model [f(x,y);g(x,y)] and also returns the stab currents in network refernce frame
 * @param [in] DYNStabModel dynstab - the dynamic stab model
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
extern PetscErrorCode DYNStabModelDAERHSFunction(DYNStabModel,PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*);
/**
 * @brief Computes the RHS Jacobian of the DAE equations
 * @param [in] DYNStabModel dynstab     - the dynamic stab model
 * @param [out] Mat J - the Jacobian matrix
 * @param [in] PetscReal t - the current time
 * @param [in] PetscScalar VD - real-part of complex bus voltage
 * @param [in] PetscScalar VQ - imaginary-part of complex bus voltage
 * @param [in] PetscScalar* xdynstab    - stab variables
 * @param [in] PetscInt dynlocglob - starting location of stab variables in the global vector 
 * @param [in] PetscInt[] V_loc - global location of VD and VQ V_loc[0] = VD_loc V_loc[1] = VQ_loc
 * Notes: For each bus, the network equations are ordered as [IQ_loc;ID_loc]
 */
extern PetscErrorCode DYNStabModelDAERHSJacobian(DYNStabModel,Mat,PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscInt,PetscInt[]);
/**
 * @brief Computes the RHS Jacobian of the DAE equations w.r.t. parameters
 * @param [in] DYNStabModel dynstab     - the dynamic stab model
 * @param [in] PetscReal t - the current time
 * @param [in] const PetscScalar* xdynstab    - stab variables
 * @param [out] Mat jacP - matrix of partial derivatives of DAE equations w.r.t parameter PG
 * @param [in] PetscInt dynlocglob - starting location of stab variables in the global vector 
 * @param [in] PetscInt PGloc      - location of PG (generator real power) in the parameter vector
 */
extern PetscErrorCode DYNStabModelDAERHSJacobianP(DYNStabModel,PetscReal,const PetscScalar*,Mat,PetscInt,PetscInt);
/**
 * @brief Gets the number and indices of differential and algebraic equations for the stab model
 * @param [in] DYNStabModel dynstab - the DYNStabModel object
 * @param [out] PetscInt* ndiff - number of differential equations
 * @param [out] PetscInt* nalg  - number of algebraic equations
 * @param [out] PetscInt* an array of size nvar that has [eqtype[i] = DIFF_EQ if equation i is a differential] or [eqtype[i] = ALG_EQ if equation  i is an algebraic]
 * NOTES: The user does not need to create the array vartype
 */
extern PetscErrorCode DYNStabModelGetEquationTypes(DYNStabModel,PetscInt*,PetscInt*,PetscInt*);
/**
 * @brief Returns the model data associated with the particular implementation type
 * @param [in] DYNStabModel dynstab - the DYNStabModel object
 * @param [out] void** modeldata - the model data struct
 */
extern PetscErrorCode DYNStabModelGetModelData(DYNStabModel,void**);
/**
 * @brief Returns the generator associated with this stab
 * @param [in] DYNStabModel dynstab - the DYNStabModel object
 * @param [out] DYNGenModel* dyngen - the dyngen object
 */
extern PetscErrorCode DYNStabModelGetDynGen(DYNStabModel,DYNGenModel*);
/**
 * @brief Gets the location for the first variable for the dynstab model from the bus variables array
 * @param [in] DYNStabModel dynstab - the DYNStabModel object
 * @param [out] PetscInt* loc - the location for the first variable 
 */
extern PetscErrorCode DYNStabModelGetFirstVariableLocation(DYNStabModel,PetscInt*);

#endif
