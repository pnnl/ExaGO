/**
 * @file dyngenmodels.h
 * @brief Public header file defining the public API for dynamics generator models
 *
 */

#ifndef DYNGENMODELS_H
#define DYNGENMODELS_H

#include <petsc.h>
#include <common.h>

typedef const char* DYNGenModelType;

#define DYNGENNONE "'NONE'"

typedef struct _p_DYNGenModel *DYNGenModel;

#include <dynexcmodels.h>
#include <dynturbgovmodels.h>

/**
 * @brief Gets the bus number and generator ID for this dynamic generator model from the data file. Calls the underlying implementation
 * @param [in] DYNGenModel dyngen - the dynamic generator model object
 * @param [in] PetscInt* busnum - bus number
 * @param [in] char** genid  - the generator ID
 */
PETSC_EXTERN PetscErrorCode DYNGenModelGetBusnumID(DYNGenModel,PetscInt*,char**);
/**
 * @brief NOT IMPLIMENTED
 */
PETSC_EXTERN PetscErrorCode DYNGenModelGetBusnumIDType(char*,PetscInt*,char*,char[]);
/**
 * @brief Sets the type of the dynamic generator model
 * @param [in] DYNGenModel dyngen - the DYNGenModel object
 * @param [in] DYNGenModelType modeltype - type of the model
 */
PETSC_EXTERN PetscErrorCode DYNGenModelSetType(DYNGenModel,DYNGenModelType);
/**
 * @brief Parses the given line and reads the data for the dynamic generator model
 * @param [in] DYNGenModel dyngen - the DYNGenModel object
 * @param [in] char* line   - the file from the line
 * @param [in] PetscScalar mbase  - the machine base
 * @param [in] PetscScalar sbase  - the system base
 * Notes: 
 * The type of generator model should be set prior to calling this routine via DYNGenModelSetType()
 * The conversion from machine base to system base,if applicable, should be done in this rouine
 */
PETSC_EXTERN PetscErrorCode DYNGenModelReadData(DYNGenModel,char*,PetscScalar,PetscScalar);
/**
 * @brief Registers a dynamic generator model
 * @param [in] const char[] sname     - model name (string)
 * @param [in] PetscErrorCode (*)(DYNGenModel) createfunction  - the class constructor
 * @param [in] PetscInt sizeofstruct  - size of object (obtained with sizeof())
 */
PETSC_EXTERN PetscErrorCode DYNGenModelRegister(const char[],PetscErrorCode (*)(DYNGenModel),PetscInt);
/**
 * @brief Destroys the dyngen model object.
 * @param [in] DYNGenModel dyngen - the dynamic generator model
 * Notes: 
 * This routine is only called for dynamic simulation. It calls the destroy method on
 * the child class to free its data
 */
PETSC_EXTERN PetscErrorCode DYNGenModelDestroy(DYNGenModel);
/**
 * @brief Gets the size of the implementation type (obtained using sizeof())
 * @param [in] DYNGenModel dyngen - the DYNGenModel object
 * @param [out] PetscInt* size - the size of the derived class object
 */
PETSC_EXTERN PetscErrorCode DYNGenModelGetSizeof(DYNGenModel,PetscInt*);
/**
 * @brief Returns the number of variables for this dynamic generator model
 * @param [in] DYNGenModel dyngen - the DYNGenModel object
 * @param [out] PetscInt* Nvar - number of variables (dofs) for this generator model
 */
PETSC_EXTERN PetscErrorCode DYNGenModelGetNvar(DYNGenModel,PetscInt*);
/**
 * @brief Sets the initial conditions for this generator model
 * @param [in] DYNGenModel dyngen - the DYNGenModel object
 * @param [in] PetscScalar Pg     - the real generator power output
 * @param [in] PetscScalar Qg     - the reactive generator power output
 * @param [in] PetscScalar VD     - real component of the bus voltage
 * @param [in] PetscScalar VQ     - imaginary component of the bus voltage
 * @param [out] PetscScalar* x   - the array of the variables for the bus on which this generator is incident.
 * Notes: 
 * The initial conditions are populated in the array x
 * Other constants, parameters can be also set during this function call in
 * its implementation struct (for e.g. setting mechanical power Pm etc.)
 * The locations to insert the variables for this generator should be obtained from
 * DYNGenModelGetFirstVariableLocation()   
 * x contains all the variables for the generator bus
 */
PETSC_EXTERN PetscErrorCode DYNGenModelSetInitialConditions(DYNGenModel,PetscScalar,PetscScalar,PetscScalar,PetscScalar,PetscScalar*);
PETSC_EXTERN PetscErrorCode DYNGenModelSetInitialConditionsP(DYNGenModel,PetscScalar,PetscScalar,PetscScalar,PetscScalar,PetscInt,PetscInt,PetscInt,PetscInt,Mat,PetscInt);
PETSC_EXTERN PetscErrorCode DYNGenModelSetInitialFieldVoltageDiff(DYNGenModel,PetscScalar,PetscScalar,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*,PetscScalar*,PetscScalar*);
/**Computes the RHS of the DAE function for the generator model [f(x,y);g(x,y)] and also returns the generator currents in network refernce frame
 * @param [in] DYNGenModel dyngen - the DYNGenModel object
 * @param [in] PetscReal t      - the current time
 * @param [in] PetscScalar VD         - real-part of complex bus voltage
 * @param [in] PetscScalar VQ         - imaginary-part of complex bus voltage
 * @param [in] PetscScalar* x      - array of the variables for the bus on which this generator is incident
 * @param [out] PetscScalar* f      - array of rhs of DAE equations for this generator
 * @param [out] PetscScalar* IGD    - real-part of generator current in network reference frame
 * @param [out] PetscScalar* IGQ    - imaginary-part of generator current in network reference frame
 * Notes: 
 * The locations for the variables (and the corresponding locations to insert entries in f)
 * for this generator model should be obtained by DYNGenModelGetFirstVariableLocation()
 */
PETSC_EXTERN PetscErrorCode DYNGenModelDAERHSFunction(DYNGenModel,PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*,PetscScalar*,PetscScalar*);
/**
 * @brief Computes the RHS Jacobian of the DAE equations
 * @param [in] DYNGenModel dyngen - the DYNGenModel object
 * @param [out] Mat J          - the Jacobian matrix
 * @param [in] PetscReal t          - the current time
 * @param [in] PetscScalar VD         - real-part of complex bus voltage
 * @param [in] PetscScalar VQ         - imaginary-part of complex bus voltage
 * @param [in] PetscScalar* xdyngen    - generator variables
 * @param [in] PetscInt dynlocglob - starting location of generator variables in the global vector 
 * @param [in] PetscInt[] V_loc      - global location of VD and VQ V_loc[0] = VD_loc V_loc[1] = VQ_loc
 * @param [in] PetscInt[] I_loc      - global location of ID and IQ I_loc[0] = ID_loc I_loc[1] = IQ_loc
 * Notes: 
 * For each bus, the network equations are ordered as [IQ_loc;ID_loc]
 */
PETSC_EXTERN PetscErrorCode DYNGenModelDAERHSJacobian(DYNGenModel,Mat,PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscInt,PetscInt[],PetscInt[]);
/**
 * @brief Gets the number and indices of differential and algebraic equations for the generator model
 * @param [in] DYNGenModel dyngen - the DYNGenModel object
 * @param [out] PetscInt* ndiff - number of differential equations
 * @param [out] PetscInt* nalg  - number of algebraic equations
 * @param [out] PetscInt* eqtype - an array of size nvar that has [eqtype[i] = DIFF_EQ if equation i is a differential] or [eqtype[i] = ALG_EQ if equation  i is an algebraic]
 * Notes: 
 * The user does not need to create the array vartype
 */
PETSC_EXTERN PetscErrorCode DYNGenModelGetEquationTypes(DYNGenModel,PetscInt*,PetscInt*,PetscInt*);
/**
 * @brief Returns the model data associated with the particular implementation type
 * @param [in] DYNGenModel dyngen - the DYNGenModel object
 * @param [out] void** modeldata - the model data struct
 */
PETSC_EXTERN PetscErrorCode DYNGenModelGetModelData(DYNGenModel,void**);
/**
 * @brief Returns the field voltage at t=t0
 * @param [in] DYNGenModel dyngen - the DYNGenModel object
 * @param [out] PetscScalar* Efd0 - Initial field voltage at t=t0
 */
PETSC_EXTERN PetscErrorCode DYNGenModelGetInitialFieldVoltage(DYNGenModel,PetscScalar*);
/**
 * @brief Returns the mechanical power at t=t0
 * @param [in] DYNGenModel dyngen - the DYNGenModel object
 * @param [out] PetscScalar* Pmech0 - Initial mechanical power at t=t0
 */
PETSC_EXTERN PetscErrorCode DYNGenModelGetInitialMechanicalPower(DYNGenModel,PetscScalar*);
/**
 * @brief Returns the field current Ifd
 * @param [in] DYNGenModel dyngen - the DYNGenModel object
 * @param [out] PetscInt* Ifd - Field current Ifd
 */
PETSC_EXTERN PetscErrorCode DYNGenModelGetFieldCurrent(DYNGenModel,PetscScalar*,PetscScalar*);
/**
 * @brief Gets the location for the first variable for the dyngen model from the bus variables array
 * @param [in] DYNGenModel dyngen - the DYNGenModel object
 * @param [out] PetscInt* loc - the location for the first variable
 */
PETSC_EXTERN PetscErrorCode DYNGenModelGetFirstVariableLocation(DYNGenModel,PetscInt*);
/**
 * @brief Returns the exciter associated with this generator
 * @param [in] DYNGenModel dyngen - the DYNGenModel object
 * @param [out] DYNExcModel* dynexc - the dynexc object
 */
PETSC_EXTERN PetscErrorCode DYNGenModelGetDynExc(DYNGenModel,DYNExcModel*);
/**
 * @brief Returns the turbine governor associated with this generator
 * @param [in] DYNGenModel dyngen - the DYNGenModel object
 * @param [out] DYNTurbgovModel* dynturbgov - the dynturbgov object
 */
PETSC_EXTERN PetscErrorCode DYNGenModelGetDynTurbgov(DYNGenModel,DYNTurbgovModel*);
/**
 * @brief Returns the speed deviation and its time derivative, (optional) of the generator at the current time.
 * @param [in] DYNGenModel dyngen - the DYNGenModel object
 * @param [in] PetscReal t      - the current time
 * @param [in] const PetscScalar* xdyn   - the state variables for this generator at time t
 * @param [out] PetscScalar* spd     - the generator speed deviation
 * @param [out] PetscScalar* dspd_dt - time-derivative of generator speed deviation (set to NULL if not required)
 */
PETSC_EXTERN PetscErrorCode DYNGenModelGetSpeedDeviation(DYNGenModel,PetscReal,const PetscScalar*,PetscScalar*,PetscScalar*);
/**
 * @brief Gets the location of the speed deviation variable for the generator model relative to the bus variable starting location
 * @param [in] DYNGenModel dyngen - the DYNGenModel object
 * @param [out] PetscInt* loc - the location of the speed deviation
 * Notes:
 * The location of the speed deviation is relative to the bus variables. So, one should first obtain
 * the first variable location for the generator model implementation using DYNGenModelGetFirstVariableLocation() 
 * and then add the location of the dw variable to return loc. For example, if dw is the 3rd variable for this
 * generator model implementation then the following code should be used in the implementation
 * ierr = DYNGenModelGetFirstVariable(dyngen,&firstvarloc);
 * *loc = firstvarloc + 2;
 */
PETSC_EXTERN PetscErrorCode DYNGenModelGetSpeedDeviationLocation(DYNGenModel,PetscInt*);
/**
 * @brief Gets the global location (location in the global vector) of the field voltage for the exciter model
 * @param [in] DYNGenModel dyngen - the DYNGenModel object
 * @param [out] PetscInt* locglob - the global location of the speed deviation variable
 */
PETSC_EXTERN PetscErrorCode DYNGenModelGetSpeedDeviationGlobalLocation(DYNGenModel,PetscInt*);
/**
 * @brief Returns the frequency of the generator at the current time.
 * @param [in] DYNGenModel dyngen - the DYNGenModel object
 * @param [in] PetscReal t      - the current time
 * @param [in] const PetscScalar* xdyn   - the state variables for this generator at time t
 * @param [out] PetscScalar* frequency   - the generator frequency in Hz
 */
PETSC_EXTERN PetscErrorCode DYNGenModelGetFrequency(DYNGenModel,PetscReal,const PetscScalar*,PetscScalar*);
/**
 * @brief Set the frequency of the generator at the current time.
 * @param [in] DYNGenModel dyngen - the DYNGenModel object
 * @param [in] PetscReal t      - the current time
 * @param [in] const PetscScalar* xdyn   - the state variables for this generator at time t
 * @param [in] PetscScalar* frequency   - the value to be set
 */
PETSC_EXTERN PetscErrorCode DYNGenModelSetFrequency(DYNGenModel,PetscScalar*,PetscScalar);
PETSC_EXTERN PetscErrorCode DYNGenModelGetFrequencyLimits(DYNGenModel,PetscScalar*,PetscScalar*);
/**
 * @brief Returns the partial derivative of the generator frequency w.r.t the machine state that governs it
 * @param [in] DYNGenModel dyngen - the DYNGenModel object
 * @param [in] PetscReal t      - the current time
 * @param [in] const PetscScalar* xdyn   - the state variables for this generator at time t
 * @param [out] PetscScalar* dfreqdstate   - dfreq_dstate
 * @param [out] PetscInt* stateloc      - location w.r.t. to the bus variables
 */
PETSC_EXTERN PetscErrorCode DYNGenModelGetdFreqdState(DYNGenModel,PetscReal,const PetscScalar*,PetscScalar*,PetscInt*);
/**
 * @brief Sets the Jacobian of the DAE w.r.t. parameters
 * @param [in] DYNGenModel dyngen - the DYNGenModel object
 * @param [in] PetscReal t          - the current time
 * @param [in] const PetscScalar* x      - state variables for this bus
 * @param [out] Mat jacP - Matrix of partial derivatives of DYN DAE equations w.r.t. parameters
 * @param [in] PetscInt dynlocglob - starting location of the variables for this dyngen model
 * @param [in] PetscInt Valoc  - location of parameter voltage angle (Va)
 * @param [in] PetscInt Vmloc  - location of parameter voltage magnitude (Vm)
 * @param [in] PetscInt Pgloc  - location of parameter generator real power output (Pg)
 * @param [in] PetscInt Qgloc  - location of parameter generator reactive power output (Qg)
 */
PETSC_EXTERN PetscErrorCode DYNGenModelDAERHSJacobianP(DYNGenModel,PetscReal,const PetscScalar*,Mat,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt);
/**
 * @brief Sets the Jacobian of the DAE w.r.t. parameters
 * @param [in] DYNGenModel dyngen - the DYNGenModel object
 * @param [in] PetscReal t          - the current time
 * @param [in] const PetscScalar* x      - state variables for this bus
 * @param [out] Vec* jacP - Matrix of partial derivatives of DYN DAE equations w.r.t. parameters
 * @param [in] PetscInt dynlocglob - starting location of the variables for this dyngen model
 * @param [in] PetscInt Valoc  - location of parameter voltage angle (Va)
 * @param [in] PetscInt Vmloc  - location of parameter voltage magnitude (Vm)
 * @param [in] PetscInt Pgloc  - location of parameter generator real power output (Pg)
 * @param [in] PetscInt Qgloc  - location of parameter generator reactive power output (Qg)
 */
PETSC_EXTERN PetscErrorCode DYNGenModelDAEFWDRHSJacobianP(DYNGenModel,PetscReal,const PetscScalar*,Vec*,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt);
PETSC_EXTERN PetscErrorCode DYNGenModelSetStatus(DYNGenModel,PetscInt);

PETSC_EXTERN PetscErrorCode DYNGenModelGetCurrent(DYNGenModel,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*,PetscScalar*);
PETSC_EXTERN PetscErrorCode DYNGenModelSetVoltage(DYNGenModel,PetscScalar,PetscScalar);
#endif
