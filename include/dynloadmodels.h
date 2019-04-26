/**
 * @file dynloadmodels.h
 * @brief Public header file defining the public API for dynamics load models
 *
 */

#ifndef DYNLOADMODELS_H
#define DYNLOADMODELS_H

#include <petsc.h>
#include <common.h>

typedef const char* DYNLoadModelType;

#define DYNLOADNONE "'NONE'"

typedef struct _p_DYNLoadModel *DYNLoadModel;

/**
 * @brief Gets the bus number and load ID for this dynamic load model from the data file. Calls the underlying implementation
 * @param [in] DYNLoadModel dynload - the dynamic load model object
 * @param [in] PetscInt* busnum - bus number
 * @param [in] char** loadid  - the load ID
 */
extern PetscErrorCode DYNLoadModelGetBusnumID(DYNLoadModel,PetscInt*,char**);
/**
 * @brief NOT IMPLIMENTED
 */
extern PetscErrorCode DYNLoadModelGetBusnumIDType(char*,PetscInt*,char*,char[]);
/**
 * @brief Sets the type of the dynamic load model
 * @param [in] DYNLoadModel dynload - the DYNLoadModel object
 * @param [in] DYNLoadModelType modeltype - type of the model
 */
extern PetscErrorCode DYNLoadModelSetType(DYNLoadModel,DYNLoadModelType);
/**
 * @brief Parses the given line and reads the data for the dynamic load model
 * @param [in] DYNLoadModel dynload - the DYNLoadModel object
 * @param [in] char* line   - the file from the line
 * @param [in] PetscScalar mbase  - the machine base
 * @param [in] PetscScalar sbase  - the system base
 * Notes: 
 * The type of load model should be set prior to calling this routine via DYNLoadModelSetType()
 * The conversion from machine base to system base,if applicable, should be done in this rouine
 */
extern PetscErrorCode DYNLoadModelReadData(DYNLoadModel,char*,PetscScalar,PetscScalar);
/**
 * @brief Registers a dynamic load model
 * @param [in] const char[] sname     - model name (string)
 * @param [in] PetscErrorCode (*)(DYNLoadModel) createfunction  - the class constructor
 * @param [in] PetscInt sizeofstruct  - size of object (obtained with sizeof())
 */
extern PetscErrorCode DYNLoadModelRegister(const char[],PetscErrorCode (*)(DYNLoadModel),PetscInt);
/**
 * @brief Destroys the dynload model object.
 * @param [in] DYNLoadModel dynload - the dynamic load model
 * Notes: 
 * This routine is only called for dynamic simulation. It calls the destroy method on
 * the child class to free its data
 */
extern PetscErrorCode DYNLoadModelDestroy(DYNLoadModel);
/**
 * @brief Gets the size of the implementation type (obtained using sizeof())
 * @param [in] DYNLoadModel dynload - the DYNLoadModel object
 * @param [out] PetscInt* size - the size of the derived class object
 */
extern PetscErrorCode DYNLoadModelGetSizeof(DYNLoadModel,PetscInt*);
/**
 * @brief Returns the number of variables for this dynamic loaderator model
 * @param [in] DYNLoadModel dynload - the DYNLoadModel object
 * @param [out] PetscInt* Nvar - number of variables (dofs) for this load model
 */
extern PetscErrorCode DYNLoadModelGetNvar(DYNLoadModel,PetscInt*);
/**
 * @brief Sets the initial conditions for this load model
 * @param [in] DYNLoadModel dynload - the DYNLoadModel object
 * @param [in] PetscScalar Pd     - the real load power input
 * @param [in] PetscScalar Qd     - the reactive load power output
 * @param [in] PetscScalar VD     - real component of the bus voltage
 * @param [in] PetscScalar VQ     - imaginary component of the bus voltage
 * @param [out] PetscScalar* x   - the array of the variables for the bus on which this load is incident.
 * Notes: 
 * The initial conditions are populated in the array x
 * Other constants, parameters can be also set during this function call in
 * its implementation struct (for e.g. setting mechanical power Pm etc.)
 * The locations to insert the variables for this load should be obtained from
 * DYNLoadModelGetFirstVariableLocation()   
 * x contains all the variables for the load bus
 */
extern PetscErrorCode DYNLoadModelSetInitialConditions(DYNLoadModel,PetscScalar,PetscScalar,PetscScalar,PetscScalar,PetscScalar*);
/**Computes the RHS of the DAE function for the load model [f(x,y);g(x,y)] and also returns the load currents in network refernce frame
 * @param [in] DYNLoadModel dynload - the DYNLoadModel object
 * @param [in] PetscReal t      - the current time
 * @param [in] PetscScalar VD         - real-part of complex bus voltage
 * @param [in] PetscScalar VQ         - imaginary-part of complex bus voltage
 * @param [in] PetscScalar* x      - array of the variables for the bus on which this load is incident
 * @param [out] PetscScalar* f      - array of rhs of DAE equations for this load
 * @param [out] PetscScalar* ILD    - real-part of load current in network reference frame
 * @param [out] PetscScalar* ILQ    - imaginary-part of load current in network reference frame
 * Notes: 
 * The locations for the variables (and the corresponding locations to insert entries in f)
 * for this load model should be obtained by DYNLoadModelGetFirstVariableLocation()
 */
extern PetscErrorCode DYNLoadModelDAERHSFunction(DYNLoadModel,PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscScalar*,PetscScalar*,PetscScalar*);
/**
 * @brief Computes the RHS Jacobian of the DAE equations
 * @param [in] DYNLoadModel dynload - the DYNLoadModel object
 * @param [out] Mat J          - the Jacobian matrix
 * @param [in] PetscReal t          - the current time
 * @param [in] PetscScalar VD         - real-part of complex bus voltage
 * @param [in] PetscScalar VQ         - imaginary-part of complex bus voltage
 * @param [in] PetscScalar* xdynload    - load variables
 * @param [in] PetscInt dynlocglob - starting location of load variables in the global vector 
 * @param [in] PetscInt[] V_loc      - global location of VD and VQ V_loc[0] = VD_loc V_loc[1] = VQ_loc
 * @param [in] PetscInt[] I_loc      - global location of ID and IQ I_loc[0] = ID_loc I_loc[1] = IQ_loc
 * Notes: 
 * For each bus, the network equations are ordered as [IQ_loc;ID_loc]
 */
extern PetscErrorCode DYNLoadModelDAERHSJacobian(DYNLoadModel,Mat,PetscReal,PetscScalar,PetscScalar,PetscScalar*,PetscInt,PetscInt[],PetscInt[]);
/**
 * @brief Gets the number and indices of differential and algebraic equations for the load model
 * @param [in] DYNLoadModel dynload - the DYNLoadModel object
 * @param [out] PetscInt* ndiff - number of differential equations
 * @param [out] PetscInt* nalg  - number of algebraic equations
 * @param [out] PetscInt* eqtype - an array of size nvar that has [eqtype[i] = DIFF_EQ if equation i is a differential] or [eqtype[i] = ALG_EQ if equation  i is an algebraic]
 * Notes: 
 * The user does not need to create the array vartype
 */
extern PetscErrorCode DYNLoadModelGetEquationTypes(DYNLoadModel,PetscInt*,PetscInt*,PetscInt*);
/**
 * @brief Returns the model data associated with the particular implementation type
 * @param [in] DYNLoadModel dynload - the DYNLoadModel object
 * @param [out] void** modeldata - the model data struct
 */
extern PetscErrorCode DYNLoadModelGetModelData(DYNLoadModel,void**);
/**
 * @brief Returns the mechanical power at t=t0
 * @param [in] DYNLoadModel dynload - the DYNLoadModel object
 * @param [out] PetscScalar* Pmech0 - Initial mechanical power at t=t0
 */
extern PetscErrorCode DYNLoadModelGetInitialMechanicalPower(DYNLoadModel,PetscScalar*);
/**
 * @brief Gets the location for the first variable for the dynload model from the bus variables array
 * @param [in] DYNLoadModel dynload - the DYNLoadModel object
 * @param [out] PetscInt* loc - the location for the first variable
 */
extern PetscErrorCode DYNLoadModelGetFirstVariableLocation(DYNLoadModel,PetscInt*);
/**
 * @brief Returns the speed deviation and its time derivative, (optional) of the load at the current time.
 * @param [in] DYNLoadModel dynload - the DYNLoadModel object
 * @param [in] PetscReal t      - the current time
 * @param [in] const PetscScalar* xdyn   - the state variables for this load at time t
 * @param [out] PetscScalar* spd     - the load speed deviation
 * @param [out] PetscScalar* dspd_dt - time-derivative of load speed deviation (set to NULL if not required)
 */
extern PetscErrorCode DYNLoadModelGetSpeedDeviation(DYNLoadModel,PetscReal,const PetscScalar*,PetscScalar*,PetscScalar*);
/**
 * @brief Gets the location of the speed deviation variable for the load model relative to the bus variable starting location
 * @param [in] DYNLoadModel dynload - the DYNLoadModel object
 * @param [out] PetscInt* loc - the location of the speed deviation
 * Notes:
 * The location of the speed deviation is relative to the bus variables. So, one should first obtain
 * the first variable location for the load model implementation using DYNLoadModelGetFirstVariableLocation() 
 * and then add the location of the dw variable to return loc. For example, if dw is the 3rd variable for this
 * load model implementation then the following code should be used in the implementation
 * ierr = DYNLoadModelGetFirstVariable(dynload,&firstvarloc);
 * *loc = firstvarloc + 2;
 */
extern PetscErrorCode DYNLoadModelGetSpeedDeviationLocation(DYNLoadModel,PetscInt*);
/**
 * @brief Gets the global location (location in the global vector) of the field voltage for the exciter model
 * @param [in] DYNLoadModel dynload - the DYNLoadModel object
 * @param [out] PetscInt* locglob - the global location of the speed deviation variable
 */
extern PetscErrorCode DYNLoadModelGetSpeedDeviationGlobalLocation(DYNLoadModel,PetscInt*);
/**
 * @brief Returns the frequency of the load at the current time.
 * @param [in] DYNLoadModel dynload - the DYNLoadModel object
 * @param [in] PetscReal t      - the current time
 * @param [in] const PetscScalar* xdyn   - the state variables for this load at time t
 * @param [out] PetscScalar* frequency   - the load frequency in Hz
 */
extern PetscErrorCode DYNLoadModelGetFrequency(DYNLoadModel,PetscReal,const PetscScalar*,PetscScalar*);
/**
 * @brief Set the frequency of the load at the current time.
 * @param [in] DYNLoadModel dynload - the DYNLoadModel object
 * @param [in] PetscReal t      - the current time
 * @param [in] const PetscScalar* xdyn   - the state variables for this load at time t
 * @param [in] PetscScalar* frequency   - the value to be set
 */
extern PetscErrorCode DYNLoadModelSetFrequency(DYNLoadModel,PetscScalar*,PetscScalar);
extern PetscErrorCode DYNLoadModelGetFrequencyLimits(DYNLoadModel,PetscScalar*,PetscScalar*);
/**
 * @brief Returns the partial derivative of the load frequency w.r.t the machine state that governs it
 * @param [in] DYNLoadModel dynload - the DYNLoadModel object
 * @param [in] PetscReal t      - the current time
 * @param [in] const PetscScalar* xdyn   - the state variables for this load at time t
 * @param [out] PetscScalar* dfreqdstate   - dfreq_dstate
 * @param [out] PetscInt* stateloc      - location w.r.t. to the bus variables
 */
extern PetscErrorCode DYNLoadModelGetdFreqdState(DYNLoadModel,PetscReal,const PetscScalar*,PetscScalar*,PetscInt*);
/**
 * @brief Sets the Jacobian of the DAE w.r.t. parameters
 * @param [in] DYNLoadModel dynload - the DYNLoadModel object
 * @param [in] PetscReal t          - the current time
 * @param [in] const PetscScalar* x      - state variables for this bus
 * @param [out] Mat jacP - Matrix of partial derivatives of DYN DAE equations w.r.t. parameters
 * @param [in] PetscInt dynlocglob - starting location of the variables for this dynload model
 * @param [in] PetscInt Valoc  - location of parameter voltage angle (Va)
 * @param [in] PetscInt Vmloc  - location of parameter voltage magnitude (Vm)
 * @param [in] PetscInt Pgloc  - location of parameter load real power output (Pg)
 * @param [in] PetscInt Qgloc  - location of parameter load reactive power output (Qg)
 */
extern PetscErrorCode DYNLoadModelDAERHSJacobianP(DYNLoadModel,PetscReal,const PetscScalar*,Mat,PetscInt,PetscInt,PetscInt);
/**
 * @brief Sets the Jacobian of the DAE w.r.t. parameters
 * @param [in] DYNLoadModel dynload - the DYNLoadModel object
 * @param [in] PetscReal t          - the current time
 * @param [in] const PetscScalar* x      - state variables for this bus
 * @param [out] Vec* jacP - Matrix of partial derivatives of DYN DAE equations w.r.t. parameters
 * @param [in] PetscInt dynlocglob - starting location of the variables for this dynload model
 * @param [in] PetscInt Valoc  - location of parameter voltage angle (Va)
 * @param [in] PetscInt Vmloc  - location of parameter voltage magnitude (Vm)
 * @param [in] PetscInt Pgloc  - location of parameter load real power output (Pg)
 * @param [in] PetscInt Qgloc  - location of parameter load reactive power output (Qg)
 */
extern PetscErrorCode DYNLoadModelDAEFWDRHSJacobianP(DYNLoadModel,PetscReal,const PetscScalar*,Vec*,PetscInt,PetscInt,PetscInt);

extern PetscErrorCode DYNLoadModelSetStatus(DYNLoadModel,PetscInt);

extern PetscErrorCode DYNLoadModelSetConstantPowerLoad(DYNLoadModel,PetscScalar,PetscScalar);

#endif
