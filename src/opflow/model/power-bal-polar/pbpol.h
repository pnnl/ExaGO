#include <opflow.h>

typedef struct _p_FormPBPOL *PBPOL;

struct _p_FormPBPOL
{
  PetscInt Nx;
  PetscErrorCode ModelDestroy(OPFLOW opflow);
  PetscErrorCode SetVariableBounds(OPFLOW opflow,Vec Xl,Vec Xu);
  PetscErrorCode SetConstraintBounds(OPFLOW opflow,Vec Gl,Vec Gu);
  PetscErrorCode SetVariableandConstraintBounds(OPFLOW opflow,Vec Xl,Vec Xu, Vec Gl, Vec Gu);
  PetscErrorCode SetInitialGuess(OPFLOW opflow,Vec X);
  PetscErrorCode ComputeEqualityConstraints(OPFLOW opflow,Vec X,Vec Ge);
  PetscErrorCode ComputeEqualityConstraintJacobian(OPFLOW opflow,Vec X,Mat Je);
  PetscErrorCode ComputeInequalityConstraints(OPFLOW opflow,Vec X,Vec Gi);
  PetscErrorCode ComputeInequalityConstraintJacobian(OPFLOW opflow,Vec X,Mat Ji);
  PetscErrorCode ComputeConstraints(OPFLOW opflow,Vec X,Vec G);
  PetscErrorCode ComputeObjective(OPFLOW opflow,Vec X,PetscScalar *obj);
  PetscErrorCode ComputeGradient(OPFLOW opflow,Vec X,Vec grad);
  PetscErrorCode ComputeObjandGradient(OPFLOW opflow,Vec X,PetscScalar *obj,Vec Grad);
  PetscErrorCode ModelSetNumVariables(OPFLOW opflow,PetscInt *busnvar,PetscInt *branchnvar,PetscInt *nx);
  PetscErrorCode ModelSetNumConstraints(OPFLOW opflow,PetscInt *branchnconeq,PetscInt *busnconeq,PetscInt *nconeq,PetscInt *nconineq);
  PetscErrorCode ComputeEqualityConstraintsHessian(OPFLOW opflow,Vec X,Vec Lambda,Mat H);
  PetscErrorCode ComputeInequalityConstraintsHessian(OPFLOW opflow, Vec X, Vec Lambda,Mat H);
  PetscErrorCode ComputeObjectiveHessian(OPFLOW opflow,Vec X,Mat H);
  PetscErrorCode ComputeHessian(OPFLOW opflow,Vec X,Vec Lambdae,Vec Lambdai,Mat H);
  PetscErrorCode SolutionToPS(OPFLOW opflow);
  PetscErrorCode ModelSetUp(OPFLOW opflow);
  PetscErrorCode ModelCreate(OPFLOW opflow);
};

/** Function declarations for use with legacy OPFLOW during c -> c++ transition */
extern PetscErrorCode OPFLOWModelDestroy_PBPOL(OPFLOW opflow);
extern PetscErrorCode OPFLOWSetVariableBounds_PBPOL(OPFLOW opflow,Vec Xl,Vec Xu);
extern PetscErrorCode OPFLOWSetConstraintBounds_PBPOL(OPFLOW opflow,Vec Gl,Vec Gu);
extern PetscErrorCode OPFLOWSetVariableandConstraintBounds_PBPOL(OPFLOW opflow,Vec Xl,Vec Xu, Vec Gl, Vec Gu);
extern PetscErrorCode OPFLOWSetInitialGuess_PBPOL(OPFLOW opflow,Vec X);
extern PetscErrorCode OPFLOWComputeEqualityConstraints_PBPOL(OPFLOW opflow,Vec X,Vec Ge);
extern PetscErrorCode OPFLOWComputeEqualityConstraintJacobian_PBPOL(OPFLOW opflow,Vec X,Mat Je);
extern PetscErrorCode OPFLOWComputeInequalityConstraints_PBPOL(OPFLOW opflow,Vec X,Vec Gi);
extern PetscErrorCode OPFLOWComputeInequalityConstraintJacobian_PBPOL(OPFLOW opflow,Vec X,Mat Ji);
extern PetscErrorCode OPFLOWComputeConstraints_PBPOL(OPFLOW opflow,Vec X,Vec G);
extern PetscErrorCode OPFLOWComputeObjective_PBPOL(OPFLOW opflow,Vec X,PetscScalar *obj);
extern PetscErrorCode OPFLOWComputeGradient_PBPOL(OPFLOW opflow,Vec X,Vec grad);
extern PetscErrorCode OPFLOWComputeObjandGradient_PBPOL(OPFLOW opflow,Vec X,PetscScalar *obj,Vec Grad);
extern PetscErrorCode OPFLOWModelSetNumVariables_PBPOL(OPFLOW opflow,PetscInt *busnvar,PetscInt *branchnvar,PetscInt *nx);
extern PetscErrorCode OPFLOWModelSetNumConstraints_PBPOL(OPFLOW opflow,PetscInt *branchnconeq,PetscInt *busnconeq,PetscInt *nconeq,PetscInt *nconineq);
extern PetscErrorCode OPFLOWComputeEqualityConstraintsHessian_PBPOL(OPFLOW opflow,Vec X,Vec Lambda,Mat H);
extern PetscErrorCode OPFLOWComputeInequalityConstraintsHessian_PBPOL(OPFLOW opflow, Vec X, Vec Lambda,Mat H);
extern PetscErrorCode OPFLOWComputeObjectiveHessian_PBPOL(OPFLOW opflow,Vec X,Mat H);
extern PetscErrorCode OPFLOWComputeHessian_PBPOL(OPFLOW opflow,Vec X,Vec Lambdae,Vec Lambdai,Mat H);
extern PetscErrorCode OPFLOWSolutionToPS_PBPOL(OPFLOW opflow);
extern PetscErrorCode OPFLOWModelSetUp_PBPOL(OPFLOW opflow);
extern PetscErrorCode OPFLOWModelCreate_PBPOL(OPFLOW opflow);
