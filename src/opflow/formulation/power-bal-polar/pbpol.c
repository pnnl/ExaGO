#include <private/opflowimpl.h>
#include "pbpol.h"

PetscErrorCode OPFLOWFormulationDestroy_PBPOL(OPFLOW opflow)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = PetscFree(opflow->formulation);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSetVariableBounds_PBPOL(OPFLOW opflow,Vec Xl,Vec Xu)
{
  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSetConstraintBounds_PBPOL(OPFLOW opflow,Vec Gl,Vec Gu)
{
  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSetVariableandConstraintBounds_PBPOL(OPFLOW opflow,Vec Xl,Vec Xu, Vec Gl, Vec Gu)
{
  PetscErrorCode ierr;
  PS             ps=opflow->ps;
  PetscScalar    *xl,*xu,*gl,*gu;
  PetscInt       i;
  PSLINE         line;
  PSBUS          bus;
  PetscInt       loc,gloc=0;

  PetscFunctionBegin;

  /* Get array pointers */
  ierr = VecGetArray(Xl,&xl);CHKERRQ(ierr);
  ierr = VecGetArray(Xu,&xu);CHKERRQ(ierr);
  ierr = VecGetArray(Gl,&gl);CHKERRQ(ierr);
  ierr = VecGetArray(Gu,&gu);CHKERRQ(ierr);

  for(i=0; i < ps->nbus; i++) {
    PetscInt k;

    bus = &ps->bus[i];

    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);

    /* Bounds on voltage angles and bounds on real power mismatch equality constraints */
    xl[loc] = PETSC_NINFINITY; xu[loc] = PETSC_INFINITY;
    gl[gloc] = 0.0;   gu[gloc] = 0.0;

    /* Bounds on voltage magnitudes and bounds on reactive power mismatch equality constraints */
    xl[loc+1] = bus->Vmin; xu[loc+1] = bus->Vmax;
    gl[gloc+1] = 0.0;       gu[gloc+1] = 0.0;

    if(bus->ide == REF_BUS || bus->ide == ISOLATED_BUS) xl[loc] = xu[loc] = bus->va*PETSC_PI/180.0;
    if(bus->ide == ISOLATED_BUS) xl[loc+1] = xu[loc+1] = bus->vm;
    
    for(k=0; k < bus->ngen; k++) {
      PSGEN gen;

      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      loc = loc+2;
      if(!gen->status) xl[loc] = xu[loc] = xl[loc+1] = xu[loc+1] = 0.0;
      else {
	xl[loc] = gen->pb; /* PGmin */
	xu[loc] = gen->pt; /* PGmax */
	xl[loc+1] = gen->qb; /* QGmin */
	xu[loc+1] = gen->qt; /* QGmax */
	/* pb, pt, qb, qt are converted in p.u. in ps.c */
      }
    }
    gloc += 2;
  }
  
  for(i=0; i < ps->Nbranch; i++) {
    line = &ps->line[i];

    /* Line flow inequality constraints */
    if(!line->status) gl[gloc] = gu[gloc] = gl[gloc+1] = gu[gloc+1] = 0.0;
    else {
      gl[gloc] = gl[gloc+1] = 0.0; 
      gu[gloc] = gu[gloc+1] = (line->rateA/ps->MVAbase)*(line->rateA/ps->MVAbase);
    }    
    gloc += 2;
  }

  ierr = VecRestoreArray(Xl,&xl);CHKERRQ(ierr);
  ierr = VecRestoreArray(Xu,&xu);CHKERRQ(ierr);
  ierr = VecRestoreArray(Gl,&gl);CHKERRQ(ierr);
  ierr = VecRestoreArray(Gu,&gu);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSetInitialGuess_PBPOL(OPFLOW opflow,Vec X)
{
  PetscErrorCode ierr;
  PS             ps=opflow->ps;
  const PetscScalar    *xl,*xu;
  PetscScalar    *x;
  PetscInt       i;
  PSBUS          bus;
  PetscInt       loc;

  PetscFunctionBegin;

  /* Get array pointers */
  ierr = VecGetArray(X,&x);CHKERRQ(ierr);
  ierr = VecGetArrayRead(opflow->Xl,&xl);CHKERRQ(ierr);
  ierr = VecGetArrayRead(opflow->Xu,&xu);CHKERRQ(ierr);
  
  for(i=0; i < ps->nbus; i++) {
    PetscInt k;

    bus = &ps->bus[i];

    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);

    /* Initial guess for voltage angles and bounds on voltage magnitudes */
    x[loc]   = (xl[loc] + xu[loc])/2.0;
    x[loc+1] = (xl[loc+1] + xu[loc+1])/2.0;

    for(k=0; k < bus->ngen; k++) {
      PSGEN gen;

      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      loc = loc+2;

      x[loc]   = (xl[loc] + xu[loc])/2.0;   /* Initial guess for Pg */
      x[loc+1] = (xl[loc+1] + xu[loc+1])/2.0; /* Initial guess for Qg */
    }
  }

  ierr = VecRestoreArray(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(opflow->Xl,&xl);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(opflow->Xu,&xu);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeEqualityConstraints_PBPOL(OPFLOW opflow,Vec X,Vec Ge)
{
  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeInequalityConstraints_PBPOL(OPFLOW opflow,Vec X,Vec Gi)
{
  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeConstraints_PBPOL(OPFLOW opflow,Vec X,Vec G)
{
  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeObjandGradient_PBPOL(OPFLOW opflow,Vec X,PetscScalar *obj,Vec Grad)
{
  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeObjective_PBPOL(OPFLOW opflow,Vec X,PetscScalar *obj)
{
  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeGradient_PBPOL(OPFLOW opflow,Vec X,Vec Grad)
{
  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWFormulationSetNumVariables_PBPOL(OPFLOW opflow,PetscInt *busnvar,PetscInt *branchnvar,PetscInt *nx)
{
  PetscInt i,ngen;
  PS       ps=opflow->ps;
  PSBUS    bus;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  
  *nx = 0;
  /* No variables for the branches */
  for(i=0; i < ps->nbranch; i++) {
    branchnvar[i] = 0;
    *nx += branchnvar[i];
  }

  /* Variables for the buses */
  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSGetNGen(bus,&ngen);CHKERRQ(ierr);
    /* Number of variables = 2 + 2*ngen (2 variables for voltage + Pg, Qg for each gen) */
    busnvar[i] = 2 + 2*ngen;
    *nx += busnvar[i];
  }

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWFormulationSetNumConstraints_PBPOL(OPFLOW opflow,PetscInt *nconeq,PetscInt *nconineq)
{
  PS  ps = opflow->ps;

  PetscFunctionBegin;
  *nconeq = 2*ps->nbus;
  *nconineq = 2*ps->nbranch;

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeEqualityConstraintJacobian(OPFLOW opflow,Vec X,Mat Jac_Ge)
{
  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeInequalityConstraintJacobian(OPFLOW opflow,Vec X,Mat Jac_Gi)
{
  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeHessian(OPFLOW opflow,Vec X,Vec Lambda,Mat H)
{
  PetscFunctionBegin;

  PetscFunctionReturn(0);
}


PetscErrorCode OPFLOWFormulationCreate_PBPOL(OPFLOW opflow)
{
  PBPOL pbpol;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  
  ierr = PetscCalloc1(1,&pbpol);CHKERRQ(ierr);

  opflow->formulation = pbpol;

  /* Inherit Ops */
  opflow->formops.destroy = OPFLOWFormulationDestroy_PBPOL;
  opflow->formops.setnumvariables = OPFLOWFormulationSetNumVariables_PBPOL;
  opflow->formops.setnumconstraints = OPFLOWFormulationSetNumConstraints_PBPOL;
  opflow->formops.setvariablebounds = OPFLOWSetVariableBounds_PBPOL;
  opflow->formops.setconstraintbounds = OPFLOWSetConstraintBounds_PBPOL;
  opflow->formops.setvariableandconstraintbounds = OPFLOWSetVariableandConstraintBounds_PBPOL;
  opflow->formops.setinitialguess = OPFLOWSetInitialGuess_PBPOL;
  opflow->formops.computeequalityconstraints = OPFLOWComputeEqualityConstraints_PBPOL;
  opflow->formops.computeinequalityconstraints = OPFLOWComputeInequalityConstraints_PBPOL;
  opflow->formops.computeequalityconstraintjacobian = OPFLOWComputeEqualityConstraintJacobian;
  opflow->formops.computeinequalityconstraintjacobian = OPFLOWComputeInequalityConstraintJacobian;
  opflow->formops.computehessian = OPFLOWComputeHessian;
  opflow->formops.computeobjandgradient = OPFLOWComputeObjandGradient_PBPOL;
  opflow->formops.computeobjective = OPFLOWComputeObjective_PBPOL;
  opflow->formops.computegradient  = OPFLOWComputeGradient_PBPOL;
  
  PetscFunctionReturn(0);
}
