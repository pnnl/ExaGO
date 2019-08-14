#include <private/psimpl.h>
#include <private/opflowimpl.h>
#include <private/scopflowpipsimpl.h>
#include <math.h>

/* Implements functions for the objective function and gradient of
   objective function 
*/

/*
  OPFLOWObjectiveFunction - The objective function for the optimal power flow

  Input Parameters:
+ opflow - the OPFLOW object
. X      - the current iterate

  Output Parameters:
. obj - the objective function value (scalar)
*/
PetscErrorCode OPFLOWObjectiveFunction(OPFLOW opflow,Vec X, PetscScalar* obj)
{
  PetscErrorCode ierr;
  const PetscScalar *x;
  PS             ps=opflow->ps;
  PetscInt       i;
  PSBUS          bus;
  PetscInt       loc;

  PetscFunctionBegin;
  *obj = 0.0;
  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);

  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
    
    PetscInt k;
    PSGEN    gen;
    PetscScalar Pg;
    for(k=0; k < bus->ngen; k++) {
      loc = loc+2;
      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      if(!gen->status) continue;
      Pg = x[loc]*ps->MVAbase;
      *obj += gen->cost_alpha*Pg*Pg + gen->cost_beta*Pg + gen->cost_gamma;
    }
  }
  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


/*
  OPFLOWObjGradientFunction - The gradient of the objective function for the optimal power flow

  Input Parameters:
+ opflow - the OPFLOW object
. X      - the current iterate

  Output Parameters:
. grad - the objective function gradient
*/
PetscErrorCode OPFLOWObjGradientFunction(OPFLOW opflow,Vec X, Vec grad)
{
  PetscErrorCode ierr;
  const PetscScalar *x;
  PetscScalar    *df;
  PS             ps=opflow->ps;
  PetscInt       i;
  PSBUS          bus;
  PetscInt       loc;

  PetscFunctionBegin;
  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecSet(grad,0.0);CHKERRQ(ierr);
  ierr = VecGetArray(grad,&df);CHKERRQ(ierr);

  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
    
    PetscInt k;
    PSGEN    gen;
    PetscScalar Pg;
    for(k=0; k < bus->ngen; k++) {
      loc = loc+2;
      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      if(!gen->status) continue;
      Pg = x[loc]*ps->MVAbase;
      df[loc] = ps->MVAbase*(2*gen->cost_alpha*Pg + gen->cost_beta);
    }
  }
  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(grad,&df);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

int str_eval_f(double* x0, double* x1, double* obj, CallBackDataPtr cbd) 
{
  PetscErrorCode ierr;
  int row = cbd->row_node_id;
  SCOPFLOW scopflow=(SCOPFLOW)cbd->prob;
  OPFLOW   opflow=scopflow->opflows[row];
  
  if(row == 0 ) {
    ierr = VecPlaceArray(opflow->X,x0);CHKERRQ(ierr);
    ierr = OPFLOWObjectiveFunction(opflow,opflow->X,obj);CHKERRQ(ierr);
    ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
  } else {
    ierr = VecPlaceArray(opflow->X,x1);CHKERRQ(ierr);
    ierr = OPFLOWObjectiveFunction(opflow,opflow->X,obj);CHKERRQ(ierr);
    ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
  }

  return 1;
}

int str_eval_grad_f(double* x0, double* x1, double* grad, CallBackDataPtr cbd) 
{
  PetscErrorCode ierr;
  int row = cbd->row_node_id;
  int col = cbd->col_node_id;
  SCOPFLOW scopflow=(SCOPFLOW)cbd->prob;
  OPFLOW   opflow=scopflow->opflows[row];
  double   *x;
  
  if(row == col) {
    if(row == 0) x = x0;
    else x = x1;
      
    ierr = VecPlaceArray(opflow->X,x);CHKERRQ(ierr);
    ierr = VecPlaceArray(opflow->gradobj,grad);CHKERRQ(ierr);
    ierr = OPFLOWObjGradientFunction(opflow,opflow->X,opflow->gradobj);CHKERRQ(ierr);
    ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
    ierr = VecResetArray(opflow->gradobj);CHKERRQ(ierr);
  }
    
  return 1;
}
