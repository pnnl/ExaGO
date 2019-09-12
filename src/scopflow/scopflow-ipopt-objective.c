#include <private/psimpl.h>
#include <private/opflowimpl.h>
#include <private/scopflowipoptimpl.h>
#include <math.h>

/* Implements functions for the objective (objective, gradient, and Hessian)
*/

/*
  SCOPFLOWObjectiveFunction - The objective function for the security constrained optimal power flow

  Input Parameters:
+ opflow - the OPFLOW object
. row    - scenario number
- X      - the current iterate

  Output Parameters:
. obj - the objective function value (scalar)
*/
PetscErrorCode SCOPFLOWObjectiveFunction(SCOPFLOW scopflow,PetscInt row,Vec X, PetscScalar* obj)
{
  PetscErrorCode ierr;
  const PetscScalar *x;
  OPFLOW         opflow=scopflow->opflows[row];
  PS             ps=opflow->ps;
  PetscInt       i;
  PSBUS          bus;
  PetscInt       loc;

  PetscFunctionBegin;

  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);

  *obj = 0.0;
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
  SCOPFLOWObjGradientFunction - The gradient of the objective function for the 
                                security constrained optimal power flow

  Input Parameters:
+ opflow - the OPFLOW object
. row    - scenario number
- X      - the current iterate

  Output Parameters:
. grad - the objective function gradient
*/
PetscErrorCode SCOPFLOWObjGradientFunction(SCOPFLOW scopflow,PetscInt row,Vec X, Vec grad)
{
  PetscErrorCode ierr;
  const PetscScalar *x;
  PetscScalar    *df;
  OPFLOW         opflow=scopflow->opflows[row];
  PS             ps=opflow->ps;
  PetscInt       i;
  PSBUS          bus;
  PetscInt       loc;

  PetscFunctionBegin;
  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);
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

/*
  SCOPFLOWComputeObjectiveHessian - Computes the Hessian for the objective function part
  
  Input Parameters:
+ scopflow - the SCOPFLOW object
. scenario - the scenario number
- X        - solution vecto X

  Output Parameters:
. H - the Hessian part for the objective function

*/
PetscErrorCode SCOPFLOWComputeObjectiveHessian(SCOPFLOW scopflow,PetscInt scenario,Vec X,Mat H) 
{
  PetscErrorCode ierr;
  OPFLOW         opflow=scopflow->opflows[scenario];
  PS             ps=opflow->ps;
  PetscInt       i,k;
  PSBUS          bus;
  PSGEN          gen;
  PetscInt       xloc;
  const PetscScalar *x;
  PetscInt       row[2],col[2];
  PetscScalar    val[2];
  PetscScalar    obj_factor = scopflow->obj_factor;


  PetscFunctionBegin;

  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);

  // for the part of objective
  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    ierr = PSBUSGetVariableLocation(bus,&xloc);CHKERRQ(ierr);
   
    for(k=0; k < bus->ngen; k++) {
      xloc = xloc+2;
      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      if(!gen->status) continue;
      row[0] = xloc;
      col[0] = xloc;
      val[0] = obj_factor*2.0*gen->cost_alpha*ps->MVAbase*ps->MVAbase;
      ierr = MatSetValues(H,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);
    }
  }

  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


int str_eval_f(double* x0, double* x1, double* obj, CallBackDataPtr cbd) 
{
  PetscErrorCode ierr;
  int row = cbd->row_node_id;
  SCOPFLOW scopflow=(SCOPFLOW)cbd->prob;
  OPFLOW   opflow=scopflow->opflows[row];
  
  *obj = 0.0;
  if(row == 0) {
    ierr = VecPlaceArray(opflow->X,x0);CHKERRQ(ierr);
    ierr = SCOPFLOWObjectiveFunction(scopflow,row,opflow->X,obj);CHKERRQ(ierr);
    ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
  } else {
    if(!scopflow->first_stage_gen_cost_only) {
      ierr = VecPlaceArray(opflow->X,x1);CHKERRQ(ierr);
      ierr = SCOPFLOWObjectiveFunction(scopflow,row,opflow->X,obj);CHKERRQ(ierr);
      ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
    }
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
    if(row == 0) {
      x = x0;
      ierr = VecPlaceArray(opflow->X,x);CHKERRQ(ierr);
      ierr = VecPlaceArray(opflow->gradobj,grad);CHKERRQ(ierr);
      ierr = VecSet(opflow->gradobj,0.0);CHKERRQ(ierr);
      ierr = SCOPFLOWObjGradientFunction(scopflow,row,opflow->X,opflow->gradobj);CHKERRQ(ierr);
      ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
      ierr = VecResetArray(opflow->gradobj);CHKERRQ(ierr);
    } else {
      if(scopflow->first_stage_gen_cost_only) {
	ierr = VecPlaceArray(opflow->gradobj,grad);CHKERRQ(ierr);
	ierr = VecSet(opflow->gradobj,0.0);CHKERRQ(ierr);
	ierr = VecResetArray(opflow->gradobj);CHKERRQ(ierr);
	return 1;
      }
      x = x1;
      ierr = VecPlaceArray(opflow->X,x);CHKERRQ(ierr);
      ierr = VecPlaceArray(opflow->gradobj,grad);CHKERRQ(ierr);
      ierr = VecSet(opflow->gradobj,0.0);CHKERRQ(ierr);
      ierr = SCOPFLOWObjGradientFunction(scopflow,row,opflow->X,opflow->gradobj);CHKERRQ(ierr);
      ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
      ierr = VecResetArray(opflow->gradobj);CHKERRQ(ierr);
    }
  }
    
  return 1;
}
