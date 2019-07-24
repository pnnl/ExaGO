#include <private/opflowimpl.h>
#include <../src/tao/constrained/impls/ipm/ipm.h> /*I "ipm.h" I*/

/****************************************************
  Designed to be a side file dedicated to
  the functions dealing with the Objective function
  Order to file:
    OPFLOWObjectiveandGradientFunction
    OPFLOWCreateObjectiveHessian
    OPFLOWHessian
*****************************************************/

/*
  OPFLOWObjectiveandGradientFunction - The objective and gradient function for the optimal power flow

  Input Parameters:
+ nlp - the TAO object
. X      - the current iterate

  Output Parameters:
+ obj - the objective function value (scalar)
- grad - the gradient vector

  Objective function:
    f(x) = SUM(a*(Pg*MVAbase)^2 + b*(Pg*MVAbase) + c); for each xi
  Gradient of Objective:
    grad = f_x st. grad[i] = 2*a*(Pg*MVAbase)*MVAbase + b*MVAbase
*/
PetscErrorCode OPFLOWObjectiveandGradientFunction(Tao nlp,Vec X,PetscScalar* obj,Vec grad,void* ctx)
{
  PetscErrorCode ierr;
  PetscInt       i,k;
  PetscInt       loc;
  PetscScalar    *df,Pg,sobj;
  OPFLOW         opflow=(OPFLOW)ctx;
  PS             ps=opflow->ps;
  Vec            localgrad;
  Vec            localX=opflow->localX;
  PSBUS          bus;
  PSGEN          gen;
  const PetscScalar *x;

  PetscFunctionBegin;
  *obj = 0.0;

  ierr = VecDuplicate(localX,&localgrad);CHKERRQ(ierr);

  ierr = DMGlobalToLocalBegin(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = VecSet(localgrad,0.0);CHKERRQ(ierr);

  ierr = VecGetArrayRead(localX,&x);CHKERRQ(ierr);
  ierr = VecGetArray(localgrad,&df);CHKERRQ(ierr);

  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);

    for(k=0; k < bus->ngen; k++) {
      loc = loc+2;
      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      if(!gen->status) continue;
      Pg = x[loc]*ps->MVAbase;
      sobj += gen->cost_alpha*Pg*Pg + gen->cost_beta*Pg + gen->cost_gamma;
      df[loc] = ps->MVAbase*(2*gen->cost_alpha*Pg + gen->cost_beta);
    }
  }
  ierr = VecRestoreArrayRead(localX,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(localgrad,&df);CHKERRQ(ierr);

  ierr = DMLocalToGlobalBegin(ps->networkdm,localgrad,INSERT_VALUES,grad);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(ps->networkdm,localgrad,INSERT_VALUES,grad);CHKERRQ(ierr);

  ierr = VecDestroy(&localgrad);CHKERRQ(ierr);
  ierr = MPI_Allreduce(&sobj,obj,1,MPI_DOUBLE,MPI_SUM,opflow->comm->type);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  OPFLOWCreateObjectiveHessian - Creates the constant Hessian matrix for the objective function f(x)

  Input Parameters
. the optimal power flow application object

  Output Parameters
  N/A - sets directly into tao->hessian

  Objective function:
    f(x) = SUM(a*(Pg*MVAbase)^2 + b*(Pg*MVAbase) + c); for each xi
  Hessian function:
    H(x) = f_xx = constant diagonal matrix: 2*a*MVAbase*MVAbase;
 */
PetscErrorCode OPFLOWCreateObjectiveHessian(OPFLOW opflow)
{
  PetscErrorCode ierr;
  PetscInt       i,k,n,loc;
  PetscScalar    val;
  PS             ps=opflow->ps;
  PSBUS          bus;
  PSGEN          gen;
  Tao            tao=opflow->nlp;
  Mat            H=tao->hessian;

  PetscFunctionBegin;
  ierr = VecGetLocalSize(opflow->X,&n);CHKERRQ(ierr);
  ierr = MatSetSizes(H,n,n,PETSC_DECIDE,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = MatSetType(H,MATAIJ);CHKERRQ(ierr);
  ierr = MatSeqAIJSetPreallocation(H,1,NULL);CHKERRQ(ierr);
  ierr = MatMPIAIJSetPreallocation(H,1,NULL,0,NULL);CHKERRQ(ierr);

  for (i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSGetVariableGlobalLocation(bus,&loc);CHKERRQ(ierr);
    for (k=0; k < bus->ngen; k++) {
      loc = loc + 2;
      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      if (!gen->status) continue;
      val = 2*ps->MVAbase*ps->MVAbase*gen->cost_alpha;
      ierr = MatSetValue(H,loc,loc,val,INSERT_VALUES);CHKERRQ(ierr);
    }
  }
  ierr = MatAssemblyBegin(H,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(H,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  //ierr = MatView(H,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  OPFLOWHessian - Tao resets hessian each itteration,
      this function prevents wasted recalculations

  Input Parameters
  tao - Tao solver object
  X   - the current iterate
  ctx - application data set

  Output Parameters (both same as entered)
  H - Hessain matrix
  H_pre - Preconditioner for Hessian
 */
PetscErrorCode OPFLOWHessian(Tao tao, Vec X, Mat H, Mat H_pre, void* ctx)
{
  PetscFunctionBegin;
  H     = tao->hessian;
  H_pre = tao->hessian;
  PetscFunctionReturn(0);
}
