#include <exago_config.h>
#if defined(EXAGO_ENABLE_IPOPT)

#include <private/opflowimpl.h>
#include <private/tcopflowimpl.h>
#include <private/scopflowimpl.h>
#include "scopflow-ipopt.h"

/* IPOPT callback functions */
Bool eval_scopflow_f(PetscInt n, PetscScalar* x, Bool new_x,
            PetscScalar* obj_value, UserDataPtr user_data)
{
  PetscErrorCode ierr;
  SCOPFLOW  scopflow=(SCOPFLOW)user_data;

  *obj_value = 0.0;

  ierr = VecPlaceArray(scopflow->X,x);CHKERRQ(ierr);
  ierr = (*scopflow->modelops.computeobjective)(scopflow,scopflow->X,obj_value);CHKERRQ(ierr);
  ierr = VecResetArray(scopflow->X);CHKERRQ(ierr);
  
  return TRUE;
}

Bool eval_scopflow_grad_f(PetscInt n, PetscScalar* x, Bool new_x,
                 PetscScalar* grad_f, UserDataPtr user_data)
{

  PetscErrorCode ierr;
  SCOPFLOW  scopflow=(SCOPFLOW)user_data;

  ierr = VecPlaceArray(scopflow->X,x);CHKERRQ(ierr);
  ierr = VecPlaceArray(scopflow->gradobj,grad_f);CHKERRQ(ierr);
  ierr = (*scopflow->modelops.computegradient)(scopflow,scopflow->X,scopflow->gradobj);CHKERRQ(ierr);
  ierr = VecResetArray(scopflow->X);CHKERRQ(ierr);
  ierr = VecResetArray(scopflow->gradobj);CHKERRQ(ierr);
  
  return TRUE;
}

Bool eval_scopflow_g(PetscInt n, PetscScalar* x, Bool new_x,
            PetscInt m, PetscScalar* g, UserDataPtr user_data)
{
  PetscErrorCode ierr;
  SCOPFLOW  scopflow=(SCOPFLOW)user_data;

  ierr = VecPlaceArray(scopflow->X,x);CHKERRQ(ierr);
  ierr = VecPlaceArray(scopflow->G,g);CHKERRQ(ierr);
  ierr = (*scopflow->modelops.computeconstraints)(scopflow,scopflow->X,scopflow->G);CHKERRQ(ierr);
  ierr = VecResetArray(scopflow->X);CHKERRQ(ierr);
  ierr = VecResetArray(scopflow->G);CHKERRQ(ierr);

  return TRUE;
}

Bool eval_scopflow_jac_g(PetscInt n, PetscScalar *x, Bool new_x,
                PetscInt m, PetscInt nele_jac,
                PetscInt *iRow, PetscInt *jCol, PetscScalar *values,
                UserDataPtr user_data)
{

  PetscErrorCode ierr;
  SCOPFLOW       scopflow=(SCOPFLOW)user_data;
  PetscInt       *iRowstart = iRow,*jColstart=jCol;
  PetscInt       roffset,coffset;
  PetscInt       nrow,ncol;
  PetscInt       nvals;
  const PetscInt *cols;
  const PetscScalar *vals;
  PetscInt i,j;

  if(values == NULL) {
    /* Set locations only */

    roffset = 0;
    coffset = 0;

    /* Compute Jacobian */
    ierr = (*scopflow->modelops.computejacobian)(scopflow,scopflow->X,scopflow->Jac);CHKERRQ(ierr);

    ierr = MatGetSize(scopflow->Jac,&nrow,&ncol);CHKERRQ(ierr);
    
    /* Copy over locations to triplet format */
    for(i=0; i < nrow; i++) {
      ierr = MatGetRow(scopflow->Jac,i,&nvals,&cols,&vals);CHKERRQ(ierr);
      for(j=0; j < nvals; j++) {
	iRowstart[j] = roffset + i;
	jColstart[j] = coffset + cols[j];
      }
      /* Increment iRow,jCol pointers */
      iRowstart += nvals;
      jColstart += nvals;
      ierr = MatRestoreRow(scopflow->Jac,i,&nvals,&cols,&vals);CHKERRQ(ierr);
    }
  } else {
    ierr = VecPlaceArray(scopflow->X,x);CHKERRQ(ierr);
    /* Compute jacobian */
    ierr = (*scopflow->modelops.computejacobian)(scopflow,scopflow->X,scopflow->Jac);CHKERRQ(ierr);

    ierr = MatGetSize(scopflow->Jac,&nrow,&ncol);CHKERRQ(ierr);
    /* Copy over values */
    for(i=0; i < nrow; i++) {
      ierr = MatGetRow(scopflow->Jac,i,&nvals,&cols,&vals);CHKERRQ(ierr);
      for(j=0; j < nvals; j++) {
	values[j] = vals[j];
      }
      values += nvals;
      ierr = MatRestoreRow(scopflow->Jac,i,&nvals,&cols,&vals);CHKERRQ(ierr);
    }
    ierr = VecResetArray(scopflow->X);CHKERRQ(ierr);
  }

  return TRUE;
}


Bool eval_scopflow_h(PetscInt n, PetscScalar *x, Bool new_x, PetscScalar obj_factor,
            PetscInt m, PetscScalar *lambda, Bool new_lambda,
            PetscInt nele_hess, PetscInt *iRow, PetscInt *jCol,
            PetscScalar *values, UserDataPtr user_data)
{
  PetscErrorCode ierr;
  PetscInt       nrow;
  SCOPFLOW         scopflow=(SCOPFLOW)user_data;
  PetscInt       j,k;
  const PetscInt *cols;
  const PetscScalar *vals;
  PetscInt ctr=0;
  PetscInt nvals;

  scopflow->obj_factor = obj_factor;

  ierr = MatGetSize(scopflow->Hes,&nrow,NULL);CHKERRQ(ierr);

  if(values == NULL) {
    /* Set locations only */

    /* Compute Hessian */
    ierr = (*scopflow->modelops.computehessian)(scopflow,scopflow->X,scopflow->Lambda,scopflow->Hes);CHKERRQ(ierr);

    /* Copy over locations to triplet format */
    /* Note that IPOPT
       requires a lower diagonal Hessian (see note https://coin-or.github.io/Ipopt/INTERFACES.html#INTERFACE_CODE)
       Hence, we only add lower diagonal locations
    */
    for(j=0; j < nrow; j++) {
      ierr = MatGetRow(scopflow->Hes,j,&nvals,&cols,&vals);CHKERRQ(ierr);
      ctr = 0;
      for(k=0; k < nvals; k++) {
	if(cols[k] >= j) { /* upper triangle */
	  /* save as lower triangle locations */
	    iRow[ctr] = cols[k];
	    jCol[ctr] = j;
	    ctr++;
	}
      }
      iRow += ctr;
      jCol += ctr;
      ierr = MatRestoreRow(scopflow->Hes,j,&nvals,&cols,&vals);CHKERRQ(ierr);
    }
  } else {
    /* Copy values */
    ierr = VecPlaceArray(scopflow->X,x);CHKERRQ(ierr);
    ierr = VecPlaceArray(scopflow->Lambda,lambda);CHKERRQ(ierr);

    /* Compute Hessian */
    ierr = (*scopflow->modelops.computehessian)(scopflow,scopflow->X,scopflow->Lambda,scopflow->Hes);CHKERRQ(ierr);

    /* Copy over values to triplet format */
    /* Note that IPOPT
       requires a lower diagonal Hessian (see note https://coin-or.github.io/Ipopt/INTERFACES.html#INTERFACE_CODE)
       Hence, we only add lower diagonal locations
    */
    for(j=0; j < nrow; j++) {
      ierr = MatGetRow(scopflow->Hes,j,&nvals,&cols,&vals);CHKERRQ(ierr);
      ctr = 0;
      for(k=0; k < nvals; k++) {
	if(cols[k] >= j) { /* upper triangle */
	  /* save as lower triangle locations */
	  values[ctr] = vals[k];
	  ctr++;
	}
      }
      values += ctr;
      ierr = MatRestoreRow(scopflow->Hes,j,&nvals,&cols,&vals);CHKERRQ(ierr);
    }

    ierr = VecResetArray(scopflow->X);CHKERRQ(ierr);
    ierr = VecResetArray(scopflow->Lambda);CHKERRQ(ierr);
  }

  return TRUE;
}

Bool SCOPFLOWSolverMonitor_IPOPT(Index alg_mod,Index iter_count,Number obj_value,Number inf_pr,
			       Number inf_du,Number mu,Number d_norm,Number regularization_size,
			       Number alpha_du,Number alpha_pr,Index ls_trials,
                               UserDataPtr user_data)
{
  SCOPFLOW  scopflow=(SCOPFLOW)user_data;
  scopflow->numiter = iter_count;
  return 1;
}

PetscErrorCode SCOPFLOWSolverSolve_IPOPT(SCOPFLOW scopflow)
{
  PetscErrorCode     ierr;
  SCOPFLOWSolver_IPOPT scopflowipopt = (SCOPFLOWSolver_IPOPT)scopflow->solver;
  MatInfo            info_jac,info_hes;
  PetscScalar        *x,*g,*xl,*xu,*gl,*gu,*lam;

  PetscFunctionBegin;

  scopflowipopt->nnz_jac_g = scopflowipopt->nnz_hes = 0;

  /* Compute nonzeros for the Jacobian */
  ierr = (*scopflow->modelops.computejacobian)(scopflow,scopflow->X,scopflow->Jac);CHKERRQ(ierr);
  ierr = MatSetOption(scopflow->Jac,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE);CHKERRQ(ierr);

  ierr = MatGetInfo(scopflow->Jac,MAT_LOCAL,&info_jac);CHKERRQ(ierr);
  scopflowipopt->nnz_jac_g = info_jac.nz_used;

  /* Compute non-zeros for Hessian */
  ierr = (*scopflow->modelops.computehessian)(scopflow,scopflow->X,scopflow->Lambda,scopflow->Hes);CHKERRQ(ierr);
  ierr = MatSetOption(scopflow->Hes,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE);CHKERRQ(ierr);

  ierr  = MatGetInfo(scopflow->Hes,MAT_LOCAL,&info_hes);CHKERRQ(ierr);
  scopflowipopt->nnz_hes = (info_hes.nz_used  - scopflow->Nx)/2 + scopflow->Nx;

  ierr = VecGetArray(scopflow->Xl,&xl);CHKERRQ(ierr);
  ierr = VecGetArray(scopflow->Xu,&xu);CHKERRQ(ierr);
  ierr = VecGetArray(scopflow->Gl,&gl);CHKERRQ(ierr);
  ierr = VecGetArray(scopflow->Gu,&gu);CHKERRQ(ierr);

  /* Create IPOPT problem */
  scopflowipopt->nlp = CreateIpoptProblem(scopflow->Nx,xl,xu,scopflow->Ncon,gl,gu,scopflowipopt->nnz_jac_g,scopflowipopt->nnz_hes,0,&eval_scopflow_f,
				   &eval_scopflow_g,&eval_scopflow_grad_f,
				   &eval_scopflow_jac_g,&eval_scopflow_h);

  ierr = VecRestoreArray(scopflow->Xl,&xl);CHKERRQ(ierr);
  ierr = VecRestoreArray(scopflow->Xu,&xu);CHKERRQ(ierr);
  ierr = VecRestoreArray(scopflow->Gl,&gl);CHKERRQ(ierr);
  ierr = VecRestoreArray(scopflow->Gu,&gu);CHKERRQ(ierr);

  ierr = VecGetArray(scopflow->X,&x);CHKERRQ(ierr);
  ierr = VecGetArray(scopflow->G,&g);CHKERRQ(ierr);
  ierr = VecGetArray(scopflow->Lambda,&lam);CHKERRQ(ierr);

  /* IPOPT solver options */
  AddIpoptNumOption(scopflowipopt->nlp, (char*)"tol", scopflow->tolerance);
  AddIpoptIntOption(scopflowipopt->nlp, (char*)"max_iter", 5000);
  AddIpoptStrOption(scopflowipopt->nlp, (char*)"mu_strategy",(char*)"monotone");
  AddIpoptStrOption(scopflowipopt->nlp, (char*)"fixed_variable_treatment", (char*)"relax_bounds");
  AddIpoptStrOption(scopflowipopt->nlp,(char*)"inf_pr_output",(char*)"internal");
  AddIpoptNumOption(scopflowipopt->nlp,(char*)"constr_mult_init_max", 0.0);
  AddIpoptNumOption(scopflowipopt->nlp,(char*)"residual_ratio_max", 1e3);
  AddIpoptNumOption(scopflowipopt->nlp,(char*)"residual_ratio_singular",1e4);


  /** Add intermediate callback to get solver info.
   * Called by IPOPT each iteration. 
  */

  SetIntermediateCallback(scopflowipopt->nlp,SCOPFLOWSolverMonitor_IPOPT);

  /* Solve */
  scopflowipopt->solve_status = IpoptSolve(scopflowipopt->nlp,x,g,&scopflow->obj,lam,NULL,NULL,scopflow);

  ierr = VecRestoreArray(scopflow->X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(scopflow->G,&g);CHKERRQ(ierr);
  ierr = VecRestoreArray(scopflow->Lambda,&lam);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverDestroy_IPOPT(SCOPFLOW scopflow)
{
  PetscErrorCode     ierr;
  SCOPFLOWSolver_IPOPT ipopt = (SCOPFLOWSolver_IPOPT)scopflow->solver;

  PetscFunctionBegin;

  if(ipopt->nlp) {
    FreeIpoptProblem(ipopt->nlp);
    ipopt->nlp = NULL;
  }

  ierr = PetscFree(ipopt);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverSetUp_IPOPT(SCOPFLOW scopflow)
{
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverGetObjective_IPOPT(SCOPFLOW scopflow,PetscReal *obj)
{
  PetscFunctionBegin;
  *obj = scopflow->obj;
  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverGetSolution_IPOPT(SCOPFLOW scopflow,PetscInt cont_num,Vec *X)
{
  PetscErrorCode ierr;
  TCOPFLOW       tcopflow;
  OPFLOW         opflow;
  Vec            Xi;
  PetscInt       nxi;
  PetscScalar    *xi,*x;
  PetscInt       ix=scopflow->xstarti[cont_num];

  PetscFunctionBegin;
  if(!scopflow->ismultiperiod) {
    opflow = scopflow->opflows[cont_num];
    nxi = opflow->nx;
    Xi  = opflow->X;
  } else {
    tcopflow = scopflow->tcopflows[cont_num];
    nxi = tcopflow->Nx;
    Xi  = tcopflow->X;
  }
    
  ierr = VecGetArray(Xi,&xi);CHKERRQ(ierr);
  ierr = VecGetArray(scopflow->X,&x);CHKERRQ(ierr);

  ierr = PetscArraycpy(xi,x+ix,nxi);CHKERRQ(ierr);

  ierr = VecRestoreArray(Xi,&xi);CHKERRQ(ierr);
  ierr = VecRestoreArray(scopflow->X,&x);CHKERRQ(ierr);

  *X = Xi;
  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverGetConstraints_IPOPT(SCOPFLOW scopflow,PetscInt cont_num,Vec *G)
{
  PetscErrorCode ierr;
  TCOPFLOW       tcopflow;
  OPFLOW         opflow;
  Vec            Gi;
  PetscInt       ngi;
  PetscScalar    *gi,*g;
  PetscInt       ig=scopflow->gstarti[cont_num];

  PetscFunctionBegin;
  if(!scopflow->ismultiperiod) {
    opflow = scopflow->opflows[cont_num];
    ngi    = opflow->ncon;
    Gi = opflow->G;
  } else {
    tcopflow = scopflow->tcopflows[cont_num];
    Gi  = tcopflow->G;
    ngi      = tcopflow->Ncon;
  }

  ierr = VecGetArray(Gi,&gi);CHKERRQ(ierr);
  ierr = VecGetArray(scopflow->G,&g);CHKERRQ(ierr);

  ierr = PetscArraycpy(gi,g+ig,ngi);CHKERRQ(ierr);

  ierr = VecRestoreArray(Gi,&gi);CHKERRQ(ierr);
  ierr = VecRestoreArray(scopflow->G,&g);CHKERRQ(ierr);

  *G = Gi;

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverGetConstraintMultipliers_IPOPT(SCOPFLOW scopflow,PetscInt cont_num,Vec *Lambda)
{
  PetscErrorCode ierr;
  TCOPFLOW       tcopflow;
  OPFLOW         opflow;
  Vec            Lambdai;
  PetscInt       ngi;
  PetscScalar    *lambdai,*lambda;
  PetscInt       ig=scopflow->gstarti[cont_num];

  PetscFunctionBegin;
  if(!scopflow->ismultiperiod) {
    opflow = scopflow->opflows[cont_num];
    ngi    = opflow->ncon;
    Lambdai = opflow->Lambda;
  } else {
    tcopflow = scopflow->tcopflows[cont_num];
    Lambdai  = tcopflow->Lambda;
    ngi      = tcopflow->Ncon;
  }

  ierr = VecGetArray(Lambdai,&lambdai);CHKERRQ(ierr);
  ierr = VecGetArray(scopflow->Lambda,&lambda);CHKERRQ(ierr);

  ierr = PetscArraycpy(lambdai,lambda+ig,ngi);CHKERRQ(ierr);

  ierr = VecRestoreArray(Lambdai,&lambdai);CHKERRQ(ierr);
  ierr = VecRestoreArray(scopflow->Lambda,&lambda);CHKERRQ(ierr);

  *Lambda = Lambdai;

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverGetConvergenceStatus_IPOPT(SCOPFLOW scopflow,PetscBool *status)
{
  SCOPFLOWSolver_IPOPT ipopt = (SCOPFLOWSolver_IPOPT)scopflow->solver;

  PetscFunctionBegin;
  if(ipopt->solve_status < 2) *status = PETSC_TRUE; /* See IpReturnCodes_inc.h in IPOPT. The first two denote convergence */
  else *status = PETSC_FALSE;

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverCreate_IPOPT(SCOPFLOW scopflow)
{
  PetscErrorCode ierr;
  SCOPFLOWSolver_IPOPT ipopt;
  
  PetscFunctionBegin;

  if(scopflow->comm->size > 1) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"IPOPT solver does not support execution in parallel\n",scopflow->comm->size); 
  ierr = PetscCalloc1(1,&ipopt);CHKERRQ(ierr);

  ipopt->nlp = NULL;
  ipopt->nnz_jac_g = 0;
  ipopt->nnz_hes = 0;
  scopflow->solver = ipopt;

  scopflow->solverops.setup = SCOPFLOWSolverSetUp_IPOPT;
  scopflow->solverops.solve = SCOPFLOWSolverSolve_IPOPT;
  scopflow->solverops.destroy = SCOPFLOWSolverDestroy_IPOPT;
  scopflow->solverops.getobjective = SCOPFLOWSolverGetObjective_IPOPT;
  scopflow->solverops.getsolution  = SCOPFLOWSolverGetSolution_IPOPT;
  scopflow->solverops.getconvergencestatus = SCOPFLOWSolverGetConvergenceStatus_IPOPT;
  scopflow->solverops.getconstraints = SCOPFLOWSolverGetConstraints_IPOPT;
  scopflow->solverops.getconstraintmultipliers = SCOPFLOWSolverGetConstraintMultipliers_IPOPT;

  PetscFunctionReturn(0);
}

#endif
