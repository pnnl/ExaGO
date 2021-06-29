#include <exago_config.h>
#if defined(EXAGO_ENABLE_IPOPT)

#include <private/opflowimpl.h>
#include <private/scopflowimpl.h>
#include <private/sopflowimpl.h>
#include "sopflow-ipopt.h"

/* IPOPT callback functions */
Bool eval_sopflow_f(PetscInt n, PetscScalar* x, Bool new_x,
            PetscScalar* obj_value, UserDataPtr user_data)
{
  PetscErrorCode ierr;
  SOPFLOW  sopflow=(SOPFLOW)user_data;

  *obj_value = 0.0;

  ierr = VecPlaceArray(sopflow->X,x);CHKERRQ(ierr);
  ierr = (*sopflow->modelops.computeobjective)(sopflow,sopflow->X,obj_value);CHKERRQ(ierr);
  ierr = VecResetArray(sopflow->X);CHKERRQ(ierr);
  
  return TRUE;
}

Bool eval_sopflow_grad_f(PetscInt n, PetscScalar* x, Bool new_x,
                 PetscScalar* grad_f, UserDataPtr user_data)
{

  PetscErrorCode ierr;
  SOPFLOW  sopflow=(SOPFLOW)user_data;

  ierr = VecPlaceArray(sopflow->X,x);CHKERRQ(ierr);
  ierr = VecPlaceArray(sopflow->gradobj,grad_f);CHKERRQ(ierr);
  ierr = (*sopflow->modelops.computegradient)(sopflow,sopflow->X,sopflow->gradobj);CHKERRQ(ierr);
  ierr = VecResetArray(sopflow->X);CHKERRQ(ierr);
  ierr = VecResetArray(sopflow->gradobj);CHKERRQ(ierr);
  
  return TRUE;
}

Bool eval_sopflow_g(PetscInt n, PetscScalar* x, Bool new_x,
            PetscInt m, PetscScalar* g, UserDataPtr user_data)
{
  PetscErrorCode ierr;
  SOPFLOW  sopflow=(SOPFLOW)user_data;

  ierr = VecPlaceArray(sopflow->X,x);CHKERRQ(ierr);
  ierr = VecPlaceArray(sopflow->G,g);CHKERRQ(ierr);
  ierr = (*sopflow->modelops.computeconstraints)(sopflow,sopflow->X,sopflow->G);CHKERRQ(ierr);
  ierr = VecResetArray(sopflow->X);CHKERRQ(ierr);
  ierr = VecResetArray(sopflow->G);CHKERRQ(ierr);

  return TRUE;
}

Bool eval_sopflow_jac_g(PetscInt n, PetscScalar *x, Bool new_x,
                PetscInt m, PetscInt nele_jac,
                PetscInt *iRow, PetscInt *jCol, PetscScalar *values,
                UserDataPtr user_data)
{

  PetscErrorCode ierr;
  SOPFLOW       sopflow=(SOPFLOW)user_data;
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
    ierr = (*sopflow->modelops.computejacobian)(sopflow,sopflow->X,sopflow->Jac);CHKERRQ(ierr);

    ierr = MatGetSize(sopflow->Jac,&nrow,&ncol);CHKERRQ(ierr);
    
    /* Copy over locations to triplet format */
    for(i=0; i < nrow; i++) {
      ierr = MatGetRow(sopflow->Jac,i,&nvals,&cols,&vals);CHKERRQ(ierr);
      for(j=0; j < nvals; j++) {
	iRowstart[j] = roffset + i;
	jColstart[j] = coffset + cols[j];
      }
      /* Increment iRow,jCol pointers */
      iRowstart += nvals;
      jColstart += nvals;
      ierr = MatRestoreRow(sopflow->Jac,i,&nvals,&cols,&vals);CHKERRQ(ierr);
    }
  } else {
    ierr = VecPlaceArray(sopflow->X,x);CHKERRQ(ierr);
    /* Compute jacobian */
    ierr = (*sopflow->modelops.computejacobian)(sopflow,sopflow->X,sopflow->Jac);CHKERRQ(ierr);

    ierr = MatGetSize(sopflow->Jac,&nrow,&ncol);CHKERRQ(ierr);
    /* Copy over values */
    for(i=0; i < nrow; i++) {
      ierr = MatGetRow(sopflow->Jac,i,&nvals,&cols,&vals);CHKERRQ(ierr);
      for(j=0; j < nvals; j++) {
	values[j] = vals[j];
      }
      values += nvals;
      ierr = MatRestoreRow(sopflow->Jac,i,&nvals,&cols,&vals);CHKERRQ(ierr);
    }
    ierr = VecResetArray(sopflow->X);CHKERRQ(ierr);
  }

  return TRUE;
}


Bool eval_sopflow_h(PetscInt n, PetscScalar *x, Bool new_x, PetscScalar obj_factor,
            PetscInt m, PetscScalar *lambda, Bool new_lambda,
            PetscInt nele_hess, PetscInt *iRow, PetscInt *jCol,
            PetscScalar *values, UserDataPtr user_data)
{
  PetscErrorCode ierr;
  PetscInt       nrow;
  SOPFLOW         sopflow=(SOPFLOW)user_data;
  PetscInt       j,k;
  const PetscInt *cols;
  const PetscScalar *vals;
  PetscInt ctr=0;
  PetscInt nvals;

  sopflow->obj_factor = obj_factor;

  ierr = MatGetSize(sopflow->Hes,&nrow,NULL);CHKERRQ(ierr);

  if(values == NULL) {
    /* Set locations only */

    /* Compute Hessian */
    ierr = (*sopflow->modelops.computehessian)(sopflow,sopflow->X,sopflow->Lambda,sopflow->Hes);CHKERRQ(ierr);

    /* Copy over locations to triplet format */
    /* Note that IPOPT
       requires a lower diagonal Hessian (see note https://coin-or.github.io/Ipopt/INTERFACES.html#INTERFACE_CODE)
       Hence, we only add lower diagonal locations
    */
    for(j=0; j < nrow; j++) {
      ierr = MatGetRow(sopflow->Hes,j,&nvals,&cols,&vals);CHKERRQ(ierr);
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
      ierr = MatRestoreRow(sopflow->Hes,j,&nvals,&cols,&vals);CHKERRQ(ierr);
    }
  } else {
    /* Copy values */
    ierr = VecPlaceArray(sopflow->X,x);CHKERRQ(ierr);
    ierr = VecPlaceArray(sopflow->Lambda,lambda);CHKERRQ(ierr);

    /* Compute Hessian */
    ierr = (*sopflow->modelops.computehessian)(sopflow,sopflow->X,sopflow->Lambda,sopflow->Hes);CHKERRQ(ierr);

    /* Copy over values to triplet format */
    /* Note that IPOPT
       requires a lower diagonal Hessian (see note https://coin-or.github.io/Ipopt/INTERFACES.html#INTERFACE_CODE)
       Hence, we only add lower diagonal locations
    */
    for(j=0; j < nrow; j++) {
      ierr = MatGetRow(sopflow->Hes,j,&nvals,&cols,&vals);CHKERRQ(ierr);
      ctr = 0;
      for(k=0; k < nvals; k++) {
	if(cols[k] >= j) { /* upper triangle */
	  /* save as lower triangle locations */
	  values[ctr] = vals[k];
	  ctr++;
	}
      }
      values += ctr;
      ierr = MatRestoreRow(sopflow->Hes,j,&nvals,&cols,&vals);CHKERRQ(ierr);
    }

    ierr = VecResetArray(sopflow->X);CHKERRQ(ierr);
    ierr = VecResetArray(sopflow->Lambda);CHKERRQ(ierr);
  }

  return TRUE;
}

Bool SOPFLOWSolverMonitor_IPOPT(Index alg_mod,Index iter_count,Number obj_value,Number inf_pr,
			       Number inf_du,Number mu,Number d_norm,Number regularization_size,
			       Number alpha_du,Number alpha_pr,Index ls_trials,
                               UserDataPtr user_data)
{
  SOPFLOW  sopflow=(SOPFLOW)user_data;
  sopflow->numiter = iter_count;
  return 1;
}

PetscErrorCode SOPFLOWSolverSolve_IPOPT(SOPFLOW sopflow)
{
  PetscErrorCode     ierr;
  SOPFLOWSolver_IPOPT sopflowipopt = (SOPFLOWSolver_IPOPT)sopflow->solver;
  MatInfo            info_jac,info_hes;
  PetscScalar        *x,*g,*xl,*xu,*gl,*gu,*lam;

  PetscFunctionBegin;

  sopflowipopt->nnz_jac_g = sopflowipopt->nnz_hes = 0;

  /* Compute nonzeros for the Jacobian */
  ierr = (*sopflow->modelops.computejacobian)(sopflow,sopflow->X,sopflow->Jac);CHKERRQ(ierr);
  ierr = MatSetOption(sopflow->Jac,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE);CHKERRQ(ierr);

  ierr = MatGetInfo(sopflow->Jac,MAT_LOCAL,&info_jac);CHKERRQ(ierr);
  sopflowipopt->nnz_jac_g = info_jac.nz_used;

  /* Compute non-zeros for Hessian */
  ierr = (*sopflow->modelops.computehessian)(sopflow,sopflow->X,sopflow->Lambda,sopflow->Hes);CHKERRQ(ierr);
  ierr = MatSetOption(sopflow->Hes,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE);CHKERRQ(ierr);

  ierr  = MatGetInfo(sopflow->Hes,MAT_LOCAL,&info_hes);CHKERRQ(ierr);
  sopflowipopt->nnz_hes = (info_hes.nz_used  - sopflow->Nx)/2 + sopflow->Nx;

  ierr = VecGetArray(sopflow->Xl,&xl);CHKERRQ(ierr);
  ierr = VecGetArray(sopflow->Xu,&xu);CHKERRQ(ierr);
  ierr = VecGetArray(sopflow->Gl,&gl);CHKERRQ(ierr);
  ierr = VecGetArray(sopflow->Gu,&gu);CHKERRQ(ierr);

  /* Create IPOPT problem */
  sopflowipopt->nlp = CreateIpoptProblem(sopflow->Nx,xl,xu,sopflow->Ncon,gl,gu,sopflowipopt->nnz_jac_g,sopflowipopt->nnz_hes,0,&eval_sopflow_f,
				   &eval_sopflow_g,&eval_sopflow_grad_f,
				   &eval_sopflow_jac_g,&eval_sopflow_h);

  ierr = VecRestoreArray(sopflow->Xl,&xl);CHKERRQ(ierr);
  ierr = VecRestoreArray(sopflow->Xu,&xu);CHKERRQ(ierr);
  ierr = VecRestoreArray(sopflow->Gl,&gl);CHKERRQ(ierr);
  ierr = VecRestoreArray(sopflow->Gu,&gu);CHKERRQ(ierr);

  /* IPOPT tolerance */
  {
    static char tol[] = "tol";
    AddIpoptNumOption(sopflowipopt->nlp,tol,sopflow->tolerance);
  }

  ierr = VecGetArray(sopflow->X,&x);CHKERRQ(ierr);
  ierr = VecGetArray(sopflow->G,&g);CHKERRQ(ierr);
  ierr = VecGetArray(sopflow->Lambda,&lam);CHKERRQ(ierr);

  /* Add intermediate callback to get solver info
     This is called by IPOPT each iteration
  */
  SetIntermediateCallback(sopflowipopt->nlp,SOPFLOWSolverMonitor_IPOPT);

  /* Solve */
  sopflowipopt->solve_status = IpoptSolve(sopflowipopt->nlp,x,g,&sopflow->obj,lam,NULL,NULL,sopflow);

  ierr = VecRestoreArray(sopflow->X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(sopflow->G,&g);CHKERRQ(ierr);
  ierr = VecRestoreArray(sopflow->Lambda,&lam);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSolverDestroy_IPOPT(SOPFLOW sopflow)
{
  PetscErrorCode     ierr;
  SOPFLOWSolver_IPOPT ipopt = (SOPFLOWSolver_IPOPT)sopflow->solver;

  PetscFunctionBegin;

  if(ipopt->nlp) {
    FreeIpoptProblem(ipopt->nlp);
    ipopt->nlp = NULL;
  }

  ierr = PetscFree(ipopt);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSolverSetUp_IPOPT(SOPFLOW sopflow)
{
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSolverGetObjective_IPOPT(SOPFLOW sopflow,PetscReal *obj)
{
  PetscFunctionBegin;
  *obj = sopflow->obj;
  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSolverGetSolution_IPOPT(SOPFLOW sopflow,PetscInt scen_num,Vec *X)
{
  PetscErrorCode ierr;
  SCOPFLOW       scopflow;
  OPFLOW         opflow;
  Vec            Xi;
  PetscInt       nxi;
  PetscScalar    *xi,*x;
  PetscInt       ix=sopflow->xstarti[scen_num];

  PetscFunctionBegin;
  if(!sopflow->ismulticontingency) {
    opflow = sopflow->opflows[scen_num];
    nxi = opflow->nx;
    Xi  = opflow->X;
  } else {
    scopflow = sopflow->scopflows[scen_num];
    nxi = scopflow->nx;
    Xi  = scopflow->X;
  }

  ierr = VecGetArray(Xi,&xi);CHKERRQ(ierr);
  ierr = VecGetArray(sopflow->X,&x);CHKERRQ(ierr);

  ierr = PetscArraycpy(xi,x+ix,nxi);CHKERRQ(ierr);

  ierr = VecRestoreArray(Xi,&xi);CHKERRQ(ierr);
  ierr = VecRestoreArray(sopflow->X,&x);CHKERRQ(ierr);

  *X = Xi;
  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSolverGetConstraints_IPOPT(SOPFLOW sopflow,PetscInt scen_num,Vec *G)
{
  PetscErrorCode ierr;
  SCOPFLOW       scopflow;
  OPFLOW         opflow;
  Vec            Gi;
  PetscInt       ngi;
  PetscScalar    *gi,*g;
  PetscInt       ig=sopflow->gstarti[scen_num];;

  PetscFunctionBegin;
  if(!sopflow->ismulticontingency) {
    opflow = sopflow->opflows[scen_num];
    ngi    = opflow->ncon;
    Gi = opflow->G;
  } else {
    scopflow = sopflow->scopflows[scen_num];
    Gi  = scopflow->G;
    ngi      = scopflow->ncon;
  }

  ierr = VecGetArray(Gi,&gi);CHKERRQ(ierr);
  ierr = VecGetArray(sopflow->G,&g);CHKERRQ(ierr);

  ierr = PetscArraycpy(gi,g+ig,ngi);CHKERRQ(ierr);

  ierr = VecRestoreArray(Gi,&gi);CHKERRQ(ierr);
  ierr = VecRestoreArray(sopflow->G,&g);CHKERRQ(ierr);

  *G = Gi;

  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSolverGetConstraintMultipliers_IPOPT(SOPFLOW sopflow,PetscInt scen_num,Vec *Lambda)
{
  PetscErrorCode ierr;
  SCOPFLOW       scopflow;
  OPFLOW         opflow;
  Vec            Lambdai;
  PetscInt       ngi;
  PetscScalar    *lambdai,*lambda;
  PetscInt       ig=sopflow->gstarti[scen_num];

  PetscFunctionBegin;
  if(!sopflow->ismulticontingency) {
    opflow = sopflow->opflows[scen_num];
    ngi    = opflow->ncon;
    Lambdai = opflow->Lambda;
  } else {
    scopflow = sopflow->scopflows[scen_num];
    Lambdai  = scopflow->Lambda;
    ngi      = scopflow->ncon;
  }

  ierr = VecGetArray(Lambdai,&lambdai);CHKERRQ(ierr);
  ierr = VecGetArray(sopflow->Lambda,&lambda);CHKERRQ(ierr);

  ierr = PetscArraycpy(lambdai,lambda+ig,ngi);CHKERRQ(ierr);

  ierr = VecRestoreArray(Lambdai,&lambdai);CHKERRQ(ierr);
  ierr = VecRestoreArray(sopflow->Lambda,&lambda);CHKERRQ(ierr);

  *Lambda = Lambdai;

  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSolverGetConvergenceStatus_IPOPT(SOPFLOW sopflow,PetscBool *status)
{
  SOPFLOWSolver_IPOPT ipopt = (SOPFLOWSolver_IPOPT)sopflow->solver;

  PetscFunctionBegin;
  if(ipopt->solve_status < 2) *status = PETSC_TRUE; /* See IpReturnCodes_inc.h in IPOPT. The first two denote convergence */
  else *status = PETSC_FALSE;

  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSolverCreate_IPOPT(SOPFLOW sopflow)
{
  PetscErrorCode ierr;
  SOPFLOWSolver_IPOPT ipopt;
  
  PetscFunctionBegin;

  if(sopflow->comm->size > 1) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"IPOPT solver does not support execution in parallel\n",sopflow->comm->size); 
  ierr = PetscCalloc1(1,&ipopt);CHKERRQ(ierr);

  ipopt->nlp = NULL;
  ipopt->nnz_jac_g = 0;
  ipopt->nnz_hes = 0;
  sopflow->solver = ipopt;

  sopflow->solverops.setup = SOPFLOWSolverSetUp_IPOPT;
  sopflow->solverops.solve = SOPFLOWSolverSolve_IPOPT;
  sopflow->solverops.destroy = SOPFLOWSolverDestroy_IPOPT;
  sopflow->solverops.getobjective = SOPFLOWSolverGetObjective_IPOPT;
  sopflow->solverops.getsolution  = SOPFLOWSolverGetSolution_IPOPT;
  sopflow->solverops.getconvergencestatus = SOPFLOWSolverGetConvergenceStatus_IPOPT;
  sopflow->solverops.getconstraints = SOPFLOWSolverGetConstraints_IPOPT;
  sopflow->solverops.getconstraintmultipliers = SOPFLOWSolverGetConstraintMultipliers_IPOPT;

  PetscFunctionReturn(0);
}

#endif
