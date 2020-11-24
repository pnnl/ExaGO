#include <exago_config.h>
#if defined(EXAGO_HAVE_IPOPT)

#include <private/opflowimpl.h>
#include <private/scopflowimpl.h>
#include "scopflow-ipopt2.h"

/* IPOPT callback functions */
Bool evalnew_scopflow_f(PetscInt n, PetscScalar* x, Bool new_x,
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

Bool evalnew_scopflow_grad_f(PetscInt n, PetscScalar* x, Bool new_x,
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

Bool evalnew_scopflow_g(PetscInt n, PetscScalar* x, Bool new_x,
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

Bool evalnew_scopflow_jac_g(PetscInt n, PetscScalar *x, Bool new_x,
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


Bool evalnew_scopflow_h(PetscInt n, PetscScalar *x, Bool new_x, PetscScalar obj_factor,
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

PetscErrorCode SCOPFLOWSolverSolve_IPOPTNEW(SCOPFLOW scopflow)
{
  PetscErrorCode     ierr;
  SCOPFLOWSolver_IPOPTNEW scopflowipoptnew = (SCOPFLOWSolver_IPOPTNEW)scopflow->solver;
  MatInfo            info_jac,info_hes;
  PetscScalar        *x,*g,*xl,*xu,*gl,*gu,*lam;

  PetscFunctionBegin;

  scopflowipoptnew->nnz_jac_g = scopflowipoptnew->nnz_hes = 0;

  /* Compute nonzeros for the Jacobian */
  ierr = (*scopflow->modelops.computejacobian)(scopflow,scopflow->X,scopflow->Jac);CHKERRQ(ierr);
  ierr = MatSetOption(scopflow->Jac,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE);CHKERRQ(ierr);

  ierr = MatGetInfo(scopflow->Jac,MAT_LOCAL,&info_jac);CHKERRQ(ierr);
  scopflowipoptnew->nnz_jac_g = info_jac.nz_used;

  /* Compute non-zeros for Hessian */
  ierr = (*scopflow->modelops.computehessian)(scopflow,scopflow->X,scopflow->Lambda,scopflow->Hes);CHKERRQ(ierr);
  ierr = MatSetOption(scopflow->Hes,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE);CHKERRQ(ierr);

  ierr  = MatGetInfo(scopflow->Hes,MAT_LOCAL,&info_hes);CHKERRQ(ierr);
  scopflowipoptnew->nnz_hes = (info_hes.nz_used  - scopflow->Nx)/2 + scopflow->Nx;

  ierr = VecGetArray(scopflow->Xl,&xl);CHKERRQ(ierr);
  ierr = VecGetArray(scopflow->Xu,&xu);CHKERRQ(ierr);
  ierr = VecGetArray(scopflow->Gl,&gl);CHKERRQ(ierr);
  ierr = VecGetArray(scopflow->Gu,&gu);CHKERRQ(ierr);

  /* Create IPOPT problem */
  scopflowipoptnew->nlp = CreateIpoptProblem(scopflow->Nx,xl,xu,scopflow->Ncon,gl,gu,scopflowipoptnew->nnz_jac_g,scopflowipoptnew->nnz_hes,0,&evalnew_scopflow_f,
				   &evalnew_scopflow_g,&evalnew_scopflow_grad_f,
				   &evalnew_scopflow_jac_g,&evalnew_scopflow_h);

  ierr = VecRestoreArray(scopflow->Xl,&xl);CHKERRQ(ierr);
  ierr = VecRestoreArray(scopflow->Xu,&xu);CHKERRQ(ierr);
  ierr = VecRestoreArray(scopflow->Gl,&gl);CHKERRQ(ierr);
  ierr = VecRestoreArray(scopflow->Gu,&gu);CHKERRQ(ierr);

  ierr = VecGetArray(scopflow->X,&x);CHKERRQ(ierr);
  ierr = VecGetArray(scopflow->G,&g);CHKERRQ(ierr);
  ierr = VecGetArray(scopflow->Lambda,&lam);CHKERRQ(ierr);

  /* Solve */
  scopflowipoptnew->solve_status = IpoptSolve(scopflowipoptnew->nlp,x,g,&scopflow->obj,lam,NULL,NULL,scopflow);

  ierr = VecRestoreArray(scopflow->X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(scopflow->G,&g);CHKERRQ(ierr);
  ierr = VecRestoreArray(scopflow->Lambda,&lam);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverDestroy_IPOPTNEW(SCOPFLOW scopflow)
{
  PetscErrorCode     ierr;
  SCOPFLOWSolver_IPOPTNEW ipoptnew = (SCOPFLOWSolver_IPOPTNEW)scopflow->solver;

  PetscFunctionBegin;

  if(ipoptnew->nlp) {
    FreeIpoptProblem(ipoptnew->nlp);
    ipoptnew->nlp = NULL;
  }

  ierr = PetscFree(ipoptnew);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverSetUp_IPOPTNEW(SCOPFLOW scopflow)
{
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverGetObjective_IPOPTNEW(SCOPFLOW scopflow,PetscReal *obj)
{
  PetscFunctionBegin;
  *obj = scopflow->obj;
  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverGetSolution_IPOPTNEW(SCOPFLOW scopflow,PetscInt cont_num,Vec *X)
{
  PetscErrorCode ierr;
  OPFLOW         opflow=scopflow->opflows[cont_num];
  Vec            Xi=opflow->X;
  PetscInt       nxi=opflow->nx;
  PetscScalar    *xi,*x;
  PetscInt       ix=scopflow->xstarti[cont_num];

  PetscFunctionBegin;
  ierr = VecGetArray(Xi,&xi);CHKERRQ(ierr);
  ierr = VecGetArray(scopflow->X,&x);CHKERRQ(ierr);

  ierr = PetscArraycpy(xi,x+ix,nxi);CHKERRQ(ierr);

  ierr = VecRestoreArray(Xi,&xi);CHKERRQ(ierr);
  ierr = VecRestoreArray(scopflow->X,&x);CHKERRQ(ierr);

  *X = Xi;
  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverGetConstraints_IPOPTNEW(SCOPFLOW scopflow,PetscInt cont_num,Vec *G)
{
  PetscErrorCode ierr;
  OPFLOW         opflow=scopflow->opflows[cont_num];
  Vec            Gi=opflow->G;
  PetscInt       ngi=opflow->ncon;
  PetscScalar    *gi,*g;
  PetscInt       ig=scopflow->gstarti[cont_num];

  PetscFunctionBegin;
  ierr = VecGetArray(Gi,&gi);CHKERRQ(ierr);
  ierr = VecGetArray(scopflow->G,&g);CHKERRQ(ierr);

  ierr = PetscArraycpy(gi,g+ig,ngi);CHKERRQ(ierr);

  ierr = VecRestoreArray(Gi,&gi);CHKERRQ(ierr);
  ierr = VecRestoreArray(scopflow->G,&g);CHKERRQ(ierr);

  *G = Gi;

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverGetConstraintMultipliers_IPOPTNEW(SCOPFLOW scopflow,PetscInt cont_num,Vec *Lambda)
{
  PetscErrorCode ierr;
  OPFLOW         opflow=scopflow->opflows[cont_num];
  Vec            Lambdai=opflow->Lambda;
  PetscInt       ngi=opflow->ncon;
  PetscScalar    *lambdai,*lambda;
  PetscInt       ig=scopflow->gstarti[cont_num];

  PetscFunctionBegin;
  ierr = VecGetArray(Lambdai,&lambdai);CHKERRQ(ierr);
  ierr = VecGetArray(scopflow->Lambda,&lambda);CHKERRQ(ierr);

  ierr = PetscArraycpy(lambdai,lambda+ig,ngi);CHKERRQ(ierr);

  ierr = VecRestoreArray(Lambdai,&lambdai);CHKERRQ(ierr);
  ierr = VecRestoreArray(scopflow->Lambda,&lambda);CHKERRQ(ierr);

  *Lambda = Lambdai;

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverGetConvergenceStatus_IPOPTNEW(SCOPFLOW scopflow,PetscBool *status)
{
  SCOPFLOWSolver_IPOPTNEW ipoptnew = (SCOPFLOWSolver_IPOPTNEW)scopflow->solver;

  PetscFunctionBegin;
  if(ipoptnew->solve_status < 2) *status = PETSC_TRUE; /* See IpReturnCodes_inc.h in IPOPT. The first two denote convergence */
  else *status = PETSC_FALSE;

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverCreate_IPOPTNEW(SCOPFLOW scopflow)
{
  PetscErrorCode ierr;
  SCOPFLOWSolver_IPOPTNEW ipoptnew;
  
  PetscFunctionBegin;

  if(scopflow->comm->size > 1) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"IPOPT solver does not support execution in parallel\n",scopflow->comm->size); 
  ierr = PetscCalloc1(1,&ipoptnew);CHKERRQ(ierr);

  ipoptnew->nlp = NULL;
  ipoptnew->nnz_jac_g = 0;
  ipoptnew->nnz_hes = 0;
  scopflow->solver = ipoptnew;

  scopflow->solverops.setup = SCOPFLOWSolverSetUp_IPOPTNEW;
  scopflow->solverops.solve = SCOPFLOWSolverSolve_IPOPTNEW;
  scopflow->solverops.destroy = SCOPFLOWSolverDestroy_IPOPTNEW;
  scopflow->solverops.getobjective = SCOPFLOWSolverGetObjective_IPOPTNEW;
  scopflow->solverops.getsolution  = SCOPFLOWSolverGetSolution_IPOPTNEW;
  scopflow->solverops.getconvergencestatus = SCOPFLOWSolverGetConvergenceStatus_IPOPTNEW;
  scopflow->solverops.getconstraints = SCOPFLOWSolverGetConstraints_IPOPTNEW;
  scopflow->solverops.getconstraintmultipliers = SCOPFLOWSolverGetConstraintMultipliers_IPOPTNEW;

  PetscFunctionReturn(0);
}

#endif
