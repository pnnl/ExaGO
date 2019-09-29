#if defined(SCOPFLOW_HAVE_IPOPT)

#include <private/opflowimpl.h>
#include "opflow-ipopt.h"

/* IPOPT callback functions */
Bool eval_opflow_f(PetscInt n, PetscScalar* x, Bool new_x,
            PetscScalar* obj_value, UserDataPtr user_data);

Bool eval_opflow_grad_f(PetscInt n, PetscScalar* x, Bool new_x,
                 PetscScalar* grad_f, UserDataPtr user_data);

Bool eval_opflow_g(PetscInt n, PetscScalar* x, Bool new_x,
            PetscInt m, PetscScalar* g, UserDataPtr user_data);

Bool eval_opflow_jac_g(PetscInt n, PetscScalar *x, Bool new_x,
                PetscInt m, PetscInt nele_jac,
                PetscInt *iRow, PetscInt *jCol, PetscScalar *values,
                UserDataPtr user_data);

Bool eval_opflow_h(PetscInt n, PetscScalar *x, Bool new_x, PetscScalar obj_factor,
            PetscInt m, PetscScalar *lambda, Bool new_lambda,
            PetscInt nele_hess, PetscInt *iRow, PetscInt *jCol,
            PetscScalar *values, UserDataPtr user_data);


static int CCMatrixToMatrixMarketValuesOnly(CCMatrix ccmatrix,PetscInt nz,PetscScalar *values)
{
  PetscErrorCode ierr;

  ierr = PetscMemcpy(values,ccmatrix->values,nz*sizeof(PetscScalar));

  return 0;
}

static int CCMatrixToMatrixMarketLocationsOnly(CCMatrix ccmatrix,PetscInt ncol,PetscInt *iRow,PetscInt *jCol,PetscInt roffset,PetscInt coffset,PetscInt nval)
{
  PetscInt *rowidx;
  PetscInt *colptr;
  PetscInt j,k,ctr=0;
  
  rowidx = ccmatrix->rowidx;
  colptr = ccmatrix->colptr;

  /* Copy from compressed column to (row,col,val) format */
  for(j=0; j < ncol; j++) {
    for(k=colptr[j]; k < colptr[j+1]; k++) {
      iRow[ctr] = rowidx[k] + roffset;
      jCol[ctr] = j + coffset;
      ctr++;
    }
  }
  if(ctr != nval) SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,"Incorrect number of entries ctr = %d given = %d\n",ctr,nval);

  return 0;
}

PetscErrorCode OPFLOWSolverSolve_IPOPT(OPFLOW opflow)
{
  OPFLOWSolver_IPOPT ipopt=opflow->solver;

  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverDestroy_IPOPT(OPFLOW opflow)
{
  PetscErrorCode     ierr;
  OPFLOWSolver_IPOPT ipopt=opflow->solver;

  PetscFunctionBegin;

  if(ipopt->nlp) {
    FreeIpoptProblem(ipopt->nlp);
  }
  ierr = PetscFree(ipopt);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverSetUp_IPOPT(OPFLOW opflow)
{
  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolverCreate_IPOPT(OPFLOW opflow)
{
  PetscErrorCode ierr;
  OPFLOWSolver_IPOPT ipopt;
  
  PetscFunctionBegin;
  ierr = PetscCalloc1(1,&ipopt);CHKERRQ(ierr);

  ipopt->nlp = NULL;
  ipopt->nnz_jac_g = 0;
  ipopt->nnz_hes = 0;
  opflow->solver = ipopt;

  opflow->solverops.setup = OPFLOWSolverSetUp_IPOPT;
  opflow->solverops.solve = OPFLOWSolverSolve_IPOPT;
  opflow->solverops.destroy = OPFLOWSolverDestroy_IPOPT;


  PetscFunctionReturn(0);
}

#endif
