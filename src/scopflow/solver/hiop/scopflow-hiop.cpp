#include <exago_config.h>

#if defined(EXAGO_ENABLE_HIOP)

#include <private/opflowimpl.h>
#include <private/tcopflowimpl.h>
#include <private/scopflowimpl.h>
#include "scopflow-hiop.hpp"

PetscErrorCode SCOPFLOWBaseAuxObjectiveFunction(OPFLOW opflow,const double* x,double* obj,void* ctx)
{
  PetscErrorCode      ierr;
  SCOPFLOW            scopflow=(SCOPFLOW)ctx;
  SCOPFLOWSolver_HIOP hiop = (SCOPFLOWSolver_HIOP)scopflow->solver;

  PetscFunctionBegin;

  if(hiop->pridecompprob->include_r_) {
    hiop->pridecompprob->rec_evaluator->eval_f(opflow->nx,x,false,*obj);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWBaseAuxGradientFunction(OPFLOW opflow,const double* x,double* grad,void* ctx)
{
  PetscErrorCode ierr;
  SCOPFLOW       scopflow=(SCOPFLOW)ctx;
  SCOPFLOWSolver_HIOP hiop = (SCOPFLOWSolver_HIOP)scopflow->solver;

  PetscFunctionBegin;

  if(hiop->pridecompprob->include_r_) {
    hiop->pridecompprob->rec_evaluator->eval_grad(opflow->nx,x,false,grad);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWBaseAuxHessianFunction(OPFLOW opflow,const double* x,Mat Hess,void* ctx)
{
  PetscErrorCode      ierr;
  SCOPFLOW            scopflow=(SCOPFLOW)ctx;
  SCOPFLOWSolver_HIOP hiop = (SCOPFLOWSolver_HIOP)scopflow->solver;
  PetscInt            row,col;
  PetscScalar         val;
  PetscInt            i;

  PetscFunctionBegin;
  if(hiop->pridecompprob->include_r_) {
    for(i=0; i < hiop->pridecompprob->nxcoup; i++) {
      row = col = hiop->pridecompprob->loc_xcoup[i];
      val = hiop->pridecompprob->rec_evaluator->get_rhess()[i];
      val *= opflow->obj_factor;
      ierr = MatSetValues(Hess,1,&row,1,&col,&val,ADD_VALUES);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

SCOPFLOWHIOPInterface::~SCOPFLOWHIOPInterface()
{
  //    printf("Exiting application\n");
}

SCOPFLOWHIOPInterface::SCOPFLOWHIOPInterface(SCOPFLOW scopflowin)
{
  int  i,k,j=0;
  PS    ps;
  PSBUS bus;
  PSGEN gen;
  OPFLOW opflowbase;

  scopflow = scopflowin;

  rec_evaluator = NULL;

  opflowbase = scopflow->opflow0;
  ps = opflowbase->ps;

  /* Get the number of coupling variables */
  nxcoup = 0;

  for(i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    for(k=0; k < bus->ngen; k++) {
      PSBUSGetGen(bus,k,&gen);
      if(gen->status) nxcoup++;  
    }
  }
  /* Set the coupling indices */
  for(i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    for(k=0; k < bus->ngen; k++) {
      PSBUSGetGen(bus,k,&gen);
      if(gen->status) {
	loc_xcoup.push_back(gen->startxpowloc); /* Location for Pg */
	j++;
      }
    }
  }
  assert(j == nxcoup);
}

hiop::hiopSolveStatus SCOPFLOWHIOPInterface::solve_master(double* x,
                                     const bool& include_r,
                                     const double& rval, 
                                     const double* grad,
                                     const double*hess)
{

  PetscErrorCode ierr;
  double         obj;
  OPFLOW         opflow;
  Vec            X;
  const PetscScalar *xsol;
  include_r_ = include_r;

  opflow = scopflow->opflows[0];

  //  printf("Rank[%d]:Enter solve_master\n",scopflow->comm->rank);
  ierr = OPFLOWSolve(opflow);
  ierr = OPFLOWGetObjective(opflow,&obj);
  ierr = OPFLOWGetSolution(opflow,&X);
  ierr = VecGetArrayRead(X,&xsol);
  ierr = PetscMemcpy(x,xsol,opflow->nx*sizeof(double));
  ierr = VecRestoreArrayRead(X,&xsol);

  printf("Rank[%d]:Exit solve_master\n",scopflow->comm->rank);
  
  return hiop::Solve_Success;
}

bool SCOPFLOWHIOPInterface::set_recourse_approx_evaluator(const int n,hiopInterfacePriDecProblem::RecourseApproxEvaluator* evaluator)
{
  assert(n == nxcoup);
  if(rec_evaluator == NULL) {
    rec_evaluator = new hiopInterfacePriDecProblem::
      RecourseApproxEvaluator(n, evaluator->get_S(), loc_xcoup,
			      evaluator->get_rval(), evaluator->get_rgrad(),
			      evaluator->get_rhess(), evaluator->get_x0());
  }
  
  assert(rec_evaluator->get_rgrad() != NULL);
  rec_evaluator->set_rval(evaluator->get_rval());
  rec_evaluator->set_rgrad(n,evaluator->get_rgrad());
  rec_evaluator->set_rhess(n,evaluator->get_rhess());
  rec_evaluator->set_x0(n,evaluator->get_x0());

  return true;
}

/* Note: x only holds the coupled variables which, in this case, are the generator real power variables
   for the base-case
*/
bool SCOPFLOWHIOPInterface::eval_f_rterm(size_t idx, const int& n, const double* x, double& rval)
{
  PetscErrorCode ierr;
  OPFLOW opflow;
  OPFLOW opflow0; /* base case OPFLOW */
  PS     ps,ps0;
  PSBUS  bus,bus0;
  PSGEN  gen,gen0;
  PetscInt i,k,j=0;
  PetscInt cont_num=idx+1;
  
  //  printf("[Rank %d] contingency %d Came in recourse objective function\n",scopflow->comm->rank,cont_num);

  if(scopflow->cstart <= cont_num && cont_num < scopflow->cend) {
    opflow0 = scopflow->opflow0;
    opflow = scopflow->opflows[cont_num-scopflow->cstart];

    /* Update generator set-points */
    ps = opflow->ps;
    ps0 = opflow0->ps; /* ps0 is only used for accessing locations */ 
    for(i=0; i < ps->nbus; i++) {
      bus = &ps->bus[i];
      bus0 = &ps0->bus[i];
      for(k=0; k < bus->ngen; k++) {
	ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
	ierr = PSBUSGetGen(bus0,k,&gen0);CHKERRQ(ierr);
	if(gen0->status) {
	  gen->pgs = x[j++];
	}
      }
    }
    assert(j == n);
    /* Solve */
    ierr = OPFLOWSolve(opflow);
    ierr = OPFLOWGetObjective(opflow,&rval);
  }
  return true;
}


bool SCOPFLOWHIOPInterface::eval_grad_rterm(size_t idx, const int& n, double* x, double* grad)
{
  PetscErrorCode ierr;
  OPFLOW opflow;
  OPFLOW opflow0; /* base case OPFLOW */
  PS     ps,ps0;
  PSBUS  bus,bus0;
  PSGEN  gen,gen0;
  PetscInt i,k,j=0;
  const PetscScalar *lam,*lameq;
  Vec    Lambda;
  PetscInt cont_num=idx+1;
  
  //  printf("[Rank %d] contingency %d Came in recourse gradient function\n",scopflow->comm->rank,cont_num);
  if(scopflow->cstart <= cont_num && cont_num < scopflow->cend) {
    opflow0 = scopflow->opflow0;
    opflow = scopflow->opflows[cont_num-scopflow->cstart];

    ierr = OPFLOWGetConstraintMultipliers(opflow,&Lambda);
    ierr = VecGetArrayRead(Lambda,&lam);
    lameq = lam;
    /* Update generator set-points */
    ps = opflow->ps;
    ps0 = opflow0->ps; 
    for(i=0; i < ps->nbus; i++) {
      bus = &ps->bus[i];
      bus0 = &ps0->bus[i];
      for(k=0; k < bus->ngen; k++) {
	ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
	ierr = PSBUSGetGen(bus0,k,&gen0);CHKERRQ(ierr);
	if(gen0->status && gen->status) {
	  /* Get the lagrange multiplier for the generator set-point equality constraint x_i - x_0 
	   gradient is the partial derivative for it (note that it is negative) */
	  if(opflow->has_gensetpoint) {
	    grad[j++] = -lameq[gen->starteqloc+1];
	  }
	}
      }
    }
    ierr = VecRestoreArrayRead(Lambda,&lam);
    return true;
  }
}

size_t SCOPFLOWHIOPInterface::get_num_rterms() const
{
  return scopflow->Nc-1;
}

size_t SCOPFLOWHIOPInterface::get_num_vars() const
{
  return scopflow->opflow0->nx;
}

void SCOPFLOWHIOPInterface::get_solution(double* x) const
{
  OPFLOW opflow;
  Vec    X;
  const double *xarr;
  PetscInt i;

  if(scopflow->cstart == 0) {
    opflow = scopflow->opflows[0];
    OPFLOWGetSolution(opflow,&X);
    VecGetArrayRead(X,&xarr);
    for(i=0; i < opflow->nx; i++) x[i] = xarr[i];
    VecRestoreArrayRead(X,&xarr);
  }
}
 
double SCOPFLOWHIOPInterface::get_objective()
{
  OPFLOW opflow;
  double obj;
  PetscErrorCode ierr;

  if(scopflow->cstart == 0) {
    opflow = scopflow->opflows[0];
    ierr = OPFLOWGetObjective(opflow,&obj);CHKERRQ(ierr);
    return obj;
  }
}

PetscErrorCode SCOPFLOWSolverSolve_HIOP(SCOPFLOW scopflow)
{
  PetscErrorCode      ierr;
  PetscInt            c;
  SCOPFLOWSolver_HIOP hiop = (SCOPFLOWSolver_HIOP)scopflow->solver;

  PetscFunctionBegin;

  auto status = hiop->pridecsolver->run();
  /*  if(!scopflow->comm->rank) {
    ierr = VecCopy(scopflow->opflow0->X,scopflow->opflows[0]->X);CHKERRQ(ierr);
    ierr = VecCopy(scopflow->opflow0->G,scopflow->opflows[0]->G);CHKERRQ(ierr);
    ierr = VecCopy(scopflow->opflow0->Lambda,scopflow->opflows[0]->Lambda);CHKERRQ(ierr);
    scopflow->opflows[0]->obj = scopflow->opflow0->obj;
  }
  */
  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverDestroy_HIOP(SCOPFLOW scopflow)
{
  PetscErrorCode     ierr;
  SCOPFLOWSolver_HIOP hiop = (SCOPFLOWSolver_HIOP)scopflow->solver;

  PetscFunctionBegin;

  delete hiop->pridecompprob;
  delete hiop->pridecsolver;

  ierr = PetscFree(hiop);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverGetObjective_HIOP(SCOPFLOW scopflow,PetscReal *obj)
{
  PetscErrorCode ierr;
  PetscReal      temp=0.0;
  PetscFunctionBegin;
  if(!scopflow->comm->rank) {
    if(!scopflow->ismultiperiod) {
      for(int i=0; i < scopflow->nc; i++) {
	temp += scopflow->opflows[i]->obj;
      }
    } else {
      temp = scopflow->tcopflows[0]->obj;
    }
  }
  // ierr = MPI_Bcast(&temp,1,MPI_REAL,0,scopflow->comm->type);CHKERRQ(ierr);
  ierr = MPI_Bcast(&temp,1,MPI_DOUBLE,0,scopflow->comm->type);CHKERRQ(ierr);
  scopflow->obj = temp;
  *obj = scopflow->obj;
  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverGetSolution_HIOP(SCOPFLOW scopflow,PetscInt cont_num,Vec *X)
{
  TCOPFLOW       tcopflow;
  OPFLOW         opflow;

  PetscFunctionBegin;

  if(scopflow->cstart <= cont_num && cont_num < scopflow->cend) {
    if(!scopflow->ismultiperiod) {
      opflow = scopflow->opflows[cont_num-scopflow->cstart];
      *X = opflow->X; 
    } else {
      tcopflow = scopflow->tcopflows[cont_num-scopflow->cstart];
      *X = tcopflow->X; 
    }
  } else {
    *X = NULL;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverGetConstraints_HIOP(SCOPFLOW scopflow,PetscInt cont_num,Vec *G)
{
  TCOPFLOW            tcopflow;
  OPFLOW              opflow;

  PetscFunctionBegin;

  if(scopflow->cstart <= cont_num && cont_num < scopflow->cend) {
    if(!scopflow->ismultiperiod) {
      opflow = scopflow->opflows[cont_num-scopflow->cstart];
      *G = opflow->G; 
    } else {
      tcopflow = scopflow->tcopflows[cont_num-scopflow->cstart];
      *G = tcopflow->G; 
    }      
  } else {
    *G = NULL;
  }

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverGetConstraintMultipliers_HIOP(SCOPFLOW scopflow,PetscInt cont_num,Vec *Lambda)
{
  TCOPFLOW            tcopflow;
  OPFLOW              opflow;

  PetscFunctionBegin;
  if(scopflow->cstart <= cont_num && cont_num < scopflow->cend) {
    if(!scopflow->ismultiperiod) {
      opflow = scopflow->opflows[cont_num-scopflow->cstart];
      *Lambda = opflow->Lambda; 
    } else {
      tcopflow = scopflow->tcopflows[cont_num-scopflow->cstart];
      *Lambda = tcopflow->Lambda; 
    }      
  } else {
    *Lambda = NULL;
  }

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverGetConvergenceStatus_HIOP(SCOPFLOW scopflow,PetscBool *status)
{
  PetscFunctionBegin;

  *status = PETSC_TRUE;

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverSetUp_HIOP(SCOPFLOW scopflow)
{
  SCOPFLOWSolver_HIOP hiop=(SCOPFLOWSolver_HIOP)scopflow->solver;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  hiop->pridecompprob = new SCOPFLOWHIOPInterface(scopflow);

  hiop->pridecsolver = new hiop::hiopAlgPrimalDecomposition(hiop->pridecompprob,hiop->pridecompprob->nxcoup,hiop->pridecompprob->loc_xcoup,scopflow->comm->type);

  /* Add auxillary functions to base case */
  if(scopflow->cstart == 0) {
    ierr = OPFLOWSetAuxillaryObjective(scopflow->opflows[0],SCOPFLOWBaseAuxObjectiveFunction,SCOPFLOWBaseAuxGradientFunction,SCOPFLOWBaseAuxHessianFunction,scopflow);CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

extern "C" {

PetscErrorCode SCOPFLOWSolverCreate_HIOP(SCOPFLOW scopflow)
{
  PetscErrorCode ierr;
  SCOPFLOWSolver_HIOP hiop;
  
  PetscFunctionBegin;

  ierr = PetscCalloc1(1,&hiop);CHKERRQ(ierr);

  scopflow->solver = hiop;

  scopflow->solverops.setup = SCOPFLOWSolverSetUp_HIOP;
  scopflow->solverops.solve = SCOPFLOWSolverSolve_HIOP;
  scopflow->solverops.destroy = SCOPFLOWSolverDestroy_HIOP;
  scopflow->solverops.getobjective = SCOPFLOWSolverGetObjective_HIOP;
  scopflow->solverops.getsolution  = SCOPFLOWSolverGetSolution_HIOP;
  scopflow->solverops.getconvergencestatus = SCOPFLOWSolverGetConvergenceStatus_HIOP;
  scopflow->solverops.getconstraints = SCOPFLOWSolverGetConstraints_HIOP;
  scopflow->solverops.getconstraintmultipliers = SCOPFLOWSolverGetConstraintMultipliers_HIOP;

  PetscFunctionReturn(0);
}

} // end of extern "C"

#endif
