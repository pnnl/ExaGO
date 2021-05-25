#include <exago_config.h>

#if defined(EXAGO_ENABLE_HIOP)
#if defined(EXAGO_ENABLE_HIOP_DISTRIBUTED)

#include <private/opflowimpl.h>
#include <private/scopflowimpl.h>
#include <private/sopflowimpl.h>
#include "sopflow-hiop.hpp"

PetscErrorCode SOPFLOWBaseAuxObjectiveFunction(OPFLOW opflow,const double* x,double* obj,void* ctx)
{
  PetscErrorCode      ierr;
  SOPFLOW            sopflow=(SOPFLOW)ctx;
  SOPFLOWSolver_HIOP hiop = (SOPFLOWSolver_HIOP)sopflow->solver;

  PetscFunctionBegin;

  if(hiop->pridecompprob->include_r_) {
    hiop->pridecompprob->rec_evaluator->eval_f(opflow->nx,x,false,*obj);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWBaseAuxGradientFunction(OPFLOW opflow,const double* x,double* grad,void* ctx)
{
  PetscErrorCode ierr;
  SOPFLOW       sopflow=(SOPFLOW)ctx;
  SOPFLOWSolver_HIOP hiop = (SOPFLOWSolver_HIOP)sopflow->solver;

  PetscFunctionBegin;

  if(hiop->pridecompprob->include_r_) {
    hiop->pridecompprob->rec_evaluator->eval_grad(opflow->nx,x,false,grad);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWBaseAuxHessianFunction(OPFLOW opflow,const double* x,Mat Hess,void* ctx)
{
  PetscErrorCode      ierr;
  SOPFLOW            sopflow=(SOPFLOW)ctx;
  SOPFLOWSolver_HIOP hiop = (SOPFLOWSolver_HIOP)sopflow->solver;
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

PetscErrorCode SOPFLOWBaseAuxSCOPFLOWObjectiveFunction(SCOPFLOW scopflow,const double* x,double* obj,void* ctx)
{
  PetscErrorCode      ierr;
  SOPFLOW            sopflow=(SOPFLOW)ctx;
  SOPFLOWSolver_HIOP hiop = (SOPFLOWSolver_HIOP)sopflow->solver;

  PetscFunctionBegin;

  if(hiop->pridecompprob->include_r_) {
    hiop->pridecompprob->rec_evaluator->eval_f(scopflow->opflows[0]->nx,x,false,*obj);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWBaseAuxSCOPFLOWGradientFunction(SCOPFLOW scopflow,const double* x,double* grad,void* ctx)
{
  PetscErrorCode ierr;
  SOPFLOW       sopflow=(SOPFLOW)ctx;
  SOPFLOWSolver_HIOP hiop = (SOPFLOWSolver_HIOP)sopflow->solver;

  PetscFunctionBegin;

  if(hiop->pridecompprob->include_r_) {
    hiop->pridecompprob->rec_evaluator->eval_grad(scopflow->opflows[0]->nx,x,false,grad);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWBaseAuxSCOPFLOWHessianFunction(SCOPFLOW scopflow,const double* x,Mat Hess,void* ctx)
{
  PetscErrorCode      ierr;
  SOPFLOW            sopflow=(SOPFLOW)ctx;
  SOPFLOWSolver_HIOP hiop = (SOPFLOWSolver_HIOP)sopflow->solver;
  OPFLOW              opflow=scopflow->opflows[0];
  PetscInt            row,col;
  PetscScalar         val;
  PetscInt            i;

  PetscFunctionBegin;
  if(hiop->pridecompprob->include_r_) {
    for(i=0; i < hiop->pridecompprob->nxcoup; i++) {
      row = col = hiop->pridecompprob->loc_xcoup[i];
      val = hiop->pridecompprob->rec_evaluator->get_rhess()[i];
      val *= scopflow->obj_factor;
      ierr = MatSetValues(Hess,1,&row,1,&col,&val,ADD_VALUES);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

SOPFLOWHIOPInterface::~SOPFLOWHIOPInterface()
{
  //    printf("Exiting application\n");
}

SOPFLOWHIOPInterface::SOPFLOWHIOPInterface(SOPFLOW sopflowin)
{
  int  i,k,j=0;
  PS    ps;
  PSBUS bus;
  PSGEN gen;
  OPFLOW opflowbase;

  sopflow = sopflowin;

  rec_evaluator = NULL;

  if(!sopflow->ismulticontingency) {
    opflowbase = sopflow->opflow0;
    ps = opflowbase->ps;
  } else {
    opflowbase = sopflow->scopflows[0]->opflows[0];
    ps = opflowbase->ps;
  }

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
	loc_xcoup.push_back(gen->startxpsetloc); /* Location for Pgset */
	j++;
      }
    }
  }
  assert(j == nxcoup);
}

hiop::hiopSolveStatus SOPFLOWHIOPInterface::solve_master(double* x,
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

  if(!sopflow->ismulticontingency) {
    opflow = sopflow->opflow0;
    
    //  printf("Rank[%d]:Enter solve_master\n",sopflow->comm->rank);
    ierr = OPFLOWSolve(opflow);
    ierr = OPFLOWGetObjective(opflow,&obj);
    ierr = OPFLOWGetSolution(opflow,&X);
    ierr = VecGetArrayRead(X,&xsol);
    ierr = PetscMemcpy(x,xsol,opflow->nx*sizeof(double));
    ierr = VecRestoreArrayRead(X,&xsol);
  } else {
    SCOPFLOW scopflow = sopflow->scopflows[0];
    ierr = SCOPFLOWSolve(scopflow);
    ierr = SCOPFLOWGetObjective(scopflow,&obj);
    ierr = SCOPFLOWGetSolution(scopflow,0,&X);
    ierr = VecGetArrayRead(X,&xsol);
    ierr = PetscMemcpy(x,xsol,scopflow->opflows[0]->nx*sizeof(double));
    ierr = VecRestoreArrayRead(X,&xsol);
  }
  //  printf("Rank[%d]:Exit solve_master\n",sopflow->comm->rank);
  
  return hiop::Solve_Success;
}

bool SOPFLOWHIOPInterface::set_recourse_approx_evaluator(const int n,hiopInterfacePriDecProblem::RecourseApproxEvaluator* evaluator)
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

extern PetscErrorCode SOPFLOWUpdateOPFLOWVariableBounds(OPFLOW,Vec,Vec,void*);

/* Note: x only holds the coupled variables which, in this case, are the generator real power variables
   for the base-case
*/
bool SOPFLOWHIOPInterface::eval_f_rterm(size_t idx, const int& n, const double* x, double& rval)
{
  PetscErrorCode ierr;
  OPFLOW opflow0; /* base case OPFLOW */
  PS     ps,ps0;
  PSBUS  bus,bus0;
  PSGEN  gen,gen0;
  PetscInt i,k,j,g=0;
  PetscInt scen_num=idx+1;
  
  //  printf("[Rank %d] contingency %d Came in recourse objective function\n",sopflow->comm->rank,scen_num);

  if(!sopflow->ismulticontingency) {
    opflow0 = sopflow->opflow0;
    
    ierr = OPFLOWCreate(PETSC_COMM_SELF,&opflowscen);CHKERRQ(ierr);
    ierr = OPFLOWSetModel(opflowscen,OPFLOWMODEL_PBPOL);CHKERRQ(ierr);
    
    ierr = OPFLOWReadMatPowerData(opflowscen,sopflow->netfile);CHKERRQ(ierr);
    /* Set up the PS object for opflow */
    ps = opflowscen->ps;
    ierr = PSSetUp(ps);CHKERRQ(ierr);
    
    if(sopflow->scenfileset) {
      ierr = PSApplyScenario(ps,sopflow->scenlist.scen[scen_num]);
    }
    
    ierr = OPFLOWHasGenSetPoint(opflowscen,PETSC_TRUE);CHKERRQ(ierr); /* Activates ramping variables */
    //  ierr = OPFLOWSetObjectiveType(opflowscen,NO_OBJ);CHKERRQ(ierr);
    ierr = OPFLOWSetUpdateVariableBoundsFunction(opflowscen,SOPFLOWUpdateOPFLOWVariableBounds,(void*)sopflow);
    
    ierr = OPFLOWSetUp(opflowscen);CHKERRQ(ierr);
    
    /* Update generator set-points */
    ps = opflowscen->ps;
    ps0 = opflow0->ps; 
    for(i=0; i < ps->nbus; i++) {
      bus = &ps->bus[i];
      bus0 = &ps0->bus[i];
      for(k=0; k < bus->ngen; k++) {
	ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
	ierr = PSBUSGetGen(bus0,k,&gen0);CHKERRQ(ierr);
	if(gen0->status) {
	  gen0->pgs = gen->pgs = x[g++];
	}
      }
    }
    assert(g == n);
    /* Solve */
    ierr = OPFLOWSolve(opflowscen);
    ierr = OPFLOWGetObjective(opflowscen,&rval);
  } else {

    ierr = SCOPFLOWCreate(PETSC_COMM_SELF,&scopflowscen);CHKERRQ(ierr);
    ierr = SCOPFLOWSetNetworkData(scopflowscen,sopflow->netfile);CHKERRQ(ierr);
    ierr = SCOPFLOWSetContingencyData(scopflowscen,NATIVE,sopflow->scopflows[0]->ctgcfile);CHKERRQ(ierr);

    if(sopflow->scenfileset) {
      ierr = SCOPFLOWSetScenario(scopflowscen,&sopflow->scenlist.scen[scen_num]);CHKERRQ(ierr);
    }

    ierr = SCOPFLOWSetUp(scopflowscen);CHKERRQ(ierr);

    /* Update generator set-points */
    ps = scopflowscen->opflows[0]->ps;
    for(i=0; i < ps->nbus; i++) {
      bus = &ps->bus[i];
      for(k=0; k < bus->ngen; k++) {
	ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
	if(gen->status) {
	  gen->pgs = x[g++];
	}
      }
    }

    assert(g == n);
    /* Solve */
    ierr = SCOPFLOWSolve(scopflowscen);
    ierr = SCOPFLOWGetObjective(scopflowscen,&rval);
  }

  return true;
}

bool SOPFLOWHIOPInterface::eval_grad_rterm(size_t idx, const int& n, double* x, double* grad)
{
  PetscErrorCode ierr;
  OPFLOW opflow;
  OPFLOW opflow0; /* base case OPFLOW */
  PS     ps,ps0;
  PSBUS  bus,bus0;
  PSGEN  gen,gen0;
  PetscInt i,k,j,g=0;
  const PetscScalar *lam,*lameq;
  Vec    Lambda;
  PetscInt scen_num=idx+1;
  
  // printf("[Rank %d] contingency %d Came in recourse gradient function\n",sopflow->comm->rank,scen_num);

  if(!sopflow->ismulticontingency) {
    opflow0 = sopflow->opflow0;
    opflow = opflowscen;
    ierr = OPFLOWGetConstraintMultipliers(opflow,&Lambda);

    ierr = VecGetArrayRead(Lambda,&lam);
    lameq = lam;

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
	    grad[g++] = -lameq[gen->starteqloc+1];
	  }
	}
      }
    }
    ierr = VecRestoreArrayRead(Lambda,&lam);

    ierr = OPFLOWDestroy(&opflow);CHKERRQ(ierr);

  } else {
    opflow = scopflowscen->opflows[0];
    ierr = SCOPFLOWGetConstraintMultipliers(scopflowscen,0,&Lambda);
    ierr = VecGetArrayRead(Lambda,&lam);
    lameq = lam;

    ps = opflow->ps;
    for(i=0; i < ps->nbus; i++) {
      bus = &ps->bus[i];
      for(k=0; k < bus->ngen; k++) {
	ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
	if(gen->status) {
	  /* Get the lagrange multiplier for the generator set-point equality constraint x_i - x_0 
	     gradient is the partial derivative for it (note that it is negative) */
	  if(opflow->has_gensetpoint) {
	    grad[g++] = -lameq[gen->starteqloc+1];
	  }
	}
      }
    }
    ierr = VecRestoreArrayRead(Lambda,&lam);

    ierr = SCOPFLOWDestroy(&scopflowscen);CHKERRQ(ierr);
  }

  return true;
}

size_t SOPFLOWHIOPInterface::get_num_rterms() const
{
  return sopflow->Ns-1;
}

size_t SOPFLOWHIOPInterface::get_num_vars() const
{
  if(!sopflow->ismulticontingency) return sopflow->opflow0->nx;
  else return sopflow->scopflows[0]->opflows[0]->nx;
}

void SOPFLOWHIOPInterface::get_solution(double* x) const
{
  OPFLOW opflow;
  Vec    X;
  const double *xarr;
  PetscInt i;

  if(sopflow->sstart == 0) {
    if(!sopflow->ismulticontingency) {
      opflow = sopflow->opflow0;
      OPFLOWGetSolution(opflow,&X);
      VecGetArrayRead(X,&xarr);
      for(i=0; i < opflow->nx; i++) x[i] = xarr[i];
      VecRestoreArrayRead(X,&xarr);
    } else {
      SCOPFLOW scopflow=sopflow->scopflows[0];
      SCOPFLOWGetSolution(scopflow,0,&X);
      VecGetArrayRead(X,&xarr);
      for(i=0; i < scopflow->opflows[0]->nx; i++) x[i] = xarr[i];
      VecRestoreArrayRead(X,&xarr);
    }
  }
}
 
double SOPFLOWHIOPInterface::get_objective()
{
  OPFLOW opflow;
  double obj;
  PetscErrorCode ierr;

  if(sopflow->sstart == 0) {
    if(!sopflow->ismulticontingency) {
      opflow = sopflow->opflow0;
      ierr = OPFLOWGetObjective(opflow,&obj);CHKERRQ(ierr);
    } else {
      ierr = SCOPFLOWGetObjective(sopflow->scopflows[0],&obj);CHKERRQ(ierr);
    }
    return obj;
  }
}

PetscErrorCode SOPFLOWSolverSolve_HIOP(SOPFLOW sopflow)
{
  PetscErrorCode      ierr;
  PetscInt            c;
  SOPFLOWSolver_HIOP hiop = (SOPFLOWSolver_HIOP)sopflow->solver;

  PetscFunctionBegin;

  hiop->status = hiop->pridecsolver->run();

  /* Reset callbacks */
  if(!sopflow->ismulticontingency) {
    ierr = OPFLOWSetAuxillaryObjective(sopflow->opflow0,NULL,NULL,NULL,sopflow);CHKERRQ(ierr);
  } else {
    ierr = SCOPFLOWSetAuxillaryObjective(sopflow->scopflows[0],NULL,NULL,NULL,sopflow);CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSolverDestroy_HIOP(SOPFLOW sopflow)
{
  PetscErrorCode     ierr;
  SOPFLOWSolver_HIOP hiop = (SOPFLOWSolver_HIOP)sopflow->solver;

  PetscFunctionBegin;

  delete hiop->pridecompprob;
  delete hiop->pridecsolver;

  ierr = PetscFree(hiop);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSolverGetObjective_HIOP(SOPFLOW sopflow,PetscReal *obj)
{
  PetscErrorCode ierr;
  PetscReal      temp=0.0;
  PetscFunctionBegin;
  if(!sopflow->comm->rank) {
    if(!sopflow->ismulticontingency) {
      for(int i=0; i < sopflow->ns; i++) {
	temp += sopflow->opflows[i]->obj;
      }
    } else {
      for(int i=0; i < sopflow->ns; i++) {
	temp += sopflow->scopflows[i]->obj;
      }
    }
  }
  // ierr = MPI_Bcast(&temp,1,MPI_REAL,0,sopflow->comm->type);CHKERRQ(ierr);
  ierr = MPI_Bcast(&temp,1,MPI_DOUBLE,0,sopflow->comm->type);CHKERRQ(ierr);
  sopflow->obj = temp;
  *obj = sopflow->obj;
  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSolverGetSolution_HIOP(SOPFLOW sopflow,PetscInt scen_num,Vec *X)
{
  OPFLOW         opflow,opflow0;
  PS             ps,ps0;
  PSBUS          bus,bus0;
  PSGEN          gen,gen0;
  PetscInt       i,k;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  if(sopflow->sstart <= scen_num && scen_num < sopflow->send) {
    if(!sopflow->ismulticontingency) {
      opflow = sopflow->opflows[scen_num-sopflow->sstart];
      opflow0 = sopflow->opflow0;

      /* Update generator set-points */
      ps = opflow->ps;
      ps0 = opflow0->ps; 
      for(i=0; i < ps->nbus; i++) {
	bus = &ps->bus[i];
	bus0 = &ps0->bus[i];
	for(k=0; k < bus->ngen; k++) {
	  ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
	  ierr = PSBUSGetGen(bus0,k,&gen0);CHKERRQ(ierr);
	  if(gen0->status) {
	    gen->pgs = gen0->pgs;
	  }
	}
      }
      
      /* Solve */
      ierr = OPFLOWSolve(opflow);
      *X = opflow->X; 
    } else {
      SCOPFLOW scopflow = sopflow->scopflows[scen_num-sopflow->sstart];
      OPFLOW   opflow0 = sopflow->scopflows[0]->opflow0;
      OPFLOW   opflow = scopflow->opflows[0];

      /* Update generator set-points */
      ps = opflow->ps;
      ps0 = opflow0->ps; 
      for(i=0; i < ps->nbus; i++) {
	bus = &ps->bus[i];
	bus0 = &ps0->bus[i];
	for(k=0; k < bus->ngen; k++) {
	  ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
	  ierr = PSBUSGetGen(bus0,k,&gen0);CHKERRQ(ierr);
	  if(gen0->status) {
	    gen->pgs = gen0->pgs;
	  }
	}
      }
      
      /* Solve */
      ierr = SCOPFLOWSolve(scopflow);
      *X = scopflow->X;
    }
  } else {
    *X = NULL;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSolverGetConstraints_HIOP(SOPFLOW sopflow,PetscInt scen_num,Vec *G)
{
  OPFLOW              opflow;

  PetscFunctionBegin;

  if(sopflow->sstart <= scen_num && scen_num < sopflow->send) {
    if(!sopflow->ismulticontingency) {
      opflow = sopflow->opflows[scen_num-sopflow->sstart];
      *G = opflow->G;
    } 
   } else {
    *G = NULL;
  }

  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSolverGetConstraintMultipliers_HIOP(SOPFLOW sopflow,PetscInt scen_num,Vec *Lambda)
{
  OPFLOW              opflow;

  PetscFunctionBegin;
  if(sopflow->sstart <= scen_num && scen_num < sopflow->send) {
    if(!sopflow->ismulticontingency) {
      opflow = sopflow->opflows[scen_num-sopflow->sstart];
      *Lambda = opflow->Lambda; 
    }
  } else {
    *Lambda = NULL;
  }

  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSolverGetConvergenceStatus_HIOP(SOPFLOW sopflow,PetscBool *status)
{
  SOPFLOWSolver_HIOP hiop=(SOPFLOWSolver_HIOP)sopflow->solver;
  PetscFunctionBegin;

  if(hiop->status == hiop::Solve_Success) *status = PETSC_TRUE;
  else *status = PETSC_FALSE;

  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSolverSetUp_HIOP(SOPFLOW sopflow)
{
  SOPFLOWSolver_HIOP hiop=(SOPFLOWSolver_HIOP)sopflow->solver;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  hiop->pridecompprob = new SOPFLOWHIOPInterface(sopflow);

  hiop->pridecsolver = new hiop::hiopAlgPrimalDecomposition(hiop->pridecompprob,hiop->pridecompprob->nxcoup,hiop->pridecompprob->loc_xcoup,sopflow->comm->type);

  /* Add auxillary functions to base case */
  if(sopflow->sstart == 0) {
    if(!sopflow->ismulticontingency) {
      ierr = OPFLOWSetAuxillaryObjective(sopflow->opflow0,SOPFLOWBaseAuxObjectiveFunction,SOPFLOWBaseAuxGradientFunction,SOPFLOWBaseAuxHessianFunction,sopflow);CHKERRQ(ierr);
    } else {
      ierr = SCOPFLOWSetAuxillaryObjective(sopflow->scopflows[0],SOPFLOWBaseAuxSCOPFLOWObjectiveFunction,SOPFLOWBaseAuxSCOPFLOWGradientFunction,SOPFLOWBaseAuxSCOPFLOWHessianFunction,sopflow);CHKERRQ(ierr);
    }
  }

  PetscFunctionReturn(0);
}


PetscErrorCode SOPFLOWSolverCreate_HIOP(SOPFLOW sopflow)
{
  PetscErrorCode ierr;
  SOPFLOWSolver_HIOP hiop;
  
  PetscFunctionBegin;

  ierr = PetscCalloc1(1,&hiop);CHKERRQ(ierr);

  sopflow->solver = hiop;

  sopflow->solverops.setup = SOPFLOWSolverSetUp_HIOP;
  sopflow->solverops.solve = SOPFLOWSolverSolve_HIOP;
  sopflow->solverops.destroy = SOPFLOWSolverDestroy_HIOP;
  sopflow->solverops.getobjective = SOPFLOWSolverGetObjective_HIOP;
  sopflow->solverops.getsolution  = SOPFLOWSolverGetSolution_HIOP;
  sopflow->solverops.getconvergencestatus = SOPFLOWSolverGetConvergenceStatus_HIOP;
  sopflow->solverops.getconstraints = SOPFLOWSolverGetConstraints_HIOP;
  sopflow->solverops.getconstraintmultipliers = SOPFLOWSolverGetConstraintMultipliers_HIOP;

  PetscFunctionReturn(0);
}

#endif
#endif
