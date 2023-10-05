#include <exago_config.h>

#if defined(EXAGO_ENABLE_HIOP)

#include "scopflow_hiop.h"
#include <private/opflowimpl.h>
#include <private/scopflowimpl.h>
#include <private/tcopflowimpl.h>

extern const char *HIOPMemSpaceChoices[];

PetscErrorCode SCOPFLOWBaseAuxObjectiveFunction(OPFLOW opflow, const double *x,
                                                double *obj, void *ctx) {
  SCOPFLOW scopflow = (SCOPFLOW)ctx;
  SCOPFLOWSolver_HIOP hiop = (SCOPFLOWSolver_HIOP)scopflow->solver;

  PetscFunctionBegin;

  if (hiop->pridecompprob->include_r_) {
    hiop->pridecompprob->rec_evaluator->eval_f(opflow->nx, x, false, *obj);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWBaseAuxGradientFunction(OPFLOW opflow, const double *x,
                                               double *grad, void *ctx) {
  SCOPFLOW scopflow = (SCOPFLOW)ctx;
  SCOPFLOWSolver_HIOP hiop = (SCOPFLOWSolver_HIOP)scopflow->solver;

  PetscFunctionBegin;

  if (hiop->pridecompprob->include_r_) {
    hiop->pridecompprob->rec_evaluator->eval_grad(opflow->nx, x, false, grad);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWBaseAuxHessianFunction(OPFLOW opflow, const double *x,
                                              Mat Hess, void *ctx) {
  (void)x;
  PetscErrorCode ierr;
  SCOPFLOW scopflow = (SCOPFLOW)ctx;
  SCOPFLOWSolver_HIOP hiop = (SCOPFLOWSolver_HIOP)scopflow->solver;
  PetscInt row, col;
  PetscScalar val;
  PetscInt i;

  PetscFunctionBegin;
  if (hiop->pridecompprob->include_r_) {
    for (i = 0; i < hiop->pridecompprob->nxcoup; i++) {
      row = col = hiop->pridecompprob->loc_xcoup[i];
      val = hiop->pridecompprob->rec_evaluator->get_rhess()
                ->local_data_const()[i];
      val *= opflow->obj_factor;
      ierr = MatSetValues(Hess, 1, &row, 1, &col, &val, ADD_VALUES);
      CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

SCOPFLOWHIOPInterface::~SCOPFLOWHIOPInterface() {
  //    printf("Exiting application\n");
}

SCOPFLOWHIOPInterface::SCOPFLOWHIOPInterface(SCOPFLOW scopflowin) {
  int i, k, j = 0;
  PS ps0;
  PSBUS bus0;
  PSGEN gen0;
  OPFLOW opflowbase;

  scopflow = scopflowin;

  rec_evaluator = NULL;

  opflowbase = scopflow->opflow0;
  ps0 = opflowbase->ps;

  /* Get the number of coupling variables */
  nxcoup = 0;

  for (i = 0; i < ps0->nbus; i++) {
    bus0 = &ps0->bus[i];
    for (k = 0; k < bus0->ngen; k++) {
      PSBUSGetGen(bus0, k, &gen0);
      if (gen0->status && !gen0->isrenewable)
        nxcoup++;
    }
  }

  PetscCalloc1(nxcoup, &loc_xcoup);
  /* Set the coupling indices */
  for (i = 0; i < ps0->nbus; i++) {
    bus0 = &ps0->bus[i];
    for (k = 0; k < bus0->ngen; k++) {
      PSBUSGetGen(bus0, k, &gen0);
      if (gen0->status && !gen0->isrenewable) {
        loc_xcoup[j] = gen0->startxpowloc; /* Location for Pg */
        j++;
      }
    }
  }
  assert(j == nxcoup);
}

hiop::hiopSolveStatus SCOPFLOWHIOPInterface::solve_master(
    hiop::hiopVector &xvec, const bool &include_r, const double &rval,
    const double *grad, const double *hess, const char *master_options_file) {
  (void)rval;
  (void)grad;
  (void)hess;
  (void)master_options_file;
  double *x = xvec.local_data();
  double obj;
  PetscErrorCode ierr;
  OPFLOW opflow;
  Vec X;
  const PetscScalar *xsol;
  include_r_ = include_r;

  opflow = scopflow->opflow0;

  //  printf("Rank[%d]:Enter solve_master\n",scopflow->comm->rank);
  ierr = OPFLOWSolve(opflow);
  ExaGOCheckError(ierr);
  ierr = OPFLOWGetObjective(opflow, &obj);
  ExaGOCheckError(ierr);
  ierr = OPFLOWGetSolution(opflow, &X);
  ExaGOCheckError(ierr);
  ierr = VecGetArrayRead(X, &xsol);
  ExaGOCheckError(ierr);
  ierr = PetscMemcpy(x, xsol, opflow->nx * sizeof(double));
  ExaGOCheckError(ierr);
  ierr = VecRestoreArrayRead(X, &xsol);
  ExaGOCheckError(ierr);
  //  printf("Rank[%d]:Exit solve_master\n",scopflow->comm->rank);

  return hiop::Solve_Success;
}

bool SCOPFLOWHIOPInterface::set_recourse_approx_evaluator(
    const int n,
    hiopInterfacePriDecProblem::RecourseApproxEvaluator *evaluator) {
  assert(n == nxcoup);
  rec_evaluator = evaluator;

  return true;
}

extern PetscErrorCode SCOPFLOWUpdateOPFLOWVariableBounds(OPFLOW, Vec, Vec,
                                                         void *);

/* Note: x only holds the coupled variables which, in this case, are the
   generator real power variables for the base-case
*/
bool SCOPFLOWHIOPInterface::eval_f_rterm(hiop::size_type idx, const int &n,
                                         const double *x, double &rval) {
  PetscErrorCode ierr;
  OPFLOW opflow0; /* base case OPFLOW */
  PS ps, ps0;
  PSBUS bus, bus0;
  PSGEN gen, gen0;
  PetscInt i, k, g = 0;
  PetscInt cont_num = idx + 1;

  //  printf("[Rank %d] contingency %d Came in recourse objective
  //  function\n",scopflow->comm->rank,cont_num);

  opflow0 = scopflow->opflow0;

  ierr = OPFLOWCreate(PETSC_COMM_SELF, &opflowctgc);
  CHKERRQ(ierr);
  ierr = OPFLOWSetModel(opflowctgc, scopflow->subproblem_model);
  CHKERRQ(ierr);
  ierr = OPFLOWSetSolver(opflowctgc, scopflow->subproblem_solver);
  CHKERRQ(ierr);
  if (scopflow->subproblem_solver == "HIOP") {
    ierr = OPFLOWSetHIOPComputeMode(opflowctgc, scopflow->compute_mode);
    CHKERRQ(ierr);
    ierr = OPFLOWSetHIOPVerbosityLevel(opflowctgc, scopflow->verbosity_level);
    CHKERRQ(ierr);
  }
  ierr = OPFLOWSetInitializationType(opflowctgc, scopflow->type);
  CHKERRQ(ierr);
  ierr = OPFLOWSetGenBusVoltageType(opflowctgc, scopflow->genbusvoltagetype);
  CHKERRQ(ierr);
  ierr = OPFLOWHasBusPowerImbalance(opflowctgc,
                                    scopflow->enable_powerimbalance_variables);
  CHKERRQ(ierr);
  ierr = OPFLOWIgnoreLineflowConstraints(opflowctgc,
                                         scopflow->ignore_lineflow_constraints);
  CHKERRQ(ierr);

  ierr = OPFLOWReadMatPowerData(opflowctgc, scopflow->netfile);
  CHKERRQ(ierr);
  /* Set up the PS object for opflow */
  ps = opflowctgc->ps;
  ierr = PSSetUp(ps);
  CHKERRQ(ierr);

  /* Set contingencies */
  if (scopflow->ctgcfileset) {
    Contingency ctgc = scopflow->ctgclist->cont[cont_num];
    ierr = PSApplyContingency(ps, ctgc);
    CHKERRQ(ierr);
  }

  /* Update generator set-points */
  ps = opflowctgc->ps;
  ps0 = opflow0->ps;
  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    bus0 = &ps0->bus[i];
    for (k = 0; k < bus->ngen; k++) {
      ierr = PSBUSGetGen(bus, k, &gen);
      CHKERRQ(ierr);
      ierr = PSBUSGetGen(bus0, k, &gen0);
      CHKERRQ(ierr);
      if (gen0->status && !gen0->isrenewable) {
        gen0->pgs = gen->pgs = gen->pg = x[g++];
      }
    }
  }
  assert(g == n);

  ierr = OPFLOWHasGenSetPoint(opflowctgc, PETSC_TRUE);
  CHKERRQ(ierr); /* Activates ramping variables */
  ierr = OPFLOWSetObjectiveType(opflowctgc, NO_OBJ);
  CHKERRQ(ierr);
  ierr = OPFLOWSetUpdateVariableBoundsFunction(
      opflowctgc, SCOPFLOWUpdateOPFLOWVariableBounds, (void *)scopflow);

  ierr = OPFLOWSetUp(opflowctgc);
  CHKERRQ(ierr);

  /* Solve */
  ierr = OPFLOWSolve(opflowctgc);
  ierr = OPFLOWSolutionToPS(opflowctgc);
  CHKERRQ(ierr);
  ierr = OPFLOWGetObjective(opflowctgc, &rval);
  //  ierr = OPFLOWPrintSolution(opflowctgc);

  return true;
}

bool SCOPFLOWHIOPInterface::eval_grad_rterm(hiop::size_type idx, const int &n,
                                            double *x,
                                            hiop::hiopVector &gradvec) {
  (void)idx;
  (void)n;
  (void)x;

  PetscErrorCode ierr;
  double *grad = gradvec.local_data();
  OPFLOW opflow;
  OPFLOW opflow0; /* base case OPFLOW */
  PS ps, ps0;
  PSBUS bus, bus0;
  PSGEN gen, gen0;
  PetscInt i, k, g = 0;
  const PetscScalar *lam, *lameq;
  Vec Lambda;

  // PetscInt cont_num = idx + 1;
  // printf("[Rank %d] contingency %d Came in recourse gradient
  // function\n",scopflow->comm->rank,cont_num);

  opflow0 = scopflow->opflow0;
  opflow = opflowctgc;

  ierr = OPFLOWGetConstraintMultipliers(opflow, &Lambda);
  ierr = VecGetArrayRead(Lambda, &lam);
  lameq = lam;
  /* Update generator set-points */
  ps = opflow->ps;
  ps0 = opflow0->ps;
  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    bus0 = &ps0->bus[i];
    for (k = 0; k < bus->ngen; k++) {
      ierr = PSBUSGetGen(bus, k, &gen);
      CHKERRQ(ierr);
      ierr = PSBUSGetGen(bus0, k, &gen0);
      CHKERRQ(ierr);
      if (!gen0->status || gen0->isrenewable)
        continue;

      if (gen->status) {
        /* Get the lagrange multiplier for the generator set-point equality
           constraint x_i - x_0
           gradient is the partial derivative for it (note that it is negative)
         */

        /*	  ierr = PetscPrintf(PETSC_COMM_SELF,"Gen[%d]: Pg = %lf Pb = %lf
         * Pt = %lf Pgs = %lf Ramp rate = %lf lambda =
         * %lf\n",gen->bus_i,gen->pg,
         * gen->pb,gen->pt,gen->pgs,gen->ramp_rate_30min,lameq[gen->starteqloc+1]);CHKERRQ(ierr);
         */
        grad[g++] = -lameq[gen->starteqloc + 1];
      } else {
        grad[g++] = 0.0;
      }
    }
  }
  ierr = VecRestoreArrayRead(Lambda, &lam);

  ierr = OPFLOWDestroy(&opflow);
  CHKERRQ(ierr);
  return true;
}

hiop::size_type SCOPFLOWHIOPInterface::get_num_rterms() const {
  return (hiop::size_type)scopflow->Nc - 1;
}

hiop::size_type SCOPFLOWHIOPInterface::get_num_vars() const {
  return (hiop::size_type)scopflow->opflow0->nx;
}

void SCOPFLOWHIOPInterface::get_solution(double *x) const {
  OPFLOW opflow;
  Vec X;
  const double *xarr;
  PetscInt i;

  if (scopflow->cstart == 0) {
    opflow = scopflow->opflow0;
    OPFLOWGetSolution(opflow, &X);
    VecGetArrayRead(X, &xarr);
    for (i = 0; i < opflow->nx; i++)
      x[i] = xarr[i];
    VecRestoreArrayRead(X, &xarr);
  }
}

double SCOPFLOWHIOPInterface::get_objective() {
  OPFLOW opflow;
  double obj;
  PetscErrorCode ierr;

  if (scopflow->cstart == 0) {
    opflow = scopflow->opflow0;
    ierr = OPFLOWGetObjective(opflow, &obj);
    CHKERRQ(ierr);
    return obj;
  }

  throw ExaGOError("scopflow->cstart != 0 so no obj value to return");
}

PetscErrorCode SCOPFLOWSolverSolve_HIOP(SCOPFLOW scopflow) {
  PetscErrorCode ierr;
  SCOPFLOWSolver_HIOP hiop = (SCOPFLOWSolver_HIOP)scopflow->solver;

  PetscFunctionBegin;

  hiop->status = hiop->pridecsolver->run();

  /* Reset callbacks */
  ierr = OPFLOWSetAuxillaryObjective(scopflow->opflow0, NULL, NULL, NULL,
                                     scopflow);
  CHKERRQ(ierr);

  /* Get objective value */
  ierr = SCOPFLOWGetBaseObjective(scopflow, &scopflow->objbase);
  CHKERRQ(ierr);

  /* Get number of iterations */
  scopflow->numiter = hiop->pridecsolver->getNumIterations();

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverDestroy_HIOP(SCOPFLOW scopflow) {
  PetscErrorCode ierr;
  SCOPFLOWSolver_HIOP hiop = (SCOPFLOWSolver_HIOP)scopflow->solver;

  PetscFunctionBegin;

  ierr = PetscFree(hiop->pridecompprob->loc_xcoup);
  CHKERRQ(ierr);
  delete hiop->pridecompprob;
  delete hiop->pridecsolver;

  ierr = PetscFree(hiop);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverGetBaseObjective_HIOP(SCOPFLOW scopflow,
                                                   PetscReal *obj) {
  PetscErrorCode ierr;
  PetscReal temp = 0.0;
  PetscFunctionBegin;
  if (!scopflow->comm->rank) {
    // Only the objective for the base problem is returned
    if (!scopflow->ismultiperiod) {
      ierr = OPFLOWGetObjective(scopflow->opflow0, &temp);
    } else {
      temp = scopflow->tcopflows[0]->obj;
    }
  }
  // ierr = MPI_Bcast(&temp,1,MPI_REAL,0,scopflow->comm->type);CHKERRQ(ierr);
  ierr = MPI_Bcast(&temp, 1, MPI_DOUBLE, 0, scopflow->comm->type);
  CHKERRQ(ierr);
  scopflow->objbase = temp;
  *obj = scopflow->objbase;
  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverGetTotalObjective_HIOP(SCOPFLOW scopflow,
                                                    PetscReal *obj) {
  PetscErrorCode ierr;
  PetscReal temp;
  PetscFunctionBegin;
  if (!scopflow->ismultiperiod) {
    for (int i = 0; i < scopflow->nc; i++) {
      temp = scopflow->opflows[0]->obj;
    }
  } else {
    temp = scopflow->tcopflows[0]->obj;
  }

  ierr = MPI_Allreduce(&temp, &scopflow->objtot, 1, MPI_DOUBLE, MPI_SUM,
                       scopflow->comm->type);
  CHKERRQ(ierr);
  *obj = scopflow->objtot;
  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverGetSolution_HIOP(SCOPFLOW scopflow,
                                              PetscInt cont_num, Vec *X) {
  TCOPFLOW tcopflow;
  OPFLOW opflow, opflow0;
  PS ps, ps0;
  PSBUS bus, bus0;
  PSGEN gen, gen0;
  PetscInt i, k;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  if (scopflow->cstart <= cont_num && cont_num < scopflow->cend) {
    if (!scopflow->ismultiperiod) {
      opflow = scopflow->opflows[cont_num - scopflow->cstart];
      opflow0 = scopflow->opflow0;

      /* Update generator set-points */
      ps = opflow->ps;
      ps0 = opflow0->ps;
      for (i = 0; i < ps->nbus; i++) {
        bus = &ps->bus[i];
        bus0 = &ps0->bus[i];
        for (k = 0; k < bus->ngen; k++) {
          ierr = PSBUSGetGen(bus, k, &gen);
          CHKERRQ(ierr);
          ierr = PSBUSGetGen(bus0, k, &gen0);
          CHKERRQ(ierr);
          if (gen0->status && !gen0->isrenewable) {
            gen->pgs = gen0->pgs;
          }
        }
      }

      /* Solve */
      ierr = OPFLOWSolve(opflow);
      *X = opflow->X;
    } else {
      tcopflow = scopflow->tcopflows[cont_num - scopflow->cstart];
      *X = tcopflow->X;
    }
  } else {
    *X = NULL;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverGetConstraints_HIOP(SCOPFLOW scopflow,
                                                 PetscInt cont_num, Vec *G) {
  TCOPFLOW tcopflow;
  OPFLOW opflow;

  PetscFunctionBegin;

  if (scopflow->cstart <= cont_num && cont_num < scopflow->cend) {
    if (!scopflow->ismultiperiod) {
      opflow = scopflow->opflows[cont_num - scopflow->cstart];
      *G = opflow->G;
    } else {
      tcopflow = scopflow->tcopflows[cont_num - scopflow->cstart];
      *G = tcopflow->G;
    }
  } else {
    *G = NULL;
  }

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverGetConstraintMultipliers_HIOP(SCOPFLOW scopflow,
                                                           PetscInt cont_num,
                                                           Vec *Lambda) {
  TCOPFLOW tcopflow;
  OPFLOW opflow;

  PetscFunctionBegin;
  if (scopflow->cstart <= cont_num && cont_num < scopflow->cend) {
    if (!scopflow->ismultiperiod) {
      opflow = scopflow->opflows[cont_num - scopflow->cstart];
      *Lambda = opflow->Lambda;
    } else {
      tcopflow = scopflow->tcopflows[cont_num - scopflow->cstart];
      *Lambda = tcopflow->Lambda;
    }
  } else {
    *Lambda = NULL;
  }

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverGetConvergenceStatus_HIOP(SCOPFLOW scopflow,
                                                       PetscBool *status) {
  SCOPFLOWSolver_HIOP hiop = (SCOPFLOWSolver_HIOP)scopflow->solver;
  PetscFunctionBegin;

  if (hiop->status == hiop::Solve_Success)
    *status = PETSC_TRUE;
  else
    *status = PETSC_FALSE;

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverSetUp_HIOP(SCOPFLOW scopflow) {
  SCOPFLOWSolver_HIOP hiop = (SCOPFLOWSolver_HIOP)scopflow->solver;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  hiop->pridecompprob = new SCOPFLOWHIOPInterface(scopflow);

  hiop->pridecsolver = new hiop::hiopAlgPrimalDecomposition(
      hiop->pridecompprob, hiop->pridecompprob->nxcoup,
      hiop->pridecompprob->loc_xcoup, scopflow->comm->type);

  /* Add auxillary functions to base case */
  if (scopflow->cstart == 0) {
    ierr = OPFLOWSetAuxillaryObjective(
        scopflow->opflow0, SCOPFLOWBaseAuxObjectiveFunction,
        SCOPFLOWBaseAuxGradientFunction, SCOPFLOWBaseAuxHessianFunction,
        scopflow);
    CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSolverCreate_HIOP(SCOPFLOW scopflow) {
  PetscErrorCode ierr;
  SCOPFLOWSolver_HIOP hiop;

  PetscFunctionBegin;

  ierr = PetscCalloc1(1, &hiop);
  CHKERRQ(ierr);

  scopflow->solver = hiop;

  scopflow->solverops.setup = SCOPFLOWSolverSetUp_HIOP;
  scopflow->solverops.solve = SCOPFLOWSolverSolve_HIOP;
  scopflow->solverops.destroy = SCOPFLOWSolverDestroy_HIOP;
  scopflow->solverops.getbaseobjective = SCOPFLOWSolverGetBaseObjective_HIOP;
  scopflow->solverops.gettotalobjective = SCOPFLOWSolverGetTotalObjective_HIOP;
  scopflow->solverops.getsolution = SCOPFLOWSolverGetSolution_HIOP;
  scopflow->solverops.getconvergencestatus =
      SCOPFLOWSolverGetConvergenceStatus_HIOP;
  scopflow->solverops.getconstraints = SCOPFLOWSolverGetConstraints_HIOP;
  scopflow->solverops.getconstraintmultipliers =
      SCOPFLOWSolverGetConstraintMultipliers_HIOP;

  PetscFunctionReturn(0);
}

#endif
