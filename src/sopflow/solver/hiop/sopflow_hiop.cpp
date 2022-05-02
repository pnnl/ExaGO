#include <exago_config.h>

#if defined(EXAGO_ENABLE_HIOP)

#include "sopflow_hiop.h"
#include <private/opflowimpl.h>
#include <private/sopflowimpl.h>

PetscErrorCode SOPFLOWBaseAuxObjectiveFunction(OPFLOW opflow, const double *x,
                                               double *obj, void *ctx) {
  PetscErrorCode ierr;
  SOPFLOW sopflow = (SOPFLOW)ctx;
  SOPFLOWSolver_HIOP hiop = (SOPFLOWSolver_HIOP)sopflow->solver;

  PetscFunctionBegin;

  if (hiop->pridecompprob->include_r_) {
    hiop->pridecompprob->rec_evaluator->eval_f(opflow->nx, x, false, *obj);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWBaseAuxGradientFunction(OPFLOW opflow, const double *x,
                                              double *grad, void *ctx) {
  PetscErrorCode ierr;
  SOPFLOW sopflow = (SOPFLOW)ctx;
  SOPFLOWSolver_HIOP hiop = (SOPFLOWSolver_HIOP)sopflow->solver;

  PetscFunctionBegin;

  if (hiop->pridecompprob->include_r_) {
    hiop->pridecompprob->rec_evaluator->eval_grad(opflow->nx, x, false, grad);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWBaseAuxHessianFunction(OPFLOW opflow, const double *x,
                                             Mat Hess, void *ctx) {
  PetscErrorCode ierr;
  SOPFLOW sopflow = (SOPFLOW)ctx;
  SOPFLOWSolver_HIOP hiop = (SOPFLOWSolver_HIOP)sopflow->solver;
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

SOPFLOWHIOPInterface::~SOPFLOWHIOPInterface() {
  PetscFree(loc_xcoup);
  //    printf("Exiting application\n");
}

SOPFLOWHIOPInterface::SOPFLOWHIOPInterface(SOPFLOW sopflowin) {
  int i, k, j = 0;
  PS ps0;
  PSBUS bus0;
  PSGEN gen0;
  OPFLOW opflowbase;

  sopflow = sopflowin;

  rec_evaluator = NULL;

  opflowbase = sopflow->opflow0;
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

hiop::hiopSolveStatus SOPFLOWHIOPInterface::solve_master(
    hiop::hiopVector &xvec, const bool &include_r, const double &rval,
    const double *grad, const double *hess, const char *master_options_file) {

  PetscErrorCode ierr;
  double *x = xvec.local_data();
  double obj;
  OPFLOW opflow;
  Vec X;
  const PetscScalar *xsol;
  include_r_ = include_r;

  opflow = sopflow->opflow0;

  //  printf("Rank[%d]:Enter solve_master\n",sopflow->comm->rank);
  ierr = OPFLOWSolve(opflow);
  ierr = OPFLOWSolutionToPS(opflow);
  ierr = OPFLOWGetObjective(opflow, &obj);
  ierr = OPFLOWGetSolution(opflow, &X);
  //  ierr = OPFLOWPrintSolution(opflow);

  ierr = VecGetArrayRead(X, &xsol);
  ierr = PetscMemcpy(x, xsol, opflow->nx * sizeof(double));
  ierr = VecRestoreArrayRead(X, &xsol);

  //  printf("Rank[%d]:Exit solve_master\n",sopflow->comm->rank);

  return hiop::Solve_Success;
}

bool SOPFLOWHIOPInterface::set_recourse_approx_evaluator(
    const int n,
    hiopInterfacePriDecProblem::RecourseApproxEvaluator *evaluator) {
  assert(n == nxcoup);
  rec_evaluator = evaluator;

  return true;
}

extern PetscErrorCode SOPFLOWUpdateOPFLOWVariableBounds(OPFLOW, Vec, Vec,
                                                        void *);

/* Note: x only holds the coupled variables which, in this case, are the
   generator real power variables for the base-case
*/
bool SOPFLOWHIOPInterface::eval_f_rterm(size_t idx, const int &n,
                                        const double *x, double &rval) {
  PetscErrorCode ierr;
  OPFLOW opflow0; /* base case OPFLOW */
  PS ps, ps0;
  PSBUS bus, bus0;
  PSGEN gen, gen0;
  PetscInt i, k, j, g = 0;
  PetscInt s = idx + 1;
  PetscInt scen_num, cont_num = 0;
  PetscBool issubproblemsolver_hiop = PETSC_FALSE;

  //  printf("[Rank %d] contingency %d Came in recourse objective
  //  function\n",sopflow->comm->rank,s);

  opflow0 = sopflow->opflow0;

  ierr = OPFLOWCreate(PETSC_COMM_SELF, &opflowscen);
  CHKERRQ(ierr);
  ierr = OPFLOWSetModel(opflowscen, sopflow->subproblem_model);
  CHKERRQ(ierr);
  ierr = OPFLOWSetInitializationType(opflowscen, sopflow->initialization_type);
  CHKERRQ(ierr);
  ierr = OPFLOWSetSolver(opflowscen, sopflow->subproblem_solver);
  CHKERRQ(ierr);
  ierr =
      PetscStrcmp(sopflow->subproblem_solver, "HIOP", &issubproblemsolver_hiop);
  CHKERRQ(ierr);
  if (issubproblemsolver_hiop) {
    ierr = OPFLOWSetHIOPComputeMode(opflowscen, sopflow->compute_mode);
    CHKERRQ(ierr);
    ierr = OPFLOWSetHIOPVerbosityLevel(opflowscen, sopflow->verbosity_level);
    CHKERRQ(ierr);
  }

  ierr = OPFLOWReadMatPowerData(opflowscen, sopflow->netfile);
  CHKERRQ(ierr);
  /* Set up the PS object for opflow */
  ps = opflowscen->ps;
  ierr = PSSetUp(ps);
  CHKERRQ(ierr);

  scen_num = s;
  cont_num = 0;
  if (!sopflow->flatten_contingencies) {
    scen_num = s;
    if (sopflow->scenfileset) {
      ierr = PSApplyScenario(ps, sopflow->scenlist.scen[s]);
      CHKERRQ(ierr);
    }
  } else {
    if (sopflow->scenfileset) {
      scen_num = s / sopflow->Nc;
      cont_num = s % sopflow->Nc;
      ierr = PSApplyScenario(ps, sopflow->scenlist.scen[scen_num]);
      CHKERRQ(ierr);

      if (cont_num != 0) {
        ierr = PSApplyContingency(ps, sopflow->ctgclist->cont[cont_num]);
        CHKERRQ(ierr);
        //      ierr = OPFLOWSetObjectiveType(opflowscen,NO_OBJ);CHKERRQ(ierr);
      }
    }
  }

  ierr = OPFLOWSetWeight(opflowscen, sopflow->scenlist.scen[scen_num].prob);
  CHKERRQ(ierr);

  ierr = OPFLOWHasGenSetPoint(opflowscen, PETSC_TRUE);
  CHKERRQ(ierr); /* Activates ramping variables */

  /* Update generator set-points */
  ps = opflowscen->ps;
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
        gen->pgs = x[g++];
      }
    }
  }

  ierr = OPFLOWSetUpdateVariableBoundsFunction(
      opflowscen, SOPFLOWUpdateOPFLOWVariableBounds, (void *)sopflow);

  ierr = OPFLOWSetUp(opflowscen);
  CHKERRQ(ierr);

  //  assert(g == n);
  /* Solve */
  /*
  ierr = PetscPrintf(PETSC_COMM_SELF,
                     "\n*******************************************\n",
                     scen_num + 1);
  CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF,
                     "SC = %d, SCENARIO NUMBER = %d, CONTINGENCY NUMBER = %d\n",
                     s, scen_num + 1, cont_num);
  CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF,
                     "\n*******************************************\n",
                     scen_num + 1);
  CHKERRQ(ierr);
  */

  ierr = OPFLOWSolve(opflowscen);
  ierr = OPFLOWGetObjective(opflowscen, &rval);
  ierr = OPFLOWSolutionToPS(opflowscen);
  CHKERRQ(ierr);
  //  ierr = OPFLOWPrintSolution(opflowscen);
  return true;
}

bool SOPFLOWHIOPInterface::eval_grad_rterm(size_t idx, const int &n, double *x,
                                           hiop::hiopVector &gradvec) {
  PetscErrorCode ierr;
  double *grad = gradvec.local_data();
  OPFLOW opflow;
  OPFLOW opflow0; /* base case OPFLOW */
  PS ps, ps0;
  PSBUS bus, bus0;
  PSGEN gen, gen0;
  PetscInt i, k, j, g = 0;
  const PetscScalar *lam, *lameq;
  Vec Lambda;
  PetscInt s = idx + 1;

  // printf("[Rank %d] contingency %d Came in recourse gradient
  // function\n",sopflow->comm->rank,scen_num);
  opflow0 = sopflow->opflow0;
  opflow = opflowscen;

  ierr = OPFLOWGetConstraintMultipliers(opflow, &Lambda);
  ierr = VecGetArrayRead(Lambda, &lam);
  lameq = lam;

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

        /*
          ierr = PetscPrintf(PETSC_COMM_SELF,"Gen[%d]: Pg = %lf Pb = %lf Pt =
          %lf Pgs = %lf Ramp rate = %lf lambda = %lf\n",gen->bus_i,gen->pg,
          gen->pb,gen->pt,gen->pgs,gen->ramp_rate_30min,lameq[gen->starteqloc+1]);CHKERRQ(ierr);
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

size_t SOPFLOWHIOPInterface::get_num_rterms() const {
  return (sopflow->Ns * sopflow->Nc) - 1;
}

size_t SOPFLOWHIOPInterface::get_num_vars() const {
  return sopflow->opflow0->nx;
}

void SOPFLOWHIOPInterface::get_solution(double *x) const {
  OPFLOW opflow;
  Vec X;
  const double *xarr;
  PetscInt i;

  if (sopflow->sstart == 0) {
    opflow = sopflow->opflow0;
    OPFLOWGetSolution(opflow, &X);
    VecGetArrayRead(X, &xarr);
    for (i = 0; i < opflow->nx; i++)
      x[i] = xarr[i];
    VecRestoreArrayRead(X, &xarr);
  }
}

double SOPFLOWHIOPInterface::get_objective() {
  OPFLOW opflow;
  double obj;
  PetscErrorCode ierr;

  if (sopflow->sstart == 0) {
    opflow = sopflow->opflow0;
    ierr = OPFLOWGetObjective(opflow, &obj);
    CHKERRQ(ierr);
    return obj;
  }
}

PetscErrorCode SOPFLOWSolverSolve_HIOP(SOPFLOW sopflow) {
  PetscErrorCode ierr;
  PetscInt c;
  SOPFLOWSolver_HIOP hiop = (SOPFLOWSolver_HIOP)sopflow->solver;

  PetscFunctionBegin;

  hiop->status = hiop->pridecsolver->run();

  /* Reset callbacks */
  ierr =
      OPFLOWSetAuxillaryObjective(sopflow->opflow0, NULL, NULL, NULL, sopflow);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSolverDestroy_HIOP(SOPFLOW sopflow) {
  PetscErrorCode ierr;
  SOPFLOWSolver_HIOP hiop = (SOPFLOWSolver_HIOP)sopflow->solver;

  PetscFunctionBegin;

  delete hiop->pridecompprob;
  delete hiop->pridecsolver;

  ierr = PetscFree(hiop);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSolverGetBaseObjective_HIOP(SOPFLOW sopflow,
                                                  PetscReal *obj) {
  PetscErrorCode ierr;
  PetscReal temp = 0.0;
  PetscFunctionBegin;
  if (!sopflow->comm->rank) {
    if (!sopflow->ismulticontingency) {
      temp = sopflow->opflows[0]->obj;
    }
  }
  // ierr = MPI_Bcast(&temp,1,MPI_REAL,0,sopflow->comm->type);CHKERRQ(ierr);
  ierr = MPI_Bcast(&temp, 1, MPI_DOUBLE, 0, sopflow->comm->type);
  CHKERRQ(ierr);
  sopflow->objbase = temp;
  *obj = sopflow->objbase;
  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSolverGetTotalObjective_HIOP(SOPFLOW sopflow,
                                                   PetscReal *obj) {
  PetscErrorCode ierr;
  PetscReal temp = 0.0;
  PetscFunctionBegin;
  if (!sopflow->ismulticontingency) {
    for (int i = 0; i < sopflow->ns; i++) {
      temp += sopflow->opflows[i]->obj;
    }
  } else {
    //    temp = sopflow->scopflows[0]->obj;
  }

  // ierr = MPI_Bcast(&temp,1,MPI_REAL,0,sopflow->comm->type);CHKERRQ(ierr);
  ierr = MPI_Allreduce(&temp, &sopflow->objtot, 1, MPI_DOUBLE, MPI_SUM,
                       sopflow->comm->type);
  CHKERRQ(ierr);
  *obj = sopflow->objtot;
  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSolverGetSolution_HIOP(SOPFLOW sopflow, PetscInt scen_num,
                                             Vec *X) {
  OPFLOW opflow, opflow0;
  PS ps, ps0;
  PSBUS bus, bus0;
  PSGEN gen, gen0;
  PetscInt i, k;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  if (sopflow->sstart <= scen_num && scen_num < sopflow->send) {
    if (!sopflow->ismulticontingency) {
      opflow = sopflow->opflows[scen_num - sopflow->sstart];
      opflow0 = sopflow->opflow0;

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
      CHKERRQ(ierr);
      *X = opflow->X;
    }
  } else {
    *X = NULL;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSolverGetConstraints_HIOP(SOPFLOW sopflow,
                                                PetscInt scen_num, Vec *G) {
  OPFLOW opflow;

  PetscFunctionBegin;

  if (sopflow->sstart <= scen_num && scen_num < sopflow->send) {
    if (!sopflow->ismulticontingency) {
      opflow = sopflow->opflows[scen_num - sopflow->sstart];
      *G = opflow->G;
    }
  } else {
    *G = NULL;
  }

  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSolverGetConstraintMultipliers_HIOP(SOPFLOW sopflow,
                                                          PetscInt scen_num,
                                                          Vec *Lambda) {
  OPFLOW opflow;

  PetscFunctionBegin;
  if (sopflow->sstart <= scen_num && scen_num < sopflow->send) {
    if (!sopflow->ismulticontingency) {
      opflow = sopflow->opflows[scen_num - sopflow->sstart];
      *Lambda = opflow->Lambda;
    }
  } else {
    *Lambda = NULL;
  }

  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSolverGetConvergenceStatus_HIOP(SOPFLOW sopflow,
                                                      PetscBool *status) {
  SOPFLOWSolver_HIOP hiop = (SOPFLOWSolver_HIOP)sopflow->solver;
  PetscFunctionBegin;

  if (hiop->status == hiop::Solve_Success)
    *status = PETSC_TRUE;
  else
    *status = PETSC_FALSE;

  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSolverSetUp_HIOP(SOPFLOW sopflow) {
  SOPFLOWSolver_HIOP hiop = (SOPFLOWSolver_HIOP)sopflow->solver;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  hiop->pridecompprob = new SOPFLOWHIOPInterface(sopflow);

  hiop->pridecsolver = new hiop::hiopAlgPrimalDecomposition(
      hiop->pridecompprob, hiop->pridecompprob->nxcoup,
      hiop->pridecompprob->loc_xcoup, sopflow->comm->type);

  /* Add auxillary functions to base case */
  if (sopflow->sstart == 0) {
    ierr = OPFLOWSetAuxillaryObjective(
        sopflow->opflow0, SOPFLOWBaseAuxObjectiveFunction,
        SOPFLOWBaseAuxGradientFunction, SOPFLOWBaseAuxHessianFunction, sopflow);
    CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSolverCreate_HIOP(SOPFLOW sopflow) {
  PetscErrorCode ierr;
  SOPFLOWSolver_HIOP hiop;

  PetscFunctionBegin;

  ierr = PetscCalloc1(1, &hiop);
  CHKERRQ(ierr);

  sopflow->solver = hiop;

  sopflow->solverops.setup = SOPFLOWSolverSetUp_HIOP;
  sopflow->solverops.solve = SOPFLOWSolverSolve_HIOP;
  sopflow->solverops.destroy = SOPFLOWSolverDestroy_HIOP;
  sopflow->solverops.getbaseobjective = SOPFLOWSolverGetBaseObjective_HIOP;
  sopflow->solverops.gettotalobjective = SOPFLOWSolverGetTotalObjective_HIOP;
  sopflow->solverops.getsolution = SOPFLOWSolverGetSolution_HIOP;
  sopflow->solverops.getconvergencestatus =
      SOPFLOWSolverGetConvergenceStatus_HIOP;
  sopflow->solverops.getconstraints = SOPFLOWSolverGetConstraints_HIOP;
  sopflow->solverops.getconstraintmultipliers =
      SOPFLOWSolverGetConstraintMultipliers_HIOP;

  PetscFunctionReturn(0);
}

#endif
