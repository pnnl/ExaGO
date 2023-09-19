#include <exago_config.h>

#if defined(EXAGO_ENABLE_RAJA)
#if defined(EXAGO_ENABLE_HIOP_SPARSE)

#include <private/opflowimpl.h>
#include "pbpolrajahiopsparsekernels.hpp"

/* Initialization is done on the host through this function. Copying over values
 * to the device is done in OPFLOWSetInitialGuessArray_PBPOLRAJAHIOPSPARSE
 */
extern PetscErrorCode OPFLOWSetInitialGuess_PBPOL(OPFLOW, Vec, Vec);

PetscErrorCode OPFLOWSetInitialGuess_PBPOLRAJAHIOPSPARSE(OPFLOW opflow, Vec X,
                                                         Vec Lambda) {
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = OPFLOWSetInitialGuess_PBPOL(opflow, X, Lambda);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

extern PetscErrorCode OPFLOWSetConstraintBounds_PBPOL(OPFLOW, Vec, Vec);

/* The constraint bounds are also calculated on the host.
 */
PetscErrorCode OPFLOWSetConstraintBounds_PBPOLRAJAHIOPSPARSE(OPFLOW opflow,
                                                             Vec Gl, Vec Gu) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = OPFLOWSetConstraintBounds_PBPOL(opflow, Gl, Gu);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

extern PetscErrorCode OPFLOWSetVariableBounds_PBPOL(OPFLOW, Vec, Vec);

/* The variable bounds are also calculated on the host.
 */
PetscErrorCode OPFLOWSetVariableBounds_PBPOLRAJAHIOPSPARSE(OPFLOW opflow,
                                                           Vec Xl, Vec Xu) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = OPFLOWSetVariableBounds_PBPOL(opflow, Xl, Xu);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolutionToPS_PBPOLRAJAHIOPSPARSE(OPFLOW opflow) {
  PetscErrorCode ierr;
  PS ps = (PS)opflow->ps;
  PetscInt i, k;
  Vec X, Lambda;
  PSBUS bus;
  PSGEN gen;
  PSLOAD load;
  PSLINE line;
  const PetscScalar *x, *lambda, *lambdae, *lambdai;
  PetscInt loc, gloc = 0;
  PetscScalar Gff, Bff, Gft, Bft, Gtf, Btf, Gtt, Btt;
  PetscScalar Vmf, Vmt, thetaf, thetat, thetaft, thetatf;
  PetscScalar Pf, Qf, Pt, Qt;
  PSBUS busf, bust;
  const PSBUS *connbuses;
  PetscInt xlocf, xloct;

  PetscFunctionBegin;

  ierr = OPFLOWGetSolution(opflow, &X);
  CHKERRQ(ierr);
  ierr = OPFLOWGetConstraintMultipliers(opflow, &Lambda);
  CHKERRQ(ierr);

  ierr = VecGetArrayRead(X, &x);
  CHKERRQ(ierr);
  ierr = VecGetArrayRead(Lambda, &lambda);
  CHKERRQ(ierr);
  lambdae = lambda;
  if (opflow->Nconineq) {
    lambdai = lambdae + opflow->nconeq;
  }

  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    loc = bus->startxVloc;

    bus->va = x[loc];
    bus->vm = x[loc + 1];

    gloc = bus->starteqloc;
    bus->mult_pmis = lambdae[gloc];
    bus->mult_qmis = lambdae[gloc + 1];

    if (opflow->include_powerimbalance_variables) {
      loc = bus->startxpimbloc;
      bus->pimb = x[loc];
      bus->qimb = x[loc + 1];
    }

    for (k = 0; k < bus->ngen; k++) {
      ierr = PSBUSGetGen(bus, k, &gen);
      CHKERRQ(ierr);
      if (!gen->status) {
        gen->pg = gen->qg = 0.0;
        continue;
      }
      loc = gen->startxpowloc;

      gen->pg = x[loc];
      gen->qg = x[loc + 1];

      if (opflow->has_gensetpoint) {
        gloc += gen->nconeq;
      }
    }

    if (opflow->include_loadloss_variables) {
      for (k = 0; k < bus->nload; k++) {
        ierr = PSBUSGetLoad(bus, k, &load);
        CHKERRQ(ierr);
        loc = load->startxloadlossloc;
        load->pl = load->pl - x[loc];
        load->ql = load->ql - x[loc + 1];
      }
    }
  }

  if (!opflow->ignore_lineflow_constraints) {
    for (i = 0; i < ps->nline; i++) {
      line = &ps->line[i];
      if (!line->status) {
        line->mult_sf = line->mult_st = 0.0;
        continue;
      }

      Gff = line->yff[0];
      Bff = line->yff[1];
      Gft = line->yft[0];
      Bft = line->yft[1];
      Gtf = line->ytf[0];
      Btf = line->ytf[1];
      Gtt = line->ytt[0];
      Btt = line->ytt[1];

      ierr = PSLINEGetConnectedBuses(line, &connbuses);
      CHKERRQ(ierr);
      busf = connbuses[0];
      bust = connbuses[1];

      xlocf = busf->startxVloc;
      xloct = bust->startxVloc;

      thetaf = x[xlocf];
      Vmf = x[xlocf + 1];
      thetat = x[xloct];
      Vmt = x[xloct + 1];
      thetaft = thetaf - thetat;
      thetatf = thetat - thetaf;

      Pf = Gff * Vmf * Vmf +
           Vmf * Vmt * (Gft * cos(thetaft) + Bft * sin(thetaft));
      Qf = -Bff * Vmf * Vmf +
           Vmf * Vmt * (-Bft * cos(thetaft) + Gft * sin(thetaft));

      Pt = Gtt * Vmt * Vmt +
           Vmt * Vmf * (Gtf * cos(thetatf) + Btf * sin(thetatf));
      Qt = -Btt * Vmt * Vmt +
           Vmt * Vmf * (-Btf * cos(thetatf) + Gtf * sin(thetatf));

      line->pf = Pf;
      line->qf = Qf;
      line->pt = Pt;
      line->qt = Qt;
      line->sf = PetscSqrtScalar(Pf * Pf + Qf * Qf);
      line->st = PetscSqrtScalar(Pt * Pt + Qt * Qt);

      if (line->rateA > 1e5) {
        line->mult_sf = line->mult_st = 0.0;
      } else {
        gloc = line->startineqloc;
        line->mult_sf = lambdai[gloc];
        line->mult_st = lambdai[gloc + 1];
      }
    }
  }

  ierr = VecRestoreArrayRead(X, &x);
  CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(Lambda, &lambda);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/* Reuse PBPOL model set up for obtaining locations */
extern PetscErrorCode OPFLOWModelSetUp_PBPOL(OPFLOW);

PetscErrorCode OPFLOWModelSetUp_PBPOLRAJAHIOPSPARSE(OPFLOW opflow) {
  PetscErrorCode ierr;
  PbpolModelRajaHiop *pbpolrajahiopsparse =
      reinterpret_cast<PbpolModelRajaHiop *>(opflow->model);

  PetscFunctionBegin;

  ierr = OPFLOWModelSetUp_PBPOL(opflow);
  CHKERRQ(ierr);

  ierr = pbpolrajahiopsparse->busparams.allocate(opflow);
  ierr = pbpolrajahiopsparse->genparams.allocate(opflow);
  ierr = pbpolrajahiopsparse->lineparams.allocate(opflow);
  ierr = pbpolrajahiopsparse->loadparams.allocate(opflow);

  BUSParamsRajaHiop *busparams = &pbpolrajahiopsparse->busparams;
  GENParamsRajaHiop *genparams = &pbpolrajahiopsparse->genparams;
  LOADParamsRajaHiop *loadparams = &pbpolrajahiopsparse->loadparams;
  LINEParamsRajaHiop *lineparams = &pbpolrajahiopsparse->lineparams;

  /* Need to compute the number of nonzeros in equality, inequality constraint
   * Jacobians and Hessian */
  int nnz_eqjac = 0;

  // Find nonzero entries in equality constraint Jacobian by row. Using
  // OPFLOWComputeEqualityConstraintJacobian_PBPOL() as a guide.

  PS ps = (PS)opflow->ps;

  for (int ibus = 0; ibus < ps->nbus; ++ibus) {

    PSBUS bus = &(ps->bus[ibus]);

    // Nonzero entries used by each *bus* starts here 

    // no matter what, each bus uses 2 rows and 2 columns
    // row 1 = real, row2 = reactive

    busparams->jacsp_idx[ibus] = nnz_eqjac;
    nnz_eqjac += 2;
    busparams->jacsq_idx[ibus] = nnz_eqjac;
    nnz_eqjac += 2;
    
    if (bus->ide == ISOLATED_BUS) {
      continue;
    }
      
    if (opflow->include_powerimbalance_variables) {
      // 2 more entries on both real and reactive
      nnz_eqjac += 4;
    }

    for (int igen = 0; igen < bus->ngen; ++igen) {
      PSGEN gen;
      ierr = PSBUSGetGen(bus, igen, &gen);
      CHKERRQ(ierr);
      if (!gen->status)
        continue;
      // each active generator uses 1 real and reactive entry on each bus
      genparams->eqjacspbus_idx[igen] = nnz_eqjac;
      nnz_eqjac += 2;
    }
      
    if (opflow->include_loadloss_variables) {
      // each load adds one real and reactive entry on each bus row
      for (int iload = 0; iload < bus->nload; ++iload) {
        loadparams->jacsp_idx[iload] = nnz_eqjac;
        nnz_eqjac += 2;
      }
    }

    const PSLINE *connlines;
    int nconnlines;
    ierr = PSBUSGetSupportingLines(bus, &nconnlines, &connlines);
    CHKERRQ(ierr);
    
    for (int iconn = 0; iconn < nconnlines; iconn++) {
      // each *active* connected line uses 4 entries total in each bus row
      PSLINE line = connlines[iconn];
      if (!line->status)
        continue;

      // each line adds 4 entries for the to bus and 4 entries for the
      // from bus. The current bus is one of these and those entries
      // have already been counted.
      nnz_eqjac += 4;
    }
      
    if (opflow->has_gensetpoint) {
      for (int igen = 0; igen < bus->ngen; ++igen) {
        PSGEN gen;
        ierr = PSBUSGetGen(bus, igen, &gen);
        CHKERRQ(ierr);
        
        if (!gen->status || gen->isrenewable)
          continue;

        // each generator uses 2 rows, 3 columns real, 1 column reactive
        genparams->eqjacspbus_idx[igen] = nnz_eqjac;
        nnz_eqjac += 4;
      }
    }
  }

  printf("Equality Jacobian nonzero count: %d vs %d\n",
         opflow->nnz_eqjacsp, nnz_eqjac);

  ierr = busparams->copy(opflow);
  ierr = genparams->copy(opflow);
  ierr = lineparams->copy(opflow);
  ierr = loadparams->copy(opflow);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWModelDestroy_PBPOLRAJAHIOPSPARSE(OPFLOW opflow) {
  PbpolModelRajaHiop *pbpolrajahiopsparse =
      reinterpret_cast<PbpolModelRajaHiop *>(opflow->model);

  PetscFunctionBegin;
  pbpolrajahiopsparse->destroy(opflow);
  delete pbpolrajahiopsparse;
  pbpolrajahiopsparse = nullptr;

  PetscFunctionReturn(0);
}

/* reuse numvariables and numconstraints functions from PBPOL model */
extern PetscErrorCode OPFLOWModelSetNumVariables_PBPOL(OPFLOW, PetscInt *,
                                                       PetscInt *, PetscInt *);
extern PetscErrorCode OPFLOWModelSetNumConstraints_PBPOL(OPFLOW, PetscInt *,
                                                         PetscInt *, PetscInt *,
                                                         PetscInt *);
extern PetscErrorCode OPFLOWComputeEqualityConstraintJacobian_PBPOL(OPFLOW, Vec,
                                                                    Mat);
extern PetscErrorCode OPFLOWComputeInequalityConstraintJacobian_PBPOL(OPFLOW,
                                                                      Vec, Mat);
extern PetscErrorCode OPFLOWComputeHessian_PBPOL(OPFLOW, Vec, Vec, Vec, Mat);

extern PetscErrorCode OPFLOWSolutionCallback_PBPOLRAJAHIOPSPARSE(
    OPFLOW, const double *, const double *, const double *, const double *,
    const double *, double);

PetscErrorCode OPFLOWModelCreate_PBPOLRAJAHIOPSPARSE(OPFLOW opflow) {

  PetscFunctionBegin;

  PbpolModelRajaHiop *pbpolrajahiopsparse = new PbpolModelRajaHiop();

  opflow->model = pbpolrajahiopsparse;

  /* PBPOLRAJAHIOPSPARSE models only support VARIABLE_WITHIN_BOUNDS
   * opflow->genbusvoltagetype
   */
  opflow->genbusvoltagetype = VARIABLE_WITHIN_BOUNDS;

  opflow->spdnordering = PETSC_FALSE;

  /* Inherit Ops */
  opflow->modelops.destroy = OPFLOWModelDestroy_PBPOLRAJAHIOPSPARSE;
  opflow->modelops.setnumvariables = OPFLOWModelSetNumVariables_PBPOL;
  opflow->modelops.setnumconstraints = OPFLOWModelSetNumConstraints_PBPOL;
  opflow->modelops.setvariablebounds =
      OPFLOWSetVariableBounds_PBPOLRAJAHIOPSPARSE;
  opflow->modelops.setvariableboundsarray =
      OPFLOWSetVariableBoundsArray_PBPOLRAJAHIOPSPARSE;
  opflow->modelops.setconstraintbounds =
      OPFLOWSetConstraintBounds_PBPOLRAJAHIOPSPARSE;
  opflow->modelops.setconstraintboundsarray =
      OPFLOWSetConstraintBoundsArray_PBPOLRAJAHIOPSPARSE;
  opflow->modelops.setinitialguess = OPFLOWSetInitialGuess_PBPOLRAJAHIOPSPARSE;
  opflow->modelops.setinitialguessarray =
      OPFLOWSetInitialGuessArray_PBPOLRAJAHIOPSPARSE;
  opflow->modelops.computeequalityconstraintsarray =
      OPFLOWComputeEqualityConstraintsArray_PBPOLRAJAHIOPSPARSE;
  opflow->modelops.computeinequalityconstraintsarray =
      OPFLOWComputeInequalityConstraintsArray_PBPOLRAJAHIOPSPARSE;
  opflow->modelops.computeobjectivearray =
      OPFLOWComputeObjectiveArray_PBPOLRAJAHIOPSPARSE;
  opflow->modelops.computegradientarray =
      OPFLOWComputeGradientArray_PBPOLRAJAHIOPSPARSE;
  opflow->modelops.solutiontops = OPFLOWSolutionToPS_PBPOLRAJAHIOPSPARSE;
  opflow->modelops.setup = OPFLOWModelSetUp_PBPOLRAJAHIOPSPARSE;
  opflow->modelops.computeequalityconstraintjacobian =
      OPFLOWComputeEqualityConstraintJacobian_PBPOL;
  opflow->modelops.computesparseequalityconstraintjacobianhiop =
      OPFLOWComputeSparseEqualityConstraintJacobian_PBPOLRAJAHIOPSPARSE;
  opflow->modelops.computeinequalityconstraintjacobian =
      OPFLOWComputeInequalityConstraintJacobian_PBPOL;
  opflow->modelops.computesparseinequalityconstraintjacobianhiop =
      OPFLOWComputeSparseInequalityConstraintJacobian_PBPOLRAJAHIOPSPARSE;
  opflow->modelops.computehessian = OPFLOWComputeHessian_PBPOL;
  opflow->modelops.computesparsehessianhiop =
      OPFLOWComputeSparseHessian_PBPOLRAJAHIOPSPARSE;
  opflow->modelops.solutioncallbackhiop =
      OPFLOWSolutionCallback_PBPOLRAJAHIOPSPARSE;

  PetscFunctionReturn(0);
}

#endif
#endif
