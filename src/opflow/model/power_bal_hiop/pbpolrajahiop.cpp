#include <exago_config.h>

#if defined(EXAGO_ENABLE_RAJA)

#include "pbpolrajahiopkernels.h"
#include <private/opflowimpl.h>

/* Initialization is done on the host through this function. Copying over values
 * to the device is done in OPFLOWSetInitialGuessArray_PBPOLRAJAHIOP
 */
extern PetscErrorCode OPFLOWSetInitialGuessArray_PBPOLHIOP(OPFLOW, double *);

PetscErrorCode OPFLOWSetInitialGuess_PBPOLRAJAHIOP(OPFLOW opflow, Vec X,
                                                   Vec Lambda) {
  PetscErrorCode ierr;
  double *x;

  PetscFunctionBegin;

  /* Get array pointers */
  ierr = VecGetArray(X, &x);
  CHKERRQ(ierr);

  ierr = OPFLOWSetInitialGuessArray_PBPOLHIOP(opflow, x);
  CHKERRQ(ierr);

  ierr = VecRestoreArray(X, &x);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

extern PetscErrorCode OPFLOWSetConstraintBounds_PBPOL(OPFLOW, Vec, Vec);

/* The constraint bounds are also calculated on the host.
 */
PetscErrorCode OPFLOWSetConstraintBounds_PBPOLRAJAHIOP(OPFLOW opflow, Vec Gl,
                                                       Vec Gu) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = OPFLOWSetConstraintBounds_PBPOL(opflow, Gl, Gu);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

extern PetscErrorCode OPFLOWSetVariableBounds_PBPOL(OPFLOW, Vec, Vec);

/* The variable bounds are also calculated on the host.
 */
PetscErrorCode OPFLOWSetVariableBounds_PBPOLRAJAHIOP(OPFLOW opflow, Vec Xl,
                                                     Vec Xu) {
  PetscErrorCode ierr;
  double *xl, *xu, *xlt, *xut;
  Vec Xlt, Xut;
  int i;

  PetscFunctionBegin;
  ierr = VecDuplicate(opflow->Xl, &Xlt);
  CHKERRQ(ierr);
  ierr = VecDuplicate(opflow->Xu, &Xut);
  CHKERRQ(ierr);
  ierr = OPFLOWSetVariableBounds_PBPOL(opflow, Xlt, Xut);
  CHKERRQ(ierr);

  ierr = VecGetArray(Xl, &xl);
  CHKERRQ(ierr);
  ierr = VecGetArray(Xu, &xu);
  CHKERRQ(ierr);
  ierr = VecGetArray(Xlt, &xlt);
  CHKERRQ(ierr);
  ierr = VecGetArray(Xut, &xut);
  CHKERRQ(ierr);

  /* Copy values from natural to sparse-dense format */
  for (i = 0; i < opflow->nx; i++) {
    xl[opflow->idxn2sd_map[i]] = xlt[i];
    xu[opflow->idxn2sd_map[i]] = xut[i];
  }

  ierr = VecRestoreArray(Xlt, &xlt);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(Xut, &xut);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(Xl, &xl);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(Xu, &xu);
  CHKERRQ(ierr);

  ierr = VecDestroy(&Xlt);
  CHKERRQ(ierr);
  ierr = VecDestroy(&Xut);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolutionToPS_PBPOLRAJAHIOP(OPFLOW opflow) {
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
      bus->pimb = x[loc] - x[loc + 1];
      bus->qimb = x[loc + 2] - x[loc + 3];
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
        load->pl_loss = x[loc];
        load->ql_loss = x[loc + 1];
      }
    }
  }

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

    Pf =
        Gff * Vmf * Vmf + Vmf * Vmt * (Gft * cos(thetaft) + Bft * sin(thetaft));
    Qf = -Bff * Vmf * Vmf +
         Vmf * Vmt * (-Bft * cos(thetaft) + Gft * sin(thetaft));

    Pt =
        Gtt * Vmt * Vmt + Vmt * Vmf * (Gtf * cos(thetatf) + Btf * sin(thetatf));
    Qt = -Btt * Vmt * Vmt +
         Vmt * Vmf * (-Btf * cos(thetatf) + Gtf * sin(thetatf));

    line->pf = Pf;
    line->qf = Qf;
    line->pt = Pt;
    line->qt = Qt;
    line->sf = PetscSqrtScalar(Pf * Pf + Qf * Qf);
    line->st = PetscSqrtScalar(Pt * Pt + Qt * Qt);

    if (opflow->ignore_lineflow_constraints || line->rateA > 1e5) {
      line->mult_sf = line->mult_st = 0.0;
    } else {
      gloc = line->startineqloc;
      line->mult_sf = lambdai[gloc];
      line->mult_st = lambdai[gloc + 1];
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

PetscErrorCode OPFLOWModelSetUp_PBPOLRAJAHIOP(OPFLOW opflow) {
  PetscErrorCode ierr;
  PbpolModelRajaHiop *pbpolrajahiop =
      reinterpret_cast<PbpolModelRajaHiop *>(opflow->model);
  PetscInt *idxn2sd_map;
  int ngen, nxsparse = 0, nxdense = 0;
  PS ps = (PS)opflow->ps;
  PSBUS bus;
  PSGEN gen;
  PSLOAD load;
  PSLINE line;
  PetscInt i, k;
  PetscInt loc, locglob;
  PetscInt eqloc = 0, ineqloc = 0;

  PetscFunctionBegin;

  ierr = OPFLOWModelSetUp_PBPOL(opflow);
  CHKERRQ(ierr);

  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    nxsparse += bus->nxshunt + bus->nxpimb;
    nxdense += bus->nxV;

    /* gen */
    for (k = 0; k < bus->ngen; k++) {
      ierr = PSBUSGetGen(bus, k, &gen);
      CHKERRQ(ierr);
      if (!gen->status)
        continue;
      nxsparse += gen->nx;
    }

    if (opflow->include_loadloss_variables) {
      for (k = 0; k < bus->nload; k++) {
        ierr = PSBUSGetLoad(bus, k, &load);
        CHKERRQ(ierr);
        if (!load->status)
          continue;
        nxsparse += load->nx;
      }
    }

    if (opflow->use_agc) {
      if (bus->ide == REF_BUS) {
        nxdense += ps->nx;
      }
    }
  }

  opflow->nxsparse = nxsparse;
  opflow->nxdense = nxdense;
  /* Create natural to sparse dense variable mapping */

  idxn2sd_map = opflow->idxn2sd_map;

  ps = opflow->ps;

  int spct = 0, dnct = 0;

  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    loc = bus->startxVloc;

    idxn2sd_map[loc] = nxsparse + dnct;
    idxn2sd_map[loc + 1] = nxsparse + dnct + 1;

    dnct += bus->nxV;

    if (opflow->include_powerimbalance_variables) { /* Power imbalance variables
                                                     */
      loc = bus->startxpimbloc;
      idxn2sd_map[loc] = spct;
      idxn2sd_map[loc + 1] = spct + 1;
      idxn2sd_map[loc + 2] = spct + 2;
      idxn2sd_map[loc + 3] = spct + 3;
      spct += bus->nxpimb;
    }

    PSBUSGetNGen(bus, &ngen);
    for (k = 0; k < ngen; k++) {
      PSBUSGetGen(bus, k, &gen);
      if (!gen->status)
        continue;

      loc = gen->startxpowloc;
      idxn2sd_map[loc] = spct;
      idxn2sd_map[loc + 1] = spct + 1;

      spct += gen->nxpow;

      if (opflow->has_gensetpoint) {
        loc = gen->startxpdevloc;
        idxn2sd_map[loc] = spct;

        spct += gen->nxpdev;

        loc = gen->startxpsetloc;

        idxn2sd_map[loc] = spct;

        spct += gen->nxpset;
      }
    }

    if (opflow->include_loadloss_variables) {
      for (k = 0; k < bus->nload; k++) {
        ierr = PSBUSGetLoad(bus, k, &load);
        CHKERRQ(ierr);
        if (!load->status)
          continue;

        loc = load->startxloadlossloc;
        idxn2sd_map[loc] = spct;
        idxn2sd_map[loc + 1] = spct + 1;

        spct += load->nxloadloss;
      }
    }

    if (opflow->use_agc) {
      if (bus->ide == REF_BUS) {
        idxn2sd_map[loc] = spct;

        spct += ps->nx;
        loc += ps->nx;
      }
    }
  }

  ierr = pbpolrajahiop->busparams.allocate(opflow);
  ierr = pbpolrajahiop->genparams.allocate(opflow);
  ierr = pbpolrajahiop->lineparams.allocate(opflow);
  ierr = pbpolrajahiop->loadparams.allocate(opflow);

  BUSParamsRajaHiop *busparams = &pbpolrajahiop->busparams;
  GENParamsRajaHiop *genparams = &pbpolrajahiop->genparams;
  LOADParamsRajaHiop *loadparams = &pbpolrajahiop->loadparams;
  LINEParamsRajaHiop *lineparams = &pbpolrajahiop->lineparams;
  int geni = 0, loadi = 0, gi;
  int nnz_eqjacsp = 0;
  int nnz_ineqjacsp = 0;
  int nnz_hesssp = 0;
  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    // We need to go over each row since HiOp sparse Jacobian needs the values
    // to be inserted row-wise (!)
    // So we first go over the partial derivatives w.r.t. real power balance
    // equations and then over the second one for reactive power balance
    // equations.

    if (opflow->include_powerimbalance_variables) {
      busparams->jacsp_idx[i] = nnz_eqjacsp;
      nnz_eqjacsp += 2;
    }

    gi = 0;
    for (k = 0; k < bus->ngen; k++) {
      ierr = PSBUSGetGen(bus, k, &gen);
      CHKERRQ(ierr);
      if (!gen->status)
        continue;
      genparams->eqjacspbus_idx[geni + gi] = nnz_eqjacsp;
      genparams->hesssp_idx[geni + gi] = nnz_hesssp;
      nnz_eqjacsp += 1;
      nnz_hesssp += 1;
      gi++;
    }

    if (opflow->include_loadloss_variables) {
      for (k = 0; k < bus->nload; k++) {
        ierr = PSBUSGetLoad(bus, k, &load);
        CHKERRQ(ierr);
        loadparams->jacsp_idx[loadi + k] = nnz_eqjacsp;
        loadparams->hesssp_idx[loadi + k] = nnz_hesssp;
        nnz_eqjacsp += 1;
        nnz_hesssp += 2;
      }
    }

    /* Jacobian elements for reactive power contributions */

    if (opflow->include_powerimbalance_variables) {
      busparams->jacsq_idx[i] = nnz_eqjacsp;
      nnz_eqjacsp += 2;
    }

    gi = 0;
    for (k = 0; k < bus->ngen; k++) {
      ierr = PSBUSGetGen(bus, k, &gen);
      CHKERRQ(ierr);
      if (!gen->status)
        continue;
      genparams->eqjacsqbus_idx[geni + gi] = nnz_eqjacsp;
      nnz_eqjacsp += 1;
      gi++;
    }

    if (opflow->include_loadloss_variables) {
      for (k = 0; k < bus->nload; k++) {
        ierr = PSBUSGetLoad(bus, k, &load);
        CHKERRQ(ierr);
        loadparams->jacsq_idx[loadi + k] = nnz_eqjacsp;
        nnz_eqjacsp += 1;
      }
    }

    gi = 0;
    if (opflow->has_gensetpoint) {
      /* Jacobian contributions for generator set-point formulation */
      for (k = 0; k < bus->ngen; k++) {
        ierr = PSBUSGetGen(bus, k, &gen);
        CHKERRQ(ierr);
        if (!gen->status)
          continue;
        genparams->eqjacspgen_idx[geni + gi] = nnz_eqjacsp;
        nnz_eqjacsp += 4; /* 4 Jacobian elements contributed to equality
                             constrained Jacobian */
        genparams->ineqjacspgen_idx[geni + gi] = nnz_ineqjacsp;
        nnz_ineqjacsp += 0; /* 0 Jacobian elements contributed to inequality
                               constrained Jacobian */
        gi++;
        nnz_hesssp += 0; /* 3 Hessian elements from inequality constraints */
      }
    }

    geni += bus->ngenON;
    loadi += bus->nload;
  }

  /* Save the number of nonzeros so we can use it for HIOP blocks info */
  opflow->nnz_eqjacsp = nnz_eqjacsp;
  opflow->nnz_ineqjacsp = nnz_ineqjacsp;
  opflow->nnz_hesssp = nnz_hesssp;

  ierr = pbpolrajahiop->busparams.copy(opflow);
  ierr = pbpolrajahiop->genparams.copy(opflow);
  ierr = pbpolrajahiop->lineparams.copy(opflow);
  ierr = pbpolrajahiop->loadparams.copy(opflow);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWModelDestroy_PBPOLRAJAHIOP(OPFLOW opflow) {
  PetscErrorCode ierr;
  PbpolModelRajaHiop *pbpolrajahiop =
      reinterpret_cast<PbpolModelRajaHiop *>(opflow->model);

  PetscFunctionBegin;
  pbpolrajahiop->destroy(opflow);
  delete pbpolrajahiop;
  pbpolrajahiop = nullptr;

  PetscFunctionReturn(0);
}

/* reuse numvariables and numconstraints functions from PBPOL model */
extern PetscErrorCode OPFLOWModelSetNumVariables_PBPOL(OPFLOW, PetscInt *,
                                                       PetscInt *, PetscInt *);
extern PetscErrorCode OPFLOWModelSetNumConstraints_PBPOL(OPFLOW, PetscInt *,
                                                         PetscInt *, PetscInt *,
                                                         PetscInt *);

PetscErrorCode OPFLOWModelCreate_PBPOLRAJAHIOP(OPFLOW opflow) {
  PetscErrorCode ierr;

  PetscFunctionBegin;

  PbpolModelRajaHiop *pbpol = new PbpolModelRajaHiop(opflow);

  opflow->model = pbpol;

  /* HIOP models only support VARIABLE_WITHIN_BOUNDS opflow->genbusvoltagetype
   */
  opflow->genbusvoltagetype = VARIABLE_WITHIN_BOUNDS;

  opflow->spdnordering = PETSC_TRUE;

  /* Inherit Ops */
  opflow->modelops.destroy = OPFLOWModelDestroy_PBPOLRAJAHIOP;
  opflow->modelops.setnumvariables = OPFLOWModelSetNumVariables_PBPOL;
  opflow->modelops.setnumconstraints = OPFLOWModelSetNumConstraints_PBPOL;
  opflow->modelops.setvariablebounds = OPFLOWSetVariableBounds_PBPOLRAJAHIOP;
  opflow->modelops.setvariableboundsarray =
      OPFLOWSetVariableBoundsArray_PBPOLRAJAHIOP;
  opflow->modelops.setconstraintbounds =
      OPFLOWSetConstraintBounds_PBPOLRAJAHIOP;
  opflow->modelops.setconstraintboundsarray =
      OPFLOWSetConstraintBoundsArray_PBPOLRAJAHIOP;
  opflow->modelops.setinitialguess = OPFLOWSetInitialGuess_PBPOLRAJAHIOP;
  opflow->modelops.setinitialguessarray =
      OPFLOWSetInitialGuessArray_PBPOLRAJAHIOP;
  opflow->modelops.computeequalityconstraintsarray =
      OPFLOWComputeEqualityConstraintsArray_PBPOLRAJAHIOP;
  opflow->modelops.computeinequalityconstraintsarray =
      OPFLOWComputeInequalityConstraintsArray_PBPOLRAJAHIOP;
  opflow->modelops.computeobjectivearray =
      OPFLOWComputeObjectiveArray_PBPOLRAJAHIOP;
  opflow->modelops.computegradientarray =
      OPFLOWComputeGradientArray_PBPOLRAJAHIOP;
  opflow->modelops.solutiontops = OPFLOWSolutionToPS_PBPOLRAJAHIOP;
  opflow->modelops.setup = OPFLOWModelSetUp_PBPOLRAJAHIOP;
  opflow->modelops.computesparseequalityconstraintjacobianhiop =
      OPFLOWComputeSparseEqualityConstraintJacobian_PBPOLRAJAHIOP;
  opflow->modelops.computesparseinequalityconstraintjacobianhiop =
      OPFLOWComputeSparseInequalityConstraintJacobian_PBPOLRAJAHIOP;
  opflow->modelops.computesparsehessianhiop =
      OPFLOWComputeSparseHessian_PBPOLRAJAHIOP;
  opflow->modelops.computedenseequalityconstraintjacobianhiop =
      OPFLOWComputeDenseEqualityConstraintJacobian_PBPOLRAJAHIOP;
  opflow->modelops.computedenseinequalityconstraintjacobianhiop =
      OPFLOWComputeDenseInequalityConstraintJacobian_PBPOLRAJAHIOP;
  opflow->modelops.computedensehessianhiop =
      OPFLOWComputeDenseHessian_PBPOLRAJAHIOP;

  PetscFunctionReturn(0);
}

#endif
