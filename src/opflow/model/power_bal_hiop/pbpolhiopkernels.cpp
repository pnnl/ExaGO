#include <exago_config.h>

#if defined(EXAGO_ENABLE_HIOP)

#include "pbpolhiop.h"
#include <private/opflowimpl.h>
#include <private/psimpl.h>

/** INITIAL GUESS **/
PetscErrorCode OPFLOWSetInitialGuessArray_PBPOLHIOP(OPFLOW opflow, double *x) {
  PetscErrorCode ierr;
  PS ps = opflow->ps;
  const PetscScalar *xl, *xu;
  PetscInt i;
  PSBUS bus;
  PetscInt loc, loc_nat;
  int *idxn2sd_map = opflow->idxn2sd_map;

  PetscFunctionBegin;

  /* Get array pointers */
  ierr = VecGetArrayRead(opflow->Xl, &xl);
  CHKERRQ(ierr);
  ierr = VecGetArrayRead(opflow->Xu, &xu);
  CHKERRQ(ierr);

  for (i = 0; i < ps->nbus; i++) {
    PetscInt k;

    bus = &ps->bus[i];

    loc_nat = bus->startxVloc;
    loc = idxn2sd_map[loc_nat];

    if (bus->ide == ISOLATED_BUS) {
      x[loc] = bus->va * PETSC_PI / 180.0;
      x[loc + 1] = bus->vm;
    } else {
      if (opflow->initializationtype == OPFLOWINIT_MIDPOINT) {
        /* Initial guess for voltage angles and bounds on voltage magnitudes */
        x[loc] = (xl[loc] + xu[loc]) / 2.0;
        x[loc + 1] = (xl[loc + 1] + xu[loc + 1]) / 2.0;
      } else if (opflow->initializationtype == OPFLOWINIT_FROMFILE ||
                 opflow->initializationtype == OPFLOWINIT_ACPF) {
        x[loc] = bus->va * PETSC_PI / 180.0;
        x[loc + 1] = PetscMax(bus->Vmin, PetscMin(bus->vm, bus->Vmax));
      } else if (opflow->initializationtype == OPFLOWINIT_FLATSTART) {
        x[loc] = 0.0;
        x[loc + 1] = 1.0;
      }
    }

    if (opflow->include_powerimbalance_variables) {
      loc_nat = bus->startxpimbloc;
      loc = idxn2sd_map[loc_nat];
      x[loc] = x[loc + 1] = x[loc + 2] = x[loc + 3] = 0.0;
    }

    for (k = 0; k < bus->ngen; k++) {
      PSGEN gen;

      ierr = PSBUSGetGen(bus, k, &gen);
      CHKERRQ(ierr);
      if (!gen->status)
        continue;

      loc_nat = gen->startxpowloc;
      loc = idxn2sd_map[loc_nat];

      if (opflow->initializationtype == OPFLOWINIT_MIDPOINT ||
          opflow->initializationtype == OPFLOWINIT_FLATSTART) {
        x[loc] = 0.5 * (xl[loc] + xu[loc]);
        x[loc + 1] = 0.5 * (xl[loc + 1] + xu[loc + 1]);
      } else if (opflow->initializationtype == OPFLOWINIT_FROMFILE ||
                 opflow->initializationtype == OPFLOWINIT_ACPF) {
        x[loc] = PetscMax(gen->pb, PetscMin(gen->pg, gen->pt));
        x[loc + 1] = PetscMax(gen->qb, PetscMin(gen->qg, gen->qt));
      }

      if (opflow->has_gensetpoint) {
        loc_nat = gen->startxpdevloc;
        loc = idxn2sd_map[loc_nat];
        x[loc] = 0.0;
        loc_nat = gen->startxpsetloc;
        loc = idxn2sd_map[loc_nat];
        x[loc] = gen->pgs;
      }
    }

    if (opflow->include_loadloss_variables) {
      for (k = 0; k < bus->nload; k++) {
        PSLOAD load;
        ierr = PSBUSGetLoad(bus, k, &load);
        CHKERRQ(ierr);
        if (!load->status)
          continue;

        loc_nat = load->startxloadlossloc;
        loc = idxn2sd_map[loc_nat];

        /* Initial value for real and reactive power load loss */
        x[loc] = 0.0;
        x[loc + 1] = 0.0;
      }
    }
    if (opflow->use_agc) {
      if (bus->ide == REF_BUS) {
        loc_nat = ps->startxloc;
        loc = idxn2sd_map[loc_nat];
        x[loc] = 0.0;
      }
    }
  }

  ierr = VecRestoreArrayRead(opflow->Xl, &xl);
  CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(opflow->Xu, &xu);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/** CONSTRAINT BOUNDS  **/
PetscErrorCode OPFLOWSetConstraintBoundsArray_PBPOLHIOP(OPFLOW opflow,
                                                        double *gl,
                                                        double *gu) {
  PBPOLHIOP pbpolhiop = (PBPOLHIOP)opflow->model;
  BUSParams *busparams = &pbpolhiop->busparams;
  LINEParams *lineparams = &pbpolhiop->lineparams;
  GENParams *genparams = &pbpolhiop->genparams;
  PetscInt i;
  PS ps = opflow->ps;

  PetscFunctionBegin;

  /* Equality constraints (all zeros) */
  for (i = 0; i < busparams->nbus; i++) {
    gl[busparams->gidx[i]] = 0.0;
    gu[busparams->gidx[i]] = 0.0;

    gl[busparams->gidx[i] + 1] = 0.0;
    gu[busparams->gidx[i] + 1] = 0.0;
  }

  if (opflow->has_gensetpoint) {
    for (i = 0; i < genparams->ngenON; i++) {
      gl[genparams->geqidxgen[i]] = gu[genparams->geqidxgen[i]] = 0.0;
      gl[genparams->geqidxgen[i] + 1] = gu[genparams->geqidxgen[i] + 1] = 0.0;
    }
  }

  /* Inequality constraints */
  for (i = 0; i < lineparams->nlinelim; i++) {
    int j = lineparams->linelimidx[i];
    gl[lineparams->gbineqidx[i]] = 0.0;
    gu[lineparams->gbineqidx[i]] = (lineparams->rateA[j] / ps->MVAbase) *
                                   (lineparams->rateA[j] / ps->MVAbase);
    gl[lineparams->gbineqidx[i] + 1] = 0.0;
    gu[lineparams->gbineqidx[i] + 1] = (lineparams->rateA[j] / ps->MVAbase) *
                                       (lineparams->rateA[j] / ps->MVAbase);
  }

  PetscFunctionReturn(0);
}

/** EQUALITY CONSTRAINTS */
PetscErrorCode OPFLOWComputeEqualityConstraintsArray_PBPOLHIOP(OPFLOW opflow,
                                                               const double *x,
                                                               double *ge) {
  PBPOLHIOP pbpolhiop = (PBPOLHIOP)opflow->model;
  BUSParams *busparams = &pbpolhiop->busparams;
  GENParams *genparams = &pbpolhiop->genparams;
  LOADParams *loadparams = &pbpolhiop->loadparams;
  LINEParams *lineparams = &pbpolhiop->lineparams;
  PetscInt i;
  PetscInt flps = 0;
  PetscErrorCode ierr;
  double Pg, delPg, Pgset;

  PetscFunctionBegin;

  for (i = 0; i < opflow->nconeq; i++)
    ge[i] = 0.0;

  /* Generator contributions */
  for (i = 0; i < genparams->ngenON; i++) {
    ge[genparams->gidxbus[i]] -= x[genparams->xidx[i]];
    ge[genparams->gidxbus[i] + 1] -= x[genparams->xidx[i] + 1];

    flps += 2;
    if (opflow->has_gensetpoint) {
      Pg = x[genparams->xidx[i]];
      delPg = x[genparams->xidx[i] + 2];
      Pgset = x[genparams->xidx[i] + 3];

      ge[genparams->geqidxgen[i]] = Pgset + delPg - Pg;
      ge[genparams->geqidxgen[i] + 1] = Pgset - genparams->pgs[i];

      flps += 3;
    }
  }

  /* Load contributions */
  for (i = 0; i < loadparams->nload; i++) {
    if (opflow->include_loadloss_variables) {
      ge[loadparams->gidx[i]] += loadparams->pl[i] - x[loadparams->xidx[i]];
      ge[loadparams->gidx[i] + 1] +=
          loadparams->ql[i] - x[loadparams->xidx[i] + 1];
    } else {
      ge[loadparams->gidx[i]] += loadparams->pl[i];
      ge[loadparams->gidx[i] + 1] += loadparams->ql[i];
    }
  }
  flps += loadparams->nload * 2;

  /* Bus contributions */
  for (i = 0; i < busparams->nbus; i++) {
    double theta = x[busparams->xidx[i]];
    double Vm = x[busparams->xidx[i] + 1];
    ge[busparams->gidx[i]] +=
        busparams->isisolated[i] *
            (theta - busparams->va[i] * PETSC_PI / 180.0) +
        busparams->ispvpq[i] * Vm * Vm * busparams->gl[i];
    ge[busparams->gidx[i] + 1] +=
        busparams->isisolated[i] * (Vm - busparams->vm[i]) -
        busparams->ispvpq[i] * Vm * Vm * busparams->bl[i];

    /* Power imbalance addition (Second Bus Variable)*/
    if (opflow->include_powerimbalance_variables) {
      double Pimbplus, Pimbminus, Qimbplus, Qimbminus, Pimb, Qimb;
      Pimbplus = x[busparams->xidxpimb[i]];
      Pimbminus = x[busparams->xidxpimb[i] + 1];
      Qimbplus = x[busparams->xidxpimb[i] + 2];
      Qimbminus = x[busparams->xidxpimb[i] + 3];

      Pimb = Pimbplus - Pimbminus;
      Qimb = Qimbplus - Qimbminus;
      ge[busparams->gidx[i]] += Pimb;
      ge[busparams->gidx[i] + 1] += Qimb;
    }
  }
  flps += busparams->nbus * 14;

  /* Line contributions */
  for (i = 0; i < lineparams->nlineON; i++) {
    double Pf, Qf, Pt, Qt;
    double thetaf = x[lineparams->xidxf[i]], Vmf = x[lineparams->xidxf[i] + 1];
    double thetat = x[lineparams->xidxt[i]], Vmt = x[lineparams->xidxt[i] + 1];
    double thetaft = thetaf - thetat;
    double thetatf = thetat - thetaf;

    Pf = lineparams->Gff[i] * Vmf * Vmf +
         Vmf * Vmt *
             (lineparams->Gft[i] * cos(thetaft) +
              lineparams->Bft[i] * sin(thetaft));
    Qf = -lineparams->Bff[i] * Vmf * Vmf +
         Vmf * Vmt *
             (-lineparams->Bft[i] * cos(thetaft) +
              lineparams->Gft[i] * sin(thetaft));
    Pt = lineparams->Gtt[i] * Vmt * Vmt +
         Vmt * Vmf *
             (lineparams->Gtf[i] * cos(thetatf) +
              lineparams->Btf[i] * sin(thetatf));
    Qt = -lineparams->Btt[i] * Vmt * Vmt +
         Vmt * Vmf *
             (-lineparams->Btf[i] * cos(thetatf) +
              lineparams->Gtf[i] * sin(thetatf));

    /* Atomic operation */
    ge[lineparams->geqidxf[i]] += Pf;
    ge[lineparams->geqidxf[i] + 1] += Qf;
    ge[lineparams->geqidxt[i]] += Pt;
    ge[lineparams->geqidxt[i] + 1] += Qt;
  }
  flps += lineparams->nlineON *
          (46 + (4 * EXAGO_FLOPS_COSOP) + (4 * EXAGO_FLOPS_SINOP));

  ierr = PetscLogFlops(flps);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/** INEQUALITY CONSTRAINTS **/
PetscErrorCode
OPFLOWComputeInequalityConstraintsArray_PBPOLHIOP(OPFLOW opflow,
                                                  const double *x, double *gi) {
  PBPOLHIOP pbpolhiop = (PBPOLHIOP)opflow->model;
  GENParams *genparams = &pbpolhiop->genparams;
  LINEParams *lineparams = &pbpolhiop->lineparams;
  PetscInt i;
  PetscInt flps = 0;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  for (i = 0; i < opflow->nconineq; i++)
    gi[i] = 0.0;

  /* Line contributions */
  for (i = 0; i < lineparams->nlinelim; i++) {
    int j = lineparams->linelimidx[i];
    double Pf, Qf, Pt, Qt, Sf2, St2;
    double thetaf = x[lineparams->xidxf[j]], Vmf = x[lineparams->xidxf[j] + 1];
    double thetat = x[lineparams->xidxt[j]], Vmt = x[lineparams->xidxt[j] + 1];
    double thetaft = thetaf - thetat;
    double thetatf = thetat - thetaf;

    Pf = lineparams->Gff[j] * Vmf * Vmf +
         Vmf * Vmt *
             (lineparams->Gft[j] * cos(thetaft) +
              lineparams->Bft[j] * sin(thetaft));
    Qf = -lineparams->Bff[j] * Vmf * Vmf +
         Vmf * Vmt *
             (-lineparams->Bft[j] * cos(thetaft) +
              lineparams->Gft[j] * sin(thetaft));
    Pt = lineparams->Gtt[j] * Vmt * Vmt +
         Vmt * Vmf *
             (lineparams->Gtf[j] * cos(thetatf) +
              lineparams->Btf[j] * sin(thetatf));
    Qt = -lineparams->Btt[j] * Vmt * Vmt +
         Vmt * Vmf *
             (-lineparams->Btf[j] * cos(thetatf) +
              lineparams->Gtf[j] * sin(thetatf));

    Sf2 = Pf * Pf + Qf * Qf;
    St2 = Pt * Pt + Qt * Qt;

    gi[lineparams->gineqidx[i]] = Sf2;
    gi[lineparams->gineqidx[i] + 1] = St2;
  }
  flps += lineparams->nlinelim *
          (72 + (4 * EXAGO_FLOPS_COSOP) + (4 * EXAGO_FLOPS_SINOP));

  ierr = PetscLogFlops(flps);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/** OBJECTIVE FUNCTION **/
PetscErrorCode OPFLOWComputeObjectiveArray_PBPOLHIOP(OPFLOW opflow,
                                                     const double *x,
                                                     double *obj) {
  PBPOLHIOP pbpolhiop = (PBPOLHIOP)opflow->model;
  GENParams *genparams = &pbpolhiop->genparams;
  LOADParams *loadparams = &pbpolhiop->loadparams;
  BUSParams *busparams = &pbpolhiop->busparams;
  PetscInt i;
  PetscErrorCode ierr;
  PS ps = opflow->ps;
  double obj_val = 0.0;
  int isobj_gencost = opflow->obj_gencost;
  double MVAbase = ps->MVAbase;
  double powerimbal_penalty = opflow->powerimbalance_penalty;

  PetscFunctionBegin;

  /* Power Imbalance Contributions (BUS) */
  if (opflow->include_powerimbalance_variables) {
    for (int i = 0; i < busparams->nbus; i++) {
      double Pimbplus, Pimbminus, Qimbplus, Qimbminus;
      Pimbplus = x[busparams->xidxpimb[i]];
      Pimbminus = x[busparams->xidxpimb[i] + 1];
      Qimbplus = x[busparams->xidxpimb[i] + 2];
      Qimbminus = x[busparams->xidxpimb[i] + 3];

      obj_val += busparams->powerimbalance_penalty[i] * ps->MVAbase *
                 (Pimbplus + Pimbminus + Qimbplus + Qimbminus);
    }
  }

  if (opflow->objectivetype == MIN_GEN_COST) {
    /* Generator objective function contributions */
    for (i = 0; i < genparams->ngenON; i++) {
      double Pg = x[genparams->xidx[i]] * MVAbase;
      obj_val += isobj_gencost *
                 (genparams->cost_alpha[i] * Pg * Pg +
                  genparams->cost_beta[i] * Pg + genparams->cost_gamma[i]);
    }
  }

  /* LoadLoss objective function contribution if present */
  if (opflow->include_loadloss_variables) {
    PetscScalar Pdloss, Qdloss;
    for (i = 0; i < loadparams->nload; i++) {
      Pdloss = x[loadparams->xidx[i]];
      Qdloss = x[loadparams->xidx[i] + 1];
      obj_val +=
          loadparams->loadloss_penalty[i] * ps->MVAbase * (Pdloss + Qdloss);
    }
  }
  *obj = obj_val;
  ierr = PetscLogFlops(genparams->ngenON * 8.0);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/** GRADIENT **/
PetscErrorCode OPFLOWComputeGradientArray_PBPOLHIOP(OPFLOW opflow,
                                                    const double *x,
                                                    double *grad) {
  PBPOLHIOP pbpolhiop = (PBPOLHIOP)opflow->model;
  GENParams *genparams = &pbpolhiop->genparams;
  LOADParams *loadparams = &pbpolhiop->loadparams;
  BUSParams *busparams = &pbpolhiop->busparams;
  PetscInt i;
  PetscErrorCode ierr;
  PS ps = opflow->ps;
  int isobj_gencost = opflow->obj_gencost;
  double MVAbase = ps->MVAbase;
  double powerimbal_penalty = opflow->powerimbalance_penalty;

  PetscFunctionBegin;

  for (i = 0; i < opflow->nx; i++)
    grad[i] = 0.0;

  /* Power Imbalance Contributions (Second BUS Variable) */
  if (opflow->include_powerimbalance_variables) {
    for (int i = 0; i < busparams->nbus; i++) {
      grad[busparams->xidxpimb[i]] =
          busparams->powerimbalance_penalty[i] * ps->MVAbase;
      grad[busparams->xidxpimb[i] + 1] =
          busparams->powerimbalance_penalty[i] * ps->MVAbase * 1.0;
      grad[busparams->xidxpimb[i] + 2] =
          busparams->powerimbalance_penalty[i] * ps->MVAbase;
      grad[busparams->xidxpimb[i] + 3] =
          busparams->powerimbalance_penalty[i] * ps->MVAbase;
    }
  }

  if (opflow->objectivetype == MIN_GEN_COST) {
    /* Generator gradient contributions */
    for (i = 0; i < genparams->ngenON; i++) {
      double Pg = x[genparams->xidx[i]] * MVAbase;
      grad[genparams->xidx[i]] =
          isobj_gencost * MVAbase *
          (2.0 * genparams->cost_alpha[i] * Pg + genparams->cost_beta[i]);
    }
  }

  /* Loadloss gradient contributions */
  if (opflow->include_loadloss_variables) {
    PetscScalar Pdloss, Qdloss;
    for (i = 0; i < loadparams->nload; i++) {
      Pdloss = x[loadparams->xidx[i]];
      Qdloss = x[loadparams->xidx[i] + 1];
      grad[loadparams->xidx[i]] = loadparams->loadloss_penalty[i] * ps->MVAbase;
      grad[loadparams->xidx[i] + 1] =
          loadparams->loadloss_penalty[i] * ps->MVAbase;
    }
  }

  ierr = PetscLogFlops(genparams->ngenON * 6);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/** VARIABLE BOUNDS **/
PetscErrorCode OPFLOWSetVariableBoundsArray_PBPOLHIOP(OPFLOW opflow, double *xl,
                                                      double *xu) {
  PBPOLHIOP pbpolhiop = (PBPOLHIOP)opflow->model;
  BUSParams *busparams = &pbpolhiop->busparams;
  GENParams *genparams = &pbpolhiop->genparams;
  LOADParams *loadparams = &pbpolhiop->loadparams;
  PetscInt i;

  PetscFunctionBegin;

  /* Bounds for bus voltages */
  for (i = 0; i < busparams->nbus; i++) {
    xl[busparams->xidx[i]] =
        busparams->ispvpq[i] * PETSC_NINFINITY +
        busparams->isisolated[i] * busparams->va[i] +
        busparams->isref[i] * busparams->va[i] * PETSC_PI / 180.0;
    xu[busparams->xidx[i]] =
        busparams->ispvpq[i] * PETSC_INFINITY +
        busparams->isisolated[i] * busparams->va[i] +
        busparams->isref[i] * busparams->va[i] * PETSC_PI / 180.0;

    xl[busparams->xidx[i] + 1] = busparams->isref[i] * busparams->vmin[i] +
                                 busparams->ispvpq[i] * busparams->vmin[i] +
                                 busparams->isisolated[i] * busparams->vm[i];
    xu[busparams->xidx[i] + 1] = busparams->isref[i] * busparams->vmax[i] +
                                 busparams->ispvpq[i] * busparams->vmax[i] +
                                 busparams->isisolated[i] * busparams->vm[i];

    /* Bounds for Power Imbalance Variables (second bus variables) */
    if (opflow->include_powerimbalance_variables) {
      xl[busparams->xidxpimb[i]] = xl[busparams->xidxpimb[i] + 1] =
          xl[busparams->xidxpimb[i] + 2] = xl[busparams->xidxpimb[i] + 3] = 0.0;
      xu[busparams->xidxpimb[i]] = xu[busparams->xidxpimb[i] + 1] =
          xu[busparams->xidxpimb[i] + 2] = xu[busparams->xidxpimb[i] + 3] =
              PETSC_INFINITY;
    }
  }

  /* Generator lower and upper bounds on variables */
  for (i = 0; i < genparams->ngenON; i++) {
    xl[genparams->xidx[i]] = genparams->pb[i];
    xu[genparams->xidx[i]] = genparams->pt[i];
    xl[genparams->xidx[i] + 1] = genparams->qb[i];
    xu[genparams->xidx[i] + 1] = genparams->qt[i];

    if (opflow->has_gensetpoint) {
      xl[genparams->xidx[i] + 2] = genparams->pb[i] - genparams->pt[i];
      xu[genparams->xidx[i] + 2] = genparams->pt[i] - genparams->pb[i];
      xl[genparams->xidx[i] + 3] = genparams->pb[i];
      xu[genparams->xidx[i] + 3] =
          10000.0; // genparams->pt[i]; This is a temporary fix for now. The
                   // proper fix is to check the generator fuel type and set it
                   // o gen->pt if it is non-renewable and 10000.0 otherwise.
    }
  }

  /* Load loss lower and upper bounds */
  if (opflow->include_loadloss_variables) {
    for (i = 0; i < loadparams->nload; i++) {
      xl[loadparams->xidx[i]] = PetscMin(0.0, loadparams->pl[i]);
      xu[loadparams->xidx[i]] = PetscMax(0.0, loadparams->pl[i]);
      xl[loadparams->xidx[i] + 1] = PetscMin(0.0, loadparams->ql[i]);
      xu[loadparams->xidx[i] + 1] = PetscMax(0.0, loadparams->ql[i]);
    }
  }

  PetscFunctionReturn(0);
}

/** Custom routines that work with HIOP interface only */
PetscErrorCode OPFLOWComputeSparseEqualityConstraintJacobian_PBPOLHIOP(
    OPFLOW opflow, const double *x, int *iJacS, int *jJacS, double *MJacS) {
  PBPOLHIOP pbpolhiop = (PBPOLHIOP)opflow->model;
  GENParams *genparams = &pbpolhiop->genparams;
  LOADParams *loadparams = &pbpolhiop->loadparams;
  BUSParams *busparams = &pbpolhiop->busparams;

  int i;

  if (iJacS != NULL && jJacS != NULL) {

    /* Power Imbalance Contributions (BUS Second Variables?) */
    if (opflow->include_powerimbalance_variables) {
      for (int i = 0; i < busparams->nbus; i++) {
        iJacS[busparams->jacsp_idx[i]] = busparams->gidx[i];
        jJacS[busparams->jacsp_idx[i]] = busparams->xidxpimb[i];
        iJacS[busparams->jacsp_idx[i] + 1] = busparams->gidx[i];
        jJacS[busparams->jacsp_idx[i] + 1] = busparams->xidxpimb[i] + 1;

        iJacS[busparams->jacsq_idx[i]] = busparams->gidx[i] + 1;
        jJacS[busparams->jacsq_idx[i]] = busparams->xidxpimb[i] + 2;
        iJacS[busparams->jacsq_idx[i] + 1] = busparams->gidx[i] + 1;
        jJacS[busparams->jacsq_idx[i] + 1] = busparams->xidxpimb[i] + 3;
      }
    }

    /* Generator contributions */
    for (i = 0; i < genparams->ngenON; i++) {
      iJacS[genparams->eqjacspbus_idx[i]] = genparams->gidxbus[i];
      jJacS[genparams->eqjacspbus_idx[i]] = genparams->xidx[i];

      iJacS[genparams->eqjacsqbus_idx[i]] = genparams->gidxbus[i] + 1;
      jJacS[genparams->eqjacsqbus_idx[i]] = genparams->xidx[i] + 1;

      if (opflow->has_gensetpoint) {
        iJacS[genparams->eqjacspgen_idx[i]] = genparams->geqidxgen[i];
        jJacS[genparams->eqjacspgen_idx[i]] = genparams->xidx[i]; // Pg

        iJacS[genparams->eqjacspgen_idx[i] + 1] = genparams->geqidxgen[i];
        jJacS[genparams->eqjacspgen_idx[i] + 1] =
            genparams->xidx[i] + 2; // delPg

        iJacS[genparams->eqjacspgen_idx[i] + 2] = genparams->geqidxgen[i];
        jJacS[genparams->eqjacspgen_idx[i] + 2] =
            genparams->xidx[i] + 3; // Pgset

        iJacS[genparams->eqjacspgen_idx[i] + 3] = genparams->geqidxgen[i] + 1;
        jJacS[genparams->eqjacspgen_idx[i] + 3] =
            genparams->xidx[i] + 3; // Pgset
      }
    }

    /* Loadloss contributions */
    if (opflow->include_loadloss_variables) {
      for (i = 0; i < loadparams->nload; i++) {
        iJacS[loadparams->jacsp_idx[i]] = loadparams->gidx[i];
        jJacS[loadparams->jacsp_idx[i]] = loadparams->xidx[i];
        iJacS[loadparams->jacsq_idx[i]] = loadparams->gidx[i] + 1;
        jJacS[loadparams->jacsq_idx[i]] = loadparams->xidx[i] + 1;
      }
    }
  }

  if (MJacS != NULL) {
    /* Power Imbalance Contributions (BUS) */
    if (opflow->include_powerimbalance_variables) {
      for (int i = 0; i < busparams->nbus; i++) {
        MJacS[busparams->jacsp_idx[i]] = 1.0;
        MJacS[busparams->jacsp_idx[i] + 1] = -1.0;
        MJacS[busparams->jacsq_idx[i]] = 1.0;
        MJacS[busparams->jacsq_idx[i] + 1] = -1.0;
      }
    }

    /* Generator contributions */
    for (i = 0; i < genparams->ngenON; i++) {
      MJacS[genparams->eqjacspbus_idx[i]] = -1.0;
      MJacS[genparams->eqjacsqbus_idx[i]] = -1.0;

      if (opflow->has_gensetpoint) {
        MJacS[genparams->eqjacspgen_idx[i]] = -1.0; // Pg

        MJacS[genparams->eqjacspgen_idx[i] + 1] = 1.0; // delPg

        MJacS[genparams->eqjacspgen_idx[i] + 2] = 1.0; // Pgset

        MJacS[genparams->eqjacspgen_idx[i] + 3] = 1.0; // Pgset
      }
    }

    /* Jacobian from loadloss contribution */
    if (opflow->include_loadloss_variables) {
      for (i = 0; i < loadparams->nload; i++) {
        MJacS[loadparams->jacsp_idx[i]] = -1;
        MJacS[loadparams->jacsq_idx[i]] = -1;
      }
    }
  }

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeSparseInequalityConstraintJacobian_PBPOLHIOP(
    OPFLOW opflow, const double *x, int *iJacS, int *jJacS, double *MJacS) {
  PBPOLHIOP pbpolhiop = (PBPOLHIOP)opflow->model;
  GENParams *genparams = &pbpolhiop->genparams;
  LOADParams *loadparams = &pbpolhiop->loadparams;

  int i;

  if (iJacS != NULL && jJacS != NULL) {
  }

  if (MJacS != NULL) {
  }

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeSparseHessian_PBPOLHIOP(OPFLOW opflow,
                                                    const double *x,
                                                    const double *lambda,
                                                    int *iHSS, int *jHSS,
                                                    double *MHSS) {
  PBPOLHIOP pbpolhiop = (PBPOLHIOP)opflow->model;
  PetscErrorCode ierr;
  GENParams *genparams = &pbpolhiop->genparams;
  LOADParams *loadparams = &pbpolhiop->loadparams;
  BUSParams *busparams = &pbpolhiop->busparams;
  PS ps = opflow->ps;
  int i;
  double obj_factor = opflow->obj_factor;
  int isobj_gencost = opflow->obj_gencost;
  double MVAbase = ps->MVAbase;
  double powerimbal_penalty = opflow->powerimbalance_penalty;
  PetscInt flps = 0;
  int loc;

  if (iHSS != NULL && jHSS != NULL) {

    /* Generator contributions */
    for (i = 0; i < genparams->ngenON; i++) {
      loc = genparams->hesssp_idx[i];
      iHSS[loc] = genparams->xidx[i];
      jHSS[loc] = genparams->xidx[i];
    }
    /* Loadloss contributions */
    if (opflow->include_loadloss_variables) {
      for (i = 0; i < loadparams->nload; i++) {
        loc = loadparams->hesssp_idx[i];
        iHSS[loc] = loadparams->xidx[i];
        jHSS[loc] = loadparams->xidx[i];
        iHSS[loc + 1] = loadparams->xidx[i] + 1;
        jHSS[loc + 1] = loadparams->xidx[i] + 1;
      }
    }
  }

  if (MHSS != NULL) {

    if (opflow->objectivetype == MIN_GEN_COST) {
      /* Generator contributions */
      for (i = 0; i < genparams->ngenON; i++) {
        loc = genparams->hesssp_idx[i];
        MHSS[loc] = isobj_gencost * obj_factor * 2.0 *
                    genparams->cost_alpha[i] * MVAbase * MVAbase;
      }
    } else if (opflow->objectivetype == NO_OBJ) {
      for (i = 0; i < genparams->ngenON; i++) {
        loc = genparams->hesssp_idx[i];
        MHSS[loc] = 0.0;
      }
    }

    /* Loadloss contributions */
    if (opflow->include_loadloss_variables) {
      for (i = 0; i < loadparams->nload; i++) {
        loc = loadparams->hesssp_idx[i];
        MHSS[loc] = 0.0;
        MHSS[loc + 1] = 0.0;
      }
    }
    flps += 5 * genparams->ngenON;
  }

  ierr = PetscLogFlops(flps);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeDenseEqualityConstraintJacobian_PBPOLHIOP(
    OPFLOW opflow, const double *x, double *JacD) {
  int i, j, row[2], col[4];
  PBPOLHIOP pbpolhiop = (PBPOLHIOP)opflow->model;
  BUSParams *busparams = &pbpolhiop->busparams;
  LINEParams *lineparams = &pbpolhiop->lineparams;
  double val[8];
  int nxsparse = opflow->nxsparse;
  int nxdense = opflow->nxdense;
  const int JacDnrows = opflow->nconeq;
  PetscInt flps = 0;
  PetscErrorCode ierr;
  PetscFunctionBegin;

  if (JacD == NULL)
    PetscFunctionReturn(0);

  /* Zero out JacD */
  for (i = 0; i < JacDnrows; i++) {
    for (j = 0; j < nxdense; j++) {
      JacD[(i * nxdense) + j] = 0.0;
    }
  }

  /* Jacobian from bus contributions */
  for (i = 0; i < busparams->nbus; i++) {
    double Vm = x[busparams->xidx[i] + 1];
    row[0] = busparams->gidx[i];
    row[1] = busparams->gidx[i] + 1;

    col[0] = busparams->xidx[i] - nxsparse;
    col[1] = busparams->xidx[i] + 1 - nxsparse;

    val[0] = busparams->isisolated[i] * 1.0 + busparams->ispvpq[i] * 0.0;
    val[1] = busparams->isisolated[i] * 0.0 +
             busparams->ispvpq[i] * 2 * Vm * busparams->gl[i];
    val[2] = 0.0;
    val[3] = busparams->isisolated[i] * 1.0 +
             busparams->ispvpq[i] * -2 * Vm * busparams->bl[i];

    JacD[(nxdense * row[0]) + col[0]] += val[0];
    JacD[(nxdense * row[0]) + col[1]] += val[1];
    JacD[(nxdense * row[1]) + col[0]] += val[2];
    JacD[(nxdense * row[1]) + col[1]] += val[3];
  }
  flps += busparams->nbus * 15;

  /* Jacobian from line contributions */
  for (i = 0; i < lineparams->nlineON; i++) {
    double thetaf = x[lineparams->xidxf[i]], Vmf = x[lineparams->xidxf[i] + 1];
    double thetat = x[lineparams->xidxt[i]], Vmt = x[lineparams->xidxt[i] + 1];
    double thetaft = thetaf - thetat;
    double thetatf = thetat - thetaf;

    row[0] = lineparams->geqidxf[i];
    row[1] = lineparams->geqidxf[i] + 1;

    col[0] = lineparams->xidxf[i] - nxsparse;
    col[1] = lineparams->xidxf[i] + 1 - nxsparse;
    col[2] = lineparams->xidxt[i] - nxsparse;
    col[3] = lineparams->xidxt[i] + 1 - nxsparse;

    /* dPf_dthetaf */
    val[0] = Vmf * Vmt *
             (-lineparams->Gft[i] * sin(thetaft) +
              lineparams->Bft[i] * cos(thetaft));
    /*dPf_dVmf */
    val[1] = 2 * lineparams->Gff[i] * Vmf +
             Vmt * (lineparams->Gft[i] * cos(thetaft) +
                    lineparams->Bft[i] * sin(thetaft));
    /*dPf_dthetat */
    val[2] =
        Vmf * Vmt *
        (lineparams->Gft[i] * sin(thetaft) - lineparams->Bft[i] * cos(thetaft));
    /* dPf_dVmt */
    val[3] = Vmf * (lineparams->Gft[i] * cos(thetaft) +
                    lineparams->Bft[i] * sin(thetaft));

    JacD[(nxdense * row[0]) + col[0]] += val[0];
    JacD[(nxdense * row[0]) + col[1]] += val[1];
    JacD[(nxdense * row[0]) + col[2]] += val[2];
    JacD[(nxdense * row[0]) + col[3]] += val[3];

    /* dQf_dthetaf */
    val[4] =
        Vmf * Vmt *
        (lineparams->Bft[i] * sin(thetaft) + lineparams->Gft[i] * cos(thetaft));
    /* dQf_dVmf */
    val[5] = -2 * lineparams->Bff[i] * Vmf +
             Vmt * (-lineparams->Bft[i] * cos(thetaft) +
                    lineparams->Gft[i] * sin(thetaft));
    /* dQf_dthetat */
    val[6] = Vmf * Vmt *
             (-lineparams->Bft[i] * sin(thetaft) -
              lineparams->Gft[i] * cos(thetaft));
    /* dQf_dVmt */
    val[7] = Vmf * (-lineparams->Bft[i] * cos(thetaft) +
                    lineparams->Gft[i] * sin(thetaft));

    JacD[(nxdense * row[1]) + col[0]] += val[4];
    JacD[(nxdense * row[1]) + col[1]] += val[5];
    JacD[(nxdense * row[1]) + col[2]] += val[6];
    JacD[(nxdense * row[1]) + col[3]] += val[7];

    row[0] = lineparams->geqidxt[i];
    row[1] = lineparams->geqidxt[i] + 1;

    col[0] = lineparams->xidxt[i] - nxsparse;
    col[1] = lineparams->xidxt[i] + 1 - nxsparse;
    col[2] = lineparams->xidxf[i] - nxsparse;
    col[3] = lineparams->xidxf[i] + 1 - nxsparse;

    /* dPt_dthetat */
    val[0] = Vmt * Vmf *
             (-lineparams->Gtf[i] * sin(thetatf) +
              lineparams->Btf[i] * cos(thetatf));
    /* dPt_dVmt */
    val[1] = 2 * lineparams->Gtt[i] * Vmt +
             Vmf * (lineparams->Gtf[i] * cos(thetatf) +
                    lineparams->Btf[i] * sin(thetatf));
    /* dPt_dthetaf */
    val[2] =
        Vmt * Vmf *
        (lineparams->Gtf[i] * sin(thetatf) - lineparams->Btf[i] * cos(thetatf));
    /* dPt_dVmf */
    val[3] = Vmt * (lineparams->Gtf[i] * cos(thetatf) +
                    lineparams->Btf[i] * sin(thetatf));

    JacD[(nxdense * row[0]) + col[0]] += val[0];
    JacD[(nxdense * row[0]) + col[1]] += val[1];
    JacD[(nxdense * row[0]) + col[2]] += val[2];
    JacD[(nxdense * row[0]) + col[3]] += val[3];

    /* dQt_dthetat */
    val[4] =
        Vmt * Vmf *
        (lineparams->Btf[i] * sin(thetatf) + lineparams->Gtf[i] * cos(thetatf));
    /* dQt_dVmt */
    val[5] = -2 * lineparams->Btt[i] * Vmt +
             Vmf * (-lineparams->Btf[i] * cos(thetatf) +
                    lineparams->Gtf[i] * sin(thetatf));
    /* dQt_dthetaf */
    val[6] = Vmt * Vmf *
             (-lineparams->Btf[i] * sin(thetatf) -
              lineparams->Gtf[i] * cos(thetatf));
    /* dQt_dVmf */
    val[7] = Vmt * (-lineparams->Btf[i] * cos(thetatf) +
                    lineparams->Gtf[i] * sin(thetatf));

    JacD[(nxdense * row[1]) + col[0]] += val[4];
    JacD[(nxdense * row[1]) + col[1]] += val[5];
    JacD[(nxdense * row[1]) + col[2]] += val[6];
    JacD[(nxdense * row[1]) + col[3]] += val[7];
  }
  flps += (188 + (16 * EXAGO_FLOPS_COSOP) + (16 * EXAGO_FLOPS_SINOP)) *
          lineparams->nlineON;

  ierr = PetscLogFlops(flps);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeDenseInequalityConstraintJacobian_PBPOLHIOP(
    OPFLOW opflow, const double *x, double *JacD) {
  PetscErrorCode ierr;
  int i, k;
  PBPOLHIOP pbpolhiop = (PBPOLHIOP)opflow->model;
  LINEParams *lineparams = &pbpolhiop->lineparams;
  int row[2], col[4];
  double val[4];
  int nxsparse = opflow->nxsparse;
  int nxdense = opflow->nxdense;
  int JacDnrows = opflow->nconineq;
  PetscInt flps = 0;

  PetscFunctionBegin;

  if (JacD == NULL)
    PetscFunctionReturn(0);

  for (i = 0; i < JacDnrows; i++) {
    for (k = 0; k < nxdense; k++) {
      JacD[(i * nxdense) + k] = 0.0;
    }
  }

  for (i = 0; i < lineparams->nlinelim; i++) {
    int j = lineparams->linelimidx[i];

    double Pf, Qf, Pt, Qt;
    double thetaf = x[lineparams->xidxf[j]], Vmf = x[lineparams->xidxf[j] + 1];
    double thetat = x[lineparams->xidxt[j]], Vmt = x[lineparams->xidxt[j] + 1];
    double thetaft = thetaf - thetat;
    double thetatf = thetat - thetaf;
    double dSf2_dPf, dSf2_dQf, dSt2_dPt, dSt2_dQt;
    double dPf_dthetaf, dPf_dVmf, dPf_dthetat, dPf_dVmt;
    double dQf_dthetaf, dQf_dVmf, dQf_dthetat, dQf_dVmt;
    double dPt_dthetaf, dPt_dVmf, dPt_dthetat, dPt_dVmt;
    double dQt_dthetaf, dQt_dVmf, dQt_dthetat, dQt_dVmt;
    double dSf2_dthetaf, dSf2_dVmf, dSf2_dthetat, dSf2_dVmt;
    double dSt2_dthetaf, dSt2_dVmf, dSt2_dthetat, dSt2_dVmt;
    double Gff = lineparams->Gff[j], Bff = lineparams->Bff[j];
    double Gft = lineparams->Gft[j], Bft = lineparams->Bft[j];
    double Gtf = lineparams->Gtf[j], Btf = lineparams->Btf[j];
    double Gtt = lineparams->Gtt[j], Btt = lineparams->Btt[j];

    Pf =
        Gff * Vmf * Vmf + Vmf * Vmt * (Gft * cos(thetaft) + Bft * sin(thetaft));
    Qf = -Bff * Vmf * Vmf +
         Vmf * Vmt * (-Bft * cos(thetaft) + Gft * sin(thetaft));
    Pt =
        Gtt * Vmt * Vmt + Vmt * Vmf * (Gtf * cos(thetatf) + Btf * sin(thetatf));
    Qt = -Btt * Vmt * Vmt +
         Vmt * Vmf * (-Btf * cos(thetatf) + Gtf * sin(thetatf));

    dSf2_dPf = 2 * Pf;
    dSf2_dQf = 2 * Qf;
    dSt2_dPt = 2 * Pt;
    dSt2_dQt = 2 * Qt;

    dPf_dthetaf = Vmf * Vmt * (-Gft * sin(thetaft) + Bft * cos(thetaft));
    dPf_dVmf = 2 * Gff * Vmf + Vmt * (Gft * cos(thetaft) + Bft * sin(thetaft));
    dPf_dthetat = Vmf * Vmt * (Gft * sin(thetaft) - Bft * cos(thetaft));
    dPf_dVmt = Vmf * (Gft * cos(thetaft) + Bft * sin(thetaft));

    dQf_dthetaf = Vmf * Vmt * (Bft * sin(thetaft) + Gft * cos(thetaft));
    dQf_dVmf =
        -2 * Bff * Vmf + Vmt * (-Bft * cos(thetaft) + Gft * sin(thetaft));
    dQf_dthetat = Vmf * Vmt * (-Bft * sin(thetaft) - Gft * cos(thetaft));
    dQf_dVmt = Vmf * (-Bft * cos(thetaft) + Gft * sin(thetaft));

    dPt_dthetat = Vmt * Vmf * (-Gtf * sin(thetatf) + Btf * cos(thetatf));
    dPt_dVmt = 2 * Gtt * Vmt + Vmf * (Gtf * cos(thetatf) + Btf * sin(thetatf));
    dPt_dthetaf = Vmt * Vmf * (Gtf * sin(thetatf) - Btf * cos(thetatf));
    dPt_dVmf = Vmt * (Gtf * cos(thetatf) + Btf * sin(thetatf));

    dQt_dthetat = Vmt * Vmf * (Btf * sin(thetatf) + Gtf * cos(thetatf));
    dQt_dVmt =
        -2 * Btt * Vmt + Vmf * (-Btf * cos(thetatf) + Gtf * sin(thetatf));
    dQt_dthetaf = Vmt * Vmf * (-Btf * sin(thetatf) - Gtf * cos(thetatf));
    dQt_dVmf = Vmt * (-Btf * cos(thetatf) + Gtf * sin(thetatf));

    dSf2_dthetaf = dSf2_dPf * dPf_dthetaf + dSf2_dQf * dQf_dthetaf;
    dSf2_dthetat = dSf2_dPf * dPf_dthetat + dSf2_dQf * dQf_dthetat;
    dSf2_dVmf = dSf2_dPf * dPf_dVmf + dSf2_dQf * dQf_dVmf;
    dSf2_dVmt = dSf2_dPf * dPf_dVmt + dSf2_dQf * dQf_dVmt;

    row[0] = lineparams->gineqidx[i];

    col[0] = lineparams->xidxf[j] - nxsparse;
    col[1] = lineparams->xidxf[j] + 1 - nxsparse;
    col[2] = lineparams->xidxt[j] - nxsparse;
    col[3] = lineparams->xidxt[j] + 1 - nxsparse;

    val[0] = dSf2_dthetaf;
    val[1] = dSf2_dVmf;
    val[2] = dSf2_dthetat;
    val[3] = dSf2_dVmt;

    JacD[(nxdense * row[0]) + col[0]] += val[0];
    JacD[(nxdense * row[0]) + col[1]] += val[1];
    JacD[(nxdense * row[0]) + col[2]] += val[2];
    JacD[(nxdense * row[0]) + col[3]] += val[3];

    dSt2_dthetaf = dSt2_dPt * dPt_dthetaf + dSt2_dQt * dQt_dthetaf;
    dSt2_dthetat = dSt2_dPt * dPt_dthetat + dSt2_dQt * dQt_dthetat;
    dSt2_dVmf = dSt2_dPt * dPt_dVmf + dSt2_dQt * dQt_dVmf;
    dSt2_dVmt = dSt2_dPt * dPt_dVmt + dSt2_dQt * dQt_dVmt;

    row[0] = lineparams->gineqidx[i] + 1;

    col[0] = lineparams->xidxt[j] - nxsparse;
    col[1] = lineparams->xidxt[j] + 1 - nxsparse;
    col[2] = lineparams->xidxf[j] - nxsparse;
    col[3] = lineparams->xidxf[j] + 1 - nxsparse;

    val[0] = dSt2_dthetat;
    val[1] = dSt2_dVmt;
    val[2] = dSt2_dthetaf;
    val[3] = dSt2_dVmf;

    JacD[(nxdense * row[0]) + col[0]] += val[0];
    JacD[(nxdense * row[0]) + col[1]] += val[1];
    JacD[(nxdense * row[0]) + col[2]] += val[2];
    JacD[(nxdense * row[0]) + col[3]] += val[3];
  }

  flps += (206 + (20 * EXAGO_FLOPS_COSOP) + (20 * EXAGO_FLOPS_SINOP)) *
          lineparams->nlinelim;
  ierr = PetscLogFlops(flps);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeDenseEqualityConstraintHessian_PBPOLHIOP(
    OPFLOW opflow, const double *x, const double *lambda, double *HDD) {
  PetscErrorCode ierr;
  PBPOLHIOP pbpolhiop = (PBPOLHIOP)opflow->model;
  BUSParams *busparams = &pbpolhiop->busparams;
  LINEParams *lineparams = &pbpolhiop->lineparams;
  PS ps = opflow->ps;
  int i;
  int row[16], col[16];
  double val[16];
  int nxsparse = opflow->nxsparse;
  int nxdense = opflow->nxdense;
  double obj_factor = opflow->obj_factor;
  double MVAbase = ps->MVAbase;
  PetscInt flps = 0;

  PetscFunctionBegin;

  /* Hessian from bus contributions */
  for (i = 0; i < busparams->nbus; i++) {
    row[0] = busparams->xidx[i] + 1 - nxsparse;
    col[0] = row[0];
    val[0] = busparams->ispvpq[i] *
             (lambda[busparams->gidx[i]] * 2 * busparams->gl[i] +
              lambda[busparams->gidx[i] + 1] * (-2 * busparams->bl[i]));
    HDD[(nxdense * row[0]) + col[0]] += val[0];
  }
  flps += 10 * busparams->nbus;

  /* Hessian from line contributions */
  for (i = 0; i < lineparams->nlineON; i++) {
    int gloc;
    double Gff, Bff, Gft, Bft, Gtf, Btf, Gtt, Btt;
    Gff = lineparams->Gff[i];
    Bff = lineparams->Bff[i];
    Gft = lineparams->Gft[i];
    Bft = lineparams->Bft[i];
    Gtf = lineparams->Gtf[i];
    Btf = lineparams->Btf[i];
    Gtt = lineparams->Gtt[i];
    Btt = lineparams->Btt[i];

    double thetaf = x[lineparams->xidxf[i]], Vmf = x[lineparams->xidxf[i] + 1];
    double thetat = x[lineparams->xidxt[i]], Vmt = x[lineparams->xidxt[i] + 1];
    double thetaft = thetaf - thetat;
    double thetatf = thetat - thetaf;

    double dPf_dthetaf_dthetaf, dPf_dthetaf_dVmf, dPf_dthetaf_dthetat,
        dPf_dthetaf_dVmt;
    double dPf_dVmf_dthetaf, dPf_dVmf_dVmf, dPf_dVmf_dthetat, dPf_dVmf_dVmt;
    double dPf_dthetat_dthetaf, dPf_dthetat_dVmf, dPf_dthetat_dthetat,
        dPf_dthetat_dVmt;
    double dPf_dVmt_dthetaf, dPf_dVmt_dVmf, dPf_dVmt_dthetat, dPf_dVmt_dVmt;

    /* dPf_dthetaf = Vmf*Vmt*(-Gft*sin(thetaft) + Bft*cos(thetaft)); */
    dPf_dthetaf_dthetaf =
        -Vmf * Vmt * (Gft * cos(thetaft) + Bft * sin(thetaft));
    dPf_dthetaf_dVmf = Vmt * (-Gft * sin(thetaft) + Bft * cos(thetaft));
    dPf_dthetaf_dthetat = Vmf * Vmt * (Gft * cos(thetaft) + Bft * sin(thetaft));
    dPf_dthetaf_dVmt = Vmf * (-Gft * sin(thetaft) + Bft * cos(thetaft));

    /* dPf_Vmf  = 2*Gff*Vmf + Vmt*(Gft*cos(thetaft) + Bft*sin(thetaft)); */
    dPf_dVmf_dthetaf = Vmt * (-Gft * sin(thetaft) + Bft * cos(thetaft));
    dPf_dVmf_dVmf = 2 * Gff;
    dPf_dVmf_dthetat = Vmt * (Gft * sin(thetaft) - Bft * cos(thetaft));
    dPf_dVmf_dVmt = (Gft * cos(thetaft) + Bft * sin(thetaft));

    /* dPf_dthetat = Vmf*Vmt*(Gft*sin(thetaft) - Bft*cos(thetaft)); */
    dPf_dthetat_dthetaf = Vmf * Vmt * (Gft * cos(thetaft) + Bft * sin(thetaft));
    dPf_dthetat_dVmf = Vmt * (Gft * sin(thetaft) - Bft * cos(thetaft));
    dPf_dthetat_dthetat =
        Vmf * Vmt * (-Gft * cos(thetaft) - Bft * sin(thetaft));
    dPf_dthetat_dVmt = Vmf * (Gft * sin(thetaft) - Bft * cos(thetaft));

    /* dPf_dVmt = Vmf*(Gft*cos(thetaft) + Bft*sin(thetaft)); */
    dPf_dVmt_dthetaf = Vmf * (-Gft * sin(thetaft) + Bft * cos(thetaft));
    dPf_dVmt_dVmf = (Gft * cos(thetaft) + Bft * sin(thetaft));
    dPf_dVmt_dthetat = Vmf * (Gft * sin(thetaft) - Bft * cos(thetaft));
    dPf_dVmt_dVmt = 0.0;

    double dQf_dthetaf_dthetaf, dQf_dthetaf_dVmf, dQf_dthetaf_dthetat,
        dQf_dthetaf_dVmt;
    double dQf_dVmf_dthetaf, dQf_dVmf_dVmf, dQf_dVmf_dthetat, dQf_dVmf_dVmt;
    double dQf_dthetat_dthetaf, dQf_dthetat_dVmf, dQf_dthetat_dthetat,
        dQf_dthetat_dVmt;
    double dQf_dVmt_dthetaf, dQf_dVmt_dVmf, dQf_dVmt_dthetat, dQf_dVmt_dVmt;

    /* dQf_dthetaf = Vmf*Vmt*(Bft*sin(thetaft) + Gft*cos(thetaft)); */
    dQf_dthetaf_dthetaf = Vmf * Vmt * (Bft * cos(thetaft) - Gft * sin(thetaft));
    dQf_dthetaf_dVmf = Vmt * (Bft * sin(thetaft) + Gft * cos(thetaft));
    dQf_dthetaf_dthetat =
        Vmf * Vmt * (-Bft * cos(thetaft) + Gft * sin(thetaft));
    dQf_dthetaf_dVmt = Vmf * (Bft * sin(thetaft) + Gft * cos(thetaft));

    /* dQf_dVmf = -2*Bff*Vmf + Vmt*(-Bft*cos(thetaft) + Gft*sin(thetaft)); */
    dQf_dVmf_dthetaf = Vmt * (Bft * sin(thetaft) + Gft * cos(thetaft));
    dQf_dVmf_dVmf = -2 * Bff;
    dQf_dVmf_dthetat = Vmt * (-Bft * sin(thetaft) - Gft * cos(thetaft));
    dQf_dVmf_dVmt = (-Bft * cos(thetaft) + Gft * sin(thetaft));

    /* dQf_dthetat = Vmf*Vmt*(-Bft*sin(thetaft) - Gft*cos(thetaft)); */
    dQf_dthetat_dthetaf =
        Vmf * Vmt * (-Bft * cos(thetaft) + Gft * sin(thetaft));
    dQf_dthetat_dVmf = Vmt * (-Bft * sin(thetaft) - Gft * cos(thetaft));
    dQf_dthetat_dthetat = Vmf * Vmt * (Bft * cos(thetaft) - Gft * sin(thetaft));
    dQf_dthetat_dVmt = Vmf * (-Bft * sin(thetaft) - Gft * cos(thetaft));

    /* dQf_dVmt = Vmf*(-Bft*cos(thetaft) + Gft*sin(thetaft)); */
    dQf_dVmt_dthetaf = Vmf * (Bft * sin(thetaft) + Gft * cos(thetaft));
    dQf_dVmt_dVmf = (-Bft * cos(thetaft) + Gft * sin(thetaft));
    dQf_dVmt_dthetat = Vmf * (-Bft * sin(thetaft) - Gft * cos(thetaft));
    dQf_dVmt_dVmt = 0.0;

    row[0] = lineparams->xidxf[i] - nxsparse;
    row[1] = lineparams->xidxf[i] + 1 - nxsparse;
    col[0] = lineparams->xidxf[i] - nxsparse;
    col[1] = lineparams->xidxf[i] + 1 - nxsparse;
    col[2] = lineparams->xidxt[i] - nxsparse;
    col[3] = lineparams->xidxt[i] + 1 - nxsparse;

    gloc = lineparams->geqidxf[i];

    val[0] = val[1] = val[2] = val[3] = val[4] = val[5] = val[6] = val[7] = 0.0;

    val[0] = lambda[gloc] * dPf_dthetaf_dthetaf +
             lambda[gloc + 1] * dQf_dthetaf_dthetaf;
    val[1] =
        lambda[gloc] * dPf_dthetaf_dVmf + lambda[gloc + 1] * dQf_dthetaf_dVmf;
    val[2] = lambda[gloc] * dPf_dthetaf_dthetat +
             lambda[gloc + 1] * dQf_dthetaf_dthetat;
    val[3] =
        lambda[gloc] * dPf_dthetaf_dVmt + lambda[gloc + 1] * dQf_dthetaf_dVmt;

    HDD[(nxdense * row[0]) + col[0]] += val[0];
    HDD[(nxdense * row[0]) + col[1]] += val[1];
    HDD[(nxdense * row[0]) + col[2]] += val[2];
    HDD[(nxdense * row[0]) + col[3]] += val[3];

    val[4] =
        lambda[gloc] * dPf_dVmf_dthetaf + lambda[gloc + 1] * dQf_dVmf_dthetaf;
    val[5] = lambda[gloc] * dPf_dVmf_dVmf + lambda[gloc + 1] * dQf_dVmf_dVmf;
    val[6] =
        lambda[gloc] * dPf_dVmf_dthetat + lambda[gloc + 1] * dQf_dVmf_dthetat;
    val[7] = lambda[gloc] * dPf_dVmf_dVmt + lambda[gloc + 1] * dQf_dVmf_dVmt;

    HDD[(nxdense * row[1]) + col[0]] += val[4];
    HDD[(nxdense * row[1]) + col[1]] += val[5];
    HDD[(nxdense * row[1]) + col[2]] += val[6];
    HDD[(nxdense * row[1]) + col[3]] += val[7];

    row[0] = lineparams->xidxt[i] - nxsparse;
    row[1] = lineparams->xidxt[i] + 1 - nxsparse;

    col[0] = lineparams->xidxf[i] - nxsparse;
    col[1] = lineparams->xidxf[i] + 1 - nxsparse;
    col[2] = lineparams->xidxt[i] - nxsparse;
    col[3] = lineparams->xidxt[i] + 1 - nxsparse;

    val[0] = val[1] = val[2] = val[3] = val[4] = val[5] = val[6] = val[7] = 0.0;

    val[0] = lambda[gloc] * dPf_dthetat_dthetaf +
             lambda[gloc + 1] * dQf_dthetat_dthetaf;
    val[1] =
        lambda[gloc] * dPf_dthetat_dVmf + lambda[gloc + 1] * dQf_dthetat_dVmf;
    val[2] = lambda[gloc] * dPf_dthetat_dthetat +
             lambda[gloc + 1] * dQf_dthetat_dthetat;
    val[3] =
        lambda[gloc] * dPf_dthetat_dVmt + lambda[gloc + 1] * dQf_dthetat_dVmt;

    HDD[(nxdense * row[0]) + col[0]] += val[0];
    HDD[(nxdense * row[0]) + col[1]] += val[1];
    HDD[(nxdense * row[0]) + col[2]] += val[2];
    HDD[(nxdense * row[0]) + col[3]] += val[3];

    val[4] =
        lambda[gloc] * dPf_dVmt_dthetaf + lambda[gloc + 1] * dQf_dVmt_dthetaf;
    val[5] = lambda[gloc] * dPf_dVmt_dVmf + lambda[gloc + 1] * dQf_dVmt_dVmf;
    val[6] =
        lambda[gloc] * dPf_dVmt_dthetat + lambda[gloc + 1] * dQf_dVmt_dthetat;
    val[7] = lambda[gloc] * dPf_dVmt_dVmt + lambda[gloc + 1] * dQf_dVmt_dVmt;

    HDD[(nxdense * row[1]) + col[0]] += val[4];
    HDD[(nxdense * row[1]) + col[1]] += val[5];
    HDD[(nxdense * row[1]) + col[2]] += val[6];
    HDD[(nxdense * row[1]) + col[3]] += val[7];

    double dPt_dthetat_dthetat, dPt_dthetat_dVmt, dPt_dthetat_dthetaf,
        dPt_dthetat_dVmf;
    double dPt_dVmt_dthetat, dPt_dVmt_dVmt, dPt_dVmt_dthetaf, dPt_dVmt_dVmf;
    double dPt_dthetaf_dthetat, dPt_dthetaf_dVmt, dPt_dthetaf_dthetaf,
        dPt_dthetaf_dVmf;
    double dPt_dVmf_dthetat, dPt_dVmf_dVmt, dPt_dVmf_dthetaf, dPt_dVmf_dVmf;

    /* dPt_dthetat = Vmf*Vmt*(-Gtf*sin(thetatf) + Btf*cos(thetatf)); */
    dPt_dthetat_dthetat =
        Vmf * Vmt * (-Gtf * cos(thetatf) - Btf * sin(thetatf));
    dPt_dthetat_dVmt = Vmf * (-Gtf * sin(thetatf) + Btf * cos(thetatf));
    dPt_dthetat_dthetaf = Vmf * Vmt * (Gtf * cos(thetatf) + Btf * sin(thetatf));
    dPt_dthetat_dVmf = Vmt * (-Gtf * sin(thetatf) + Btf * cos(thetatf));

    /* dPt_Vmt  = 2*Gtt*Vmt + Vmf*(Gtf*cos(thetatf) + Btf*sin(thetatf)); */
    dPt_dVmt_dthetat = Vmf * (-Gtf * sin(thetatf) + Bft * cos(thetatf));
    dPt_dVmt_dVmt = 2 * Gtt;
    dPt_dVmt_dthetaf = Vmf * (Gtf * sin(thetatf) - Btf * cos(thetatf));
    dPt_dVmt_dVmf = (Gtf * cos(thetatf) + Btf * sin(thetatf));

    /* dPt_dthetaf = Vmf*Vmt*(Gtf*sin(thetatf) - Btf*cos(thetatf)); */
    dPt_dthetaf_dthetat = Vmf * Vmt * (Gtf * cos(thetatf) + Btf * sin(thetatf));
    dPt_dthetaf_dVmt = Vmf * (Gtf * sin(thetatf) - Btf * cos(thetatf));
    dPt_dthetaf_dthetaf =
        Vmf * Vmt * (-Gtf * cos(thetatf) - Btf * sin(thetatf));
    dPt_dthetaf_dVmf = Vmt * (Gtf * sin(thetatf) - Btf * cos(thetatf));

    /* dPt_dVmf = Vmt*(Gtf*cos(thetatf) + Btf*sin(thetatf)); */
    dPt_dVmf_dthetat = Vmt * (-Gtf * sin(thetatf) + Btf * cos(thetatf));
    dPt_dVmf_dVmt = (Gtf * cos(thetatf) + Btf * sin(thetatf));
    dPt_dVmf_dthetaf = Vmt * (Gtf * sin(thetatf) - Btf * cos(thetatf));
    dPt_dVmf_dVmf = 0.0;

    double dQt_dthetaf_dthetaf, dQt_dthetaf_dVmf, dQt_dthetaf_dthetat,
        dQt_dthetaf_dVmt;
    double dQt_dVmf_dthetaf, dQt_dVmf_dVmf, dQt_dVmf_dthetat, dQt_dVmf_dVmt;
    double dQt_dthetat_dthetaf, dQt_dthetat_dVmf, dQt_dthetat_dthetat,
        dQt_dthetat_dVmt;
    double dQt_dVmt_dthetaf, dQt_dVmt_dVmf, dQt_dVmt_dthetat, dQt_dVmt_dVmt;

    /* dQt_dthetat = Vmf*Vmt*(Btf*sin(thetatf) + Gtf*cos(thetatf)); */
    dQt_dthetat_dthetat = Vmf * Vmt * (Btf * cos(thetatf) - Gtf * sin(thetatf));
    dQt_dthetat_dVmt = Vmf * (Btf * sin(thetatf) + Gtf * cos(thetatf));
    dQt_dthetat_dthetaf =
        Vmf * Vmt * (-Btf * cos(thetatf) + Gtf * sin(thetatf));
    dQt_dthetat_dVmf = Vmf * (Btf * sin(thetatf) + Gtf * cos(thetatf));

    /* dQt_dVmt = -2*Btt*Vmt + Vmf*(-Btf*cos(thetatf) + Gtf*sin(thetatf)); */
    dQt_dVmt_dthetat = Vmf * (Btf * sin(thetatf) + Gtf * cos(thetatf));
    dQt_dVmt_dVmt = -2 * Btt;
    dQt_dVmt_dthetaf = Vmf * (-Btf * sin(thetatf) + Gtf * cos(thetatf));
    dQt_dVmt_dVmf = (-Btf * cos(thetatf) + Gtf * sin(thetatf));

    /* dQt_dthetaf = Vmf*Vmt*(-Btf*sin(thetatf) - Gtf*cos(thetatf)); */
    dQt_dthetaf_dthetat =
        Vmf * Vmt * (-Btf * cos(thetatf) + Gtf * sin(thetatf));
    dQt_dthetaf_dVmt = Vmf * (-Btf * sin(thetatf) - Gtf * cos(thetatf));
    dQt_dthetaf_dthetaf = Vmf * Vmt * (Btf * cos(thetatf) - Gtf * sin(thetatf));
    dQt_dthetaf_dVmf = Vmt * (-Btf * sin(thetatf) - Gtf * cos(thetatf));

    /* dQt_dVmf = Vmt*(-Btf*cos(thetatf) + Gtf*sin(thetatf)); */
    dQt_dVmf_dthetat = Vmt * (Btf * sin(thetatf) + Gtf * cos(thetatf));
    dQt_dVmf_dVmt = (-Btf * cos(thetatf) + Gtf * sin(thetatf));
    dQt_dVmf_dthetaf = Vmt * (-Btf * sin(thetatf) - Gtf * cos(thetatf));
    dQt_dVmf_dVmf = 0.0;

    row[0] = lineparams->xidxt[i] - nxsparse;
    row[1] = lineparams->xidxt[i] + 1 - nxsparse;
    col[0] = lineparams->xidxt[i] - nxsparse;
    col[1] = lineparams->xidxt[i] + 1 - nxsparse;
    col[2] = lineparams->xidxf[i] - nxsparse;
    col[3] = lineparams->xidxf[i] + 1 - nxsparse;

    gloc = lineparams->geqidxt[i];

    val[0] = val[1] = val[2] = val[3] = val[4] = val[5] = val[6] = val[7] = 0.0;

    val[0] = lambda[gloc] * dPt_dthetat_dthetat +
             lambda[gloc + 1] * dQt_dthetat_dthetat;
    val[1] =
        lambda[gloc] * dPt_dthetat_dVmt + lambda[gloc + 1] * dQt_dthetat_dVmt;
    val[2] = lambda[gloc] * dPt_dthetat_dthetaf +
             lambda[gloc + 1] * dQt_dthetat_dthetaf;
    val[3] =
        lambda[gloc] * dPt_dthetat_dVmf + lambda[gloc + 1] * dQt_dthetat_dVmf;

    HDD[(nxdense * row[0]) + col[0]] += val[0];
    HDD[(nxdense * row[0]) + col[1]] += val[1];
    HDD[(nxdense * row[0]) + col[2]] += val[2];
    HDD[(nxdense * row[0]) + col[3]] += val[3];

    val[4] =
        lambda[gloc] * dPt_dVmt_dthetat + lambda[gloc + 1] * dQt_dVmt_dthetat;
    val[5] = lambda[gloc] * dPt_dVmt_dVmt + lambda[gloc + 1] * dQt_dVmt_dVmt;
    val[6] =
        lambda[gloc] * dPt_dVmt_dthetaf + lambda[gloc + 1] * dQt_dVmt_dthetaf;
    val[7] = lambda[gloc] * dPt_dVmt_dVmf + lambda[gloc + 1] * dQt_dVmt_dVmf;

    HDD[(nxdense * row[1]) + col[0]] += val[4];
    HDD[(nxdense * row[1]) + col[1]] += val[5];
    HDD[(nxdense * row[1]) + col[2]] += val[6];
    HDD[(nxdense * row[1]) + col[3]] += val[7];

    row[0] = lineparams->xidxf[i] - nxsparse;
    row[1] = lineparams->xidxf[i] + 1 - nxsparse;
    col[0] = lineparams->xidxt[i] - nxsparse;
    col[1] = lineparams->xidxt[i] + 1 - nxsparse;
    col[2] = lineparams->xidxf[i] - nxsparse;
    col[3] = lineparams->xidxf[i] + 1 - nxsparse;

    val[0] = val[1] = val[2] = val[3] = val[4] = val[5] = val[6] = val[7] = 0.0;

    val[0] = lambda[gloc] * dPt_dthetaf_dthetat +
             lambda[gloc + 1] * dQt_dthetaf_dthetat;
    val[1] =
        lambda[gloc] * dPt_dthetaf_dVmt + lambda[gloc + 1] * dQt_dthetaf_dVmt;
    val[2] = lambda[gloc] * dPt_dthetaf_dthetaf +
             lambda[gloc + 1] * dQt_dthetaf_dthetaf;
    val[3] =
        lambda[gloc] * dPt_dthetaf_dVmf + lambda[gloc + 1] * dQt_dthetaf_dVmf;

    HDD[(nxdense * row[0]) + col[0]] += val[0];
    HDD[(nxdense * row[0]) + col[1]] += val[1];
    HDD[(nxdense * row[0]) + col[2]] += val[2];
    HDD[(nxdense * row[0]) + col[3]] += val[3];

    val[4] =
        lambda[gloc] * dPt_dVmf_dthetat + lambda[gloc + 1] * dQt_dVmf_dthetat;
    val[5] = lambda[gloc] * dPt_dVmf_dVmt + lambda[gloc + 1] * dQt_dVmf_dVmt;
    val[6] =
        lambda[gloc] * dPt_dVmf_dthetaf + lambda[gloc + 1] * dQt_dVmf_dthetaf;
    val[7] = lambda[gloc] * dPt_dVmf_dVmf + lambda[gloc + 1] * dQt_dVmf_dVmf;

    HDD[(nxdense * row[1]) + col[0]] += val[4];
    HDD[(nxdense * row[1]) + col[1]] += val[5];
    HDD[(nxdense * row[1]) + col[2]] += val[6];
    HDD[(nxdense * row[1]) + col[3]] += val[7];
  }
  flps += (56 * (EXAGO_FLOPS_SINOP + EXAGO_FLOPS_SINOP) + 462) *
          lineparams->nlineON;
  ierr = PetscLogFlops(flps);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeDenseInequalityConstraintHessian_PBPOLHIOP(
    OPFLOW opflow, const double *x, const double *lambda, double *HDD) {
  int i;
  PBPOLHIOP pbpolhiop = (PBPOLHIOP)opflow->model;
  LINEParams *lineparams = &pbpolhiop->lineparams;
  int row[12], col[12];
  double val[12];
  int nxsparse = opflow->nxsparse;
  int nxdense = opflow->nxdense;
  PetscErrorCode ierr;
  PetscInt flps = 0;

  PetscFunctionBegin;

  // Hessian from line contributions
  for (i = 0; i < lineparams->nlinelim; i++) {
    int j = lineparams->linelimidx[i];
    int gloc;

    double Pf, Qf, Pt, Qt;
    double thetaf = x[lineparams->xidxf[j]], Vmf = x[lineparams->xidxf[j] + 1];
    double thetat = x[lineparams->xidxt[j]], Vmt = x[lineparams->xidxt[j] + 1];
    double thetaft = thetaf - thetat;
    double thetatf = thetat - thetaf;
    double Gff = lineparams->Gff[j], Bff = lineparams->Bff[j];
    double Gft = lineparams->Gft[j], Bft = lineparams->Bft[j];
    double Gtf = lineparams->Gtf[j], Btf = lineparams->Btf[j];
    double Gtt = lineparams->Gtt[j], Btt = lineparams->Btt[j];

    Pf =
        Gff * Vmf * Vmf + Vmf * Vmt * (Gft * cos(thetaft) + Bft * sin(thetaft));
    Qf = -Bff * Vmf * Vmf +
         Vmf * Vmt * (-Bft * cos(thetaft) + Gft * sin(thetaft));

    Pt =
        Gtt * Vmt * Vmt + Vmt * Vmf * (Gtf * cos(thetatf) + Btf * sin(thetatf));
    Qt = -Btt * Vmt * Vmt +
         Vmt * Vmf * (-Btf * cos(thetatf) + Gtf * sin(thetatf));

    double dSf2_dPf, dSf2_dQf, dSt2_dPt, dSt2_dQt;

    dSf2_dPf = 2. * Pf;
    dSf2_dQf = 2. * Qf;
    dSt2_dPt = 2. * Pt;
    dSt2_dQt = 2. * Qt;

    double dPf_dthetaf, dPf_dVmf, dPf_dthetat, dPf_dVmt;
    double dQf_dthetaf, dQf_dVmf, dQf_dthetat, dQf_dVmt;
    double dPt_dthetaf, dPt_dVmf, dPt_dthetat, dPt_dVmt;
    double dQt_dthetaf, dQt_dVmf, dQt_dthetat, dQt_dVmt;

    dPf_dthetaf = Vmf * Vmt * (-Gft * sin(thetaft) + Bft * cos(thetaft));
    dPf_dVmf = 2. * Gff * Vmf + Vmt * (Gft * cos(thetaft) + Bft * sin(thetaft));
    dPf_dthetat = Vmf * Vmt * (Gft * sin(thetaft) - Bft * cos(thetaft));
    dPf_dVmt = Vmf * (Gft * cos(thetaft) + Bft * sin(thetaft));

    dQf_dthetaf = Vmf * Vmt * (Bft * sin(thetaft) + Gft * cos(thetaft));
    dQf_dVmf =
        -2. * Bff * Vmf + Vmt * (-Bft * cos(thetaft) + Gft * sin(thetaft));
    dQf_dthetat = Vmf * Vmt * (-Bft * sin(thetaft) - Gft * cos(thetaft));
    dQf_dVmt = Vmf * (-Bft * cos(thetaft) + Gft * sin(thetaft));

    dPt_dthetat = Vmt * Vmf * (-Gtf * sin(thetatf) + Btf * cos(thetatf));
    dPt_dVmt = 2. * Gtt * Vmt + Vmf * (Gtf * cos(thetatf) + Btf * sin(thetatf));
    dPt_dthetaf = Vmt * Vmf * (Gtf * sin(thetatf) - Btf * cos(thetatf));
    dPt_dVmf = Vmt * (Gtf * cos(thetatf) + Btf * sin(thetatf));

    dQt_dthetat = Vmt * Vmf * (Btf * sin(thetatf) + Gtf * cos(thetatf));
    dQt_dVmt =
        -2. * Btt * Vmt + Vmf * (-Btf * cos(thetatf) + Gtf * sin(thetatf));
    dQt_dthetaf = Vmt * Vmf * (-Btf * sin(thetatf) - Gtf * cos(thetatf));
    dQt_dVmf = Vmt * (-Btf * cos(thetatf) + Gtf * sin(thetatf));

    double d2Pf_dthetaf_dthetaf, d2Pf_dthetaf_dVmf, d2Pf_dthetaf_dthetat,
        d2Pf_dthetaf_dVmt;
    double d2Pf_dVmf_dthetaf, d2Pf_dVmf_dVmf, d2Pf_dVmf_dthetat, d2Pf_dVmf_dVmt;
    double d2Pf_dthetat_dthetaf, d2Pf_dthetat_dVmf, d2Pf_dthetat_dthetat,
        d2Pf_dthetat_dVmt;
    double d2Pf_dVmt_dthetaf, d2Pf_dVmt_dVmf, d2Pf_dVmt_dthetat, d2Pf_dVmt_dVmt;

    /* dPf_dthetaf = Vmf*Vmt*(-Gft*sin(thetaft) + Bft*cos(thetaft)); */
    d2Pf_dthetaf_dthetaf =
        -Vmf * Vmt * (Gft * cos(thetaft) + Bft * sin(thetaft));
    d2Pf_dthetaf_dVmf = Vmt * (-Gft * sin(thetaft) + Bft * cos(thetaft));
    d2Pf_dthetaf_dthetat =
        Vmf * Vmt * (Gft * cos(thetaft) + Bft * sin(thetaft));
    d2Pf_dthetaf_dVmt = Vmf * (-Gft * sin(thetaft) + Bft * cos(thetaft));

    /* dPf_Vmf  = 2*Gff*Vmf + Vmt*(Gft*cos(thetaft) + Bft*sin(thetaft)); */
    d2Pf_dVmf_dthetaf = Vmt * (-Gft * sin(thetaft) + Bft * cos(thetaft));
    d2Pf_dVmf_dVmf = 2 * Gff;
    d2Pf_dVmf_dthetat = Vmt * (Gft * sin(thetaft) - Bft * cos(thetaft));
    d2Pf_dVmf_dVmt = (Gft * cos(thetaft) + Bft * sin(thetaft));

    /* dPf_dthetat = Vmf*Vmt*(Gft*sin(thetaft) - Bft*cos(thetaft)); */
    d2Pf_dthetat_dthetaf =
        Vmf * Vmt * (Gft * cos(thetaft) + Bft * sin(thetaft));
    d2Pf_dthetat_dVmf = Vmt * (Gft * sin(thetaft) - Bft * cos(thetaft));
    d2Pf_dthetat_dthetat =
        Vmf * Vmt * (-Gft * cos(thetaft) - Bft * sin(thetaft));
    d2Pf_dthetat_dVmt = Vmf * (Gft * sin(thetaft) - Bft * cos(thetaft));

    /* dPf_dVmt = Vmf*(Gft*cos(thetaft) + Bft*sin(thetaft)); */
    d2Pf_dVmt_dthetaf = Vmf * (-Gft * sin(thetaft) + Bft * cos(thetaft));
    d2Pf_dVmt_dVmf = (Gft * cos(thetaft) + Bft * sin(thetaft));
    d2Pf_dVmt_dthetat = Vmf * (Gft * sin(thetaft) - Bft * cos(thetaft));
    d2Pf_dVmt_dVmt = 0.0;

    double d2Qf_dthetaf_dthetaf, d2Qf_dthetaf_dVmf, d2Qf_dthetaf_dthetat,
        d2Qf_dthetaf_dVmt;
    double d2Qf_dVmf_dthetaf, d2Qf_dVmf_dVmf, d2Qf_dVmf_dthetat, d2Qf_dVmf_dVmt;
    double d2Qf_dthetat_dthetaf, d2Qf_dthetat_dVmf, d2Qf_dthetat_dthetat,
        d2Qf_dthetat_dVmt;
    double d2Qf_dVmt_dthetaf, d2Qf_dVmt_dVmf, d2Qf_dVmt_dthetat, d2Qf_dVmt_dVmt;

    /* dQf_dthetaf = Vmf*Vmt*(Bft*sin(thetaft) + Gft*cos(thetaft)); */
    d2Qf_dthetaf_dthetaf =
        Vmf * Vmt * (Bft * cos(thetaft) - Gft * sin(thetaft));
    d2Qf_dthetaf_dVmf = Vmt * (Bft * sin(thetaft) + Gft * cos(thetaft));
    d2Qf_dthetaf_dthetat =
        Vmf * Vmt * (-Bft * cos(thetaft) + Gft * sin(thetaft));
    d2Qf_dthetaf_dVmt = Vmf * (Bft * sin(thetaft) + Gft * cos(thetaft));

    /* dQf_dVmf = -2*Bff*Vmf + Vmt*(-Bft*cos(thetaft) + Gft*sin(thetaft)); */
    d2Qf_dVmf_dthetaf = Vmt * (Bft * sin(thetaft) + Gft * cos(thetaft));
    d2Qf_dVmf_dVmf = -2 * Bff;
    d2Qf_dVmf_dthetat = Vmt * (-Bft * sin(thetaft) - Gft * cos(thetaft));
    d2Qf_dVmf_dVmt = (-Bft * cos(thetaft) + Gft * sin(thetaft));

    /* dQf_dthetat = Vmf*Vmt*(-Bft*sin(thetaft) - Gft*cos(thetaft)); */
    d2Qf_dthetat_dthetaf =
        Vmf * Vmt * (-Bft * cos(thetaft) + Gft * sin(thetaft));
    d2Qf_dthetat_dVmf = Vmt * (-Bft * sin(thetaft) - Gft * cos(thetaft));
    d2Qf_dthetat_dthetat =
        Vmf * Vmt * (Bft * cos(thetaft) - Gft * sin(thetaft));
    d2Qf_dthetat_dVmt = Vmf * (-Bft * sin(thetaft) - Gft * cos(thetaft));

    /* dQf_dVmt = Vmf*(-Bft*cos(thetaft) + Gft*sin(thetaft)); */
    d2Qf_dVmt_dthetaf = Vmf * (Bft * sin(thetaft) + Gft * cos(thetaft));
    d2Qf_dVmt_dVmf = (-Bft * cos(thetaft) + Gft * sin(thetaft));
    d2Qf_dVmt_dthetat = Vmf * (-Bft * sin(thetaft) - Gft * cos(thetaft));
    d2Qf_dVmt_dVmt = 0.0;

    double d2Pt_dthetat_dthetat, d2Pt_dthetat_dVmt, d2Pt_dthetat_dthetaf,
        d2Pt_dthetat_dVmf;
    double d2Pt_dVmt_dthetat, d2Pt_dVmt_dVmt, d2Pt_dVmt_dthetaf, d2Pt_dVmt_dVmf;
    double d2Pt_dthetaf_dthetat, d2Pt_dthetaf_dVmt, d2Pt_dthetaf_dthetaf,
        d2Pt_dthetaf_dVmf;
    double d2Pt_dVmf_dthetat, d2Pt_dVmf_dVmt, d2Pt_dVmf_dthetaf, d2Pt_dVmf_dVmf;

    /* dPt_dthetat = Vmf*Vmt*(-Gtf*sin(thetatf) + Btf*cos(thetatf)); */
    d2Pt_dthetat_dthetat =
        Vmf * Vmt * (-Gtf * cos(thetatf) - Btf * sin(thetatf));
    d2Pt_dthetat_dVmt = Vmf * (-Gtf * sin(thetatf) + Btf * cos(thetatf));
    d2Pt_dthetat_dthetaf =
        Vmf * Vmt * (Gtf * cos(thetatf) + Btf * sin(thetatf));
    d2Pt_dthetat_dVmf = Vmt * (-Gtf * sin(thetatf) + Btf * cos(thetatf));

    /* dPt_Vmt  = 2*Gtt*Vmt + Vmf*(Gtf*cos(thetatf) + Btf*sin(thetatf)); */
    d2Pt_dVmt_dthetat = Vmf * (-Gtf * sin(thetatf) + Bft * cos(thetatf));
    d2Pt_dVmt_dVmt = 2 * Gtt;
    d2Pt_dVmt_dthetaf = Vmf * (Gtf * sin(thetatf) - Btf * cos(thetatf));
    d2Pt_dVmt_dVmf = (Gtf * cos(thetatf) + Btf * sin(thetatf));

    /* dPt_dthetaf = Vmf*Vmt*(Gtf*sin(thetatf) - Btf*cos(thetatf)); */
    d2Pt_dthetaf_dthetat =
        Vmf * Vmt * (Gtf * cos(thetatf) + Btf * sin(thetatf));
    d2Pt_dthetaf_dVmt = Vmf * (Gtf * sin(thetatf) - Btf * cos(thetatf));
    d2Pt_dthetaf_dthetaf =
        Vmf * Vmt * (-Gtf * cos(thetatf) - Btf * sin(thetatf));
    d2Pt_dthetaf_dVmf = Vmt * (Gtf * sin(thetatf) - Btf * cos(thetatf));

    /* dPt_dVmf = Vmt*(Gtf*cos(thetatf) + Btf*sin(thetatf)); */
    d2Pt_dVmf_dthetat = Vmt * (-Gtf * sin(thetatf) + Btf * cos(thetatf));
    d2Pt_dVmf_dVmt = (Gtf * cos(thetatf) + Btf * sin(thetatf));
    d2Pt_dVmf_dthetaf = Vmt * (Gtf * sin(thetatf) - Btf * cos(thetatf));
    d2Pt_dVmf_dVmf = 0.0;

    double d2Qt_dthetaf_dthetaf, d2Qt_dthetaf_dVmf, d2Qt_dthetaf_dthetat,
        d2Qt_dthetaf_dVmt;
    double d2Qt_dVmf_dthetaf, d2Qt_dVmf_dVmf, d2Qt_dVmf_dthetat, d2Qt_dVmf_dVmt;
    double d2Qt_dthetat_dthetaf, d2Qt_dthetat_dVmf, d2Qt_dthetat_dthetat,
        d2Qt_dthetat_dVmt;
    double d2Qt_dVmt_dthetaf, d2Qt_dVmt_dVmf, d2Qt_dVmt_dthetat, d2Qt_dVmt_dVmt;

    /* dQt_dthetat = Vmf*Vmt*(Btf*sin(thetatf) + Gtf*cos(thetatf)); */
    d2Qt_dthetat_dthetat =
        Vmf * Vmt * (Btf * cos(thetatf) - Gtf * sin(thetatf));
    d2Qt_dthetat_dVmt = Vmf * (Btf * sin(thetatf) + Gtf * cos(thetatf));
    d2Qt_dthetat_dthetaf =
        Vmf * Vmt * (-Btf * cos(thetatf) + Gtf * sin(thetatf));
    d2Qt_dthetat_dVmf = Vmf * (Btf * sin(thetatf) + Gtf * cos(thetatf));

    /* dQt_dVmt = -2*Btt*Vmt + Vmf*(-Btf*cos(thetatf) + Gtf*sin(thetatf)); */
    d2Qt_dVmt_dthetat = Vmf * (Btf * sin(thetatf) + Gtf * cos(thetatf));
    d2Qt_dVmt_dVmt = -2 * Btt;
    d2Qt_dVmt_dthetaf = Vmf * (-Btf * sin(thetatf) + Gtf * cos(thetatf));
    d2Qt_dVmt_dVmf = (-Btf * cos(thetatf) + Gtf * sin(thetatf));

    /* dQt_dthetaf = Vmf*Vmt*(-Btf*sin(thetatf) - Gtf*cos(thetatf)); */
    d2Qt_dthetaf_dthetat =
        Vmf * Vmt * (-Btf * cos(thetatf) + Gtf * sin(thetatf));
    d2Qt_dthetaf_dVmt = Vmf * (-Btf * sin(thetatf) - Gtf * cos(thetatf));
    d2Qt_dthetaf_dthetaf =
        Vmf * Vmt * (Btf * cos(thetatf) - Gtf * sin(thetatf));
    d2Qt_dthetaf_dVmf = Vmt * (-Btf * sin(thetatf) - Gtf * cos(thetatf));

    /* dQt_dVmf = Vmt*(-Btf*cos(thetatf) + Gtf*sin(thetatf)); */
    d2Qt_dVmf_dthetat = Vmt * (Btf * sin(thetatf) + Gtf * cos(thetatf));
    d2Qt_dVmf_dVmt = (-Btf * cos(thetatf) + Gtf * sin(thetatf));
    d2Qt_dVmf_dthetaf = Vmt * (-Btf * sin(thetatf) - Gtf * cos(thetatf));
    d2Qt_dVmf_dVmf = 0.0;

    double d2Sf2_dthetaf_dthetaf = 0.0, d2Sf2_dthetaf_dVmf = 0.0,
           d2Sf2_dthetaf_dthetat = 0.0, d2Sf2_dthetaf_dVmt = 0.0;
    double d2St2_dthetaf_dthetaf = 0.0, d2St2_dthetaf_dVmf = 0.0,
           d2St2_dthetaf_dthetat = 0.0, d2St2_dthetaf_dVmt = 0.0;

    d2Sf2_dthetaf_dthetaf =
        2 * dPf_dthetaf * dPf_dthetaf + dSf2_dPf * d2Pf_dthetaf_dthetaf +
        2 * dQf_dthetaf * dQf_dthetaf + dSf2_dQf * d2Qf_dthetaf_dthetaf;
    d2Sf2_dthetaf_dVmf =
        2 * dPf_dVmf * dPf_dthetaf + dSf2_dPf * d2Pf_dthetaf_dVmf +
        2 * dQf_dVmf * dQf_dthetaf + dSf2_dQf * d2Qf_dthetaf_dVmf;
    d2Sf2_dthetaf_dthetat =
        2 * dPf_dthetat * dPf_dthetaf + dSf2_dPf * d2Pf_dthetaf_dthetat +
        2 * dQf_dthetat * dQf_dthetaf + dSf2_dQf * d2Qf_dthetaf_dthetat;
    d2Sf2_dthetaf_dVmt =
        2 * dPf_dVmt * dPf_dthetaf + dSf2_dPf * d2Pf_dthetaf_dVmt +
        2 * dQf_dVmt * dQf_dthetaf + dSf2_dQf * d2Qf_dthetaf_dVmt;

    d2St2_dthetaf_dthetaf =
        2 * dPt_dthetaf * dPt_dthetaf + dSt2_dPt * d2Pt_dthetaf_dthetaf +
        2 * dQt_dthetaf * dQt_dthetaf + dSt2_dQt * d2Qt_dthetaf_dthetaf;
    d2St2_dthetaf_dVmf =
        2 * dPt_dVmf * dPt_dthetaf + dSt2_dPt * d2Pt_dthetaf_dVmf +
        2 * dQt_dVmf * dQt_dthetaf + dSt2_dQt * d2Qt_dthetaf_dVmf;
    d2St2_dthetaf_dthetat =
        2 * dPt_dthetat * dPt_dthetaf + dSt2_dPt * d2Pt_dthetaf_dthetat +
        2 * dQt_dthetat * dQt_dthetaf + dSt2_dQt * d2Qt_dthetaf_dthetat;
    d2St2_dthetaf_dVmt =
        2 * dPt_dVmt * dPt_dthetaf + dSt2_dPt * d2Pt_dthetaf_dVmt +
        2 * dQt_dVmt * dQt_dthetaf + dSt2_dQt * d2Qt_dthetaf_dVmt;

    val[0] = val[1] = val[2] = val[3] = 0.0;

    row[0] = lineparams->xidxf[j] - nxsparse;
    col[0] = lineparams->xidxf[j] - nxsparse;
    col[1] = lineparams->xidxf[j] + 1 - nxsparse;
    col[2] = lineparams->xidxt[j] - nxsparse;
    col[3] = lineparams->xidxt[j] + 1 - nxsparse;

    gloc = lineparams->gineqidx[i];

    val[0] = lambda[gloc] * d2Sf2_dthetaf_dthetaf +
             lambda[gloc + 1] * d2St2_dthetaf_dthetaf;
    val[1] = lambda[gloc] * d2Sf2_dthetaf_dVmf +
             lambda[gloc + 1] * d2St2_dthetaf_dVmf;
    val[2] = lambda[gloc] * d2Sf2_dthetaf_dthetat +
             lambda[gloc + 1] * d2St2_dthetaf_dthetat;
    val[3] = lambda[gloc] * d2Sf2_dthetaf_dVmt +
             lambda[gloc + 1] * d2St2_dthetaf_dVmt;

    HDD[(nxdense * row[0]) + col[0]] += val[0];
    HDD[(nxdense * row[0]) + col[1]] += val[1];
    HDD[(nxdense * row[0]) + col[2]] += val[2];
    HDD[(nxdense * row[0]) + col[3]] += val[3];

    double d2Sf2_dVmf_dthetaf, d2Sf2_dVmf_dVmf, d2Sf2_dVmf_dthetat,
        d2Sf2_dVmf_dVmt;
    double d2St2_dVmf_dthetaf, d2St2_dVmf_dVmf, d2St2_dVmf_dthetat,
        d2St2_dVmf_dVmt;

    d2Sf2_dVmf_dthetaf =
        2 * dPf_dthetaf * dPf_dVmf + dSf2_dPf * d2Pf_dVmf_dthetaf +
        2 * dQf_dthetaf * dQf_dVmf + dSf2_dQf * d2Qf_dVmf_dthetaf;
    d2Sf2_dVmf_dVmf = 2 * dPf_dVmf * dPf_dVmf + dSf2_dPf * d2Pf_dVmf_dVmf +
                      2 * dQf_dVmf * dQf_dVmf + dSf2_dQf * d2Qf_dVmf_dVmf;
    d2Sf2_dVmf_dthetat =
        2 * dPf_dthetat * dPf_dVmf + dSf2_dPf * d2Pf_dVmf_dthetat +
        2 * dQf_dthetat * dQf_dVmf + dSf2_dQf * d2Qf_dVmf_dthetat;
    d2Sf2_dVmf_dVmt = 2 * dPf_dVmt * dPf_dVmf + dSf2_dPf * d2Pf_dVmf_dVmt +
                      2 * dQf_dVmt * dQf_dVmf + dSf2_dQf * d2Qf_dVmf_dVmt;

    d2St2_dVmf_dthetaf =
        2 * dPt_dthetaf * dPt_dVmf + dSt2_dPt * d2Pt_dVmf_dthetaf +
        2 * dQt_dthetaf * dQt_dVmf + dSt2_dQt * d2Qt_dVmf_dthetaf;
    d2St2_dVmf_dVmf = 2 * dPt_dVmf * dPt_dVmf + dSt2_dPt * d2Pt_dVmf_dVmf +
                      2 * dQt_dVmf * dQt_dVmf + dSt2_dQt * d2Qt_dVmf_dVmf;
    d2St2_dVmf_dthetat =
        2 * dPt_dthetat * dPt_dVmf + dSt2_dPt * d2Pt_dVmf_dthetat +
        2 * dQt_dthetat * dQt_dVmf + dSt2_dQt * d2Qt_dVmf_dthetat;
    d2St2_dVmf_dVmt = 2 * dPt_dVmt * dPt_dVmf + dSt2_dPt * d2Pt_dVmf_dVmt +
                      2 * dQt_dVmt * dQt_dVmf + dSt2_dQt * d2Qt_dVmf_dVmt;

    val[0] = val[1] = val[2] = val[3] = 0.0;
    col[0] = lineparams->xidxf[j] - nxsparse;
    col[1] = lineparams->xidxf[j] + 1 - nxsparse;
    col[2] = lineparams->xidxt[j] - nxsparse;
    col[3] = lineparams->xidxt[j] + 1 - nxsparse;

    row[0] = lineparams->xidxf[j] + 1 - nxsparse;

    val[0] = lambda[gloc] * d2Sf2_dVmf_dthetaf +
             lambda[gloc + 1] * d2St2_dVmf_dthetaf;
    val[1] =
        lambda[gloc] * d2Sf2_dVmf_dVmf + lambda[gloc + 1] * d2St2_dVmf_dVmf;
    val[2] = lambda[gloc] * d2Sf2_dVmf_dthetat +
             lambda[gloc + 1] * d2St2_dVmf_dthetat;
    val[3] =
        lambda[gloc] * d2Sf2_dVmf_dVmt + lambda[gloc + 1] * d2St2_dVmf_dVmt;

    HDD[(nxdense * row[0]) + col[0]] += val[0];
    HDD[(nxdense * row[0]) + col[1]] += val[1];
    HDD[(nxdense * row[0]) + col[2]] += val[2];
    HDD[(nxdense * row[0]) + col[3]] += val[3];

    double d2Sf2_dthetat_dthetaf, d2Sf2_dthetat_dVmf, d2Sf2_dthetat_dthetat,
        d2Sf2_dthetat_dVmt;
    double d2St2_dthetat_dthetaf, d2St2_dthetat_dVmf, d2St2_dthetat_dthetat,
        d2St2_dthetat_dVmt;

    d2Sf2_dthetat_dthetaf =
        2 * dPf_dthetaf * dPf_dthetat + dSf2_dPf * d2Pf_dthetat_dthetaf +
        2 * dQf_dthetat * dQf_dthetaf + dSf2_dQf * d2Qf_dthetat_dthetaf;
    d2Sf2_dthetat_dVmf =
        2 * dPf_dVmf * dPf_dthetat + dSf2_dPf * d2Pf_dthetat_dVmf +
        2 * dQf_dthetat * dQf_dVmf + dSf2_dQf * d2Qf_dthetat_dVmf;
    d2Sf2_dthetat_dthetat =
        2 * dPf_dthetat * dPf_dthetat + dSf2_dPf * d2Pf_dthetat_dthetat +
        2 * dQf_dthetat * dQf_dthetat + dSf2_dQf * d2Qf_dthetat_dthetat;
    d2Sf2_dthetat_dVmt =
        2 * dPf_dVmt * dPf_dthetat + dSf2_dPf * d2Pf_dthetat_dVmt +
        2 * dQf_dthetat * dQf_dVmt + dSf2_dQf * d2Qf_dthetat_dVmt;

    d2St2_dthetat_dthetaf =
        2 * dPt_dthetaf * dPt_dthetat + dSt2_dPt * d2Pt_dthetat_dthetaf +
        2 * dQt_dthetaf * dQt_dthetat + dSt2_dQt * d2Qt_dthetat_dthetaf;
    d2St2_dthetat_dVmf =
        2 * dPt_dVmf * dPt_dthetat + dSt2_dPt * d2Pt_dthetat_dVmf +
        2 * dQt_dVmf * dQt_dthetat + dSt2_dQt * d2Qt_dthetat_dVmf;
    d2St2_dthetat_dthetat =
        2 * dPt_dthetat * dPt_dthetat + dSt2_dPt * d2Pt_dthetat_dthetat +
        2 * dQt_dthetat * dQt_dthetat + dSt2_dQt * d2Qt_dthetat_dthetat;
    d2St2_dthetat_dVmt =
        2 * dPt_dVmt * dPt_dthetat + dSt2_dPt * d2Pt_dthetat_dVmt +
        2 * dQt_dVmt * dQt_dthetat + dSt2_dQt * d2Qt_dthetat_dVmt;

    val[0] = val[1] = val[2] = val[3] = 0.0;

    col[0] = lineparams->xidxf[j] - nxsparse;
    col[1] = lineparams->xidxf[j] + 1 - nxsparse;
    col[2] = lineparams->xidxt[j] - nxsparse;
    col[3] = lineparams->xidxt[j] + 1 - nxsparse;

    row[0] = lineparams->xidxt[j] - nxsparse;

    val[0] = lambda[gloc] * d2Sf2_dthetat_dthetaf +
             lambda[gloc + 1] * d2St2_dthetat_dthetaf;
    val[1] = lambda[gloc] * d2Sf2_dthetat_dVmf +
             lambda[gloc + 1] * d2St2_dthetat_dVmf;
    val[2] = lambda[gloc] * d2Sf2_dthetat_dthetat +
             lambda[gloc + 1] * d2St2_dthetat_dthetat;
    val[3] = lambda[gloc] * d2Sf2_dthetat_dVmt +
             lambda[gloc + 1] * d2St2_dthetat_dVmt;

    HDD[(nxdense * row[0]) + col[0]] += val[0];
    HDD[(nxdense * row[0]) + col[1]] += val[1];
    HDD[(nxdense * row[0]) + col[2]] += val[2];
    HDD[(nxdense * row[0]) + col[3]] += val[3];

    double d2Sf2_dVmt_dthetaf, d2Sf2_dVmt_dVmf, d2Sf2_dVmt_dthetat,
        d2Sf2_dVmt_dVmt;
    double d2St2_dVmt_dthetaf, d2St2_dVmt_dVmf, d2St2_dVmt_dthetat,
        d2St2_dVmt_dVmt;

    d2Sf2_dVmt_dthetaf =
        2 * dPf_dthetaf * dPf_dVmt + dSf2_dPf * d2Pf_dVmt_dthetaf +
        2 * dQf_dthetaf * dQf_dVmt + dSf2_dQf * d2Qf_dVmt_dthetaf;
    d2Sf2_dVmt_dVmf = 2 * dPf_dVmf * dPf_dVmt + dSf2_dPf * d2Pf_dVmt_dVmf +
                      2 * dQf_dVmf * dQf_dVmt + dSf2_dQf * d2Qf_dVmt_dVmf;
    d2Sf2_dVmt_dthetat =
        2 * dPf_dthetat * dPf_dVmt + dSf2_dPf * d2Pf_dVmt_dthetat +
        2 * dQf_dthetat * dQf_dVmt + dSf2_dQf * d2Qf_dVmt_dthetat;
    d2Sf2_dVmt_dVmt = 2 * dPf_dVmt * dPf_dVmt + dSf2_dPf * d2Pf_dVmt_dVmt +
                      2 * dQf_dVmt * dQf_dVmt + dSf2_dQf * d2Qf_dVmt_dVmt;

    d2St2_dVmt_dthetaf =
        2 * dPt_dthetaf * dPt_dVmt + dSt2_dPt * d2Pt_dVmt_dthetaf +
        2 * dQt_dthetaf * dQt_dVmt + dSt2_dQt * d2Qt_dVmt_dthetaf;
    d2St2_dVmt_dVmf = 2 * dPt_dVmf * dPt_dVmt + dSt2_dPt * d2Pt_dVmt_dVmf +
                      2 * dQt_dVmf * dQt_dVmt + dSt2_dQt * d2Qt_dVmt_dVmf;
    d2St2_dVmt_dthetat =
        2 * dPt_dthetat * dPt_dVmt + dSt2_dPt * d2Pt_dVmt_dthetat +
        2 * dQt_dthetat * dQt_dVmt + dSt2_dQt * d2Qt_dVmt_dthetat;
    d2St2_dVmt_dVmt = 2 * dPt_dVmt * dPt_dVmt + dSt2_dPt * d2Pt_dVmt_dVmt +
                      2 * dQt_dVmt * dQt_dVmt + dSt2_dQt * d2Qt_dVmt_dVmt;

    val[0] = val[1] = val[2] = val[3] = 0.0;

    row[0] = lineparams->xidxt[j] + 1 - nxsparse;
    col[0] = lineparams->xidxf[j] - nxsparse;
    col[1] = lineparams->xidxf[j] + 1 - nxsparse;
    col[2] = lineparams->xidxt[j] - nxsparse;
    col[3] = lineparams->xidxt[j] + 1 - nxsparse;

    val[0] = lambda[gloc] * d2Sf2_dVmt_dthetaf +
             lambda[gloc + 1] * d2St2_dVmt_dthetaf;
    val[1] =
        lambda[gloc] * d2Sf2_dVmt_dVmf + lambda[gloc + 1] * d2St2_dVmt_dVmf;
    val[2] = lambda[gloc] * d2Sf2_dVmt_dthetat +
             lambda[gloc + 1] * d2St2_dVmt_dthetat;
    val[3] =
        lambda[gloc] * d2Sf2_dVmt_dVmt + lambda[gloc + 1] * d2St2_dVmt_dVmt;

    HDD[(nxdense * row[0]) + col[0]] += val[0];
    HDD[(nxdense * row[0]) + col[1]] += val[1];
    HDD[(nxdense * row[0]) + col[2]] += val[2];
    HDD[(nxdense * row[0]) + col[3]] += val[3];
  }
  flps += (972 + (92 * EXAGO_FLOPS_COSOP) + (92 * EXAGO_FLOPS_SINOP)) *
          lineparams->nlinelim;
  ierr = PetscLogFlops(flps);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeDenseHessian_PBPOLHIOP(OPFLOW opflow,
                                                   const double *x,
                                                   const double *lambda,
                                                   double *HDD) {
  PetscErrorCode ierr;
  int i, j;
  int nxdense = opflow->nxdense;

  if (!HDD)
    PetscFunctionReturn(0);

  for (i = 0; i < nxdense * nxdense; i++)
    HDD[i] = 0.0;

  /* Equality constraint Hessian */
  ierr = OPFLOWComputeDenseEqualityConstraintHessian_PBPOLHIOP(opflow, x,
                                                               lambda, HDD);
  CHKERRQ(ierr);

  if (opflow->nconineq) {
    ierr = OPFLOWComputeDenseInequalityConstraintHessian_PBPOLHIOP(
        opflow, x, lambda + opflow->nconeq, HDD);
    CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}
#endif
