#include "dcopf.h"
#include "exago_config.h"
#include <private/opflowimpl.h>

PetscErrorCode OPFLOWModelDestroy_DCOPF(OPFLOW opflow) {
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = PetscFree(opflow->model);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSetVariableBounds_DCOPF(OPFLOW opflow, Vec Xl, Vec Xu) {
  PetscErrorCode ierr;
  PS ps = opflow->ps;
  PetscScalar *xl, *xu;
  PetscInt i;
  PSBUS bus;
  PSGEN gen;
  PetscInt loc;

  PetscFunctionBegin;

  /* Get array pointers */
  ierr = VecGetArray(Xl, &xl);
  CHKERRQ(ierr);
  ierr = VecGetArray(Xu, &xu);
  CHKERRQ(ierr);

  for (i = 0; i < ps->nbus; i++) {
    PetscInt k;

    bus = &ps->bus[i];

    /* Bounds on bus variables */
    loc = bus->startxVloc;

    xl[loc] = PETSC_NINFINITY;
    xu[loc] = PETSC_INFINITY;

    if (bus->ide == REF_BUS || bus->ide == ISOLATED_BUS)
      xl[loc] = xu[loc] = bus->va;

    if (opflow->include_powerimbalance_variables) {
      loc = bus->startxpimbloc;

      xl[loc] = xl[loc + 1] = 0.0;
      xu[loc] = xu[loc + 1] = PETSC_INFINITY;
    }

    /* Bounds on generator variables */

    for (k = 0; k < bus->ngen; k++) {
      ierr = PSBUSGetGen(bus, k, &gen);
      CHKERRQ(ierr);
      if (!gen->status)
        continue;

      loc = gen->startxpowloc;

      xl[loc] = gen->pb; /* PGmin */
      xu[loc] = gen->pt; /* PGmax */
      /* pb, pt are converted to p.u. in ps.c */

      if (opflow->has_gensetpoint && !gen->isrenewable) {
        loc = gen->startxpdevloc;

        xl[loc] = gen->pb - gen->pt;
        xu[loc] = gen->pt - gen->pb;

        loc = gen->startxpsetloc;

        xl[loc] = gen->pb;
        xu[loc] = gen->pt;
      }
    }

    if (opflow->include_loadloss_variables) {
      for (k = 0; k < bus->nload; k++) {
        PSLOAD load;
        ierr = PSBUSGetLoad(bus, k, &load);
        CHKERRQ(ierr);
        loc = load->startxloadlossloc;
        if (!load->status)
          xl[loc] = xu[loc] = 0.0;
        else {
          xl[loc] = PetscMin(0.0, load->pl);
          xu[loc] = PetscMax(0.0, load->pl);
        }
      }
    }

    if (opflow->use_agc) {
      if (bus->ide == REF_BUS) {
        loc = ps->startxloc;
        xl[loc] = PETSC_NINFINITY;
        xu[loc] = PETSC_INFINITY;
      }
    }
  }

  ierr = VecRestoreArray(Xl, &xl);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(Xu, &xu);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSetConstraintBounds_DCOPF(OPFLOW opflow, Vec Gl, Vec Gu) {
  PetscErrorCode ierr;
  PS ps = opflow->ps;
  PetscScalar *gl, *gu;
  PetscInt i, k, ngen;
  PSBUS bus;
  PSLINE line;
  PSGEN gen;
  PetscInt gloc;

  PetscFunctionBegin;

  /* Get array pointers */
  ierr = VecGetArray(Gl, &gl);
  CHKERRQ(ierr);
  ierr = VecGetArray(Gu, &gu);
  CHKERRQ(ierr);

  /** EQUALITY CONSTRAINTS **/
  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    gloc = bus->starteqloc;
    /* Equality constraint bounds for real power mismatch at the
     * bus */
    gl[gloc] = 0.0;
    gu[gloc] = 0.0;

    if (opflow->has_gensetpoint) {
      ierr = PSBUSGetNGen(bus, &ngen);
      CHKERRQ(ierr);
      for (k = 0; k < ngen; k++) {
        ierr = PSBUSGetGen(bus, k, &gen);
        CHKERRQ(ierr);
        if (!gen->status || gen->isrenewable)
          continue;
        gloc = gen->starteqloc;
        gl[gloc] = gu[gloc] = 0.0;
        gl[gloc + 1] = gu[gloc + 1] = 0.0;
      }
    }
  }

  /** INEQUALITY CONSTRAINTS **/
  gl += opflow->nconeq;
  gu += opflow->nconeq;

  if (opflow->has_gensetpoint) {
    for (i = 0; i < ps->nbus; i++) {
      bus = &ps->bus[i];
      for (k = 0; k < bus->ngen; k++) {
        ierr = PSBUSGetGen(bus, k, &gen);
        CHKERRQ(ierr);
        if (!gen->status || gen->isrenewable)
          continue;
        gloc = gen->startineqloc;
        if (opflow->use_agc) {
          gl[gloc] = gl[gloc + 1] = 0.0;
          gu[gloc] = gu[gloc + 1] = PETSC_INFINITY;
        }
      }
    }
  }
  /* Line flow constraint bounds */
  if (!opflow->ignore_lineflow_constraints) {
    for (i = 0; i < ps->nline; i++) {
      line = &ps->line[i];
      if (!line->status || line->rateA > 1e5)
        continue;

      gloc = line->startineqloc;
      /* Line flow inequality constraints */
      gl[gloc] = gl[gloc + 1] = -(line->rateA / ps->MVAbase);
      gu[gloc] = gu[gloc + 1] = (line->rateA / ps->MVAbase);
    }
  }

  ierr = VecRestoreArray(Gl, &gl);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(Gu, &gu);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSetVariableandConstraintBounds_DCOPF(OPFLOW opflow, Vec Xl,
                                                          Vec Xu, Vec Gl,
                                                          Vec Gu) {
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = OPFLOWSetVariableBounds_DCOPF(opflow, Xl, Xu);
  CHKERRQ(ierr);
  ierr = OPFLOWSetConstraintBounds_DCOPF(opflow, Gl, Gu);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSetInitialGuess_DCOPF(OPFLOW opflow, Vec X, Vec Lambda) {
  PetscErrorCode ierr;
  PS ps = opflow->ps;
  const PetscScalar *xl, *xu;
  PetscScalar *x;
  PetscInt i;
  PSBUS bus;
  PetscInt loc;

  PetscFunctionBegin;

  /* Get array pointers */
  ierr = VecGetArray(X, &x);
  CHKERRQ(ierr);
  ierr = VecGetArrayRead(opflow->Xl, &xl);
  CHKERRQ(ierr);
  ierr = VecGetArrayRead(opflow->Xu, &xu);
  CHKERRQ(ierr);

  for (i = 0; i < ps->nbus; i++) {
    PetscInt k;

    bus = &ps->bus[i];

    loc = bus->startxVloc;

    if (bus->ide == ISOLATED_BUS) {
      x[loc] = bus->va;
    } else {
      if (opflow->initializationtype == OPFLOWINIT_MIDPOINT) {
        /* Initial guess for voltage angles and bounds on voltage magnitudes */
        x[loc] = (xl[loc] + xu[loc]) / 2.0;
      } else if (opflow->initializationtype == OPFLOWINIT_FROMFILE ||
                 opflow->initializationtype == OPFLOWINIT_ACPF) {
        x[loc] = bus->va;
      } else if (opflow->initializationtype == OPFLOWINIT_FLATSTART) {
        x[loc] = 0.0;
      }
    }

    if (opflow->include_powerimbalance_variables) {
      loc = bus->startxpimbloc;
      x[loc] = x[loc + 1] = 0.0;
    }

    for (k = 0; k < bus->ngen; k++) {
      PSGEN gen;

      ierr = PSBUSGetGen(bus, k, &gen);
      CHKERRQ(ierr);
      if (!gen->status)
        continue;

      loc = gen->startxpowloc;

      if (opflow->initializationtype == OPFLOWINIT_MIDPOINT ||
          opflow->initializationtype == OPFLOWINIT_FLATSTART) {
        x[loc] = 0.5 * (xl[loc] + xu[loc]);
      } else if (opflow->initializationtype == OPFLOWINIT_FROMFILE ||
                 opflow->initializationtype == OPFLOWINIT_ACPF) {
        x[loc] = PetscMax(gen->pb, PetscMin(gen->pg, gen->pt));
      }

      if (opflow->has_gensetpoint && !gen->isrenewable) {
        loc = gen->startxpdevloc;
        x[loc] = 0.0;
        loc = gen->startxpsetloc;
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
        loc = load->startxloadlossloc;
        /* Initial value for real power load loss */
        x[loc] = 0.0;
      }
    }
    if (opflow->use_agc) {
      if (bus->ide == REF_BUS) {
        loc = ps->startxloc;
        x[loc] = 0.0;
      }
    }
  }

  ierr = VecRestoreArray(X, &x);
  CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(opflow->Xl, &xl);
  CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(opflow->Xu, &xu);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeEqualityConstraints_DCOPF(OPFLOW opflow, Vec X,
                                                      Vec Ge) {
  PetscErrorCode ierr;
  PetscInt i, k, nconnlines;
  PetscInt gloc, row[2];
  PetscInt xloc, xlocf, xloct;
  PetscScalar val[2];
  PetscScalar Pg, Pd;
  PetscScalar Bdc, Pshift;
  PetscScalar thetaf, thetat, thetaft, thetatf;
  PetscScalar Pf, Pt;
  PetscScalar theta;
  PS ps = opflow->ps;
  PSLOAD load;
  PSLINE line;
  PSBUS bus, busf, bust;
  PSGEN gen;
  const PSBUS *connbuses;
  const PSLINE *connlines;
  const PetscScalar *x;
  double flps = 0.0;
  PetscScalar delP = 0.0;

  PetscFunctionBegin;
  ierr = VecSet(Ge, 0.0);
  CHKERRQ(ierr);

  ierr = VecGetArrayRead(X, &x);
  CHKERRQ(ierr);

  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    gloc = bus->starteqloc;

    row[0] = gloc;

    xloc = bus->startxVloc;

    theta = x[xloc];

    if (bus->ide == ISOLATED_BUS) {
      row[0] = gloc;
      val[0] = theta - bus->va;
      ierr = VecSetValues(Ge, 1, row, val, ADD_VALUES);
      CHKERRQ(ierr);
      continue;
    }

    /* Power imbalance addition */
    if (opflow->include_powerimbalance_variables) {
      PetscScalar Pimbplus, Pimbminus, Pimb;
      xloc = bus->startxpimbloc;
      Pimbplus = x[xloc];
      Pimbminus = x[xloc + 1];

      Pimb = Pimbplus - Pimbminus;
      val[0] = Pimb;

      ierr = VecSetValues(Ge, 1, row, val, ADD_VALUES);
      CHKERRQ(ierr);

      flps += 1.0;
    }

    for (k = 0; k < bus->ngen; k++) {
      ierr = PSBUSGetGen(bus, k, &gen);
      CHKERRQ(ierr);
      if (!gen->status)
        continue;

      xloc = gen->startxpowloc;

      Pg = x[xloc];

      val[0] = -Pg;

      ierr = VecSetValues(Ge, 1, row, val, ADD_VALUES);
      CHKERRQ(ierr);

      flps += 2.0;
    }

    for (k = 0; k < bus->nload; k++) {
      ierr = PSBUSGetLoad(bus, k, &load);
      CHKERRQ(ierr);
      xloc = load->startxloadlossloc;
      if (opflow->include_loadloss_variables) {
        Pd = load->pl - x[xloc];
      } else {
        Pd = load->pl;
      }

      val[0] = Pd;

      ierr = VecSetValues(Ge, 1, row, val, ADD_VALUES);
      CHKERRQ(ierr);
      flps += 1.0;
    }

    ierr = PSBUSGetSupportingLines(bus, &nconnlines, &connlines);
    CHKERRQ(ierr);
    for (k = 0; k < nconnlines; k++) {
      line = connlines[k];
      if (!line->status)
        continue;

      Bdc = line->bdc;
      Pshift = line->pshift;

      ierr = PSLINEGetConnectedBuses(line, &connbuses);
      CHKERRQ(ierr);
      busf = connbuses[0];
      bust = connbuses[1];

      xlocf = busf->startxVloc;
      xloct = bust->startxVloc;

      thetaf = x[xlocf];
      thetat = x[xloct];
      thetaft = thetaf - thetat;
      thetatf = thetat - thetaf;

      if (bus == busf) {
        Pf = Bdc * thetaft + Pshift;

        val[0] = Pf;

        ierr = VecSetValues(Ge, 1, row, val, ADD_VALUES);
        CHKERRQ(ierr);

        flps += 1.0;
      } else {
        Pt = Bdc * thetatf - Pshift;

        val[0] = Pt;

        ierr = VecSetValues(Ge, 1, row, val, ADD_VALUES);
        CHKERRQ(ierr);
        flps += 1.0;
      }
    }

    if (opflow->has_gensetpoint) {
      for (k = 0; k < bus->ngen; k++) {
        PetscScalar delPg; /* Ramp */
        PetscScalar Pgset;
        ierr = PSBUSGetGen(bus, k, &gen);
        CHKERRQ(ierr);
        if (!gen->status || gen->isrenewable)
          continue;

        Pg = x[gen->startxpowloc];

        delPg = x[gen->startxpdevloc];
        Pgset = x[gen->startxpsetloc];
        gloc = gen->starteqloc;

        row[0] = gloc;
        val[0] = Pgset + delPg - Pg;

        ierr = VecSetValues(Ge, 1, row, val, ADD_VALUES);
        CHKERRQ(ierr);

        row[0] = gloc + 1;
        val[0] = Pgset - gen->pgs;

        ierr = VecSetValues(Ge, 1, row, val, ADD_VALUES);
        CHKERRQ(ierr);

        flps += 3.0;
      }
    }
  }
  ierr = VecAssemblyBegin(Ge);
  CHKERRQ(ierr);
  ierr = VecAssemblyEnd(Ge);
  CHKERRQ(ierr);

  ierr = VecRestoreArrayRead(X, &x);
  CHKERRQ(ierr);

  PetscLogFlops(flps);
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeEqualityConstraintJacobian_DCOPF(OPFLOW opflow,
                                                             Vec X, Mat Je) {
  PetscErrorCode ierr;
  PetscInt i, k, row[2], col[4], gidx, flps = 0;
  PetscInt nconnlines, locglob, loc, locglobf, locglobt, locf, loct;
  PetscScalar val[8], Bdc, Pshift;
  PetscScalar thetaf, thetat, thetaft, thetatf;
  PS ps = opflow->ps;
  PSBUS bus;
  PSLINE line;
  PSBUS busf, bust;
  PSGEN gen;
  PSLOAD load;
  DM networkdm = ps->networkdm;
  const PSLINE *connlines;
  const PSBUS *connbuses;
  const PetscScalar *xarr;

  PetscFunctionBegin;
  ierr = MatZeroEntries(Je);
  CHKERRQ(ierr);

  ierr = VecGetArrayRead(X, &xarr);
  CHKERRQ(ierr);

  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    gidx = bus->starteqloc;

    row[0] = gidx;

    loc = bus->startxVloc;
    locglob = bus->startxVlocglob;

    col[0] = locglob;
    col[1] = locglob + 1;
    /* Isolated and reference bus */
    if (bus->ide == ISOLATED_BUS) {
      val[0] = 1.0;
      ierr = MatSetValues(Je, 1, row, 1, col, val, ADD_VALUES);
      CHKERRQ(ierr);
      continue;
    }

    /* Power imbalance Jacobian terms */
    if (opflow->include_powerimbalance_variables) {
      locglob = bus->startxpimblocglob;

      val[0] = 1;
      val[1] = -1;
      col[0] = locglob;
      col[1] = locglob + 1;
      ierr = MatSetValues(Je, 1, row, 2, col, val, ADD_VALUES);
      CHKERRQ(ierr);
    }

    for (k = 0; k < bus->ngen; k++) {
      ierr = PSBUSGetGen(bus, k, &gen);
      CHKERRQ(ierr);
      if (!gen->status)
        continue;

      locglob = gen->startxpowlocglob;

      val[0] = -1;
      col[0] = locglob;
      ierr = MatSetValues(Je, 1, row, 1, col, val, ADD_VALUES);
      CHKERRQ(ierr);
    }

    if (opflow->include_loadloss_variables) {
      for (k = 0; k < bus->nload; k++) {
        ierr = PSBUSGetLoad(bus, k, &load);
        CHKERRQ(ierr);

        locglob = load->startxloadlosslocglob;

        val[0] = -1;
        col[0] = locglob;
        ierr = MatSetValues(Je, 1, row, 1, col, val, ADD_VALUES);
        CHKERRQ(ierr);
      }
    }

    /* Partial derivatives of network equations */
    /* Get the lines supporting the bus */
    ierr = PSBUSGetSupportingLines(bus, &nconnlines, &connlines);
    CHKERRQ(ierr);

    for (k = 0; k < nconnlines; k++) {
      line = connlines[k];
      if (!line->status)
        continue;
      Bdc = line->bdc;
      Pshift = line->pshift;

      /* Get the connected buses to this line */
      ierr = PSLINEGetConnectedBuses(line, &connbuses);
      CHKERRQ(ierr);
      busf = connbuses[0];
      bust = connbuses[1];

      locf = busf->startxVloc;
      loct = bust->startxVloc;

      locglobf = busf->startxVlocglob;
      locglobt = bust->startxVlocglob;

      thetaf = xarr[locf];
      thetat = xarr[loct];
      thetaft = thetaf - thetat;
      thetatf = thetat - thetaf;

      if (bus == busf) {
        col[0] = locglobf;
        col[1] = locglobt;
        /* dPf_dthetaf */
        val[0] = Bdc;
        /*dPf_dthetat */
        val[1] = -Bdc;

        ierr = MatSetValues(Je, 1, row, 2, col, val, ADD_VALUES);
        CHKERRQ(ierr);
      } else {
        col[0] = locglobt;
        col[1] = locglobf;

        /* dPt_dthetat */
        val[0] = Bdc;
        /* dPt_dthetaf */
        val[1] = -Bdc;

        ierr = MatSetValues(Je, 1, row, 2, col, val, ADD_VALUES);
        CHKERRQ(ierr);
      }
    }
    flps += nconnlines * 2;

    if (opflow->has_gensetpoint) {
      ierr = PSBUSGetVariableGlobalLocation(bus, &locglob);
      CHKERRQ(ierr);

      for (k = 0; k < bus->ngen; k++) {
        ierr = PSBUSGetGen(bus, k, &gen);
        CHKERRQ(ierr);
        if (!gen->status || gen->isrenewable)
          continue;

        locglob = gen->startxpowlocglob;

        gidx = gen->starteqloc;
        row[0] = gidx;

        col[0] = locglob;
        val[0] = -1.0;

        locglob = gen->startxpdevlocglob;

        col[1] = locglob;
        val[1] = 1.0;

        locglob = gen->startxpsetlocglob;
        col[2] = locglob;
        val[2] = 1.0;

        ierr = MatSetValues(Je, 1, row, 3, col, val, ADD_VALUES);
        CHKERRQ(ierr);

        gidx = gen->starteqloc + 1;
        locglob = gen->startxpsetlocglob;
        row[0] = gidx;
        col[0] = locglob;
        val[0] = 1.0;

        ierr = MatSetValues(Je, 1, row, 1, col, val, ADD_VALUES);
        CHKERRQ(ierr);
      }
    }
  }
  ierr = VecRestoreArrayRead(X, &xarr);
  CHKERRQ(ierr);

  ierr = MatAssemblyBegin(Je, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Je, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);

  PetscLogFlops(flps);
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeInequalityConstraints_DCOPF(OPFLOW opflow, Vec X,
                                                        Vec Gi) {
  PetscErrorCode ierr;
  PetscInt i, k;
  PetscInt gloc;
  PetscInt xlocf, xloct;
  PetscScalar *g;
  PetscScalar Bdc, Pshift;
  PetscScalar thetaf, thetat, thetaft, thetatf;
  PetscScalar Pf, Pt;
  PS ps = opflow->ps;
  PSLINE line;
  PSBUS bus;
  PSBUS busf, bust;
  const PSBUS *connbuses;
  PSGEN gen;
  const PetscScalar *x;
  double flps = 0.0;
  PetscScalar Pg, delP, delPg, Pgs;
  PetscInt xloc;

  PetscFunctionBegin;
  ierr = VecSet(Gi, 0.0);
  CHKERRQ(ierr);

  ierr = VecGetArrayRead(X, &x);
  CHKERRQ(ierr);
  ierr = VecGetArray(Gi, &g);
  CHKERRQ(ierr);

  if (opflow->has_gensetpoint) {
    for (i = 0; i < ps->nbus; i++) {
      bus = &ps->bus[i];
      for (k = 0; k < bus->ngen; k++) {
        ierr = PSBUSGetGen(bus, k, &gen);
        CHKERRQ(ierr);
        if (!gen->status || gen->isrenewable)
          continue;

        gloc = gen->startineqloc;

        if (opflow->use_agc) {
          Pg = x[gen->startxpowloc];
          delPg = x[gen->startxpdevloc];
          delP = x[ps->startxloc];

          g[gloc] = (gen->apf * delP - delPg) * (Pg - gen->pt);
          g[gloc + 1] = (delPg - gen->apf * delP) * (gen->pb - Pg);

          flps += 8.0;
        }
      }
    }
  }

  if (!opflow->ignore_lineflow_constraints) {
    for (i = 0; i < ps->nline; i++) {
      line = &ps->line[i];
      if (!line->status || line->rateA > 1e5)
        continue;

      gloc = line->startineqloc;

      Bdc = line->bdc;
      Pshift = line->pshift;

      ierr = PSLINEGetConnectedBuses(line, &connbuses);
      CHKERRQ(ierr);
      busf = connbuses[0];
      bust = connbuses[1];

      xlocf = busf->startxVloc;
      xloct = bust->startxVloc;

      thetaf = x[xlocf];
      thetat = x[xloct];
      thetaft = thetaf - thetat;
      thetatf = thetat - thetaf;

      Pf = Bdc * thetaft + Pshift;
      Pt = Bdc * thetatf - Pshift;

      g[gloc] = Pf;
      g[gloc + 1] = Pt;

      flps += 4.0;
    }
  }

  ierr = VecRestoreArrayRead(X, &x);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(Gi, &g);
  CHKERRQ(ierr);

  PetscLogFlops(flps);
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeInequalityConstraintJacobian_DCOPF(OPFLOW opflow,
                                                               Vec X, Mat Ji) {
  PetscErrorCode ierr;
  PetscInt i, k, flps = 0;
  PetscInt row[2], col[4];
  PetscInt rstart, rend;
  PetscInt gloc, xlocf, xloct;
  PetscScalar val[4];
  PetscScalar Bdc, Pshift;
  PetscScalar thetaf, thetat, thetaft, thetatf;
  PetscScalar Pf, Pt;
  PetscScalar dPf_dthetaf, dPf_dthetat;
  PetscScalar dPt_dthetaf, dPt_dthetat;

  PS ps = opflow->ps;
  MPI_Comm comm = opflow->comm->type;
  PSLINE line;
  PSBUS busf, bust;
  const PSBUS *connbuses;
  PSBUS bus;
  PSGEN gen;
  const PetscScalar *x;
  PetscScalar Pg, delPg, delP;
  PetscInt xloc, loc;

  PetscFunctionBegin;
  ierr = MatZeroEntries(Ji);
  CHKERRQ(ierr);

  ierr = VecGetArrayRead(X, &x);
  CHKERRQ(ierr);

  ierr = MatGetOwnershipRange(Ji, &rstart, &rend);
  CHKERRQ(ierr);

  if (opflow->has_gensetpoint) {
    for (i = 0; i < ps->nbus; i++) {
      bus = &ps->bus[i];
      for (k = 0; k < bus->ngen; k++) {
        ierr = PSBUSGetGen(bus, k, &gen);
        CHKERRQ(ierr);
        if (!gen->status || gen->isrenewable)
          continue;

        gloc = gen->startineqloc;

        Pg = x[gen->startxpowloc];

        if (opflow->use_agc) {
          Pg = x[gen->startxpowloc];
          delP = x[ps->startxloc];
          delPg = x[gen->startxpdevloc];

          //	  g[gloc] = (gen->apf*delP - delPg)*(Pg - gen->pt);
          row[0] = gloc;
          col[0] = gen->startxpowloc;
          col[1] = gen->startxpdevloc;
          col[2] = ps->startxloc;
          val[0] = gen->apf * delP - delPg;
          val[1] = -(Pg - gen->pt);
          val[2] = gen->apf * (Pg - gen->pt);
          ierr = MatSetValues(Ji, 1, row, 3, col, val, ADD_VALUES);
          CHKERRQ(ierr);

          // g[gloc+1] = (delPg - gen->apf*delP)*(gen->pb - Pg);
          row[0] = gloc + 1;
          col[0] = gen->startxpowloc;
          col[1] = gen->startxpdevloc;
          col[2] = ps->startxloc;
          val[0] = gen->apf * delP - delPg;
          val[1] = gen->pb - Pg;
          val[2] = -gen->apf * (gen->pb - Pg);
          ierr = MatSetValues(Ji, 1, row, 3, col, val, ADD_VALUES);
          CHKERRQ(ierr);
        }
      }
    }
  }

  if (!opflow->ignore_lineflow_constraints) {
    for (i = 0; i < ps->nline; i++) {
      line = &ps->line[i];
      if (!line->status || line->rateA > 1e5)
        continue;

      gloc = line->startineqloc;

      Bdc = line->bdc;
      Pshift = line->pshift;

      ierr = PSLINEGetConnectedBuses(line, &connbuses);
      CHKERRQ(ierr);
      busf = connbuses[0];
      bust = connbuses[1];

      xlocf = busf->startxVloc;
      xloct = bust->startxVloc;

      thetaf = x[xlocf];
      thetat = x[xloct];
      thetaft = thetaf - thetat;
      thetatf = thetat - thetaf;

      dPf_dthetaf = Bdc;
      dPf_dthetat = -Bdc;

      dPt_dthetat = Bdc;
      dPt_dthetaf = -Bdc;

      row[0] = gloc;
      col[0] = xlocf;
      col[1] = xloct;

      val[0] = dPf_dthetaf;
      val[1] = dPf_dthetat;

      ierr = MatSetValues(Ji, 1, row, 2, col, val, ADD_VALUES);
      CHKERRQ(ierr);

      row[0] = gloc + 1;
      col[0] = xloct;
      col[1] = xlocf;

      val[0] = dPt_dthetat;
      val[1] = dPt_dthetaf;

      ierr = MatSetValues(Ji, 1, row, 2, col, val, ADD_VALUES);
      CHKERRQ(ierr);
    }
    flps += ps->nline * 2;
  }
  ierr = VecRestoreArrayRead(X, &x);
  CHKERRQ(ierr);

  ierr = MatAssemblyBegin(Ji, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Ji, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);

  PetscLogFlops(flps);
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeConstraints_DCOPF(OPFLOW opflow, Vec X, Vec G) {
  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeObjective_DCOPF(OPFLOW opflow, Vec X,
                                            PetscScalar *obj) {
  PetscErrorCode ierr;
  const PetscScalar *x;
  PS ps = opflow->ps;
  PetscInt i;
  PSBUS bus;
  PSGEN gen;
  PetscInt loc;
  PetscInt k;
  PetscScalar Pg;
  double flps = 0.0;

  PetscFunctionBegin;

  ierr = VecGetArrayRead(X, &x);
  CHKERRQ(ierr);

  *obj = 0.0;
  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    if (opflow->include_powerimbalance_variables) {
      loc = bus->startxpimbloc;
      PetscScalar Pimbplus, Pimbminus;
      loc = bus->startxpimbloc;
      Pimbplus = x[loc];
      Pimbminus = x[loc + 1];

      *obj +=
          opflow->powerimbalance_penalty * ps->MVAbase * (Pimbplus + Pimbminus);
      flps += 5.0;
    }

    for (k = 0; k < bus->ngen; k++) {
      ierr = PSBUSGetGen(bus, k, &gen);
      CHKERRQ(ierr);
      if (!gen->status)
        continue;

      if (opflow->objectivetype == MIN_GEN_COST) {
        loc = gen->startxpowloc;
        Pg = x[loc] * ps->MVAbase;
        *obj +=
            gen->cost_alpha * Pg * Pg + gen->cost_beta * Pg + gen->cost_gamma;
        flps += 7.0;
      } else if (opflow->objectivetype == MIN_GENSETPOINT_DEVIATION) {
        loc = gen->startxpdevloc;
        PetscScalar delPg;
        delPg = x[loc];

        *obj += delPg * delPg;
        flps += 3.0;
      }
    }

    if (opflow->include_loadloss_variables) {
      PSLOAD load;
      PetscScalar Pdloss;
      for (k = 0; k < bus->nload; k++) {
        ierr = PSBUSGetLoad(bus, k, &load);
        CHKERRQ(ierr);
        if (!load->status)
          continue;
        loc = load->startxloadlossloc;
        Pdloss = x[loc];
        *obj += opflow->loadloss_penalty * ps->MVAbase * Pdloss;
        flps += 3.0;
      }
    }
  }
  ierr = VecRestoreArrayRead(X, &x);
  CHKERRQ(ierr);

  PetscLogFlops(flps);
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeGradient_DCOPF(OPFLOW opflow, Vec X, Vec grad) {
  PetscErrorCode ierr;
  const PetscScalar *x;
  PetscScalar *df;
  PS ps = opflow->ps;
  PetscInt i;
  PSBUS bus;
  PSGEN gen;
  PSLOAD load;
  PetscInt loc;
  PetscInt k;
  PetscScalar Pg;
  double flps = 0.0;

  PetscFunctionBegin;
  ierr = VecGetArrayRead(X, &x);
  CHKERRQ(ierr);
  ierr = VecGetArray(grad, &df);
  CHKERRQ(ierr);

  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    if (opflow->include_powerimbalance_variables) {
      PetscScalar Pimb, Qimb;
      loc = bus->startxpimbloc;

      df[loc] = opflow->powerimbalance_penalty * ps->MVAbase;
      df[loc + 1] = opflow->powerimbalance_penalty * ps->MVAbase;

      flps += 2.0;
    }

    for (k = 0; k < bus->ngen; k++) {
      ierr = PSBUSGetGen(bus, k, &gen);
      CHKERRQ(ierr);
      if (!gen->status)
        continue;

      if (opflow->objectivetype == MIN_GEN_COST) {
        loc = gen->startxpowloc;
        Pg = x[loc] * ps->MVAbase;
        df[loc] = ps->MVAbase * (2 * gen->cost_alpha * Pg + gen->cost_beta);
        flps += 5.0;
      } else if (opflow->objectivetype == MIN_GENSETPOINT_DEVIATION) {
        PetscScalar delPg;
        loc = gen->startxpdevloc;
        delPg = x[loc];

        df[loc] = 2 * delPg;
        flps += 1.0;
      }
    }

    if (opflow->include_loadloss_variables) {
      PSLOAD load;
      PetscScalar Pdloss;
      for (k = 0; k < bus->nload; k++) {
        ierr = PSBUSGetLoad(bus, k, &load);
        CHKERRQ(ierr);
        if (!load->status)
          continue;
        loc = load->startxloadlossloc;
        Pdloss = x[loc];
        df[loc] = opflow->loadloss_penalty * ps->MVAbase;
        flps += 1.0;
      }
    }
  }
  ierr = VecRestoreArrayRead(X, &x);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(grad, &df);
  CHKERRQ(ierr);

  PetscLogFlops(flps);
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeObjandGradient_DCOPF(OPFLOW opflow, Vec X,
                                                 PetscScalar *obj, Vec Grad) {
  PetscErrorCode ierr;
  PetscFunctionBegin;

  ierr = OPFLOWComputeObjective_DCOPF(opflow, X, obj);
  CHKERRQ(ierr);
  ierr = OPFLOWComputeGradient_DCOPF(opflow, X, Grad);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/* Note: This function also sets the number of variables for each component,
 * i.e., bus, gen, load, and line */
PetscErrorCode OPFLOWModelSetNumVariables_DCOPF(OPFLOW opflow,
                                                PetscInt *busnvar,
                                                PetscInt *branchnvar,
                                                PetscInt *nx) {
  PetscInt i, ngen, nload, k;
  PS ps = opflow->ps;
  PSBUS bus;
  PSGEN gen;
  PSLOAD load;
  PSLINE line;
  PetscErrorCode ierr;
  PetscBool isghost;

  PetscFunctionBegin;

  *nx = 0;
  /* No variables for the branches */
  for (i = 0; i < ps->nline; i++) {
    line = &ps->line[i];
    branchnvar[i] = line->nx = 0;
    *nx += branchnvar[i];
  }

  /* Variables for the buses */
  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSIsGhosted(bus, &isghost);
    CHKERRQ(ierr);
    //    if(isghost) continue;
    busnvar[i] = 1; /* 1 variable for the bus (voltage angle) */
    bus->nxV = busnvar[i];

    if (opflow->include_powerimbalance_variables) {
      busnvar[i] += 2;
      bus->nxpimb = 2;
    }

    bus->nxshunt = 0;

    bus->nx = bus->nxV + bus->nxshunt + bus->nxpimb;

    if (opflow->use_agc) {
      if (bus->ide == REF_BUS) {
        ps->nx = 1; /* System deltaP used */
        busnvar[i] += 1;
      }
    }

    ierr = PSBUSGetNGen(bus, &ngen);
    CHKERRQ(ierr);
    for (k = 0; k < ngen; k++) {
      ierr = PSBUSGetGen(bus, k, &gen);
      CHKERRQ(ierr);
      if (!gen->status)
        continue;
      busnvar[i] += 1; /* (1 variable for generator power Pg) */
      gen->nxpow = 1;
      if (opflow->has_gensetpoint && !gen->isrenewable) {
        gen->nxpset = 1;
        gen->nxpdev = 1;
        busnvar[i] += 2;
      }
      gen->nx = gen->nxpow + gen->nxpdev;
    }

    if (opflow->include_loadloss_variables) {
      for (k = 0; k < bus->nload; k++) {
        ierr = PSBUSGetLoad(bus, k, &load);
        CHKERRQ(ierr);
        if (!load->status)
          continue;
        /* Load loss variables */
        busnvar[i] += 1;
        load->nxloadloss = 1;
        load->nx = load->nxloadloss;
      }
    }

    if (!isghost)
      *nx += busnvar[i];
  }

  PetscFunctionReturn(0);
}

/* Note: This function also sets the number of constraints for each component,
 * i.e., bus, gen, load, and line */
PetscErrorCode OPFLOWModelSetNumConstraints_DCOPF(OPFLOW opflow,
                                                  PetscInt *branchnconeq,
                                                  PetscInt *busnconeq,
                                                  PetscInt *nconeq,
                                                  PetscInt *nconineq) {
  PetscInt i, k, ngen;
  PS ps = opflow->ps;
  PSBUS bus;
  PSLINE line;
  PSGEN gen;
  PetscErrorCode ierr;
  PetscBool isghost;

  PetscFunctionBegin;
  *nconeq = *nconineq = 0;

  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSIsGhosted(bus, &isghost);
    CHKERRQ(ierr);
    if (isghost)
      continue;
    *nconeq += 1;
    bus->nconeq = 1;
    bus->nconineq = 0;

    if (opflow->use_agc) {
      if (bus->ide == REF_BUS) {
        ps->nconeq = ps->nconineq = 0;
      }
    }

    ierr = PSBUSGetNGen(bus, &ngen);
    CHKERRQ(ierr);
    for (k = 0; k < ngen; k++) {
      ierr = PSBUSGetGen(bus, k, &gen);
      CHKERRQ(ierr);
      if (!gen->status)
        continue;
      if (opflow->has_gensetpoint && !gen->isrenewable) {
        *nconeq += 2;
        gen->nconeq = 2;
        *nconineq += 0;
        gen->nconineq = 0;
        if (opflow->use_agc) {
          *nconineq += 2;
          gen->nconineq += 2;
        }
      }
    }
  }

  if (!opflow->ignore_lineflow_constraints) {
    for (i = 0; i < ps->nline; i++) {
      line = &ps->line[i];
      if (line->status && line->rateA < 1e5) {
        *nconineq += 2; /* Real power flow line constraint (from and to bus) */
        line->nconineq = 2;
        line->nconeq = 0;
      }
    }
  }

  PetscFunctionReturn(0);
}

/*
  OPFLOWComputeEqualityConstraintsHessian - Computes the Hessian for the
equality constraints function part

  Input Parameters:
+ opflow   - the OPFLOW object
. X        - solution vector X
- Lambda   - Lagrangian multiplier vector

  Output Parameters:
. H - the Hessian part for the equality constraints

*/
PetscErrorCode OPFLOWComputeEqualityConstraintsHessian_DCOPF(OPFLOW opflow,
                                                             Vec X, Vec Lambda,
                                                             Mat H) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}

/*
  OPFLOWComputeInequalityConstraintsHessian - Computes the Inequality
Constraints Hessian

  Input Parameters:
+ opflow   - the OPFLOW object
. X        - the solution vector
- Lambda   - Lagrangian multipler vector

  Output Parameters:
+ H   - the Hessian matrix

*/
PetscErrorCode OPFLOWComputeInequalityConstraintsHessian_DCOPF(OPFLOW opflow,
                                                               Vec X,
                                                               Vec Lambda,
                                                               Mat H) {
  PetscErrorCode ierr;
  PS ps = opflow->ps;
  PetscInt i, k, flps = 0;
  PSLINE line;
  const PSBUS *connbuses;
  PetscInt xloc;
  PetscInt xlocf, xloct, xlocglobf, xlocglobt;
  PSBUS busf, bust;
  const PetscScalar *x;
  const PetscScalar *lambda;
  PetscInt gloc;
  PetscInt row[12], col[12];
  PetscScalar val[12];
  PSBUS bus;
  PSGEN gen;
  PetscScalar Pg, delPg, delP;
  PetscInt xlocglob, loc;

  PetscFunctionBegin;

  ierr = VecGetArrayRead(X, &x);
  CHKERRQ(ierr);
  ierr = VecGetArrayRead(Lambda, &lambda);
  CHKERRQ(ierr);

  if (opflow->has_gensetpoint) {
    for (i = 0; i < ps->nbus; i++) {
      bus = &ps->bus[i];
      for (k = 0; k < bus->ngen; k++) {
        ierr = PSBUSGetGen(bus, k, &gen);
        CHKERRQ(ierr);
        if (!gen->status || gen->isrenewable)
          continue;

        gloc = gen->startineqloc;

        if (opflow->use_agc) {
          Pg = x[gen->startxpowloc];
          delPg = x[gen->startxpdevloc];
          delP = x[ps->startxloc];

          //	  df1_dPg = gen->apf*delP - delPg;
          // 	  df2_dPg = gen->apf*delP - delPg;

          row[0] = gen->startxpowloc;

          col[0] = gen->startxpowloc;
          col[1] = gen->startxpdevloc;
          col[2] = ps->startxloc;

          val[0] = 0.0;
          val[1] = -lambda[gloc] - lambda[gloc + 1];
          val[2] = gen->apf * (lambda[gloc] + lambda[gloc + 1]);

          ierr = MatSetValues(H, 1, row, 3, col, val, ADD_VALUES);

          //	  df1_ddelPg = -(Pg - gen->pt);
          //	  df2_ddelPg = gen->pb - Pg;
          row[0] = gen->startxpdevloc;
          val[0] = -lambda[gloc] - lambda[gloc + 1];
          ierr = MatSetValues(H, 1, row, 1, col, val, ADD_VALUES);
          CHKERRQ(ierr);

          //	  df1_ddelP = gen->apf*(Pg - gen->pt);
          //	  df2_ddelP = -gen->apf*(gen->pb - Pg);
          row[0] = ps->startxloc;
          val[0] = gen->apf * (lambda[gloc] + lambda[gloc + 1]);
          ierr = MatSetValues(H, 1, row, 1, col, val, ADD_VALUES);
          CHKERRQ(ierr);
        }
      }
    }
  }

  ierr = VecRestoreArrayRead(X, &x);
  CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(Lambda, &lambda);
  CHKERRQ(ierr);

  PetscLogFlops(flps);
  PetscFunctionReturn(0);
}

/*
  OPFLOWComputeObjectiveHessian - Computes the Hessian for the objective
function part

  Input Parameters:
+ opflow - the OPFLOW object
- X        - solution vecto X

  Output Parameters:
. H - the Hessian part for the objective function

*/
PetscErrorCode OPFLOWComputeObjectiveHessian_DCOPF(OPFLOW opflow, Vec X,
                                                   Mat H) {
  PetscErrorCode ierr;
  PS ps = opflow->ps;
  PetscInt i, k;
  PSBUS bus;
  PSGEN gen;
  PetscInt xlocglob;
  const PetscScalar *x;
  PetscInt row[2], col[2];
  PetscScalar val[2];
  PetscScalar obj_factor = opflow->obj_factor;
  PetscInt flps = 0;

  PetscFunctionBegin;

  ierr = VecGetArrayRead(X, &x);
  CHKERRQ(ierr);

  // for the part of objective
  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    for (k = 0; k < bus->ngen; k++) {
      ierr = PSBUSGetGen(bus, k, &gen);
      CHKERRQ(ierr);
      if (!gen->status)
        continue;

      if (opflow->objectivetype == MIN_GEN_COST) {
        xlocglob = gen->startxpowlocglob;

        row[0] = xlocglob;
        col[0] = xlocglob;

        val[0] = obj_factor * 2.0 * gen->cost_alpha * ps->MVAbase * ps->MVAbase;
        ierr = MatSetValues(H, 1, row, 1, col, val, ADD_VALUES);
        CHKERRQ(ierr);
        flps += 4;
      } else if (opflow->objectivetype == MIN_GENSETPOINT_DEVIATION) {
        xlocglob = gen->startxpdevlocglob;
        row[0] = xlocglob;
        col[0] = xlocglob;
        val[0] = obj_factor * 2.0;
        ierr = MatSetValues(H, 1, row, 1, col, val, ADD_VALUES);
        CHKERRQ(ierr);

        flps += 1;
      }
    }
  }

  ierr = VecRestoreArrayRead(X, &x);
  CHKERRQ(ierr);
  PetscLogFlops(flps);
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeHessian_DCOPF(OPFLOW opflow, Vec X, Vec Lambdae,
                                          Vec Lambdai, Mat H) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = MatZeroEntries(H);
  CHKERRQ(ierr);

  /* Objective function Hessian */
  ierr = OPFLOWComputeObjectiveHessian_DCOPF(opflow, X, H);
  CHKERRQ(ierr);

  if (opflow->nconineq) {
    ierr =
        OPFLOWComputeInequalityConstraintsHessian_DCOPF(opflow, X, Lambdai, H);
    CHKERRQ(ierr);
  }

  ierr = MatAssemblyBegin(H, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = MatAssemblyEnd(H, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolutionToPS_DCOPF(OPFLOW opflow) {
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
  PetscScalar Bdc, Pshift;
  PetscScalar thetaf, thetat, thetaft, thetatf;
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

    gloc = bus->starteqloc;
    bus->mult_pmis = lambdae[gloc];
    bus->mult_qmis = 0.0;

    if (opflow->include_powerimbalance_variables) {
      loc = bus->startxpimbloc;
      PetscScalar Pimbplus, Pimbminus;
      Pimbplus = x[loc];
      Pimbminus = x[loc + 1];
      bus->pimb = Pimbplus - Pimbminus;
      bus->qimb = 0.0;
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
      gen->qg = 0.0;
    }

    if (opflow->include_loadloss_variables) {
      for (k = 0; k < bus->nload; k++) {
        ierr = PSBUSGetLoad(bus, k, &load);
        CHKERRQ(ierr);
        loc = load->startxloadlossloc;
        load->pl = load->pl - x[loc];
      }
    }
  }

  for (i = 0; i < ps->nline; i++) {
    line = &ps->line[i];
    if (!line->status) {
      line->mult_sf = line->mult_st = 0.0;
      continue;
    }

    Bdc = line->bdc;
    Pshift = line->pshift;

    ierr = PSLINEGetConnectedBuses(line, &connbuses);
    CHKERRQ(ierr);
    busf = connbuses[0];
    bust = connbuses[1];

    xlocf = busf->startxVloc;
    xloct = bust->startxVloc;

    thetaf = x[xlocf];
    thetat = x[xloct];
    thetaft = thetaf - thetat;
    thetatf = thetat - thetaf;

    Pf = Bdc * thetaft + Pshift;
    Pt = Bdc * thetatf - Pshift;

    line->pf = Pf;
    line->qf = 0.0;
    line->pt = Pt;
    line->qt = 0.0;
    line->sf = Pf;
    line->st = Pt;

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

PetscErrorCode OPFLOWModelSetUp_DCOPF(OPFLOW opflow) {
  PetscErrorCode ierr;
  PS ps = (PS)opflow->ps;
  PSBUS bus;
  PSGEN gen;
  PSLOAD load;
  PSLINE line;
  PetscInt i, k;
  PetscInt loc, locglob;
  PetscInt eqloc = 0, ineqloc = 0;

  PetscFunctionBegin;

  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSGetVariableLocation(bus, &loc);
    CHKERRQ(ierr);
    ierr = PSBUSGetVariableGlobalLocation(bus, &locglob);
    CHKERRQ(ierr);

    bus->startxVloc = loc;
    bus->startxVlocglob = locglob;

    bus->startxshuntloc = bus->startxVloc + bus->nxV;
    bus->startxshuntlocglob = bus->startxVlocglob + bus->nxV;

    bus->startxpimbloc = bus->startxshuntloc + bus->nxshunt;
    bus->startxpimblocglob = bus->startxshuntlocglob + bus->nxshunt;

    bus->nx = bus->nxV + bus->nxshunt + bus->nxpimb;

    loc += bus->nx;
    locglob += bus->nx;

    bus->starteqloc = eqloc;
    eqloc += bus->nconeq;

    if (opflow->genbusvoltagetype == FIXED_WITHIN_QBOUNDS) {
      /* Inequality constraints */
      if (bus->ide == PV_BUS || bus->ide == REF_BUS) {
        bus->startineqloc = ineqloc;
        ineqloc += bus->nconineq;
      }
    }

    /* gen */
    for (k = 0; k < bus->ngen; k++) {
      ierr = PSBUSGetGen(bus, k, &gen);
      CHKERRQ(ierr);
      if (!gen->status)
        continue;

      gen->startxpowloc = loc;
      gen->startxpowlocglob = locglob;

      if (opflow->has_gensetpoint && !gen->isrenewable) {
        gen->startxpdevloc = gen->startxpowloc + gen->nxpow;
        gen->startxpdevlocglob = gen->startxpowlocglob + gen->nxpow;

        gen->startxpsetloc = gen->startxpdevloc + gen->nxpdev;
        gen->startxpsetlocglob = gen->startxpdevlocglob + gen->nxpdev;

        gen->starteqloc = eqloc;
        eqloc += gen->nconeq;
        gen->startineqloc = ineqloc;
        ineqloc += gen->nconineq;
      }

      gen->nx = gen->nxpow + gen->nxpdev + gen->nxpset;
      loc += gen->nx;
      locglob += gen->nx;
    }

    for (k = 0; k < bus->nload; k++) {
      ierr = PSBUSGetLoad(bus, k, &load);
      CHKERRQ(ierr);
      if (!load->status)
        continue;
      if (opflow->include_loadloss_variables) {
        load->startxloadlossloc = loc;
        load->startxloadlosslocglob = locglob;

        load->nx = load->nxloadloss;

        loc += load->nx;
        locglob += load->nx;
      }
    }

    if (opflow->use_agc) {
      if (bus->ide == REF_BUS) {
        ps->startxloc = loc;
        ps->startxlocglob = locglob;
        loc += ps->nx;
        locglob += ps->nx;
      }
    }
  }

  for (i = 0; i < ps->nline; i++) {
    line = &ps->line[i];
    if (line->status && line->rateA < 1e5) {
      line->startineqloc = ineqloc;
      ineqloc += line->nconineq;
    }
  }

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWModelCreate_DCOPF(OPFLOW opflow) {
  DCOPF dcopf;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = PetscCalloc1(1, &dcopf);
  CHKERRQ(ierr);

  opflow->model = dcopf;

  /* Inherit Ops */
  opflow->modelops.destroy = OPFLOWModelDestroy_DCOPF;
  opflow->modelops.setnumvariables = OPFLOWModelSetNumVariables_DCOPF;
  opflow->modelops.setnumconstraints = OPFLOWModelSetNumConstraints_DCOPF;
  opflow->modelops.setvariablebounds = OPFLOWSetVariableBounds_DCOPF;
  opflow->modelops.setconstraintbounds = OPFLOWSetConstraintBounds_DCOPF;
  opflow->modelops.setvariableandconstraintbounds =
      OPFLOWSetVariableandConstraintBounds_DCOPF;
  opflow->modelops.setinitialguess = OPFLOWSetInitialGuess_DCOPF;
  opflow->modelops.setup = OPFLOWModelSetUp_DCOPF;
  opflow->modelops.computeequalityconstraints =
      OPFLOWComputeEqualityConstraints_DCOPF;
  opflow->modelops.computeinequalityconstraints =
      OPFLOWComputeInequalityConstraints_DCOPF;
  opflow->modelops.computeequalityconstraintjacobian =
      OPFLOWComputeEqualityConstraintJacobian_DCOPF;
  opflow->modelops.computeinequalityconstraintjacobian =
      OPFLOWComputeInequalityConstraintJacobian_DCOPF;
  opflow->modelops.computehessian = OPFLOWComputeHessian_DCOPF;
  opflow->modelops.computeobjandgradient = OPFLOWComputeObjandGradient_DCOPF;
  opflow->modelops.computeobjective = OPFLOWComputeObjective_DCOPF;
  opflow->modelops.computegradient = OPFLOWComputeGradient_DCOPF;
  opflow->modelops.solutiontops = OPFLOWSolutionToPS_DCOPF;

  PetscFunctionReturn(0);
}
