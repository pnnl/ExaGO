#include "exago_config.h"
#include "ibcar.h"
#include <private/opflowimpl.h>

/**
    This is a revised model of the current-balance-cartesian (ibcar) model with
    the following differences:
    -- The line flows considered are based on "current", not the "MVA flow"
    -- Additional variables for the "current flow" on the lines (at to and from
ends) are introduced.
**/

PetscErrorCode OPFLOWModelDestroy_IBCAR2(OPFLOW opflow) {
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = PetscFree(opflow->model);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSetVariableBounds_IBCAR2(OPFLOW opflow, Vec Xl, Vec Xu) {
  PetscErrorCode ierr;
  PS ps = opflow->ps;
  PetscScalar *xl, *xu;
  PetscInt i;
  PSBUS bus;
  PSLINE line;
  PSGEN gen;
  PetscInt loc;
  PetscInt k;

  PetscFunctionBegin;

  /* Get array pointers */
  ierr = VecGetArray(Xl, &xl);
  CHKERRQ(ierr);
  ierr = VecGetArray(Xu, &xu);
  CHKERRQ(ierr);

  for (i = 0; i < ps->nline; i++) {
    line = &ps->line[i];

    ierr = PSLINEGetVariableLocation(line, &loc);
    CHKERRQ(ierr);
    /* Bounds on real and imaginary parts of line from and to current flows */
    xl[loc] = PETSC_NINFINITY;
    xu[loc] = PETSC_INFINITY;
    xl[loc + 1] = PETSC_NINFINITY;
    xu[loc + 1] = PETSC_INFINITY;
    xl[loc + 2] = PETSC_NINFINITY;
    xu[loc + 2] = PETSC_INFINITY;
    xl[loc + 3] = PETSC_NINFINITY;
    xu[loc + 3] = PETSC_INFINITY;
  }

  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    ierr = PSBUSGetVariableLocation(bus, &loc);
    CHKERRQ(ierr);

    /* Bounds on real and imaginary part of voltages */
    xl[loc] = -bus->Vmax;
    xu[loc] = bus->Vmax;
    xl[loc + 1] = -bus->Vmax;
    xu[loc + 1] = bus->Vmax;

    if (bus->ide == ISOLATED_BUS) {
      xl[loc] = xu[loc] = bus->vm * PetscCosScalar(bus->va);
      xl[loc + 1] = xu[loc + 1] = bus->vm * PetscSinScalar(bus->va);
    }

    for (k = 0; k < bus->ngen; k++) {
      ierr = PSBUSGetGen(bus, k, &gen);
      CHKERRQ(ierr);
      if (!gen->status)
        continue;
      loc = loc + 2;
      xl[loc] = gen->pb;     /* PGmin */
      xu[loc] = gen->pt;     /* PGmax */
      xl[loc + 1] = gen->qb; /* QGmin */
      xu[loc + 1] = gen->qt; /* QGmax */
      /* pb, pt, qb, qt are converted in p.u. in ps.c */
    }

    if (opflow->include_loadloss_variables) {
      for (k = 0; k < bus->nload; k++) {
        PSLOAD load;
        ierr = PSBUSGetLoad(bus, k, &load);
        CHKERRQ(ierr);
        loc += 2;
        if (!load->status)
          xl[loc] = xu[loc] = xl[loc + 1] = xu[loc + 1] = 0.0;
        else {
          xl[loc] = PetscMin(0.0, load->pl);
          xu[loc] = PetscMax(0.0, load->pl);
          xl[loc + 1] = PetscMin(0.0, load->ql);
          xu[loc + 1] = PetscMax(0.0, load->ql);
        }
      }
    }

    if (opflow->include_powerimbalance_variables) {
      loc += 2;
      xl[loc] = xl[loc + 1] = PETSC_NINFINITY;
      xu[loc] = xu[loc + 1] = PETSC_INFINITY;
    }
  }

  ierr = VecRestoreArray(Xl, &xl);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(Xu, &xu);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSetConstraintBounds_IBCAR2(OPFLOW opflow, Vec Gl, Vec Gu) {
  PetscErrorCode ierr;
  PS ps = opflow->ps;
  PetscScalar *gl, *gu;
  PetscInt i;
  PSLINE line;
  PSBUS bus;
  PetscInt gloc = 0;

  PetscFunctionBegin;

  /* Get array pointers */
  ierr = VecGetArray(Gl, &gl);
  CHKERRQ(ierr);
  ierr = VecGetArray(Gu, &gu);
  CHKERRQ(ierr);

  /* Equality constraint bounds for line flow balance equations */
  for (i = 0; i < ps->nline; i++) {
    line = &ps->line[i];
    if (!line->status)
      continue;
    gl[gloc] = gl[gloc + 1] = gl[gloc + 2] = gl[gloc + 3] = 0.0;
    gu[gloc] = gu[gloc + 1] = gu[gloc + 2] = gu[gloc + 3] = 0.0;

    gloc += 4;
  }
  /* Equality constraint bounds for balance equations and ref. bus angle */
  for (i = 0; i < ps->nbus; i++) {

    bus = &ps->bus[i];

    gl[gloc] = 0.0;
    gu[gloc] = 0.0;
    gl[gloc + 1] = 0.0;
    gu[gloc + 1] = 0.0;

    gloc += 2;
    if (bus->ide == REF_BUS) {
      /* Equality constraint on reference bus angle */
      /* V_I -V_R*tan(Va) */
      gl[gloc] = gu[gloc] = 0;
      gloc++;
    }
  }

  /* Inequality constraints on voltage magnitude */
  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    gl[gloc] = bus->Vmin * bus->Vmin;
    gu[gloc] = bus->Vmax * bus->Vmax;
    if (bus->ide == ISOLATED_BUS) {
      gl[gloc] = PETSC_NINFINITY;
      gu[gloc] = PETSC_INFINITY;
    }
    gloc++;
  }

  if (!opflow->ignore_lineflow_constraints) {
    /* Bounds on inequality constraints for line flows */
    for (i = 0; i < ps->nline; i++) {
      line = &ps->line[i];
      if (!line->status || line->rateA > 1e5)
        continue;
      gl[gloc] = gl[gloc + 1] = 0.0;
      gu[gloc] = gu[gloc + 1] =
          (line->rateA / ps->MVAbase) * (line->rateA / ps->MVAbase);
      gloc += 2;
    }
  }

  ierr = VecRestoreArray(Gl, &gl);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(Gu, &gu);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSetVariableandConstraintBounds_IBCAR2(OPFLOW opflow,
                                                           Vec Xl, Vec Xu,
                                                           Vec Gl, Vec Gu) {
  PetscErrorCode ierr;
  PetscFunctionBegin;

  ierr = OPFLOWSetVariableBounds_IBCAR2(opflow, Xl, Xu);
  CHKERRQ(ierr);
  ierr = OPFLOWSetConstraintBounds_IBCAR2(opflow, Gl, Gu);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSetInitialGuess_IBCAR2(OPFLOW opflow, Vec X, Vec Lambda) {
  PetscErrorCode ierr;
  PS ps = opflow->ps;
  const PetscScalar *xl, *xu;
  PetscScalar *x;
  PetscInt i;
  PSBUS bus;
  PSGEN gen;
  PSLINE line;
  PetscInt loc;
  PetscInt k;

  PetscFunctionBegin;

  /* Get array pointers */
  ierr = VecGetArray(X, &x);
  CHKERRQ(ierr);
  ierr = VecGetArrayRead(opflow->Xl, &xl);
  CHKERRQ(ierr);
  ierr = VecGetArrayRead(opflow->Xu, &xu);
  CHKERRQ(ierr);

  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    ierr = PSBUSGetVariableLocation(bus, &loc);
    CHKERRQ(ierr);

    if (bus->ide == ISOLATED_BUS) {
      x[loc] = bus->vm * PetscCosScalar(bus->va);
      x[loc + 1] = bus->vm * PetscSinScalar(bus->va);
    } else {
      if (opflow->initializationtype == OPFLOWINIT_MIDPOINT) {
        x[loc] = 0.5 * (bus->Vmin + bus->Vmax) * PetscCosScalar(bus->va);
        x[loc + 1] = 0; // 0.5*(bus->Vmin +
                        // bus->Vmax)*PetscCosScalar(bus->va);
      } else if (opflow->initializationtype == OPFLOWINIT_FROMFILE ||
                 opflow->initializationtype == OPFLOWINIT_ACPF) {
        x[loc] = PetscMax(bus->Vmin, PetscMin(bus->vm, bus->Vmax)) *
                 PetscCosScalar(bus->va);
        x[loc + 1] = PetscMax(bus->Vmin, PetscMin(bus->vm, bus->Vmax)) *
                     PetscSinScalar(bus->va);
      } else if (opflow->initializationtype == OPFLOWINIT_FLATSTART) {
        x[loc] = 1.0;
        x[loc + 1] = 0.0;
      }
    }

    for (k = 0; k < bus->ngen; k++) {
      ierr = PSBUSGetGen(bus, k, &gen);
      CHKERRQ(ierr);
      if (!gen->status)
        continue;
      loc = loc + 2;

      if (opflow->initializationtype == OPFLOWINIT_MIDPOINT ||
          opflow->initializationtype == OPFLOWINIT_FLATSTART) {
        x[loc] = 0.5 * (xl[loc] + xu[loc]);
        x[loc + 1] = 0.5 * (xl[loc + 1] + xu[loc + 1]);
      } else if (opflow->initializationtype == OPFLOWINIT_FROMFILE ||
                 opflow->initializationtype == OPFLOWINIT_ACPF) {
        x[loc] = PetscMax(gen->pb, PetscMin(gen->pg, gen->pt));
        x[loc + 1] = PetscMax(gen->qb, PetscMin(gen->qg, gen->qt));
      }
    }

    if (opflow->include_loadloss_variables) {
      for (k = 0; k < bus->nload; k++) {
        PSLOAD load;

        ierr = PSBUSGetLoad(bus, k, &load);
        CHKERRQ(ierr);
        loc += 2;
        /* Initial value for real and reactive power load loss */
        x[loc] = 0.0;
        x[loc + 1] = 0.0;
      }
    }

    if (opflow->include_powerimbalance_variables) {
      loc += 2;
      x[loc] = x[loc + 1] = 0.0;
    }
  }

  /* Initialize branch currents */
  for (i = 0; i < ps->nline; i++) {
    line = &ps->line[i];
    if (!line->status)
      continue;
    ierr = PSLINEGetVariableLocation(line, &loc);
    CHKERRQ(ierr);
    PetscScalar Gff, Bff, Gft, Bft, Gtf, Btf, Gtt, Btt;
    const PSBUS *connbuses;
    PSBUS busf, bust;
    PetscInt xlocf, xloct;
    PetscScalar Vrf, Vif, Vrt, Vit;

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

    ierr = PSBUSGetVariableLocation(busf, &xlocf);
    CHKERRQ(ierr);
    ierr = PSBUSGetVariableLocation(bust, &xloct);
    CHKERRQ(ierr);

    Vrf = x[xlocf];
    Vif = x[xlocf + 1];
    Vrt = x[xloct];
    Vit = x[xloct + 1];

    x[loc] = Gff * Vrf - Bff * Vif + Gft * Vrt - Bft * Vit;
    x[loc + 1] = Bff * Vrf + Gff * Vif + Bft * Vrt + Gft * Vit;

    x[loc + 2] = Gtt * Vrt - Btt * Vit + Gtf * Vrf - Btf * Vif;
    x[loc + 3] = Btt * Vrt + Gtt * Vit + Btf * Vrf + Gtf * Vif;
  }

  ierr = VecRestoreArray(X, &x);
  CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(opflow->Xl, &xl);
  CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(opflow->Xu, &xu);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeEqualityConstraints_IBCAR2(OPFLOW opflow, Vec X,
                                                       Vec Ge) {
  PetscErrorCode ierr;
  PetscInt i, k, nconnlines;
  PetscInt gloc = 0, row[6];
  PetscInt xloc, xlocf, xloct, xlocl;
  PetscScalar val[6];
  PetscScalar Pg, Qg, Pd, Qd;
  PetscScalar Gff, Bff, Gft, Bft, Gtf, Btf, Gtt, Btt;
  PetscScalar Vrf, Vif, Vrt, Vit;
  PetscScalar Irf, Iif, Irt, Iit;
  PetscScalar Vr, Vi, Vm, Vm2;
  PS ps = opflow->ps;
  PSLOAD load;
  PSLINE line;
  PSBUS bus, busf, bust;
  PSGEN gen;
  const PSBUS *connbuses;
  const PSLINE *connlines;
  const PetscScalar *x;
  PetscInt rstart;
  PetscInt flps = 0;

  PetscFunctionBegin;
  ierr = VecSet(Ge, 0.0);
  CHKERRQ(ierr);
  ierr = VecGetOwnershipRange(Ge, &rstart, NULL);
  CHKERRQ(ierr);

  ierr = VecGetArrayRead(X, &x);
  CHKERRQ(ierr);

  for (i = 0; i < ps->nline; i++) {
    line = &ps->line[i];
    if (!line->status)
      continue;

    ierr = PSLINEGetVariableLocation(line, &xloc);
    CHKERRQ(ierr);
    Irf = x[xloc];
    Iif = x[xloc + 1];
    Irt = x[xloc + 2];
    Iit = x[xloc + 3];

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

    ierr = PSBUSGetVariableLocation(busf, &xlocf);
    CHKERRQ(ierr);
    ierr = PSBUSGetVariableLocation(bust, &xloct);
    CHKERRQ(ierr);

    Vrf = x[xlocf];
    Vif = x[xlocf + 1];
    Vrt = x[xloct];
    Vit = x[xloct + 1];

    val[0] = Gff * Vrf - Bff * Vif + Gft * Vrt - Bft * Vit - Irf;
    val[1] = Bff * Vrf + Gff * Vif + Bft * Vrt + Gft * Vit - Iif;
    val[2] = Gtt * Vrt - Btt * Vit + Gtf * Vrf - Btf * Vif - Irt;
    val[3] = Btt * Vrt + Gtt * Vit + Btf * Vrf + Gtf * Vif - Iit;

    row[0] = rstart + gloc;
    row[1] = row[0] + 1;
    row[2] = row[1] + 1;
    row[3] = row[2] + 1;
    ierr = VecSetValues(Ge, 4, row, val, ADD_VALUES);
    CHKERRQ(ierr);

    gloc += 4;
  }
  flps += 32 * ps->nline;

  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    row[0] = rstart + gloc;
    row[1] = row[0] + 1;

    ierr = PSBUSGetVariableLocation(bus, &xloc);
    CHKERRQ(ierr);

    /* Real and imaginary voltages for the bus */
    Vr = x[xloc];
    Vi = x[xloc + 1];
    Vm = PetscSqrtScalar(Vr * Vr + Vi * Vi);
    Vm2 = Vm * Vm;
    flps += 4 + EXAGO_FLOPS_SQRTOP;

    if (bus->ide == ISOLATED_BUS) {
      val[0] = Vr - bus->vm * PetscCosScalar(bus->va);
      val[1] = Vi - bus->vm * PetscSinScalar(bus->va);

      ierr = VecSetValues(Ge, 2, row, val, ADD_VALUES);
      CHKERRQ(ierr);
      gloc += 2;
      flps += 8 + EXAGO_FLOPS_SINOP + EXAGO_FLOPS_COSOP;
      continue;
    }

    /* Shunt current injections */
    val[0] = bus->gl * Vr - bus->bl * Vi;
    val[1] = bus->bl * Vr + bus->gl * Vi;
    ierr = VecSetValues(Ge, 2, row, val, ADD_VALUES);
    CHKERRQ(ierr);
    flps += 6;

    /* Generation current injection */
    for (k = 0; k < bus->ngen; k++) {
      ierr = PSBUSGetGen(bus, k, &gen);
      CHKERRQ(ierr);
      if (!gen->status)
        continue;
      xloc += 2;
      Pg = x[xloc];
      Qg = x[xloc + 1];

      val[0] = -(Pg * Vr + Qg * Vi) / Vm2;
      val[1] = -(-Qg * Vr + Pg * Vi) / Vm2;
      ierr = VecSetValues(Ge, 2, row, val, ADD_VALUES);
      CHKERRQ(ierr);
    }
    flps += 8 * bus->ngen;

    /* Load current injection */
    for (k = 0; k < bus->nload; k++) {
      ierr = PSBUSGetLoad(bus, k, &load);
      CHKERRQ(ierr);
      if (opflow->include_loadloss_variables) {
        xloc += 2;
        Pd = load->pl - x[xloc];
        Qd = load->ql - x[xloc + 1];
        flps += 2;
      } else {
        Pd = load->pl;
        Qd = load->ql;
      }
      val[0] = (Pd * Vr + Qd * Vi) / Vm2;
      val[1] = (-Qd * Vr + Pd * Vi) / Vm2;
      flps += 8;

      ierr = VecSetValues(Ge, 2, row, val, ADD_VALUES);
      CHKERRQ(ierr);
    }

    /* Power imbalance current addition */
    if (opflow->include_powerimbalance_variables) {
      PetscScalar Pimb, Qimb;
      xloc += 2;
      Pimb = x[xloc];
      Qimb = x[xloc + 1];
      val[0] = (Pimb * Vr + Qimb * Vi) / Vm2;
      val[1] = (-Qimb * Vr + Pimb * Vi) / Vm2;
      ierr = VecSetValues(Ge, 2, row, val, ADD_VALUES);
      CHKERRQ(ierr);
      flps += 8;
    }

    /* Branch flow injections */
    ierr = PSBUSGetSupportingLines(bus, &nconnlines, &connlines);
    CHKERRQ(ierr);
    for (k = 0; k < nconnlines; k++) {
      line = connlines[k];
      if (!line->status)
        continue;

      ierr = PSLINEGetConnectedBuses(line, &connbuses);
      CHKERRQ(ierr);
      busf = connbuses[0];
      bust = connbuses[1];

      ierr = PSLINEGetVariableLocation(line, &xlocl);
      CHKERRQ(ierr);
      if (bus == busf) {
        Irf = x[xlocl];
        Iif = x[xlocl + 1];

        val[0] = Irf;
        val[1] = Iif;
        ierr = VecSetValues(Ge, 2, row, val, ADD_VALUES);
        CHKERRQ(ierr);
      } else {
        Irt = x[xlocl + 2];
        Iit = x[xlocl + 3];

        val[0] = Irt;
        val[1] = Iit;
        ierr = VecSetValues(Ge, 2, row, val, ADD_VALUES);
        CHKERRQ(ierr);
      }
    }
    gloc += 2;
    if (bus->ide == REF_BUS) {
      /* Equality constraint for angle */
      row[2] = row[1] + 1;
      val[2] = Vi - Vr * PetscTanScalar(bus->va);
      ierr = VecSetValues(Ge, 1, &row[2], &val[2], ADD_VALUES);
      CHKERRQ(ierr);
      gloc += 1;
      flps += 4 + EXAGO_FLOPS_TANOP;
    }
  }
  ierr = VecAssemblyBegin(Ge);
  CHKERRQ(ierr);
  ierr = VecAssemblyEnd(Ge);
  CHKERRQ(ierr);

  ierr = VecRestoreArrayRead(X, &x);
  CHKERRQ(ierr);
  ierr = PetscLogFlops(flps);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeEqualityConstraintJacobian_IBCAR2(OPFLOW opflow,
                                                              Vec X, Mat Je) {
  PetscErrorCode ierr;
  PetscInt i, k, row[3], col[4], genctr, gloc = 0, loadctr;
  PetscInt nconnlines, locglob, loc, locglobf, locglobt;
  PetscInt locl, locglobl;
  PetscScalar Vr, Vi, Vm, Vm2, Vm4, val[8], Gff, Bff, Gft, Bft, Gtf, Btf, Gtt,
      Btt;
  PS ps = opflow->ps;
  PSBUS bus;
  PSGEN gen;
  PSLINE line;
  PSBUS busf, bust;
  const PSLINE *connlines;
  const PSBUS *connbuses;
  const PetscScalar *x;
  PetscInt rstart;
  PetscInt flps = 0;

  PetscFunctionBegin;
  ierr = MatZeroEntries(Je);
  CHKERRQ(ierr);

  ierr = MatGetOwnershipRange(Je, &rstart, NULL);
  CHKERRQ(ierr);

  for (i = 0; i < ps->nline; i++) {
    line = &ps->line[i];
    if (!line->status)
      continue;

    ierr = PSLINEGetVariableLocation(line, &locl);
    CHKERRQ(ierr);
    ierr = PSLINEGetVariableGlobalLocation(line, &locglobl);
    CHKERRQ(ierr);

    Gff = line->yff[0];
    Bff = line->yff[1];
    Gft = line->yft[0];
    Bft = line->yft[1];
    Gtf = line->ytf[0];
    Btf = line->ytf[1];
    Gtt = line->ytt[0];
    Btt = line->ytt[1];

    /* Get the connected buses to this line */
    ierr = PSLINEGetConnectedBuses(line, &connbuses);
    CHKERRQ(ierr);
    busf = connbuses[0];
    bust = connbuses[1];

    ierr = PSBUSGetVariableGlobalLocation(busf, &locglobf);
    CHKERRQ(ierr);
    ierr = PSBUSGetVariableGlobalLocation(bust, &locglobt);
    CHKERRQ(ierr);

    row[0] = locglobl;
    row[1] = locglobl + 1;
    col[0] = locglobf;
    col[1] = locglobf + 1;
    col[2] = locglobt;
    col[3] = locglobt + 1;
    /* dIrf_dVrf */
    val[0] = Gff;
    /* dIrf_dVif */
    val[1] = -Bff;
    /* dIrf_dVrt */
    val[2] = Gft;
    /* dIrf_dVit */
    val[3] = -Bft;

    /* dIif_dVrf */
    val[4] = Bff;
    /* dIif_dVif */
    val[5] = Gff;
    /* dIif_dVrt */
    val[6] = Bft;
    /* dIif_dVit */
    val[7] = Gft;

    ierr = MatSetValues(Je, 2, row, 4, col, val, ADD_VALUES);
    CHKERRQ(ierr);

    col[0] = locglobl;
    col[1] = locglobl + 1;
    /* dIrf_dIrf */
    val[0] = -1.0;
    /* dIrf_dIif */
    val[1] = 0.0;
    /* dIif_dIrf */
    val[2] = 0.0;
    /* dIif_dIif */
    val[3] = -1.0;

    ierr = MatSetValues(Je, 2, row, 2, col, val, ADD_VALUES);
    CHKERRQ(ierr);

    row[0] = locglobl + 2;
    row[1] = locglobl + 3;
    col[0] = locglobt;
    col[1] = locglobt + 1;
    col[2] = locglobf;
    col[3] = locglobf + 1;

    /* dIrt_dVrt */
    val[0] = Gtt;
    /* dIrt_dVit */
    val[1] = -Btt;
    /* dIrt_dVrf */
    val[2] = Gtf;
    /* dIrt_dVif */
    val[3] = -Btf;

    /* dIit_dVrt */
    val[4] = Btt;
    /* dIit_dVit */
    val[5] = Gtt;
    /* dIit_dVrf */
    val[6] = Btf;
    /* dIit_dVif */
    val[7] = Gtf;

    ierr = MatSetValues(Je, 2, row, 4, col, val, ADD_VALUES);
    CHKERRQ(ierr);

    col[0] = locglobl + 2;
    col[1] = locglobl + 3;

    /* dIrt_dIrt */
    val[0] = -1.0;
    /* dIrt_dIit */
    val[1] = 0.0;
    /* dIit_dIrt */
    val[2] = 0.0;
    /* dIit_dIit */
    val[3] = -1.0;

    ierr = MatSetValues(Je, 2, row, 2, col, val, ADD_VALUES);
    CHKERRQ(ierr);

    gloc += 4;
  }

  ierr = VecGetArrayRead(X, &x);
  CHKERRQ(ierr);
  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    row[0] = rstart + gloc;
    row[1] = row[0] + 1;

    ierr = PSBUSGetVariableLocation(bus, &loc);
    CHKERRQ(ierr);
    ierr = PSBUSGetVariableGlobalLocation(bus, &locglob);
    CHKERRQ(ierr);

    Vr = x[loc];
    Vi = x[loc + 1];
    Vm = PetscSqrtScalar(Vr * Vr + Vi * Vi);
    Vm2 = Vm * Vm;
    Vm4 = Vm2 * Vm2;
    flps += 6 + EXAGO_FLOPS_SQRTOP;

    col[0] = locglob;
    col[1] = locglob + 1;
    /* Isolated and reference bus */
    if (bus->ide == ISOLATED_BUS) {
      val[0] = val[3] = 1.0;
      val[1] = val[2] = 0.0;
      ierr = MatSetValues(Je, 2, row, 2, col, val, ADD_VALUES);
      CHKERRQ(ierr);
      gloc += 2;
      continue;
    }
    /* Partial derivative for shunt contribution */
    val[0] = bus->gl;
    val[1] = -bus->bl;
    val[2] = bus->bl;
    val[3] = bus->gl;
    ierr = MatSetValues(Je, 2, row, 2, col, val, ADD_VALUES);
    CHKERRQ(ierr);

    genctr = 0;
    for (k = 0; k < bus->ngen; k++) {
      ierr = PSBUSGetGen(bus, k, &gen);
      CHKERRQ(ierr);
      if (!gen->status)
        continue;
      PetscScalar Pg, Qg;
      loc += 2;
      Pg = x[loc];
      Qg = x[loc + 1];

      col[0] = locglob;
      col[1] = locglob + 1;

      /* dIrg_dVr */
      val[0] = (Pg * Vr * Vr + 2 * Qg * Vr * Vi - Pg * Vi * Vi) / Vm4;
      /* dIrg_dVi */
      val[1] = (Qg * Vi * Vi + 2 * Pg * Vr * Vi - Qg * Vr * Vr) / Vm4;
      /* dIig_dVr */
      val[2] = (-Qg * Vr * Vr + 2 * Pg * Vr * Vi + Qg * Vi * Vi) / Vm4;
      /* dIig_dVi */
      val[3] = (Pg * Vi * Vi - 2 * Qg * Vr * Vi - Pg * Vr * Vr) / Vm4;

      ierr = MatSetValues(Je, 2, row, 2, col, val, ADD_VALUES);
      CHKERRQ(ierr);

      col[0] = locglob + 2 + genctr;
      col[1] = locglob + 2 + genctr + 1;

      /* dIrg_dPg */
      val[0] = -Vr / Vm2;
      /* dIrg_dQg */
      val[1] = -Vi / Vm2;
      /* dIig_dPg */
      val[2] = -Vi / Vm2;
      /* dIig_dQg */
      val[3] = Vr / Vm2;

      ierr = MatSetValues(Je, 2, row, 2, col, val, ADD_VALUES);
      CHKERRQ(ierr);

      genctr += 2;
    }
    flps += 44 * bus->ngen;

    loadctr = 0;
    for (k = 0; k < bus->nload; k++) {
      PSLOAD load;
      PetscScalar Pd, Qd, Pdloss, Qdloss;
      ierr = PSBUSGetLoad(bus, k, &load);
      CHKERRQ(ierr);
      if (opflow->include_loadloss_variables) {
        loc += 2;
        Pdloss = x[loc];
        Qdloss = x[loc + 1];

        Pd = load->pl - Pdloss;
        Qd = load->ql - Qdloss;
        flps += 2;
      } else {
        Pd = load->pl;
        Qd = load->ql;
      }

      col[0] = locglob;
      col[1] = locglob + 1;

      /* dIrd_dVr */
      val[0] = (-Pd * Vr * Vr - 2 * Qd * Vr * Vi + Pd * Vi * Vi) / Vm4;
      /* dIrd_dVi */
      val[1] = (-Qd * Vi * Vi - 2 * Pd * Vr * Vi + Qd * Vr * Vr) / Vm4;
      /* dIid_dVr */
      val[2] = (Qd * Vr * Vr - 2 * Pd * Vr * Vi - Qd * Vi * Vi) / Vm4;
      /* dIid_dVi */
      val[3] = (-Pd * Vi * Vi + 2 * Qd * Vr * Vi + Pd * Vr * Vr) / Vm4;
      flps += 40;

      ierr = MatSetValues(Je, 2, row, 2, col, val, ADD_VALUES);
      CHKERRQ(ierr);

      if (opflow->include_loadloss_variables) {

        col[0] = locglob + 2 + 2 * bus->ngen + loadctr;
        col[1] = locglob + 2 + 2 * bus->ngen + loadctr + 1;

        /* dIrd_dPdloss */
        val[0] = -Vr / Vm2;
        /* dIrd_dQdloss */
        val[1] = -Vi / Vm2;
        /* dIid_dPdloss */
        val[2] = -Vi / Vm2;
        /* dIid_dQdloss */
        val[3] = Vr / Vm2;

        ierr = MatSetValues(Je, 2, row, 2, col, val, ADD_VALUES);
        CHKERRQ(ierr);
        loadctr += 2;
        flps += 4;
      }
    }

    /* Power imbalance Jacobian terms */
    if (opflow->include_powerimbalance_variables) {
      PetscScalar Pimb, Qimb;
      loc += 2;
      Pimb = x[loc];
      Qimb = x[loc + 1];

      col[0] = locglob;
      col[1] = locglob + 1;

      /* dIrg_dVr */
      val[0] = (-Pimb * Vr * Vr - 2 * Qimb * Vr * Vi + Pimb * Vi * Vi) / Vm4;
      /* dIrg_dVi */
      val[1] = (-Qimb * Vi * Vi - 2 * Pimb * Vr * Vi + Qimb * Vr * Vr) / Vm4;
      /* dIig_dVr */
      val[2] = (Qimb * Vr * Vr - 2 * Pimb * Vr * Vi - Qimb * Vi * Vi) / Vm4;
      /* dIig_dVi */
      val[3] = (-Pimb * Vi * Vi + 2 * Qimb * Vr * Vi + Pimb * Vr * Vr) / Vm4;

      ierr = MatSetValues(Je, 2, row, 2, col, val, ADD_VALUES);
      CHKERRQ(ierr);

      col[0] = locglob + 2 + 2 * bus->ngen +
               opflow->include_loadloss_variables * 2 * bus->nload;
      col[1] = locglob + 2 + 2 * bus->ngen +
               opflow->include_loadloss_variables * 2 * bus->nload + 1;

      /* dIrmb_dPimb */
      val[0] = Vr / Vm2;
      /* dIrmb_dQimb */
      val[1] = Vi / Vm2;
      /* dIiimb_dPimb */
      val[2] = Vi / Vm2;
      /* dIiimb_dQimb */
      val[3] = -Vr / Vm2;
      flps += 44;

      ierr = MatSetValues(Je, 2, row, 2, col, val, ADD_VALUES);
      CHKERRQ(ierr);
    }

    /* Partial derivatives of network equations */
    /* Get the lines supporting the bus */
    ierr = PSBUSGetSupportingLines(bus, &nconnlines, &connlines);
    CHKERRQ(ierr);

    for (k = 0; k < nconnlines; k++) {
      line = connlines[k];
      if (!line->status)
        continue;

      /* Get the connected buses to this line */
      ierr = PSLINEGetConnectedBuses(line, &connbuses);
      CHKERRQ(ierr);
      busf = connbuses[0];
      bust = connbuses[1];

      ierr = PSLINEGetVariableGlobalLocation(line, &locglobl);
      CHKERRQ(ierr);

      if (bus == busf) {
        col[0] = locglobl;
        col[1] = locglobl + 1;
      } else {
        col[0] = locglobl + 2;
        col[1] = locglobl + 3;
      }
      val[0] = val[3] = 1.0;
      val[1] = val[2] = 0.0;
      ierr = MatSetValues(Je, 2, row, 2, col, val, ADD_VALUES);
      CHKERRQ(ierr);
    }
    gloc += 2;
    if (bus->ide == REF_BUS) {
      /* Partial of equality constraint for ref. bus angle */
      row[2] = row[1] + 1;
      col[0] = locglob;
      col[1] = locglob + 1;
      val[0] = -PetscTanScalar(bus->va);
      val[1] = 1.0;
      ierr = MatSetValues(Je, 1, row + 2, 2, col, val, ADD_VALUES);
      CHKERRQ(ierr);
      gloc += 1;
      flps += 2 + EXAGO_FLOPS_TANOP;
    }
  }
  ierr = VecRestoreArrayRead(X, &x);
  CHKERRQ(ierr);

  ierr = MatAssemblyBegin(Je, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Je, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = PetscLogFlops(flps);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeInequalityConstraints_IBCAR2(OPFLOW opflow, Vec X,
                                                         Vec Gi) {
  PetscErrorCode ierr;
  PetscInt i;
  PetscInt gloc = 0;
  PetscInt xloc;
  PetscScalar *g;
  PetscScalar Irf, Iif, Irt, Iit;
  PetscScalar Vr, Vi;
  PS ps = opflow->ps;
  PSBUS bus;
  PSLINE line;
  const PetscScalar *x;
  PetscInt flps = 0;

  PetscFunctionBegin;
  ierr = VecSet(Gi, 0.0);
  CHKERRQ(ierr);

  ierr = VecGetArrayRead(X, &x);
  CHKERRQ(ierr);
  ierr = VecGetArray(Gi, &g);
  CHKERRQ(ierr);

  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSGetVariableLocation(bus, &xloc);
    CHKERRQ(ierr);

    Vr = x[xloc];
    Vi = x[xloc + 1];

    g[gloc] = Vr * Vr + Vi * Vi;
    gloc++;
  }
  flps += 3 * ps->nbus;

  if (!opflow->ignore_lineflow_constraints) {
    for (i = 0; i < ps->nline; i++) {
      line = &ps->line[i];
      if (!line->status || line->rateA > 1e5)
        continue;

      ierr = PSLINEGetVariableLocation(line, &xloc);
      CHKERRQ(ierr);

      Irf = x[xloc];
      Iif = x[xloc + 1];
      Irt = x[xloc + 2];
      Iit = x[xloc + 3];

      g[gloc] = Irf * Irf + Iif * Iif;
      g[gloc + 1] = Irt * Irt + Iit * Iit;

      gloc += 2;
    }
    flps += 6 * ps->nline;
  }

  ierr = VecRestoreArrayRead(X, &x);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(Gi, &g);
  CHKERRQ(ierr);
  ierr = PetscLogFlops(flps);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeInequalityConstraintJacobian_IBCAR2(OPFLOW opflow,
                                                                Vec X, Mat Ji) {
  PetscErrorCode ierr;
  PetscInt i;
  PetscInt row[2], col[4];
  PetscInt rstart, rend;
  PetscInt gloc = 0, xloc, xlocglob;
  PetscScalar val[4];
  PetscScalar Irf, Iif, Irt, Iit;
  PetscScalar Vr, Vi;
  PS ps = opflow->ps;
  MPI_Comm comm = opflow->comm->type;
  PSLINE line;
  PSBUS bus;
  const PetscScalar *x;
  PetscInt flps = 0;

  PetscFunctionBegin;
  ierr = MatZeroEntries(Ji);
  CHKERRQ(ierr);

  ierr = VecGetArrayRead(X, &x);
  CHKERRQ(ierr);

  ierr = MatGetOwnershipRange(Ji, &rstart, &rend);
  CHKERRQ(ierr);
  ierr = MPI_Barrier(comm);
  CHKERRQ(ierr);

  gloc = rstart;
  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    row[0] = gloc;

    ierr = PSBUSGetVariableLocation(bus, &xloc);
    CHKERRQ(ierr);
    ierr = PSBUSGetVariableGlobalLocation(bus, &xlocglob);
    CHKERRQ(ierr);

    Vr = x[xloc];
    Vi = x[xloc + 1];

    col[0] = xlocglob;
    col[1] = xlocglob + 1;
    val[0] = 2 * Vr;
    val[1] = 2 * Vi;

    ierr = MatSetValues(Ji, 1, row, 2, col, val, ADD_VALUES);
    CHKERRQ(ierr);

    gloc += 1;
  }
  flps += 2 * ps->nbus;

  if (opflow->ignore_lineflow_constraints) {
    ierr = VecRestoreArrayRead(X, &x);
    CHKERRQ(ierr);

    ierr = MatAssemblyBegin(Ji, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(Ji, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);

    PetscFunctionReturn(0);
  }

  for (i = 0; i < ps->nline; i++) {
    line = &ps->line[i];
    if (!line->status || line->rateA > 1e5)
      continue;

    ierr = PSLINEGetVariableLocation(line, &xloc);
    CHKERRQ(ierr);
    ierr = PSLINEGetVariableGlobalLocation(line, &xlocglob);
    CHKERRQ(ierr);

    row[0] = gloc;
    col[0] = xlocglob;
    col[1] = xlocglob + 1;

    Irf = x[xloc];
    Iif = x[xloc + 1];

    /* dIf2_dIrf */
    val[0] = 2 * Irf;
    /* dIf_dIif */
    val[1] = 2 * Iif;

    ierr = MatSetValues(Ji, 1, row, 2, col, val, ADD_VALUES);
    CHKERRQ(ierr);

    row[0] = gloc + 1;
    col[0] = xlocglob + 2;
    col[1] = xlocglob + 3;

    Irt = x[xloc + 2];
    Iit = x[xloc + 3];

    /* dIt2_dIrt */
    val[0] = 2 * Irt;
    /* dIt2_dIit */
    val[1] = 2 * Iit;

    ierr = MatSetValues(Ji, 1, row, 2, col, val, ADD_VALUES);
    CHKERRQ(ierr);

    gloc += 2;
  }
  flps += 4 * ps->nline;

  ierr = VecRestoreArrayRead(X, &x);
  CHKERRQ(ierr);

  ierr = MatAssemblyBegin(Ji, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Ji, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = PetscLogFlops(flps);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeConstraints_IBCAR2(OPFLOW opflow, Vec X, Vec G) {
  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeObjective_IBCAR2(OPFLOW opflow, Vec X,
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
  PetscInt flps = 0;

  PetscFunctionBegin;

  ierr = VecGetArrayRead(X, &x);
  CHKERRQ(ierr);

  *obj = 0.0;
  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    ierr = PSBUSGetVariableLocation(bus, &loc);
    CHKERRQ(ierr);

    for (k = 0; k < bus->ngen; k++) {
      ierr = PSBUSGetGen(bus, k, &gen);
      CHKERRQ(ierr);
      if (!gen->status)
        continue;
      loc = loc + 2;
      if (opflow->obj_gencost) {
        Pg = x[loc] * ps->MVAbase;
        *obj +=
            gen->cost_alpha * Pg * Pg + gen->cost_beta * Pg + gen->cost_gamma;
        flps += 7;
      }
    }

    if (opflow->include_loadloss_variables) {
      PSLOAD load;
      PetscScalar Pdloss, Qdloss;
      for (k = 0; k < bus->nload; k++) {
        loc += 2;
        ierr = PSBUSGetLoad(bus, k, &load);
        CHKERRQ(ierr);
        if (!load->status)
          continue;
        Pdloss = x[loc];
        Qdloss = x[loc + 1];
        *obj += opflow->loadloss_penalty * ps->MVAbase * ps->MVAbase *
                (Pdloss * Pdloss + Qdloss * Qdloss);
        flps += 7;
      }
    }

    if (opflow->include_powerimbalance_variables) {
      PetscScalar Pimb, Qimb;
      loc += 2;
      Pimb = x[loc];
      Qimb = x[loc + 1];
      *obj += opflow->powerimbalance_penalty * ps->MVAbase * ps->MVAbase *
              (Pimb * Pimb + Qimb * Qimb);
      flps += 7;
    }
  }
  ierr = VecRestoreArrayRead(X, &x);
  CHKERRQ(ierr);
  ierr = PetscLogFlops(flps);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeGradient_IBCAR2(OPFLOW opflow, Vec X, Vec grad) {
  PetscErrorCode ierr;
  const PetscScalar *x;
  PetscScalar *df;
  PS ps = opflow->ps;
  PetscInt i;
  PSBUS bus;
  PSGEN gen;
  PetscInt loc;
  PetscInt k;
  PetscScalar Pg;
  PetscInt flps = 0;

  PetscFunctionBegin;
  ierr = VecGetArrayRead(X, &x);
  CHKERRQ(ierr);
  ierr = VecGetArray(grad, &df);
  CHKERRQ(ierr);

  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    ierr = PSBUSGetVariableLocation(bus, &loc);
    CHKERRQ(ierr);

    for (k = 0; k < bus->ngen; k++) {
      ierr = PSBUSGetGen(bus, k, &gen);
      CHKERRQ(ierr);
      if (!gen->status)
        continue;
      loc = loc + 2;
      if (opflow->obj_gencost) {
        Pg = x[loc] * ps->MVAbase;
        df[loc] = ps->MVAbase * (2 * gen->cost_alpha * Pg + gen->cost_beta);
        flps += 5;
      }
    }

    if (opflow->include_loadloss_variables) {
      PSLOAD load;
      PetscScalar Pdloss, Qdloss;
      for (k = 0; k < bus->nload; k++) {
        loc += 2;
        ierr = PSBUSGetLoad(bus, k, &load);
        CHKERRQ(ierr);
        if (!load->status)
          continue;
        Pdloss = x[loc];
        Qdloss = x[loc + 1];
        df[loc] =
            opflow->loadloss_penalty * ps->MVAbase * ps->MVAbase * 2 * Pdloss;
        df[loc + 1] =
            opflow->loadloss_penalty * ps->MVAbase * ps->MVAbase * 2 * Qdloss;
        flps += 8;
      }
    }

    if (opflow->include_powerimbalance_variables) {
      PetscScalar Pimb, Qimb;
      loc += 2;
      Pimb = x[loc];
      Qimb = x[loc + 1];
      df[loc] =
          opflow->powerimbalance_penalty * ps->MVAbase * ps->MVAbase * 2 * Pimb;
      df[loc + 1] =
          opflow->powerimbalance_penalty * ps->MVAbase * ps->MVAbase * 2 * Qimb;
      flps += 8;
    }
  }
  ierr = VecRestoreArrayRead(X, &x);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(grad, &df);
  CHKERRQ(ierr);
  ierr = PetscLogFlops(flps);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeObjandGradient_IBCAR2(OPFLOW opflow, Vec X,
                                                  PetscScalar *obj, Vec Grad) {
  PetscErrorCode ierr;
  PetscFunctionBegin;

  ierr = OPFLOWComputeObjective_IBCAR2(opflow, X, obj);
  CHKERRQ(ierr);
  ierr = OPFLOWComputeGradient_IBCAR2(opflow, X, Grad);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWModelSetNumVariables_IBCAR2(OPFLOW opflow,
                                                 PetscInt *busnvar,
                                                 PetscInt *branchnvar,
                                                 PetscInt *nx) {
  PetscInt i, ngen, nload, k;
  PS ps = opflow->ps;
  PSBUS bus;
  PSLINE line;
  PSGEN gen;
  PetscErrorCode ierr;
  PetscBool isghost;

  PetscFunctionBegin;

  *nx = 0;
  /* No variables for the branches */
  for (i = 0; i < ps->nline; i++) {
    line = &ps->line[i];
    if (!line->status)
      branchnvar[i] = 0;
    else
      branchnvar[i] = 4; /* Ifr, Ifi, Itr, Iti From and to bus currents (real
                            and imaginary parts)*/
    *nx += branchnvar[i];
  }

  /* Variables for the buses */
  /* Real and imaginary part of voltages at each bus + Pg,Qg for each gen */
  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSIsGhosted(bus, &isghost);
    CHKERRQ(ierr);
    busnvar[i] = 2; /* 2 variables for the voltages */
    ierr = PSBUSGetNGen(bus, &ngen);
    CHKERRQ(ierr);
    for (k = 0; k < bus->ngen; k++) {
      ierr = PSBUSGetGen(bus, k, &gen);
      CHKERRQ(ierr);
      if (!gen->status)
        continue;
      busnvar[i] += 2; /* (2 variables Pg, Qg for each gen) */
    }

    if (opflow->include_loadloss_variables) {
      ierr = PSBUSGetNLoad(bus, &nload);
      CHKERRQ(ierr);
      /* Load loss variables..Real and imaginary part of the load loss */
      busnvar[i] += 2 * nload;
    }

    if (opflow->include_powerimbalance_variables)
      busnvar[i] += 2;
    if (!isghost)
      *nx += busnvar[i];
  }

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWModelSetNumConstraints_IBCAR2(OPFLOW opflow,
                                                   PetscInt *branchnconeq,
                                                   PetscInt *busnconeq,
                                                   PetscInt *nconeq,
                                                   PetscInt *nconineq) {
  PetscInt i;
  PS ps = opflow->ps;
  PSBUS bus;
  PSLINE line;
  PetscErrorCode ierr;
  PetscBool isghost;

  PetscFunctionBegin;

  *nconeq = 0;
  *nconineq = 0;

  for (i = 0; i < ps->nline; i++) {
    line = &ps->line[i];
    if (line->status)
      *nconeq += 4;
  }

  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSIsGhosted(bus, &isghost);
    CHKERRQ(ierr);
    if (isghost)
      continue;
    *nconeq += 2; /* Power balance constraints */
    if (bus->ide == REF_BUS)
      *nconeq += 1; /* Reference angle constraint */
    *nconineq += 1; /* Voltage magnitude constraint */
  }

  if (!opflow->ignore_lineflow_constraints) {
    for (i = 0; i < ps->nline; i++) {
      line = &ps->line[i];
      if (line->status && line->rateA < 1e5)
        *nconineq += 2; /* Line flow constraints */
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
PetscErrorCode OPFLOWComputeEqualityConstraintsHessian_IBCAR2(OPFLOW opflow,
                                                              Vec X, Vec Lambda,
                                                              Mat H) {
  PetscErrorCode ierr;
  PS ps = opflow->ps;
  PetscInt i, k;
  PSBUS bus;
  PSGEN gen;
  PSLINE line;
  PetscInt xloc, xlocglob;
  const PetscScalar *x;
  const PetscScalar *lambda;
  PetscInt gloc = 0;
  PetscInt row[12], col[12];
  PetscScalar val[12];
  PetscScalar Vr, Vi, Vm2, Vm4, Vm6;
  PetscInt flps = 0;

  PetscFunctionBegin;

  ierr = VecGetArrayRead(X, &x);
  CHKERRQ(ierr);
  ierr = VecGetArrayRead(Lambda, &lambda);
  CHKERRQ(ierr);

  for (i = 0; i < ps->nline; i++) {
    line = &ps->line[i];
    if (!line->status)
      continue;
    gloc += 4; /* No Hessian contribution from line flow constraints, only
                  increment gloc counter */
  }

  // For equality constraints (power flow) */
  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    ierr = PSBUSGetVariableLocation(bus, &xloc);
    CHKERRQ(ierr);
    ierr = PSBUSGetVariableGlobalLocation(bus, &xlocglob);
    CHKERRQ(ierr);

    Vr = x[xloc];
    Vi = x[xloc + 1];
    Vm2 = Vr * Vr + Vi * Vi;
    Vm4 = Vm2 * Vm2;
    Vm6 = Vm2 * Vm2 * Vm2;
    flps += 6;

    PetscInt genctr = 0;
    for (k = 0; k < bus->ngen; k++) {
      ierr = PSBUSGetGen(bus, k, &gen);
      CHKERRQ(ierr);
      if (!gen->status)
        continue;

      PetscScalar Pg, Qg;
      PetscScalar d2Irg_dVr2, d2Irg_dVr_dVi, d2Irg_dVr_dPg, d2Irg_dVr_dQg;
      PetscScalar d2Irg_dVi_dVr, d2Irg_dVi2, d2Irg_dVi_dPg, d2Irg_dVi_dQg;
      PetscScalar d2Iig_dVr2, d2Iig_dVr_dVi, d2Iig_dVr_dPg, d2Iig_dVr_dQg;
      PetscScalar d2Iig_dVi_dVr, d2Iig_dVi2, d2Iig_dVi_dPg, d2Iig_dVi_dQg;

      xloc += 2;
      Pg = x[xloc];
      Qg = x[xloc + 1];

      d2Irg_dVr2 = -2 *
                   (Pg * Vr * Vr * Vr + 3 * Qg * Vr * Vr * Vi -
                    3 * Pg * Vr * Vi * Vi - Qg * Vi * Vi * Vi) /
                   Vm6;
      d2Irg_dVr_dVi = 2 *
                      (Pg * Vi * Vi * Vi - 3 * Qg * Vr * Vi * Vi -
                       3 * Pg * Vr * Vr * Vi + Qg * Vr * Vr * Vr) /
                      Vm6;
      d2Irg_dVr_dPg = (Vr * Vr - Vi * Vi) / Vm4;
      d2Irg_dVr_dQg = 2 * Vr * Vi / Vm4;

      d2Irg_dVi_dVr = 2 *
                      (Qg * Vr * Vr * Vr - 3 * Pg * Vi * Vr * Vr -
                       3 * Qg * Vi * Vi * Vr + Pg * Vi * Vi * Vi) /
                      Vm6;
      d2Irg_dVi2 = -2 *
                   (Qg * Vi * Vi * Vi + 3 * Pg * Vr * Vi * Vi -
                    3 * Qg * Vr * Vr * Vi - Pg * Vr * Vr * Vr) /
                   Vm6;
      d2Irg_dVi_dPg = 2 * Vr * Vi / Vm4;
      d2Irg_dVi_dQg = (Vi * Vi - Vr * Vr) / Vm4;

      d2Iig_dVr2 = 2 *
                   (Qg * Vr * Vr * Vr - 3 * Pg * Vi * Vr * Vr -
                    3 * Qg * Vi * Vi * Vr + Pg * Vi * Vi * Vi) /
                   Vm6;
      d2Iig_dVr_dVi = -2 *
                      (Qg * Vi * Vi * Vi + 3 * Pg * Vr * Vi * Vi -
                       3 * Qg * Vr * Vr * Vi - Pg * Vr * Vr * Vr) /
                      Vm6;
      d2Iig_dVr_dPg = 2 * Vr * Vi / Vm4;
      d2Iig_dVr_dQg = (Vi * Vi - Vr * Vr) / Vm4;

      d2Iig_dVi_dVr = 2 *
                      (Pg * Vr * Vr * Vr + 3 * Qg * Vr * Vr * Vi -
                       3 * Pg * Vr * Vi * Vi - Qg * Vi * Vi * Vi) /
                      Vm6;
      d2Iig_dVi2 = -2 *
                   (Pg * Vi * Vi * Vi - 3 * Qg * Vr * Vi * Vi -
                    3 * Pg * Vr * Vr * Vi + Qg * Vr * Vr * Vr) /
                   Vm6;
      d2Iig_dVi_dPg = (Vi * Vi - Vr * Vr) / Vm4;
      d2Iig_dVi_dQg = -2 * Vr * Vi / Vm4;

      row[0] = xlocglob;
      row[1] = xlocglob + 1;
      col[0] = xlocglob;
      col[1] = xlocglob + 1;
      col[2] = xlocglob + 2 + genctr;
      col[3] = xlocglob + 2 + genctr + 1;

      val[0] = val[1] = val[2] = val[3] = val[4] = val[5] = val[6] = val[7] =
          0.0;

      val[0] = lambda[gloc] * d2Irg_dVr2 + lambda[gloc + 1] * d2Iig_dVr2;
      val[1] = lambda[gloc] * d2Irg_dVr_dVi + lambda[gloc + 1] * d2Iig_dVr_dVi;
      val[2] = lambda[gloc] * d2Irg_dVr_dPg + lambda[gloc + 1] * d2Iig_dVr_dPg;
      val[3] = lambda[gloc] * d2Irg_dVr_dQg + lambda[gloc + 1] * d2Iig_dVr_dQg;

      val[4] = lambda[gloc] * d2Irg_dVi_dVr + lambda[gloc + 1] * d2Iig_dVi_dVr;
      val[5] = lambda[gloc] * d2Irg_dVi2 + lambda[gloc + 1] * d2Iig_dVi2;
      val[6] = lambda[gloc] * d2Irg_dVi_dPg + lambda[gloc + 1] * d2Iig_dVi_dPg;
      val[7] = lambda[gloc] * d2Irg_dVi_dQg + lambda[gloc + 1] * d2Iig_dVi_dQg;

      ierr = MatSetValues(H, 2, row, 4, col, val, ADD_VALUES);
      CHKERRQ(ierr);

      PetscScalar d2Irg_dPg_dVr, d2Irg_dPg_dVi, d2Irg_dQg_dVr, d2Irg_dQg_dVi;
      PetscScalar d2Iig_dPg_dVr, d2Iig_dPg_dVi, d2Iig_dQg_dVr, d2Iig_dQg_dVi;

      row[0] = xlocglob + 2 + genctr;
      row[1] = xlocglob + 2 + genctr + 1;
      col[0] = xlocglob;
      col[1] = xlocglob + 1;

      d2Irg_dPg_dVr = (Vr * Vr - Vi * Vi) / Vm4;
      d2Irg_dPg_dVi = 2 * Vr * Vi / Vm4;
      d2Irg_dQg_dVr = 2 * Vr * Vi / Vm4;
      d2Irg_dQg_dVi = (Vi * Vi - Vr * Vr) / Vm4;

      d2Iig_dPg_dVr = 2 * Vr * Vi / Vm4;
      d2Iig_dPg_dVi = (Vi * Vi - Vr * Vr) / Vm4;

      d2Iig_dQg_dVr = (Vi * Vi - Vr * Vr) / Vm4;
      d2Iig_dQg_dVi = -2 * Vr * Vi / Vm4;

      val[0] = val[1] = val[2] = val[3] = 0.0;

      val[0] = lambda[gloc] * d2Irg_dPg_dVr + lambda[gloc + 1] * d2Iig_dPg_dVr;
      val[1] = lambda[gloc] * d2Irg_dPg_dVi + lambda[gloc + 1] * d2Iig_dPg_dVi;
      val[2] = lambda[gloc] * d2Irg_dQg_dVr + lambda[gloc + 1] * d2Iig_dQg_dVr;
      val[3] = lambda[gloc] * d2Irg_dQg_dVi + lambda[gloc + 1] * d2Iig_dQg_dVi;

      ierr = MatSetValues(H, 2, row, 2, col, val, ADD_VALUES);
      CHKERRQ(ierr);

      genctr += 2;
    }
    flps += 245 * bus->ngen;

    PetscInt loadctr = 0;
    for (k = 0; k < bus->nload; k++) {
      PSLOAD load;
      PetscScalar Pd, Qd, Pdloss, Qdloss;
      PetscScalar d2Ird_dVr2, d2Ird_dVr_dVi, d2Ird_dVr_dPdloss,
          d2Ird_dVr_dQdloss;
      PetscScalar d2Ird_dVi_dVr, d2Ird_dVi2, d2Ird_dVi_dPdloss,
          d2Ird_dVi_dQdloss;
      PetscScalar d2Iid_dVr2, d2Iid_dVr_dVi, d2Iid_dVr_dPdloss,
          d2Iid_dVr_dQdloss;
      PetscScalar d2Iid_dVi_dVr, d2Iid_dVi2, d2Iid_dVi_dPdloss,
          d2Iid_dVi_dQdloss;

      ierr = PSBUSGetLoad(bus, k, &load);
      CHKERRQ(ierr);
      if (opflow->include_loadloss_variables) {
        xloc += 2;
        Pdloss = x[xloc];
        Qdloss = x[xloc + 1];
        Pd = load->pl - Pdloss;
        Qd = load->ql - Qdloss;
        flps += 2;
      } else {
        Pd = load->pl;
        Qd = load->ql;
      }

      d2Ird_dVr2 = 2 *
                   (Pd * Vr * Vr * Vr + 3 * Qd * Vr * Vr * Vi -
                    3 * Pd * Vr * Vi * Vi - Qd * Vi * Vi * Vi) /
                   Vm6;
      d2Ird_dVr_dVi = -2 *
                      (Pd * Vi * Vi * Vi - 3 * Qd * Vr * Vi * Vi -
                       3 * Pd * Vr * Vr * Vi + Qd * Vr * Vr * Vr) /
                      Vm6;

      d2Ird_dVr_dPdloss = (Vr * Vr - Vi * Vi) / Vm4;
      d2Ird_dVr_dQdloss = 2 * Vr * Vi / Vm4;

      d2Ird_dVi_dVr = -2 *
                      (Qd * Vr * Vr * Vr - 3 * Pd * Vi * Vr * Vr -
                       3 * Qd * Vi * Vi * Vr + Pd * Vi * Vi * Vi) /
                      Vm6;
      d2Ird_dVi2 = 2 *
                   (Qd * Vi * Vi * Vi + 3 * Pd * Vr * Vi * Vi -
                    3 * Qd * Vr * Vr * Vi - Pd * Vr * Vr * Vr) /
                   Vm6;

      d2Ird_dVi_dPdloss = 2 * Vr * Vi / Vm4;
      d2Ird_dVi_dQdloss = (Vi * Vi - Vr * Vr) / Vm4;

      d2Iid_dVr2 = -2 *
                   (Qd * Vr * Vr * Vr - 3 * Pd * Vi * Vr * Vr -
                    3 * Qd * Vi * Vi * Vr + Pd * Vi * Vi * Vi) /
                   Vm6;
      d2Iid_dVr_dVi = 2 *
                      (Qd * Vi * Vi * Vi + 3 * Pd * Vr * Vi * Vi -
                       3 * Qd * Vr * Vr * Vi - Pd * Vr * Vr * Vr) /
                      Vm6;

      d2Iid_dVr_dPdloss = 2 * Vr * Vi / Vm4;
      d2Iid_dVr_dQdloss = (Vi * Vi - Vr * Vr) / Vm4;

      d2Iid_dVi_dVr = -2 *
                      (Pd * Vr * Vr * Vr + 3 * Qd * Vr * Vr * Vi -
                       3 * Pd * Vr * Vi * Vi - Qd * Vi * Vi * Vi) /
                      Vm6;
      d2Iid_dVi2 = 2 *
                   (Pd * Vi * Vi * Vi - 3 * Qd * Vr * Vi * Vi -
                    3 * Pd * Vr * Vr * Vi + Qd * Vr * Vr * Vr) /
                   Vm6;

      d2Iid_dVi_dPdloss = (Vi * Vi - Vr * Vr) / Vm4;
      d2Iid_dVi_dQdloss = -2 * Vr * Vi / Vm4;

      row[0] = xlocglob;
      row[1] = xlocglob + 1;
      col[0] = xlocglob;
      col[1] = xlocglob + 1;

      val[0] = val[1] = val[2] = val[3] = 0.0;
      val[0] = lambda[gloc] * d2Ird_dVr2 + lambda[gloc + 1] * d2Iid_dVr2;
      val[1] = lambda[gloc] * d2Ird_dVr_dVi + lambda[gloc + 1] * d2Iid_dVr_dVi;
      val[2] = lambda[gloc] * d2Ird_dVi_dVr + lambda[gloc + 1] * d2Iid_dVi_dVr;
      val[3] = lambda[gloc] * d2Ird_dVi2 + lambda[gloc + 1] * d2Iid_dVi2;
      flps += 192;

      ierr = MatSetValues(H, 2, row, 2, col, val, ADD_VALUES);
      CHKERRQ(ierr);

      if (opflow->include_loadloss_variables) {

        col[0] = xlocglob + 2 + 2 * bus->ngen + loadctr;
        col[1] = xlocglob + 2 + 2 * bus->ngen + loadctr + 1;

        val[0] = val[1] = val[2] = val[3] = 0.0;
        val[0] = lambda[gloc] * d2Ird_dVr_dPdloss +
                 lambda[gloc + 1] * d2Iid_dVr_dPdloss;
        val[1] = lambda[gloc] * d2Ird_dVr_dQdloss +
                 lambda[gloc + 1] * d2Iid_dVr_dQdloss;
        val[2] = lambda[gloc] * d2Ird_dVi_dPdloss +
                 lambda[gloc + 1] * d2Iid_dVi_dPdloss;
        val[3] = lambda[gloc] * d2Ird_dVi_dQdloss +
                 lambda[gloc + 1] * d2Iid_dVi_dQdloss;

        ierr = MatSetValues(H, 2, row, 2, col, val, ADD_VALUES);
        CHKERRQ(ierr);

        PetscScalar d2Ird_dPdloss_dVr, d2Ird_dPdloss_dVi, d2Ird_dQdloss_dVr,
            d2Ird_dQdloss_dVi;
        PetscScalar d2Iid_dPdloss_dVr, d2Iid_dPdloss_dVi, d2Iid_dQdloss_dVr,
            d2Iid_dQdloss_dVi;

        row[0] = xlocglob + 2 + genctr;
        row[1] = xlocglob + 2 + genctr + 1;
        col[0] = xlocglob;
        col[1] = xlocglob + 1;

        d2Ird_dPdloss_dVr = (Vr * Vr - Vi * Vi) / Vm4;
        d2Ird_dPdloss_dVi = 2 * Vr * Vi / Vm4;
        d2Ird_dQdloss_dVr = 2 * Vr * Vi / Vm4;
        d2Ird_dQdloss_dVi = (Vi * Vi - Vr * Vr) / Vm4;

        d2Iid_dPdloss_dVr = 2 * Vr * Vi / Vm4;
        d2Iid_dPdloss_dVi = (Vi * Vi - Vr * Vr) / Vm4;

        d2Iid_dQdloss_dVr = (Vi * Vi - Vr * Vr) / Vm4;
        d2Iid_dQdloss_dVi = -2 * Vr * Vi / Vm4;

        val[0] = val[1] = val[2] = val[3] = 0.0;

        val[0] = lambda[gloc] * d2Ird_dPdloss_dVr +
                 lambda[gloc + 1] * d2Iid_dPdloss_dVr;
        val[1] = lambda[gloc] * d2Ird_dPdloss_dVi +
                 lambda[gloc + 1] * d2Iid_dPdloss_dVi;
        val[2] = lambda[gloc] * d2Ird_dQdloss_dVr +
                 lambda[gloc + 1] * d2Iid_dQdloss_dVr;
        val[3] = lambda[gloc] * d2Ird_dQdloss_dVi +
                 lambda[gloc + 1] * d2Iid_dQdloss_dVi;
        flps += 52;

        ierr = MatSetValues(H, 2, row, 2, col, val, ADD_VALUES);
        CHKERRQ(ierr);
      }
      loadctr += 2;
    }

    if (opflow->include_powerimbalance_variables) {
      PetscScalar Pimb, Qimb;
      PetscScalar d2Irimb_dVr2, d2Irimb_dVr_dVi, d2Irimb_dVr_dPimb,
          d2Irimb_dVr_dQimb;
      PetscScalar d2Irimb_dVi_dVr, d2Irimb_dVi2, d2Irimb_dVi_dPimb,
          d2Irimb_dVi_dQimb;
      PetscScalar d2Iiimb_dVr2, d2Iiimb_dVr_dVi, d2Iiimb_dVr_dPimb,
          d2Iiimb_dVr_dQimb;
      PetscScalar d2Iiimb_dVi_dVr, d2Iiimb_dVi2, d2Iiimb_dVi_dPimb,
          d2Iiimb_dVi_dQimb;

      xloc += 2;
      Pimb = x[xloc];
      Qimb = x[xloc + 1];

      d2Irimb_dVr2 = 2 *
                     (Pimb * Vr * Vr * Vr + 3 * Qimb * Vr * Vr * Vi -
                      3 * Pimb * Vr * Vi * Vi - Qimb * Vi * Vi * Vi) /
                     Vm6;
      d2Irimb_dVr_dVi = -2 *
                        (Pimb * Vi * Vi * Vi - 3 * Qimb * Vr * Vi * Vi -
                         3 * Pimb * Vr * Vr * Vi + Qimb * Vr * Vr * Vr) /
                        Vm6;
      d2Irimb_dVr_dPimb = (Vi * Vi - Vr * Vr) / Vm4;
      d2Irimb_dVr_dQimb = -2 * Vr * Vi / Vm4;

      d2Irimb_dVi_dVr = -2 *
                        (Qimb * Vr * Vr * Vr - 3 * Pimb * Vi * Vr * Vr -
                         3 * Qimb * Vi * Vi * Vr + Pimb * Vi * Vi * Vi) /
                        Vm6;
      d2Irimb_dVi2 = 2 *
                     (Qimb * Vi * Vi * Vi + 3 * Pimb * Vr * Vi * Vi -
                      3 * Qimb * Vr * Vr * Vi - Pimb * Vr * Vr * Vr) /
                     Vm6;
      d2Irimb_dVi_dPimb = -2 * Vr * Vi / Vm4;
      d2Irimb_dVi_dQimb = (Vr * Vr - Vi * Vi) / Vm4;

      d2Iiimb_dVr2 = -2 *
                     (Qimb * Vr * Vr * Vr - 3 * Pimb * Vi * Vr * Vr -
                      3 * Qimb * Vi * Vi * Vr + Pimb * Vi * Vi * Vi) /
                     Vm6;
      d2Iiimb_dVr_dVi = 2 *
                        (Qimb * Vi * Vi * Vi + 3 * Pimb * Vr * Vi * Vi -
                         3 * Qimb * Vr * Vr * Vi - Pimb * Vr * Vr * Vr) /
                        Vm6;
      d2Iiimb_dVr_dPimb = -2 * Vr * Vi / Vm4;
      d2Iiimb_dVr_dQimb = (Vr * Vr - Vi * Vi) / Vm4;

      d2Iiimb_dVi_dVr = -2 *
                        (Pimb * Vr * Vr * Vr + 3 * Qimb * Vr * Vr * Vi -
                         3 * Pimb * Vr * Vi * Vi - Qimb * Vi * Vi * Vi) /
                        Vm6;
      d2Iiimb_dVi2 = 2 *
                     (Pimb * Vi * Vi * Vi - 3 * Qimb * Vr * Vi * Vi -
                      3 * Pimb * Vr * Vr * Vi + Qimb * Vr * Vr * Vr) /
                     Vm6;
      d2Iiimb_dVi_dPimb = (Vr * Vr - Vi * Vi) / Vm4;
      d2Iiimb_dVi_dQimb = 2 * Vr * Vi / Vm4;

      row[0] = xlocglob;
      row[1] = xlocglob + 1;
      col[0] = xlocglob;
      col[1] = xlocglob + 1;
      col[2] = xlocglob + 2 + 2 * bus->ngen +
               2 * opflow->include_loadloss_variables * bus->nload;
      col[3] = xlocglob + 2 + 2 * bus->ngen +
               2 * opflow->include_loadloss_variables * bus->nload + 1;

      val[0] = val[1] = val[2] = val[3] = val[4] = val[5] = val[6] = val[7] =
          0.0;

      val[0] = lambda[gloc] * d2Irimb_dVr2 + lambda[gloc + 1] * d2Iiimb_dVr2;
      val[1] =
          lambda[gloc] * d2Irimb_dVr_dVi + lambda[gloc + 1] * d2Iiimb_dVr_dVi;
      val[2] = lambda[gloc] * d2Irimb_dVr_dPimb +
               lambda[gloc + 1] * d2Iiimb_dVr_dPimb;
      val[3] = lambda[gloc] * d2Irimb_dVr_dQimb +
               lambda[gloc + 1] * d2Iiimb_dVr_dQimb;

      val[4] =
          lambda[gloc] * d2Irimb_dVi_dVr + lambda[gloc + 1] * d2Iiimb_dVi_dVr;
      val[5] = lambda[gloc] * d2Irimb_dVi2 + lambda[gloc + 1] * d2Iiimb_dVi2;
      val[6] = lambda[gloc] * d2Irimb_dVi_dPimb +
               lambda[gloc + 1] * d2Iiimb_dVi_dPimb;
      val[7] = lambda[gloc] * d2Irimb_dVi_dQimb +
               lambda[gloc + 1] * d2Iiimb_dVi_dQimb;

      ierr = MatSetValues(H, 2, row, 4, col, val, ADD_VALUES);
      CHKERRQ(ierr);

      PetscScalar d2Irimb_dPimb_dVr, d2Irimb_dPimb_dVi, d2Irimb_dQimb_dVr,
          d2Irimb_dQimb_dVi;
      PetscScalar d2Iiimb_dPimb_dVr, d2Iiimb_dPimb_dVi, d2Iiimb_dQimb_dVr,
          d2Iiimb_dQimb_dVi;

      row[0] = xlocglob + 2 + 2 * bus->ngen +
               2 * opflow->include_loadloss_variables * bus->nload;
      row[1] = xlocglob + 2 + 2 * bus->ngen +
               2 * opflow->include_loadloss_variables * bus->nload + 1;

      col[0] = xlocglob;
      col[1] = xlocglob + 1;

      d2Irimb_dPimb_dVr = (Vi * Vi - Vr * Vr) / Vm4;
      d2Irimb_dPimb_dVi = -2 * Vr * Vi / Vm4;
      d2Irimb_dQimb_dVr = -2 * Vr * Vi / Vm4;
      d2Irimb_dQimb_dVi = (Vr * Vr - Vi * Vi) / Vm4;

      d2Iiimb_dPimb_dVr = -2 * Vr * Vi / Vm4;
      d2Iiimb_dPimb_dVi = (Vr * Vr - Vi * Vi) / Vm4;

      d2Iiimb_dQimb_dVr = (Vr * Vr - Vi * Vi) / Vm4;
      d2Iiimb_dQimb_dVi = 2 * Vr * Vi / Vm4;

      val[0] = val[1] = val[2] = val[3] = 0.0;

      val[0] = lambda[gloc] * d2Irimb_dPimb_dVr +
               lambda[gloc + 1] * d2Iiimb_dPimb_dVr;
      val[1] = lambda[gloc] * d2Irimb_dPimb_dVi +
               lambda[gloc + 1] * d2Iiimb_dPimb_dVi;
      val[2] = lambda[gloc] * d2Irimb_dQimb_dVr +
               lambda[gloc + 1] * d2Iiimb_dQimb_dVr;
      val[3] = lambda[gloc] * d2Irimb_dQimb_dVi +
               lambda[gloc + 1] * d2Iiimb_dQimb_dVi;
      flps += 246;

      ierr = MatSetValues(H, 2, row, 2, col, val, ADD_VALUES);
      CHKERRQ(ierr);
    }

    gloc += 2;

    if (bus->ide == REF_BUS)
      gloc++;
  }

  ierr = VecRestoreArrayRead(X, &x);
  CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(Lambda, &lambda);
  CHKERRQ(ierr);
  ierr = PetscLogFlops(flps);
  CHKERRQ(ierr);

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
PetscErrorCode OPFLOWComputeInequalityConstraintsHessian_IBCAR2(OPFLOW opflow,
                                                                Vec X,
                                                                Vec Lambda,
                                                                Mat H) {
  PetscErrorCode ierr;
  PS ps = opflow->ps;
  PetscInt i;
  PSLINE line;
  PetscInt xloc, xlocglob;
  PSBUS bus;
  const PetscScalar *x;
  const PetscScalar *lambda;
  PetscInt gloc = 0;
  PetscInt row[12], col[12];
  PetscScalar val[12];
  PetscScalar dVm2_dVr2, dVm2_dVi2;
  PetscInt flps = 0;

  PetscFunctionBegin;

  ierr = VecGetArrayRead(X, &x);
  CHKERRQ(ierr);
  ierr = VecGetArrayRead(Lambda, &lambda);
  CHKERRQ(ierr);

  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSGetVariableLocation(bus, &xloc);
    CHKERRQ(ierr);
    ierr = PSBUSGetVariableGlobalLocation(bus, &xlocglob);
    CHKERRQ(ierr);

    dVm2_dVr2 = 2.0;
    dVm2_dVi2 = 2.0;

    row[0] = xlocglob;
    row[1] = xlocglob + 1;
    col[0] = xlocglob;
    col[1] = xlocglob + 1;
    val[0] = lambda[gloc] * dVm2_dVr2;
    val[1] = val[2] = 0.0;
    val[3] = lambda[gloc] * dVm2_dVi2;

    ierr = MatSetValues(H, 2, row, 2, col, val, ADD_VALUES);
    CHKERRQ(ierr);
    gloc++;
  }
  flps += 2 * ps->nbus;

  if (opflow->ignore_lineflow_constraints) {
    ierr = VecRestoreArrayRead(X, &x);
    CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(Lambda, &lambda);
    CHKERRQ(ierr);

    PetscFunctionReturn(0);
  }

  val[0] = val[1] = val[2] = val[3] = 0.0;
  // for the part of line constraints
  for (i = 0; i < ps->nline; i++) {
    line = &ps->line[i];
    if (!line->status || line->rateA > 1e5)
      continue;

    ierr = PSLINEGetVariableGlobalLocation(line, &xlocglob);
    CHKERRQ(ierr);

    PetscScalar dIf2_dIrf_dIrf, dIf2_dIif_dIif;
    PetscScalar dIt2_dIrt_dIrt, dIt2_dIit_dIit;

    row[0] = xlocglob;
    row[1] = xlocglob + 1;
    col[0] = xlocglob;
    col[1] = xlocglob + 1;
    dIf2_dIrf_dIrf = 2.0;
    dIf2_dIif_dIif = 2.0;
    dIt2_dIrt_dIrt = 2.0;
    dIt2_dIit_dIit = 2.0;

    val[0] = lambda[gloc] * dIf2_dIrf_dIrf;
    val[1] = val[2] = 0.0;
    val[3] = lambda[gloc] * dIf2_dIif_dIif;

    ierr = MatSetValues(H, 2, row, 2, col, val, ADD_VALUES);
    CHKERRQ(ierr);

    gloc++;

    row[0] = xlocglob + 2;
    row[1] = xlocglob + 3;
    col[0] = xlocglob + 2;
    col[1] = xlocglob + 3;

    val[0] = lambda[gloc] * dIt2_dIrt_dIrt;
    val[1] = val[2] = 0.0;
    val[3] = lambda[gloc] * dIt2_dIit_dIit;

    ierr = MatSetValues(H, 2, row, 2, col, val, ADD_VALUES);
    CHKERRQ(ierr);

    gloc++;
  }
  flps += 4 * ps->nline;

  ierr = VecRestoreArrayRead(X, &x);
  CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(Lambda, &lambda);
  CHKERRQ(ierr);
  ierr = PetscLogFlops(flps);
  CHKERRQ(ierr);

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
PetscErrorCode OPFLOWComputeObjectiveHessian_IBCAR2(OPFLOW opflow, Vec X,
                                                    Mat H) {
  PetscErrorCode ierr;
  PS ps = opflow->ps;
  PetscInt i, k;
  PSBUS bus;
  PSGEN gen;
  PetscInt xloc, xlocglob;
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

    ierr = PSBUSGetVariableLocation(bus, &xloc);
    CHKERRQ(ierr);
    ierr = PSBUSGetVariableGlobalLocation(bus, &xlocglob);
    CHKERRQ(ierr);

    for (k = 0; k < bus->ngen; k++) {
      ierr = PSBUSGetGen(bus, k, &gen);
      CHKERRQ(ierr);
      if (!gen->status)
        continue;
      xlocglob = xlocglob + 2;
      row[0] = xlocglob;
      col[0] = xlocglob;
      if (opflow->obj_gencost) {
        val[0] = obj_factor * 2.0 * gen->cost_alpha * ps->MVAbase * ps->MVAbase;
        ierr = MatSetValues(H, 1, row, 1, col, val, ADD_VALUES);
        CHKERRQ(ierr);
        flps += 4.0;
      }
    }

    if (opflow->include_loadloss_variables) {
      PSLOAD load;
      for (k = 0; k < bus->nload; k++) {
        xlocglob = xlocglob + 2;
        ierr = PSBUSGetLoad(bus, k, &load);
        CHKERRQ(ierr);
        if (!load->status)
          continue;
        row[0] = xlocglob;
        col[0] = xlocglob;
        val[0] = obj_factor * 2.0 * opflow->loadloss_penalty * ps->MVAbase *
                 ps->MVAbase;
        ierr = MatSetValues(H, 1, row, 1, col, val, ADD_VALUES);
        CHKERRQ(ierr);

        row[0] = xlocglob + 1;
        col[0] = xlocglob + 1;
        val[0] = obj_factor * 2.0 * opflow->loadloss_penalty * ps->MVAbase *
                 ps->MVAbase;
        ierr = MatSetValues(H, 1, row, 1, col, val, ADD_VALUES);
        CHKERRQ(ierr);
      }
      flps += 8 * bus->nload;
    }

    if (opflow->include_powerimbalance_variables) {
      xlocglob = xlocglob + 2;
      row[0] = xlocglob;
      col[0] = xlocglob;
      val[0] = obj_factor * 2.0 * opflow->powerimbalance_penalty * ps->MVAbase *
               ps->MVAbase;
      ierr = MatSetValues(H, 1, row, 1, col, val, ADD_VALUES);
      CHKERRQ(ierr);

      row[0] = xlocglob + 1;
      col[0] = xlocglob + 1;
      val[0] = obj_factor * 2.0 * opflow->powerimbalance_penalty * ps->MVAbase *
               ps->MVAbase;
      flps += 8;
      ierr = MatSetValues(H, 1, row, 1, col, val, ADD_VALUES);
      CHKERRQ(ierr);
    }
  }

  ierr = VecRestoreArrayRead(X, &x);
  CHKERRQ(ierr);
  ierr = PetscLogFlops(flps);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeHessian_IBCAR2(OPFLOW opflow, Vec X, Vec Lambdae,
                                           Vec Lambdai, Mat H) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = MatZeroEntries(H);
  CHKERRQ(ierr);

  /* Objective function Hessian */
  ierr = OPFLOWComputeObjectiveHessian_IBCAR2(opflow, X, H);
  CHKERRQ(ierr);

  /* Equality constraints Hessian */
  ierr = OPFLOWComputeEqualityConstraintsHessian_IBCAR2(opflow, X, Lambdae, H);
  CHKERRQ(ierr);

  /* Inequality constraints Hessian */
  if (opflow->nconineq) {
    ierr =
        OPFLOWComputeInequalityConstraintsHessian_IBCAR2(opflow, X, Lambdai, H);
    CHKERRQ(ierr);
  }

  ierr = MatAssemblyBegin(H, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = MatAssemblyEnd(H, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWModelCreate_IBCAR2(OPFLOW opflow) {
  IBCAR ibcar2;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = PetscCalloc1(1, &ibcar2);
  CHKERRQ(ierr);

  opflow->model = ibcar2;

  /* Inherit Ops */
  opflow->modelops.destroy = OPFLOWModelDestroy_IBCAR2;
  opflow->modelops.setnumvariables = OPFLOWModelSetNumVariables_IBCAR2;
  opflow->modelops.setnumconstraints = OPFLOWModelSetNumConstraints_IBCAR2;
  opflow->modelops.setvariablebounds = OPFLOWSetVariableBounds_IBCAR2;
  opflow->modelops.setconstraintbounds = OPFLOWSetConstraintBounds_IBCAR2;
  opflow->modelops.setvariableandconstraintbounds =
      OPFLOWSetVariableandConstraintBounds_IBCAR2;
  opflow->modelops.setinitialguess = OPFLOWSetInitialGuess_IBCAR2;
  opflow->modelops.computeequalityconstraints =
      OPFLOWComputeEqualityConstraints_IBCAR2;
  opflow->modelops.computeinequalityconstraints =
      OPFLOWComputeInequalityConstraints_IBCAR2;
  opflow->modelops.computeequalityconstraintjacobian =
      OPFLOWComputeEqualityConstraintJacobian_IBCAR2;
  opflow->modelops.computeinequalityconstraintjacobian =
      OPFLOWComputeInequalityConstraintJacobian_IBCAR2;
  opflow->modelops.computehessian = OPFLOWComputeHessian_IBCAR2;
  opflow->modelops.computeobjandgradient = OPFLOWComputeObjandGradient_IBCAR2;
  opflow->modelops.computeobjective = OPFLOWComputeObjective_IBCAR2;
  opflow->modelops.computegradient = OPFLOWComputeGradient_IBCAR2;

  PetscFunctionReturn(0);
}
