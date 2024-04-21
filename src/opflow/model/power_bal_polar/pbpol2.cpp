 #include "pbpol.h"
#include "exago_config.h"
#include <private/opflowimpl.h>

PetscErrorCode OPFLOWModelDestroy_PBPOL2(OPFLOW opflow) {
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = PetscFree(opflow->model);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSetVariableBounds_PBPOL2(OPFLOW opflow, Vec Xl, Vec Xu) {
  PetscErrorCode ierr;
  PS ps = opflow->ps;
  PetscScalar *xl, *xu;
  PetscInt i;
  PSBUS bus;
  PSGEN gen;
  PSLINE line;
  PetscInt loc;

  PetscFunctionBegin;

  /* Get array pointers */
  ierr = VecGetArray(Xl, &xl);
  CHKERRQ(ierr);
  ierr = VecGetArray(Xu, &xu);
  CHKERRQ(ierr);

  for (i = 0; i < ps->nline; i++) {
    line = &ps->line[i];
    if(!line->status) continue;
    loc = line->startxflowloc;
    
    xl[loc] = xl[loc + 1] = xl[loc + 2] = xl[loc + 3] = PETSC_NINFINITY;
    xu[loc] = xu[loc + 1] = xu[loc + 2] = xu[loc + 3] = PETSC_INFINITY;

    if (!opflow->ignore_lineflow_constraints && opflow->allow_lineflow_violation) {
      /* Bounds on slack variables */
      loc = line->startxslackloc;
      xl[loc] = xl[loc + 1] = 0.0;
      xu[loc] = xu[loc + 1] = PETSC_INFINITY;
    }
  }

  for (i = 0; i < ps->nbus; i++) {
    PetscInt k;

    bus = &ps->bus[i];

    /* Bounds on bus variables */
    loc = bus->startxVloc;

    xl[loc] = PETSC_NINFINITY;
    xu[loc] = PETSC_INFINITY;

    if (opflow->genbusvoltagetype == VARIABLE_WITHIN_BOUNDS) {
      xl[loc + 1] = bus->Vmin;
      xu[loc + 1] = bus->Vmax;
    } else if (opflow->genbusvoltagetype == FIXED_WITHIN_QBOUNDS) {
      if (bus->ide == REF_BUS || bus->ide == PV_BUS) {
        xl[loc + 1] = 0.0;
        xu[loc + 1] = 2.0;
      } else {
        xl[loc + 1] = bus->Vmin;
        xu[loc + 1] = bus->Vmax; /* PQ buses */
      }
    } else if (opflow->genbusvoltagetype == FIXED_AT_SETPOINT) {
      if (bus->ide == REF_BUS || bus->ide == PV_BUS) {
        xl[loc + 1] = bus->vm;
        xu[loc + 1] = bus->vm;
      } else {
        xl[loc + 1] = bus->Vmin;
        xu[loc + 1] = bus->Vmax; /* PQ buses */
      }
    }

    if (bus->ide == REF_BUS || bus->ide == ISOLATED_BUS)
      xl[loc] = xu[loc] = bus->va;
    if (bus->ide == ISOLATED_BUS)
      xl[loc + 1] = xu[loc + 1] = bus->vm;

    if (opflow->include_powerimbalance_variables) {
      loc = bus->startxpimbloc;

      xl[loc] = xl[loc + 1] = xl[loc + 2] = xl[loc + 3] = 0.0;
      xu[loc] = xu[loc + 1] = xu[loc + 2] = xu[loc + 3] = PETSC_INFINITY;

      loc += bus->nxpimb;
    }

    /* Bounds on generator variables */

    for (k = 0; k < bus->ngen; k++) {
      ierr = PSBUSGetGen(bus, k, &gen);
      CHKERRQ(ierr);
      if (!gen->status)
        continue;

      loc = gen->startxpowloc;

      xl[loc] = gen->pb;     /* PGmin */
      xu[loc] = gen->pt;     /* PGmax */
      xl[loc + 1] = gen->qb; /* QGmin */
      xu[loc + 1] = gen->qt; /* QGmax */
      /* pb, pt, qb, qt are converted in p.u. in ps.c */

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
          xl[loc] = xu[loc] = xl[loc + 1] = xu[loc + 1] = 0.0;
        else {
          xl[loc] = PetscMin(0.0, load->pl * load->loss_frac);
          xu[loc] = PetscMax(0.0, load->pl * load->loss_frac);
          xl[loc + 1] = PetscMin(0.0, load->ql * load->loss_frac);
          xu[loc + 1] = PetscMax(0.0, load->ql * load->loss_frac);
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

PetscErrorCode OPFLOWSetConstraintBounds_PBPOL2(OPFLOW opflow, Vec Gl, Vec Gu) {
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
    /* Equality constraint bounds for real and reactive power mismatch at the
     * bus */
    gl[gloc] = 0.0;
    gu[gloc] = 0.0;
    gl[gloc + 1] = 0.0;
    gu[gloc + 1] = 0.0;

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

  /* Line flow equality constraints */
  for (i = 0; i < ps->nline; i++) {
    line = &ps->line[i];
    
    if(!line->status) continue;
    /* Equality constraints */
    gloc = line->starteqloc;
    gl[gloc] = gl[gloc + 1] = gl[gloc + 2] = gl[gloc + 3] = 0.0;
    gu[gloc] = gu[gloc + 1] = gu[gloc + 2] = gu[gloc + 3] = 0.0;
  }

  /** INEQUALITY CONSTRAINTS **/
  gl += opflow->nconeq;
  gu += opflow->nconeq;

  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    if (opflow->genbusvoltagetype == FIXED_WITHIN_QBOUNDS) {
      /* Inequality constraint bounds for voltages */
      if (bus->ide == PV_BUS || bus->ide == REF_BUS) {
        gloc = bus->startineqloc;
        gl[gloc] = gl[gloc + 1] = 0;
        gu[gloc] = gu[gloc + 1] = PETSC_INFINITY;
      }
    }
  }

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
  for (i = 0; i < ps->nline; i++) {
    line = &ps->line[i];

    if(!opflow->ignore_lineflow_constraints && line->monitor_flow) {
      gloc = line->startineqloc;
      /* Line flow inequality constraints */
      gl[gloc] = gl[gloc + 1] = 0.0;
      gu[gloc] = gu[gloc + 1] =
        (line->rateA / ps->MVAbase) * (line->rateA / ps->MVAbase);
    }
  }

  ierr = VecRestoreArray(Gl, &gl);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(Gu, &gu);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSetVariableandConstraintBounds_PBPOL2(OPFLOW opflow, Vec Xl,
                                                          Vec Xu, Vec Gl,
                                                          Vec Gu) {
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = OPFLOWSetVariableBounds_PBPOL2(opflow, Xl, Xu);
  CHKERRQ(ierr);
  ierr = OPFLOWSetConstraintBounds_PBPOL2(opflow, Gl, Gu);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSetInitialGuess_PBPOL2(OPFLOW opflow, Vec X, Vec Lambda) {
  PetscErrorCode ierr;
  PS ps = opflow->ps;
  const PetscScalar *xl, *xu;
  PetscScalar *x, *lambda;
  PetscInt i;
  PSBUS bus;
  PSLINE line;
  PetscInt loc, gloc;

  PetscFunctionBegin;

  /* Get array pointers */
  ierr = VecGetArray(X, &x);
  CHKERRQ(ierr);
  ierr = VecGetArray(Lambda, &lambda);
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
      x[loc + 1] = bus->vm;
    } else {
      if (opflow->initializationtype == OPFLOWINIT_MIDPOINT) {
        /* Initial guess for voltage angles and bounds on voltage magnitudes */
        x[loc] = (xl[loc] + xu[loc]) / 2.0;
        x[loc + 1] = (xl[loc + 1] + xu[loc + 1]) / 2.0;
      } else if (opflow->initializationtype == OPFLOWINIT_FROMFILE ||
                 opflow->initializationtype == OPFLOWINIT_ACPF ||
                 opflow->initializationtype == OPFLOWINIT_DCOPF) {
        x[loc] = bus->va;
        x[loc + 1] = PetscMax(bus->Vmin, PetscMin(bus->vm, bus->Vmax));
      } else if (opflow->initializationtype == OPFLOWINIT_FLATSTART) {
        x[loc] = 0.0;
        x[loc + 1] = 1.0;
      }
    }

    if (opflow->include_powerimbalance_variables) {
      loc = bus->startxpimbloc;
      x[loc] = x[loc + 1] = x[loc + 2] = x[loc + 3] = 0.0;
      if (bus->pimb > 0) {
        x[loc] = bus->pimb;
        x[loc + 1] = 0.0;
      } else {
        x[loc] = 0.0;
        x[loc + 1] = -bus->pimb;
      }
      if (bus->qimb > 0) {
        x[loc + 2] = bus->qimb;
        x[loc + 3] = 0.0;
      } else {
        x[loc + 2] = 0.0;
        x[loc + 3] = -bus->qimb;
      }
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
        x[loc + 1] = 0.5 * (xl[loc + 1] + xu[loc + 1]);
      } else if (opflow->initializationtype == OPFLOWINIT_FROMFILE ||
                 opflow->initializationtype == OPFLOWINIT_ACPF ||
                 opflow->initializationtype == OPFLOWINIT_DCOPF) {
        x[loc] = PetscMax(gen->pb, PetscMin(gen->pg, gen->pt));
        x[loc + 1] = PetscMax(gen->qb, PetscMin(gen->qg, gen->qt));
      }
      loc += gen->nxpow;

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
        /* Initial value for real and reactive power load loss */
        x[loc] = load->pl_loss;
        x[loc + 1] = load->ql_loss;
      }
    }
    if (opflow->use_agc) {
      if (bus->ide == REF_BUS) {
        loc = ps->startxloc;
        x[loc] = 0.0;
      }
    }

    gloc = bus->starteqloc;
    lambda[gloc] = bus->mult_pmis;
    lambda[gloc + 1] = bus->mult_qmis;
  }

  PetscScalar *lambdai = lambda + opflow->nconeq;

  for (i = 0; i < ps->nline; i++) {
    line = &ps->line[i];
    if(!line->status) continue;

    loc = line->startxflowloc;

    if (opflow->initializationtype == OPFLOWINIT_MIDPOINT ||
	opflow->initializationtype == OPFLOWINIT_FLATSTART) {
        x[loc] = 0.5 * (xl[loc] + xu[loc]);
        x[loc + 1] = 0.5 * (xl[loc + 1] + xu[loc + 1]);
	x[loc + 2] = 0.5 * (xl[loc + 2] + xu[loc + 2]);
        x[loc + 3] = 0.5 * (xl[loc + 3] + xu[loc + 3]);
    } else if (opflow->initializationtype == OPFLOWINIT_FROMFILE ||
                 opflow->initializationtype == OPFLOWINIT_ACPF ||
                 opflow->initializationtype == OPFLOWINIT_DCOPF) {
      x[loc] = line->pf;     x[loc + 1] = line->qf;
      x[loc + 2] = line->pt; x[loc + 3] = line->qt;
    }
    if (!opflow->ignore_lineflow_constraints && opflow->allow_lineflow_violation && line->monitor_flow) {
      loc = line->startxslackloc;
      /* Initialize slacks for line flow violations */
      x[loc] = x[loc + 1] = 0.0;
    }

    if(!opflow->ignore_lineflow_constraints && line->monitor_flow) {
      gloc = line->startineqloc;
      lambdai[gloc] = line->mult_sf;
      lambdai[gloc + 1] = line->mult_st;
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

PetscErrorCode OPFLOWComputeEqualityConstraints_PBPOL2(OPFLOW opflow, Vec X,
                                                      Vec Ge) {
  PetscErrorCode ierr;
  PetscInt i, k, nconnlines;
  PetscInt gloc, row[4];
  PetscInt xloc, xlocf, xloct;
  PetscScalar val[4];
  PetscScalar Pg, Qg, Pd, Qd;
  PetscScalar Gff, Bff, Gft, Bft, Gtf, Btf, Gtt, Btt;
  PetscScalar Vmf, Vmt, thetaf, thetat, thetaft, thetatf;
  PetscScalar Pf, Qf, Pt, Qt;
  PetscScalar theta, Vm;
  PS ps = opflow->ps;
  PSLOAD load;
  PSLINE line;
  PSBUS bus, busf, bust;
  PSGEN gen;
  const PSBUS *connbuses;
  const PSLINE *connlines;
  const PetscScalar *x;
  double flps = 0.0;

  PetscFunctionBegin;
  ierr = VecSet(Ge, 0.0);
  CHKERRQ(ierr);

  ierr = VecGetArrayRead(X, &x);
  CHKERRQ(ierr);

  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    gloc = bus->starteqloc;

    row[0] = gloc;
    row[1] = row[0] + 1;

    xloc = bus->startxVloc;

    theta = x[xloc];
    Vm = x[xloc + 1];

    if (bus->ide == ISOLATED_BUS) {
      row[0] = gloc;
      row[1] = row[0] + 1;
      val[0] = theta - bus->va;
      val[1] = Vm - bus->vm;
      ierr = VecSetValues(Ge, 2, row, val, ADD_VALUES);
      CHKERRQ(ierr);
      continue;
    }

    /* Shunt injections */
    val[0] = Vm * Vm * bus->gl;
    val[1] = -Vm * Vm * bus->bl;
    ierr = VecSetValues(Ge, 2, row, val, ADD_VALUES);
    CHKERRQ(ierr);

    flps += 5.0;

    /* Power imbalance addition */
    if (opflow->include_powerimbalance_variables) {
      PetscScalar Pimbplus, Pimbminus, Qimbplus, Qimbminus, Pimb, Qimb;
      xloc = bus->startxpimbloc;
      Pimbplus = x[xloc];
      Pimbminus = x[xloc + 1];
      Qimbplus = x[xloc + 2];
      Qimbminus = x[xloc + 3];

      Pimb = Pimbplus - Pimbminus;
      Qimb = Qimbplus - Qimbminus;
      val[0] = Pimb;
      val[1] = Qimb;
      ierr = VecSetValues(Ge, 2, row, val, ADD_VALUES);
      CHKERRQ(ierr);

      flps += 2.0;
    }

    for (k = 0; k < bus->ngen; k++) {
      ierr = PSBUSGetGen(bus, k, &gen);
      CHKERRQ(ierr);
      if (!gen->status)
        continue;

      xloc = gen->startxpowloc;

      Pg = x[xloc];
      Qg = x[xloc + 1];

      val[0] = -Pg;
      val[1] = -Qg;
      ierr = VecSetValues(Ge, 2, row, val, ADD_VALUES);
      CHKERRQ(ierr);

      flps += 2.0;
    }

    for (k = 0; k < bus->nload; k++) {
      ierr = PSBUSGetLoad(bus, k, &load);
      CHKERRQ(ierr);
      xloc = load->startxloadlossloc;
      if (opflow->include_loadloss_variables) {
        Pd = load->pl - x[xloc];
        Qd = load->ql - x[xloc + 1];
        xloc += load->nxloadloss;
      } else {
        Pd = load->pl;
        Qd = load->ql;
      }

      val[0] = Pd;
      val[1] = Qd;
      ierr = VecSetValues(Ge, 2, row, val, ADD_VALUES);
      CHKERRQ(ierr);
      flps += 2.0;
    }

    ierr = PSBUSGetSupportingLines(bus, &nconnlines, &connlines);
    CHKERRQ(ierr);
    for (k = 0; k < nconnlines; k++) {
      line = connlines[k];
      if (!line->status)
        continue;

      xloc = line->startxflowloc;

      Pf = x[line->startxflowloc];
      Qf = x[line->startxflowloc + 1];
      Pt = x[line->startxflowloc + 2];
      Qt = x[line->startxflowloc + 3];
      
      ierr = PSLINEGetConnectedBuses(line, &connbuses);
      CHKERRQ(ierr);
      busf = connbuses[0];
      bust = connbuses[1];

      if (bus == busf) {
        val[0] = Pf;
        val[1] = Qf;
        ierr = VecSetValues(Ge, 2, row, val, ADD_VALUES);
        CHKERRQ(ierr);
      } else {
        val[0] = Pt;
        val[1] = Qt;
        ierr = VecSetValues(Ge, 2, row, val, ADD_VALUES);
        CHKERRQ(ierr);
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

  for(i = 0; i < ps->nline; i++) {
    line = &ps->line[i];
    if(!line->status) continue;

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

    Pf = x[line->startxflowloc];
    Qf = x[line->startxflowloc + 1];
    Pt = x[line->startxflowloc + 2];
    Qt = x[line->startxflowloc + 3];

    /* From and to real and reactive flow balance */
    val[0] = Gff * Vmf * Vmf +
      Vmf * Vmt * (Gft * cos(thetaft) + Bft * sin(thetaft)) - Pf;
    val[1] = -Bff * Vmf * Vmf +
      Vmf * Vmt * (-Bft * cos(thetaft) + Gft * sin(thetaft)) - Qf;
    val[2] = Gtt * Vmt * Vmt +
      Vmt * Vmf * (Gtf * cos(thetatf) + Btf * sin(thetatf)) - Pt;
    val[3] = -Btt * Vmt * Vmt +
      Vmt * Vmf * (-Btf * cos(thetatf) + Gtf * sin(thetatf)) - Qt;

    row[0] = line->starteqloc;
    row[1] = row[0] + 1;
    row[2] = row[1] + 1;
    row[3] = row[2] + 1;

    ierr = VecSetValues(Ge, 4, row, val, ADD_VALUES);
    CHKERRQ(ierr);    
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

PetscErrorCode OPFLOWComputeEqualityConstraintJacobian_PBPOL2(OPFLOW opflow,
                                                             Vec X, Mat Je) {
  PetscErrorCode ierr;
  PetscInt i, k, row[2], col[4], gidx, flps = 0;
  PetscInt nconnlines, locglob, loc, locglobf, locglobt, locf, loct;
  PetscScalar Vm, val[8], Gff, Bff, Gft, Bft, Gtf, Btf, Gtt, Btt;
  PetscScalar thetaf, thetat, Vmf, Vmt, thetaft, thetatf;
  PS ps = opflow->ps;
  PSBUS bus;
  PSLINE line;
  PSBUS busf, bust;
  PSGEN gen;
  PSLOAD load;
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
    row[1] = row[0] + 1;

    loc = bus->startxVloc;
    locglob = bus->startxVlocglob;

    Vm = xarr[loc + 1];

    col[0] = locglob;
    col[1] = locglob + 1;
    /* Isolated and reference bus */
    if (bus->ide == ISOLATED_BUS) {
      val[0] = val[3] = 1.0;
      val[1] = val[2] = 0.0;
      ierr = MatSetValues(Je, 2, row, 2, col, val, ADD_VALUES);
      CHKERRQ(ierr);
      continue;
    }
    /* Shunt injections */
    val[0] = 0.0;
    val[1] = 2 * Vm * bus->gl;
    val[2] = 0.0;
    val[3] = -2 * Vm * bus->bl; /* Partial derivative for shunt contribution */
    ierr = MatSetValues(Je, 2, row, 2, col, val, ADD_VALUES);
    CHKERRQ(ierr);

    /* Power imbalance Jacobian terms */
    if (opflow->include_powerimbalance_variables) {
      locglob = bus->startxpimblocglob;

      val[0] = 1;
      val[1] = -1;
      col[0] = locglob;
      col[1] = locglob + 1;
      ierr = MatSetValues(Je, 1, row, 2, col, val, ADD_VALUES);
      CHKERRQ(ierr);

      col[0] = locglob + 2;
      col[1] = locglob + 3;
      ierr = MatSetValues(Je, 1, row + 1, 2, col, val, ADD_VALUES);
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

      col[0] = locglob + 1;
      ierr = MatSetValues(Je, 1, row + 1, 1, col, val, ADD_VALUES);
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

        col[0] = locglob + 1;
        ierr = MatSetValues(Je, 1, row + 1, 1, col, val, ADD_VALUES);
        CHKERRQ(ierr);
      }
    }

    /* Partial derivatives of network equations */
    /* Get the lines supporting the bus */
    ierr = PSBUSGetSupportingLines(bus, &nconnlines, &connlines);
    CHKERRQ(ierr);

    for (k = 0; k < nconnlines; k++) {
      line = connlines[k];
      if (!line->status) continue;

      /* Get the connected buses to this line */
      ierr = PSLINEGetConnectedBuses(line, &connbuses);
      CHKERRQ(ierr);
      busf = connbuses[0];
      bust = connbuses[1];

      if (bus == busf) {
        col[0] = line->startxflowloc;
        col[1] = line->startxflowloc + 1;

	val[0] = 1.0; val[1] = 1.0;
	
        ierr = MatSetValues(Je, 1, row, 1, col, val, ADD_VALUES);
	CHKERRQ(ierr);
	ierr = MatSetValues(Je, 1, row+1, 1, col+1, val+1, ADD_VALUES);
	CHKERRQ(ierr);
      } else {
	col[0] = line->startxflowloc + 2;
	col[1] = line->startxflowloc + 3;

	val[0] = 1.0; val[1] = 1.0;

	ierr = MatSetValues(Je, 1, row, 1, col, val, ADD_VALUES);
	CHKERRQ(ierr);
	ierr = MatSetValues(Je, 1, row+1, 1, col+1, val+1, ADD_VALUES);
	CHKERRQ(ierr);
      }
    }

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

  for (i = 0; i < ps->nline; i++) {
    line = &ps->line[i];
    if (!line->status) continue;
    
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

    locf = busf->startxVloc;
    loct = bust->startxVloc;

    locglobf = busf->startxVlocglob;
    locglobt = bust->startxVlocglob;

    thetaf = xarr[locf];
    Vmf = xarr[locf + 1];
    thetat = xarr[loct];
    Vmt = xarr[loct + 1];
    thetaft = thetaf - thetat;
    thetatf = thetat - thetaf;

    col[0] = locglobf;
    col[1] = locglobf + 1;
    col[2] = locglobt;
    col[3] = locglobt + 1;
    /* dPf_dthetaf */
    val[0] = Vmf * Vmt * (-Gft * sin(thetaft) + Bft * cos(thetaft));
    /*dPf_dVmf */
    val[1] =
      2 * Gff * Vmf + Vmt * (Gft * cos(thetaft) + Bft * sin(thetaft));
    /*dPf_dthetat */
    val[2] = Vmf * Vmt * (Gft * sin(thetaft) - Bft * cos(thetaft));
    /* dPf_dVmt */
    val[3] = Vmf * (Gft * cos(thetaft) + Bft * sin(thetaft));
    
    /* dQf_dthetaf */
    val[4] = Vmf * Vmt * (Bft * sin(thetaft) + Gft * cos(thetaft));
    /* dQf_dVmf */
    val[5] =
      -2 * Bff * Vmf + Vmt * (-Bft * cos(thetaft) + Gft * sin(thetaft));
    /* dQf_dthetat */
    val[6] = Vmf * Vmt * (-Bft * sin(thetaft) - Gft * cos(thetaft));
    /* dQf_dVmt */
    val[7] = Vmf * (-Bft * cos(thetaft) + Gft * sin(thetaft));

    row[0] = line->starteqloc;
    row[1] = line->starteqloc + 1;
    ierr = MatSetValues(Je, 2, row, 4, col, val, ADD_VALUES);
    CHKERRQ(ierr);

    col[0] = line->startxflowloc;
    val[0] = -1.0;
    ierr = MatSetValues(Je, 1, row, 1, col, val, ADD_VALUES);
    CHKERRQ(ierr);

    col[0] = line->startxflowloc + 1;
    val[0] = -1.0;
    ierr = MatSetValues(Je, 1, row+1, 1, col, val, ADD_VALUES);
    CHKERRQ(ierr);

    col[0] = locglobt;
    col[1] = locglobt + 1;
    col[2] = locglobf;
    col[3] = locglobf + 1;
    /* dPt_dthetat */
    val[0] = Vmt * Vmf * (-Gtf * sin(thetatf) + Btf * cos(thetatf));
    /* dPt_dVmt */
    val[1] =
      2 * Gtt * Vmt + Vmf * (Gtf * cos(thetatf) + Btf * sin(thetatf));
    /* dPt_dthetaf */
    val[2] = Vmt * Vmf * (Gtf * sin(thetatf) - Btf * cos(thetatf));
    /* dPt_dVmf */
    val[3] = Vmt * (Gtf * cos(thetatf) + Btf * sin(thetatf));
    
    /* dQt_dthetat */
    val[4] = Vmt * Vmf * (Btf * sin(thetatf) + Gtf * cos(thetatf));
    /* dQt_dVmt */
    val[5] =
      -2 * Btt * Vmt + Vmf * (-Btf * cos(thetatf) + Gtf * sin(thetatf));
    /* dQt_dthetaf */
    val[6] = Vmt * Vmf * (-Btf * sin(thetatf) - Gtf * cos(thetatf));
    /* dQt_dVmf */
    val[7] = Vmt * (-Btf * cos(thetatf) + Gtf * sin(thetatf));

    row[0] = line->starteqloc + 2;
    row[1] = line->starteqloc + 3;
    
    ierr = MatSetValues(Je, 2, row, 4, col, val, ADD_VALUES);
    CHKERRQ(ierr);
    
    col[0] = line->startxflowloc + 2;
    val[0] = -1.0;
    ierr = MatSetValues(Je, 1, row, 1, col, val, ADD_VALUES);
    CHKERRQ(ierr);

    col[0] = line->startxflowloc + 3;
    val[0] = -1.0;
    ierr = MatSetValues(Je, 1, row+1, 1, col, val, ADD_VALUES);
    CHKERRQ(ierr);

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

PetscErrorCode OPFLOWComputeInequalityConstraints_PBPOL2(OPFLOW opflow, Vec X,
                                                        Vec Gi) {
  PetscErrorCode ierr;
  PetscInt i, k;
  PetscInt gloc;
  PetscScalar *g;
  PetscScalar Pf, Qf, Pt, Qt, Sf2, St2;
  PS ps = opflow->ps;
  PSLINE line;
  PSBUS bus;
  PSGEN gen;
  const PetscScalar *x;
  double flps = 0.0;
  PetscScalar Pg, delP, delPg;
  PetscInt xloc;
  PetscScalar V, Vset, Q, Qmax, Qmin, Qg;

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

  if (opflow->genbusvoltagetype == FIXED_WITHIN_QBOUNDS) {
    for (i = 0; i < ps->nbus; i++) {
      bus = &ps->bus[i];
      if (bus->ide == PV_BUS || bus->ide == REF_BUS) {
        gloc = bus->startineqloc;
        xloc = bus->startxVloc;
        V = x[xloc + 1];

        /*to get the correct reserves and condition when PV should become PQ we
         * need to loop through all gen connected to a bus and sum reactive
         * reserves and limits*/
        Q = 0;
        Qmax = 0;
        Qmin = 0;
        for (k = 0; k < bus->ngen; k++) {
          ierr = PSBUSGetGen(bus, k, &gen);
          CHKERRQ(ierr);
          if (!gen->status)
            continue;
          xloc = gen->startxpowloc;
          Qg = x[xloc + 1];
          Q = Q + Qg;
          Qmax = Qmax + gen->qt;
          Qmin = Qmin + gen->qb;
          Vset = gen->vs;
        }

        /* Inequality constraints */
        g[gloc] = (Q - Qmax) * (Vset - V);
        g[gloc + 1] = (Qmin - Q) * (V - Vset);
      }
    }
  }

  if (!opflow->ignore_lineflow_constraints) {
    for (i = 0; i < ps->nline; i++) {
      line = &ps->line[i];
      if(!line->monitor_flow) continue;
     
      xloc = line->startxflowloc;

      Pf = x[xloc];     Qf = x[xloc + 1];
      Pt = x[xloc + 2]; Qt = x[xloc + 3];

      gloc = line->startineqloc;

      Sf2 = Pf * Pf + Qf * Qf;
      St2 = Pt * Pt + Qt * Qt;

      g[gloc] = Sf2;
      g[gloc + 1] = St2;

      if (opflow->allow_lineflow_violation) {
        PetscInt loc;
        PetscScalar xsft_slack = 0.0, xstf_slack = 0.0;
        loc = line->startxslackloc;
        // Slack variables for from and to side
        xsft_slack = x[loc];
        xstf_slack = x[loc + 1];

        g[gloc] -= xsft_slack;
        g[gloc + 1] -= xstf_slack;
      }

      flps += 160.0;
    }
  }

  ierr = VecRestoreArrayRead(X, &x);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(Gi, &g);
  CHKERRQ(ierr);

  PetscLogFlops(flps);
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeInequalityConstraintJacobian_PBPOL2(OPFLOW opflow,
                                                               Vec X, Mat Ji) {
  PetscErrorCode ierr;
  PetscInt i, k, flps = 0;
  PetscInt row[2], col[4];
  PetscInt rstart, rend;
  PetscInt loc, gloc;
  PetscScalar val[4];
  PetscScalar Pf, Qf, Pt, Qt;
  PetscScalar dSf2_dPf, dSf2_dQf, dSt2_dPt, dSt2_dQt;
  PS ps = opflow->ps;
  PSLINE line;
  PSBUS bus;
  PSGEN gen;
  const PetscScalar *x;
  PetscScalar Pg, delPg, delP;
  PetscInt xloc;
  PetscScalar V, Vset, Q, Qmax, Qmin, Qg;

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

  if (opflow->genbusvoltagetype == FIXED_WITHIN_QBOUNDS) {
    /*adding inequality Constraint Jacobian calc */

    for (i = 0; i < ps->nbus; i++) {
      bus = &ps->bus[i];
      if (bus->ide == PV_BUS || bus->ide == REF_BUS) {
        gloc = bus->startineqloc;
        xloc = bus->startxVloc;
        V = x[xloc + 1];

        Q = 0;
        Qmax = 0;
        Qmin = 0;
        for (k = 0; k < bus->ngen; k++) {
          ierr = PSBUSGetGen(bus, k, &gen);
          CHKERRQ(ierr);
          if (!gen->status)
            continue;
          loc = gen->startxpowloc;
          Qg = x[loc + 1];
          Q = Q + Qg;
          Qmax = Qmax + gen->qt;
          Qmin = Qmin + gen->qb;
          Vset = gen->vs;
          /* Partial derivative of eq1 and eq2 w.r.t. Qg */
          row[0] = gloc;
          row[1] = gloc + 1;
          col[0] = loc + 1;
          val[0] = Vset - V; // deq1_dQg
          val[1] = Vset - V; // deq2_dQg
          ierr = MatSetValues(Ji, 2, row, 1, col, val, ADD_VALUES);
          CHKERRQ(ierr);
        }

        val[0] = Qmax - Q; // deq1_dV
        val[1] = Qmin - Q; // deq2_dV
        row[0] = gloc;
        row[1] = gloc + 1;
        col[0] = xloc + 1;
        ierr = MatSetValues(Ji, 2, row, 1, col, val, ADD_VALUES);
        CHKERRQ(ierr);
      }
    }
  }

  if (!opflow->ignore_lineflow_constraints) {
    for (i = 0; i < ps->nline; i++) {
      line = &ps->line[i];
      if(!line->monitor_flow) continue;

      gloc = line->startineqloc;

      xloc = line->startxflowloc;
      
      Pf = x[xloc];     Qf = x[xloc + 1];
      Pt = x[xloc + 2]; Qt = x[xloc + 3];


      dSf2_dPf = 2 * Pf;
      dSf2_dQf = 2 * Qf;
      dSt2_dPt = 2 * Pt;
      dSt2_dQt = 2 * Qt;

      row[0] = gloc;


      col[0] = line->startxflowloc;
      col[1] = line->startxflowloc + 1;
      
      val[0] = dSf2_dPf;
      val[1] = dSf2_dQf;
      ierr = MatSetValues(Ji, 1, row, 2, col, val, ADD_VALUES);
      CHKERRQ(ierr);

      row[0] = gloc + 1;
      col[0] = line->startxflowloc + 2;
      col[1] = line->startxflowloc + 3;
      val[0] = dSt2_dPt;
      val[1] = dSt2_dQt;
      ierr = MatSetValues(Ji, 1, row, 2, col, val, ADD_VALUES);
      CHKERRQ(ierr);

      if (opflow->allow_lineflow_violation) {
        loc = line->startxslackloc;
        row[0] = gloc;
        col[0] = loc;
        val[0] = -1.0;

        ierr = MatSetValues(Ji, 1, row, 1, col, val, ADD_VALUES);
        CHKERRQ(ierr);

        row[0] = gloc + 1;
        col[0] = loc + 1;
        val[0] = -1.0;

        ierr = MatSetValues(Ji, 1, row, 1, col, val, ADD_VALUES);
        CHKERRQ(ierr);
      }
    }
    flps += opflow->nlinesmon *
      (160 + (20 * EXAGO_FLOPS_SINOP) + (20 * EXAGO_FLOPS_COSOP));
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

PetscErrorCode OPFLOWComputeConstraints_PBPOL2(OPFLOW opflow, Vec X, Vec G) {
  // empty function
  (void)opflow;
  (void)X;
  (void)G;
  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeObjective_PBPOL2(OPFLOW opflow, Vec X,
                                            PetscScalar *obj) {
  PetscErrorCode ierr;
  const PetscScalar *x;
  PS ps = opflow->ps;
  PetscInt i;
  PSBUS bus;
  PSGEN gen;
  PSLINE line;
  PetscInt loc;
  PetscInt k;
  PetscScalar Pg;
  double flps = 0.0;

  PetscFunctionBegin;

  ierr = VecGetArrayRead(X, &x);
  CHKERRQ(ierr);

  *obj = 0.0;

  for (i = 0; i < ps->nline; i++) {
    line = &ps->line[i];
    if (!opflow->ignore_lineflow_constraints && opflow->allow_lineflow_violation && line->monitor_flow) {
      loc = line->startxslackloc;
      // Slack variables for from and to side
      PetscScalar xsft_slack = 0.0, xstf_slack = 0.0;
      xsft_slack = x[loc];
      xstf_slack = x[loc + 1];

      *obj += opflow->lineflowviolation_penalty * (xsft_slack + xstf_slack);
      flps += 2.0;
    }
  }

  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    if (opflow->include_powerimbalance_variables) {
      loc = bus->startxpimbloc;
      PetscScalar Pimbplus, Pimbminus, Qimbplus, Qimbminus;
      loc = bus->startxpimbloc;
      Pimbplus = x[loc];
      Pimbminus = x[loc + 1];
      Qimbplus = x[loc + 2];
      Qimbminus = x[loc + 3];

      *obj += opflow->powerimbalance_penalty * ps->MVAbase *
              (Pimbplus + Pimbminus + Qimbplus + Qimbminus);
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
      PetscScalar Pdloss, Qdloss;
      for (k = 0; k < bus->nload; k++) {
        ierr = PSBUSGetLoad(bus, k, &load);
        CHKERRQ(ierr);
        if (!load->status)
          continue;
        loc = load->startxloadlossloc;
        Pdloss = x[loc];
        Qdloss = x[loc + 1];
        *obj += load->loss_cost * ps->MVAbase * (Pdloss + Qdloss);
        flps += 4.0;
      }
    }
  }
  ierr = VecRestoreArrayRead(X, &x);
  CHKERRQ(ierr);

  /* Multiply by the weight */
  *obj *= opflow->weight;

  PetscLogFlops(flps);
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeGradient_PBPOL2(OPFLOW opflow, Vec X, Vec grad) {
  PetscErrorCode ierr;
  const PetscScalar *x;
  PetscScalar *df;
  PS ps = opflow->ps;
  PetscInt i;
  PSBUS bus;
  PSGEN gen;
  PSLINE line;
  PetscInt loc;
  PetscInt k;
  PetscScalar Pg;
  double flps = 0.0;

  PetscFunctionBegin;
  ierr = VecGetArrayRead(X, &x);
  CHKERRQ(ierr);
  ierr = VecGetArray(grad, &df);
  CHKERRQ(ierr);

  for (i = 0; i < ps->nline; i++) {
    line = &ps->line[i];
    if (!opflow->ignore_lineflow_constraints && opflow->allow_lineflow_violation && line->monitor_flow) {
      loc = line->startxslackloc;

      df[loc] = opflow->lineflowviolation_penalty;
      df[loc + 1] = opflow->lineflowviolation_penalty;

      flps += 0.0;
    }
  }

  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    if (opflow->include_powerimbalance_variables) {
      loc = bus->startxpimbloc;

      df[loc] = opflow->powerimbalance_penalty * ps->MVAbase;
      df[loc + 1] = opflow->powerimbalance_penalty * ps->MVAbase;
      df[loc + 2] = opflow->powerimbalance_penalty * ps->MVAbase;
      df[loc + 3] = opflow->powerimbalance_penalty * ps->MVAbase;

      flps += 4.0;
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
      for (k = 0; k < bus->nload; k++) {
        ierr = PSBUSGetLoad(bus, k, &load);
        CHKERRQ(ierr);
        if (!load->status)
          continue;
        loc = load->startxloadlossloc;
        df[loc] = load->loss_cost * ps->MVAbase;
        df[loc + 1] = load->loss_cost * ps->MVAbase;
        flps += 2.0;
      }
    }
  }
  ierr = VecRestoreArrayRead(X, &x);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(grad, &df);
  CHKERRQ(ierr);

  /* Scale the objective by the weight (probability) of the problem */
  ierr = VecScale(grad, opflow->weight);
  CHKERRQ(ierr);

  PetscLogFlops(flps);
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeObjandGradient_PBPOL2(OPFLOW opflow, Vec X,
                                                 PetscScalar *obj, Vec Grad) {
  PetscErrorCode ierr;
  PetscFunctionBegin;

  ierr = OPFLOWComputeObjective_PBPOL2(opflow, X, obj);
  CHKERRQ(ierr);
  ierr = OPFLOWComputeGradient_PBPOL2(opflow, X, Grad);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/* Note: This function also sets the number of variables for each component,
 * i.e., bus, gen, load, and line */
PetscErrorCode OPFLOWModelSetNumVariables_PBPOL2(OPFLOW opflow,
                                                PetscInt *busnvar,
                                                PetscInt *branchnvar,
                                                PetscInt *nx) {
  PetscInt i, ngen, k;
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
    if(!line->status) continue;
    
    branchnvar[i] = line->nx = 0;
    /* 4 variables for the flow */
    branchnvar[i] += 4;
    line->nxflow = 4;
    if (!opflow->ignore_lineflow_constraints && opflow->allow_lineflow_violation && line->monitor_flow) {
      /* Two variables for line flow slacks
         - From side flow (Sft) and To side flow (Stf)
      */
      branchnvar[i] += 2;
      line->nxslack += 2;
    }
    line->nx = line->nxslack + line->nxflow;
    *nx += branchnvar[i];
  }

  /* Variables for the buses */
  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSIsGhosted(bus, &isghost);
    CHKERRQ(ierr);
    //    if(isghost) continue;
    busnvar[i] = 2; /* 2 variables for the bus voltages */
    bus->nxV = busnvar[i];

    if (opflow->include_powerimbalance_variables) {
      busnvar[i] += 4;
      bus->nxpimb = 4;
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
      busnvar[i] += 2; /* (2 variables for voltage + Pg, Qg for each gen) */
      gen->nxpow = 2;
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
        /* Load loss variables..Real and imaginary part of the load loss */
        busnvar[i] += 2;
        load->nxloadloss = 2;
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
PetscErrorCode OPFLOWModelSetNumConstraints_PBPOL2(OPFLOW opflow,
                                                  PetscInt *branchnconeq,
                                                  PetscInt *busnconeq,
                                                  PetscInt *nconeq,
                                                  PetscInt *nconineq) {
  (void)branchnconeq;
  (void)busnconeq;
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
    *nconeq += 2;
    bus->nconeq = 2;
    bus->nconineq = 0;

    if (opflow->use_agc) {
      if (bus->ide == REF_BUS) {
        ps->nconeq = ps->nconineq = 0;
      }
    }

    if (opflow->genbusvoltagetype == FIXED_WITHIN_QBOUNDS) {
      if (bus->ide == PV_BUS || bus->ide == REF_BUS) {
        *nconineq += 2;
        bus->nconineq = 2;
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

  for (i = 0; i < ps->nline; i++) {
    line = &ps->line[i];
    
    if(!line->status) continue;
    
    line->nconeq = 4; /* 4 equality constraints for line */
    *nconeq += 4;
    line->nconineq = 0;
    if (!opflow->ignore_lineflow_constraints && line->monitor_flow) {
      *nconineq += 2; /* Number of line flow inequality constraints */
      line->nconineq = 2;
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
PetscErrorCode OPFLOWComputeEqualityConstraintsHessian_PBPOL2(OPFLOW opflow,
                                                             Vec X, Vec Lambda,
                                                             Mat H) {
  PetscErrorCode ierr;
  PS ps = opflow->ps;
  PetscInt i, flps = 0;
  PSLINE line;
  PSBUS bus;
  PetscInt nconnlines;
  const PSLINE *connlines;
  const PSBUS *connbuses;
  PetscInt xloc, xlocf, xloct;
  PSBUS busf, bust;
  const PetscScalar *x;
  const PetscScalar *lambda;
  PetscInt gloc;
  PetscInt row[16], col[16];
  PetscScalar val[16];

  PetscFunctionBegin;

  ierr = VecGetArrayRead(X, &x);
  CHKERRQ(ierr);
  ierr = VecGetArrayRead(Lambda, &lambda);
  CHKERRQ(ierr);

  // For equality constraints (power flow) */
  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    gloc = bus->starteqloc;

    xloc = bus->startxVloc;

    row[0] = xloc + 1;
    col[0] = xloc + 1;

    val[0] = lambda[gloc] * 2 * bus->gl + lambda[gloc + 1] * (-2 * bus->bl);
    ierr = MatSetValues(H, 1, row, 1, col, val, ADD_VALUES);
    CHKERRQ(ierr);

    ierr = PSBUSGetSupportingLines(bus, &nconnlines, &connlines);
    CHKERRQ(ierr);
  }

  for (i = 0; i < ps->nline; i++) {
    line = &ps->line[i];
    if (!line->status) continue;

    PetscScalar Gff, Bff, Gft, Bft, Gtf, Btf, Gtt, Btt;
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

    PetscScalar Vmf, Vmt, thetaf, thetat, thetaft, thetatf;

    thetaf = x[xlocf];
    Vmf = x[xlocf + 1];
    thetat = x[xloct];
    Vmt = x[xloct + 1];
    thetaft = thetaf - thetat;
    thetatf = thetat - thetaf;
    
    PetscScalar dPf_dthetaf_dthetaf, dPf_dthetaf_dVmf, dPf_dthetaf_dthetat,
      dPf_dthetaf_dVmt;
    PetscScalar dPf_dVmf_dthetaf, dPf_dVmf_dVmf, dPf_dVmf_dthetat,
      dPf_dVmf_dVmt;
    PetscScalar dPf_dthetat_dthetaf, dPf_dthetat_dVmf, dPf_dthetat_dthetat,
      dPf_dthetat_dVmt;
    PetscScalar dPf_dVmt_dthetaf, dPf_dVmt_dVmf, dPf_dVmt_dthetat,
      dPf_dVmt_dVmt;
    
    /* dPf_dthetaf = Vmf*Vmt*(-Gft*sin(thetaft) + Bft*cos(thetaft)); */
    dPf_dthetaf_dthetaf =
      -Vmf * Vmt * (Gft * cos(thetaft) + Bft * sin(thetaft));
    dPf_dthetaf_dVmf = Vmt * (-Gft * sin(thetaft) + Bft * cos(thetaft));
    dPf_dthetaf_dthetat =
      Vmf * Vmt * (Gft * cos(thetaft) + Bft * sin(thetaft));
    dPf_dthetaf_dVmt = Vmf * (-Gft * sin(thetaft) + Bft * cos(thetaft));
    
    /* dPf_Vmf  = 2*Gff*Vmf + Vmt*(Gft*cos(thetaft) + Bft*sin(thetaft)); */
    dPf_dVmf_dthetaf = Vmt * (-Gft * sin(thetaft) + Bft * cos(thetaft));
    dPf_dVmf_dVmf = 2 * Gff;
    dPf_dVmf_dthetat = Vmt * (Gft * sin(thetaft) - Bft * cos(thetaft));
    dPf_dVmf_dVmt = (Gft * cos(thetaft) + Bft * sin(thetaft));
    
    /* dPf_dthetat = Vmf*Vmt*(Gft*sin(thetaft) - Bft*cos(thetaft)); */
    dPf_dthetat_dthetaf =
      Vmf * Vmt * (Gft * cos(thetaft) + Bft * sin(thetaft));
    dPf_dthetat_dVmf = Vmt * (Gft * sin(thetaft) - Bft * cos(thetaft));
    dPf_dthetat_dthetat =
      Vmf * Vmt * (-Gft * cos(thetaft) - Bft * sin(thetaft));
    dPf_dthetat_dVmt = Vmf * (Gft * sin(thetaft) - Bft * cos(thetaft));
    
    /* dPf_dVmt = Vmf*(Gft*cos(thetaft) + Bft*sin(thetaft)); */
    dPf_dVmt_dthetaf = Vmf * (-Gft * sin(thetaft) + Bft * cos(thetaft));
    dPf_dVmt_dVmf = (Gft * cos(thetaft) + Bft * sin(thetaft));
    dPf_dVmt_dthetat = Vmf * (Gft * sin(thetaft) - Bft * cos(thetaft));
    dPf_dVmt_dVmt = 0.0;
    
    PetscScalar dQf_dthetaf_dthetaf, dQf_dthetaf_dVmf, dQf_dthetaf_dthetat,
      dQf_dthetaf_dVmt;
    PetscScalar dQf_dVmf_dthetaf, dQf_dVmf_dVmf, dQf_dVmf_dthetat,
      dQf_dVmf_dVmt;
    PetscScalar dQf_dthetat_dthetaf, dQf_dthetat_dVmf, dQf_dthetat_dthetat,
      dQf_dthetat_dVmt;
    PetscScalar dQf_dVmt_dthetaf, dQf_dVmt_dVmf, dQf_dVmt_dthetat,
      dQf_dVmt_dVmt;
    
    /* dQf_dthetaf = Vmf*Vmt*(Bft*sin(thetaft) + Gft*cos(thetaft)); */
    dQf_dthetaf_dthetaf =
      Vmf * Vmt * (Bft * cos(thetaft) - Gft * sin(thetaft));
    dQf_dthetaf_dVmf = Vmt * (Bft * sin(thetaft) + Gft * cos(thetaft));
    dQf_dthetaf_dthetat =
      Vmf * Vmt * (-Bft * cos(thetaft) + Gft * sin(thetaft));
    dQf_dthetaf_dVmt = Vmf * (Bft * sin(thetaft) + Gft * cos(thetaft));
    
    /* dQf_dVmf = -2*Bff*Vmf + Vmt*(-Bft*cos(thetaft) + Gft*sin(thetaft));
     */
    dQf_dVmf_dthetaf = Vmt * (Bft * sin(thetaft) + Gft * cos(thetaft));
    dQf_dVmf_dVmf = -2 * Bff;
    dQf_dVmf_dthetat = Vmt * (-Bft * sin(thetaft) - Gft * cos(thetaft));
    dQf_dVmf_dVmt = (-Bft * cos(thetaft) + Gft * sin(thetaft));
    
    /* dQf_dthetat = Vmf*Vmt*(-Bft*sin(thetaft) - Gft*cos(thetaft)); */
    dQf_dthetat_dthetaf =
      Vmf * Vmt * (-Bft * cos(thetaft) + Gft * sin(thetaft));
    dQf_dthetat_dVmf = Vmt * (-Bft * sin(thetaft) - Gft * cos(thetaft));
    dQf_dthetat_dthetat =
      Vmf * Vmt * (Bft * cos(thetaft) - Gft * sin(thetaft));
    dQf_dthetat_dVmt = Vmf * (-Bft * sin(thetaft) - Gft * cos(thetaft));
    
    /* dQf_dVmt = Vmf*(-Bft*cos(thetaft) + Gft*sin(thetaft)); */
    dQf_dVmt_dthetaf = Vmf * (Bft * sin(thetaft) + Gft * cos(thetaft));
    dQf_dVmt_dVmf = (-Bft * cos(thetaft) + Gft * sin(thetaft));
    dQf_dVmt_dthetat = Vmf * (-Bft * sin(thetaft) - Gft * cos(thetaft));
    dQf_dVmt_dVmt = 0.0;

        PetscScalar dPt_dthetat_dthetat, dPt_dthetat_dVmt, dPt_dthetat_dthetaf,
      dPt_dthetat_dVmf;
    PetscScalar dPt_dVmt_dthetat, dPt_dVmt_dVmt, dPt_dVmt_dthetaf,
      dPt_dVmt_dVmf;
    PetscScalar dPt_dthetaf_dthetat, dPt_dthetaf_dVmt, dPt_dthetaf_dthetaf,
      dPt_dthetaf_dVmf;
    PetscScalar dPt_dVmf_dthetat, dPt_dVmf_dVmt, dPt_dVmf_dthetaf,
      dPt_dVmf_dVmf;
    
    /* dPt_dthetat = Vmf*Vmt*(-Gtf*sin(thetatf) + Btf*cos(thetatf)); */
    dPt_dthetat_dthetat =
      Vmf * Vmt * (-Gtf * cos(thetatf) - Btf * sin(thetatf));
    dPt_dthetat_dVmt = Vmf * (-Gtf * sin(thetatf) + Btf * cos(thetatf));
    dPt_dthetat_dthetaf =
      Vmf * Vmt * (Gtf * cos(thetatf) + Btf * sin(thetatf));
    dPt_dthetat_dVmf = Vmt * (-Gtf * sin(thetatf) + Btf * cos(thetatf));
    
    /* dPt_Vmt  = 2*Gtt*Vmt + Vmf*(Gtf*cos(thetatf) + Btf*sin(thetatf)); */
    dPt_dVmt_dthetat = Vmf * (-Gtf * sin(thetatf) + Bft * cos(thetatf));
    dPt_dVmt_dVmt = 2 * Gtt;
    dPt_dVmt_dthetaf = Vmf * (Gtf * sin(thetatf) - Btf * cos(thetatf));
    dPt_dVmt_dVmf = (Gtf * cos(thetatf) + Btf * sin(thetatf));
    
    /* dPt_dthetaf = Vmf*Vmt*(Gtf*sin(thetatf) - Btf*cos(thetatf)); */
    dPt_dthetaf_dthetat =
      Vmf * Vmt * (Gtf * cos(thetatf) + Btf * sin(thetatf));
    dPt_dthetaf_dVmt = Vmf * (Gtf * sin(thetatf) - Btf * cos(thetatf));
    dPt_dthetaf_dthetaf =
      Vmf * Vmt * (-Gtf * cos(thetatf) - Btf * sin(thetatf));
    dPt_dthetaf_dVmf = Vmt * (Gtf * sin(thetatf) - Btf * cos(thetatf));
    
    /* dPt_dVmf = Vmt*(Gtf*cos(thetatf) + Btf*sin(thetatf)); */
    dPt_dVmf_dthetat = Vmt * (-Gtf * sin(thetatf) + Btf * cos(thetatf));
    dPt_dVmf_dVmt = (Gtf * cos(thetatf) + Btf * sin(thetatf));
    dPt_dVmf_dthetaf = Vmt * (Gtf * sin(thetatf) - Btf * cos(thetatf));
    dPt_dVmf_dVmf = 0.0;
    
    PetscScalar dQt_dthetaf_dthetaf, dQt_dthetaf_dVmf, dQt_dthetaf_dthetat,
      dQt_dthetaf_dVmt;
    PetscScalar dQt_dVmf_dthetaf, dQt_dVmf_dVmf, dQt_dVmf_dthetat,
            dQt_dVmf_dVmt;
    PetscScalar dQt_dthetat_dthetaf, dQt_dthetat_dVmf, dQt_dthetat_dthetat,
      dQt_dthetat_dVmt;
    PetscScalar dQt_dVmt_dthetaf, dQt_dVmt_dVmf, dQt_dVmt_dthetat,
      dQt_dVmt_dVmt;
    
    /* dQt_dthetat = Vmf*Vmt*(Btf*sin(thetatf) + Gtf*cos(thetatf)); */
    dQt_dthetat_dthetat =
      Vmf * Vmt * (Btf * cos(thetatf) - Gtf * sin(thetatf));
    dQt_dthetat_dVmt = Vmf * (Btf * sin(thetatf) + Gtf * cos(thetatf));
    dQt_dthetat_dthetaf =
      Vmf * Vmt * (-Btf * cos(thetatf) + Gtf * sin(thetatf));
    dQt_dthetat_dVmf = Vmf * (Btf * sin(thetatf) + Gtf * cos(thetatf));
    
    /* dQt_dVmt = -2*Btt*Vmt + Vmf*(-Btf*cos(thetatf) + Gtf*sin(thetatf));
     */
    dQt_dVmt_dthetat = Vmf * (Btf * sin(thetatf) + Gtf * cos(thetatf));
    dQt_dVmt_dVmt = -2 * Btt;
    dQt_dVmt_dthetaf = Vmf * (-Btf * sin(thetatf) + Gtf * cos(thetatf));
    dQt_dVmt_dVmf = (-Btf * cos(thetatf) + Gtf * sin(thetatf));
    
    /* dQt_dthetaf = Vmf*Vmt*(-Btf*sin(thetatf) - Gtf*cos(thetatf)); */
    dQt_dthetaf_dthetat =
      Vmf * Vmt * (-Btf * cos(thetatf) + Gtf * sin(thetatf));
    dQt_dthetaf_dVmt = Vmf * (-Btf * sin(thetatf) - Gtf * cos(thetatf));
    dQt_dthetaf_dthetaf =
      Vmf * Vmt * (Btf * cos(thetatf) - Gtf * sin(thetatf));
    dQt_dthetaf_dVmf = Vmt * (-Btf * sin(thetatf) - Gtf * cos(thetatf));
    
    /* dQt_dVmf = Vmt*(-Btf*cos(thetatf) + Gtf*sin(thetatf)); */
    dQt_dVmf_dthetat = Vmt * (Btf * sin(thetatf) + Gtf * cos(thetatf));
    dQt_dVmf_dVmt = (-Btf * cos(thetatf) + Gtf * sin(thetatf));
    dQt_dVmf_dthetaf = Vmt * (-Btf * sin(thetatf) - Gtf * cos(thetatf));
    dQt_dVmf_dVmf = 0.0;

    
    row[0] = xlocf;
    row[1] = xlocf + 1;
    
    col[0] = xlocf;
    col[1] = xlocf + 1;
    col[2] = xloct;
    col[3] = xloct + 1;

    gloc = line->starteqloc;
    
    val[0] = val[1] = val[2] = val[3] = val[4] = val[5] = val[6] = val[7] =
      0.0;
    
    val[0] = lambda[gloc]     * dPf_dthetaf_dthetaf +
             lambda[gloc + 1] * dQf_dthetaf_dthetaf +
             lambda[gloc + 2] * dPt_dthetaf_dthetaf +
             lambda[gloc + 3] * dQt_dthetaf_dthetaf;
    
    val[1] = lambda[gloc]     * dPf_dthetaf_dVmf +
             lambda[gloc + 1] * dQf_dthetaf_dVmf +
             lambda[gloc + 2] * dPt_dthetaf_dVmf +
             lambda[gloc + 3] * dQt_dthetaf_dVmf;

    val[2] = lambda[gloc]     * dPf_dthetaf_dthetat +
             lambda[gloc + 1] * dQf_dthetaf_dthetat +
             lambda[gloc + 2] * dPt_dthetaf_dthetat +
             lambda[gloc + 3] * dQt_dthetaf_dthetat;

    val[3] = lambda[gloc]     * dPf_dthetaf_dVmt +
             lambda[gloc + 1] * dQf_dthetaf_dVmt +
             lambda[gloc + 2] * dPt_dthetaf_dVmt +
             lambda[gloc + 3] * dQt_dthetaf_dVmt;

    val[4] = lambda[gloc]     * dPf_dVmf_dthetaf +
             lambda[gloc + 1] * dQf_dVmf_dthetaf +
             lambda[gloc + 2] * dPt_dVmf_dthetaf +
             lambda[gloc + 3] * dQt_dVmf_dthetaf;

    val[5] = lambda[gloc]     * dPf_dVmf_dVmf +
             lambda[gloc + 1] * dQf_dVmf_dVmf +
             lambda[gloc + 2] * dPt_dVmf_dVmf +
             lambda[gloc + 3] * dQt_dVmf_dVmf;

    val[6] = lambda[gloc]     * dPf_dVmf_dthetat +
             lambda[gloc + 1] * dQf_dVmf_dthetat +
             lambda[gloc + 2] * dPt_dVmf_dthetat +
             lambda[gloc + 3] * dQt_dVmf_dthetat;

    val[7] = lambda[gloc]     * dPf_dVmf_dVmt +
             lambda[gloc + 1] * dQf_dVmf_dVmt +
             lambda[gloc + 2] * dPt_dVmf_dVmt +
             lambda[gloc + 3] * dQt_dVmf_dVmt;
    
    ierr = MatSetValues(H, 2, row, 4, col, val, ADD_VALUES);
    CHKERRQ(ierr);

    row[0] = xloct;
    row[1] = xloct + 1;
    
    col[0] = xlocf;
    col[1] = xlocf + 1;
    col[2] = xloct;
    col[3] = xloct + 1;

    
    val[0] = lambda[gloc]     * dPf_dthetat_dthetaf +
             lambda[gloc + 1] * dQf_dthetat_dthetaf +
             lambda[gloc + 2] * dPt_dthetat_dthetaf +
             lambda[gloc + 3] * dQt_dthetat_dthetaf;
    
    val[1] = lambda[gloc]     * dPf_dthetat_dVmf +
             lambda[gloc + 1] * dQf_dthetat_dVmf +
             lambda[gloc + 2] * dPt_dthetat_dVmf +
             lambda[gloc + 3] * dQt_dthetat_dVmf;

    val[2] = lambda[gloc]     * dPf_dthetat_dthetat +
             lambda[gloc + 1] * dQf_dthetat_dthetat +
             lambda[gloc + 2] * dPt_dthetat_dthetat +
             lambda[gloc + 3] * dQt_dthetat_dthetat;

    val[3] = lambda[gloc]     * dPf_dthetat_dVmt +
             lambda[gloc + 1] * dQf_dthetat_dVmt +
             lambda[gloc + 2] * dPt_dthetat_dVmt +
             lambda[gloc + 3] * dQt_dthetat_dVmt;

    val[4] = lambda[gloc]     * dPf_dVmt_dthetaf +
             lambda[gloc + 1] * dQf_dVmt_dthetaf +
             lambda[gloc + 2] * dPt_dVmt_dthetaf +
             lambda[gloc + 3] * dQt_dVmt_dthetaf;

    val[5] = lambda[gloc]     * dPf_dVmt_dVmf +
             lambda[gloc + 1] * dQf_dVmt_dVmf +
             lambda[gloc + 2] * dPt_dVmt_dVmf +
             lambda[gloc + 3] * dQt_dVmt_dVmf;

    val[6] = lambda[gloc]     * dPf_dVmt_dthetat +
             lambda[gloc + 1] * dQf_dVmt_dthetat +
             lambda[gloc + 2] * dPt_dVmt_dthetat +
             lambda[gloc + 3] * dQt_dVmt_dthetat;

    val[7] = lambda[gloc]     * dPf_dVmt_dVmt +
             lambda[gloc + 1] * dQf_dVmt_dVmt +
             lambda[gloc + 2] * dPt_dVmt_dVmt +
             lambda[gloc + 3] * dQt_dVmt_dVmt;
    
    ierr = MatSetValues(H, 2, row, 4, col, val, ADD_VALUES);
    CHKERRQ(ierr);
  }

  ierr = VecRestoreArrayRead(X, &x);
  CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(Lambda, &lambda);
  CHKERRQ(ierr);

  PetscLogFlops(flps);
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
PetscErrorCode OPFLOWComputeInequalityConstraintsHessian_PBPOL2(OPFLOW opflow,
                                                               Vec X,
                                                               Vec Lambda,
                                                               Mat H) {
  PetscErrorCode ierr;
  PS ps = opflow->ps;
  PetscInt i, k, flps = 0;
  PSLINE line;
  PetscInt xloc;
  const PetscScalar *x;
  const PetscScalar *lambda;
  PetscInt gloc;
  PetscInt row[12], col[12];
  PetscScalar val[12];
  PSBUS bus;
  PSGEN gen;
  PetscInt loc;
  PetscScalar Q, Qmax, Qmin, Qg;

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
          // Pg = x[gen->startxpowloc];
          // delPg = x[gen->startxpdevloc];
          // delP = x[ps->startxloc];

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

  if (opflow->genbusvoltagetype == FIXED_WITHIN_QBOUNDS) {
    for (i = 0; i < ps->nbus; i++) {
      bus = &ps->bus[i];
      if (bus->ide == PV_BUS || bus->ide == REF_BUS) {
        gloc = bus->startineqloc;
        xloc = bus->startxVloc;
        // xlocglob = bus->startxVlocglob;
        // V = x[xloc + 1];

        Q = 0;
        Qmax = 0;
        Qmin = 0;
        for (k = 0; k < bus->ngen; k++) {
          ierr = PSBUSGetGen(bus, k, &gen);
          CHKERRQ(ierr);
          if (!gen->status)
            continue;
          loc = gen->startxpowloc;
          Qg = x[loc + 1];
          Q = Q + Qg;
          Qmax = Qmax + gen->qt;
          Qmin = Qmin + gen->qb;
          // Vset = gen->vs;

          row[0] = loc + 1;
          col[0] = xloc + 1;
          val[0] = -(
              lambda[gloc] +
              lambda[gloc + 1]); // lam_eq1*d2eq1_dQg_dV + lam_eq2*d2eq2_dQg_dV
          ierr = MatSetValues(H, 1, row, 1, col, val, ADD_VALUES);
          CHKERRQ(ierr);
          row[0] = xloc + 1;
          col[0] = loc + 1;
          val[0] = -(
              lambda[gloc] +
              lambda[gloc + 1]); // lam_eq1* d2eq1_dQg_dV + lam_eq2*d2eq2_dV_dQg
          ierr = MatSetValues(H, 1, row, 1, col, val, ADD_VALUES);
          CHKERRQ(ierr);
        }
      }
    }
  }

  // for the part of line constraints
  if (!opflow->ignore_lineflow_constraints) {
    for (i = 0; i < ps->nline; i++) {
      line = &ps->line[i];
      if(!line->monitor_flow) continue;

      gloc = line->startineqloc;
      xloc = line->startxflowloc;

      row[0] = xloc;
      col[0] = xloc;
      val[0] = lambda[gloc]*2.0;
      ierr = MatSetValues(H, 1, row, 1, col, val, ADD_VALUES);
      CHKERRQ(ierr);

      row[0] = xloc+1;
      col[0] = xloc+1;
      val[0] = lambda[gloc]*2.0;
      ierr = MatSetValues(H, 1, row, 1, col, val, ADD_VALUES);
      CHKERRQ(ierr);

      row[0] = xloc + 2;
      col[0] = xloc + 2;
      val[0] = lambda[gloc + 1]*2.0;
      ierr = MatSetValues(H, 1, row, 1, col, val, ADD_VALUES);
      CHKERRQ(ierr);

      row[0] = xloc + 3;
      col[0] = xloc + 3;
      val[0] = lambda[gloc + 1]*2.0;
      ierr = MatSetValues(H, 1, row, 1, col, val, ADD_VALUES);
      CHKERRQ(ierr);
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
PetscErrorCode OPFLOWComputeObjectiveHessian_PBPOL2(OPFLOW opflow, Vec X,
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
  PetscScalar weight = opflow->weight;
  PetscInt flps = 0;

  PetscFunctionBegin;

  ierr = VecGetArrayRead(X, &x);
  CHKERRQ(ierr);

  // for the part of objective
  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    ierr = PSBUSGetVariableGlobalLocation(bus, &xlocglob);
    CHKERRQ(ierr);

    if (opflow->include_powerimbalance_variables) {
      xlocglob = bus->startxpimblocglob;

      row[0] = xlocglob;
      col[0] = xlocglob;
      val[0] = 0.0;

      ierr = MatSetValues(H, 1, row, 1, col, val, ADD_VALUES);
      CHKERRQ(ierr);

      row[0] = xlocglob + 1;
      col[0] = xlocglob + 1;
      val[0] = 0.0;

      ierr = MatSetValues(H, 1, row, 1, col, val, ADD_VALUES);
      CHKERRQ(ierr);
    }

    for (k = 0; k < bus->ngen; k++) {
      ierr = PSBUSGetGen(bus, k, &gen);
      CHKERRQ(ierr);
      if (!gen->status)
        continue;

      if (opflow->objectivetype == MIN_GEN_COST) {
        xlocglob = gen->startxpowlocglob;

        row[0] = xlocglob;
        col[0] = xlocglob;

        val[0] = weight * obj_factor * 2.0 * gen->cost_alpha * ps->MVAbase *
                 ps->MVAbase;
        ierr = MatSetValues(H, 1, row, 1, col, val, ADD_VALUES);
        CHKERRQ(ierr);
        flps += 4;
      } else if (opflow->objectivetype == MIN_GENSETPOINT_DEVIATION) {
        xlocglob = gen->startxpdevlocglob;
        row[0] = xlocglob;
        col[0] = xlocglob;
        val[0] = weight * obj_factor * 2.0;
        ierr = MatSetValues(H, 1, row, 1, col, val, ADD_VALUES);
        CHKERRQ(ierr);

        flps += 1;
      }
    }

    if (opflow->include_loadloss_variables) {
      PSLOAD load;
      for (k = 0; k < bus->nload; k++) {
        ierr = PSBUSGetLoad(bus, k, &load);
        CHKERRQ(ierr);
        if (!load->status)
          continue;
        xlocglob = load->startxloadlosslocglob;
        row[0] = xlocglob;
        col[0] = xlocglob;
        val[0] = 0.0;
        ierr = MatSetValues(H, 1, row, 1, col, val, ADD_VALUES);
        CHKERRQ(ierr);

        row[0] = xlocglob + 1;
        col[0] = xlocglob + 1;
        val[0] = 0.0;
        ierr = MatSetValues(H, 1, row, 1, col, val, ADD_VALUES);
        CHKERRQ(ierr);
      }
    }
  }

  ierr = VecRestoreArrayRead(X, &x);
  CHKERRQ(ierr);
  PetscLogFlops(flps);
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeHessian_PBPOL2(OPFLOW opflow, Vec X, Vec Lambdae,
                                          Vec Lambdai, Mat H) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = MatZeroEntries(H);
  CHKERRQ(ierr);

  /* Objective function Hessian */
  ierr = OPFLOWComputeObjectiveHessian_PBPOL2(opflow, X, H);
  CHKERRQ(ierr);

  /* Equality constraints Hessian */
  ierr = OPFLOWComputeEqualityConstraintsHessian_PBPOL2(opflow, X, Lambdae, H);
  CHKERRQ(ierr);

  /* Inequality constraints Hessian */
  if (opflow->nconineq) {
    ierr =
        OPFLOWComputeInequalityConstraintsHessian_PBPOL2(opflow, X, Lambdai, H);
    CHKERRQ(ierr);
  }

  ierr = MatAssemblyBegin(H, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = MatAssemblyEnd(H, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolutionToPS_PBPOL2(OPFLOW opflow) {
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
      PetscScalar Pimbplus, Pimbminus, Qimbplus, Qimbminus;
      Pimbplus = x[loc];
      Pimbminus = x[loc + 1];
      Qimbplus = x[loc + 2];
      Qimbminus = x[loc + 3];
      bus->pimb = Pimbplus - Pimbminus;
      bus->qimb = Qimbplus - Qimbminus;
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
  }

  for (i = 0; i < opflow->nlinesmon; i++) {
    line = &ps->line[opflow->linesmon[i]];
    gloc = line->startineqloc;
    line->mult_sf = lambdai[gloc];
    line->mult_st = lambdai[gloc + 1];
  }

  ierr = VecRestoreArrayRead(X, &x);
  CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(Lambda, &lambda);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWModelSetUp_PBPOL2(OPFLOW opflow) {
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
    if(!line->status) continue;

    ierr = PSLINEGetVariableLocation(line, &loc);
    line->startxflowloc = loc;
    CHKERRQ(ierr);
    /* Set starting location for slack variable */
    if (!opflow->ignore_lineflow_constraints) {
      if(opflow->allow_lineflow_violation) {
	line->startxslackloc = line->startxflowloc + line->nxflow;
      }
      line->startineqloc = ineqloc;
      ineqloc += line->nconineq;
    }
    line->starteqloc = eqloc;
    eqloc += line->nconeq;
  }

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWModelCreate_PBPOL2(OPFLOW opflow) {
  PBPOL pbpol;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = PetscCalloc1(1, &pbpol);
  CHKERRQ(ierr);

  opflow->model = pbpol;

  /* Inherit Ops */
  opflow->modelops.destroy = OPFLOWModelDestroy_PBPOL2;
  opflow->modelops.setnumvariables = OPFLOWModelSetNumVariables_PBPOL2;
  opflow->modelops.setnumconstraints = OPFLOWModelSetNumConstraints_PBPOL2;
  opflow->modelops.setvariablebounds = OPFLOWSetVariableBounds_PBPOL2;
  opflow->modelops.setconstraintbounds = OPFLOWSetConstraintBounds_PBPOL2;
  opflow->modelops.setvariableandconstraintbounds =
      OPFLOWSetVariableandConstraintBounds_PBPOL2;
  opflow->modelops.setinitialguess = OPFLOWSetInitialGuess_PBPOL2;
  opflow->modelops.setup = OPFLOWModelSetUp_PBPOL2;
  opflow->modelops.computeequalityconstraints =
      OPFLOWComputeEqualityConstraints_PBPOL2;
  opflow->modelops.computeinequalityconstraints =
      OPFLOWComputeInequalityConstraints_PBPOL2;
  opflow->modelops.computeequalityconstraintjacobian =
      OPFLOWComputeEqualityConstraintJacobian_PBPOL2;
  opflow->modelops.computeinequalityconstraintjacobian =
      OPFLOWComputeInequalityConstraintJacobian_PBPOL2;
  opflow->modelops.computehessian = OPFLOWComputeHessian_PBPOL2;
  opflow->modelops.computeobjandgradient = OPFLOWComputeObjandGradient_PBPOL2;
  opflow->modelops.computeobjective = OPFLOWComputeObjective_PBPOL2;
  opflow->modelops.computegradient = OPFLOWComputeGradient_PBPOL2;
  opflow->modelops.solutiontops = OPFLOWSolutionToPS_PBPOL2;

  PetscFunctionReturn(0);
}
