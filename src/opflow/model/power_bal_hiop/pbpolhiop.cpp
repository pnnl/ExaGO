#include <exago_config.h>

#if defined(EXAGO_ENABLE_HIOP)

#include "pbpolhiop.h"
#include <private/opflowimpl.h>
#include <private/psimpl.h>

/* Functions to create and destroy data arrays for different
   component classes
*/
PetscErrorCode DestroyBusParams(OPFLOW opflow, BUSParams *busparams) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = PetscFree(busparams->isref);
  CHKERRQ(ierr);
  ierr = PetscFree(busparams->isisolated);
  CHKERRQ(ierr);
  ierr = PetscFree(busparams->ispvpq);
  CHKERRQ(ierr);
  ierr = PetscFree(busparams->vmin);
  CHKERRQ(ierr);
  ierr = PetscFree(busparams->vmax);
  CHKERRQ(ierr);
  ierr = PetscFree(busparams->va);
  CHKERRQ(ierr);
  ierr = PetscFree(busparams->vm);
  CHKERRQ(ierr);
  ierr = PetscFree(busparams->gl);
  CHKERRQ(ierr);
  ierr = PetscFree(busparams->bl);
  CHKERRQ(ierr);
  ierr = PetscFree(busparams->xidx);
  CHKERRQ(ierr);
  ierr = PetscFree(busparams->gidx);
  CHKERRQ(ierr);
  if (opflow->include_powerimbalance_variables) {
    ierr = PetscFree(busparams->xidxpimb);
    CHKERRQ(ierr);
    ierr = PetscFree(busparams->jacsp_idx);
    CHKERRQ(ierr);
    ierr = PetscFree(busparams->jacsq_idx);
    CHKERRQ(ierr);
    ierr = PetscFree(busparams->hesssp_idx);
    CHKERRQ(ierr);
    ierr = PetscFree(busparams->powerimbalance_penalty);
    CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

/* Create data for buses that is used in different computations */
PetscErrorCode CreateBusParams(OPFLOW opflow, BUSParams *busparams) {
  PS ps = opflow->ps;
  PetscInt loc, gloc = 0;
  PSBUS bus;
  PetscInt i;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  busparams->nbus = ps->nbus;

  /* Allocate the arrays */
  ierr = PetscCalloc1(busparams->nbus, &busparams->isref);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(busparams->nbus, &busparams->isisolated);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(busparams->nbus, &busparams->ispvpq);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(busparams->nbus, &busparams->vmin);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(busparams->nbus, &busparams->vmax);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(busparams->nbus, &busparams->va);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(busparams->nbus, &busparams->vm);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(busparams->nbus, &busparams->gl);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(busparams->nbus, &busparams->bl);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(busparams->nbus, &busparams->xidx);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(busparams->nbus, &busparams->gidx);
  CHKERRQ(ierr);
  if (opflow->include_powerimbalance_variables) {
    ierr = PetscCalloc1(busparams->nbus, &busparams->powerimbalance_penalty);
    CHKERRQ(ierr);
    ierr = PetscCalloc1(busparams->nbus, &busparams->xidxpimb);
    CHKERRQ(ierr);
    ierr = PetscCalloc1(busparams->nbus, &busparams->jacsp_idx);
    CHKERRQ(ierr);
    ierr = PetscCalloc1(busparams->nbus, &busparams->jacsq_idx);
    CHKERRQ(ierr);
    ierr = PetscCalloc1(busparams->nbus, &busparams->hesssp_idx);
    CHKERRQ(ierr);
  }

  /* Populate the arrays */
  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    loc = bus->startxVloc;

    busparams->xidx[i] = opflow->idxn2sd_map[loc];
    busparams->gidx[i] = bus->starteqloc;

    if (bus->ide == REF_BUS)
      busparams->isref[i] = 1;
    else if (bus->ide == ISOLATED_BUS)
      busparams->isisolated[i] = 1;
    else
      busparams->ispvpq[i] = 1;

    if (opflow->genbusvoltagetype == FIXED_AT_SETPOINT) {
      if (bus->ide == REF_BUS || bus->ide == PV_BUS) {
        /* Hold voltage at reference and PV buses */
        busparams->vmin[i] = bus->vm;
        busparams->vmax[i] = bus->vm;
      } else {
        busparams->vmin[i] = bus->Vmin;
        busparams->vmax[i] = bus->Vmax;
      }
    } else {
      busparams->vmin[i] = bus->Vmin;
      busparams->vmax[i] = bus->Vmax;
    }
    busparams->vm[i] = bus->vm;
    busparams->va[i] = bus->va;
    busparams->gl[i] = bus->gl;
    busparams->bl[i] = bus->bl;

    if (opflow->include_powerimbalance_variables) {
      busparams->powerimbalance_penalty[i] = opflow->powerimbalance_penalty;
      loc = bus->startxpimbloc;
      busparams->xidxpimb[i] = opflow->idxn2sd_map[loc];
      // busparams->xidxpimb[i+1] = opflow->idxn2sd_map[loc+1];
    }
  }

  PetscFunctionReturn(0);
}

PetscErrorCode DestroyLineParams(OPFLOW opflow, LINEParams *lineparams) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = PetscFree(lineparams->Gff);
  CHKERRQ(ierr);
  ierr = PetscFree(lineparams->Bff);
  CHKERRQ(ierr);
  ierr = PetscFree(lineparams->Gft);
  CHKERRQ(ierr);
  ierr = PetscFree(lineparams->Bft);
  CHKERRQ(ierr);
  ierr = PetscFree(lineparams->Gtf);
  CHKERRQ(ierr);
  ierr = PetscFree(lineparams->Btf);
  CHKERRQ(ierr);
  ierr = PetscFree(lineparams->Gtt);
  CHKERRQ(ierr);
  ierr = PetscFree(lineparams->Btt);
  CHKERRQ(ierr);
  ierr = PetscFree(lineparams->rateA);
  CHKERRQ(ierr);

  ierr = PetscFree(lineparams->xidxf);
  CHKERRQ(ierr);
  ierr = PetscFree(lineparams->xidxt);
  CHKERRQ(ierr);

  ierr = PetscFree(lineparams->geqidxf);
  CHKERRQ(ierr);
  ierr = PetscFree(lineparams->geqidxt);
  CHKERRQ(ierr);
  ierr = PetscFree(lineparams->gineqidx);
  CHKERRQ(ierr);
  ierr = PetscFree(lineparams->gbineqidx);
  CHKERRQ(ierr);

  ierr = PetscFree(lineparams->linelimidx);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/* Create data for lines that is used in different computations */
PetscErrorCode CreateLineParams(OPFLOW opflow, LINEParams *lineparams) {
  PS ps = opflow->ps;
  PetscInt linei = 0, linelimi = 0;
  PSLINE line;
  PetscInt i;
  PetscErrorCode ierr;
  const PSBUS *connbuses;
  PSBUS busf, bust;

  PetscFunctionBegin;

  ierr = PSGetNumActiveLines(ps, &lineparams->nlineON, NULL);
  CHKERRQ(ierr);

  lineparams->nlinelim = 0;
  /* Get the number of lines that are active and have finite limits. These lines
     will be only considered in inequality constraints */
  if (opflow->nconineq) {
    for (i = 0; i < ps->nline; i++) {
      line = &ps->line[i];

      if (!line->status || line->rateA > 1e5)
        continue;
      lineparams->nlinelim++;
    }
  }
  /* Allocate arrays */
  ierr = PetscCalloc1(lineparams->nlineON, &lineparams->Gff);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(lineparams->nlineON, &lineparams->Bff);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(lineparams->nlineON, &lineparams->Gft);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(lineparams->nlineON, &lineparams->Bft);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(lineparams->nlineON, &lineparams->Gtf);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(lineparams->nlineON, &lineparams->Btf);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(lineparams->nlineON, &lineparams->Gtt);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(lineparams->nlineON, &lineparams->Btt);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(lineparams->nlineON, &lineparams->rateA);
  CHKERRQ(ierr);

  ierr = PetscCalloc1(lineparams->nlineON, &lineparams->xidxf);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(lineparams->nlineON, &lineparams->xidxt);
  CHKERRQ(ierr);

  ierr = PetscCalloc1(lineparams->nlineON, &lineparams->geqidxf);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(lineparams->nlineON, &lineparams->geqidxt);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(lineparams->nlineON, &lineparams->gineqidx);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(lineparams->nlineON, &lineparams->gbineqidx);
  CHKERRQ(ierr);

  if (opflow->nconineq) {
    ierr = PetscCalloc1(lineparams->nlinelim, &lineparams->linelimidx);
    CHKERRQ(ierr);
  }

  /* Populate arrays */
  for (i = 0; i < ps->nline; i++) {
    line = &ps->line[i];

    if (!line->status)
      continue;

    lineparams->Gff[linei] = line->yff[0];
    lineparams->Bff[linei] = line->yff[1];
    lineparams->Gft[linei] = line->yft[0];
    lineparams->Bft[linei] = line->yft[1];
    lineparams->Gtf[linei] = line->ytf[0];
    lineparams->Btf[linei] = line->ytf[1];
    lineparams->Gtt[linei] = line->ytt[0];
    lineparams->Btt[linei] = line->ytt[1];
    lineparams->rateA[linei] = line->rateA;

    ierr = PSLINEGetConnectedBuses(line, &connbuses);
    CHKERRQ(ierr);
    busf = connbuses[0];
    bust = connbuses[1];

    int xidxf, xidxt;

    xidxf = busf->startxVloc;
    xidxt = bust->startxVloc;

    lineparams->xidxf[linei] = opflow->idxn2sd_map[xidxf];
    lineparams->xidxt[linei] = opflow->idxn2sd_map[xidxt];

    /*
       Each bus has two equality (balance) constraints, hence the use of
       coefficient 2 to map the location of the equality constraint for the bus
    */
    lineparams->geqidxf[linei] = busf->starteqloc;
    lineparams->geqidxt[linei] = bust->starteqloc;

    if (opflow->nconineq) {
      if (line->rateA < 1e5) {
        lineparams->gbineqidx[linelimi] = opflow->nconeq + line->startineqloc;
        lineparams->gineqidx[linelimi] = line->startineqloc;
        lineparams->linelimidx[linelimi] = linei;
        linelimi++;
      }
    }
    linei++;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode DestroyLoadParams(OPFLOW opflow, LOADParams *loadparams) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscFree(loadparams->pl);
  CHKERRQ(ierr);
  ierr = PetscFree(loadparams->ql);
  CHKERRQ(ierr);
  ierr = PetscFree(loadparams->loadloss_penalty);
  CHKERRQ(ierr);
  ierr = PetscFree(loadparams->xidx);
  CHKERRQ(ierr);
  ierr = PetscFree(loadparams->gidx);
  CHKERRQ(ierr);
  if (opflow->include_loadloss_variables) {
    ierr = PetscFree(loadparams->jacsp_idx);
    CHKERRQ(ierr);
    ierr = PetscFree(loadparams->jacsq_idx);
    CHKERRQ(ierr);
    ierr = PetscFree(loadparams->hesssp_idx);
    CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

/* Create data for loads that is used in different computations */
PetscErrorCode CreateLoadParams(OPFLOW opflow, LOADParams *loadparams) {
  PS ps = opflow->ps;
  PetscInt loc, loadi = 0;
  PSLOAD load;
  PSBUS bus;
  PetscInt i, j;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  /* Get the number of active generators (STATUS ON) */
  ierr = PSGetNumLoads(ps, &loadparams->nload, NULL);
  CHKERRQ(ierr);

  /* Allocate arrays */
  ierr = PetscCalloc1(loadparams->nload, &loadparams->pl);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(loadparams->nload, &loadparams->ql);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(loadparams->nload, &loadparams->loadloss_penalty);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(loadparams->nload, &loadparams->xidx);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(loadparams->nload, &loadparams->gidx);
  CHKERRQ(ierr);
  if (opflow->include_loadloss_variables) {
    ierr = PetscCalloc1(loadparams->nload, &loadparams->jacsp_idx);
    CHKERRQ(ierr);
    ierr = PetscCalloc1(loadparams->nload, &loadparams->jacsq_idx);
    CHKERRQ(ierr);
    ierr = PetscCalloc1(loadparams->nload, &loadparams->hesssp_idx);
    CHKERRQ(ierr);
  }

  /* Insert data in loadparams */
  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSGetVariableLocation(bus, &loc);
    CHKERRQ(ierr);
    for (j = 0; j < bus->nload; j++) {
      ierr = PSBUSGetLoad(bus, j, &load);
      CHKERRQ(ierr);

      loc = load->startxloadlossloc;

      loadparams->pl[loadi] = load->pl;
      loadparams->ql[loadi] = load->ql;
      loadparams->loadloss_penalty[loadi] = opflow->loadloss_penalty;

      loadparams->xidx[loadi] = opflow->idxn2sd_map[loc];
      loadparams->gidx[loadi] = bus->starteqloc;
      loadi++;
    }
  }

  PetscFunctionReturn(0);
}

PetscErrorCode DestroyGenParams(OPFLOW opflow, GENParams *genparams) {
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = PetscFree(genparams->cost_alpha);
  CHKERRQ(ierr);
  ierr = PetscFree(genparams->cost_beta);
  CHKERRQ(ierr);
  ierr = PetscFree(genparams->cost_gamma);
  CHKERRQ(ierr);

  ierr = PetscFree(genparams->pt);
  CHKERRQ(ierr);
  ierr = PetscFree(genparams->pb);
  CHKERRQ(ierr);
  ierr = PetscFree(genparams->qt);
  CHKERRQ(ierr);
  ierr = PetscFree(genparams->qb);
  CHKERRQ(ierr);

  if (opflow->has_gensetpoint) {
    ierr = PetscFree(genparams->geqidxgen);
    CHKERRQ(ierr);
    ierr = PetscFree(genparams->gineqidxgen);
    CHKERRQ(ierr);
    ierr = PetscFree(genparams->gbineqidxgen);
    CHKERRQ(ierr);
    ierr = PetscFree(genparams->pgs);
    CHKERRQ(ierr);
    ierr = PetscFree(genparams->eqjacspgen_idx);
    CHKERRQ(ierr);
    ierr = PetscFree(genparams->ineqjacspgen_idx);
    CHKERRQ(ierr);
  }

  ierr = PetscFree(genparams->xidx);
  CHKERRQ(ierr);
  ierr = PetscFree(genparams->gidxbus);
  CHKERRQ(ierr);

  ierr = PetscFree(genparams->eqjacspbus_idx);
  CHKERRQ(ierr);
  ierr = PetscFree(genparams->eqjacsqbus_idx);
  CHKERRQ(ierr);
  ierr = PetscFree(genparams->hesssp_idx);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
/* Create data for generators that is used in different computations */
PetscErrorCode CreateGenParams(OPFLOW opflow, GENParams *genparams) {
  PS ps = opflow->ps;
  PetscInt loc, gloc = 0, geni = 0;
  PSGEN gen;
  PSBUS bus;
  PetscInt i, j;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  /* Get the number of active generators (STATUS ON) */
  ierr = PSGetNumActiveGenerators(ps, &genparams->ngenON, NULL);
  CHKERRQ(ierr);

  /* Allocate the arrays */
  ierr = PetscCalloc1(genparams->ngenON, &genparams->cost_alpha);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(genparams->ngenON, &genparams->cost_beta);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(genparams->ngenON, &genparams->cost_gamma);
  CHKERRQ(ierr);

  ierr = PetscCalloc1(genparams->ngenON, &genparams->pt);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(genparams->ngenON, &genparams->pb);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(genparams->ngenON, &genparams->qt);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(genparams->ngenON, &genparams->qb);
  CHKERRQ(ierr);

  ierr = PetscCalloc1(genparams->ngenON, &genparams->hesssp_idx);
  CHKERRQ(ierr);

  if (opflow->has_gensetpoint) {
    ierr = PetscCalloc1(genparams->ngenON, &genparams->geqidxgen);
    CHKERRQ(ierr);
    ierr = PetscCalloc1(genparams->ngenON, &genparams->gineqidxgen);
    CHKERRQ(ierr);
    ierr = PetscCalloc1(genparams->ngenON, &genparams->gbineqidxgen);
    CHKERRQ(ierr);
    ierr = PetscCalloc1(genparams->ngenON, &genparams->pgs);
    CHKERRQ(ierr);

    ierr = PetscCalloc1(genparams->ngenON, &genparams->eqjacspgen_idx);
    CHKERRQ(ierr);
    ierr = PetscCalloc1(genparams->ngenON, &genparams->ineqjacspgen_idx);
    CHKERRQ(ierr);
  }

  ierr = PetscCalloc1(genparams->ngenON, &genparams->xidx);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(genparams->ngenON, &genparams->gidxbus);
  CHKERRQ(ierr);

  ierr = PetscCalloc1(genparams->ngenON, &genparams->eqjacspbus_idx);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(genparams->ngenON, &genparams->eqjacsqbus_idx);
  CHKERRQ(ierr);

  /* Insert data in genparams */
  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    gloc = bus->starteqloc;

    for (j = 0; j < bus->ngen; j++) {
      ierr = PSBUSGetGen(bus, j, &gen);
      CHKERRQ(ierr);
      if (!gen->status)
        continue;

      loc = gen->startxpowloc;

      genparams->cost_alpha[geni] = gen->cost_alpha;
      genparams->cost_beta[geni] = gen->cost_beta;
      genparams->cost_gamma[geni] = gen->cost_gamma;
      genparams->pt[geni] = gen->pt;
      genparams->pb[geni] = gen->pb;
      genparams->qt[geni] = gen->qt;
      genparams->qb[geni] = gen->qb;
      if (opflow->has_gensetpoint) {
        genparams->pgs[geni] = gen->pgs;
      }

      genparams->xidx[geni] = opflow->idxn2sd_map[loc];
      genparams->gidxbus[geni] = gloc;
      if (opflow->has_gensetpoint) {
        genparams->geqidxgen[geni] = gen->starteqloc;
        genparams->gineqidxgen[geni] = gen->startineqloc;
        genparams->gbineqidxgen[geni] = opflow->nconeq + gen->startineqloc;
      }

      geni++;
    }
  }

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolutionToPS_PBPOLHIOP(OPFLOW opflow) {
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

PetscErrorCode OPFLOWModelSetUp_PBPOLHIOP(OPFLOW opflow) {
  PetscErrorCode ierr;
  PBPOLHIOP pbpolhiop = (PBPOLHIOP)opflow->model;
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

  ierr = CreateBusParams(opflow, &pbpolhiop->busparams);
  ierr = CreateGenParams(opflow, &pbpolhiop->genparams);
  ierr = CreateLineParams(opflow, &pbpolhiop->lineparams);
  ierr = CreateLoadParams(opflow, &pbpolhiop->loadparams);

  BUSParams *busparams = &pbpolhiop->busparams;
  GENParams *genparams = &pbpolhiop->genparams;
  LOADParams *loadparams = &pbpolhiop->loadparams;
  LINEParams *lineparams = &pbpolhiop->lineparams;

  int geni = 0, loadi = 0, gi;
  int nnz_eqjacsp = 0;
  int nnz_ineqjacsp = 0;
  int nnz_hesssp = 0;
  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    /* We need to go over each row since HiOp sparse Jacobian needs the values
     to be inserted row-wise (!)
     So we first go over the partial derivatives w.r.t. real power balance
     equations and then over the second one for reactive power balance
     equations.
    */

    if (opflow->include_powerimbalance_variables) {
      busparams->jacsp_idx[i] = nnz_eqjacsp;
      busparams->hesssp_idx[i] = nnz_hesssp;
      nnz_eqjacsp += 1;
      nnz_hesssp += 2; /* Includes contributions for real and imaginary part */
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
      nnz_eqjacsp += 1;
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

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWModelDestroy_PBPOLHIOP(OPFLOW opflow) {
  PetscErrorCode ierr;
  PBPOLHIOP pbpolhiop = (PBPOLHIOP)opflow->model;

  PetscFunctionBegin;

  ierr = DestroyBusParams(opflow, &pbpolhiop->busparams);
  ierr = DestroyGenParams(opflow, &pbpolhiop->genparams);
  ierr = DestroyLineParams(opflow, &pbpolhiop->lineparams);
  ierr = DestroyLoadParams(opflow, &pbpolhiop->loadparams);

  ierr = PetscFree(opflow->model);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* reuse numvariables and numconstraints functions from PBPOL model */
extern PetscErrorCode OPFLOWModelSetNumVariables_PBPOL(OPFLOW, PetscInt *,
                                                       PetscInt *, PetscInt *);
extern PetscErrorCode OPFLOWModelSetNumConstraints_PBPOL(OPFLOW, PetscInt *,
                                                         PetscInt *, PetscInt *,
                                                         PetscInt *);

PetscErrorCode OPFLOWModelCreate_PBPOLHIOP(OPFLOW opflow) {
  PBPOLHIOP pbpol;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = PetscCalloc1(1, &pbpol);
  CHKERRQ(ierr);

  opflow->model = pbpol;

  /* HIOP models only support VARIABLE_WITHIN_BOUNDS opflow->genbusvoltagetype
   */
  opflow->genbusvoltagetype = VARIABLE_WITHIN_BOUNDS;

  opflow->spdnordering = PETSC_TRUE;

  /* Inherit Ops */
  opflow->modelops.destroy = OPFLOWModelDestroy_PBPOLHIOP;
  opflow->modelops.setnumvariables = OPFLOWModelSetNumVariables_PBPOL;
  opflow->modelops.setnumconstraints = OPFLOWModelSetNumConstraints_PBPOL;
  opflow->modelops.setvariableboundsarray =
      OPFLOWSetVariableBoundsArray_PBPOLHIOP;
  opflow->modelops.setconstraintboundsarray =
      OPFLOWSetConstraintBoundsArray_PBPOLHIOP;
  opflow->modelops.setinitialguessarray = OPFLOWSetInitialGuessArray_PBPOLHIOP;
  opflow->modelops.computeequalityconstraintsarray =
      OPFLOWComputeEqualityConstraintsArray_PBPOLHIOP;
  opflow->modelops.computeinequalityconstraintsarray =
      OPFLOWComputeInequalityConstraintsArray_PBPOLHIOP;
  opflow->modelops.computeobjectivearray =
      OPFLOWComputeObjectiveArray_PBPOLHIOP;
  opflow->modelops.computegradientarray = OPFLOWComputeGradientArray_PBPOLHIOP;
  opflow->modelops.solutiontops = OPFLOWSolutionToPS_PBPOLHIOP;
  opflow->modelops.setup = OPFLOWModelSetUp_PBPOLHIOP;
  opflow->modelops.computesparseequalityconstraintjacobianhiop =
      OPFLOWComputeSparseEqualityConstraintJacobian_PBPOLHIOP;
  opflow->modelops.computesparseinequalityconstraintjacobianhiop =
      OPFLOWComputeSparseInequalityConstraintJacobian_PBPOLHIOP;
  opflow->modelops.computesparsehessianhiop =
      OPFLOWComputeSparseHessian_PBPOLHIOP;
  opflow->modelops.computedenseequalityconstraintjacobianhiop =
      OPFLOWComputeDenseEqualityConstraintJacobian_PBPOLHIOP;
  opflow->modelops.computedenseinequalityconstraintjacobianhiop =
      OPFLOWComputeDenseInequalityConstraintJacobian_PBPOLHIOP;
  opflow->modelops.computedensehessianhiop =
      OPFLOWComputeDenseHessian_PBPOLHIOP;

  PetscFunctionReturn(0);
}

#endif
