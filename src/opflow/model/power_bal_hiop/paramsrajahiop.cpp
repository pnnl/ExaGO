
#include "paramsrajahiop.h"
#include <private/opflowimpl.h>

template <typename T, typename I> T *paramAlloc(umpire::Allocator &a, I N) {
  return static_cast<T *>(a.allocate(N * sizeof(T)));
}

/* Functions to create and destroy data arrays for different
   component classes
*/
int BUSParamsRajaHiop::destroy(OPFLOW opflow) {
  h_allocator_.deallocate(isref);
  h_allocator_.deallocate(isisolated);
  h_allocator_.deallocate(ispvpq);
  h_allocator_.deallocate(vmin);
  h_allocator_.deallocate(vmax);
  h_allocator_.deallocate(va);
  h_allocator_.deallocate(vm);
  h_allocator_.deallocate(gl);
  h_allocator_.deallocate(bl);
  h_allocator_.deallocate(xidx);
  h_allocator_.deallocate(gidx);
  h_allocator_.deallocate(jacsp_idx);
  h_allocator_.deallocate(jacsq_idx);
  if (opflow->include_powerimbalance_variables) {
    h_allocator_.deallocate(xidxpimb);
    h_allocator_.deallocate(powerimbalance_penalty);
  }

#ifdef EXAGO_ENABLE_GPU
  d_allocator_.deallocate(isref_dev_);
  d_allocator_.deallocate(isisolated_dev_);
  d_allocator_.deallocate(ispvpq_dev_);
  d_allocator_.deallocate(vmin_dev_);
  d_allocator_.deallocate(vmax_dev_);
  d_allocator_.deallocate(va_dev_);
  d_allocator_.deallocate(vm_dev_);
  d_allocator_.deallocate(gl_dev_);
  d_allocator_.deallocate(bl_dev_);
  d_allocator_.deallocate(xidx_dev_);
  d_allocator_.deallocate(gidx_dev_);
  d_allocator_.deallocate(jacsp_idx_dev_);
  d_allocator_.deallocate(jacsq_idx_dev_);
  if (opflow->include_powerimbalance_variables) {
    d_allocator_.deallocate(xidxpimb_dev_);
    d_allocator_.deallocate(powerimbalance_penalty_dev_);
  }
#endif

  return 0;
}

/* Copy data from host to device */
int BUSParamsRajaHiop::copy(OPFLOW opflow) {
  /* Allocate the arrays */
  auto &resmgr = umpire::ResourceManager::getInstance();
  h_allocator_ = resmgr.getAllocator("HOST");

#ifdef EXAGO_ENABLE_GPU
  d_allocator_ = resmgr.getAllocator("DEVICE");

  // Copy host data from the host to the device
  resmgr.copy(isref_dev_, isref);
  resmgr.copy(isisolated_dev_, isisolated);
  resmgr.copy(ispvpq_dev_, ispvpq);

  resmgr.copy(vmin_dev_, vmin);
  resmgr.copy(vmax_dev_, vmax);
  resmgr.copy(va_dev_, va);
  resmgr.copy(vm_dev_, vm);
  resmgr.copy(gl_dev_, gl);
  resmgr.copy(bl_dev_, bl);

  resmgr.copy(xidx_dev_, xidx);
  resmgr.copy(gidx_dev_, gidx);
  resmgr.copy(jacsp_idx_dev_, jacsp_idx);
  resmgr.copy(jacsq_idx_dev_, jacsq_idx);
  if (opflow->include_powerimbalance_variables) {
    resmgr.copy(xidxpimb_dev_, xidxpimb);
    resmgr.copy(powerimbalance_penalty_dev_, powerimbalance_penalty);
  }
#else
  isref_dev_ = isref;
  isisolated_dev_ = isisolated;
  ispvpq_dev_ = ispvpq;
  vmin_dev_ = vmin;
  vmax_dev_ = vmax;
  va_dev_ = va;
  vm_dev_ = vm;
  gl_dev_ = gl;
  bl_dev_ = bl;
  xidx_dev_ = xidx;
  xidxpimb_dev_ = xidxpimb;
  gidx_dev_ = gidx;
  jacsp_idx_dev_ = jacsp_idx;
  jacsq_idx_dev_ = jacsq_idx;
  powerimbalance_penalty_dev_ = powerimbalance_penalty;
#endif
  return 0;
}

/* Create data for buses that is used in different computations */
int BUSParamsRajaHiop::allocate(OPFLOW opflow) {
  PS ps = opflow->ps;
  PetscInt loc;
  PSBUS bus;

  nbus = ps->nbus;
  /* Allocate the arrays */
  auto &resmgr = umpire::ResourceManager::getInstance();
  h_allocator_ = resmgr.getAllocator("HOST");

  // Allocate data on the host
  isref = paramAlloc<int>(h_allocator_, nbus);
  isisolated = paramAlloc<int>(h_allocator_, nbus);
  ispvpq = paramAlloc<int>(h_allocator_, nbus);

  vmin = paramAlloc<double>(h_allocator_, nbus);
  vmax = paramAlloc<double>(h_allocator_, nbus);
  va = paramAlloc<double>(h_allocator_, nbus);
  vm = paramAlloc<double>(h_allocator_, nbus);
  gl = paramAlloc<double>(h_allocator_, nbus);
  bl = paramAlloc<double>(h_allocator_, nbus);

  xidx = paramAlloc<int>(h_allocator_, nbus);
  gidx = paramAlloc<int>(h_allocator_, nbus);

  jacsp_idx = paramAlloc<int>(h_allocator_, nbus);
  jacsq_idx = paramAlloc<int>(h_allocator_, nbus);
  if (opflow->include_powerimbalance_variables) {
    xidxpimb = paramAlloc<int>(h_allocator_, nbus);
    powerimbalance_penalty = paramAlloc<double>(h_allocator_, nbus);
  }

  /* Memzero arrays */
  resmgr.memset(isref, 0, nbus * sizeof(int));
  resmgr.memset(ispvpq, 0, nbus * sizeof(int));
  resmgr.memset(isisolated, 0, nbus * sizeof(int));

  for (int i = 0; i < nbus; i++) {
    bus = &ps->bus[i];
    loc = bus->startxVloc;

    xidx[i] = opflow->idxn2sd_map[loc];
    gidx[i] = bus->starteqloc;

    if (bus->ide == REF_BUS)
      isref[i] = 1;
    else if (bus->ide == ISOLATED_BUS)
      isisolated[i] = 1;
    else
      ispvpq[i] = 1;

    if (opflow->genbusvoltagetype == FIXED_AT_SETPOINT) {
      if (bus->ide == REF_BUS || bus->ide == PV_BUS) {
        /* Hold voltage at reference and PV buses */
        vmin[i] = bus->vm;
        vmax[i] = bus->vm;
      } else {
        vmin[i] = bus->Vmin;
        vmax[i] = bus->Vmax;
      }
    } else {
      vmin[i] = bus->Vmin;
      vmax[i] = bus->Vmax;
    }
    vm[i] = bus->vm;
    va[i] = bus->va;
    gl[i] = bus->gl;
    bl[i] = bus->bl;

    if (opflow->include_powerimbalance_variables) {
      loc = bus->startxpimbloc;
      xidxpimb[i] = opflow->idxn2sd_map[loc];
      powerimbalance_penalty[i] = opflow->powerimbalance_penalty;
    }
  }

#ifdef EXAGO_ENABLE_GPU
  d_allocator_ = resmgr.getAllocator("DEVICE");
  // Allocate data on the device
  isref_dev_ = paramAlloc<int>(d_allocator_, nbus);
  isisolated_dev_ = paramAlloc<int>(d_allocator_, nbus);
  ispvpq_dev_ = paramAlloc<int>(d_allocator_, nbus);

  vmin_dev_ = paramAlloc<double>(d_allocator_, nbus);
  vmax_dev_ = paramAlloc<double>(d_allocator_, nbus);
  va_dev_ = paramAlloc<double>(d_allocator_, nbus);
  vm_dev_ = paramAlloc<double>(d_allocator_, nbus);
  gl_dev_ = paramAlloc<double>(d_allocator_, nbus);
  bl_dev_ = paramAlloc<double>(d_allocator_, nbus);

  xidx_dev_ = paramAlloc<int>(d_allocator_, nbus);
  gidx_dev_ = paramAlloc<int>(d_allocator_, nbus);

  jacsp_idx_dev_ = paramAlloc<int>(d_allocator_, nbus);
  jacsq_idx_dev_ = paramAlloc<int>(d_allocator_, nbus);
  if (opflow->include_powerimbalance_variables) {
    xidxpimb_dev_ = paramAlloc<int>(d_allocator_, nbus);
    powerimbalance_penalty_dev_ = paramAlloc<double>(d_allocator_, nbus);
  }
#endif
  return 0;
}

int LINEParamsRajaHiop::copy(OPFLOW opflow) {
  /* Allocate the arrays */
  auto &resmgr = umpire::ResourceManager::getInstance();
  h_allocator_ = resmgr.getAllocator("HOST");

#ifdef EXAGO_ENABLE_GPU
  // Copy data from the host to the device
  resmgr.copy(Gff_dev_, Gff);
  resmgr.copy(Bff_dev_, Bff);
  resmgr.copy(Gft_dev_, Gft);
  resmgr.copy(Bft_dev_, Bft);
  resmgr.copy(Gtf_dev_, Gtf);
  resmgr.copy(Btf_dev_, Btf);
  resmgr.copy(Gtt_dev_, Gtt);
  resmgr.copy(Btt_dev_, Btt);
  resmgr.copy(rateA_dev_, rateA);

  resmgr.copy(xidxf_dev_, xidxf);
  resmgr.copy(xidxt_dev_, xidxt);

  resmgr.copy(geqidxf_dev_, geqidxf);
  resmgr.copy(geqidxt_dev_, geqidxt);
  resmgr.copy(busf_idx_dev_, busf_idx);
  resmgr.copy(bust_idx_dev_, bust_idx);
  resmgr.copy(jacf_idx_dev_, jacf_idx);
  resmgr.copy(jact_idx_dev_, jact_idx);

  if (opflow->nlinesmon) {
    resmgr.copy(gineqidx_dev_, gineqidx);
    resmgr.copy(gbineqidx_dev_, gbineqidx);
    resmgr.copy(linelimidx_dev_, linelimidx);
    resmgr.copy(jac_ieq_idx_dev_, jac_ieq_idx);
  }
#else
  Gff_dev_ = Gff;
  Bff_dev_ = Bff;
  Gft_dev_ = Gft;
  Bft_dev_ = Bft;
  Gtf_dev_ = Gtf;
  Btf_dev_ = Btf;
  Gtt_dev_ = Gtt;
  Btt_dev_ = Btt;
  rateA_dev_ = rateA;
  xidxf_dev_ = xidxf;
  xidxt_dev_ = xidxt;
  geqidxf_dev_ = geqidxf;
  geqidxt_dev_ = geqidxt;
  busf_idx_dev_ = busf_idx;
  bust_idx_dev_ = bust_idx;
  jacf_idx_dev_ = jacf_idx;
  jact_idx_dev_ = jact_idx;
  if (opflow->nlinesmon) {
    gineqidx_dev_ = gineqidx;
    gbineqidx_dev_ = gbineqidx;
    linelimidx_dev_ = linelimidx;
    jac_ieq_idx_dev_ = jac_ieq_idx;
  }
#endif
  return 0;
}

int LINEParamsRajaHiop::destroy(OPFLOW opflow) {
  // Destroy parameter arrays on the host
  h_allocator_.deallocate(Gff);
  h_allocator_.deallocate(Bff);
  h_allocator_.deallocate(Gft);
  h_allocator_.deallocate(Bft);
  h_allocator_.deallocate(Gtf);
  h_allocator_.deallocate(Btf);
  h_allocator_.deallocate(Gtt);
  h_allocator_.deallocate(Btt);
  h_allocator_.deallocate(rateA);

  h_allocator_.deallocate(xidxf);
  h_allocator_.deallocate(xidxt);

  h_allocator_.deallocate(geqidxf);
  h_allocator_.deallocate(geqidxt);
  h_allocator_.deallocate(busf_idx);
  h_allocator_.deallocate(bust_idx);
  h_allocator_.deallocate(jacf_idx);
  h_allocator_.deallocate(jact_idx);

  if (opflow->nlinesmon) {
    h_allocator_.deallocate(gineqidx);
    h_allocator_.deallocate(gbineqidx);
    h_allocator_.deallocate(linelimidx);
    h_allocator_.deallocate(jac_ieq_idx);
  }

#ifdef EXAGO_ENABLE_GPU
  // Destroy parameter arrays on the device
  d_allocator_.deallocate(Gff_dev_);
  d_allocator_.deallocate(Bff_dev_);
  d_allocator_.deallocate(Gft_dev_);
  d_allocator_.deallocate(Bft_dev_);
  d_allocator_.deallocate(Gtf_dev_);
  d_allocator_.deallocate(Btf_dev_);
  d_allocator_.deallocate(Gtt_dev_);
  d_allocator_.deallocate(Btt_dev_);
  d_allocator_.deallocate(rateA_dev_);

  d_allocator_.deallocate(xidxf_dev_);
  d_allocator_.deallocate(xidxt_dev_);

  d_allocator_.deallocate(geqidxf_dev_);
  d_allocator_.deallocate(geqidxt_dev_);
  d_allocator_.deallocate(busf_idx_dev_);
  d_allocator_.deallocate(bust_idx_dev_);
  d_allocator_.deallocate(jacf_idx_dev_);
  d_allocator_.deallocate(jact_idx_dev_);

  if (opflow->nlinesmon) {
    d_allocator_.deallocate(gineqidx_dev_);
    d_allocator_.deallocate(gbineqidx_dev_);
    d_allocator_.deallocate(linelimidx_dev_);
    d_allocator_.deallocate(jac_ieq_idx_dev_);
  }
#endif

  return 0;
}

/* Create data for lines that is used in different computations */
int LINEParamsRajaHiop::allocate(OPFLOW opflow) {
  PS ps = opflow->ps;
  PetscInt linei = 0;
  PSLINE line;
  PetscInt i;
  PetscErrorCode ierr;
  const PSBUS *connbuses;
  PSBUS busf, bust;

  PetscFunctionBegin;
  ierr = PSGetNumActiveLines(ps, &nlineON, NULL);
  CHKERRQ(ierr);

  nlinelim = opflow->nlinesmon;

  /* Allocate data arrays */
  auto &resmgr = umpire::ResourceManager::getInstance();
  h_allocator_ = resmgr.getAllocator("HOST");

  // Allocate data on the host
  Gff = paramAlloc<double>(h_allocator_, nlineON);
  Bff = paramAlloc<double>(h_allocator_, nlineON);
  Gft = paramAlloc<double>(h_allocator_, nlineON);
  Bft = paramAlloc<double>(h_allocator_, nlineON);
  Gtf = paramAlloc<double>(h_allocator_, nlineON);
  Btf = paramAlloc<double>(h_allocator_, nlineON);
  Gtt = paramAlloc<double>(h_allocator_, nlineON);
  Btt = paramAlloc<double>(h_allocator_, nlineON);
  rateA = paramAlloc<double>(h_allocator_, nlineON);

  xidxf = paramAlloc<int>(h_allocator_, nlineON);
  xidxt = paramAlloc<int>(h_allocator_, nlineON);

  geqidxf = paramAlloc<int>(h_allocator_, nlineON);
  geqidxt = paramAlloc<int>(h_allocator_, nlineON);

  busf_idx = paramAlloc<int>(h_allocator_, nlineON);
  bust_idx = paramAlloc<int>(h_allocator_, nlineON);
  jacf_idx = paramAlloc<int>(h_allocator_, nlineON);
  jact_idx = paramAlloc<int>(h_allocator_, nlineON);

  if (opflow->nlinesmon) {
    linelimidx = paramAlloc<int>(h_allocator_, nlinelim);
    gineqidx = paramAlloc<int>(h_allocator_, nlinelim);
    gbineqidx = paramAlloc<int>(h_allocator_, nlinelim);
    jac_ieq_idx = paramAlloc<int>(h_allocator_, nlinelim);
  }

  PetscInt j = 0;
  /* Populate arrays */
  for (i = 0; i < ps->nline; i++) {
    line = &ps->line[i];

    if (!line->status)
      continue;

    Gff[linei] = line->yff[0];
    Bff[linei] = line->yff[1];
    Gft[linei] = line->yft[0];
    Bft[linei] = line->yft[1];
    Gtf[linei] = line->ytf[0];
    Btf[linei] = line->ytf[1];
    Gtt[linei] = line->ytt[0];
    Btt[linei] = line->ytt[1];
    rateA[linei] = line->rateA;

    ierr = PSLINEGetConnectedBuses(line, &connbuses);
    CHKERRQ(ierr);
    busf = connbuses[0];
    bust = connbuses[1];

    int xidxfi, xidxti;
    xidxfi = busf->startxVloc;
    xidxti = bust->startxVloc;

    xidxf[linei] = opflow->idxn2sd_map[xidxfi];
    xidxt[linei] = opflow->idxn2sd_map[xidxti];

    /*
       Each bus has two equality (balance) constraints, hence the use of
       coefficient 2 to map the location of the equality constraint for the bus
    */
    geqidxf[linei] = busf->starteqloc;
    geqidxt[linei] = bust->starteqloc;
    busf_idx[linei] = ps->busext2intmap[line->fbus];
    bust_idx[linei] = ps->busext2intmap[line->tbus];
    jacf_idx[linei] = 0;
    jact_idx[linei] = 0;

    if (j < opflow->nlinesmon && opflow->linesmon[j] == i) {
      gbineqidx[j] = opflow->nconeq + line->startineqloc;
      gineqidx[j] = line->startineqloc;
      linelimidx[j] = linei;
      jac_ieq_idx[j] = 0;
      j++;
    }

    linei++;
  }

#ifdef EXAGO_ENABLE_GPU
  d_allocator_ = resmgr.getAllocator("DEVICE");
  // Allocate data on the device
  Gff_dev_ = paramAlloc<double>(d_allocator_, nlineON);
  Bff_dev_ = paramAlloc<double>(d_allocator_, nlineON);
  Gft_dev_ = paramAlloc<double>(d_allocator_, nlineON);
  Bft_dev_ = paramAlloc<double>(d_allocator_, nlineON);
  Gtf_dev_ = paramAlloc<double>(d_allocator_, nlineON);
  Btf_dev_ = paramAlloc<double>(d_allocator_, nlineON);
  Gtt_dev_ = paramAlloc<double>(d_allocator_, nlineON);
  Btt_dev_ = paramAlloc<double>(d_allocator_, nlineON);
  rateA_dev_ = paramAlloc<double>(d_allocator_, nlineON);

  xidxf_dev_ = paramAlloc<int>(d_allocator_, nlineON);
  xidxt_dev_ = paramAlloc<int>(d_allocator_, nlineON);

  geqidxf_dev_ = paramAlloc<int>(d_allocator_, nlineON);
  geqidxt_dev_ = paramAlloc<int>(d_allocator_, nlineON);

  busf_idx_dev_ = paramAlloc<int>(d_allocator_, nlineON);
  bust_idx_dev_ = paramAlloc<int>(d_allocator_, nlineON);
  jacf_idx_dev_ = paramAlloc<int>(d_allocator_, nlineON);
  jact_idx_dev_ = paramAlloc<int>(d_allocator_, nlineON);

  if (opflow->nconineq) {
    gineqidx_dev_ = paramAlloc<int>(d_allocator_, nlinelim);
    gbineqidx_dev_ = paramAlloc<int>(d_allocator_, nlinelim);
    linelimidx_dev_ = paramAlloc<int>(d_allocator_, nlinelim);
    jac_ieq_idx_dev_ = paramAlloc<int>(d_allocator_, nlinelim);
  }
#endif
  return 0;
}

int LOADParamsRajaHiop::copy(OPFLOW opflow) {
  /* Allocate the arrays */
  auto &resmgr = umpire::ResourceManager::getInstance();
  h_allocator_ = resmgr.getAllocator("HOST");

#ifdef EXAGO_ENABLE_GPU
  // Copy data from host to device
  resmgr.copy(pl_dev_, pl);
  resmgr.copy(ql_dev_, ql);
  resmgr.copy(xidx_dev_, xidx);
  resmgr.copy(gidx_dev_, gidx);
  if (opflow->include_loadloss_variables) {
    resmgr.copy(jacsp_idx_dev_, jacsp_idx);
    resmgr.copy(jacsq_idx_dev_, jacsq_idx);
    resmgr.copy(hesssp_idx_dev_, hesssp_idx);
    resmgr.copy(loadloss_penalty_dev_, loadloss_penalty);
  }
#else
  pl_dev_ = pl;
  ql_dev_ = ql;
  xidx_dev_ = xidx;
  gidx_dev_ = gidx;
  jacsp_idx_dev_ = jacsp_idx;
  jacsq_idx_dev_ = jacsq_idx;
  hesssp_idx_dev_ = hesssp_idx;
  loadloss_penalty_dev_ = loadloss_penalty;
#endif

  return 0;
}

int LOADParamsRajaHiop::destroy(OPFLOW opflow) {
  h_allocator_.deallocate(pl);
  h_allocator_.deallocate(ql);
  h_allocator_.deallocate(xidx);
  h_allocator_.deallocate(gidx);
  if (opflow->include_loadloss_variables) {
    h_allocator_.deallocate(loadloss_penalty);
    h_allocator_.deallocate(jacsp_idx);
    h_allocator_.deallocate(jacsq_idx);
    h_allocator_.deallocate(hesssp_idx);
  }
#ifdef EXAGO_ENABLE_GPU
  d_allocator_.deallocate(pl_dev_);
  d_allocator_.deallocate(ql_dev_);
  d_allocator_.deallocate(xidx_dev_);
  d_allocator_.deallocate(gidx_dev_);
  if (opflow->include_loadloss_variables) {
    d_allocator_.deallocate(loadloss_penalty_dev_);
    d_allocator_.deallocate(jacsp_idx_dev_);
    d_allocator_.deallocate(jacsq_idx_dev_);
    d_allocator_.deallocate(hesssp_idx_dev_);
  }
#endif
  return 0;
}

/* Create data for loads that is used in different computations */
int LOADParamsRajaHiop::allocate(OPFLOW opflow) {
  PS ps = opflow->ps;
  PetscInt loc, loadi = 0;
  PSLOAD load;
  PSBUS bus;
  PetscInt i, j;
  PetscErrorCode ierr;

  /* Get the number of active generators (STATUS ON) */
  ierr = PSGetNumLoads(ps, &nload, NULL);
  CHKERRQ(ierr);

  /* Allocate arrays */
  auto &resmgr = umpire::ResourceManager::getInstance();
  h_allocator_ = resmgr.getAllocator("HOST");

  // Allocate data on the host
  pl = paramAlloc<double>(h_allocator_, nload);
  ql = paramAlloc<double>(h_allocator_, nload);
  xidx = paramAlloc<int>(h_allocator_, nload);
  gidx = paramAlloc<int>(h_allocator_, nload);
  if (opflow->include_loadloss_variables) {
    loadloss_penalty = paramAlloc<double>(h_allocator_, nload);
    jacsp_idx = paramAlloc<int>(h_allocator_, nload);
    jacsq_idx = paramAlloc<int>(h_allocator_, nload);
    hesssp_idx = paramAlloc<int>(h_allocator_, nload);
  }
  /* Insert data in loadparams */
  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSGetVariableLocation(bus, &loc);
    CHKERRQ(ierr);
    for (j = 0; j < bus->nload; j++) {
      ierr = PSBUSGetLoad(bus, j, &load);
      CHKERRQ(ierr);
      pl[loadi] = load->pl;
      ql[loadi] = load->ql;
      if (opflow->include_loadloss_variables) {
        loc = load->startxloadlossloc;
        loadloss_penalty[loadi] = load->loss_cost;
        xidx[loadi] = opflow->idxn2sd_map[loc];
      }
      gidx[loadi] = bus->starteqloc;
      loadi++;
    }
  }
#ifdef EXAGO_ENABLE_GPU
  d_allocator_ = resmgr.getAllocator("DEVICE");
  // Allocate data on the device
  pl_dev_ = paramAlloc<double>(d_allocator_, nload);
  ql_dev_ = paramAlloc<double>(d_allocator_, nload);
  xidx_dev_ = paramAlloc<int>(d_allocator_, nload);
  gidx_dev_ = paramAlloc<int>(d_allocator_, nload);
  if (opflow->include_loadloss_variables) {
    loadloss_penalty_dev_ = paramAlloc<double>(d_allocator_, nload);
    jacsp_idx_dev_ = paramAlloc<int>(d_allocator_, nload);
    jacsq_idx_dev_ = paramAlloc<int>(d_allocator_, nload);
    hesssp_idx_dev_ = paramAlloc<int>(d_allocator_, nload);
  }

#endif

  return (0);
}

int GENParamsRajaHiop::destroy(OPFLOW opflow) {
  // Free arrays on the host
  h_allocator_.deallocate(cost_alpha);
  h_allocator_.deallocate(cost_beta);
  h_allocator_.deallocate(cost_gamma);
  h_allocator_.deallocate(pt);
  h_allocator_.deallocate(pb);
  h_allocator_.deallocate(qt);
  h_allocator_.deallocate(qb);
  h_allocator_.deallocate(isrenewable);
  h_allocator_.deallocate(xidx);
  h_allocator_.deallocate(gidxbus);
  h_allocator_.deallocate(eqjacspbus_idx);
  h_allocator_.deallocate(eqjacsqbus_idx);
  h_allocator_.deallocate(hesssp_idx);
  if (opflow->has_gensetpoint) {
    h_allocator_.deallocate(geqidxgen);
    h_allocator_.deallocate(gineqidxgen);
    h_allocator_.deallocate(gbineqidxgen);
    h_allocator_.deallocate(pgs);
    h_allocator_.deallocate(eqjacspgen_idx);
    h_allocator_.deallocate(ineqjacspgen_idx);
  }
#ifdef EXAGO_ENABLE_GPU
  // Free arrays on the device
  d_allocator_.deallocate(cost_alpha_dev_);
  d_allocator_.deallocate(cost_beta_dev_);
  d_allocator_.deallocate(cost_gamma_dev_);
  d_allocator_.deallocate(pt_dev_);
  d_allocator_.deallocate(pb_dev_);
  d_allocator_.deallocate(qt_dev_);
  d_allocator_.deallocate(qb_dev_);
  d_allocator_.deallocate(isrenewable_dev_);
  d_allocator_.deallocate(xidx_dev_);
  d_allocator_.deallocate(gidxbus_dev_);
  d_allocator_.deallocate(eqjacspbus_idx_dev_);
  d_allocator_.deallocate(eqjacsqbus_idx_dev_);
  d_allocator_.deallocate(hesssp_idx_dev_);
  if (opflow->has_gensetpoint) {
    d_allocator_.deallocate(geqidxgen_dev_);
    d_allocator_.deallocate(gineqidxgen_dev_);
    d_allocator_.deallocate(gbineqidxgen_dev_);
    d_allocator_.deallocate(pgs_dev_);
    d_allocator_.deallocate(eqjacspgen_idx_dev_);
    d_allocator_.deallocate(ineqjacspgen_idx_dev_);
  }

#endif
  return 0;
}

int GENParamsRajaHiop::copy(OPFLOW opflow) {
  /* Allocate arrays on the host */
  auto &resmgr = umpire::ResourceManager::getInstance();
  h_allocator_ = resmgr.getAllocator("HOST");

#ifdef EXAGO_ENABLE_GPU
  // Copy host data to the device
  resmgr.copy(cost_alpha_dev_, cost_alpha);
  resmgr.copy(cost_beta_dev_, cost_beta);
  resmgr.copy(cost_gamma_dev_, cost_gamma);

  resmgr.copy(pt_dev_, pt);
  resmgr.copy(pb_dev_, pb);
  resmgr.copy(qt_dev_, qt);
  resmgr.copy(qb_dev_, qb);
  resmgr.copy(isrenewable_dev_, isrenewable);

  resmgr.copy(xidx_dev_, xidx);
  resmgr.copy(gidxbus_dev_, gidxbus);

  resmgr.copy(eqjacspbus_idx_dev_, eqjacspbus_idx);
  resmgr.copy(eqjacsqbus_idx_dev_, eqjacsqbus_idx);
  resmgr.copy(hesssp_idx_dev_, hesssp_idx);
  if (opflow->has_gensetpoint) {
    resmgr.copy(geqidxgen_dev_, geqidxgen);
    resmgr.copy(gineqidxgen_dev_, gineqidxgen);
    resmgr.copy(gbineqidxgen_dev_, gbineqidxgen);
    resmgr.copy(eqjacspgen_idx_dev_, eqjacspgen_idx);
    resmgr.copy(ineqjacspgen_idx_dev_, ineqjacspgen_idx);
    resmgr.copy(pgs_dev_, pgs);
  }
#else
  cost_alpha_dev_ = cost_alpha;
  cost_beta_dev_ = cost_beta;
  cost_gamma_dev_ = cost_gamma;
  pt_dev_ = pt;
  pb_dev_ = pb;
  qt_dev_ = qt;
  qb_dev_ = qb;
  isrenewable_dev_ = isrenewable;
  xidx_dev_ = xidx;
  gidxbus_dev_ = gidxbus;
  eqjacspbus_idx_dev_ = eqjacspbus_idx;
  eqjacsqbus_idx_dev_ = eqjacsqbus_idx;
  hesssp_idx_dev_ = hesssp_idx;
  geqidxgen_dev_ = geqidxgen;
  gineqidxgen_dev_ = gineqidxgen;
  gbineqidxgen_dev_ = gbineqidxgen;
  eqjacspgen_idx_dev_ = eqjacspgen_idx;
  ineqjacspgen_idx_dev_ = ineqjacspgen_idx;
  pgs_dev_ = pgs;
#endif
  return 0;
}

/* Create data for generators that is used in different computations */
int GENParamsRajaHiop::allocate(OPFLOW opflow) {
  PS ps = opflow->ps;
  PetscInt loc, gloc = 0, geni = 0;
  PSGEN gen;
  PSBUS bus;
  PetscInt i, j;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  /* Get the number of active generators (STATUS ON) */
  ierr = PSGetNumActiveGenerators(ps, &ngenON, NULL);
  CHKERRQ(ierr);

  /* Allocate arrays on the host */
  auto &resmgr = umpire::ResourceManager::getInstance();
  h_allocator_ = resmgr.getAllocator("HOST");

  cost_alpha = paramAlloc<double>(h_allocator_, ngenON);
  cost_beta = paramAlloc<double>(h_allocator_, ngenON);
  cost_gamma = paramAlloc<double>(h_allocator_, ngenON);

  pt = paramAlloc<double>(h_allocator_, ngenON);
  pb = paramAlloc<double>(h_allocator_, ngenON);
  qt = paramAlloc<double>(h_allocator_, ngenON);
  qb = paramAlloc<double>(h_allocator_, ngenON);
  isrenewable = paramAlloc<int>(h_allocator_, ngenON);

  xidx = paramAlloc<int>(h_allocator_, ngenON);
  gidxbus = paramAlloc<int>(h_allocator_, ngenON);

  eqjacspbus_idx = paramAlloc<int>(h_allocator_, ngenON);
  eqjacsqbus_idx = paramAlloc<int>(h_allocator_, ngenON);
  hesssp_idx = paramAlloc<int>(h_allocator_, ngenON);

  if (opflow->has_gensetpoint) {
    geqidxgen = paramAlloc<int>(h_allocator_, ngenON);
    gineqidxgen = paramAlloc<int>(h_allocator_, ngenON);
    gbineqidxgen = paramAlloc<int>(h_allocator_, ngenON);
    eqjacspgen_idx = paramAlloc<int>(h_allocator_, ngenON);
    ineqjacspgen_idx = paramAlloc<int>(h_allocator_, ngenON);
    pgs = paramAlloc<double>(h_allocator_, ngenON);
  }

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

      cost_alpha[geni] = gen->cost_alpha;
      cost_beta[geni] = gen->cost_beta;
      cost_gamma[geni] = gen->cost_gamma;
      pt[geni] = gen->pt;
      pb[geni] = gen->pb;
      qt[geni] = gen->qt;
      qb[geni] = gen->qb;
      isrenewable[geni] = (int)gen->isrenewable;
      if (opflow->has_gensetpoint) {
        pgs[geni] = gen->pgs;
      }

      xidx[geni] = opflow->idxn2sd_map[loc];
      gidxbus[geni] = gloc;
      if (opflow->has_gensetpoint) {
        geqidxgen[geni] = gen->starteqloc;
        gineqidxgen[geni] = gen->startineqloc;
        gbineqidxgen[geni] = opflow->nconeq + gen->startineqloc;
      }

      geni++;
    }
  }

#ifdef EXAGO_ENABLE_GPU
  d_allocator_ = resmgr.getAllocator("DEVICE");
  // Allocate arrays on the device
  cost_alpha_dev_ = paramAlloc<double>(d_allocator_, ngenON);
  cost_beta_dev_ = paramAlloc<double>(d_allocator_, ngenON);
  cost_gamma_dev_ = paramAlloc<double>(d_allocator_, ngenON);

  pt_dev_ = paramAlloc<double>(d_allocator_, ngenON);
  pb_dev_ = paramAlloc<double>(d_allocator_, ngenON);
  qt_dev_ = paramAlloc<double>(d_allocator_, ngenON);
  qb_dev_ = paramAlloc<double>(d_allocator_, ngenON);
  isrenewable_dev_ = paramAlloc<int>(d_allocator_, ngenON);

  xidx_dev_ = paramAlloc<int>(d_allocator_, ngenON);
  gidxbus_dev_ = paramAlloc<int>(d_allocator_, ngenON);

  eqjacspbus_idx_dev_ = paramAlloc<int>(d_allocator_, ngenON);
  eqjacsqbus_idx_dev_ = paramAlloc<int>(d_allocator_, ngenON);
  hesssp_idx_dev_ = paramAlloc<int>(d_allocator_, ngenON);
  if (opflow->has_gensetpoint) {
    geqidxgen_dev_ = paramAlloc<int>(d_allocator_, ngenON);
    gineqidxgen_dev_ = paramAlloc<int>(d_allocator_, ngenON);
    gbineqidxgen_dev_ = paramAlloc<int>(d_allocator_, ngenON);
    eqjacspgen_idx_dev_ = paramAlloc<int>(d_allocator_, ngenON);
    ineqjacspgen_idx_dev_ = paramAlloc<int>(d_allocator_, ngenON);
    pgs_dev_ = paramAlloc<double>(d_allocator_, ngenON);
  }
#endif
  return 0;
}

void PbpolModelRajaHiop::destroy(OPFLOW opflow) {
  loadparams.destroy(opflow);
  lineparams.destroy(opflow);
  busparams.destroy(opflow);
  genparams.destroy(opflow);

#ifdef EXAGO_ENABLE_GPU
  if (i_jaceq != NULL) {
    // These arrays get allocated only with sparse GPU model. For other models,
    // they are not allocated. Hence, this business of checking if the array is
    // NULL

    auto &resmgr = umpire::ResourceManager::getInstance();
    umpire::Allocator h_allocator_ = resmgr.getAllocator("HOST");
    umpire::Allocator d_allocator_ = resmgr.getAllocator("DEVICE");

    h_allocator_.deallocate(i_jaceq);
    h_allocator_.deallocate(j_jaceq);
    h_allocator_.deallocate(val_jaceq);
    d_allocator_.deallocate(idx_jaceq_dev_);

    h_allocator_.deallocate(i_hess);
    h_allocator_.deallocate(j_hess);
    h_allocator_.deallocate(val_hess);

    if (opflow->nconineq) {
      h_allocator_.deallocate(i_jacineq);
      h_allocator_.deallocate(j_jacineq);
      h_allocator_.deallocate(val_jacineq);
      d_allocator_.deallocate(idx_jacineq_dev_);
    }
  }
#endif
}
