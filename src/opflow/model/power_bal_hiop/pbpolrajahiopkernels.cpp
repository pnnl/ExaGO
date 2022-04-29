#include <exago_config.h>

#if defined(EXAGO_ENABLE_RAJA)

#include <umpire/Allocator.hpp>
#include <umpire/ResourceManager.hpp>

#include <RAJA/RAJA.hpp>
#include <private/raja_exec_config.h>

#include "pbpolrajahiop.h"
#include "pbpolrajahiopkernels.h"
#include <private/opflowimpl.h>
#include <private/psimpl.h>

/************* NOTE ***********************/
/* No Load loss or power imbalance variables considered yet */
/********************************************/

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
  if (opflow->include_powerimbalance_variables) {
    h_allocator_.deallocate(xidxpimb);
    h_allocator_.deallocate(powerimbalance_penalty);
    h_allocator_.deallocate(jacsp_idx);
    h_allocator_.deallocate(jacsq_idx);
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
  if (opflow->include_powerimbalance_variables) {
    d_allocator_.deallocate(xidxpimb_dev_);
    d_allocator_.deallocate(powerimbalance_penalty_dev_);
    d_allocator_.deallocate(jacsp_idx_dev_);
    d_allocator_.deallocate(jacsq_idx_dev_);
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
  if (opflow->include_powerimbalance_variables) {
    resmgr.copy(xidxpimb_dev_, xidxpimb);
    resmgr.copy(jacsp_idx_dev_, jacsp_idx);
    resmgr.copy(jacsq_idx_dev_, jacsq_idx);
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

  if (opflow->include_powerimbalance_variables) {
    xidxpimb = paramAlloc<int>(h_allocator_, nbus);
    powerimbalance_penalty = paramAlloc<double>(h_allocator_, nbus);
    jacsp_idx = paramAlloc<int>(h_allocator_, nbus);
    jacsq_idx = paramAlloc<int>(h_allocator_, nbus);
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

  if (opflow->include_powerimbalance_variables) {
    xidxpimb_dev_ = paramAlloc<int>(d_allocator_, nbus);
    powerimbalance_penalty_dev_ = paramAlloc<double>(d_allocator_, nbus);
    jacsp_idx_dev_ = paramAlloc<int>(d_allocator_, nbus);
    jacsq_idx_dev_ = paramAlloc<int>(d_allocator_, nbus);
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

  if (opflow->nlinesmon) {
    resmgr.copy(gineqidx_dev_, gineqidx);
    resmgr.copy(gbineqidx_dev_, gbineqidx);
    resmgr.copy(linelimidx_dev_, linelimidx);
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
  if (opflow->nlinesmon) {
    gineqidx_dev_ = gineqidx;
    gbineqidx_dev_ = gbineqidx;
    linelimidx_dev_ = linelimidx;
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

  if (opflow->nlinesmon) {
    h_allocator_.deallocate(gineqidx);
    h_allocator_.deallocate(gbineqidx);
    h_allocator_.deallocate(linelimidx);
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

  if (opflow->nlinesmon) {
    d_allocator_.deallocate(gineqidx_dev_);
    d_allocator_.deallocate(gbineqidx_dev_);
    d_allocator_.deallocate(linelimidx_dev_);
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

  if (opflow->nlinesmon) {
    linelimidx = paramAlloc<int>(h_allocator_, nlinelim);
    gineqidx = paramAlloc<int>(h_allocator_, nlinelim);
    gbineqidx = paramAlloc<int>(h_allocator_, nlinelim);
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

    if (j < opflow->nlinesmon && opflow->linesmon[j] == i) {
      gbineqidx[j] = opflow->nconeq + line->startineqloc;
      gineqidx[j] = line->startineqloc;
      linelimidx[j] = linei;
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

  if (opflow->nconineq) {
    gineqidx_dev_ = paramAlloc<int>(d_allocator_, nlinelim);
    gbineqidx_dev_ = paramAlloc<int>(d_allocator_, nlinelim);
    linelimidx_dev_ = paramAlloc<int>(d_allocator_, nlinelim);
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
        loadloss_penalty[loadi] = opflow->loadloss_penalty;
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

PetscErrorCode OPFLOWSetInitialGuessArray_PBPOLRAJAHIOP(OPFLOW opflow,
                                                        double *x0_dev) {
  PetscErrorCode ierr;
  double *x;
  auto &resmgr = umpire::ResourceManager::getInstance();

  PetscFunctionBegin;
  //  ierr = PetscPrintf(MPI_COMM_SELF,"Entered
  //  OPFLOWInitialization\n");CHKERRQ(ierr);
  // Do initialization on Host
  ierr = (*opflow->modelops.setinitialguess)(opflow, opflow->X, opflow->Lambda);
  CHKERRQ(ierr);
  ierr = VecGetArray(opflow->X, &x);
  CHKERRQ(ierr);

  // Copy from host to device
  umpire::Allocator h_allocator_ = resmgr.getAllocator("HOST");
  registerWith(x, opflow->nx, resmgr, h_allocator_);
  resmgr.copy(x0_dev, x);

  ierr = VecRestoreArray(opflow->X, &x);
  CHKERRQ(ierr);

  //  ierr = PetscPrintf(MPI_COMM_SELF,"Exit Initialization\n");CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSetConstraintBoundsArray_PBPOLRAJAHIOP(OPFLOW opflow,
                                                            double *gl_dev,
                                                            double *gu_dev) {
  PbpolModelRajaHiop *pbpolrajahiop =
      reinterpret_cast<PbpolModelRajaHiop *>(opflow->model);
  LINEParamsRajaHiop *lineparams = &pbpolrajahiop->lineparams;
  PS ps = opflow->ps;
  double MVAbase = ps->MVAbase;

  PetscFunctionBegin;

  //  PetscPrintf(MPI_COMM_SELF,"Entered Constraint Bounds\n");

  /* Equality constraints (all zeros) */
  auto &resmgr = umpire::ResourceManager::getInstance();
  resmgr.memset(gl_dev, 0);
  resmgr.memset(gu_dev, 0);

  /* Inequality constraint bounds */
  int *linelimidx = lineparams->linelimidx_dev_;
  int *gbineqidx = lineparams->gbineqidx_dev_;
  double *rateA = lineparams->rateA_dev_;

  if (lineparams->nlinelim) {
    // PetscPrintf(PETSC_COMM_SELF,"nlinelim = %d ineq =
    // %d\n",lineparams->nlinelim,opflow->nconineq);
    RAJA::forall<exago_raja_exec>(
        RAJA::RangeSegment(0, lineparams->nlinelim),
        RAJA_LAMBDA(RAJA::Index_type i) {
          int j = linelimidx[i];
          gl_dev[gbineqidx[i]] = 0.0;
          gu_dev[gbineqidx[i]] = (rateA[j] / MVAbase) * (rateA[j] / MVAbase);
          gl_dev[gbineqidx[i] + 1] = 0.0;
          gu_dev[gbineqidx[i] + 1] =
              (rateA[j] / MVAbase) * (rateA[j] / MVAbase);
        });
  }

  //  PetscPrintf(MPI_COMM_SELF,"Exit Constraint Bounds\n");
  PetscFunctionReturn(0);
}

/** EQUALITY CONSTRAINTS */
PetscErrorCode OPFLOWComputeEqualityConstraintsArray_PBPOLRAJAHIOP(
    OPFLOW opflow, const double *x_dev, double *ge_dev) {
  PbpolModelRajaHiop *pbpolrajahiop =
      reinterpret_cast<PbpolModelRajaHiop *>(opflow->model);
  BUSParamsRajaHiop *busparams = &pbpolrajahiop->busparams;
  GENParamsRajaHiop *genparams = &pbpolrajahiop->genparams;
  LOADParamsRajaHiop *loadparams = &pbpolrajahiop->loadparams;
  LINEParamsRajaHiop *lineparams = &pbpolrajahiop->lineparams;
  PetscErrorCode ierr;
  PetscInt flps = 0;
  int include_loadloss_variables = opflow->include_loadloss_variables;
  int include_powerimbalance_variables =
      opflow->include_powerimbalance_variables;

  PetscFunctionBegin;
  //  PetscPrintf(MPI_COMM_SELF,"Entered Equality constraints\n");

  // Zero out array
  auto &resmgr = umpire::ResourceManager::getInstance();
  resmgr.memset(ge_dev, 0);

  /* Generator contributions */
  int *g_gidxbus = genparams->gidxbus_dev_;
  int *g_xidx = genparams->xidx_dev_;
  RAJA::forall<exago_raja_exec>(
      RAJA::RangeSegment(0, genparams->ngenON),
      RAJA_LAMBDA(RAJA::Index_type i) {
        RAJA::atomicSub<exago_raja_atomic>(&ge_dev[g_gidxbus[i]],
                                           x_dev[g_xidx[i]]);
        RAJA::atomicSub<exago_raja_atomic>(&ge_dev[g_gidxbus[i] + 1],
                                           x_dev[g_xidx[i] + 1]);
      });
  flps += genparams->ngenON * 2;

  if (opflow->has_gensetpoint) {
    int *g_geqidxgen = genparams->geqidxgen_dev_;
    double *g_pgs = genparams->pgs_dev_;
    RAJA::forall<exago_raja_exec>(
        RAJA::RangeSegment(0, genparams->ngenON),
        RAJA_LAMBDA(RAJA::Index_type i) {
          double Pg, delPg, Pgset;
          Pg = x_dev[g_xidx[i]];
          delPg = x_dev[g_xidx[i] + 2];
          Pgset = x_dev[g_xidx[i] + 3];

          ge_dev[g_geqidxgen[i]] = Pgset + delPg - Pg;
          ge_dev[g_geqidxgen[i] + 1] = Pgset - g_pgs[i];
        });
  }

  /* Load contributions */
  double *pl = loadparams->pl_dev_;
  double *ql = loadparams->ql_dev_;
  int *l_gidx = loadparams->gidx_dev_;
  int *l_xidx = loadparams->xidx_dev_;
  RAJA::forall<exago_raja_exec>(
      RAJA::RangeSegment(0, loadparams->nload),
      RAJA_LAMBDA(RAJA::Index_type i) {
        if (include_loadloss_variables) {
          RAJA::atomicAdd<exago_raja_atomic>(&ge_dev[l_gidx[i]], pl[i]);
          RAJA::atomicAdd<exago_raja_atomic>(&ge_dev[l_gidx[i] + 1], ql[i]);
          RAJA::atomicSub<exago_raja_atomic>(&ge_dev[l_gidx[i]],
                                             x_dev[l_xidx[i]]);
          RAJA::atomicSub<exago_raja_atomic>(&ge_dev[l_gidx[i] + 1],
                                             x_dev[l_xidx[i] + 1]);
        } else {
          RAJA::atomicAdd<exago_raja_atomic>(&ge_dev[l_gidx[i]], pl[i]);
          RAJA::atomicAdd<exago_raja_atomic>(&ge_dev[l_gidx[i] + 1], ql[i]);
        }
      });
  flps += loadparams->nload * 2;

  /* Bus contributions */
  int *isisolated = busparams->isisolated_dev_;
  int *ispvpq = busparams->ispvpq_dev_;
  double *gl = busparams->gl_dev_;
  double *bl = busparams->bl_dev_;
  double *vm = busparams->vm_dev_;
  double *va = busparams->va_dev_;
  int *b_xidx = busparams->xidx_dev_;
  int *b_gidx = busparams->gidx_dev_;
  int *b_xidxpimb = busparams->xidxpimb_dev_;

  RAJA::forall<exago_raja_exec>(
      RAJA::RangeSegment(0, busparams->nbus), RAJA_LAMBDA(RAJA::Index_type i) {
        double theta = x_dev[b_xidx[i]];
        double Vm = x_dev[b_xidx[i] + 1];
        RAJA::atomicAdd<exago_raja_atomic>(&ge_dev[b_gidx[i]],
                                           isisolated[i] * (theta - va[i]) +
                                               ispvpq[i] * Vm * Vm * gl[i]);

        RAJA::atomicAdd<exago_raja_atomic>(&ge_dev[b_gidx[i] + 1],
                                           isisolated[i] * (Vm - vm[i]) -
                                               ispvpq[i] * Vm * Vm * bl[i]);

        if (include_powerimbalance_variables) {
          double Pimb = x_dev[b_xidxpimb[i]] - x_dev[b_xidxpimb[i] + 1];
          double Qimb = x_dev[b_xidxpimb[i] + 2] - x_dev[b_xidxpimb[i] + 3];
          RAJA::atomicAdd<exago_raja_atomic>(&ge_dev[b_gidx[i]], Pimb);
          RAJA::atomicAdd<exago_raja_atomic>(&ge_dev[b_gidx[i] + 1], Qimb);
        }
      });
  flps += busparams->nbus * 14;

  /* Line contributions */
  double *Gff = lineparams->Gff_dev_;
  double *Gtt = lineparams->Gtt_dev_;
  double *Gft = lineparams->Gft_dev_;
  double *Gtf = lineparams->Gtf_dev_;

  double *Bff = lineparams->Bff_dev_;
  double *Btt = lineparams->Btt_dev_;
  double *Bft = lineparams->Bft_dev_;
  double *Btf = lineparams->Btf_dev_;

  int *xidxf = lineparams->xidxf_dev_;
  int *xidxt = lineparams->xidxt_dev_;
  int *geqidxf = lineparams->geqidxf_dev_;
  int *geqidxt = lineparams->geqidxt_dev_;

  RAJA::forall<exago_raja_exec>(
      RAJA::RangeSegment(0, lineparams->nlineON),
      RAJA_LAMBDA(RAJA::Index_type i) {
        double Pf, Qf, Pt, Qt;
        double thetaf = x_dev[xidxf[i]], Vmf = x_dev[xidxf[i] + 1];
        double thetat = x_dev[xidxt[i]], Vmt = x_dev[xidxt[i] + 1];
        double thetaft = thetaf - thetat;
        double thetatf = thetat - thetaf;

        Pf = Gff[i] * Vmf * Vmf +
             Vmf * Vmt * (Gft[i] * cos(thetaft) + Bft[i] * sin(thetaft));
        Qf = -Bff[i] * Vmf * Vmf +
             Vmf * Vmt * (-Bft[i] * cos(thetaft) + Gft[i] * sin(thetaft));
        Pt = Gtt[i] * Vmt * Vmt +
             Vmt * Vmf * (Gtf[i] * cos(thetatf) + Btf[i] * sin(thetatf));
        Qt = -Btt[i] * Vmt * Vmt +
             Vmt * Vmf * (-Btf[i] * cos(thetatf) + Gtf[i] * sin(thetatf));

        RAJA::atomicAdd<exago_raja_atomic>(&ge_dev[geqidxf[i]], Pf);
        RAJA::atomicAdd<exago_raja_atomic>(&ge_dev[geqidxf[i] + 1], Qf);
        RAJA::atomicAdd<exago_raja_atomic>(&ge_dev[geqidxt[i]], Pt);
        RAJA::atomicAdd<exago_raja_atomic>(&ge_dev[geqidxt[i] + 1], Qt);
      });
  flps += lineparams->nlineON *
          (46 + (4 * EXAGO_FLOPS_COSOP) + (4 * EXAGO_FLOPS_SINOP));
  //  PetscPrintf(MPI_COMM_SELF,"Exit Equality Constraints\n");

  ierr = PetscLogFlops(flps);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/** INEQUALITY CONSTRAINTS **/
PetscErrorCode OPFLOWComputeInequalityConstraintsArray_PBPOLRAJAHIOP(
    OPFLOW opflow, const double *x_dev, double *gi_dev) {
  PbpolModelRajaHiop *pbpolrajahiop =
      reinterpret_cast<PbpolModelRajaHiop *>(opflow->model);
  LINEParamsRajaHiop *lineparams = &pbpolrajahiop->lineparams;
  PetscInt flps = 0;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  //  PetscPrintf(MPI_COMM_SELF,"Entered Inequality Constraints\n");

  // Zero out array
  auto &resmgr = umpire::ResourceManager::getInstance();
  resmgr.memset(gi_dev, 0);

  /* Line contributions */
  double *Gff = lineparams->Gff_dev_;
  double *Gtt = lineparams->Gtt_dev_;
  double *Gft = lineparams->Gft_dev_;
  double *Gtf = lineparams->Gtf_dev_;

  double *Bff = lineparams->Bff_dev_;
  double *Btt = lineparams->Btt_dev_;
  double *Bft = lineparams->Bft_dev_;
  double *Btf = lineparams->Btf_dev_;

  int *linelimidx = lineparams->linelimidx_dev_;
  int *xidxf = lineparams->xidxf_dev_;
  int *xidxt = lineparams->xidxt_dev_;
  int *gineqidx = lineparams->gineqidx_dev_;
  if (lineparams->nlinelim) {
    RAJA::forall<exago_raja_exec>(
        RAJA::RangeSegment(0, lineparams->nlinelim),
        RAJA_LAMBDA(RAJA::Index_type i) {
          int j = linelimidx[i];
          double Pf, Qf, Pt, Qt, Sf2, St2;
          double thetaf = x_dev[xidxf[j]], Vmf = x_dev[xidxf[j] + 1];
          double thetat = x_dev[xidxt[j]], Vmt = x_dev[xidxt[j] + 1];
          double thetaft = thetaf - thetat;
          double thetatf = thetat - thetaf;

          Pf = Gff[j] * Vmf * Vmf +
               Vmf * Vmt * (Gft[j] * cos(thetaft) + Bft[j] * sin(thetaft));
          Qf = -Bff[j] * Vmf * Vmf +
               Vmf * Vmt * (-Bft[j] * cos(thetaft) + Gft[j] * sin(thetaft));
          Pt = Gtt[j] * Vmt * Vmt +
               Vmt * Vmf * (Gtf[j] * cos(thetatf) + Btf[j] * sin(thetatf));
          Qt = -Btt[j] * Vmt * Vmt +
               Vmt * Vmf * (-Btf[j] * cos(thetatf) + Gtf[j] * sin(thetatf));

          Sf2 = Pf * Pf + Qf * Qf;
          St2 = Pt * Pt + Qt * Qt;

          gi_dev[gineqidx[i]] = Sf2;
          gi_dev[gineqidx[i] + 1] = St2;
        });
    flps += lineparams->nlinelim *
            (72 + (4 * EXAGO_FLOPS_COSOP) + (4 * EXAGO_FLOPS_SINOP));
  }
  //  PetscPrintf(MPI_COMM_SELF,"Exit Inequality Constraints\n");
  ierr = PetscLogFlops(flps);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/** OBJECTIVE FUNCTION **/
// Note: This kernel (and all the kernels for this model assume that the data
// has been already allocated on the device. x_dev is pointer to array on the
// GPU
PetscErrorCode OPFLOWComputeObjectiveArray_PBPOLRAJAHIOP(OPFLOW opflow,
                                                         const double *x_dev,
                                                         double *obj) {
  PbpolModelRajaHiop *pbpolrajahiop =
      reinterpret_cast<PbpolModelRajaHiop *>(opflow->model);
  GENParamsRajaHiop *genparams = &pbpolrajahiop->genparams;
  LOADParamsRajaHiop *loadparams = &pbpolrajahiop->loadparams;
  BUSParamsRajaHiop *busparams = &pbpolrajahiop->busparams;
  PetscErrorCode ierr;
  PS ps = opflow->ps;
  int isobj_gencost = opflow->obj_gencost;
  double MVAbase = ps->MVAbase;
  double weight = opflow->weight;

  PetscFunctionBegin;
  //  PetscPrintf(MPI_COMM_SELF,"Entered objective function\n");

  // You need local copies of device pointers so that lambda can capture them
  double *cost_alpha = genparams->cost_alpha_dev_;
  double *cost_beta = genparams->cost_beta_dev_;
  double *cost_gamma = genparams->cost_gamma_dev_;
  int *xidx = genparams->xidx_dev_;
  int *l_xidx = loadparams->xidx_dev_;
  int *b_xidxpimb = busparams->xidxpimb_dev_;

  /* Generator objective function contributions */
  // Set up reduce sum object
  RAJA::ReduceSum<exago_raja_reduce, double> obj_val_sum(0.0);

  if (opflow->objectivetype == MIN_GEN_COST) {
    // Compute reduction on CUDA device
    RAJA::forall<exago_raja_exec>(
        RAJA::RangeSegment(0, genparams->ngenON),
        RAJA_LAMBDA(RAJA::Index_type i) {
          double Pg = x_dev[xidx[i]] * MVAbase;
          obj_val_sum +=
              weight * isobj_gencost *
              (cost_alpha[i] * Pg * Pg + cost_beta[i] * Pg + cost_gamma[i]);
        });
  }

  if (opflow->include_loadloss_variables) {
    double *loadloss_penalty_dev_ = loadparams->loadloss_penalty_dev_;
    RAJA::forall<exago_raja_exec>(
        RAJA::RangeSegment(0, loadparams->nload),
        RAJA_LAMBDA(RAJA::Index_type i) {
          double Pdloss = x_dev[l_xidx[i]];
          double Qdloss = x_dev[l_xidx[i] + 1];
          obj_val_sum +=
              weight * loadloss_penalty_dev_[i] * MVAbase * (Pdloss + Qdloss);
        });
  }

  /* Powerimbalance contributions */
  if (opflow->include_powerimbalance_variables) {
    double *powerimbalance_penalty_dev_ =
        busparams->powerimbalance_penalty_dev_;
    RAJA::forall<exago_raja_exec>(
        RAJA::RangeSegment(0, busparams->nbus),
        RAJA_LAMBDA(RAJA::Index_type i) {
          double Pimbplus = x_dev[b_xidxpimb[i]];
          double Pimbminus = x_dev[b_xidxpimb[i] + 1];
          double Qimbplus = x_dev[b_xidxpimb[i] + 2];
          double Qimbminus = x_dev[b_xidxpimb[i] + 3];
          obj_val_sum += weight * powerimbalance_penalty_dev_[i] * MVAbase *
                         (Pimbplus + Pimbminus + Qimbplus + Qimbminus);
        });
  }

  *obj = static_cast<double>(obj_val_sum.get());
  ierr = PetscLogFlops(genparams->ngenON * 8.0);
  CHKERRQ(ierr);

  //  PetscPrintf(MPI_COMM_SELF,"Exit objective function\n");
  PetscFunctionReturn(0);
}

/** GRADIENT **/
PetscErrorCode OPFLOWComputeGradientArray_PBPOLRAJAHIOP(OPFLOW opflow,
                                                        const double *x_dev,
                                                        double *grad_dev) {
  PbpolModelRajaHiop *pbpolrajahiop =
      reinterpret_cast<PbpolModelRajaHiop *>(opflow->model);
  GENParamsRajaHiop *genparams = &pbpolrajahiop->genparams;
  LOADParamsRajaHiop *loadparams = &pbpolrajahiop->loadparams;
  BUSParamsRajaHiop *busparams = &pbpolrajahiop->busparams;
  PS ps = opflow->ps;
  int isobj_gencost = opflow->obj_gencost;
  double MVAbase = ps->MVAbase;
  double weight = opflow->weight;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  //  PetscPrintf(MPI_COMM_SELF,"Entered gradient function\n");
  // Zero out array
  auto &resmgr = umpire::ResourceManager::getInstance();
  resmgr.memset(grad_dev, 0);

  /* Generator gradient contributions */
  double *cost_alpha = genparams->cost_alpha_dev_;
  double *cost_beta = genparams->cost_beta_dev_;
  int *xidx = genparams->xidx_dev_;
  int *l_xidx = loadparams->xidx_dev_;
  int *b_xidxpimb = busparams->xidxpimb_dev_;

  if (opflow->objectivetype == MIN_GEN_COST) {
    RAJA::forall<exago_raja_exec>(
        RAJA::RangeSegment(0, genparams->ngenON),
        RAJA_LAMBDA(RAJA::Index_type i) {
          double Pg = x_dev[xidx[i]] * MVAbase;
          grad_dev[xidx[i]] = weight * isobj_gencost * MVAbase *
                              (2.0 * cost_alpha[i] * Pg + cost_beta[i]);
        });
  }

  if (opflow->include_loadloss_variables) {
    double *loadloss_penalty_dev_ = loadparams->loadloss_penalty_dev_;
    RAJA::forall<exago_raja_exec>(
        RAJA::RangeSegment(0, loadparams->nload),
        RAJA_LAMBDA(RAJA::Index_type i) {
          double Pdloss = x_dev[l_xidx[i]];
          double Qdloss = x_dev[l_xidx[i] + 1];
          grad_dev[l_xidx[i]] = weight * loadloss_penalty_dev_[i] * MVAbase;
          grad_dev[l_xidx[i] + 1] = weight * loadloss_penalty_dev_[i] * MVAbase;
        });
  }

  if (opflow->include_powerimbalance_variables) {
    double *powerimbalance_penalty_dev_ =
        busparams->powerimbalance_penalty_dev_;
    RAJA::forall<exago_raja_exec>(
        RAJA::RangeSegment(0, busparams->nbus),
        RAJA_LAMBDA(RAJA::Index_type i) {
          grad_dev[b_xidxpimb[i]] =
              weight * powerimbalance_penalty_dev_[i] * MVAbase;
          grad_dev[b_xidxpimb[i] + 1] =
              weight * powerimbalance_penalty_dev_[i] * MVAbase;
          grad_dev[b_xidxpimb[i] + 2] =
              weight * powerimbalance_penalty_dev_[i] * MVAbase;
          grad_dev[b_xidxpimb[i] + 3] =
              weight * powerimbalance_penalty_dev_[i] * MVAbase;
        });
  }

  ierr = PetscLogFlops(genparams->ngenON * 6.0);
  CHKERRQ(ierr);
  //  PetscPrintf(MPI_COMM_SELF,"Exit gradient function\n");

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSetVariableBoundsArray_PBPOLRAJAHIOP(OPFLOW opflow,
                                                          double *xl_dev,
                                                          double *xu_dev) {
  PetscErrorCode ierr;
  double *xl, *xu;
  auto &resmgr = umpire::ResourceManager::getInstance();

  PetscFunctionBegin;

  // The bounds have been already computed on the host, get their pointers
  ierr = VecGetArray(opflow->Xl, &xl);
  CHKERRQ(ierr);
  ierr = VecGetArray(opflow->Xu, &xu);
  CHKERRQ(ierr);

  // Copy host to device
  umpire::Allocator h_allocator_ = resmgr.getAllocator("HOST");
  registerWith(xl, opflow->nx, resmgr, h_allocator_);
  registerWith(xu, opflow->nx, resmgr, h_allocator_);
  resmgr.copy(xl_dev, xl);
  resmgr.copy(xu_dev, xu);

  ierr = VecRestoreArray(opflow->Xl, &xl);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(opflow->Xu, &xu);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

// Note: This kernel (and all the kernels for this model assume that the data
// has been already allocated on the device. xl_dev and xu_dev are pointers to
// arrays on the GPU
PetscErrorCode OPFLOWSetVariableBoundsArray_PBPOLRAJAHIOP_old(OPFLOW opflow,
                                                              double *xl_dev,
                                                              double *xu_dev) {
  PbpolModelRajaHiop *pbpolrajahiop =
      reinterpret_cast<PbpolModelRajaHiop *>(opflow->model);
  BUSParamsRajaHiop *busparams = &pbpolrajahiop->busparams;
  GENParamsRajaHiop *genparams = &pbpolrajahiop->genparams;
  LOADParamsRajaHiop *loadparams = &pbpolrajahiop->loadparams;

  PetscFunctionBegin;
  //  PetscPrintf(MPI_COMM_SELF,"Entered variable bounds\n");
  int *xidx = busparams->xidx_dev_;
  int *ispvpq = busparams->ispvpq_dev_;
  int *isref = busparams->isref_dev_;
  int *isisolated = busparams->isisolated_dev_;
  double *va = busparams->va_dev_;
  double *vm = busparams->vm_dev_;
  double *vmin = busparams->vmin_dev_;
  double *vmax = busparams->vmax_dev_;
  int *b_xidxpimb = busparams->xidxpimb_dev_;
  int include_powerimbalance_variables =
      opflow->include_powerimbalance_variables;

  /* Bounds for bus voltages */
  RAJA::forall<exago_raja_exec>(
      RAJA::RangeSegment(0, busparams->nbus) /* index set here */,
      RAJA_LAMBDA(RAJA::Index_type i) {
        xl_dev[xidx[i]] = ispvpq[i] * PETSC_NINFINITY + isisolated[i] * va[i] +
                          isref[i] * va[i];
        xu_dev[xidx[i]] = ispvpq[i] * PETSC_INFINITY + isisolated[i] * va[i] +
                          isref[i] * va[i];

        xl_dev[xidx[i] + 1] =
            isref[i] * vmin[i] + ispvpq[i] * vmin[i] + isisolated[i] * vm[i];
        xu_dev[xidx[i] + 1] =
            isref[i] * vmax[i] + ispvpq[i] * vmax[i] + isisolated[i] * vm[i];
        /* Bounds for Power Imbalance Variables (second bus variables) */
        if (include_powerimbalance_variables) {
          xl_dev[b_xidxpimb[i]] = xl_dev[b_xidxpimb[i] + 1] =
              xl_dev[b_xidxpimb[i] + 2] = xl_dev[b_xidxpimb[i] + 3] = 0.0;
          xu_dev[b_xidxpimb[i]] = xu_dev[b_xidxpimb[i] + 1] =
              xu_dev[b_xidxpimb[i] + 2] = xu_dev[b_xidxpimb[i] + 3] =
                  PETSC_INFINITY;
        }
      });

  int *idx = genparams->xidx_dev_;
  double *pb = genparams->pb_dev_;
  double *pt = genparams->pt_dev_;
  double *qb = genparams->qb_dev_;
  double *qt = genparams->qt_dev_;

  /* Generator lower and upper bounds on variables */
  RAJA::forall<exago_raja_exec>(
      RAJA::RangeSegment(0, genparams->ngenON),
      RAJA_LAMBDA(RAJA::Index_type i) {
        xl_dev[idx[i]] = pb[i];
        xu_dev[idx[i]] = pt[i];
        xl_dev[idx[i] + 1] = qb[i];
        xu_dev[idx[i] + 1] = qt[i];
      });

  if (opflow->has_gensetpoint) {
    /* Bounds on power deviation and set-point */
    RAJA::forall<exago_raja_exec>(
        RAJA::RangeSegment(0, genparams->ngenON),
        RAJA_LAMBDA(RAJA::Index_type i) {
          xl_dev[idx[i] + 2] = pb[i] - pt[i];
          xu_dev[idx[i] + 2] = pt[i] - pb[i];
          xl_dev[idx[i] + 3] = pb[i];
          xu_dev[idx[i] + 3] = pt[i];
        });
  }

  /* Load loss lower and upper bounds */
  if (opflow->include_loadloss_variables) {
    int *l_xidx = loadparams->xidx_dev_;
    double *pl = loadparams->pl_dev_;
    double *ql = loadparams->ql_dev_;

    RAJA::forall<exago_raja_exec>(
        RAJA::RangeSegment(0, loadparams->nload),
        RAJA_LAMBDA(RAJA::Index_type i) {
          xl_dev[l_xidx[i]] = 0;
          xu_dev[l_xidx[i]] = 0;
          xl_dev[l_xidx[i] + 1] = 0;
          xu_dev[l_xidx[i] + 1] = 0;
          RAJA::atomicMin<exago_raja_atomic>(&xl_dev[l_xidx[i]], pl[i]);
          RAJA::atomicMax<exago_raja_atomic>(&xu_dev[l_xidx[i]], pl[i]);
          RAJA::atomicMin<exago_raja_atomic>(&xl_dev[l_xidx[i] + 1], ql[i]);
          RAJA::atomicMax<exago_raja_atomic>(&xu_dev[l_xidx[i] + 1], ql[i]);
        });
  }

  //  PetscPrintf(MPI_COMM_SELF,"Exit variable bounds\n");

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeSparseInequalityConstraintJacobian_PBPOLRAJAHIOP(
    OPFLOW opflow, const double *x_dev, int *iJacS_dev, int *jJacS_dev,
    double *MJacS_dev) {
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeSparseEqualityConstraintJacobian_PBPOLRAJAHIOP(
    OPFLOW opflow, const double *x_dev, int *iJacS_dev, int *jJacS_dev,
    double *MJacS_dev) {
  PbpolModelRajaHiop *pbpolrajahiop =
      reinterpret_cast<PbpolModelRajaHiop *>(opflow->model);
  GENParamsRajaHiop *genparams = &pbpolrajahiop->genparams;
  LOADParamsRajaHiop *loadparams = &pbpolrajahiop->loadparams;
  BUSParamsRajaHiop *busparams = &pbpolrajahiop->busparams;

  //  PetscPrintf(MPI_COMM_SELF,"Entered sparse jacobian\n");
  if (iJacS_dev != NULL && jJacS_dev != NULL) {

    /* Bus power imbalance contribution */
    int *b_xidxpimb = busparams->xidxpimb_dev_;
    int *b_gidx = busparams->gidx_dev_;
    int *b_jacsp_idx = busparams->jacsp_idx_dev_;
    int *b_jacsq_idx = busparams->jacsq_idx_dev_;
    if (opflow->include_powerimbalance_variables) {
      RAJA::forall<exago_raja_exec>(
          RAJA::RangeSegment(0, busparams->nbus),
          RAJA_LAMBDA(RAJA::Index_type i) {
            iJacS_dev[b_jacsp_idx[i]] = b_gidx[i];
            jJacS_dev[b_jacsp_idx[i]] = b_xidxpimb[i];
            iJacS_dev[b_jacsp_idx[i] + 1] = b_gidx[i];
            jJacS_dev[b_jacsp_idx[i] + 1] = b_xidxpimb[i] + 1;

            iJacS_dev[b_jacsq_idx[i]] = b_gidx[i] + 1;
            jJacS_dev[b_jacsq_idx[i]] = b_xidxpimb[i] + 2;
            iJacS_dev[b_jacsq_idx[i] + 1] = b_gidx[i] + 1;
            jJacS_dev[b_jacsq_idx[i] + 1] = b_xidxpimb[i] + 3;
          });
    }

    /* Generator contributions for row, col entries */
    int *g_gidxbus = genparams->gidxbus_dev_;
    int *g_xidx = genparams->xidx_dev_;
    int *eqjacspbus_idx = genparams->eqjacspbus_idx_dev_;
    int *eqjacsqbus_idx = genparams->eqjacsqbus_idx_dev_;
    RAJA::forall<exago_raja_exec>(
        RAJA::RangeSegment(0, genparams->ngenON),
        RAJA_LAMBDA(RAJA::Index_type i) {
          iJacS_dev[eqjacspbus_idx[i]] = g_gidxbus[i];
          jJacS_dev[eqjacspbus_idx[i]] = g_xidx[i];

          iJacS_dev[eqjacsqbus_idx[i]] = g_gidxbus[i] + 1;
          jJacS_dev[eqjacsqbus_idx[i]] = g_xidx[i] + 1;
        });

    if (opflow->has_gensetpoint) {
      int *eqjacspgen_idx = genparams->eqjacspgen_idx_dev_;
      int *g_geqidxgen = genparams->geqidxgen_dev_;
      int *g_xidx = genparams->xidx_dev_;

      RAJA::forall<exago_raja_exec>(
          RAJA::RangeSegment(0, genparams->ngenON),
          RAJA_LAMBDA(RAJA::Index_type i) {
            iJacS_dev[eqjacspgen_idx[i]] = g_geqidxgen[i];
            jJacS_dev[eqjacspgen_idx[i]] = g_xidx[i];

            iJacS_dev[eqjacspgen_idx[i] + 1] = g_geqidxgen[i];
            jJacS_dev[eqjacspgen_idx[i] + 1] = g_xidx[i] + 2;

            iJacS_dev[eqjacspgen_idx[i] + 2] = g_geqidxgen[i];
            jJacS_dev[eqjacspgen_idx[i] + 2] = g_xidx[i] + 3;

            iJacS_dev[eqjacspgen_idx[i] + 3] = g_geqidxgen[i] + 1;
            jJacS_dev[eqjacspgen_idx[i] + 3] = g_xidx[i] + 3;
          });
    }

    /* Loadloss contributions */
    if (opflow->include_loadloss_variables) {
      int *l_gidx = loadparams->gidx_dev_;
      int *l_xidx = loadparams->xidx_dev_;
      int *l_jacsp_idx = loadparams->jacsp_idx_dev_;

      RAJA::forall<exago_raja_exec>(
          RAJA::RangeSegment(0, loadparams->nload),
          RAJA_LAMBDA(RAJA::Index_type i) {
            iJacS_dev[l_jacsp_idx[i]] = l_gidx[i];
            jJacS_dev[l_jacsp_idx[i]] = l_xidx[i];
            iJacS_dev[l_jacsp_idx[i] + 1] = l_gidx[i] + 1;
            jJacS_dev[l_jacsp_idx[i] + 1] = l_xidx[i] + 1;
          });
    }
  }

  if (MJacS_dev != NULL) {
    /* Bus Contribution - Power imbalance */
    if (opflow->include_powerimbalance_variables) {
      int *b_jacsp_idx = busparams->jacsp_idx_dev_;
      int *b_jacsq_idx = busparams->jacsq_idx_dev_;
      RAJA::forall<exago_raja_exec>(
          RAJA::RangeSegment(0, busparams->nbus),
          RAJA_LAMBDA(RAJA::Index_type i) {
            MJacS_dev[b_jacsp_idx[i]] = 1.0;
            MJacS_dev[b_jacsp_idx[i] + 1] = -1.0;
            MJacS_dev[b_jacsq_idx[i]] = 1.0;
            MJacS_dev[b_jacsq_idx[i] + 1] = -1.0;
          });
    }

    /* Generator contributions */
    int *eqjacspbus_idx = genparams->eqjacspbus_idx_dev_;
    int *eqjacsqbus_idx = genparams->eqjacsqbus_idx_dev_;
    RAJA::forall<exago_raja_exec>(
        RAJA::RangeSegment(0, genparams->ngenON),
        RAJA_LAMBDA(RAJA::Index_type i) {
          MJacS_dev[eqjacspbus_idx[i]] = -1.0;
          MJacS_dev[eqjacsqbus_idx[i]] = -1.0;
        });

    if (opflow->has_gensetpoint) {
      int *eqjacspgen_idx = genparams->eqjacspgen_idx_dev_;

      RAJA::forall<exago_raja_exec>(
          RAJA::RangeSegment(0, genparams->ngenON),
          RAJA_LAMBDA(RAJA::Index_type i) {
            MJacS_dev[eqjacspgen_idx[i]] = -1.0;
            MJacS_dev[eqjacspgen_idx[i] + 1] = 1.0;
            MJacS_dev[eqjacspgen_idx[i] + 2] = 1.0;
            MJacS_dev[eqjacspgen_idx[i] + 3] = 1.0;
          });
    }

    /* Loadloss contributions - 2 contributions expected */
    if (opflow->include_loadloss_variables) {
      int *l_jacsp_idx = loadparams->jacsp_idx_dev_;
      RAJA::forall<exago_raja_exec>(
          RAJA::RangeSegment(0, loadparams->nload),
          RAJA_LAMBDA(RAJA::Index_type i) {
            MJacS_dev[l_jacsp_idx[i]] = -1;
            MJacS_dev[l_jacsp_idx[i] + 1] = -1;
          });
    }
  }
  //  PetscPrintf(MPI_COMM_SELF,"Exit sparse jacobian\n");
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeSparseHessian_PBPOLRAJAHIOP(
    OPFLOW opflow, const double *x_dev, const double *lambda_dev, int *iHSS_dev,
    int *jHSS_dev, double *MHSS_dev) {
  PbpolModelRajaHiop *pbpolrajahiop =
      reinterpret_cast<PbpolModelRajaHiop *>(opflow->model);
  PetscErrorCode ierr;
  GENParamsRajaHiop *genparams = &pbpolrajahiop->genparams;
  LOADParamsRajaHiop *loadparams = &pbpolrajahiop->loadparams;
  PS ps = opflow->ps;
  double obj_factor = opflow->obj_factor;
  int isobj_gencost = opflow->obj_gencost;
  double MVAbase = ps->MVAbase;
  double weight = opflow->weight;
  PetscInt flps = 0;
  //  PetscPrintf(MPI_COMM_SELF,"Entered sparse Hessian\n");

  if (iHSS_dev != NULL && jHSS_dev != NULL) {

    /* Generator contributions for row,col numbers */
    int *g_xidx = genparams->xidx_dev_;
    int *hesssp_idx = genparams->hesssp_idx_dev_;
    RAJA::forall<exago_raja_exec>(
        RAJA::RangeSegment(0, genparams->ngenON),
        RAJA_LAMBDA(RAJA::Index_type i) {
          iHSS_dev[hesssp_idx[i]] = g_xidx[i];
          jHSS_dev[hesssp_idx[i]] = g_xidx[i];
        });

    /* Loadloss contributions - two contributions*/
    if (opflow->include_loadloss_variables) {
      int *l_xidx = loadparams->xidx_dev_;
      int *l_hesssp_idx = loadparams->hesssp_idx_dev_;
      RAJA::forall<exago_raja_exec>(
          RAJA::RangeSegment(0, loadparams->nload),
          RAJA_LAMBDA(RAJA::Index_type i) {
            iHSS_dev[l_hesssp_idx[i]] = l_xidx[i];
            jHSS_dev[l_hesssp_idx[i]] = l_xidx[i];
            iHSS_dev[l_hesssp_idx[i] + 1] = l_xidx[i] + 1;
            jHSS_dev[l_hesssp_idx[i] + 1] = l_xidx[i] + 1;
          });
    }
  }

  if (MHSS_dev != NULL) {

    /* Generator contributions */
    if (opflow->objectivetype == MIN_GEN_COST) {
      int *hesssp_idx = genparams->hesssp_idx_dev_;
      double *cost_alpha = genparams->cost_alpha_dev_;

      RAJA::forall<exago_raja_exec>(
          RAJA::RangeSegment(0, genparams->ngenON),
          RAJA_LAMBDA(RAJA::Index_type i) {
            MHSS_dev[hesssp_idx[i]] = weight * isobj_gencost * obj_factor *
                                      2.0 * cost_alpha[i] * MVAbase * MVAbase;
          });
      flps += 6 * genparams->ngenON;
    } else if (opflow->objectivetype == NO_OBJ) {
      int *hesssp_idx = genparams->hesssp_idx_dev_;
      RAJA::forall<exago_raja_exec>(
          RAJA::RangeSegment(0, genparams->ngenON),
          RAJA_LAMBDA(RAJA::Index_type i) { MHSS_dev[hesssp_idx[i]] = 0.0; });
    }

    /* Loadloss contributions - 2 contributions expected */
    if (opflow->include_loadloss_variables) {
      int *l_hesssp_idx = loadparams->hesssp_idx_dev_;
      double *loadloss_penalty_dev_ = loadparams->loadloss_penalty_dev_;
      RAJA::forall<exago_raja_exec>(
          RAJA::RangeSegment(0, loadparams->nload),
          RAJA_LAMBDA(RAJA::Index_type i) {
            MHSS_dev[l_hesssp_idx[i]] = 0.0;
            MHSS_dev[l_hesssp_idx[i] + 1] = 0.0;
          });
    }
  }
  //  PetscPrintf(MPI_COMM_SELF,"Exit sparse hessian\n");

  ierr = PetscLogFlops(flps);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/**
 * @param[inout] JacD_dev Jacobian matrix with size * 2*busparams->nbus by
 * nxdense
 */
PetscErrorCode OPFLOWComputeDenseEqualityConstraintJacobian_PBPOLRAJAHIOP(
    OPFLOW opflow, const double *x_dev, double *JacD_dev) {
  PbpolModelRajaHiop *pbpolrajahiop =
      reinterpret_cast<PbpolModelRajaHiop *>(opflow->model);
  BUSParamsRajaHiop *busparams = &pbpolrajahiop->busparams;
  LINEParamsRajaHiop *lineparams = &pbpolrajahiop->lineparams;
  // LOADParamsRajaHiop     *loadparams=&pbpolrajahiop->loadparams;
  int nxsparse = opflow->nxsparse;
  const int nxdense = opflow->nxdense;
  double flps = 0.0;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  //  PetscPrintf(MPI_COMM_SELF,"Entered equality constrained dense
  //  jacobian\n");

  if (JacD_dev == NULL)
    PetscFunctionReturn(0);

  /** Although we will perform atomic operations on this view, we do not want
   * to use `RAJA::make_atomic_view` because explicit `RAJA::atomicAdd`s make
   * the atomic requirement clearer in the code */
  RAJA::View<double, RAJA::Layout<2>> JacD_view(JacD_dev, opflow->nconeq,
                                                nxdense);

  /* Zero out JacD */
  auto &resmgr = umpire::ResourceManager::getInstance();
  umpire::Allocator alloc;
#ifdef EXAGO_ENABLE_GPU
  alloc = resmgr.getAllocator("DEVICE");
#else
  alloc = resmgr.getAllocator("HOST");
#endif
  registerWith(JacD_dev, opflow->nconeq * nxdense, resmgr, alloc);
  resmgr.memset(JacD_dev, 0);

  /* Jacobian from bus contributions */
  int *isisolated = busparams->isisolated_dev_;
  int *ispvpq = busparams->ispvpq_dev_;
  double *gl = busparams->gl_dev_;
  double *bl = busparams->bl_dev_;
  int *b_xidx = busparams->xidx_dev_;
  int *b_gidx = busparams->gidx_dev_;
  RAJA::forall<exago_raja_exec>(
      RAJA::RangeSegment(0, busparams->nbus), RAJA_LAMBDA(RAJA::Index_type i) {
        double Vm = x_dev[b_xidx[i] + 1];
        int row[2], col[4];
        double val[8];

        row[0] = b_gidx[i];
        row[1] = b_gidx[i] + 1;

        col[0] = b_xidx[i] - nxsparse;
        col[1] = b_xidx[i] + 1 - nxsparse;

        val[0] = isisolated[i] * 1.0 + ispvpq[i] * 0.0;
        val[1] = isisolated[i] * 0.0 + ispvpq[i] * 2 * Vm * gl[i];
        val[2] = 0.0;
        val[3] = isisolated[i] * 1.0 + ispvpq[i] * -2 * Vm * bl[i];

        RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[0], col[0]), val[0]);
        RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[0], col[1]), val[1]);
        RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[1], col[0]), val[2]);
        RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[1], col[1]), val[3]);
      });
  flps += 14 * busparams->nbus;

  /* Jacobian from line contributions */
  double *Gff = lineparams->Gff_dev_;
  double *Gtt = lineparams->Gtt_dev_;
  double *Gft = lineparams->Gft_dev_;
  double *Gtf = lineparams->Gtf_dev_;

  double *Bff = lineparams->Bff_dev_;
  double *Btt = lineparams->Btt_dev_;
  double *Bft = lineparams->Bft_dev_;
  double *Btf = lineparams->Btf_dev_;

  int *xidxf = lineparams->xidxf_dev_;
  int *xidxt = lineparams->xidxt_dev_;
  int *geqidxf = lineparams->geqidxf_dev_;
  int *geqidxt = lineparams->geqidxt_dev_;

  RAJA::forall<exago_raja_exec>(
      RAJA::RangeSegment(0, lineparams->nlineON),
      RAJA_LAMBDA(RAJA::Index_type i) {
        int row[2], col[4];
        double val[8];
        double thetaf = x_dev[xidxf[i]], Vmf = x_dev[xidxf[i] + 1];
        double thetat = x_dev[xidxt[i]], Vmt = x_dev[xidxt[i] + 1];
        double thetaft = thetaf - thetat;
        double thetatf = thetat - thetaf;

        row[0] = geqidxf[i];
        row[1] = geqidxf[i] + 1;

        col[0] = xidxf[i] - nxsparse;
        col[1] = xidxf[i] + 1 - nxsparse;
        col[2] = xidxt[i] - nxsparse;
        col[3] = xidxt[i] + 1 - nxsparse;

        /* dPf_dthetaf */
        val[0] = Vmf * Vmt * (-Gft[i] * sin(thetaft) + Bft[i] * cos(thetaft));
        /*dPf_dVmf */
        val[1] = 2 * Gff[i] * Vmf +
                 Vmt * (Gft[i] * cos(thetaft) + Bft[i] * sin(thetaft));
        /*dPf_dthetat */
        val[2] = Vmf * Vmt * (Gft[i] * sin(thetaft) - Bft[i] * cos(thetaft));
        /* dPf_dVmt */
        val[3] = Vmf * (Gft[i] * cos(thetaft) + Bft[i] * sin(thetaft));

        RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[0], col[0]), val[0]);
        RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[0], col[1]), val[1]);
        RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[0], col[2]), val[2]);
        RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[0], col[3]), val[3]);

        /* dQf_dthetaf */
        val[4] = Vmf * Vmt * (Bft[i] * sin(thetaft) + Gft[i] * cos(thetaft));
        /* dQf_dVmf */
        val[5] = -2 * Bff[i] * Vmf +
                 Vmt * (-Bft[i] * cos(thetaft) + Gft[i] * sin(thetaft));
        /* dQf_dthetat */
        val[6] = Vmf * Vmt * (-Bft[i] * sin(thetaft) - Gft[i] * cos(thetaft));
        /* dQf_dVmt */
        val[7] = Vmf * (-Bft[i] * cos(thetaft) + Gft[i] * sin(thetaft));

        RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[1], col[0]), val[4]);
        RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[1], col[1]), val[5]);
        RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[1], col[2]), val[6]);
        RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[1], col[3]), val[7]);

        row[0] = geqidxt[i];
        row[1] = geqidxt[i] + 1;

        col[0] = xidxt[i] - nxsparse;
        col[1] = xidxt[i] + 1 - nxsparse;
        col[2] = xidxf[i] - nxsparse;
        col[3] = xidxf[i] + 1 - nxsparse;

        /* dPt_dthetat */
        val[0] = Vmt * Vmf * (-Gtf[i] * sin(thetatf) + Btf[i] * cos(thetatf));
        /* dPt_dVmt */
        val[1] = 2 * Gtt[i] * Vmt +
                 Vmf * (Gtf[i] * cos(thetatf) + Btf[i] * sin(thetatf));
        /* dPt_dthetaf */
        val[2] = Vmt * Vmf * (Gtf[i] * sin(thetatf) - Btf[i] * cos(thetatf));
        /* dPt_dVmf */
        val[3] = Vmt * (Gtf[i] * cos(thetatf) + Btf[i] * sin(thetatf));

        RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[0], col[0]), val[0]);
        RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[0], col[1]), val[1]);
        RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[0], col[2]), val[2]);
        RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[0], col[3]), val[3]);

        /* dQt_dthetat */
        val[4] = Vmt * Vmf * (Btf[i] * sin(thetatf) + Gtf[i] * cos(thetatf));
        /* dQt_dVmt */
        val[5] = -2 * Btt[i] * Vmt +
                 Vmf * (-Btf[i] * cos(thetatf) + Gtf[i] * sin(thetatf));
        /* dQt_dthetaf */
        val[6] = Vmt * Vmf * (-Btf[i] * sin(thetatf) - Gtf[i] * cos(thetatf));
        /* dQt_dVmf */
        val[7] = Vmt * (-Btf[i] * cos(thetatf) + Gtf[i] * sin(thetatf));

        RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[1], col[0]), val[4]);
        RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[1], col[1]), val[5]);
        RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[1], col[2]), val[6]);
        RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[1], col[3]), val[7]);
      });
  flps += (188 + (16 * EXAGO_FLOPS_COSOP) + (16 * EXAGO_FLOPS_SINOP)) *
          lineparams->nlineON;

  ierr = PetscLogFlops(flps);
  CHKERRQ(ierr);

  //  PetscPrintf(MPI_COMM_SELF,"Exit equality dense jacobian\n");

  PetscFunctionReturn(0);
}

/**
 * @param[inout] JacD_dev Jacobian matrix with size * 2*lineparams->nlinelim by
 * nxdense
 */
PetscErrorCode OPFLOWComputeDenseInequalityConstraintJacobian_PBPOLRAJAHIOP(
    OPFLOW opflow, const double *x_dev, double *JacD_dev) {
  PbpolModelRajaHiop *pbpolrajahiop =
      reinterpret_cast<PbpolModelRajaHiop *>(opflow->model);
  LINEParamsRajaHiop *lineparams = &pbpolrajahiop->lineparams;
  PetscErrorCode ierr;
  double flps = 0.0;
  int nxsparse = opflow->nxsparse;
  int nxdense = opflow->nxdense;
  int nx = opflow->nx;

  PetscFunctionBegin;
  //  PetscPrintf(MPI_COMM_SELF,"Enter inequality dense jacobian\n");
  /* Return if there are no inequality constraints */
  if (!lineparams->nlinelim) {
    //    PetscPrintf(MPI_COMM_SELF,"No inequality constraints. Exit inequality
    //    dense jacobian\n");
    PetscFunctionReturn(0);
  }

  if (JacD_dev == NULL)
    PetscFunctionReturn(0);

  /** Although we will perform atomic operations on this view, we do not want
   * to use `RAJA::make_atomic_view` because explicit `RAJA::atomicAdd`s make
   * the atomic requirement clearer in the code */
  RAJA::View<double, RAJA::Layout<2>> JacD_view(
      JacD_dev, 2 * lineparams->nlinelim, nx - nxsparse);

  /* Zero out JacD */
  auto &resmgr = umpire::ResourceManager::getInstance();
  umpire::Allocator alloc;
#ifdef EXAGO_ENABLE_GPU
  alloc = resmgr.getAllocator("DEVICE");
#else
  alloc = resmgr.getAllocator("HOST");
#endif
  registerWith(JacD_dev, 2 * lineparams->nlinelim * nxdense, resmgr, alloc);
  resmgr.memset(JacD_dev, 0);

  /* Line contributions */
  double *Gff_arr = lineparams->Gff_dev_;
  double *Gtt_arr = lineparams->Gtt_dev_;
  double *Gft_arr = lineparams->Gft_dev_;
  double *Gtf_arr = lineparams->Gtf_dev_;

  double *Bff_arr = lineparams->Bff_dev_;
  double *Btt_arr = lineparams->Btt_dev_;
  double *Bft_arr = lineparams->Bft_dev_;
  double *Btf_arr = lineparams->Btf_dev_;

  int *linelimidx = lineparams->linelimidx_dev_;
  int *xidxf = lineparams->xidxf_dev_;
  int *xidxt = lineparams->xidxt_dev_;
  int *gineqidx = lineparams->gineqidx_dev_;
  RAJA::forall<exago_raja_exec>(
      RAJA::RangeSegment(0, lineparams->nlinelim),
      RAJA_LAMBDA(RAJA::Index_type i) {
        int j = linelimidx[i];
        int row[2], col[4];
        double val[4];
        double Pf, Qf, Pt, Qt;
        double thetaf = x_dev[xidxf[j]], Vmf = x_dev[xidxf[j] + 1];
        double thetat = x_dev[xidxt[j]], Vmt = x_dev[xidxt[j] + 1];
        double thetaft = thetaf - thetat;
        double thetatf = thetat - thetaf;
        double dSf2_dPf, dSf2_dQf, dSt2_dPt, dSt2_dQt;
        double dPf_dthetaf, dPf_dVmf, dPf_dthetat, dPf_dVmt;
        double dQf_dthetaf, dQf_dVmf, dQf_dthetat, dQf_dVmt;
        double dPt_dthetaf, dPt_dVmf, dPt_dthetat, dPt_dVmt;
        double dQt_dthetaf, dQt_dVmf, dQt_dthetat, dQt_dVmt;
        double dSf2_dthetaf, dSf2_dVmf, dSf2_dthetat, dSf2_dVmt;
        double dSt2_dthetaf, dSt2_dVmf, dSt2_dthetat, dSt2_dVmt;
        double Gff = Gff_arr[j], Bff = Bff_arr[j];
        double Gft = Gft_arr[j], Bft = Bft_arr[j];
        double Gtf = Gtf_arr[j], Btf = Btf_arr[j];
        double Gtt = Gtt_arr[j], Btt = Btt_arr[j];

        Pf = Gff * Vmf * Vmf +
             Vmf * Vmt * (Gft * cos(thetaft) + Bft * sin(thetaft));
        Qf = -Bff * Vmf * Vmf +
             Vmf * Vmt * (-Bft * cos(thetaft) + Gft * sin(thetaft));
        Pt = Gtt * Vmt * Vmt +
             Vmt * Vmf * (Gtf * cos(thetatf) + Btf * sin(thetatf));
        Qt = -Btt * Vmt * Vmt +
             Vmt * Vmf * (-Btf * cos(thetatf) + Gtf * sin(thetatf));

        dSf2_dPf = 2 * Pf;
        dSf2_dQf = 2 * Qf;
        dSt2_dPt = 2 * Pt;
        dSt2_dQt = 2 * Qt;

        dPf_dthetaf = Vmf * Vmt * (-Gft * sin(thetaft) + Bft * cos(thetaft));
        dPf_dVmf =
            2 * Gff * Vmf + Vmt * (Gft * cos(thetaft) + Bft * sin(thetaft));
        dPf_dthetat = Vmf * Vmt * (Gft * sin(thetaft) - Bft * cos(thetaft));
        dPf_dVmt = Vmf * (Gft * cos(thetaft) + Bft * sin(thetaft));

        dQf_dthetaf = Vmf * Vmt * (Bft * sin(thetaft) + Gft * cos(thetaft));
        dQf_dVmf =
            -2 * Bff * Vmf + Vmt * (-Bft * cos(thetaft) + Gft * sin(thetaft));
        dQf_dthetat = Vmf * Vmt * (-Bft * sin(thetaft) - Gft * cos(thetaft));
        dQf_dVmt = Vmf * (-Bft * cos(thetaft) + Gft * sin(thetaft));

        dPt_dthetat = Vmt * Vmf * (-Gtf * sin(thetatf) + Btf * cos(thetatf));
        dPt_dVmt =
            2 * Gtt * Vmt + Vmf * (Gtf * cos(thetatf) + Btf * sin(thetatf));
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

        row[0] = gineqidx[i];

        col[0] = xidxf[j] - nxsparse;
        col[1] = xidxf[j] + 1 - nxsparse;
        col[2] = xidxt[j] - nxsparse;
        col[3] = xidxt[j] + 1 - nxsparse;

        val[0] = dSf2_dthetaf;
        val[1] = dSf2_dVmf;
        val[2] = dSf2_dthetat;
        val[3] = dSf2_dVmt;

        RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[0], col[0]), val[0]);
        RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[0], col[1]), val[1]);
        RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[0], col[2]), val[2]);
        RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[0], col[3]), val[3]);

        dSt2_dthetaf = dSt2_dPt * dPt_dthetaf + dSt2_dQt * dQt_dthetaf;
        dSt2_dthetat = dSt2_dPt * dPt_dthetat + dSt2_dQt * dQt_dthetat;
        dSt2_dVmf = dSt2_dPt * dPt_dVmf + dSt2_dQt * dQt_dVmf;
        dSt2_dVmt = dSt2_dPt * dPt_dVmt + dSt2_dQt * dQt_dVmt;

        row[0] = gineqidx[i] + 1;

        col[0] = xidxt[j] - nxsparse;
        col[1] = xidxt[j] + 1 - nxsparse;
        col[2] = xidxf[j] - nxsparse;
        col[3] = xidxf[j] + 1 - nxsparse;

        val[0] = dSt2_dthetat;
        val[1] = dSt2_dVmt;
        val[2] = dSt2_dthetaf;
        val[3] = dSt2_dVmf;

        RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[0], col[0]), val[0]);
        RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[0], col[1]), val[1]);
        RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[0], col[2]), val[2]);
        RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[0], col[3]), val[3]);
      });
  //  PetscPrintf(MPI_COMM_SELF,"Exit inequality dense jacobian\n");
  flps += (183 + (20 * EXAGO_FLOPS_COSOP) + (20 * EXAGO_FLOPS_SINOP)) *
          lineparams->nlinelim;
  ierr = PetscLogFlops(flps);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/**
 * @param[inout] HDD_dev Hessian matrix with size nxdense x nxdense
 */
PetscErrorCode OPFLOWComputeDenseEqualityConstraintHessian_PBPOLRAJAHIOP(
    OPFLOW opflow, const double *x_dev, const double *lambda_dev,
    double *HDD_dev) {
  PbpolModelRajaHiop *pbpolrajahiop =
      reinterpret_cast<PbpolModelRajaHiop *>(opflow->model);
  BUSParamsRajaHiop *busparams = &pbpolrajahiop->busparams;
  LINEParamsRajaHiop *lineparams = &pbpolrajahiop->lineparams;
  int nxsparse = opflow->nxsparse;
  const int nxdense = opflow->nxdense;
  PetscErrorCode ierr;
  PetscInt flps = 0;

  PetscFunctionBegin;
  //  PetscPrintf(MPI_COMM_SELF,"Enter equality dense hessian\n");

  /* Hessian from bus contributions */
  int *b_xidx = busparams->xidx_dev_;
  int *b_gidx = busparams->gidx_dev_;
  int *ispvpq = busparams->ispvpq_dev_;
  double *gl = busparams->gl_dev_;
  double *bl = busparams->bl_dev_;

  /** Although we will perform atomic operations on this view, we do not want
   * to use `RAJA::make_atomic_view` because explicit `RAJA::atomicAdd`s make
   * the atomic requirement clearer in the code */
  RAJA::View<double, RAJA::Layout<2>> HDD_view(HDD_dev, nxdense, nxdense);

  RAJA::forall<exago_raja_exec>(
      RAJA::RangeSegment(0, busparams->nbus), RAJA_LAMBDA(RAJA::Index_type i) {
        int row, col;
        double val;
        row = b_xidx[i] + 1 - nxsparse;
        col = row;
        val = ispvpq[i] * (lambda_dev[b_gidx[i]] * 2 * gl[i] +
                           lambda_dev[b_gidx[i] + 1] * (-2 * bl[i]));
        RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row, col), val);
      });
  flps += 7 * busparams->nbus;

  /* Hessian from line contributions */
  double *Gff_arr = lineparams->Gff_dev_;
  double *Gtt_arr = lineparams->Gtt_dev_;
  double *Gft_arr = lineparams->Gft_dev_;
  double *Gtf_arr = lineparams->Gtf_dev_;

  double *Bff_arr = lineparams->Bff_dev_;
  double *Btt_arr = lineparams->Btt_dev_;
  double *Bft_arr = lineparams->Bft_dev_;
  double *Btf_arr = lineparams->Btf_dev_;

  int *xidxf = lineparams->xidxf_dev_;
  int *xidxt = lineparams->xidxt_dev_;
  int *geqidxf = lineparams->geqidxf_dev_;
  int *geqidxt = lineparams->geqidxt_dev_;

  RAJA::forall<exago_raja_exec>(
      RAJA::RangeSegment(0, lineparams->nlineON),
      RAJA_LAMBDA(RAJA::Index_type i) {
        int gloc;
        int row[2], col[4];
        double val[8];
        double Gff, Bff, Gft, Bft, Gtf, Btf, Gtt, Btt;
        Gff = Gff_arr[i];
        Bff = Bff_arr[i];
        Gft = Gft_arr[i];
        Bft = Bft_arr[i];
        Gtf = Gtf_arr[i];
        Btf = Btf_arr[i];
        Gtt = Gtt_arr[i];
        Btt = Btt_arr[i];

        double thetaf = x_dev[xidxf[i]], Vmf = x_dev[xidxf[i] + 1];
        double thetat = x_dev[xidxt[i]], Vmt = x_dev[xidxt[i] + 1];
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

        double dQf_dthetaf_dthetaf, dQf_dthetaf_dVmf, dQf_dthetaf_dthetat,
            dQf_dthetaf_dVmt;
        double dQf_dVmf_dthetaf, dQf_dVmf_dVmf, dQf_dVmf_dthetat, dQf_dVmf_dVmt;
        double dQf_dthetat_dthetaf, dQf_dthetat_dVmf, dQf_dthetat_dthetat,
            dQf_dthetat_dVmt;
        double dQf_dVmt_dthetaf, dQf_dVmt_dVmf, dQf_dVmt_dthetat, dQf_dVmt_dVmt;

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

        row[0] = xidxf[i] - nxsparse;
        row[1] = xidxf[i] + 1 - nxsparse;
        col[0] = xidxf[i] - nxsparse;
        col[1] = xidxf[i] + 1 - nxsparse;
        col[2] = xidxt[i] - nxsparse;
        col[3] = xidxt[i] + 1 - nxsparse;

        gloc = geqidxf[i];

        val[0] = val[1] = val[2] = val[3] = val[4] = val[5] = val[6] = val[7] =
            0.0;

        val[0] = lambda_dev[gloc] * dPf_dthetaf_dthetaf +
                 lambda_dev[gloc + 1] * dQf_dthetaf_dthetaf;
        val[1] = lambda_dev[gloc] * dPf_dthetaf_dVmf +
                 lambda_dev[gloc + 1] * dQf_dthetaf_dVmf;
        val[2] = lambda_dev[gloc] * dPf_dthetaf_dthetat +
                 lambda_dev[gloc + 1] * dQf_dthetaf_dthetat;
        val[3] = lambda_dev[gloc] * dPf_dthetaf_dVmt +
                 lambda_dev[gloc + 1] * dQf_dthetaf_dVmt;

        val[4] = lambda_dev[gloc] * dPf_dVmf_dthetaf +
                 lambda_dev[gloc + 1] * dQf_dVmf_dthetaf;
        val[5] = lambda_dev[gloc] * dPf_dVmf_dVmf +
                 lambda_dev[gloc + 1] * dQf_dVmf_dVmf;
        val[6] = lambda_dev[gloc] * dPf_dVmf_dthetat +
                 lambda_dev[gloc + 1] * dQf_dVmf_dthetat;
        val[7] = lambda_dev[gloc] * dPf_dVmf_dVmt +
                 lambda_dev[gloc + 1] * dQf_dVmf_dVmt;

        RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[0]), val[0]);
        RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[1]), val[1]);
        RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[2]), val[2]);
        RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[3]), val[3]);

        RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[1], col[0]), val[4]);
        RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[1], col[1]), val[5]);
        RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[1], col[2]), val[6]);
        RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[1], col[3]), val[7]);

        row[0] = xidxt[i] - nxsparse;
        row[1] = xidxt[i] + 1 - nxsparse;

        col[0] = xidxf[i] - nxsparse;
        col[1] = xidxf[i] + 1 - nxsparse;
        col[2] = xidxt[i] - nxsparse;
        col[3] = xidxt[i] + 1 - nxsparse;

        val[0] = val[1] = val[2] = val[3] = val[4] = val[5] = val[6] = val[7] =
            0.0;

        val[0] = lambda_dev[gloc] * dPf_dthetat_dthetaf +
                 lambda_dev[gloc + 1] * dQf_dthetat_dthetaf;
        val[1] = lambda_dev[gloc] * dPf_dthetat_dVmf +
                 lambda_dev[gloc + 1] * dQf_dthetat_dVmf;
        val[2] = lambda_dev[gloc] * dPf_dthetat_dthetat +
                 lambda_dev[gloc + 1] * dQf_dthetat_dthetat;
        val[3] = lambda_dev[gloc] * dPf_dthetat_dVmt +
                 lambda_dev[gloc + 1] * dQf_dthetat_dVmt;

        val[4] = lambda_dev[gloc] * dPf_dVmt_dthetaf +
                 lambda_dev[gloc + 1] * dQf_dVmt_dthetaf;
        val[5] = lambda_dev[gloc] * dPf_dVmt_dVmf +
                 lambda_dev[gloc + 1] * dQf_dVmt_dVmf;
        val[6] = lambda_dev[gloc] * dPf_dVmt_dthetat +
                 lambda_dev[gloc + 1] * dQf_dVmt_dthetat;
        val[7] = lambda_dev[gloc] * dPf_dVmt_dVmt +
                 lambda_dev[gloc + 1] * dQf_dVmt_dVmt;

        RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[0]), val[0]);
        RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[1]), val[1]);
        RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[2]), val[2]);
        RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[3]), val[3]);

        RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[1], col[0]), val[4]);
        RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[1], col[1]), val[5]);
        RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[1], col[2]), val[6]);
        RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[1], col[3]), val[7]);

        //    ierr = MatSetValues(H,2,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);

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

        double dQt_dthetaf_dthetaf, dQt_dthetaf_dVmf, dQt_dthetaf_dthetat,
            dQt_dthetaf_dVmt;
        double dQt_dVmf_dthetaf, dQt_dVmf_dVmf, dQt_dVmf_dthetat, dQt_dVmf_dVmt;
        double dQt_dthetat_dthetaf, dQt_dthetat_dVmf, dQt_dthetat_dthetat,
            dQt_dthetat_dVmt;
        double dQt_dVmt_dthetaf, dQt_dVmt_dVmf, dQt_dVmt_dthetat, dQt_dVmt_dVmt;

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

        row[0] = xidxt[i] - nxsparse;
        row[1] = xidxt[i] + 1 - nxsparse;
        col[0] = xidxt[i] - nxsparse;
        col[1] = xidxt[i] + 1 - nxsparse;
        col[2] = xidxf[i] - nxsparse;
        col[3] = xidxf[i] + 1 - nxsparse;

        gloc = geqidxt[i];

        val[0] = val[1] = val[2] = val[3] = val[4] = val[5] = val[6] = val[7] =
            0.0;

        val[0] = lambda_dev[gloc] * dPt_dthetat_dthetat +
                 lambda_dev[gloc + 1] * dQt_dthetat_dthetat;
        val[1] = lambda_dev[gloc] * dPt_dthetat_dVmt +
                 lambda_dev[gloc + 1] * dQt_dthetat_dVmt;
        val[2] = lambda_dev[gloc] * dPt_dthetat_dthetaf +
                 lambda_dev[gloc + 1] * dQt_dthetat_dthetaf;
        val[3] = lambda_dev[gloc] * dPt_dthetat_dVmf +
                 lambda_dev[gloc + 1] * dQt_dthetat_dVmf;

        val[4] = lambda_dev[gloc] * dPt_dVmt_dthetat +
                 lambda_dev[gloc + 1] * dQt_dVmt_dthetat;
        val[5] = lambda_dev[gloc] * dPt_dVmt_dVmt +
                 lambda_dev[gloc + 1] * dQt_dVmt_dVmt;
        val[6] = lambda_dev[gloc] * dPt_dVmt_dthetaf +
                 lambda_dev[gloc + 1] * dQt_dVmt_dthetaf;
        val[7] = lambda_dev[gloc] * dPt_dVmt_dVmf +
                 lambda_dev[gloc + 1] * dQt_dVmt_dVmf;

        RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[0]), val[0]);
        RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[1]), val[1]);
        RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[2]), val[2]);
        RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[3]), val[3]);

        RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[1], col[0]), val[4]);
        RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[1], col[1]), val[5]);
        RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[1], col[2]), val[6]);
        RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[1], col[3]), val[7]);

        row[0] = xidxf[i] - nxsparse;
        row[1] = xidxf[i] + 1 - nxsparse;
        col[0] = xidxt[i] - nxsparse;
        col[1] = xidxt[i] + 1 - nxsparse;
        col[2] = xidxf[i] - nxsparse;
        col[3] = xidxf[i] + 1 - nxsparse;

        val[0] = val[1] = val[2] = val[3] = val[4] = val[5] = val[6] = val[7] =
            0.0;

        val[0] = lambda_dev[gloc] * dPt_dthetaf_dthetat +
                 lambda_dev[gloc + 1] * dQt_dthetaf_dthetat;
        val[1] = lambda_dev[gloc] * dPt_dthetaf_dVmt +
                 lambda_dev[gloc + 1] * dQt_dthetaf_dVmt;
        val[2] = lambda_dev[gloc] * dPt_dthetaf_dthetaf +
                 lambda_dev[gloc + 1] * dQt_dthetaf_dthetaf;
        val[3] = lambda_dev[gloc] * dPt_dthetaf_dVmf +
                 lambda_dev[gloc + 1] * dQt_dthetaf_dVmf;

        RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[0]), val[0]);
        RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[1]), val[1]);
        RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[2]), val[2]);
        RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[3]), val[3]);

        val[4] = lambda_dev[gloc] * dPt_dVmf_dthetat +
                 lambda_dev[gloc + 1] * dQt_dVmf_dthetat;
        val[5] = lambda_dev[gloc] * dPt_dVmf_dVmt +
                 lambda_dev[gloc + 1] * dQt_dVmf_dVmt;
        val[6] = lambda_dev[gloc] * dPt_dVmf_dthetaf +
                 lambda_dev[gloc + 1] * dQt_dVmf_dthetaf;
        val[7] = lambda_dev[gloc] * dPt_dVmf_dVmf +
                 lambda_dev[gloc + 1] * dQt_dVmf_dVmf;

        RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[1], col[0]), val[4]);
        RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[1], col[1]), val[5]);
        RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[1], col[2]), val[6]);
        RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[1], col[3]), val[7]);
      });

  flps += (56 * (EXAGO_FLOPS_SINOP + EXAGO_FLOPS_SINOP) + 462.0) *
          lineparams->nlineON;
  ierr = PetscLogFlops(flps);
  CHKERRQ(ierr);
  //  PetscPrintf(MPI_COMM_SELF,"Exit equality dense hessian\n");

  PetscFunctionReturn(0);
}

/**
 * @param[inout] HDD_dev Hessian matrix with size nxdense x nxdense
 */
PetscErrorCode OPFLOWComputeDenseInequalityConstraintHessian_PBPOLRAJAHIOP(
    OPFLOW opflow, const double *x_dev, const double *lambda_dev,
    double *HDD_dev) {
  PbpolModelRajaHiop *pbpolrajahiop =
      reinterpret_cast<PbpolModelRajaHiop *>(opflow->model);
  LINEParamsRajaHiop *lineparams = &pbpolrajahiop->lineparams;
  int nxsparse = opflow->nxsparse;
  const int nxdense = opflow->nxdense;
  PetscErrorCode ierr;
  PetscInt flps = 0;

  PetscFunctionBegin;
  //  PetscPrintf(MPI_COMM_SELF,"Enter inequality dense hessian\n");

  /* Return if there are no inequality constraints */
  if (!lineparams->nlinelim) {
    //    PetscPrintf(MPI_COMM_SELF,"No inequality constraints. Exit inequality
    //    dense hessian\n");
    PetscFunctionReturn(0);
  }

  /** Although we will perform atomic operations on this view, we do not want
   * to use `RAJA::make_atomic_view` because explicit `RAJA::atomicAdd`s make
   * the atomic requirement clearer in the code */
  RAJA::View<double, RAJA::Layout<2>> HDD_view(HDD_dev, nxdense, nxdense);

  // Hessian from line contributions
  double *Gff_arr = lineparams->Gff_dev_;
  double *Gtt_arr = lineparams->Gtt_dev_;
  double *Gft_arr = lineparams->Gft_dev_;
  double *Gtf_arr = lineparams->Gtf_dev_;

  double *Bff_arr = lineparams->Bff_dev_;
  double *Btt_arr = lineparams->Btt_dev_;
  double *Bft_arr = lineparams->Bft_dev_;
  double *Btf_arr = lineparams->Btf_dev_;

  int *xidxf = lineparams->xidxf_dev_;
  int *xidxt = lineparams->xidxt_dev_;
  int *gineqidx = lineparams->gineqidx_dev_;
  int *linelimidx = lineparams->linelimidx_dev_;

  RAJA::forall<exago_raja_exec>(
      RAJA::RangeSegment(0, lineparams->nlinelim),
      RAJA_LAMBDA(RAJA::Index_type i) {
        int j = linelimidx[i];
        int gloc;
        int row[2], col[4];
        double val[8];

        double Pf, Qf, Pt, Qt;
        double thetaf = x_dev[xidxf[j]], Vmf = x_dev[xidxf[j] + 1];
        double thetat = x_dev[xidxt[j]], Vmt = x_dev[xidxt[j] + 1];
        double thetaft = thetaf - thetat;
        double thetatf = thetat - thetaf;
        double Gff = Gff_arr[j], Bff = Bff_arr[j];
        double Gft = Gft_arr[j], Bft = Bft_arr[j];
        double Gtf = Gtf_arr[j], Btf = Btf_arr[j];
        double Gtt = Gtt_arr[j], Btt = Btt_arr[j];

        Pf = Gff * Vmf * Vmf +
             Vmf * Vmt * (Gft * cos(thetaft) + Bft * sin(thetaft));
        Qf = -Bff * Vmf * Vmf +
             Vmf * Vmt * (-Bft * cos(thetaft) + Gft * sin(thetaft));

        Pt = Gtt * Vmt * Vmt +
             Vmt * Vmf * (Gtf * cos(thetatf) + Btf * sin(thetatf));
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
        dPf_dVmf =
            2. * Gff * Vmf + Vmt * (Gft * cos(thetaft) + Bft * sin(thetaft));
        dPf_dthetat = Vmf * Vmt * (Gft * sin(thetaft) - Bft * cos(thetaft));
        dPf_dVmt = Vmf * (Gft * cos(thetaft) + Bft * sin(thetaft));

        dQf_dthetaf = Vmf * Vmt * (Bft * sin(thetaft) + Gft * cos(thetaft));
        dQf_dVmf =
            -2. * Bff * Vmf + Vmt * (-Bft * cos(thetaft) + Gft * sin(thetaft));
        dQf_dthetat = Vmf * Vmt * (-Bft * sin(thetaft) - Gft * cos(thetaft));
        dQf_dVmt = Vmf * (-Bft * cos(thetaft) + Gft * sin(thetaft));

        dPt_dthetat = Vmt * Vmf * (-Gtf * sin(thetatf) + Btf * cos(thetatf));
        dPt_dVmt =
            2. * Gtt * Vmt + Vmf * (Gtf * cos(thetatf) + Btf * sin(thetatf));
        dPt_dthetaf = Vmt * Vmf * (Gtf * sin(thetatf) - Btf * cos(thetatf));
        dPt_dVmf = Vmt * (Gtf * cos(thetatf) + Btf * sin(thetatf));

        dQt_dthetat = Vmt * Vmf * (Btf * sin(thetatf) + Gtf * cos(thetatf));
        dQt_dVmt =
            -2. * Btt * Vmt + Vmf * (-Btf * cos(thetatf) + Gtf * sin(thetatf));
        dQt_dthetaf = Vmt * Vmf * (-Btf * sin(thetatf) - Gtf * cos(thetatf));
        dQt_dVmf = Vmt * (-Btf * cos(thetatf) + Gtf * sin(thetatf));

        double d2Pf_dthetaf_dthetaf, d2Pf_dthetaf_dVmf, d2Pf_dthetaf_dthetat,
            d2Pf_dthetaf_dVmt;
        double d2Pf_dVmf_dthetaf, d2Pf_dVmf_dVmf, d2Pf_dVmf_dthetat,
            d2Pf_dVmf_dVmt;
        double d2Pf_dthetat_dthetaf, d2Pf_dthetat_dVmf, d2Pf_dthetat_dthetat,
            d2Pf_dthetat_dVmt;
        double d2Pf_dVmt_dthetaf, d2Pf_dVmt_dVmf, d2Pf_dVmt_dthetat,
            d2Pf_dVmt_dVmt;

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
        double d2Qf_dVmf_dthetaf, d2Qf_dVmf_dVmf, d2Qf_dVmf_dthetat,
            d2Qf_dVmf_dVmt;
        double d2Qf_dthetat_dthetaf, d2Qf_dthetat_dVmf, d2Qf_dthetat_dthetat,
            d2Qf_dthetat_dVmt;
        double d2Qf_dVmt_dthetaf, d2Qf_dVmt_dVmf, d2Qf_dVmt_dthetat,
            d2Qf_dVmt_dVmt;

        /* dQf_dthetaf = Vmf*Vmt*(Bft*sin(thetaft) + Gft*cos(thetaft)); */
        d2Qf_dthetaf_dthetaf =
            Vmf * Vmt * (Bft * cos(thetaft) - Gft * sin(thetaft));
        d2Qf_dthetaf_dVmf = Vmt * (Bft * sin(thetaft) + Gft * cos(thetaft));
        d2Qf_dthetaf_dthetat =
            Vmf * Vmt * (-Bft * cos(thetaft) + Gft * sin(thetaft));
        d2Qf_dthetaf_dVmt = Vmf * (Bft * sin(thetaft) + Gft * cos(thetaft));

        /* dQf_dVmf = -2*Bff*Vmf + Vmt*(-Bft*cos(thetaft) + Gft*sin(thetaft));
         */
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
        double d2Pt_dVmt_dthetat, d2Pt_dVmt_dVmt, d2Pt_dVmt_dthetaf,
            d2Pt_dVmt_dVmf;
        double d2Pt_dthetaf_dthetat, d2Pt_dthetaf_dVmt, d2Pt_dthetaf_dthetaf,
            d2Pt_dthetaf_dVmf;
        double d2Pt_dVmf_dthetat, d2Pt_dVmf_dVmt, d2Pt_dVmf_dthetaf,
            d2Pt_dVmf_dVmf;

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
        double d2Qt_dVmf_dthetaf, d2Qt_dVmf_dVmf, d2Qt_dVmf_dthetat,
            d2Qt_dVmf_dVmt;
        double d2Qt_dthetat_dthetaf, d2Qt_dthetat_dVmf, d2Qt_dthetat_dthetat,
            d2Qt_dthetat_dVmt;
        double d2Qt_dVmt_dthetaf, d2Qt_dVmt_dVmf, d2Qt_dVmt_dthetat,
            d2Qt_dVmt_dVmt;

        /* dQt_dthetat = Vmf*Vmt*(Btf*sin(thetatf) + Gtf*cos(thetatf)); */
        d2Qt_dthetat_dthetat =
            Vmf * Vmt * (Btf * cos(thetatf) - Gtf * sin(thetatf));
        d2Qt_dthetat_dVmt = Vmf * (Btf * sin(thetatf) + Gtf * cos(thetatf));
        d2Qt_dthetat_dthetaf =
            Vmf * Vmt * (-Btf * cos(thetatf) + Gtf * sin(thetatf));
        d2Qt_dthetat_dVmf = Vmf * (Btf * sin(thetatf) + Gtf * cos(thetatf));

        /* dQt_dVmt = -2*Btt*Vmt + Vmf*(-Btf*cos(thetatf) + Gtf*sin(thetatf));
         */
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

        row[0] = xidxf[j] - nxsparse;
        col[0] = xidxf[j] - nxsparse;
        col[1] = xidxf[j] + 1 - nxsparse;
        col[2] = xidxt[j] - nxsparse;
        col[3] = xidxt[j] + 1 - nxsparse;

        gloc = gineqidx[i];

        val[0] = lambda_dev[gloc] * d2Sf2_dthetaf_dthetaf +
                 lambda_dev[gloc + 1] * d2St2_dthetaf_dthetaf;
        val[1] = lambda_dev[gloc] * d2Sf2_dthetaf_dVmf +
                 lambda_dev[gloc + 1] * d2St2_dthetaf_dVmf;
        val[2] = lambda_dev[gloc] * d2Sf2_dthetaf_dthetat +
                 lambda_dev[gloc + 1] * d2St2_dthetaf_dthetat;
        val[3] = lambda_dev[gloc] * d2Sf2_dthetaf_dVmt +
                 lambda_dev[gloc + 1] * d2St2_dthetaf_dVmt;

        RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[0]), val[0]);
        RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[1]), val[1]);
        RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[2]), val[2]);
        RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[3]), val[3]);

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
        col[0] = xidxf[j] - nxsparse;
        col[1] = xidxf[j] + 1 - nxsparse;
        col[2] = xidxt[j] - nxsparse;
        col[3] = xidxt[j] + 1 - nxsparse;

        row[0] = xidxf[j] + 1 - nxsparse;

        val[0] = lambda_dev[gloc] * d2Sf2_dVmf_dthetaf +
                 lambda_dev[gloc + 1] * d2St2_dVmf_dthetaf;
        val[1] = lambda_dev[gloc] * d2Sf2_dVmf_dVmf +
                 lambda_dev[gloc + 1] * d2St2_dVmf_dVmf;
        val[2] = lambda_dev[gloc] * d2Sf2_dVmf_dthetat +
                 lambda_dev[gloc + 1] * d2St2_dVmf_dthetat;
        val[3] = lambda_dev[gloc] * d2Sf2_dVmf_dVmt +
                 lambda_dev[gloc + 1] * d2St2_dVmf_dVmt;

        RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[0]), val[0]);
        RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[1]), val[1]);
        RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[2]), val[2]);
        RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[3]), val[3]);

        //    ierr = MatSetValues(H,1,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);

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

        col[0] = xidxf[j] - nxsparse;
        col[1] = xidxf[j] + 1 - nxsparse;
        col[2] = xidxt[j] - nxsparse;
        col[3] = xidxt[j] + 1 - nxsparse;

        row[0] = xidxt[j] - nxsparse;

        val[0] = lambda_dev[gloc] * d2Sf2_dthetat_dthetaf +
                 lambda_dev[gloc + 1] * d2St2_dthetat_dthetaf;
        val[1] = lambda_dev[gloc] * d2Sf2_dthetat_dVmf +
                 lambda_dev[gloc + 1] * d2St2_dthetat_dVmf;
        val[2] = lambda_dev[gloc] * d2Sf2_dthetat_dthetat +
                 lambda_dev[gloc + 1] * d2St2_dthetat_dthetat;
        val[3] = lambda_dev[gloc] * d2Sf2_dthetat_dVmt +
                 lambda_dev[gloc + 1] * d2St2_dthetat_dVmt;

        RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[0]), val[0]);
        RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[1]), val[1]);
        RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[2]), val[2]);
        RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[3]), val[3]);

        //    ierr = MatSetValues(H,1,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);

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

        row[0] = xidxt[j] + 1 - nxsparse;
        col[0] = xidxf[j] - nxsparse;
        col[1] = xidxf[j] + 1 - nxsparse;
        col[2] = xidxt[j] - nxsparse;
        col[3] = xidxt[j] + 1 - nxsparse;

        val[0] = lambda_dev[gloc] * d2Sf2_dVmt_dthetaf +
                 lambda_dev[gloc + 1] * d2St2_dVmt_dthetaf;
        val[1] = lambda_dev[gloc] * d2Sf2_dVmt_dVmf +
                 lambda_dev[gloc + 1] * d2St2_dVmt_dVmf;
        val[2] = lambda_dev[gloc] * d2Sf2_dVmt_dthetat +
                 lambda_dev[gloc + 1] * d2St2_dVmt_dthetat;
        val[3] = lambda_dev[gloc] * d2Sf2_dVmt_dVmt +
                 lambda_dev[gloc + 1] * d2St2_dVmt_dVmt;

        RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[0]), val[0]);
        RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[1]), val[1]);
        RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[2]), val[2]);
        RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[3]), val[3]);
      });
  //  PetscPrintf(MPI_COMM_SELF,"Exit inequality dense hessian\n");
  flps += (972 + (92 * EXAGO_FLOPS_COSOP) + (92 * EXAGO_FLOPS_SINOP)) *
          lineparams->nlinelim;
  ierr = PetscLogFlops(flps);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeDenseHessian_PBPOLRAJAHIOP(OPFLOW opflow,
                                                       const double *x_dev,
                                                       const double *lambda_dev,
                                                       double *HDD_dev) {
  PetscErrorCode ierr;
  int nxdense = opflow->nxdense;

  if (!HDD_dev)
    PetscFunctionReturn(0);

  RAJA::forall<exago_raja_exec>(
      RAJA::RangeSegment(0, nxdense * nxdense),
      RAJA_LAMBDA(RAJA::Index_type i) { HDD_dev[i] = 0.0; });

  /* Equality constraint Hessian */
  ierr = OPFLOWComputeDenseEqualityConstraintHessian_PBPOLRAJAHIOP(
      opflow, x_dev, lambda_dev, HDD_dev);
  CHKERRQ(ierr);

  if (opflow->nconineq) {
    ierr = OPFLOWComputeDenseInequalityConstraintHessian_PBPOLRAJAHIOP(
        opflow, x_dev, lambda_dev + opflow->nconeq, HDD_dev);
    CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}
#endif
