
#include <exago_config.h>

#if defined(EXAGO_HAVE_RAJA) 

#include <umpire/Allocator.hpp>
#include <umpire/ResourceManager.hpp>

#include <RAJA/RAJA.hpp>

#include <private/psimpl.h>
#include <private/opflowimpl.h>
#include "pbpolrajahiopkernels.hpp"
#include "pbpolrajahiop.hpp"

#ifdef EXAGO_HAVE_GPU
  using exago_raja_exec = RAJA::cuda_exec<128>;
  using exago_raja_reduce = RAJA::cuda_reduce;
  using exago_raja_atomic = RAJA::cuda_atomic;
  #define RAJA_LAMBDA [=] __device__
#else
  using exago_raja_exec = RAJA::omp_parallel_for_exec;
  using exago_raja_reduce = RAJA::omp_reduce;
  using exago_raja_atomic = RAJA::omp_atomic;
  #define RAJA_LAMBDA [=]
#endif

/************* NOTE ***********************/
/* No Load loss or power imbalance variables considered yet */
/********************************************/

/* Functions to create and destroy data arrays for different
   component classes
*/
int BUSParamsRajaHiop::destroy(OPFLOW opflow)
{
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

#ifdef EXAGO_HAVE_GPU
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
#endif

  return 0;
}

/* Create data for buses that is used in different computations */
int BUSParamsRajaHiop::allocate(OPFLOW opflow)
{
  PS             ps=opflow->ps;
  PetscInt       loc,gloc=0;
  PSBUS          bus;
  PetscInt       i;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  
  nbus = ps->nbus;

  /* Allocate the arrays */
  auto& resmgr = umpire::ResourceManager::getInstance();
  h_allocator_ = resmgr.getAllocator("HOST");

  // Allocate data on the host
  isref      = paramAlloc<int>(h_allocator_, nbus);
  isisolated = paramAlloc<int>(h_allocator_, nbus);
  ispvpq     = paramAlloc<int>(h_allocator_, nbus);

  vmin = paramAlloc<double>(h_allocator_, nbus);
  vmax = paramAlloc<double>(h_allocator_, nbus);
  va   = paramAlloc<double>(h_allocator_, nbus);
  vm   = paramAlloc<double>(h_allocator_, nbus);
  gl   = paramAlloc<double>(h_allocator_, nbus);
  bl   = paramAlloc<double>(h_allocator_, nbus);

  xidx = paramAlloc<int>(h_allocator_, nbus);
  gidx = paramAlloc<int>(h_allocator_, nbus);

  /* Memzero arrays */
  resmgr.memset(isref,0,nbus*sizeof(int));
  resmgr.memset(ispvpq,0,nbus*sizeof(int));
  resmgr.memset(isisolated,0,nbus*sizeof(int));
  
  for (int i = 0; i < nbus; i++)
  {
    vmin[i] = 0.9;
    vmax[i] = 1.1;
  }
  
  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
    
    xidx[i] = opflow->idxn2sd_map[loc];
    gidx[i] = gloc;

    if(bus->ide == REF_BUS) isref[i] = 1;
    else if(bus->ide == ISOLATED_BUS) isisolated[i] = 1;
    else ispvpq[i] = 1;

    if(opflow->genbusVmfixed) {
      if(bus->ide == REF_BUS || bus->ide == PV_BUS) {
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
    vm[i]   = bus->vm;
    va[i]   = bus->va;
    gl[i]   = bus->gl;
    bl[i]   = bus->bl;

    gloc += 2;
  }

#ifdef EXAGO_HAVE_GPU
  d_allocator_ = resmgr.getAllocator("DEVICE");
  // Allocate data on the device
  isref_dev_      = paramAlloc<int>(d_allocator_, nbus);
  isisolated_dev_ = paramAlloc<int>(d_allocator_, nbus);
  ispvpq_dev_     = paramAlloc<int>(d_allocator_, nbus);

  vmin_dev_ = paramAlloc<double>(d_allocator_, nbus);
  vmax_dev_ = paramAlloc<double>(d_allocator_, nbus);
  va_dev_   = paramAlloc<double>(d_allocator_, nbus);
  vm_dev_   = paramAlloc<double>(d_allocator_, nbus);
  gl_dev_   = paramAlloc<double>(d_allocator_, nbus);
  bl_dev_   = paramAlloc<double>(d_allocator_, nbus);

  xidx_dev_ = paramAlloc<int>(d_allocator_, nbus);
  gidx_dev_ = paramAlloc<int>(d_allocator_, nbus);

  // Copy host data from the host to the device
  resmgr.copy(isref_dev_, isref);
  resmgr.copy(isisolated_dev_, isisolated);
  resmgr.copy(ispvpq_dev_, ispvpq);

  resmgr.copy(vmin_dev_, vmin);
  resmgr.copy(vmax_dev_, vmax);
  resmgr.copy(va_dev_  , va);
  resmgr.copy(vm_dev_  , vm);
  resmgr.copy(gl_dev_  , gl);
  resmgr.copy(bl_dev_  , bl);

  resmgr.copy(xidx_dev_, xidx);
  resmgr.copy(gidx_dev_, gidx);
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
  gidx_dev_ = gidx;
#endif
  PetscFunctionReturn(0);
}


int LINEParamsRajaHiop::destroy(OPFLOW opflow)
{
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
  h_allocator_.deallocate(gineqidx);
  h_allocator_.deallocate(gbineqidx);

  if(opflow->nconineq) {
    h_allocator_.deallocate(linelimidx);
  }

#ifdef EXAGO_HAVE_GPU
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
  d_allocator_.deallocate(gineqidx_dev_);
  d_allocator_.deallocate(gbineqidx_dev_);
  if(opflow->nconineq) {
    d_allocator_.deallocate(linelimidx_dev_);
  }
#endif

  return 0;
}


/* Create data for lines that is used in different computations */
int LINEParamsRajaHiop::allocate(OPFLOW opflow)
{
  PS             ps=opflow->ps;
  PetscInt       linei=0,linelimi=0;
  PSLINE         line;
  PetscInt       i;
  PetscErrorCode ierr;
  const PSBUS    *connbuses;
  PSBUS          busf,bust;
  PetscInt       gloc=0; /* offset for inequality constraint contributions */
  PetscInt       gbloc=opflow->nconeq; /* starting offset for inequality constraint bound */
  /* the above gloc, gbloc is needed because for the constraint bound calculation,
     the entire G vector is passed in, while for inequality constraints calulation only
     the inequality constraint vector is passed in 
  */

  PetscFunctionBegin;

  ierr = PSGetNumActiveLines(ps, &nlineON, NULL); CHKERRQ(ierr);

  nlinelim = 0;
  /* Get the number of lines that are active and have finite limits. These lines
     will be only considered in inequality constraints */
  if(opflow->nconineq) {
    for(i=0; i < ps->nline; i++) {
      line = &ps->line[i];

      if(!line->status || line->rateA > 1e5) continue;
      nlinelim++;
    }
  }

  /* Allocate data arrays */
  auto& resmgr = umpire::ResourceManager::getInstance();
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
  gineqidx = paramAlloc<int>(h_allocator_, nlineON);
  gbineqidx = paramAlloc<int>(h_allocator_, nlineON);

  if(opflow->nconineq) {
    linelimidx = paramAlloc<int>(h_allocator_, nlinelim);
  }

  /* Populate arrays */
  for(i=0; i < ps->nline; i++) {
    line = &ps->line[i];

    if(!line->status) continue;

    Gff[linei] = line->yff[0];
    Bff[linei] = line->yff[1];
    Gft[linei] = line->yft[0];
    Bft[linei] = line->yft[1];
    Gtf[linei] = line->ytf[0];
    Btf[linei] = line->ytf[1];
    Gtt[linei] = line->ytt[0];
    Btt[linei] = line->ytt[1];
    rateA[linei] = line->rateA;

    ierr = PSLINEGetConnectedBuses(line,&connbuses);CHKERRQ(ierr);
    busf = connbuses[0];
    bust = connbuses[1];

    int xidxfi,xidxti;
    ierr = PSBUSGetVariableLocation(busf,&xidxfi);CHKERRQ(ierr);
    ierr = PSBUSGetVariableLocation(bust,&xidxti);CHKERRQ(ierr);

    xidxf[linei] = opflow->idxn2sd_map[xidxfi];
    xidxt[linei] = opflow->idxn2sd_map[xidxti];

    /* 
       Each bus has two equality (balance) constraints, hence the use of coefficient 2
       to map the location of the equality constraint for the bus
    */
    geqidxf[linei] = 2*busf->internal_i;
    geqidxt[linei] = 2*bust->internal_i;
    
    if(opflow->nconineq) {
      if(line->status && line->rateA < 1e5) {
	gbineqidx[linelimi] = gbloc;
	gineqidx[linelimi] = gloc;
	linelimidx[linelimi] = linei;
	linelimi++;
	gbloc += 2;
	gloc += 2;
      }
    }
    linei++;
  }

#ifdef EXAGO_HAVE_GPU
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
  gineqidx_dev_  = paramAlloc<int>(d_allocator_, nlineON);
  gbineqidx_dev_ = paramAlloc<int>(d_allocator_, nlineON);

  if(opflow->nconineq) {
    linelimidx_dev_ = paramAlloc<int>(d_allocator_, nlinelim);
  }

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
  resmgr.copy(gineqidx_dev_, gineqidx);
  resmgr.copy(gbineqidx_dev_, gbineqidx);

  if(opflow->nconineq) {
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
  gineqidx_dev_ = gineqidx;
  gbineqidx_dev_ = gbineqidx;
  if(opflow->nconineq)
  {
    linelimidx_dev_ = linelimidx;
  }
#endif
  return 0;
}

int LOADParamsRajaHiop::destroy(OPFLOW opflow)
{
  h_allocator_.deallocate(pl);
  h_allocator_.deallocate(ql);
  h_allocator_.deallocate(xidx);
  h_allocator_.deallocate(gidx);
#ifdef EXAGO_HAVE_GPU
  d_allocator_.deallocate(pl_dev_);
  d_allocator_.deallocate(ql_dev_);
  d_allocator_.deallocate(xidx_dev_);
  d_allocator_.deallocate(gidx_dev_);
#endif
  return 0;
}

/* Create data for loads that is used in different computations */
int LOADParamsRajaHiop::allocate(OPFLOW opflow)
{
  PS             ps=opflow->ps;
  PetscInt       loc,gloc=0,loadi=0;
  PSLOAD         load;
  PSBUS          bus;
  PetscInt       i,j;
  PetscErrorCode ierr;
  
  /* Get the number of active generators (STATUS ON) */
  ierr = PSGetNumLoads(ps, &nload, NULL); CHKERRQ(ierr);

  /* Allocate arrays */
  auto& resmgr = umpire::ResourceManager::getInstance();
  h_allocator_ = resmgr.getAllocator("HOST");

  // Allocate data on the host
  pl   = paramAlloc<double>(h_allocator_, nload);
  ql   = paramAlloc<double>(h_allocator_, nload);
  xidx = paramAlloc<int>(h_allocator_, nload);
  gidx = paramAlloc<int>(h_allocator_, nload);

  /* Insert data in loadparams */
  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
    for(j=0; j < bus->nload; j++) {
      ierr = PSBUSGetLoad(bus,j,&load);CHKERRQ(ierr);
      
      loc += 2;

      pl[loadi] = load->pl;
      ql[loadi] = load->ql;
      xidx[loadi] = loc;
      gidx[loadi] = gloc;

      loadi++;
    }
    gloc += 2;
  }

#ifdef EXAGO_HAVE_GPU
  d_allocator_ = resmgr.getAllocator("DEVICE");
  // Allocate data on the device
  pl_dev_   = paramAlloc<double>(d_allocator_, nload);
  ql_dev_   = paramAlloc<double>(d_allocator_, nload);
  xidx_dev_ = paramAlloc<int>(d_allocator_, nload);
  gidx_dev_ = paramAlloc<int>(d_allocator_, nload);

  // Copy data from host to device
  resmgr.copy(pl_dev_, pl);
  resmgr.copy(ql_dev_, ql);
  resmgr.copy(xidx_dev_, xidx);
  resmgr.copy(gidx_dev_, gidx);
#else
  pl_dev_ = pl;
  ql_dev_ = ql;
  xidx_dev_ = xidx;
  gidx_dev_ = gidx;
#endif

  return (0);
}


int GENParamsRajaHiop::destroy(OPFLOW opflow)
{
  // Free arrays on the host
  h_allocator_.deallocate(cost_alpha);
  h_allocator_.deallocate(cost_beta);
  h_allocator_.deallocate(cost_gamma);
  h_allocator_.deallocate(pt);
  h_allocator_.deallocate(pb);
  h_allocator_.deallocate(qt);
  h_allocator_.deallocate(qb);
  h_allocator_.deallocate(xidx);
  h_allocator_.deallocate(gidx);
  h_allocator_.deallocate(jacsp_idx);
  h_allocator_.deallocate(jacsq_idx);
#ifdef EXAGO_HAVE_GPU
  // Free arrays on the device
  d_allocator_.deallocate(cost_alpha_dev_);
  d_allocator_.deallocate(cost_beta_dev_);
  d_allocator_.deallocate(cost_gamma_dev_);
  d_allocator_.deallocate(pt_dev_);
  d_allocator_.deallocate(pb_dev_);
  d_allocator_.deallocate(qt_dev_);
  d_allocator_.deallocate(qb_dev_);
  d_allocator_.deallocate(xidx_dev_);
  d_allocator_.deallocate(gidx_dev_);
  d_allocator_.deallocate(jacsp_idx_dev_);
  d_allocator_.deallocate(jacsq_idx_dev_);
#endif
  return 0;
}


/* Create data for generators that is used in different computations */
int GENParamsRajaHiop::allocate(OPFLOW opflow)
{
  PS             ps=opflow->ps;
  PetscInt       loc,gloc=0,geni=0,nnzs=0,gi;
  PSGEN          gen;
  PSBUS          bus;
  PetscInt       i,j;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  
  /* Get the number of active generators (STATUS ON) */
  ierr = PSGetNumActiveGenerators(ps,&ngenON,NULL);CHKERRQ(ierr);

  /* Allocate arrays on the host */
  auto& resmgr = umpire::ResourceManager::getInstance();
  h_allocator_ = resmgr.getAllocator("HOST");

  cost_alpha = paramAlloc<double>(h_allocator_, ngenON);
  cost_beta  = paramAlloc<double>(h_allocator_, ngenON);
  cost_gamma = paramAlloc<double>(h_allocator_, ngenON);

  pt = paramAlloc<double>(h_allocator_, ngenON);
  pb = paramAlloc<double>(h_allocator_, ngenON);
  qt = paramAlloc<double>(h_allocator_, ngenON);
  qb = paramAlloc<double>(h_allocator_, ngenON);

  xidx = paramAlloc<int>(h_allocator_, ngenON);
  gidx = paramAlloc<int>(h_allocator_, ngenON);

  jacsp_idx = paramAlloc<int>(h_allocator_,ngenON);
  jacsq_idx = paramAlloc<int>(h_allocator_,ngenON);

  /* Populate data on the host */
  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
    gi = 0;
    for(j=0; j < bus->ngen; j++) {
      ierr = PSBUSGetGen(bus,j,&gen);CHKERRQ(ierr);
      if(!gen->status) continue;
      
      loc += 2;

      cost_alpha[geni] = gen->cost_alpha;
      cost_beta[geni]  = gen->cost_beta;
      cost_gamma[geni] = gen->cost_gamma;
      pt[geni]         = gen->pt;
      pb[geni]         = gen->pb;
      qt[geni]         = gen->qt;
      qb[geni]         = gen->qb;

      xidx[geni]       = opflow->idxn2sd_map[loc];
      gidx[geni]       = gloc;
      jacsp_idx[geni]  = nnzs + gi;
      jacsq_idx[geni]  = nnzs + bus->ngenON + gi;

      geni++;
      gi++;
    }
    nnzs += 2*bus->ngenON;
    gloc += 2;
  }

#ifdef EXAGO_HAVE_GPU
  d_allocator_ = resmgr.getAllocator("DEVICE");
  // Allocate arrays on the device
  cost_alpha_dev_ = paramAlloc<double>(d_allocator_, ngenON);
  cost_beta_dev_  = paramAlloc<double>(d_allocator_, ngenON);
  cost_gamma_dev_ = paramAlloc<double>(d_allocator_, ngenON);

  pt_dev_ = paramAlloc<double>(d_allocator_, ngenON);
  pb_dev_ = paramAlloc<double>(d_allocator_, ngenON);
  qt_dev_ = paramAlloc<double>(d_allocator_, ngenON);
  qb_dev_ = paramAlloc<double>(d_allocator_, ngenON);

  xidx_dev_ = paramAlloc<int>(d_allocator_, ngenON);
  gidx_dev_ = paramAlloc<int>(d_allocator_, ngenON);

  jacsp_idx_dev_ = paramAlloc<int>(d_allocator_,ngenON);
  jacsq_idx_dev_ = paramAlloc<int>(d_allocator_,ngenON);

  // Copy host data to the device
  resmgr.copy(cost_alpha_dev_, cost_alpha);
  resmgr.copy(cost_beta_dev_ , cost_beta );
  resmgr.copy(cost_gamma_dev_, cost_gamma);

  resmgr.copy(pt_dev_, pt);
  resmgr.copy(pb_dev_, pb);
  resmgr.copy(qt_dev_, qt);
  resmgr.copy(qb_dev_, qb);

  resmgr.copy(xidx_dev_, xidx);
  resmgr.copy(gidx_dev_, gidx);

  resmgr.copy(jacsp_idx_dev_,jacsp_idx);
  resmgr.copy(jacsq_idx_dev_,jacsq_idx);
#else
  cost_alpha_dev_ = cost_alpha;
  cost_beta_dev_ = cost_beta;
  cost_gamma_dev_ = cost_gamma;
  pt_dev_ = pt;
  pb_dev_ = pb;
  qt_dev_ = qt;
  qb_dev_ = qb;
  xidx_dev_ = xidx;
  gidx_dev_ = gidx;
  jacsp_idx_dev_ = jacsp_idx;
  jacsq_idx_dev_ = jacsq_idx;
#endif
  return 0;
}


PetscErrorCode OPFLOWSetInitialGuessArray_PBPOLRAJAHIOP(OPFLOW opflow,double* x0_dev)
{
  PetscErrorCode ierr;
  double         *x;
  auto& resmgr = umpire::ResourceManager::getInstance();

  PetscFunctionBegin;
  //  ierr = PetscPrintf(MPI_COMM_SELF,"Entered OPFLOWInitialization\n");CHKERRQ(ierr);

  // Do initialization on Host
  ierr = (*opflow->modelops.setinitialguess)(opflow,opflow->X);CHKERRQ(ierr);
  ierr = VecGetArray(opflow->X,&x);CHKERRQ(ierr);

  // Copy from host to device
  umpire::Allocator h_allocator_ = resmgr.getAllocator("HOST");
  registerWith(x,opflow->nx,resmgr,h_allocator_);
  resmgr.copy(x0_dev,x);

  ierr = VecRestoreArray(opflow->X,&x);CHKERRQ(ierr);

  //  ierr = PetscPrintf(MPI_COMM_SELF,"Exit Initialization\n");CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


PetscErrorCode OPFLOWSetConstraintBoundsArray_PBPOLRAJAHIOP(OPFLOW opflow,double *gl_dev,double *gu_dev)
{
  PbpolModelRajaHiop* pbpolrajahiop = reinterpret_cast<PbpolModelRajaHiop*>(opflow->model);
  LINEParamsRajaHiop     *lineparams=&pbpolrajahiop->lineparams;
  PS             ps=opflow->ps;
  double         MVAbase = ps->MVAbase;

  PetscFunctionBegin;

  //  PetscPrintf(MPI_COMM_SELF,"Entered Constraint Bounds\n");

  /* Equallity constraints (all zeros) */
  auto& resmgr = umpire::ResourceManager::getInstance();
  resmgr.memset(gl_dev, 0);
  resmgr.memset(gu_dev, 0);

  /* Inequality constraint bounds */
  int* linelimidx = lineparams->linelimidx_dev_;
  int* gbineqidx = lineparams->gbineqidx_dev_;
  double* rateA = lineparams->rateA_dev_;

  if(lineparams->nlinelim) {
    // PetscPrintf(PETSC_COMM_SELF,"nlinelim = %d ineq = %d\n",lineparams->nlinelim,opflow->nconineq);
    RAJA::forall<exago_raja_exec>(RAJA::RangeSegment(0, lineparams->nlinelim),
    RAJA_LAMBDA (RAJA::Index_type i)
    {
      int    j=linelimidx[i];
      gl_dev[gbineqidx[i]]   = 0.0;
      gu_dev[gbineqidx[i]]   = (rateA[j]/MVAbase)*(rateA[j]/MVAbase);
      gl_dev[gbineqidx[i]+1] = 0.0;
      gu_dev[gbineqidx[i]+1] = (rateA[j]/MVAbase)*(rateA[j]/MVAbase);
    }
  );
  }

  //  PetscPrintf(MPI_COMM_SELF,"Exit Constraint Bounds\n");
  PetscFunctionReturn(0);
}


/** EQUALITY CONSTRAINTS */
PetscErrorCode OPFLOWComputeEqualityConstraintsArray_PBPOLRAJAHIOP(OPFLOW opflow,const double *x_dev, double *ge_dev)
{
  //  PBPOLRAJAHIOP       pbpolrajahiop=(PBPOLRAJAHIOP)opflow->model;
  PbpolModelRajaHiop* pbpolrajahiop = reinterpret_cast<PbpolModelRajaHiop*>(opflow->model);
  BUSParamsRajaHiop      *busparams=&pbpolrajahiop->busparams;
  GENParamsRajaHiop      *genparams=&pbpolrajahiop->genparams;
  LOADParamsRajaHiop     *loadparams=&pbpolrajahiop->loadparams;
  LINEParamsRajaHiop     *lineparams=&pbpolrajahiop->lineparams;
  PetscErrorCode ierr;
  double                 flps=0.0;
  PetscFunctionBegin;
  //  PetscPrintf(MPI_COMM_SELF,"Entered Equality constraints\n");

  // Zero out array
  auto& resmgr = umpire::ResourceManager::getInstance();
  resmgr.memset(ge_dev, 0);

  /* Generator contributions */
  int* g_gidx = genparams->gidx_dev_;
  int* g_xidx = genparams->xidx_dev_;
  RAJA::forall< exago_raja_exec >( RAJA::RangeSegment(0, genparams->ngenON),
    RAJA_LAMBDA (RAJA::Index_type i) 
    {
      RAJA::atomicSub<exago_raja_atomic>(&ge_dev[g_gidx[i]],x_dev[g_xidx[i]]);
      RAJA::atomicSub<exago_raja_atomic>(&ge_dev[g_gidx[i]+1],x_dev[g_xidx[i]+1]);
    }
  );
  flps += genparams->ngenON*2.0;

  /* Load contributions */
  double* pl = loadparams->pl_dev_;
  double* ql = loadparams->ql_dev_;
  int*    l_gidx = loadparams->gidx_dev_;
  RAJA::forall< exago_raja_exec >( RAJA::RangeSegment(0, loadparams->nload),
    RAJA_LAMBDA (RAJA::Index_type i) 
    {
      RAJA::atomicAdd<exago_raja_atomic>(&ge_dev[l_gidx[i]],pl[i]);
      RAJA::atomicAdd<exago_raja_atomic>(&ge_dev[l_gidx[i]+1],ql[i]);
    }
  );
  flps += loadparams->nload*2.0;

  /* Bus contributions */
  int* isisolated = busparams->isisolated_dev_;
  int* ispvpq = busparams->ispvpq_dev_;
  double* gl = busparams->gl_dev_;
  double* bl = busparams->bl_dev_;
  double* vm = busparams->vm_dev_;
  double* va = busparams->va_dev_;
  int* b_xidx = busparams->xidx_dev_;
  int* b_gidx = busparams->gidx_dev_;
  RAJA::forall<exago_raja_exec>(RAJA::RangeSegment(0, busparams->nbus),
    RAJA_LAMBDA (RAJA::Index_type i)
    {
      double theta= x_dev[b_xidx[i]];
      double Vm   = x_dev[b_xidx[i]+1];
      RAJA::atomicAdd<exago_raja_atomic>(&ge_dev[b_gidx[i]],
				      isisolated[i]*(theta - va[i]*PETSC_PI/180.0) + ispvpq[i]*Vm*Vm*gl[i]);

      RAJA::atomicAdd<exago_raja_atomic>(&ge_dev[b_gidx[i]+1],
				       isisolated[i]*(Vm    - vm[i]) - ispvpq[i]*Vm*Vm*bl[i]);
    }
  );
  flps += busparams->nbus*14.0;

  /* Line contributions */
  double* Gff = lineparams->Gff_dev_;
  double* Gtt = lineparams->Gtt_dev_;
  double* Gft = lineparams->Gft_dev_;
  double* Gtf = lineparams->Gtf_dev_;

  double* Bff = lineparams->Bff_dev_;
  double* Btt = lineparams->Btt_dev_;
  double* Bft = lineparams->Bft_dev_;
  double* Btf = lineparams->Btf_dev_;

  int* xidxf = lineparams->xidxf_dev_;
  int* xidxt = lineparams->xidxt_dev_;
  int* geqidxf = lineparams->geqidxf_dev_;
  int* geqidxt = lineparams->geqidxt_dev_;

  RAJA::forall<exago_raja_exec>(RAJA::RangeSegment(0, lineparams->nlineON),
    RAJA_LAMBDA (RAJA::Index_type i)
    {
      double Pf,Qf,Pt,Qt;
      double thetaf=x_dev[xidxf[i]], Vmf=x_dev[xidxf[i]+1];
      double thetat=x_dev[xidxt[i]], Vmt=x_dev[xidxt[i]+1];
      double thetaft=thetaf-thetat;
      double thetatf=thetat-thetaf;

      Pf = Gff[i]*Vmf*Vmf  + Vmf*Vmt*(Gft[i]*cos(thetaft) + Bft[i]*sin(thetaft));
      Qf = -Bff[i]*Vmf*Vmf + Vmf*Vmt*(-Bft[i]*cos(thetaft) + Gft[i]*sin(thetaft));
      Pt = Gtt[i]*Vmt*Vmt  + Vmt*Vmf*(Gtf[i]*cos(thetatf) + Btf[i]*sin(thetatf));
      Qt = -Btt[i]*Vmt*Vmt + Vmt*Vmf*(-Btf[i]*cos(thetatf) + Gtf[i]*sin(thetatf));
    
      RAJA::atomicAdd<exago_raja_atomic>(&ge_dev[geqidxf[i]]  , Pf);
      RAJA::atomicAdd<exago_raja_atomic>(&ge_dev[geqidxf[i]+1], Qf);
      RAJA::atomicAdd<exago_raja_atomic>(&ge_dev[geqidxt[i]]  , Pt);
      RAJA::atomicAdd<exago_raja_atomic>(&ge_dev[geqidxt[i]+1], Qt);
    }
  );
  flps += lineparams->nlineON*(46+(4*EXAGO_FLOPS_COSOP)+(4*EXAGO_FLOPS_SINOP));
  //  PetscPrintf(MPI_COMM_SELF,"Exit Equality Constraints\n");

  ierr = PetscLogFlops(flps);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/** INEQUALITY CONSTRAINTS **/
PetscErrorCode OPFLOWComputeInequalityConstraintsArray_PBPOLRAJAHIOP(OPFLOW opflow, const double *x_dev, double *gi_dev)
{
  PbpolModelRajaHiop* pbpolrajahiop = reinterpret_cast<PbpolModelRajaHiop*>(opflow->model);
  LINEParamsRajaHiop     *lineparams=&pbpolrajahiop->lineparams;
  double                 flps=0.0;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  //  PetscPrintf(MPI_COMM_SELF,"Entered Inequality Constraints\n");

  // Zero out array
  auto& resmgr = umpire::ResourceManager::getInstance();
  resmgr.memset(gi_dev, 0);

  /* Line contributions */
  double* Gff = lineparams->Gff_dev_;
  double* Gtt = lineparams->Gtt_dev_;
  double* Gft = lineparams->Gft_dev_;
  double* Gtf = lineparams->Gtf_dev_;

  double* Bff = lineparams->Bff_dev_;
  double* Btt = lineparams->Btt_dev_;
  double* Bft = lineparams->Bft_dev_;
  double* Btf = lineparams->Btf_dev_;

  int* linelimidx = lineparams->linelimidx_dev_;
  int* xidxf = lineparams->xidxf_dev_;
  int* xidxt = lineparams->xidxt_dev_;
  int* gineqidx = lineparams->gineqidx_dev_;
  if(lineparams->nlinelim) {
    RAJA::forall<exago_raja_exec>(RAJA::RangeSegment(0, lineparams->nlinelim),
      RAJA_LAMBDA (RAJA::Index_type i)
      {
	int    j=linelimidx[i];
	double Pf,Qf,Pt,Qt,Sf2,St2;
	double thetaf=x_dev[xidxf[j]], Vmf=x_dev[xidxf[j]+1];
	double thetat=x_dev[xidxt[j]], Vmt=x_dev[xidxt[j]+1];
	double thetaft=thetaf-thetat;
	double thetatf=thetat-thetaf;
    
	Pf = Gff[j]*Vmf*Vmf  + Vmf*Vmt*(Gft[j]*cos(thetaft) + Bft[j]*sin(thetaft));
	Qf = -Bff[j]*Vmf*Vmf + Vmf*Vmt*(-Bft[j]*cos(thetaft) + Gft[j]*sin(thetaft));
	Pt = Gtt[j]*Vmt*Vmt  + Vmt*Vmf*(Gtf[j]*cos(thetatf) + Btf[j]*sin(thetatf));
	Qt = -Btt[j]*Vmt*Vmt + Vmt*Vmf*(-Btf[j]*cos(thetatf) + Gtf[j]*sin(thetatf));
	
	Sf2 = Pf*Pf + Qf*Qf;
	St2 = Pt*Pt + Qt*Qt;
	
	gi_dev[gineqidx[i]]   = Sf2;
	gi_dev[gineqidx[i]+1] = St2;
      }
    );
    flps += lineparams->nlinelim*(72+(4*EXAGO_FLOPS_COSOP)+(4*EXAGO_FLOPS_SINOP));
  }
  //  PetscPrintf(MPI_COMM_SELF,"Exit Inequality Constraints\n");
  ierr = PetscLogFlops(flps);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/** OBJECTIVE FUNCTION **/
// Note: This kernel (and all the kernels for this model assume that the data has been already
// allocated on the device. x_dev is pointer to array on the GPU
PetscErrorCode OPFLOWComputeObjectiveArray_PBPOLRAJAHIOP(OPFLOW opflow,const double *x_dev,double *obj)
{
  //  PBPOLRAJAHIOP      pbpolrajahiop=(PBPOLRAJAHIOP)opflow->model;
  PbpolModelRajaHiop* pbpolrajahiop = reinterpret_cast<PbpolModelRajaHiop*>(opflow->model);
  GENParamsRajaHiop     *genparams=&pbpolrajahiop->genparams;
  PetscErrorCode ierr;
  PS             ps=opflow->ps;
  int            isobj_gencost=opflow->obj_gencost;
  double         MVAbase=ps->MVAbase;

  PetscFunctionBegin;

  //  PetscPrintf(MPI_COMM_SELF,"Entered objective function\n");

  // You need local copies of device pointers so that lambda can capture them
  double* cost_alpha = genparams->cost_alpha_dev_;
  double* cost_beta  = genparams->cost_beta_dev_;
  double* cost_gamma = genparams->cost_gamma_dev_;
  int*    xidx = genparams->xidx_dev_;

  /* Generator objective function contributions */
  // Set up reduce sum object
  RAJA::ReduceSum< exago_raja_reduce, double> obj_val_sum(0.0);
  // Compute reduction on CUDA device
  RAJA::forall< exago_raja_exec >( RAJA::RangeSegment(0, genparams->ngenON), 
    RAJA_LAMBDA (RAJA::Index_type i)
    {
      double Pg = x_dev[xidx[i]]*MVAbase;
      obj_val_sum += isobj_gencost*(cost_alpha[i]*Pg*Pg + cost_beta[i]*Pg + cost_gamma[i]);
    }
  );

  *obj = static_cast<double>(obj_val_sum.get());
  ierr = PetscLogFlops(genparams->ngenON*8.0);CHKERRQ(ierr);

  //  PetscPrintf(MPI_COMM_SELF,"Exit objective function\n");
  PetscFunctionReturn(0);
}

/** GRADIENT **/
PetscErrorCode OPFLOWComputeGradientArray_PBPOLRAJAHIOP(OPFLOW opflow,const double *x_dev, double* grad_dev)
{
  PbpolModelRajaHiop* pbpolrajahiop = reinterpret_cast<PbpolModelRajaHiop*>(opflow->model);
  GENParamsRajaHiop   *genparams=&pbpolrajahiop->genparams;
  PS                  ps=opflow->ps;
  int                 isobj_gencost=opflow->obj_gencost;
  double              MVAbase=ps->MVAbase;
  PetscErrorCode      ierr;

  PetscFunctionBegin;
  //  PetscPrintf(MPI_COMM_SELF,"Entered gradient function\n");

  // Zero out array
  auto& resmgr = umpire::ResourceManager::getInstance();
  resmgr.memset(grad_dev, 0);

  /* Generator gradient contributions */
  double* cost_alpha = genparams->cost_alpha_dev_;
  double* cost_beta  = genparams->cost_beta_dev_;
  int*    xidx       = genparams->xidx_dev_;
  RAJA::forall< exago_raja_exec >( RAJA::RangeSegment(0, genparams->ngenON), 
    RAJA_LAMBDA (RAJA::Index_type i)
    {
      double Pg = x_dev[xidx[i]]*MVAbase;
      grad_dev[xidx[i]] = isobj_gencost*MVAbase*(2.0*cost_alpha[i]*Pg + cost_beta[i]);
    }
  );

  ierr = PetscLogFlops(genparams->ngenON*6.0);CHKERRQ(ierr);
  //  PetscPrintf(MPI_COMM_SELF,"Exit gradient function\n");

  PetscFunctionReturn(0);
}

// Note: This kernel (and all the kernels for this model assume that the data has been already
// allocated on the device. xl_dev and xu_dev are pointers to arrays on the GPU
PetscErrorCode OPFLOWSetVariableBoundsArray_PBPOLRAJAHIOP(OPFLOW opflow,double *xl_dev,double *xu_dev)
{
  //  PBPOLRAJAHIOP      pbpolrajahiop=(PBPOLRAJAHIOP)opflow->model;
  PbpolModelRajaHiop* pbpolrajahiop = reinterpret_cast<PbpolModelRajaHiop*>(opflow->model);
  BUSParamsRajaHiop      *busparams=&pbpolrajahiop->busparams;
  GENParamsRajaHiop      *genparams=&pbpolrajahiop->genparams;

  PetscFunctionBegin;
  //  PetscPrintf(MPI_COMM_SELF,"Entered variable bounds\n");

  int* xidx = busparams->xidx_dev_;
  int* ispvpq = busparams->ispvpq_dev_;
  int* isref = busparams->isref_dev_;
  int* isisolated = busparams->isisolated_dev_;
  double* va = busparams->va_dev_;
  double* vm = busparams->vm_dev_;
  double* vmin = busparams->vmin_dev_;
  double* vmax = busparams->vmax_dev_;

  /* Bounds for bus voltages */
  RAJA::forall<exago_raja_exec>(RAJA::RangeSegment(0, busparams->nbus)/* index set here */,
    RAJA_LAMBDA (RAJA::Index_type i)
    {
      xl_dev[xidx[i]] = ispvpq[i]*PETSC_NINFINITY + isisolated[i]*va[i] + isref[i]*va[i]*PETSC_PI/180.0;
      xu_dev[xidx[i]] = ispvpq[i]*PETSC_INFINITY  + isisolated[i]*va[i] + isref[i]*va[i]*PETSC_PI/180.0;
    
      xl_dev[xidx[i]+1] = isref[i]*vmin[i]  + ispvpq[i]*vmin[i] + isisolated[i]*vm[i];
      xu_dev[xidx[i]+1] = isref[i]*vmax[i]  + ispvpq[i]*vmax[i] + isisolated[i]*vm[i];
    }
  );

  int*   idx = genparams->xidx_dev_;
  double* pb = genparams->pb_dev_;
  double* pt = genparams->pt_dev_;
  double* qb = genparams->qb_dev_;
  double* qt = genparams->qt_dev_;

  /* Generator lower and upper bounds on variables */
  RAJA::forall< exago_raja_exec >(RAJA::RangeSegment(0, genparams->ngenON),
    RAJA_LAMBDA (RAJA::Index_type i)
    {
      xl_dev[idx[i]]   = pb[i];
      xu_dev[idx[i]]   = pt[i];
      xl_dev[idx[i]+1] = qb[i];
      xu_dev[idx[i]+1] = qt[i];
    }
  );
  //  PetscPrintf(MPI_COMM_SELF,"Exit variable bounds\n");

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeSparseJacobian_PBPOLRAJAHIOP(OPFLOW opflow,int *iJacS_dev, int *jJacS_dev,double *MJacS_dev)
{
  PbpolModelRajaHiop* pbpolrajahiop = reinterpret_cast<PbpolModelRajaHiop*>(opflow->model);
  GENParamsRajaHiop     *genparams=&pbpolrajahiop->genparams;
  //  PetscPrintf(MPI_COMM_SELF,"Entered sparse jacobian\n");

  if(iJacS_dev != NULL && jJacS_dev != NULL) {
    /* Generator contributions for row, col entries */
    int* g_gidx = genparams->gidx_dev_;
    int* g_xidx = genparams->xidx_dev_;
    int* jacsp_idx = genparams->jacsp_idx_dev_;
    int* jacsq_idx = genparams->jacsq_idx_dev_;
    RAJA::forall< exago_raja_exec >( RAJA::RangeSegment(0, genparams->ngenON),
      RAJA_LAMBDA (RAJA::Index_type i) 
      {
	iJacS_dev[jacsp_idx[i]] = g_gidx[i];
	jJacS_dev[jacsp_idx[i]] = g_xidx[i];
	
	iJacS_dev[jacsq_idx[i]] = g_gidx[i]+1;
	jJacS_dev[jacsq_idx[i]] = g_xidx[i]+1;
	
      }
    );
  }

  if(MJacS_dev != NULL) {
    /* Generator contributions */
    int* jacsp_idx = genparams->jacsp_idx_dev_;
    int* jacsq_idx = genparams->jacsq_idx_dev_;
    RAJA::forall< exago_raja_exec >( RAJA::RangeSegment(0, genparams->ngenON),
      RAJA_LAMBDA (RAJA::Index_type i) 
      {
        MJacS_dev[jacsp_idx[i]] = -1.0;
        MJacS_dev[jacsq_idx[i]] = -1.0;
      }
    );
  }    
  //  PetscPrintf(MPI_COMM_SELF,"Exit sparse jacobian\n");
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeSparseHessian_PBPOLRAJAHIOP(OPFLOW opflow,const double *x_dev,int *iHSS_dev, int *jHSS_dev,double *MHSS_dev)
{
  PbpolModelRajaHiop* pbpolrajahiop = reinterpret_cast<PbpolModelRajaHiop*>(opflow->model);
  PetscErrorCode ierr;
  GENParamsRajaHiop     *genparams=&pbpolrajahiop->genparams;
  PS             ps=opflow->ps;
  double         obj_factor = opflow->obj_factor;
  int            isobj_gencost=opflow->obj_gencost;
  double         MVAbase=ps->MVAbase;
  double         flps=0.0;

  //  PetscPrintf(MPI_COMM_SELF,"Entered sparse Hessian\n");

  if(iHSS_dev != NULL && jHSS_dev != NULL) {
    /* Generator contributions for row,col numbers */
    int* g_xidx = genparams->xidx_dev_;
    int* jacsp_idx = genparams->jacsp_idx_dev_;
    int* jacsq_idx = genparams->jacsq_idx_dev_;
    RAJA::forall< exago_raja_exec >( RAJA::RangeSegment(0, genparams->ngenON),
      RAJA_LAMBDA (RAJA::Index_type i) 
      {
        iHSS_dev[jacsp_idx[i]] = g_xidx[i];
        jHSS_dev[jacsp_idx[i]] = g_xidx[i];
      
	iHSS_dev[jacsq_idx[i]] = g_xidx[i]+1;
	jHSS_dev[jacsq_idx[i]] = g_xidx[i]+1;
      }
    );
    flps += 2.0*genparams->ngenON;
  }

  if(MHSS_dev != NULL) {
    /* Generator contributions */
    int* jacsp_idx = genparams->jacsp_idx_dev_;
    int* jacsq_idx = genparams->jacsq_idx_dev_;
    double* cost_alpha = genparams->cost_alpha_dev_;

    RAJA::forall< exago_raja_exec >( RAJA::RangeSegment(0, genparams->ngenON),
      RAJA_LAMBDA (RAJA::Index_type i) 
      {
        MHSS_dev[jacsp_idx[i]] = isobj_gencost*obj_factor*2.0*cost_alpha[i]*MVAbase*MVAbase;
	MHSS_dev[jacsq_idx[i]]  = 0.0;
      }
    );
    flps += 5.0*genparams->ngenON;
  }    
  //  PetscPrintf(MPI_COMM_SELF,"Exit sparse hessian\n");

  ierr = PetscLogFlops(flps);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/**
 * @param[inout] JacD_dev Jacobian matrix with size * 2*busparams->nbus by nxdense
 */
PetscErrorCode OPFLOWComputeDenseEqualityConstraintJacobian_PBPOLRAJAHIOP(OPFLOW opflow,const double *x_dev,double *JacD_dev)
{
  PbpolModelRajaHiop* pbpolrajahiop = reinterpret_cast<PbpolModelRajaHiop*>(opflow->model);
  BUSParamsRajaHiop      *busparams=&pbpolrajahiop->busparams;
  LINEParamsRajaHiop     *lineparams=&pbpolrajahiop->lineparams;
  int            nxsparse=2*opflow->ps->ngenON;
  int            nx=opflow->nx;
  const int      nxdense=nx-nxsparse;
  double         flps = 0.0;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  //  PetscPrintf(MPI_COMM_SELF,"Entered equality constrained dense jacobian\n");

  if(JacD_dev == NULL) PetscFunctionReturn(0);

  /** Although we will perform atomic operations on this view, we do not want
   * to use `RAJA::make_atomic_view` because explicit `RAJA::atomicAdd`s make
   * the atomic requirement clearer in the code */
  RAJA::View<double, RAJA::Layout<2>> JacD_view(JacD_dev, 2*busparams->nbus, nxdense);

  /* Zero out JacD */
  auto& resmgr = umpire::ResourceManager::getInstance();
  umpire::Allocator alloc;
#ifdef EXAGO_HAVE_GPU
  alloc = resmgr.getAllocator("DEVICE");
#else
  alloc = resmgr.getAllocator("HOST");
#endif
  registerWith(JacD_dev,2*busparams->nbus * nxdense,resmgr,alloc);
  resmgr.memset(JacD_dev,0);

  /* Jacobian from bus contributions */
  int* isisolated = busparams->isisolated_dev_;
  int* ispvpq = busparams->ispvpq_dev_;
  double* gl = busparams->gl_dev_;
  double* bl = busparams->bl_dev_;
  int* b_xidx = busparams->xidx_dev_;
  int* b_gidx = busparams->gidx_dev_;
  RAJA::forall<exago_raja_exec>(RAJA::RangeSegment(0, busparams->nbus),
    RAJA_LAMBDA (RAJA::Index_type i)
    {
      double Vm    = x_dev[b_xidx[i]+1];
      int     row[2],col[4];
      double  val[8];

      row[0] = b_gidx[i];
      row[1] = b_gidx[i]+1;
    
      col[0] = b_xidx[i]   - nxsparse;
      col[1] = b_xidx[i]+1 - nxsparse;

      val[0] = isisolated[i]*1.0 + ispvpq[i]*0.0;
      val[1] = isisolated[i]*0.0 + ispvpq[i]*2*Vm*gl[i];
      val[2] = 0.0;
      val[3] = isisolated[i]*1.0 + ispvpq[i]*-2*Vm*bl[i];

      RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[0], col[0]), val[0]);
      RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[0], col[1]), val[1]);
      RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[1], col[0]), val[2]);
      RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[1], col[1]), val[3]);
    }
  );
  flps += 14.0*busparams->nbus;


  /* Jacobian from line contributions */
  double* Gff = lineparams->Gff_dev_;
  double* Gtt = lineparams->Gtt_dev_;
  double* Gft = lineparams->Gft_dev_;
  double* Gtf = lineparams->Gtf_dev_;

  double* Bff = lineparams->Bff_dev_;
  double* Btt = lineparams->Btt_dev_;
  double* Bft = lineparams->Bft_dev_;
  double* Btf = lineparams->Btf_dev_;

  int* xidxf = lineparams->xidxf_dev_;
  int* xidxt = lineparams->xidxt_dev_;
  int* geqidxf = lineparams->geqidxf_dev_;
  int* geqidxt = lineparams->geqidxt_dev_;

  RAJA::forall<exago_raja_exec>(RAJA::RangeSegment(0, lineparams->nlineON),
    RAJA_LAMBDA (RAJA::Index_type i)
    {
      int     row[2],col[4];
      double  val[8];
      double thetaf=x_dev[xidxf[i]], Vmf=x_dev[xidxf[i]+1];
      double thetat=x_dev[xidxt[i]], Vmt=x_dev[xidxt[i]+1];
      double thetaft=thetaf-thetat;
      double thetatf=thetat-thetaf;

      row[0] = geqidxf[i];
      row[1] = geqidxf[i]+1;

      col[0] = xidxf[i] - nxsparse; 
      col[1] = xidxf[i]+1 - nxsparse; 
      col[2] = xidxt[i] - nxsparse;
      col[3] = xidxt[i]+1 - nxsparse;
      
      /* dPf_dthetaf */
      val[0] = Vmf*Vmt*(-Gft[i]*sin(thetaft) + Bft[i]*cos(thetaft));
      /*dPf_dVmf */
      val[1] = 2*Gff[i]*Vmf + Vmt*(Gft[i]*cos(thetaft) + Bft[i]*sin(thetaft));
      /*dPf_dthetat */
      val[2] = Vmf*Vmt*(Gft[i]*sin(thetaft) - Bft[i]*cos(thetaft));
      /* dPf_dVmt */
      val[3] = Vmf*(Gft[i]*cos(thetaft) + Bft[i]*sin(thetaft));

      RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[0], col[0]), val[0]);
      RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[0], col[1]), val[1]);
      RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[0], col[2]), val[2]);
      RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[0], col[3]), val[3]);

      /* dQf_dthetaf */
      val[4] = Vmf*Vmt*(Bft[i]*sin(thetaft) + Gft[i]*cos(thetaft));
      /* dQf_dVmf */
      val[5] = -2*Bff[i]*Vmf + Vmt*(-Bft[i]*cos(thetaft) + Gft[i]*sin(thetaft));
      /* dQf_dthetat */
      val[6] = Vmf*Vmt*(-Bft[i]*sin(thetaft) - Gft[i]*cos(thetaft));
      /* dQf_dVmt */
      val[7] = Vmf*(-Bft[i]*cos(thetaft) + Gft[i]*sin(thetaft));

      RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[1], col[0]), val[4]);
      RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[1], col[1]), val[5]);
      RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[1], col[2]), val[6]);
      RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[1], col[3]), val[7]);

      row[0] = geqidxt[i];
      row[1] = geqidxt[i]+1;
      
      col[0] = xidxt[i] - nxsparse; 
      col[1] = xidxt[i]+1 - nxsparse;
      col[2] = xidxf[i] - nxsparse;
      col[3] = xidxf[i]+1 - nxsparse;
      
      /* dPt_dthetat */
      val[0] = Vmt*Vmf*(-Gtf[i]*sin(thetatf) + Btf[i]*cos(thetatf));
      /* dPt_dVmt */
      val[1] = 2*Gtt[i]*Vmt + Vmf*(Gtf[i]*cos(thetatf) + Btf[i]*sin(thetatf));
      /* dPt_dthetaf */
      val[2] = Vmt*Vmf*(Gtf[i]*sin(thetatf) - Btf[i]*cos(thetatf));
      /* dPt_dVmf */
      val[3] = Vmt*(Gtf[i]*cos(thetatf) + Btf[i]*sin(thetatf));
      
      RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[0], col[0]), val[0]);
      RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[0], col[1]), val[1]);
      RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[0], col[2]), val[2]);
      RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[0], col[3]), val[3]);
    
      /* dQt_dthetat */
      val[4] = Vmt*Vmf*(Btf[i]*sin(thetatf) + Gtf[i]*cos(thetatf));
      /* dQt_dVmt */
      val[5] = -2*Btt[i]*Vmt + Vmf*(-Btf[i]*cos(thetatf) + Gtf[i]*sin(thetatf));
      /* dQt_dthetaf */
      val[6] = Vmt*Vmf*(-Btf[i]*sin(thetatf) - Gtf[i]*cos(thetatf));
      /* dQt_dVmf */
      val[7] = Vmt*(-Btf[i]*cos(thetatf) + Gtf[i]*sin(thetatf));
      
      RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[1], col[0]), val[4]);
      RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[1], col[1]), val[5]);
      RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[1], col[2]), val[6]);
      RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[1], col[3]), val[7]);
    }
  );
  flps += (188+(16*EXAGO_FLOPS_COSOP)+(16*EXAGO_FLOPS_SINOP))*lineparams->nlineON;

  ierr = PetscLogFlops(flps);CHKERRQ(ierr);

  //  PetscPrintf(MPI_COMM_SELF,"Exit equality dense jacobian\n");

  PetscFunctionReturn(0);
}
    
/**
 * @param[inout] JacD_dev Jacobian matrix with size * 2*lineparams->nlinelim by nxdense
 */
PetscErrorCode OPFLOWComputeDenseInequalityConstraintJacobian_PBPOLRAJAHIOP(OPFLOW opflow,const double *x_dev,double *JacD_dev)
{
  PbpolModelRajaHiop* pbpolrajahiop = reinterpret_cast<PbpolModelRajaHiop*>(opflow->model);
  LINEParamsRajaHiop     *lineparams=&pbpolrajahiop->lineparams;
  PetscErrorCode      ierr;
  double              flps=0.0;
  int       nxsparse=2*opflow->ps->ngenON;
  int       nx=opflow->nx;

  PetscFunctionBegin;
  //  PetscPrintf(MPI_COMM_SELF,"Enter inequality dense jacobian\n");
  /* Return if there are no inequality constraints */
  if(!lineparams->nlinelim) {
    //    PetscPrintf(MPI_COMM_SELF,"No inequality constraints. Exit inequality dense jacobian\n");
    PetscFunctionReturn(0);
  }

  if(JacD_dev == NULL) PetscFunctionReturn(0);

  /** Although we will perform atomic operations on this view, we do not want
   * to use `RAJA::make_atomic_view` because explicit `RAJA::atomicAdd`s make
   * the atomic requirement clearer in the code */
  RAJA::View<double, RAJA::Layout<2>> JacD_view(JacD_dev, 2*lineparams->nlinelim, nx-nxsparse);
  
  /* Zero out JacD */
  auto& resmgr = umpire::ResourceManager::getInstance();
  umpire::Allocator alloc;
#ifdef EXAGO_HAVE_GPU
  alloc = resmgr.getAllocator("DEVICE");
#else
  alloc = resmgr.getAllocator("HOST");
#endif
  registerWith(JacD_dev, 2 * lineparams->nlinelim * (nx - nxsparse), resmgr, alloc);
  resmgr.memset(JacD_dev,0);

  /* Line contributions */
  double* Gff_arr = lineparams->Gff_dev_;
  double* Gtt_arr = lineparams->Gtt_dev_;
  double* Gft_arr = lineparams->Gft_dev_;
  double* Gtf_arr = lineparams->Gtf_dev_;

  double* Bff_arr = lineparams->Bff_dev_;
  double* Btt_arr = lineparams->Btt_dev_;
  double* Bft_arr = lineparams->Bft_dev_;
  double* Btf_arr = lineparams->Btf_dev_;

  int* linelimidx = lineparams->linelimidx_dev_;
  int* xidxf = lineparams->xidxf_dev_;
  int* xidxt = lineparams->xidxt_dev_;
  int* gineqidx = lineparams->gineqidx_dev_;
  RAJA::forall<exago_raja_exec>(RAJA::RangeSegment(0, lineparams->nlinelim),
    RAJA_LAMBDA (RAJA::Index_type i)
    {
      int    j=linelimidx[i];
      int       row[2],col[4];
      double    val[4];
      double Pf,Qf,Pt,Qt;
      double thetaf=x_dev[xidxf[j]], Vmf=x_dev[xidxf[j]+1];
      double thetat=x_dev[xidxt[j]], Vmt=x_dev[xidxt[j]+1];
      double thetaft=thetaf-thetat;
      double thetatf=thetat-thetaf;
      double dSf2_dPf, dSf2_dQf, dSt2_dPt, dSt2_dQt;
      double dPf_dthetaf,dPf_dVmf,dPf_dthetat,dPf_dVmt;
      double dQf_dthetaf,dQf_dVmf,dQf_dthetat,dQf_dVmt;
      double dPt_dthetaf,dPt_dVmf,dPt_dthetat,dPt_dVmt;
      double dQt_dthetaf,dQt_dVmf,dQt_dthetat,dQt_dVmt;
      double dSf2_dthetaf,dSf2_dVmf,dSf2_dthetat,dSf2_dVmt;
      double dSt2_dthetaf,dSt2_dVmf,dSt2_dthetat,dSt2_dVmt;
      double Gff = Gff_arr[j], Bff = Bff_arr[j];
      double Gft = Gft_arr[j], Bft = Bft_arr[j];
      double Gtf = Gtf_arr[j], Btf = Btf_arr[j];
      double Gtt = Gtt_arr[j], Btt = Btt_arr[j];
      
      Pf = Gff*Vmf*Vmf  + Vmf*Vmt*(Gft*cos(thetaft)  + Bft*sin(thetaft));
      Qf = -Bff*Vmf*Vmf + Vmf*Vmt*(-Bft*cos(thetaft) + Gft*sin(thetaft));
      Pt = Gtt*Vmt*Vmt  + Vmt*Vmf*(Gtf*cos(thetatf)  + Btf*sin(thetatf));
      Qt = -Btt*Vmt*Vmt + Vmt*Vmf*(-Btf*cos(thetatf) + Gtf*sin(thetatf));
      
      dSf2_dPf = 2*Pf;
      dSf2_dQf = 2*Qf;
      dSt2_dPt = 2*Pt;
      dSt2_dQt = 2*Qt;
      
      dPf_dthetaf = Vmf*Vmt*(-Gft*sin(thetaft) + Bft*cos(thetaft));
      dPf_dVmf    = 2*Gff*Vmf + Vmt*(Gft*cos(thetaft) + Bft*sin(thetaft));
      dPf_dthetat = Vmf*Vmt*(Gft*sin(thetaft) - Bft*cos(thetaft));
      dPf_dVmt    = Vmf*(Gft*cos(thetaft) + Bft*sin(thetaft));
      
      dQf_dthetaf = Vmf*Vmt*(Bft*sin(thetaft) + Gft*cos(thetaft));
      dQf_dVmf    = -2*Bff*Vmf + Vmt*(-Bft*cos(thetaft) + Gft*sin(thetaft));
      dQf_dthetat = Vmf*Vmt*(-Bft*sin(thetaft) - Gft*cos(thetaft));
      dQf_dVmt    = Vmf*(-Bft*cos(thetaft) + Gft*sin(thetaft));
      
      dPt_dthetat = Vmt*Vmf*(-Gtf*sin(thetatf) + Btf*cos(thetatf));
      dPt_dVmt    = 2*Gtt*Vmt + Vmf*(Gtf*cos(thetatf) + Btf*sin(thetatf));
      dPt_dthetaf = Vmt*Vmf*(Gtf*sin(thetatf) - Btf*cos(thetatf));
      dPt_dVmf    = Vmt*(Gtf*cos(thetatf) + Btf*sin(thetatf));
      
      dQt_dthetat = Vmt*Vmf*(Btf*sin(thetatf) + Gtf*cos(thetatf));
      dQt_dVmt    = -2*Btt*Vmt + Vmf*(-Btf*cos(thetatf) + Gtf*sin(thetatf));
      dQt_dthetaf = Vmt*Vmf*(-Btf*sin(thetatf) - Gtf*cos(thetatf));
      dQt_dVmf    = Vmt*(-Btf*cos(thetatf) + Gtf*sin(thetatf));
      
      dSf2_dthetaf = dSf2_dPf*dPf_dthetaf + dSf2_dQf*dQf_dthetaf;
      dSf2_dthetat = dSf2_dPf*dPf_dthetat + dSf2_dQf*dQf_dthetat;
      dSf2_dVmf    = dSf2_dPf*dPf_dVmf    + dSf2_dQf*dQf_dVmf;
      dSf2_dVmt    = dSf2_dPf*dPf_dVmt    + dSf2_dQf*dQf_dVmt;
      
      row[0] = gineqidx[i];

      col[0] = xidxf[j]   - nxsparse;
      col[1] = xidxf[j]+1 - nxsparse;
      col[2] = xidxt[j]   - nxsparse;
      col[3] = xidxt[j]+1 - nxsparse;
      
      val[0] = dSf2_dthetaf;
      val[1] = dSf2_dVmf;
      val[2] = dSf2_dthetat;
      val[3] = dSf2_dVmt;
      
      RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[0], col[0]), val[0]);
      RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[0], col[1]), val[1]);
      RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[0], col[2]), val[2]);
      RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[0], col[3]), val[3]);
      
      dSt2_dthetaf = dSt2_dPt*dPt_dthetaf + dSt2_dQt*dQt_dthetaf;
      dSt2_dthetat = dSt2_dPt*dPt_dthetat + dSt2_dQt*dQt_dthetat;
      dSt2_dVmf    = dSt2_dPt*dPt_dVmf    + dSt2_dQt*dQt_dVmf;
      dSt2_dVmt    = dSt2_dPt*dPt_dVmt    + dSt2_dQt*dQt_dVmt;
      
      row[0] = gineqidx[i]+1;
      
      col[0] = xidxt[j]   - nxsparse;
      col[1] = xidxt[j]+1 - nxsparse;
      col[2] = xidxf[j]   - nxsparse;
      col[3] = xidxf[j]+1 - nxsparse;
      
      val[0] = dSt2_dthetat;
      val[1] = dSt2_dVmt;
      val[2] = dSt2_dthetaf;
      val[3] = dSt2_dVmf;

      RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[0], col[0]), val[0]);
      RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[0], col[1]), val[1]);
      RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[0], col[2]), val[2]);
      RAJA::atomicAdd<exago_raja_atomic>(&JacD_view(row[0], col[3]), val[3]);
    }
  );
  //  PetscPrintf(MPI_COMM_SELF,"Exit inequality dense jacobian\n");
  flps += (183+(20*EXAGO_FLOPS_COSOP)+(20*EXAGO_FLOPS_SINOP))*lineparams->nlinelim;
  ierr = PetscLogFlops(flps);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/**
 * @param[inout] HDD_dev Hessian matrix with size nxdense x nxdense
 */
PetscErrorCode OPFLOWComputeDenseEqualityConstraintHessian_PBPOLRAJAHIOP(OPFLOW opflow,const double *x_dev,const double* lambda_dev, double *HDD_dev)
{
  PbpolModelRajaHiop* pbpolrajahiop = reinterpret_cast<PbpolModelRajaHiop*>(opflow->model);
  BUSParamsRajaHiop   *busparams=&pbpolrajahiop->busparams;
  LINEParamsRajaHiop  *lineparams=&pbpolrajahiop->lineparams;
  int                 nxsparse = 2*opflow->ps->ngenON;
  const int           nxdense = 2*opflow->ps->nbus;
  PetscErrorCode      ierr;
  double              flps=0.0;

  PetscFunctionBegin;
  //  PetscPrintf(MPI_COMM_SELF,"Enter equality dense hessian\n");

  /* Hessian from bus contributions */
  int* b_xidx = busparams->xidx_dev_;
  int* b_gidx = busparams->gidx_dev_;
  int* ispvpq = busparams->ispvpq_dev_;
  double* gl = busparams->gl_dev_;
  double* bl = busparams->bl_dev_;

  /** Although we will perform atomic operations on this view, we do not want
   * to use `RAJA::make_atomic_view` because explicit `RAJA::atomicAdd`s make
   * the atomic requirement clearer in the code */
  RAJA::View<double, RAJA::Layout<2>> HDD_view(HDD_dev, nxdense, nxdense);

  RAJA::forall<exago_raja_exec>(RAJA::RangeSegment(0, busparams->nbus),
    RAJA_LAMBDA (RAJA::Index_type i)
    {
      int row,col;
      double val;
      row = b_xidx[i] + 1 - nxsparse;
      col = row;
      val = ispvpq[i]*(lambda_dev[b_gidx[i]]*2*gl[i] + lambda_dev[b_gidx[i]+1]*(-2*bl[i]));
      RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row, col),val);
    }
  );
  flps += 7.0*busparams->nbus;

  /* Hessian from line contributions */
  double* Gff_arr = lineparams->Gff_dev_;
  double* Gtt_arr = lineparams->Gtt_dev_;
  double* Gft_arr = lineparams->Gft_dev_;
  double* Gtf_arr = lineparams->Gtf_dev_;

  double* Bff_arr = lineparams->Bff_dev_;
  double* Btt_arr = lineparams->Btt_dev_;
  double* Bft_arr = lineparams->Bft_dev_;
  double* Btf_arr = lineparams->Btf_dev_;

  int* xidxf = lineparams->xidxf_dev_;
  int* xidxt = lineparams->xidxt_dev_;
  int* geqidxf = lineparams->geqidxf_dev_;
  int* geqidxt = lineparams->geqidxt_dev_;

  RAJA::forall<exago_raja_exec>(RAJA::RangeSegment(0, lineparams->nlineON),
    RAJA_LAMBDA (RAJA::Index_type i)
    {
      int    gloc;
      int    row[2],col[4];
      double val[8];
      double Gff,Bff,Gft,Bft,Gtf,Btf,Gtt,Btt;
      Gff = Gff_arr[i];
      Bff = Bff_arr[i];
      Gft = Gft_arr[i];
      Bft = Bft_arr[i];
      Gtf = Gtf_arr[i];
      Btf = Btf_arr[i];
      Gtt = Gtt_arr[i];
      Btt = Btt_arr[i];
    
    double thetaf=x_dev[xidxf[i]], Vmf=x_dev[xidxf[i]+1];
    double thetat=x_dev[xidxt[i]], Vmt=x_dev[xidxt[i]+1];
    double thetaft=thetaf-thetat;
    double thetatf=thetat-thetaf;
    
    double dPf_dthetaf_dthetaf,dPf_dthetaf_dVmf,dPf_dthetaf_dthetat,dPf_dthetaf_dVmt;
    double dPf_dVmf_dthetaf,   dPf_dVmf_dVmf,   dPf_dVmf_dthetat,   dPf_dVmf_dVmt;
    double dPf_dthetat_dthetaf,dPf_dthetat_dVmf,dPf_dthetat_dthetat,dPf_dthetat_dVmt;
    double dPf_dVmt_dthetaf,   dPf_dVmt_dVmf,   dPf_dVmt_dthetat,   dPf_dVmt_dVmt;
    
    /* dPf_dthetaf = Vmf*Vmt*(-Gft*sin(thetaft) + Bft*cos(thetaft)); */
    dPf_dthetaf_dthetaf = -Vmf*Vmt*(Gft*cos(thetaft) + Bft*sin(thetaft));
    dPf_dthetaf_dVmf    =     Vmt*(-Gft*sin(thetaft) + Bft*cos(thetaft));
    dPf_dthetaf_dthetat =  Vmf*Vmt*(Gft*cos(thetaft) + Bft*sin(thetaft));
    dPf_dthetaf_dVmt    =     Vmf*(-Gft*sin(thetaft) + Bft*cos(thetaft));
    
    /* dPf_Vmf  = 2*Gff*Vmf + Vmt*(Gft*cos(thetaft) + Bft*sin(thetaft)); */
    dPf_dVmf_dthetaf    =  Vmt*(-Gft*sin(thetaft) + Bft*cos(thetaft)); 
    dPf_dVmf_dVmf       =  2*Gff;
    dPf_dVmf_dthetat    =  Vmt*(Gft*sin(thetaft) - Bft*cos(thetaft));
    dPf_dVmf_dVmt       =      (Gft*cos(thetaft) + Bft*sin(thetaft));
    
    /* dPf_dthetat = Vmf*Vmt*(Gft*sin(thetaft) - Bft*cos(thetaft)); */
    dPf_dthetat_dthetaf = Vmf*Vmt*(Gft*cos(thetaft)  + Bft*sin(thetaft));
    dPf_dthetat_dVmf    =     Vmt*(Gft*sin(thetaft)  - Bft*cos(thetaft));
    dPf_dthetat_dthetat = Vmf*Vmt*(-Gft*cos(thetaft) - Bft*sin(thetaft));
    dPf_dthetat_dVmt    =     Vmf*(Gft*sin(thetaft)  - Bft*cos(thetaft));
    
    /* dPf_dVmt = Vmf*(Gft*cos(thetaft) + Bft*sin(thetaft)); */
    dPf_dVmt_dthetaf    = Vmf*(-Gft*sin(thetaft) + Bft*cos(thetaft)); 
    dPf_dVmt_dVmf       =      (Gft*cos(thetaft) + Bft*sin(thetaft));
    dPf_dVmt_dthetat    = Vmf*(Gft*sin(thetaft) - Bft*cos(thetaft));
    dPf_dVmt_dVmt       = 0.0;
    
    double dQf_dthetaf_dthetaf,dQf_dthetaf_dVmf,dQf_dthetaf_dthetat,dQf_dthetaf_dVmt;
    double dQf_dVmf_dthetaf,   dQf_dVmf_dVmf,   dQf_dVmf_dthetat,   dQf_dVmf_dVmt;
    double dQf_dthetat_dthetaf,dQf_dthetat_dVmf,dQf_dthetat_dthetat,dQf_dthetat_dVmt;
    double dQf_dVmt_dthetaf,   dQf_dVmt_dVmf,   dQf_dVmt_dthetat,   dQf_dVmt_dVmt;
    
    /* dQf_dthetaf = Vmf*Vmt*(Bft*sin(thetaft) + Gft*cos(thetaft)); */
    dQf_dthetaf_dthetaf = Vmf*Vmt*(Bft*cos(thetaft)  - Gft*sin(thetaft));
    dQf_dthetaf_dVmf    =     Vmt*(Bft*sin(thetaft)  + Gft*cos(thetaft));
    dQf_dthetaf_dthetat = Vmf*Vmt*(-Bft*cos(thetaft) + Gft*sin(thetaft));
    dQf_dthetaf_dVmt    =     Vmf*(Bft*sin(thetaft) + Gft*cos(thetaft));
    
    /* dQf_dVmf = -2*Bff*Vmf + Vmt*(-Bft*cos(thetaft) + Gft*sin(thetaft)); */
    dQf_dVmf_dthetaf    =  Vmt*(Bft*sin(thetaft) + Gft*cos(thetaft)); 
    dQf_dVmf_dVmf       = -2*Bff;
    dQf_dVmf_dthetat    =  Vmt*(-Bft*sin(thetaft) - Gft*cos(thetaft));
    dQf_dVmf_dVmt       =      (-Bft*cos(thetaft) + Gft*sin(thetaft));
    
    /* dQf_dthetat = Vmf*Vmt*(-Bft*sin(thetaft) - Gft*cos(thetaft)); */
    dQf_dthetat_dthetaf = Vmf*Vmt*(-Bft*cos(thetaft) + Gft*sin(thetaft));
    dQf_dthetat_dVmf    =     Vmt*(-Bft*sin(thetaft) - Gft*cos(thetaft));
    dQf_dthetat_dthetat = Vmf*Vmt*( Bft*cos(thetaft) - Gft*sin(thetaft));
    dQf_dthetat_dVmt    =     Vmf*(-Bft*sin(thetaft) - Gft*cos(thetaft));
    
    /* dQf_dVmt = Vmf*(-Bft*cos(thetaft) + Gft*sin(thetaft)); */
    dQf_dVmt_dthetaf    = Vmf*(Bft*sin(thetaft) + Gft*cos(thetaft)); 
    dQf_dVmt_dVmf       =    (-Bft*cos(thetaft) + Gft*sin(thetaft));
    dQf_dVmt_dthetat    = Vmf*(-Bft*sin(thetaft) - Gft*cos(thetaft));
    dQf_dVmt_dVmt       = 0.0;
    
    row[0] = xidxf[i]     - nxsparse; 
    row[1] = xidxf[i] + 1 - nxsparse;
    col[0] = xidxf[i]     - nxsparse; 
    col[1] = xidxf[i] + 1 - nxsparse; 
    col[2] = xidxt[i]     - nxsparse;
    col[3] = xidxt[i] + 1 - nxsparse;
    
    gloc=geqidxf[i];
    
    val[0] = val[1] = val[2] = val[3] = val[4] = val[5] = val[6] = val[7] = 0.0;
    
    val[0] = lambda_dev[gloc]*dPf_dthetaf_dthetaf + lambda_dev[gloc+1]*dQf_dthetaf_dthetaf;
    val[1] = lambda_dev[gloc]*dPf_dthetaf_dVmf    + lambda_dev[gloc+1]*dQf_dthetaf_dVmf;
    val[2] = lambda_dev[gloc]*dPf_dthetaf_dthetat + lambda_dev[gloc+1]*dQf_dthetaf_dthetat;
    val[3] = lambda_dev[gloc]*dPf_dthetaf_dVmt    + lambda_dev[gloc+1]*dQf_dthetaf_dVmt;

    val[4] = lambda_dev[gloc]*dPf_dVmf_dthetaf + lambda_dev[gloc+1]*dQf_dVmf_dthetaf;
    val[5] = lambda_dev[gloc]*dPf_dVmf_dVmf    + lambda_dev[gloc+1]*dQf_dVmf_dVmf;
    val[6] = lambda_dev[gloc]*dPf_dVmf_dthetat + lambda_dev[gloc+1]*dQf_dVmf_dthetat;
    val[7] = lambda_dev[gloc]*dPf_dVmf_dVmt    + lambda_dev[gloc+1]*dQf_dVmf_dVmt;

    RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[0]), val[0]);
    RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[1]), val[1]);
    RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[2]), val[2]);
    RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[3]), val[3]);
                                                              
    RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[1], col[0]), val[4]);
    RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[1], col[1]), val[5]);
    RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[1], col[2]), val[6]);
    RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[1], col[3]), val[7]);

    row[0] = xidxt[i]     - nxsparse;
    row[1] = xidxt[i] + 1 - nxsparse;

    col[0] = xidxf[i]     - nxsparse;
    col[1] = xidxf[i] + 1 - nxsparse;
    col[2] = xidxt[i]     - nxsparse;
    col[3] = xidxt[i] + 1 - nxsparse;
    
    val[0] = val[1] = val[2] = val[3] = val[4] = val[5] = val[6] = val[7] = 0.0;
    
    val[0] = lambda_dev[gloc]*dPf_dthetat_dthetaf + lambda_dev[gloc+1]*dQf_dthetat_dthetaf;
    val[1] = lambda_dev[gloc]*dPf_dthetat_dVmf    + lambda_dev[gloc+1]*dQf_dthetat_dVmf;
    val[2] = lambda_dev[gloc]*dPf_dthetat_dthetat + lambda_dev[gloc+1]*dQf_dthetat_dthetat;
    val[3] = lambda_dev[gloc]*dPf_dthetat_dVmt    + lambda_dev[gloc+1]*dQf_dthetat_dVmt;

    val[4] = lambda_dev[gloc]*dPf_dVmt_dthetaf + lambda_dev[gloc+1]*dQf_dVmt_dthetaf;
    val[5] = lambda_dev[gloc]*dPf_dVmt_dVmf    + lambda_dev[gloc+1]*dQf_dVmt_dVmf;
    val[6] = lambda_dev[gloc]*dPf_dVmt_dthetat + lambda_dev[gloc+1]*dQf_dVmt_dthetat;
    val[7] = lambda_dev[gloc]*dPf_dVmt_dVmt    + lambda_dev[gloc+1]*dQf_dVmt_dVmt;

    RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[0]), val[0]);
    RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[1]), val[1]);
    RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[2]), val[2]);
    RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[3]), val[3]);
                                                              
    RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[1], col[0]), val[4]);
    RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[1], col[1]), val[5]);
    RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[1], col[2]), val[6]);
    RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[1], col[3]), val[7]);

    //    ierr = MatSetValues(H,2,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);
	
    double dPt_dthetat_dthetat,dPt_dthetat_dVmt,dPt_dthetat_dthetaf,dPt_dthetat_dVmf;
    double dPt_dVmt_dthetat,   dPt_dVmt_dVmt,   dPt_dVmt_dthetaf,   dPt_dVmt_dVmf;
    double dPt_dthetaf_dthetat,dPt_dthetaf_dVmt,dPt_dthetaf_dthetaf,dPt_dthetaf_dVmf;
    double dPt_dVmf_dthetat,   dPt_dVmf_dVmt,   dPt_dVmf_dthetaf,   dPt_dVmf_dVmf;

    /* dPt_dthetat = Vmf*Vmt*(-Gtf*sin(thetatf) + Btf*cos(thetatf)); */
    dPt_dthetat_dthetat = Vmf*Vmt*(-Gtf*cos(thetatf) - Btf*sin(thetatf));
    dPt_dthetat_dVmt    =     Vmf*(-Gtf*sin(thetatf) + Btf*cos(thetatf));
    dPt_dthetat_dthetaf =  Vmf*Vmt*(Gtf*cos(thetatf) + Btf*sin(thetatf));
    dPt_dthetat_dVmf    =     Vmt*(-Gtf*sin(thetatf) + Btf*cos(thetatf));
    
    /* dPt_Vmt  = 2*Gtt*Vmt + Vmf*(Gtf*cos(thetatf) + Btf*sin(thetatf)); */
    dPt_dVmt_dthetat    =  Vmf*(-Gtf*sin(thetatf) + Bft*cos(thetatf)); 
    dPt_dVmt_dVmt       =  2*Gtt;
    dPt_dVmt_dthetaf    =  Vmf*(Gtf*sin(thetatf) - Btf*cos(thetatf));
    dPt_dVmt_dVmf       =      (Gtf*cos(thetatf) + Btf*sin(thetatf));
    
    /* dPt_dthetaf = Vmf*Vmt*(Gtf*sin(thetatf) - Btf*cos(thetatf)); */
    dPt_dthetaf_dthetat = Vmf*Vmt*(Gtf*cos(thetatf)  + Btf*sin(thetatf));
    dPt_dthetaf_dVmt    =     Vmf*(Gtf*sin(thetatf)  - Btf*cos(thetatf));
    dPt_dthetaf_dthetaf = Vmf*Vmt*(-Gtf*cos(thetatf) - Btf*sin(thetatf));
    dPt_dthetaf_dVmf    =     Vmt*(Gtf*sin(thetatf)  - Btf*cos(thetatf));
    
    /* dPt_dVmf = Vmt*(Gtf*cos(thetatf) + Btf*sin(thetatf)); */
    dPt_dVmf_dthetat    = Vmt*(-Gtf*sin(thetatf) + Btf*cos(thetatf)); 
    dPt_dVmf_dVmt       =     (Gtf*cos(thetatf) + Btf*sin(thetatf));
    dPt_dVmf_dthetaf    = Vmt*(Gtf*sin(thetatf) - Btf*cos(thetatf));
    dPt_dVmf_dVmf       = 0.0;
    
    double dQt_dthetaf_dthetaf,dQt_dthetaf_dVmf,dQt_dthetaf_dthetat,dQt_dthetaf_dVmt;
    double dQt_dVmf_dthetaf,   dQt_dVmf_dVmf,   dQt_dVmf_dthetat,   dQt_dVmf_dVmt;
    double dQt_dthetat_dthetaf,dQt_dthetat_dVmf,dQt_dthetat_dthetat,dQt_dthetat_dVmt;
    double dQt_dVmt_dthetaf,   dQt_dVmt_dVmf,   dQt_dVmt_dthetat,   dQt_dVmt_dVmt;
    
    /* dQt_dthetat = Vmf*Vmt*(Btf*sin(thetatf) + Gtf*cos(thetatf)); */
    dQt_dthetat_dthetat = Vmf*Vmt*(Btf*cos(thetatf)  - Gtf*sin(thetatf));
    dQt_dthetat_dVmt    =     Vmf*(Btf*sin(thetatf)  + Gtf*cos(thetatf));
    dQt_dthetat_dthetaf = Vmf*Vmt*(-Btf*cos(thetatf) + Gtf*sin(thetatf));
    dQt_dthetat_dVmf    =     Vmf*(Btf*sin(thetatf) + Gtf*cos(thetatf));
    
    /* dQt_dVmt = -2*Btt*Vmt + Vmf*(-Btf*cos(thetatf) + Gtf*sin(thetatf)); */
    dQt_dVmt_dthetat    =  Vmf*(Btf*sin(thetatf) + Gtf*cos(thetatf)); 
    dQt_dVmt_dVmt       = -2*Btt;
    dQt_dVmt_dthetaf    =  Vmf*(-Btf*sin(thetatf) + Gtf*cos(thetatf));
    dQt_dVmt_dVmf       =      (-Btf*cos(thetatf) + Gtf*sin(thetatf));
    
    /* dQt_dthetaf = Vmf*Vmt*(-Btf*sin(thetatf) - Gtf*cos(thetatf)); */
    dQt_dthetaf_dthetat = Vmf*Vmt*(-Btf*cos(thetatf) + Gtf*sin(thetatf));
    dQt_dthetaf_dVmt    =     Vmf*(-Btf*sin(thetatf) - Gtf*cos(thetatf));
    dQt_dthetaf_dthetaf = Vmf*Vmt*( Btf*cos(thetatf) - Gtf*sin(thetatf));
    dQt_dthetaf_dVmf    =     Vmt*(-Btf*sin(thetatf) - Gtf*cos(thetatf));
    
    /* dQt_dVmf = Vmt*(-Btf*cos(thetatf) + Gtf*sin(thetatf)); */
    dQt_dVmf_dthetat    = Vmt*(Btf*sin(thetatf) + Gtf*cos(thetatf)); 
    dQt_dVmf_dVmt       =    (-Btf*cos(thetatf) + Gtf*sin(thetatf));
    dQt_dVmf_dthetaf    = Vmt*(-Btf*sin(thetatf) - Gtf*cos(thetatf));
    dQt_dVmf_dVmf       = 0.0;
    
    row[0] = xidxt[i]     - nxsparse;
    row[1] = xidxt[i] + 1 - nxsparse;
    col[0] = xidxt[i]     - nxsparse; 
    col[1] = xidxt[i] + 1 - nxsparse; 
    col[2] = xidxf[i]     - nxsparse; 
    col[3] = xidxf[i] + 1 - nxsparse;

    gloc   = geqidxt[i];

    val[0] = val[1] = val[2] = val[3] = val[4] = val[5] = val[6] = val[7] = 0.0;

    val[0] = lambda_dev[gloc]*dPt_dthetat_dthetat + lambda_dev[gloc+1]*dQt_dthetat_dthetat;
    val[1] = lambda_dev[gloc]*dPt_dthetat_dVmt    + lambda_dev[gloc+1]*dQt_dthetat_dVmt;
    val[2] = lambda_dev[gloc]*dPt_dthetat_dthetaf + lambda_dev[gloc+1]*dQt_dthetat_dthetaf;
    val[3] = lambda_dev[gloc]*dPt_dthetat_dVmf    + lambda_dev[gloc+1]*dQt_dthetat_dVmf;

    val[4] = lambda_dev[gloc]*dPt_dVmt_dthetat + lambda_dev[gloc+1]*dQt_dVmt_dthetat;
    val[5] = lambda_dev[gloc]*dPt_dVmt_dVmt    + lambda_dev[gloc+1]*dQt_dVmt_dVmt;
    val[6] = lambda_dev[gloc]*dPt_dVmt_dthetaf + lambda_dev[gloc+1]*dQt_dVmt_dthetaf;
    val[7] = lambda_dev[gloc]*dPt_dVmt_dVmf    + lambda_dev[gloc+1]*dQt_dVmt_dVmf;

    RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[0]), val[0]);
    RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[1]), val[1]);
    RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[2]), val[2]);
    RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[3]), val[3]);
                                                              
    RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[1], col[0]), val[4]);
    RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[1], col[1]), val[5]);
    RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[1], col[2]), val[6]);
    RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[1], col[3]), val[7]);

    row[0] = xidxf[i]     - nxsparse;
    row[1] = xidxf[i] + 1 - nxsparse;
    col[0] = xidxt[i]     - nxsparse; 
    col[1] = xidxt[i] + 1 - nxsparse; 
    col[2] = xidxf[i]     - nxsparse; 
    col[3] = xidxf[i] + 1 - nxsparse;
    
    val[0] = val[1] = val[2] = val[3] = val[4] = val[5] = val[6] = val[7] = 0.0;
    
    val[0] = lambda_dev[gloc]*dPt_dthetaf_dthetat + lambda_dev[gloc+1]*dQt_dthetaf_dthetat;
    val[1] = lambda_dev[gloc]*dPt_dthetaf_dVmt    + lambda_dev[gloc+1]*dQt_dthetaf_dVmt;
    val[2] = lambda_dev[gloc]*dPt_dthetaf_dthetaf + lambda_dev[gloc+1]*dQt_dthetaf_dthetaf;
    val[3] = lambda_dev[gloc]*dPt_dthetaf_dVmf    + lambda_dev[gloc+1]*dQt_dthetaf_dVmf;
    
    RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[0]), val[0]);
    RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[1]), val[1]);
    RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[2]), val[2]);
    RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[3]), val[3]);
    
    val[4] = lambda_dev[gloc]*dPt_dVmf_dthetat + lambda_dev[gloc+1]*dQt_dVmf_dthetat;
    val[5] = lambda_dev[gloc]*dPt_dVmf_dVmt    + lambda_dev[gloc+1]*dQt_dVmf_dVmt;
    val[6] = lambda_dev[gloc]*dPt_dVmf_dthetaf + lambda_dev[gloc+1]*dQt_dVmf_dthetaf;
    val[7] = lambda_dev[gloc]*dPt_dVmf_dVmf    + lambda_dev[gloc+1]*dQt_dVmf_dVmf;
                                                              
    RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[1], col[0]), val[4]);
    RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[1], col[1]), val[5]);
    RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[1], col[2]), val[6]);
    RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[1], col[3]), val[7]);
  });

  flps += ((56.0*EXAGO_FLOPS_SINOP) + (56.0*EXAGO_FLOPS_SINOP) + 462.0)*lineparams->nlineON;
  ierr = PetscLogFlops(flps);CHKERRQ(ierr);
  //  PetscPrintf(MPI_COMM_SELF,"Exit equality dense hessian\n");

  PetscFunctionReturn(0);
}

/**
 * @param[inout] HDD_dev Hessian matrix with size nxdense x nxdense
 */
PetscErrorCode OPFLOWComputeDenseInequalityConstraintHessian_PBPOLRAJAHIOP(OPFLOW opflow,const double *x_dev,const double* lambda_dev, double *HDD_dev)
{
  PbpolModelRajaHiop* pbpolrajahiop = reinterpret_cast<PbpolModelRajaHiop*>(opflow->model);
  LINEParamsRajaHiop  *lineparams=&pbpolrajahiop->lineparams;
  int                 nxsparse=2*opflow->ps->ngenON;
  const int           nxdense = 2*opflow->ps->nbus;
  PetscErrorCode      ierr;
  double              flps=0.0;

  PetscFunctionBegin;
  //  PetscPrintf(MPI_COMM_SELF,"Enter inequality dense hessian\n");

  /* Return if there are no inequality constraints */
  if(!lineparams->nlinelim) {
    //    PetscPrintf(MPI_COMM_SELF,"No inequality constraints. Exit inequality dense hessian\n");
    PetscFunctionReturn(0);
  }

  /** Although we will perform atomic operations on this view, we do not want
   * to use `RAJA::make_atomic_view` because explicit `RAJA::atomicAdd`s make
   * the atomic requirement clearer in the code */
  RAJA::View<double, RAJA::Layout<2>> HDD_view(HDD_dev, nxdense, nxdense);

  // Hessian from line contributions 
  double* Gff_arr = lineparams->Gff_dev_;
  double* Gtt_arr = lineparams->Gtt_dev_;
  double* Gft_arr = lineparams->Gft_dev_;
  double* Gtf_arr = lineparams->Gtf_dev_;

  double* Bff_arr = lineparams->Bff_dev_;
  double* Btt_arr = lineparams->Btt_dev_;
  double* Bft_arr = lineparams->Bft_dev_;
  double* Btf_arr = lineparams->Btf_dev_;

  int* xidxf = lineparams->xidxf_dev_;
  int* xidxt = lineparams->xidxt_dev_;
  int* gineqidx = lineparams->gineqidx_dev_;
  int* linelimidx = lineparams->linelimidx_dev_;

  RAJA::forall<exago_raja_exec>(RAJA::RangeSegment(0, lineparams->nlinelim),
    RAJA_LAMBDA (RAJA::Index_type i)
    {
      int j = linelimidx[i];
      int gloc;
      int    row[2],col[4];
      double val[8];

      double Pf,Qf,Pt,Qt;
      double thetaf=x_dev[xidxf[j]], Vmf=x_dev[xidxf[j]+1];
      double thetat=x_dev[xidxt[j]], Vmt=x_dev[xidxt[j]+1];
      double thetaft=thetaf-thetat;
      double thetatf=thetat-thetaf;
      double Gff = Gff_arr[j], Bff = Bff_arr[j];
      double Gft = Gft_arr[j], Bft = Bft_arr[j];
      double Gtf = Gtf_arr[j], Btf = Btf_arr[j];
      double Gtt = Gtt_arr[j], Btt = Btt_arr[j];

      Pf =  Gff*Vmf*Vmf + Vmf*Vmt*( Gft*cos(thetaft) + Bft*sin(thetaft));
      Qf = -Bff*Vmf*Vmf + Vmf*Vmt*(-Bft*cos(thetaft) + Gft*sin(thetaft));
    
      Pt =  Gtt*Vmt*Vmt + Vmt*Vmf*( Gtf*cos(thetatf) + Btf*sin(thetatf));
      Qt = -Btt*Vmt*Vmt + Vmt*Vmf*(-Btf*cos(thetatf) + Gtf*sin(thetatf));
      
      double dSf2_dPf, dSf2_dQf, dSt2_dPt, dSt2_dQt;
      
      dSf2_dPf = 2.*Pf;
      dSf2_dQf = 2.*Qf;
      dSt2_dPt = 2.*Pt;
      dSt2_dQt = 2.*Qt;
      
      double dPf_dthetaf,dPf_dVmf,dPf_dthetat,dPf_dVmt;
      double dQf_dthetaf,dQf_dVmf,dQf_dthetat,dQf_dVmt;
      double dPt_dthetaf,dPt_dVmf,dPt_dthetat,dPt_dVmt;
      double dQt_dthetaf,dQt_dVmf,dQt_dthetat,dQt_dVmt;
      
      dPf_dthetaf = 			Vmf*Vmt*(-Gft*sin(thetaft) + Bft*cos(thetaft));
      dPf_dVmf    = 2.*Gff*Vmf + 	Vmt*( Gft*cos(thetaft) + Bft*sin(thetaft));
      dPf_dthetat = 			Vmf*Vmt*( Gft*sin(thetaft) - Bft*cos(thetaft));
      dPf_dVmt    = 				Vmf*( Gft*cos(thetaft) + Bft*sin(thetaft));
      
      dQf_dthetaf = 			Vmf*Vmt*( Bft*sin(thetaft) + Gft*cos(thetaft));
      dQf_dVmf    = -2.*Bff*Vmf + 	Vmt*(-Bft*cos(thetaft) + Gft*sin(thetaft));
      dQf_dthetat = 			Vmf*Vmt*(-Bft*sin(thetaft) - Gft*cos(thetaft));
      dQf_dVmt    = 				Vmf*(-Bft*cos(thetaft) + Gft*sin(thetaft));
      
      dPt_dthetat = 			Vmt*Vmf*(-Gtf*sin(thetatf) + Btf*cos(thetatf));
      dPt_dVmt    = 2.*Gtt*Vmt + 	Vmf*( Gtf*cos(thetatf) + Btf*sin(thetatf));
      dPt_dthetaf = 			Vmt*Vmf*( Gtf*sin(thetatf) - Btf*cos(thetatf));
      dPt_dVmf    = 				Vmt*( Gtf*cos(thetatf) + Btf*sin(thetatf));
      
      dQt_dthetat = 			Vmt*Vmf*( Btf*sin(thetatf) + Gtf*cos(thetatf));
      dQt_dVmt    = -2.*Btt*Vmt + 	Vmf*(-Btf*cos(thetatf) + Gtf*sin(thetatf));
      dQt_dthetaf = 			Vmt*Vmf*(-Btf*sin(thetatf) - Gtf*cos(thetatf));
      dQt_dVmf    = 				Vmt*(-Btf*cos(thetatf) + Gtf*sin(thetatf));
      
      double d2Pf_dthetaf_dthetaf,d2Pf_dthetaf_dVmf,d2Pf_dthetaf_dthetat,d2Pf_dthetaf_dVmt;
      double d2Pf_dVmf_dthetaf,   d2Pf_dVmf_dVmf,   d2Pf_dVmf_dthetat,   d2Pf_dVmf_dVmt;
      double d2Pf_dthetat_dthetaf,d2Pf_dthetat_dVmf,d2Pf_dthetat_dthetat,d2Pf_dthetat_dVmt;
      double d2Pf_dVmt_dthetaf,   d2Pf_dVmt_dVmf,   d2Pf_dVmt_dthetat,   d2Pf_dVmt_dVmt;
      
      /* dPf_dthetaf = Vmf*Vmt*(-Gft*sin(thetaft) + Bft*cos(thetaft)); */
      d2Pf_dthetaf_dthetaf = -Vmf*Vmt*(Gft*cos(thetaft) + Bft*sin(thetaft));
      d2Pf_dthetaf_dVmf    =     Vmt*(-Gft*sin(thetaft) + Bft*cos(thetaft));
      d2Pf_dthetaf_dthetat =  Vmf*Vmt*(Gft*cos(thetaft) + Bft*sin(thetaft));
      d2Pf_dthetaf_dVmt    =     Vmf*(-Gft*sin(thetaft) + Bft*cos(thetaft));
      
      /* dPf_Vmf  = 2*Gff*Vmf + Vmt*(Gft*cos(thetaft) + Bft*sin(thetaft)); */
      d2Pf_dVmf_dthetaf    =  Vmt*(-Gft*sin(thetaft) + Bft*cos(thetaft)); 
      d2Pf_dVmf_dVmf       =  2*Gff;
      d2Pf_dVmf_dthetat    =  Vmt*(Gft*sin(thetaft) - Bft*cos(thetaft));
      d2Pf_dVmf_dVmt       =      (Gft*cos(thetaft) + Bft*sin(thetaft));
      
      /* dPf_dthetat = Vmf*Vmt*(Gft*sin(thetaft) - Bft*cos(thetaft)); */
      d2Pf_dthetat_dthetaf = Vmf*Vmt*(Gft*cos(thetaft)  + Bft*sin(thetaft));
      d2Pf_dthetat_dVmf    =     Vmt*(Gft*sin(thetaft)  - Bft*cos(thetaft));
      d2Pf_dthetat_dthetat = Vmf*Vmt*(-Gft*cos(thetaft) - Bft*sin(thetaft));
      d2Pf_dthetat_dVmt    =     Vmf*(Gft*sin(thetaft)  - Bft*cos(thetaft));
      
      /* dPf_dVmt = Vmf*(Gft*cos(thetaft) + Bft*sin(thetaft)); */
      d2Pf_dVmt_dthetaf    = Vmf*(-Gft*sin(thetaft) + Bft*cos(thetaft)); 
      d2Pf_dVmt_dVmf       =      (Gft*cos(thetaft) + Bft*sin(thetaft));
      d2Pf_dVmt_dthetat    = Vmf*(Gft*sin(thetaft) - Bft*cos(thetaft));
      d2Pf_dVmt_dVmt       = 0.0;
      
      double d2Qf_dthetaf_dthetaf,d2Qf_dthetaf_dVmf,d2Qf_dthetaf_dthetat,d2Qf_dthetaf_dVmt;
      double d2Qf_dVmf_dthetaf,   d2Qf_dVmf_dVmf,   d2Qf_dVmf_dthetat,   d2Qf_dVmf_dVmt;
      double d2Qf_dthetat_dthetaf,d2Qf_dthetat_dVmf,d2Qf_dthetat_dthetat,d2Qf_dthetat_dVmt;
      double d2Qf_dVmt_dthetaf,   d2Qf_dVmt_dVmf,   d2Qf_dVmt_dthetat,   d2Qf_dVmt_dVmt;
      
      /* dQf_dthetaf = Vmf*Vmt*(Bft*sin(thetaft) + Gft*cos(thetaft)); */
      d2Qf_dthetaf_dthetaf = Vmf*Vmt*(Bft*cos(thetaft)  - Gft*sin(thetaft));
      d2Qf_dthetaf_dVmf    =     Vmt*(Bft*sin(thetaft)  + Gft*cos(thetaft));
      d2Qf_dthetaf_dthetat = Vmf*Vmt*(-Bft*cos(thetaft) + Gft*sin(thetaft));
      d2Qf_dthetaf_dVmt    =     Vmf*(Bft*sin(thetaft) + Gft*cos(thetaft));
      
      /* dQf_dVmf = -2*Bff*Vmf + Vmt*(-Bft*cos(thetaft) + Gft*sin(thetaft)); */
      d2Qf_dVmf_dthetaf    =  Vmt*(Bft*sin(thetaft) + Gft*cos(thetaft)); 
      d2Qf_dVmf_dVmf       = -2*Bff;
      d2Qf_dVmf_dthetat    =  Vmt*(-Bft*sin(thetaft) - Gft*cos(thetaft));
      d2Qf_dVmf_dVmt       =      (-Bft*cos(thetaft) + Gft*sin(thetaft));
      
      /* dQf_dthetat = Vmf*Vmt*(-Bft*sin(thetaft) - Gft*cos(thetaft)); */
      d2Qf_dthetat_dthetaf = Vmf*Vmt*(-Bft*cos(thetaft) + Gft*sin(thetaft));
      d2Qf_dthetat_dVmf    =     Vmt*(-Bft*sin(thetaft) - Gft*cos(thetaft));
      d2Qf_dthetat_dthetat = Vmf*Vmt*( Bft*cos(thetaft) - Gft*sin(thetaft));
      d2Qf_dthetat_dVmt    =     Vmf*(-Bft*sin(thetaft) - Gft*cos(thetaft));
      
      /* dQf_dVmt = Vmf*(-Bft*cos(thetaft) + Gft*sin(thetaft)); */
      d2Qf_dVmt_dthetaf    = Vmf*(Bft*sin(thetaft) + Gft*cos(thetaft)); 
      d2Qf_dVmt_dVmf       =    (-Bft*cos(thetaft) + Gft*sin(thetaft));
      d2Qf_dVmt_dthetat    = Vmf*(-Bft*sin(thetaft) - Gft*cos(thetaft));
      d2Qf_dVmt_dVmt       = 0.0;
      
      double d2Pt_dthetat_dthetat,d2Pt_dthetat_dVmt,d2Pt_dthetat_dthetaf,d2Pt_dthetat_dVmf;
      double d2Pt_dVmt_dthetat,   d2Pt_dVmt_dVmt,   d2Pt_dVmt_dthetaf,   d2Pt_dVmt_dVmf;
      double d2Pt_dthetaf_dthetat,d2Pt_dthetaf_dVmt,d2Pt_dthetaf_dthetaf,d2Pt_dthetaf_dVmf;
      double d2Pt_dVmf_dthetat,   d2Pt_dVmf_dVmt,   d2Pt_dVmf_dthetaf,   d2Pt_dVmf_dVmf;
      
      /* dPt_dthetat = Vmf*Vmt*(-Gtf*sin(thetatf) + Btf*cos(thetatf)); */
      d2Pt_dthetat_dthetat = Vmf*Vmt*(-Gtf*cos(thetatf) - Btf*sin(thetatf));
      d2Pt_dthetat_dVmt    =     Vmf*(-Gtf*sin(thetatf) + Btf*cos(thetatf));
      d2Pt_dthetat_dthetaf =  Vmf*Vmt*(Gtf*cos(thetatf) + Btf*sin(thetatf));
      d2Pt_dthetat_dVmf    =     Vmt*(-Gtf*sin(thetatf) + Btf*cos(thetatf));
      
      /* dPt_Vmt  = 2*Gtt*Vmt + Vmf*(Gtf*cos(thetatf) + Btf*sin(thetatf)); */
      d2Pt_dVmt_dthetat    =  Vmf*(-Gtf*sin(thetatf) + Bft*cos(thetatf)); 
      d2Pt_dVmt_dVmt       =  2*Gtt;
      d2Pt_dVmt_dthetaf    =  Vmf*(Gtf*sin(thetatf) - Btf*cos(thetatf));
      d2Pt_dVmt_dVmf       =      (Gtf*cos(thetatf) + Btf*sin(thetatf));
      
      /* dPt_dthetaf = Vmf*Vmt*(Gtf*sin(thetatf) - Btf*cos(thetatf)); */
      d2Pt_dthetaf_dthetat = Vmf*Vmt*(Gtf*cos(thetatf)  + Btf*sin(thetatf));
      d2Pt_dthetaf_dVmt    =     Vmf*(Gtf*sin(thetatf)  - Btf*cos(thetatf));
      d2Pt_dthetaf_dthetaf = Vmf*Vmt*(-Gtf*cos(thetatf) - Btf*sin(thetatf));
      d2Pt_dthetaf_dVmf    =     Vmt*(Gtf*sin(thetatf)  - Btf*cos(thetatf));
      
      /* dPt_dVmf = Vmt*(Gtf*cos(thetatf) + Btf*sin(thetatf)); */
      d2Pt_dVmf_dthetat    = Vmt*(-Gtf*sin(thetatf) + Btf*cos(thetatf)); 
      d2Pt_dVmf_dVmt       =     (Gtf*cos(thetatf) + Btf*sin(thetatf));
      d2Pt_dVmf_dthetaf    = Vmt*(Gtf*sin(thetatf) - Btf*cos(thetatf));
      d2Pt_dVmf_dVmf       = 0.0;
      
      double d2Qt_dthetaf_dthetaf,d2Qt_dthetaf_dVmf,d2Qt_dthetaf_dthetat,d2Qt_dthetaf_dVmt;
      double d2Qt_dVmf_dthetaf,   d2Qt_dVmf_dVmf,   d2Qt_dVmf_dthetat,   d2Qt_dVmf_dVmt;
      double d2Qt_dthetat_dthetaf,d2Qt_dthetat_dVmf,d2Qt_dthetat_dthetat,d2Qt_dthetat_dVmt;
      double d2Qt_dVmt_dthetaf,   d2Qt_dVmt_dVmf,   d2Qt_dVmt_dthetat,   d2Qt_dVmt_dVmt;
      
      /* dQt_dthetat = Vmf*Vmt*(Btf*sin(thetatf) + Gtf*cos(thetatf)); */
      d2Qt_dthetat_dthetat = Vmf*Vmt*(Btf*cos(thetatf)  - Gtf*sin(thetatf));
      d2Qt_dthetat_dVmt    =     Vmf*(Btf*sin(thetatf)  + Gtf*cos(thetatf));
      d2Qt_dthetat_dthetaf = Vmf*Vmt*(-Btf*cos(thetatf) + Gtf*sin(thetatf));
      d2Qt_dthetat_dVmf    =     Vmf*(Btf*sin(thetatf) + Gtf*cos(thetatf));
      
      /* dQt_dVmt = -2*Btt*Vmt + Vmf*(-Btf*cos(thetatf) + Gtf*sin(thetatf)); */
      d2Qt_dVmt_dthetat    =  Vmf*(Btf*sin(thetatf) + Gtf*cos(thetatf)); 
      d2Qt_dVmt_dVmt       = -2*Btt;
      d2Qt_dVmt_dthetaf    =  Vmf*(-Btf*sin(thetatf) + Gtf*cos(thetatf));
      d2Qt_dVmt_dVmf       =      (-Btf*cos(thetatf) + Gtf*sin(thetatf));
      
      /* dQt_dthetaf = Vmf*Vmt*(-Btf*sin(thetatf) - Gtf*cos(thetatf)); */
      d2Qt_dthetaf_dthetat = Vmf*Vmt*(-Btf*cos(thetatf) + Gtf*sin(thetatf));
      d2Qt_dthetaf_dVmt    =     Vmf*(-Btf*sin(thetatf) - Gtf*cos(thetatf));
      d2Qt_dthetaf_dthetaf = Vmf*Vmt*( Btf*cos(thetatf) - Gtf*sin(thetatf));
      d2Qt_dthetaf_dVmf    =     Vmt*(-Btf*sin(thetatf) - Gtf*cos(thetatf));
      
      /* dQt_dVmf = Vmt*(-Btf*cos(thetatf) + Gtf*sin(thetatf)); */
      d2Qt_dVmf_dthetat    = Vmt*(Btf*sin(thetatf) + Gtf*cos(thetatf)); 
      d2Qt_dVmf_dVmt       =    (-Btf*cos(thetatf) + Gtf*sin(thetatf));
      d2Qt_dVmf_dthetaf    = Vmt*(-Btf*sin(thetatf) - Gtf*cos(thetatf));
      d2Qt_dVmf_dVmf       = 0.0;
      
      double d2Sf2_dthetaf_dthetaf=0.0,d2Sf2_dthetaf_dVmf=0.0,d2Sf2_dthetaf_dthetat=0.0,d2Sf2_dthetaf_dVmt=0.0;
      double d2St2_dthetaf_dthetaf=0.0,d2St2_dthetaf_dVmf=0.0,d2St2_dthetaf_dthetat=0.0,d2St2_dthetaf_dVmt=0.0;
      
      d2Sf2_dthetaf_dthetaf = 2*dPf_dthetaf*dPf_dthetaf + dSf2_dPf*d2Pf_dthetaf_dthetaf +  2*dQf_dthetaf*dQf_dthetaf + dSf2_dQf*d2Qf_dthetaf_dthetaf;
      d2Sf2_dthetaf_dVmf = 2*dPf_dVmf*dPf_dthetaf + dSf2_dPf*d2Pf_dthetaf_dVmf +  2*dQf_dVmf*dQf_dthetaf + dSf2_dQf*d2Qf_dthetaf_dVmf;
      d2Sf2_dthetaf_dthetat = 2*dPf_dthetat*dPf_dthetaf + dSf2_dPf*d2Pf_dthetaf_dthetat +  2*dQf_dthetat*dQf_dthetaf + dSf2_dQf*d2Qf_dthetaf_dthetat;
      d2Sf2_dthetaf_dVmt = 2*dPf_dVmt*dPf_dthetaf + dSf2_dPf*d2Pf_dthetaf_dVmt +  2*dQf_dVmt*dQf_dthetaf + dSf2_dQf*d2Qf_dthetaf_dVmt;
      
      d2St2_dthetaf_dthetaf = 2*dPt_dthetaf*dPt_dthetaf + dSt2_dPt*d2Pt_dthetaf_dthetaf +  2*dQt_dthetaf*dQt_dthetaf + dSt2_dQt*d2Qt_dthetaf_dthetaf;
      d2St2_dthetaf_dVmf = 2*dPt_dVmf*dPt_dthetaf + dSt2_dPt*d2Pt_dthetaf_dVmf +  2*dQt_dVmf*dQt_dthetaf + dSt2_dQt*d2Qt_dthetaf_dVmf;
      d2St2_dthetaf_dthetat = 2*dPt_dthetat*dPt_dthetaf + dSt2_dPt*d2Pt_dthetaf_dthetat +  2*dQt_dthetat*dQt_dthetaf + dSt2_dQt*d2Qt_dthetaf_dthetat;
      d2St2_dthetaf_dVmt = 2*dPt_dVmt*dPt_dthetaf + dSt2_dPt*d2Pt_dthetaf_dVmt +  2*dQt_dVmt*dQt_dthetaf + dSt2_dQt*d2Qt_dthetaf_dVmt;
      
      val[0] = val[1] = val[2] = val[3] = 0.0;

      row[0] = xidxf[j]     - nxsparse;
      col[0] = xidxf[j]     - nxsparse;
      col[1] = xidxf[j] + 1 - nxsparse;
      col[2] = xidxt[j]     - nxsparse;
      col[3] = xidxt[j] + 1 - nxsparse;
      
      gloc = gineqidx[i];

      val[0] = lambda_dev[gloc]*d2Sf2_dthetaf_dthetaf + lambda_dev[gloc+1]*d2St2_dthetaf_dthetaf;
      val[1] = lambda_dev[gloc]*d2Sf2_dthetaf_dVmf + lambda_dev[gloc+1]*d2St2_dthetaf_dVmf;
      val[2] = lambda_dev[gloc]*d2Sf2_dthetaf_dthetat + lambda_dev[gloc+1]*d2St2_dthetaf_dthetat;
      val[3] = lambda_dev[gloc]*d2Sf2_dthetaf_dVmt + lambda_dev[gloc+1]*d2St2_dthetaf_dVmt;
      
      RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[0]), val[0]);
      RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[1]), val[1]);
      RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[2]), val[2]);
      RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[3]), val[3]);

      double d2Sf2_dVmf_dthetaf,d2Sf2_dVmf_dVmf,d2Sf2_dVmf_dthetat,d2Sf2_dVmf_dVmt;
      double d2St2_dVmf_dthetaf,d2St2_dVmf_dVmf,d2St2_dVmf_dthetat,d2St2_dVmf_dVmt;
      
      d2Sf2_dVmf_dthetaf = 2*dPf_dthetaf*dPf_dVmf + dSf2_dPf*d2Pf_dVmf_dthetaf +  2*dQf_dthetaf*dQf_dVmf + dSf2_dQf*d2Qf_dVmf_dthetaf;
      d2Sf2_dVmf_dVmf = 2*dPf_dVmf*dPf_dVmf + dSf2_dPf*d2Pf_dVmf_dVmf +  2*dQf_dVmf*dQf_dVmf + dSf2_dQf*d2Qf_dVmf_dVmf;
      d2Sf2_dVmf_dthetat = 2*dPf_dthetat*dPf_dVmf + dSf2_dPf*d2Pf_dVmf_dthetat +  2*dQf_dthetat*dQf_dVmf + dSf2_dQf*d2Qf_dVmf_dthetat;
      d2Sf2_dVmf_dVmt = 2*dPf_dVmt*dPf_dVmf + dSf2_dPf*d2Pf_dVmf_dVmt +  2*dQf_dVmt*dQf_dVmf + dSf2_dQf*d2Qf_dVmf_dVmt;
      
      d2St2_dVmf_dthetaf = 2*dPt_dthetaf*dPt_dVmf + dSt2_dPt*d2Pt_dVmf_dthetaf +  2*dQt_dthetaf*dQt_dVmf + dSt2_dQt*d2Qt_dVmf_dthetaf;
      d2St2_dVmf_dVmf = 2*dPt_dVmf*dPt_dVmf + dSt2_dPt*d2Pt_dVmf_dVmf +  2*dQt_dVmf*dQt_dVmf + dSt2_dQt*d2Qt_dVmf_dVmf;
      d2St2_dVmf_dthetat = 2*dPt_dthetat*dPt_dVmf + dSt2_dPt*d2Pt_dVmf_dthetat +  2*dQt_dthetat*dQt_dVmf + dSt2_dQt*d2Qt_dVmf_dthetat;
      d2St2_dVmf_dVmt = 2*dPt_dVmt*dPt_dVmf + dSt2_dPt*d2Pt_dVmf_dVmt +  2*dQt_dVmt*dQt_dVmf + dSt2_dQt*d2Qt_dVmf_dVmt;
      
      val[0] = val[1] = val[2] = val[3] = 0.0;
      col[0] = xidxf[j]     - nxsparse;
      col[1] = xidxf[j] + 1 - nxsparse;
      col[2] = xidxt[j]     - nxsparse;
      col[3] = xidxt[j] + 1 - nxsparse;
      
      row[0] = xidxf[j] + 1 - nxsparse;

      val[0] = lambda_dev[gloc]*d2Sf2_dVmf_dthetaf + lambda_dev[gloc+1]*d2St2_dVmf_dthetaf;
      val[1] = lambda_dev[gloc]*d2Sf2_dVmf_dVmf + lambda_dev[gloc+1]*d2St2_dVmf_dVmf;
      val[2] = lambda_dev[gloc]*d2Sf2_dVmf_dthetat + lambda_dev[gloc+1]*d2St2_dVmf_dthetat;
      val[3] = lambda_dev[gloc]*d2Sf2_dVmf_dVmt + lambda_dev[gloc+1]*d2St2_dVmf_dVmt;

      RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[0]), val[0]);
      RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[1]), val[1]);
      RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[2]), val[2]);
      RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[3]), val[3]);

      //    ierr = MatSetValues(H,1,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);
      
      double d2Sf2_dthetat_dthetaf,d2Sf2_dthetat_dVmf,d2Sf2_dthetat_dthetat,d2Sf2_dthetat_dVmt;
      double d2St2_dthetat_dthetaf,d2St2_dthetat_dVmf,d2St2_dthetat_dthetat,d2St2_dthetat_dVmt;
      
      d2Sf2_dthetat_dthetaf = 2*dPf_dthetaf*dPf_dthetat + dSf2_dPf*d2Pf_dthetat_dthetaf +  2*dQf_dthetat*dQf_dthetaf + dSf2_dQf*d2Qf_dthetat_dthetaf;
      d2Sf2_dthetat_dVmf = 2*dPf_dVmf*dPf_dthetat + dSf2_dPf*d2Pf_dthetat_dVmf +  2*dQf_dthetat*dQf_dVmf + dSf2_dQf*d2Qf_dthetat_dVmf;
      d2Sf2_dthetat_dthetat = 2*dPf_dthetat*dPf_dthetat + dSf2_dPf*d2Pf_dthetat_dthetat +  2*dQf_dthetat*dQf_dthetat + dSf2_dQf*d2Qf_dthetat_dthetat;
      d2Sf2_dthetat_dVmt = 2*dPf_dVmt*dPf_dthetat + dSf2_dPf*d2Pf_dthetat_dVmt +  2*dQf_dthetat*dQf_dVmt + dSf2_dQf*d2Qf_dthetat_dVmt;
      
      d2St2_dthetat_dthetaf = 2*dPt_dthetaf*dPt_dthetat + dSt2_dPt*d2Pt_dthetat_dthetaf +  2*dQt_dthetaf*dQt_dthetat + dSt2_dQt*d2Qt_dthetat_dthetaf;
      d2St2_dthetat_dVmf = 2*dPt_dVmf*dPt_dthetat + dSt2_dPt*d2Pt_dthetat_dVmf +  2*dQt_dVmf*dQt_dthetat + dSt2_dQt*d2Qt_dthetat_dVmf;
      d2St2_dthetat_dthetat = 2*dPt_dthetat*dPt_dthetat + dSt2_dPt*d2Pt_dthetat_dthetat +  2*dQt_dthetat*dQt_dthetat + dSt2_dQt*d2Qt_dthetat_dthetat;
      d2St2_dthetat_dVmt = 2*dPt_dVmt*dPt_dthetat + dSt2_dPt*d2Pt_dthetat_dVmt +  2*dQt_dVmt*dQt_dthetat + dSt2_dQt*d2Qt_dthetat_dVmt;
      
      val[0] = val[1] = val[2] = val[3] = 0.0;
      
      col[0] = xidxf[j]     - nxsparse;
      col[1] = xidxf[j] + 1 - nxsparse;
      col[2] = xidxt[j]     - nxsparse;
      col[3] = xidxt[j] + 1 - nxsparse;
      
      row[0] = xidxt[j]     - nxsparse;
      
      val[0] = lambda_dev[gloc]*d2Sf2_dthetat_dthetaf + lambda_dev[gloc+1]*d2St2_dthetat_dthetaf;
      val[1] = lambda_dev[gloc]*d2Sf2_dthetat_dVmf + lambda_dev[gloc+1]*d2St2_dthetat_dVmf;
      val[2] = lambda_dev[gloc]*d2Sf2_dthetat_dthetat + lambda_dev[gloc+1]*d2St2_dthetat_dthetat;
      val[3] = lambda_dev[gloc]*d2Sf2_dthetat_dVmt + lambda_dev[gloc+1]*d2St2_dthetat_dVmt;
      
      RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[0]), val[0]);
      RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[1]), val[1]);
      RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[2]), val[2]);
      RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[3]), val[3]);

      //    ierr = MatSetValues(H,1,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);
      
      double d2Sf2_dVmt_dthetaf,d2Sf2_dVmt_dVmf,d2Sf2_dVmt_dthetat,d2Sf2_dVmt_dVmt;
      double d2St2_dVmt_dthetaf,d2St2_dVmt_dVmf,d2St2_dVmt_dthetat,d2St2_dVmt_dVmt;
      
      d2Sf2_dVmt_dthetaf = 2*dPf_dthetaf*dPf_dVmt + dSf2_dPf*d2Pf_dVmt_dthetaf +  2*dQf_dthetaf*dQf_dVmt + dSf2_dQf*d2Qf_dVmt_dthetaf;
      d2Sf2_dVmt_dVmf = 2*dPf_dVmf*dPf_dVmt + dSf2_dPf*d2Pf_dVmt_dVmf +  2*dQf_dVmf*dQf_dVmt + dSf2_dQf*d2Qf_dVmt_dVmf;
      d2Sf2_dVmt_dthetat = 2*dPf_dthetat*dPf_dVmt + dSf2_dPf*d2Pf_dVmt_dthetat +  2*dQf_dthetat*dQf_dVmt + dSf2_dQf*d2Qf_dVmt_dthetat;
      d2Sf2_dVmt_dVmt = 2*dPf_dVmt*dPf_dVmt + dSf2_dPf*d2Pf_dVmt_dVmt +  2*dQf_dVmt*dQf_dVmt + dSf2_dQf*d2Qf_dVmt_dVmt;
      
      d2St2_dVmt_dthetaf = 2*dPt_dthetaf*dPt_dVmt + dSt2_dPt*d2Pt_dVmt_dthetaf +  2*dQt_dthetaf*dQt_dVmt + dSt2_dQt*d2Qt_dVmt_dthetaf;
      d2St2_dVmt_dVmf = 2*dPt_dVmf*dPt_dVmt + dSt2_dPt*d2Pt_dVmt_dVmf +  2*dQt_dVmf*dQt_dVmt + dSt2_dQt*d2Qt_dVmt_dVmf;
      d2St2_dVmt_dthetat = 2*dPt_dthetat*dPt_dVmt + dSt2_dPt*d2Pt_dVmt_dthetat +  2*dQt_dthetat*dQt_dVmt + dSt2_dQt*d2Qt_dVmt_dthetat;
      d2St2_dVmt_dVmt = 2*dPt_dVmt*dPt_dVmt + dSt2_dPt*d2Pt_dVmt_dVmt +  2*dQt_dVmt*dQt_dVmt + dSt2_dQt*d2Qt_dVmt_dVmt;
      
      val[0] = val[1] = val[2] = val[3] = 0.0;

      row[0] = xidxt[j] + 1 - nxsparse;
      col[0] = xidxf[j]     - nxsparse;
      col[1] = xidxf[j] + 1 - nxsparse;
      col[2] = xidxt[j]     - nxsparse;
      col[3] = xidxt[j] + 1 - nxsparse;
      
      val[0] = lambda_dev[gloc]*d2Sf2_dVmt_dthetaf + lambda_dev[gloc+1]*d2St2_dVmt_dthetaf;
      val[1] = lambda_dev[gloc]*d2Sf2_dVmt_dVmf + lambda_dev[gloc+1]*d2St2_dVmt_dVmf;
      val[2] = lambda_dev[gloc]*d2Sf2_dVmt_dthetat + lambda_dev[gloc+1]*d2St2_dVmt_dthetat;
      val[3] = lambda_dev[gloc]*d2Sf2_dVmt_dVmt + lambda_dev[gloc+1]*d2St2_dVmt_dVmt;
    
      RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[0]), val[0]);
      RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[1]), val[1]);
      RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[2]), val[2]);
      RAJA::atomicAdd<exago_raja_atomic>(&HDD_view(row[0], col[3]), val[3]);
    }
  );
  //  PetscPrintf(MPI_COMM_SELF,"Exit inequality dense hessian\n");
  flps += (972+(92*EXAGO_FLOPS_COSOP)+(92*EXAGO_FLOPS_SINOP))*lineparams->nlinelim;
  ierr = PetscLogFlops(flps);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeDenseHessian_PBPOLRAJAHIOP(OPFLOW opflow,const double *x_dev, const double *lambda_dev,double *HDD_dev)
{
  PetscErrorCode ierr;
  int nxdense = 2*opflow->ps->nbus;

  if(!HDD_dev) PetscFunctionReturn(0);

  RAJA::forall<exago_raja_exec>(RAJA::RangeSegment(0, nxdense*nxdense),
    RAJA_LAMBDA (RAJA::Index_type i)
    {
      HDD_dev[i] = 0.0;
    }
  );

  /* Equality constraint Hessian */
  ierr = OPFLOWComputeDenseEqualityConstraintHessian_PBPOLRAJAHIOP(opflow,x_dev,lambda_dev,HDD_dev);CHKERRQ(ierr);

  if(opflow->nconineq) {
    ierr = OPFLOWComputeDenseInequalityConstraintHessian_PBPOLRAJAHIOP(opflow,x_dev,lambda_dev+opflow->nconeq,HDD_dev);CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}
#endif
