#include <exago_config.h>

#if defined(EXAGO_HAVE_RAJA) 

#include <umpire/Allocator.hpp>
#include <umpire/ResourceManager.hpp>

#include <RAJA/RAJA.hpp>

#include <private/opflowimpl.h>
#include "pbpolrajakernels.hpp"
#include "pbpolraja.hpp"

/************* NOTE ***********************/
/* No Load loss or power imbalance variables considered yet */
/********************************************/

/* Functions to create and destroy data arrays for different
   component classes
*/
int BUSParams::destroy()
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

  return 0;
}

/* Create data for buses that is used in different computations */
int BUSParams::allocate(OPFLOW opflow)
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
  d_allocator_ = resmgr.getAllocator("DEVICE");

  // Allocate data on the host
  isref      = paramAlloc<int>(h_allocator_, nbus);
  isisolated = paramAlloc<int>(h_allocator_, nbus);
  ispvpq     = paramAlloc<int>(h_allocator_, nbus);

  vmin = paramAlloc<double>(h_allocator_, nbus);
  vmax = paramAlloc<double>(h_allocator_, nbus);
  vmin = paramAlloc<double>(h_allocator_, nbus);
  va   = paramAlloc<double>(h_allocator_, nbus);
  vm   = paramAlloc<double>(h_allocator_, nbus);
  gl   = paramAlloc<double>(h_allocator_, nbus);
  bl   = paramAlloc<double>(h_allocator_, nbus);

  xidx = paramAlloc<int>(h_allocator_, nbus);
  gidx = paramAlloc<int>(h_allocator_, nbus);

  /* Populate the arrays */
  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
    
    xidx[i] = loc;
    gidx[i] = gloc;

    if(bus->ide == REF_BUS) isref[i] = 1;
    else if(bus->ide == ISOLATED_BUS) isisolated[i] = 1;
    else ispvpq[i] = 1;

    vmin[i] = bus->Vmin;
    vmax[i] = bus->Vmax;
    vm[i]   = bus->vm;
    va[i]   = bus->va;
    gl[i]   = bus->gl;
    bl[i]   = bus->bl;

    gloc += 2;
  }
    
  // Allocate data on the device
  isref_dev_      = paramAlloc<int>(d_allocator_, nbus);
  isisolated_dev_ = paramAlloc<int>(d_allocator_, nbus);
  ispvpq_dev_     = paramAlloc<int>(d_allocator_, nbus);

  vmin_dev_ = paramAlloc<double>(d_allocator_, nbus);
  vmax_dev_ = paramAlloc<double>(d_allocator_, nbus);
  vmin_dev_ = paramAlloc<double>(d_allocator_, nbus);
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
  resmgr.copy(vmin_dev_, vmin);
  resmgr.copy(va_dev_  , va);
  resmgr.copy(vm_dev_  , vm);
  resmgr.copy(gl_dev_  , gl);
  resmgr.copy(bl_dev_  , bl);

  resmgr.copy(xidx_dev_, xidx);
  resmgr.copy(gidx_dev_, gidx);

  PetscFunctionReturn(0);
}


int LINEParams::destroy()
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

  h_allocator_.deallocate(linelimidx);

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

  d_allocator_.deallocate(linelimidx_dev_);

  return 0;
}


/* Create data for lines that is used in different computations */
int LINEParams::allocate(OPFLOW opflow)
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
  for(i=0; i < ps->nline; i++) {
    line = &ps->line[i];

    if(!line->status || line->rateA > 1e5) continue;
    nlinelim++;
  }

  /* Allocate data arrays */
  auto& resmgr = umpire::ResourceManager::getInstance();
  h_allocator_ = resmgr.getAllocator("HOST");
  d_allocator_ = resmgr.getAllocator("DEVICE");

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

  linelimidx = paramAlloc<int>(h_allocator_, nlinelim);

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

    ierr = PSBUSGetVariableLocation(busf,&xidxf[linei]);CHKERRQ(ierr);
    ierr = PSBUSGetVariableLocation(bust,&xidxt[linei]);CHKERRQ(ierr);

    /* 
       Each bus has two equality (balance) constraints, hence the use of coefficient 2
       to map the location of the equality constraint for the bus
    */
    geqidxf[linei] = 2*busf->internal_i;
    geqidxt[linei] = 2*bust->internal_i;

    if(line->rateA < 1e5) {
      gbineqidx[linelimi] = gbloc;
      gineqidx[linelimi] = gloc;
      linelimidx[linelimi] = linei;
      linelimi++;
      gbloc += 2;
      gloc += 2;
    }
    linei++;
  }

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

  linelimidx_dev_ = paramAlloc<int>(d_allocator_, nlinelim);

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

  resmgr.copy(linelimidx_dev_, linelimidx);

  return 0;
}

int LOADParams::destroy()
{
  h_allocator_.deallocate(pl);
  h_allocator_.deallocate(ql);
  h_allocator_.deallocate(xidx);
  h_allocator_.deallocate(gidx);

  d_allocator_.deallocate(pl_dev_);
  d_allocator_.deallocate(ql_dev_);
  d_allocator_.deallocate(xidx_dev_);
  d_allocator_.deallocate(gidx_dev_);

  return 0;
}

/* Create data for loads that is used in different computations */
int LOADParams::allocate(OPFLOW opflow)
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
  d_allocator_ = resmgr.getAllocator("DEVICE");

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

  return (0);
}


int GENParams::destroy()
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

  return 0;
}


/* Create data for generators that is used in different computations */
int GENParams::allocate(OPFLOW opflow)
{
  PS             ps=opflow->ps;
  PetscInt       loc,gloc=0,geni=0;
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
  d_allocator_ = resmgr.getAllocator("DEVICE");

  cost_alpha = paramAlloc<double>(h_allocator_, ngenON);
  cost_beta  = paramAlloc<double>(h_allocator_, ngenON);
  cost_gamma = paramAlloc<double>(h_allocator_, ngenON);

  pt = paramAlloc<double>(h_allocator_, ngenON);
  pb = paramAlloc<double>(h_allocator_, ngenON);
  qt = paramAlloc<double>(h_allocator_, ngenON);
  qb = paramAlloc<double>(h_allocator_, ngenON);

  xidx = paramAlloc<int>(h_allocator_, ngenON);
  gidx = paramAlloc<int>(h_allocator_, ngenON);

  /* Populate data on the host */
  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
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

      xidx[geni]       = loc;
      gidx[geni]       = gloc;

      geni++;
    }
    gloc += 2;
  }

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

  return 0;
}

/* The calculations for different routines start from here */

/** CONSTRAINT BOUNDS  **/
PetscErrorCode OPFLOWSetConstraintBounds_PBPOLRAJA(OPFLOW opflow,Vec Gl, Vec Gu)
{
  PetscErrorCode ierr;
//  PBPOLRAJA         pbpolraja=(PBPOLRAJA)opflow->model;
  PbpolModelRaja* pbpolraja = reinterpret_cast<PbpolModelRaja*>(opflow->model);
  LINEParams     *lineparams=&pbpolraja->lineparams;
  PetscScalar    *gl,*gu;
  PS             ps=opflow->ps;

  PetscFunctionBegin;

  ierr = VecSet(Gl,0.0);
  ierr = VecSet(Gu,0.0);

  if(!opflow->nconineq) PetscFunctionReturn(0);

  ierr = VecGetArray(Gl,&gl);CHKERRQ(ierr);
  ierr = VecGetArray(Gu,&gu);CHKERRQ(ierr);

  auto& resmgr = umpire::ResourceManager::getInstance();
  auto hostallocator = resmgr.getAllocator("HOST");
  auto devallocator = resmgr.getAllocator("DEVICE");

  // sizes for gl == gu
  int sz; 
  VecGetSize(Gl, &sz);

  registerWith(gl, sz, resmgr, hostallocator);
  registerWith(gu, sz, resmgr, hostallocator);
  double* d_gl = paramAlloc<double>(devallocator, sz);
  double* d_gu = paramAlloc<double>(devallocator, sz);

  int* gbineqidx = lineparams->gbineqidx_dev_;
  double* rateA = lineparams->rateA_dev_;
  int* linelimidx = lineparams->linelimidx_dev_;
  // Later, we should copy this to device... 
  double MVAbase = ps->MVAbase;
  RAJA::forall< RAJA::cuda_exec<128> >(RAJA::RangeSegment(0, lineparams->nlinelim),
    [=] __device__ (RAJA::Index_type i)
    {
      int j = linelimidx[i];
      d_gl[gbineqidx[i]]   = 0.0;
      d_gu[gbineqidx[i]]   = (rateA[j]/MVAbase)*(rateA[j]/MVAbase);
      d_gl[gbineqidx[i]+1] = 0.0;
      d_gu[gbineqidx[i]+1] = (rateA[j]/MVAbase)*(rateA[j]/MVAbase);
    }
  );

  resmgr.copy(gu, d_gu);
  resmgr.copy(gl, d_gl);

  ierr = VecGetArray(Gl,&gl);CHKERRQ(ierr);
  ierr = VecGetArray(Gu,&gu);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/** EQUALITY CONSTRAINTS

Notes from RAJA Hackathon:
  - Are the four loops in this function dependent on each other?
  - Can we run each of them in parallel?
    - async raja kernels
    - potentially combine into one larger kernel wiht custom raja kernel policy

*/
PetscErrorCode OPFLOWComputeEqualityConstraints_PBPOLRAJA(OPFLOW opflow, Vec X, Vec Ge)
{
  PetscErrorCode ierr;
//  PBPOLRAJA         pbpolraja=(PBPOLRAJA)opflow->model;
  PbpolModelRaja* pbpolraja = reinterpret_cast<PbpolModelRaja*>(opflow->model);
  BUSParams      *busparams=&pbpolraja->busparams;
  GENParams      *genparams=&pbpolraja->genparams;
  LOADParams     *loadparams=&pbpolraja->loadparams;
  LINEParams     *lineparams=&pbpolraja->lineparams;
  PetscScalar    *x,*ge;

  PetscFunctionBegin;
  ierr = VecSet(Ge,0.0);CHKERRQ(ierr);
  ierr = VecGetArray(X,&x);CHKERRQ(ierr);
  ierr = VecGetArray(Ge,&ge);CHKERRQ(ierr);

  // Get sizes of resources to allocate
  int ge_sz; ierr = VecGetSize(Ge, &ge_sz);
  int x_sz; ierr = VecGetSize(X, &x_sz);

  /// TODO: handle this another way such that we do not have to dynamically
  /// create host/dev allocators all the time
  auto& resmgr = umpire::ResourceManager::getInstance();
  auto hostallocator = resmgr.getAllocator("HOST");
  auto devallocator = resmgr.getAllocator("DEVICE");

  double* d_ge = paramAlloc<double>(devallocator, ge_sz);
  double* d_x = paramAlloc<double>(devallocator, x_sz);

  // umpire's memset currently only handles bytes, so must pass in an int
  // In the future, they plan to take any T
  resmgr.memset(d_ge, 0);

  registerWith(ge, ge_sz, resmgr, hostallocator);
  registerWith(x, x_sz, resmgr, hostallocator);
  resmgr.copy(d_x, x);

  /* Equality constraints */
  /* Generator contributions */
  int* gp_gidx = genparams->gidx_dev_;
  int* gp_xidx = genparams->xidx_dev_;
  RAJA::forall< RAJA::cuda_exec<128> >( RAJA::RangeSegment(0, genparams->ngenON),
  [=] __device__ (RAJA::Index_type i) {
    RAJA::atomicSub<RAJA::cuda_atomic>(&d_ge[gp_gidx[i]], d_x[gp_xidx[i]]);
    RAJA::atomicSub<RAJA::cuda_atomic>(&d_ge[gp_gidx[i]+1], d_x[gp_xidx[i]+1]);
  });

  double* pl = loadparams->pl_dev_;
  double* ql = loadparams->ql_dev_;
  int* lp_gidx = loadparams->gidx_dev_;
  RAJA::forall< RAJA::cuda_exec<128> >( RAJA::RangeSegment(0, loadparams->nload),
    [=] __device__ (RAJA::Index_type i) {
      RAJA::atomicAdd<RAJA::cuda_atomic>(&d_ge[lp_gidx[i]], pl[i]);
      RAJA::atomicAdd<RAJA::cuda_atomic>(&d_ge[lp_gidx[i]+1], ql[i]);
    });

  int* isisolated = busparams->isisolated_dev_;
  int* ispvpq = busparams->ispvpq_dev_;
  double* bl = busparams->bl_dev_;
  double* vm = busparams->vm_dev_;
  double* va = busparams->va_dev_;
  double* gl = busparams->gl_dev_;
  int* bp_xidx = busparams->xidx_dev_;
  int* bp_gidx = busparams->gidx_dev_;
  RAJA::forall<RAJA::cuda_exec<128>>(RAJA::RangeSegment(0, busparams->nbus),
    [=] __device__ (RAJA::Index_type i)
    {
      double theta= d_x[bp_xidx[i]];
      double Vm   = d_x[bp_xidx[i]+1];

      RAJA::atomicAdd<RAJA::cuda_atomic>(
        &d_ge[bp_gidx[i]],
        isisolated[i]*(theta - va[i]*PETSC_PI/180.0) + ispvpq[i]*Vm*Vm*gl[i]);

      RAJA::atomicAdd<RAJA::cuda_atomic>(
        &d_ge[bp_gidx[i]+1],
        isisolated[i]*(Vm- vm[i]) - ispvpq[i]*Vm*Vm*bl[i]);
    });

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
  int* geqidxf = lineparams->geqidxf;
  int* geqidxt = lineparams->geqidxt;

  RAJA::forall<RAJA::cuda_exec<128>>(RAJA::RangeSegment(0, lineparams->nlineON),
    [=] __device__ (RAJA::Index_type i)
    {
      double Pf,Qf,Pt,Qt;
      double thetaf=d_x[xidxf[i]], Vmf=d_x[xidxf[i]+1];
      double thetat=d_x[xidxt[i]], Vmt=d_x[xidxt[i]+1];
      double thetaft=thetaf-thetat;
      double thetatf=thetat-thetaf;
      
      Pf = Gff[i]*Vmf*Vmf  + Vmf*Vmt*(Gft[i]*cos(thetaft) + Bft[i]*sin(thetaft));
      Qf = -Bff[i]*Vmf*Vmf + Vmf*Vmt*(-Bft[i]*cos(thetaft) + Gft[i]*sin(thetaft));
      Pt = Gtt[i]*Vmt*Vmt  + Vmt*Vmf*(Gtf[i]*cos(thetatf) + Btf[i]*sin(thetatf));
      Qt = -Btt[i]*Vmt*Vmt + Vmt*Vmf*(-Btf[i]*cos(thetatf) + Gtf[i]*sin(thetatf));

      RAJA::atomicAdd<RAJA::cuda_atomic>(&d_ge[geqidxf[i]]  , Pf);
      RAJA::atomicAdd<RAJA::cuda_atomic>(&d_ge[geqidxf[i]+1], Qf);
      RAJA::atomicAdd<RAJA::cuda_atomic>(&d_ge[geqidxt[i]]  , Pt);
      RAJA::atomicAdd<RAJA::cuda_atomic>(&d_ge[geqidxt[i]+1], Qt);
    });

  resmgr.copy(ge, d_ge);
  ierr = VecRestoreArray(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(Ge,&ge);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/** 
 * @brief INEQUALITY CONSTRAINTS 
 * 
 * @param[in/out] opflow - handle for the opflow class
 * @param[in]     X      - solution vector
 * @param[out]    Gi     - inequality constraints vector
 */
PetscErrorCode OPFLOWComputeInequalityConstraints_PBPOLRAJA(OPFLOW opflow, Vec X, Vec Gi)
{
  PetscErrorCode ierr;
//  PBPOLRAJA         pbpolraja=(PBPOLRAJA)opflow->model;
  PbpolModelRaja* pbpolraja = reinterpret_cast<PbpolModelRaja*>(opflow->model);
  LINEParams     *lineparams=&pbpolraja->lineparams;
  PetscScalar    *x,*gi;

  PetscFunctionBegin;
  ierr = VecSet(Gi,0.0);CHKERRQ(ierr);
  ierr = VecGetArray(X,&x);CHKERRQ(ierr);
  ierr = VecGetArray(Gi,&gi);CHKERRQ(ierr);

  // sizes for gl == gu
  int gi_sz; VecGetSize(Gi, &gi_sz);
  int x_sz; VecGetSize(X, &x_sz);
  auto& resmgr = umpire::ResourceManager::getInstance();
  auto hostallocator = resmgr.getAllocator("HOST");
  auto devallocator = resmgr.getAllocator("DEVICE");

  double* d_gi = paramAlloc<double>(devallocator, gi_sz);
  double* d_x =  paramAlloc<double>(devallocator, x_sz);

  registerWith(x, x_sz, resmgr, hostallocator);
  registerWith(gi, gi_sz, resmgr, hostallocator);
  resmgr.copy(d_x, x);

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
  RAJA::forall<RAJA::cuda_exec<128>>(RAJA::RangeSegment(0, lineparams->nlinelim),
    [=] __device__ (RAJA::Index_type i)
    {
      int    j=linelimidx[i];
      double Pf,Qf,Pt,Qt,Sf2,St2;
      double thetaf=d_x[xidxf[j]], Vmf=d_x[xidxf[j]+1];
      double thetat=d_x[xidxt[j]], Vmt=d_x[xidxt[j]+1];
      double thetaft=thetaf-thetat;
      double thetatf=thetat-thetaf;
      
      Pf = Gff[j]*Vmf*Vmf  + Vmf*Vmt*(Gft[j]*cos(thetaft) + Bft[j]*sin(thetaft));
      Qf = -Bff[j]*Vmf*Vmf + Vmf*Vmt*(-Bft[j]*cos(thetaft) + Gft[j]*sin(thetaft));
      Pt = Gtt[j]*Vmt*Vmt  + Vmt*Vmf*(Gtf[j]*cos(thetatf) + Btf[j]*sin(thetatf));
      Qt = -Btt[j]*Vmt*Vmt + Vmt*Vmf*(-Btf[j]*cos(thetatf) + Gtf[j]*sin(thetatf));
      
      Sf2 = Pf*Pf + Qf*Qf;
      St2 = Pt*Pt + Qt*Qt;

      /* Not atomic */
      d_gi[gineqidx[i]]   = Sf2;
      d_gi[gineqidx[i]+1] = St2;
    }
  );

  resmgr.copy(gi, d_gi);

  ierr = VecRestoreArray(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(Gi,&gi);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/** 
 * @brief OBJECTIVE FUNCTION 
 * 
 * @param[in/out] opflow - handle for the opflow class
 * @param[in]     X      - solution vector
 * @param[out]    obj    - value of the objective function
 */
PetscErrorCode OPFLOWComputeObjective_PBPOLRAJA(OPFLOW opflow,Vec X,double *obj)
{
  PetscErrorCode ierr;
//  PBPOLRAJA         pbpolraja=(PBPOLRAJA)opflow->model;
  PbpolModelRaja* pbpolraja = reinterpret_cast<PbpolModelRaja*>(opflow->model);
  GENParams     *genparams=&pbpolraja->genparams;
  PS             ps=opflow->ps;
  int            isobj_gencost=opflow->obj_gencost;
  double         MVAbase=ps->MVAbase;

  PetscFunctionBegin;

  // Access raw data from PETSc vector
  PetscScalar *x; 
  ierr = VecGetArray(X,&x); CHKERRQ(ierr);
  PetscInt Nx; ///< Size of Global Vector X
  ierr = VecGetSize(X, &Nx); CHKERRQ(ierr);

  // Allocate Umpire arrays on the host and the device
  auto& resmgr = umpire::ResourceManager::getInstance();
  umpire::Allocator hostallocator = resmgr.getAllocator("HOST");
  umpire::Allocator devallocator = resmgr.getAllocator("DEVICE");

  double* d_x = paramAlloc<double>(devallocator, Nx);
  registerWith(x, Nx, resmgr, hostallocator);
  resmgr.copy(d_x, x);

  // You need local copies of device pointers so that lambda can capture them
  double* d_cost_alpha = genparams->cost_alpha_dev_;
  double* d_cost_beta  = genparams->cost_beta_dev_;
  double* d_cost_gamma = genparams->cost_gamma_dev_;
  int* d_xidx = genparams->xidx_dev_;
  
  // Set up reduce sum object
  RAJA::ReduceSum< RAJA::cuda_reduce, double> obj_val_sum(0.0);
  // Compute reduction on CUDA device
  RAJA::forall< RAJA::cuda_exec<128> >( RAJA::RangeSegment(0, genparams->ngenON), 
    [=] __device__ (RAJA::Index_type i)
    {
      double Pg = d_x[d_xidx[i]] * MVAbase; 
      obj_val_sum += isobj_gencost*(d_cost_alpha[i]*Pg*Pg + 
                                    d_cost_beta[i]*Pg + 
                                    d_cost_gamma[i]); 
    });

  // Return objective value
  *obj = static_cast<double>(obj_val_sum.get());

  // Delete temporary copies of X on host and device
  devallocator.deallocate(d_x);

  ierr = VecRestoreArray(X,&x);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/** 
 * @brief GRADIENT 
 * 
 * @param[in/out] opflow - handle for the opflow class
 * @param[in]     X      - solution vector
 * @param[out]    Grad   - gradient vector
 */
PetscErrorCode OPFLOWComputeGradient_PBPOLRAJA(OPFLOW opflow,Vec X,Vec Grad)
{
  PetscErrorCode ierr;
//  PBPOLRAJA         pbpolraja=(PBPOLRAJA)opflow->model;
  PbpolModelRaja* pbpolraja = reinterpret_cast<PbpolModelRaja*>(opflow->model);
  GENParams     *genparams=&pbpolraja->genparams;
  PS             ps=opflow->ps;
  int            isobj_gencost=opflow->obj_gencost;
  double         MVAbase=ps->MVAbase;

  PetscFunctionBegin;

  // Access PETSc vectors' raw data
  PetscScalar *x;
  PetscScalar *grad;
  ierr = VecGetArray(X,&x);CHKERRQ(ierr);
  ierr = VecGetArray(Grad,&grad);CHKERRQ(ierr);

  // Get vector size  
  PetscInt Nx, Ngrad; // Size of Global Vector X; Global Vector Y is also of same size 
  ierr = VecGetSize(X, &Nx); CHKERRQ(ierr);
  ierr = VecGetSize(Grad, &Ngrad); CHKERRQ(ierr);
  assert(Nx == Ngrad);

  // Create allocators
  auto& resmgr = umpire::ResourceManager::getInstance();
  umpire::Allocator hostallocator = resmgr.getAllocator("HOST");
  umpire::Allocator devallocator = resmgr.getAllocator("DEVICE");

  // allocate vector data on device
  double* d_x     = paramAlloc<double>(devallocator, Nx);
  double* d_grad  = paramAlloc<double>(devallocator, Ngrad);

  registerWith(x, Nx, resmgr, hostallocator);
  registerWith(grad, Nx, resmgr, hostallocator);
  resmgr.copy(d_x, x);
  resmgr.copy(d_grad, grad);

  double* d_cost_alpha = genparams->cost_alpha_dev_;
  double* d_cost_beta  = genparams->cost_beta_dev_;
  int* d_xidx = genparams->xidx_dev_;
  RAJA::forall< RAJA::cuda_exec<128> >( RAJA::RangeSegment(0, genparams->ngenON), 
    [=] __device__ (RAJA::Index_type i)
    {
      double Pg = d_x[d_xidx[i]] * MVAbase;
      d_grad[d_xidx[i]] = isobj_gencost*MVAbase*(2.0*d_cost_alpha[i]*Pg + d_cost_beta[i]);
    });
  
  resmgr.copy(grad, d_grad);

  // Delete temporary device arrays
  devallocator.deallocate(d_x);
  devallocator.deallocate(d_grad);

  ierr = VecRestoreArray(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(Grad,&grad);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/** VARIABLE BOUNDS **/
PetscErrorCode OPFLOWSetVariableBounds_PBPOLRAJA(OPFLOW opflow,Vec Xl,Vec Xu)
{
  PetscErrorCode ierr;
//  PBPOLRAJA         pbpolraja=(PBPOLRAJA)opflow->model;
  PbpolModelRaja* pbpolraja = reinterpret_cast<PbpolModelRaja*>(opflow->model);
  BUSParams      *busparams=&pbpolraja->busparams;
  GENParams      *genparams=&pbpolraja->genparams;
  PetscScalar    *xl,*xu;

  PetscFunctionBegin;
  ierr = VecGetArray(Xl,&xl);CHKERRQ(ierr);
  ierr = VecGetArray(Xu,&xu);CHKERRQ(ierr);

  /* Bounds for bus voltages

    TODO:
      Port these two critical loops

    - create raja index set to handle indirection 
    - iterate over segment and launch
    - maximize serial memory access
    - remapping to make memory access more uniform
    - concat xl/xu into groups of 4 
      - but these vectors come from hiop...
      - Will we get these same vectors repeatedly?
      - how do you map the network/grid into memory space to minimize device traffic?
    - manual cuda kernel, profile, evaluate occupancy benchmark current impl with creating local 
      pointer to offset instead of repeating
      - gather/scatter operations
    - applies at runtime, could enhance/degredate performance depending on variables at runtime
    - partition and reorder network to make memory acces more uniform
    - **graph coloring problem**
  */

  PetscInt sz;
  ierr = VecGetSize(Xl, &sz); CHKERRQ(ierr);

  auto& resmgr = umpire::ResourceManager::getInstance();
  umpire::Allocator hostallocator = resmgr.getAllocator("HOST");
  umpire::Allocator devallocator = resmgr.getAllocator("DEVICE");

  double* d_xl  = paramAlloc<double>(devallocator, sz);
  double* d_xu  = paramAlloc<double>(devallocator, sz);
  registerWith(xl, sz, resmgr, hostallocator);
  registerWith(xu, sz, resmgr, hostallocator);

  /*
  indices may not be regularly strided
  reduce registers used in each SM since we have to dereference fewer addresses
  */
  int* d_xidx = busparams->xidx_dev_;
  int* bp_ispvpq = busparams->ispvpq_dev_;
  int* is_ref = busparams->isref_dev_;
  int* is_isolated = busparams->isisolated_dev_;
  double* va = busparams->va_dev_;
  double* vm = busparams->vm_dev_;
  double* vmin = busparams->vmin_dev_;
  double* vmax = busparams->vmax_dev_;
  RAJA::forall<RAJA::cuda_exec<128>>(RAJA::RangeSegment(0, busparams->nbus)/* index set here */,
    [=] __device__ (RAJA::Index_type i)
    {
      d_xl[d_xidx[i]] = bp_ispvpq[i]*PETSC_NINFINITY + is_isolated[i]*va[i] + is_ref[i]*va[i]*PETSC_PI/180.0;
      d_xu[d_xidx[i]] = bp_ispvpq[i]*PETSC_INFINITY  + is_isolated[i]*va[i] + is_ref[i]*va[i]*PETSC_PI/180.0;
      d_xl[d_xidx[i]+1] = is_ref[i]*vmin[i]  + bp_ispvpq[i]*vmin[i] + is_isolated[i]*vm[i];
      d_xu[d_xidx[i]+1] = is_ref[i]*vmax[i]  + bp_ispvpq[i]*vmax[i] + is_isolated[i]*vm[i];
    }
  );

  int* idx = genparams->xidx_dev_;
  double* pb = genparams->pb_dev_;
  double* pt = genparams->pt_dev_;
  double* qb = genparams->qb_dev_;
  double* qt = genparams->qt_dev_;
  RAJA::forall< RAJA::cuda_exec<128> >(RAJA::RangeSegment(0, genparams->ngenON),
    [=] __device__ (RAJA::Index_type i)
    {
      d_xl[idx[i]]   = pb[i];
      d_xu[idx[i]]   = pt[i];
      d_xl[idx[i]+1] = qb[i];
      d_xu[idx[i]+1] = qt[i];
    }
  );

  // copy back from device to petsc vec
  // assumption is that whatever we get from hiop is already on the device
  resmgr.copy(xl, d_xl);
  resmgr.copy(xu, d_xu);

  ierr = VecRestoreArray(Xl,&xl);CHKERRQ(ierr);
  ierr = VecRestoreArray(Xu,&xu);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeEqualityConstraintJacobian_PBPOLRAJA(OPFLOW opflow,Vec X,Mat Je)
{
  PetscErrorCode ierr;
  PetscInt       i,row[2],col[4];
  PbpolModelRaja* pbpolraja = reinterpret_cast<PbpolModelRaja*>(opflow->model);
  BUSParams      *busparams=&pbpolraja->busparams;
  GENParams      *genparams=&pbpolraja->genparams;
  LINEParams     *lineparams=&pbpolraja->lineparams;
  PetscScalar    val[8];
  PetscScalar    *x;

  PetscFunctionBegin;
  ierr = MatZeroEntries(Je);CHKERRQ(ierr);

  ierr = VecGetArray(X,&x);CHKERRQ(ierr);

  /* Jacobian from bus contributions */
  for(i=0; i < busparams->nbus; i++) {
    double Vm    = x[busparams->xidx[i]+1];
    row[0] = busparams->gidx[i];
    row[1] = busparams->gidx[i]+1;
    
    col[0] = busparams->xidx[i];
    col[1] = col[0]+1;

    val[0] = busparams->isisolated[i]*1.0 + busparams->ispvpq[i]*0.0;
    val[1] = busparams->isisolated[i]*0.0 + busparams->ispvpq[i]*2*Vm*busparams->gl[i];
    val[2] = 0.0;
    val[3] = busparams->isisolated[i]*1.0 + busparams->ispvpq[i]*-2*Vm*busparams->bl[i];

    ierr = MatSetValues(Je,2,row,2,col,val,ADD_VALUES);CHKERRQ(ierr);
  }

  /* Jacobian from generator contributions */
  for(i=0; i < genparams->ngenON; i++) {
    row[0] = genparams->gidx[i];
    col[0] = genparams->xidx[i];
    val[0] = -1;

    ierr = MatSetValues(Je,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);
      
    row[0] = genparams->gidx[i]+1;

    col[0] = genparams->xidx[i]+1;
    ierr = MatSetValues(Je,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);
  }

  /* Jacobian from line contributions */
  for(i=0; i < lineparams->nlineON; i++) {
    double thetaf=x[lineparams->xidxf[i]], Vmf=x[lineparams->xidxf[i]+1];
    double thetat=x[lineparams->xidxt[i]], Vmt=x[lineparams->xidxt[i]+1];
    double thetaft=thetaf-thetat;
    double thetatf=thetat-thetaf;

    row[0] = lineparams->geqidxf[i];
    row[1] = lineparams->geqidxf[i]+1;

    col[0] = lineparams->xidxf[i]; 
    col[1] = lineparams->xidxf[i]+1; 
    col[2] = lineparams->xidxt[i];
    col[3] = lineparams->xidxt[i]+1;

    /* dPf_dthetaf */
    val[0] = Vmf*Vmt*(-lineparams->Gft[i]*sin(thetaft) + lineparams->Bft[i]*cos(thetaft));
    /*dPf_dVmf */
    val[1] = 2*lineparams->Gff[i]*Vmf + Vmt*(lineparams->Gft[i]*cos(thetaft) + lineparams->Bft[i]*sin(thetaft));
    /*dPf_dthetat */
    val[2] = Vmf*Vmt*(lineparams->Gft[i]*sin(thetaft) - lineparams->Bft[i]*cos(thetaft));
    /* dPf_dVmt */
    val[3] = Vmf*(lineparams->Gft[i]*cos(thetaft) + lineparams->Bft[i]*sin(thetaft));

    /* dQf_dthetaf */
    val[4] = Vmf*Vmt*(lineparams->Bft[i]*sin(thetaft) + lineparams->Gft[i]*cos(thetaft));
    /* dQf_dVmf */
    val[5] = -2*lineparams->Bff[i]*Vmf + Vmt*(-lineparams->Bft[i]*cos(thetaft) + lineparams->Gft[i]*sin(thetaft));
    /* dQf_dthetat */
    val[6] = Vmf*Vmt*(-lineparams->Bft[i]*sin(thetaft) - lineparams->Gft[i]*cos(thetaft));
    /* dQf_dVmt */
    val[7] = Vmf*(-lineparams->Bft[i]*cos(thetaft) + lineparams->Gft[i]*sin(thetaft));

    ierr = MatSetValues(Je,2,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);

    row[0] = lineparams->geqidxt[i];
    row[1] = lineparams->geqidxt[i]+1;

    col[0] = lineparams->xidxt[i]; 
    col[1] = lineparams->xidxt[i]+1;
    col[2] = lineparams->xidxf[i];
    col[3] = lineparams->xidxf[i]+1;

    /* dPt_dthetat */
    val[0] = Vmt*Vmf*(-lineparams->Gtf[i]*sin(thetatf) + lineparams->Btf[i]*cos(thetatf));
    /* dPt_dVmt */
    val[1] = 2*lineparams->Gtt[i]*Vmt + Vmf*(lineparams->Gtf[i]*cos(thetatf) + lineparams->Btf[i]*sin(thetatf));
    /* dPt_dthetaf */
    val[2] = Vmt*Vmf*(lineparams->Gtf[i]*sin(thetatf) - lineparams->Btf[i]*cos(thetatf));
    /* dPt_dVmf */
    val[3] = Vmt*(lineparams->Gtf[i]*cos(thetatf) + lineparams->Btf[i]*sin(thetatf));
    
    /* dQt_dthetat */
    val[4] = Vmt*Vmf*(lineparams->Btf[i]*sin(thetatf) + lineparams->Gtf[i]*cos(thetatf));
    /* dQt_dVmt */
    val[5] = -2*lineparams->Btt[i]*Vmt + Vmf*(-lineparams->Btf[i]*cos(thetatf) + lineparams->Gtf[i]*sin(thetatf));
    /* dQt_dthetaf */
    val[6] = Vmt*Vmf*(-lineparams->Btf[i]*sin(thetatf) - lineparams->Gtf[i]*cos(thetatf));
    /* dQt_dVmf */
    val[7] = Vmt*(-lineparams->Btf[i]*cos(thetatf) + lineparams->Gtf[i]*sin(thetatf));
    ierr = MatSetValues(Je,2,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);
  }

  ierr = VecRestoreArray(X,&x);CHKERRQ(ierr);

  ierr = MatAssemblyBegin(Je,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Je,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeInequalityConstraintJacobian_PBPOLRAJA(OPFLOW opflow,Vec X,Mat Ji)
{
  PetscErrorCode ierr;
  PetscInt       i;
//  PBPOLRAJA         pbpolraja=(PBPOLRAJA)opflow->model;
  PbpolModelRaja* pbpolraja = reinterpret_cast<PbpolModelRaja*>(opflow->model);
  LINEParams     *lineparams=&pbpolraja->lineparams;
  PetscInt       row[2],col[4];
  PetscScalar    val[4];
  PetscScalar    *x;

  PetscFunctionBegin;
  ierr = MatZeroEntries(Ji);CHKERRQ(ierr);

  ierr = VecGetArray(X,&x);CHKERRQ(ierr);

  for(i=0; i < lineparams->nlinelim; i++) {
    int j = lineparams->linelimidx[i];

    double Pf,Qf,Pt,Qt;
    double thetaf=x[lineparams->xidxf[j]], Vmf=x[lineparams->xidxf[j]+1];
    double thetat=x[lineparams->xidxt[j]], Vmt=x[lineparams->xidxt[j]+1];
    double thetaft=thetaf-thetat;
    double thetatf=thetat-thetaf;
    double dSf2_dPf, dSf2_dQf, dSt2_dPt, dSt2_dQt;
    double dPf_dthetaf,dPf_dVmf,dPf_dthetat,dPf_dVmt;
    double dQf_dthetaf,dQf_dVmf,dQf_dthetat,dQf_dVmt;
    double dPt_dthetaf,dPt_dVmf,dPt_dthetat,dPt_dVmt;
    double dQt_dthetaf,dQt_dVmf,dQt_dthetat,dQt_dVmt;
    double dSf2_dthetaf,dSf2_dVmf,dSf2_dthetat,dSf2_dVmt;
    double dSt2_dthetaf,dSt2_dVmf,dSt2_dthetat,dSt2_dVmt;
    double Gff = lineparams->Gff[j], Bff = lineparams->Bff[j];
    double Gft = lineparams->Gft[j], Bft = lineparams->Bft[j];
    double Gtf = lineparams->Gtf[j], Btf = lineparams->Btf[j];
    double Gtt = lineparams->Gtt[j], Btt = lineparams->Btt[j];

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
      
    row[0] = lineparams->gineqidx[i];
    col[0] = lineparams->xidxf[j];
    col[1] = lineparams->xidxf[j]+1;
    col[2] = lineparams->xidxt[j];
    col[3] = lineparams->xidxt[j]+1;

    val[0] = dSf2_dthetaf;
    val[1] = dSf2_dVmf;
    val[2] = dSf2_dthetat;
    val[3] = dSf2_dVmt;
    ierr = MatSetValues(Ji,1,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);
      
    dSt2_dthetaf = dSt2_dPt*dPt_dthetaf + dSt2_dQt*dQt_dthetaf;
    dSt2_dthetat = dSt2_dPt*dPt_dthetat + dSt2_dQt*dQt_dthetat;
    dSt2_dVmf    = dSt2_dPt*dPt_dVmf    + dSt2_dQt*dQt_dVmf;
    dSt2_dVmt    = dSt2_dPt*dPt_dVmt    + dSt2_dQt*dQt_dVmt;
      
    row[0] = lineparams->gineqidx[i]+1;
    col[0] = lineparams->xidxt[j];
    col[1] = lineparams->xidxt[j]+1;
    col[2] = lineparams->xidxf[j];
    col[3] = lineparams->xidxf[j]+1;

    val[0] = dSt2_dthetat;
    val[1] = dSt2_dVmt;
    val[2] = dSt2_dthetaf;
    val[3] = dSt2_dVmf;
    ierr = MatSetValues(Ji,1,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);
      
  }
  ierr = VecRestoreArray(X,&x);CHKERRQ(ierr);

  ierr = MatAssemblyBegin(Ji,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Ji,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  OPFLOWComputeObjectiveHessian - Computes the Hessian for the objective function part
  
  Input Parameters:
+ opflow - the OPFLOW object
- X        - solution vecto X

  Output Parameters:
. H - the Hessian part for the objective function

*/
PetscErrorCode OPFLOWComputeObjectiveHessian_PBPOLRAJA(OPFLOW opflow,Vec X,Mat H) 
{
  PetscErrorCode ierr;
//  PBPOLRAJA         pbpolraja=(PBPOLRAJA)opflow->model;
  PbpolModelRaja* pbpolraja = reinterpret_cast<PbpolModelRaja*>(opflow->model);
  GENParams     *genparams=&pbpolraja->genparams;
  PS             ps=opflow->ps;
  PetscInt       i;
  PetscScalar    *x;
  PetscInt       row[2],col[2];
  PetscScalar    val[2];
  double         obj_factor = opflow->obj_factor;
  int            isobj_gencost=opflow->obj_gencost;
  double         MVAbase=ps->MVAbase;

  PetscFunctionBegin;

  ierr = VecGetArray(X,&x);CHKERRQ(ierr);

  /* Generator Objective function Hessian contributions */
  for(i=0; i < genparams->ngenON; i++) {
    row[0] = genparams->xidx[i];
    col[0] = genparams->xidx[i];
    val[0] = isobj_gencost*obj_factor*2.0*genparams->cost_alpha[i]*MVAbase*MVAbase;
    ierr = MatSetValues(H,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);
  }
    
  ierr = VecRestoreArray(X,&x);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


/*
  OPFLOWComputeEqualityConstraintsHessian - Computes the Hessian for the equality constraints function part
  
  Input Parameters:
+ opflow   - the OPFLOW object
. X        - solution vector X
- Lambda   - Lagrangian multiplier vector

  Output Parameters:
. H - the Hessian part for the equality constraints

*/
PetscErrorCode OPFLOWComputeEqualityConstraintsHessian_PBPOLRAJA(OPFLOW opflow,Vec X,Vec Lambda,Mat H) 
{
  PetscErrorCode ierr;
//  PBPOLRAJA         pbpolraja=(PBPOLRAJA)opflow->model;
  PbpolModelRaja* pbpolraja = reinterpret_cast<PbpolModelRaja*>(opflow->model);
  BUSParams      *busparams=&pbpolraja->busparams;
  LINEParams     *lineparams=&pbpolraja->lineparams;
  PetscInt       i;
  PetscInt       row[16],col[16];
  PetscScalar    val[16];
  PetscScalar    *x,*lambda;

  PetscFunctionBegin;

  ierr = VecGetArray(X,&x);CHKERRQ(ierr);
  ierr = VecGetArray(Lambda,&lambda);CHKERRQ(ierr);

  /* Hessian from bus contributions */
  for(i=0; i < busparams->nbus; i++) {
    row[0] = busparams->xidx[i]+1;
    col[0] = row[0];
    val[0] = busparams->ispvpq[i]*(lambda[busparams->gidx[i]]*2*busparams->gl[i] + lambda[busparams->gidx[i]+1]*(-2*busparams->bl[i]));
    ierr = MatSetValues(H,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);
  }

  /* Hessian from line contributions */
  for(i=0; i < lineparams->nlineON; i++) {
    int    gloc;
    double Gff,Bff,Gft,Bft,Gtf,Btf,Gtt,Btt;
    Gff = lineparams->Gff[i];
    Bff = lineparams->Bff[i];
    Gft = lineparams->Gft[i];
    Bft = lineparams->Bft[i];
    Gtf = lineparams->Gtf[i];
    Btf = lineparams->Btf[i];
    Gtt = lineparams->Gtt[i];
    Btt = lineparams->Btt[i];
    
    double thetaf=x[lineparams->xidxf[i]], Vmf=x[lineparams->xidxf[i]+1];
    double thetat=x[lineparams->xidxt[i]], Vmt=x[lineparams->xidxt[i]+1];
    double thetaft=thetaf-thetat;
    double thetatf=thetat-thetaf;
    
    double dPf_dthetaf_dthetaf,dPf_dthetaf_dVmf,dPf_dthetaf_dthetat,dPf_dthetaf_dVmt;
    double dPf_dVmf_dthetaf,   dPf_dVmf_dVmf,   dPf_dVmf_dthetat,   dPf_dVmf_dVmt;
    double dPf_dthetat_dthetaf,dPf_dthetat_dVmf,dPf_dthetat_dthetat,dPf_dthetat_dVmt;
    double dPf_dVmt_dthetaf,   dPf_dVmt_dVmf,   dPf_dVmt_dthetat,   dPf_dVmt_dVmt;
    
    /* dPf_dthetaf = Vmf*Vmt*(-Gft*sin(thetaft) + Bft*cos(thetaft)); */
    dPf_dthetaf_dthetaf = -Vmf*Vmt*(Gft*PetscCosScalar(thetaft) + Bft*PetscSinScalar(thetaft));
    dPf_dthetaf_dVmf    =     Vmt*(-Gft*PetscSinScalar(thetaft) + Bft*PetscCosScalar(thetaft));
    dPf_dthetaf_dthetat =  Vmf*Vmt*(Gft*PetscCosScalar(thetaft) + Bft*PetscSinScalar(thetaft));
    dPf_dthetaf_dVmt    =     Vmf*(-Gft*PetscSinScalar(thetaft) + Bft*PetscCosScalar(thetaft));
    
    /* dPf_Vmf  = 2*Gff*Vmf + Vmt*(Gft*cos(thetaft) + Bft*sin(thetaft)); */
    dPf_dVmf_dthetaf    =  Vmt*(-Gft*PetscSinScalar(thetaft) + Bft*PetscCosScalar(thetaft)); 
    dPf_dVmf_dVmf       =  2*Gff;
    dPf_dVmf_dthetat    =  Vmt*(Gft*PetscSinScalar(thetaft) - Bft*PetscCosScalar(thetaft));
    dPf_dVmf_dVmt       =      (Gft*PetscCosScalar(thetaft) + Bft*PetscSinScalar(thetaft));
    
    /* dPf_dthetat = Vmf*Vmt*(Gft*sin(thetaft) - Bft*cos(thetaft)); */
    dPf_dthetat_dthetaf = Vmf*Vmt*(Gft*PetscCosScalar(thetaft)  + Bft*PetscSinScalar(thetaft));
    dPf_dthetat_dVmf    =     Vmt*(Gft*PetscSinScalar(thetaft)  - Bft*PetscCosScalar(thetaft));
    dPf_dthetat_dthetat = Vmf*Vmt*(-Gft*PetscCosScalar(thetaft) - Bft*PetscSinScalar(thetaft));
    dPf_dthetat_dVmt    =     Vmf*(Gft*PetscSinScalar(thetaft)  - Bft*PetscCosScalar(thetaft));
    
    /* dPf_dVmt = Vmf*(Gft*cos(thetaft) + Bft*sin(thetaft)); */
    dPf_dVmt_dthetaf    = Vmf*(-Gft*PetscSinScalar(thetaft) + Bft*PetscCosScalar(thetaft)); 
    dPf_dVmt_dVmf       =      (Gft*PetscCosScalar(thetaft) + Bft*PetscSinScalar(thetaft));
    dPf_dVmt_dthetat    = Vmf*(Gft*PetscSinScalar(thetaft) - Bft*PetscCosScalar(thetaft));
    dPf_dVmt_dVmt       = 0.0;
    
    double dQf_dthetaf_dthetaf,dQf_dthetaf_dVmf,dQf_dthetaf_dthetat,dQf_dthetaf_dVmt;
    double dQf_dVmf_dthetaf,   dQf_dVmf_dVmf,   dQf_dVmf_dthetat,   dQf_dVmf_dVmt;
    double dQf_dthetat_dthetaf,dQf_dthetat_dVmf,dQf_dthetat_dthetat,dQf_dthetat_dVmt;
    double dQf_dVmt_dthetaf,   dQf_dVmt_dVmf,   dQf_dVmt_dthetat,   dQf_dVmt_dVmt;
    
    /* dQf_dthetaf = Vmf*Vmt*(Bft*sin(thetaft) + Gft*cos(thetaft)); */
    dQf_dthetaf_dthetaf = Vmf*Vmt*(Bft*PetscCosScalar(thetaft)  - Gft*PetscSinScalar(thetaft));
    dQf_dthetaf_dVmf    =     Vmt*(Bft*PetscSinScalar(thetaft)  + Gft*PetscCosScalar(thetaft));
    dQf_dthetaf_dthetat = Vmf*Vmt*(-Bft*PetscCosScalar(thetaft) + Gft*PetscSinScalar(thetaft));
    dQf_dthetaf_dVmt    =     Vmf*(Bft*PetscSinScalar(thetaft) + Gft*PetscCosScalar(thetaft));
    
    /* dQf_dVmf = -2*Bff*Vmf + Vmt*(-Bft*cos(thetaft) + Gft*sin(thetaft)); */
    dQf_dVmf_dthetaf    =  Vmt*(Bft*PetscSinScalar(thetaft) + Gft*PetscCosScalar(thetaft)); 
    dQf_dVmf_dVmf       = -2*Bff;
    dQf_dVmf_dthetat    =  Vmt*(-Bft*PetscSinScalar(thetaft) - Gft*PetscCosScalar(thetaft));
    dQf_dVmf_dVmt       =      (-Bft*PetscCosScalar(thetaft) + Gft*PetscSinScalar(thetaft));
    
    /* dQf_dthetat = Vmf*Vmt*(-Bft*sin(thetaft) - Gft*cos(thetaft)); */
    dQf_dthetat_dthetaf = Vmf*Vmt*(-Bft*PetscCosScalar(thetaft) + Gft*PetscSinScalar(thetaft));
    dQf_dthetat_dVmf    =     Vmt*(-Bft*PetscSinScalar(thetaft) - Gft*PetscCosScalar(thetaft));
    dQf_dthetat_dthetat = Vmf*Vmt*( Bft*PetscCosScalar(thetaft) - Gft*PetscSinScalar(thetaft));
    dQf_dthetat_dVmt    =     Vmf*(-Bft*PetscSinScalar(thetaft) - Gft*PetscCosScalar(thetaft));
    
    /* dQf_dVmt = Vmf*(-Bft*cos(thetaft) + Gft*sin(thetaft)); */
    dQf_dVmt_dthetaf    = Vmf*(Bft*PetscSinScalar(thetaft) + Gft*PetscCosScalar(thetaft)); 
    dQf_dVmt_dVmf       =    (-Bft*PetscCosScalar(thetaft) + Gft*PetscSinScalar(thetaft));
    dQf_dVmt_dthetat    = Vmf*(-Bft*PetscSinScalar(thetaft) - Gft*PetscCosScalar(thetaft));
    dQf_dVmt_dVmt       = 0.0;
    
    row[0] = lineparams->xidxf[i]; row[1] = lineparams->xidxf[i]+1;
    col[0] = lineparams->xidxf[i]; 
    col[1] = lineparams->xidxf[i]+1; 
    col[2] = lineparams->xidxt[i];
    col[3] = lineparams->xidxt[i]+1;
    
    gloc=lineparams->geqidxf[i];
    
    val[0] = val[1] = val[2] = val[3] = val[4] = val[5] = val[6] = val[7] = 0.0;
    
    val[0] = lambda[gloc]*dPf_dthetaf_dthetaf + lambda[gloc+1]*dQf_dthetaf_dthetaf;
    val[1] = lambda[gloc]*dPf_dthetaf_dVmf    + lambda[gloc+1]*dQf_dthetaf_dVmf;
    val[2] = lambda[gloc]*dPf_dthetaf_dthetat + lambda[gloc+1]*dQf_dthetaf_dthetat;
    val[3] = lambda[gloc]*dPf_dthetaf_dVmt    + lambda[gloc+1]*dQf_dthetaf_dVmt;
    
    val[4] = lambda[gloc]*dPf_dVmf_dthetaf + lambda[gloc+1]*dQf_dVmf_dthetaf;
    val[5] = lambda[gloc]*dPf_dVmf_dVmf    + lambda[gloc+1]*dQf_dVmf_dVmf;
    val[6] = lambda[gloc]*dPf_dVmf_dthetat + lambda[gloc+1]*dQf_dVmf_dthetat;
    val[7] = lambda[gloc]*dPf_dVmf_dVmt    + lambda[gloc+1]*dQf_dVmf_dVmt;
    
    ierr = MatSetValues(H,2,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);
    
    row[0] = lineparams->xidxt[i];
    row[1] = row[0]+1;
    
    col[0] = lineparams->xidxf[i];
    col[1] = col[0]+1;
    col[2] = lineparams->xidxt[i];
    col[3] = col[2]+1;
    
    val[0] = val[1] = val[2] = val[3] = val[4] = val[5] = val[6] = val[7] = 0.0;
    
    val[0] = lambda[gloc]*dPf_dthetat_dthetaf + lambda[gloc+1]*dQf_dthetat_dthetaf;
    val[1] = lambda[gloc]*dPf_dthetat_dVmf    + lambda[gloc+1]*dQf_dthetat_dVmf;
    val[2] = lambda[gloc]*dPf_dthetat_dthetat + lambda[gloc+1]*dQf_dthetat_dthetat;
    val[3] = lambda[gloc]*dPf_dthetat_dVmt    + lambda[gloc+1]*dQf_dthetat_dVmt;
    
    val[4] = lambda[gloc]*dPf_dVmt_dthetaf + lambda[gloc+1]*dQf_dVmt_dthetaf;
    val[5] = lambda[gloc]*dPf_dVmt_dVmf    + lambda[gloc+1]*dQf_dVmt_dVmf;
    val[6] = lambda[gloc]*dPf_dVmt_dthetat + lambda[gloc+1]*dQf_dVmt_dthetat;
    val[7] = lambda[gloc]*dPf_dVmt_dVmt    + lambda[gloc+1]*dQf_dVmt_dVmt;
    
    ierr = MatSetValues(H,2,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);
	
    double dPt_dthetat_dthetat,dPt_dthetat_dVmt,dPt_dthetat_dthetaf,dPt_dthetat_dVmf;
    double dPt_dVmt_dthetat,   dPt_dVmt_dVmt,   dPt_dVmt_dthetaf,   dPt_dVmt_dVmf;
    double dPt_dthetaf_dthetat,dPt_dthetaf_dVmt,dPt_dthetaf_dthetaf,dPt_dthetaf_dVmf;
    double dPt_dVmf_dthetat,   dPt_dVmf_dVmt,   dPt_dVmf_dthetaf,   dPt_dVmf_dVmf;

    /* dPt_dthetat = Vmf*Vmt*(-Gtf*sin(thetatf) + Btf*cos(thetatf)); */
    dPt_dthetat_dthetat = Vmf*Vmt*(-Gtf*PetscCosScalar(thetatf) - Btf*PetscSinScalar(thetatf));
    dPt_dthetat_dVmt    =     Vmf*(-Gtf*PetscSinScalar(thetatf) + Btf*PetscCosScalar(thetatf));
    dPt_dthetat_dthetaf =  Vmf*Vmt*(Gtf*PetscCosScalar(thetatf) + Btf*PetscSinScalar(thetatf));
    dPt_dthetat_dVmf    =     Vmt*(-Gtf*PetscSinScalar(thetatf) + Btf*PetscCosScalar(thetatf));
    
    /* dPt_Vmt  = 2*Gtt*Vmt + Vmf*(Gtf*cos(thetatf) + Btf*sin(thetatf)); */
    dPt_dVmt_dthetat    =  Vmf*(-Gtf*PetscSinScalar(thetatf) + Bft*PetscCosScalar(thetatf)); 
    dPt_dVmt_dVmt       =  2*Gtt;
    dPt_dVmt_dthetaf    =  Vmf*(Gtf*PetscSinScalar(thetatf) - Btf*PetscCosScalar(thetatf));
    dPt_dVmt_dVmf       =      (Gtf*PetscCosScalar(thetatf) + Btf*PetscSinScalar(thetatf));
    
    /* dPt_dthetaf = Vmf*Vmt*(Gtf*sin(thetatf) - Btf*cos(thetatf)); */
    dPt_dthetaf_dthetat = Vmf*Vmt*(Gtf*PetscCosScalar(thetatf)  + Btf*PetscSinScalar(thetatf));
    dPt_dthetaf_dVmt    =     Vmf*(Gtf*PetscSinScalar(thetatf)  - Btf*PetscCosScalar(thetatf));
    dPt_dthetaf_dthetaf = Vmf*Vmt*(-Gtf*PetscCosScalar(thetatf) - Btf*PetscSinScalar(thetatf));
    dPt_dthetaf_dVmf    =     Vmt*(Gtf*PetscSinScalar(thetatf)  - Btf*PetscCosScalar(thetatf));
    
    /* dPt_dVmf = Vmt*(Gtf*cos(thetatf) + Btf*sin(thetatf)); */
    dPt_dVmf_dthetat    = Vmt*(-Gtf*PetscSinScalar(thetatf) + Btf*PetscCosScalar(thetatf)); 
    dPt_dVmf_dVmt       =     (Gtf*PetscCosScalar(thetatf) + Btf*PetscSinScalar(thetatf));
    dPt_dVmf_dthetaf    = Vmt*(Gtf*PetscSinScalar(thetatf) - Btf*PetscCosScalar(thetatf));
    dPt_dVmf_dVmf       = 0.0;
    
    double dQt_dthetaf_dthetaf,dQt_dthetaf_dVmf,dQt_dthetaf_dthetat,dQt_dthetaf_dVmt;
    double dQt_dVmf_dthetaf,   dQt_dVmf_dVmf,   dQt_dVmf_dthetat,   dQt_dVmf_dVmt;
    double dQt_dthetat_dthetaf,dQt_dthetat_dVmf,dQt_dthetat_dthetat,dQt_dthetat_dVmt;
    double dQt_dVmt_dthetaf,   dQt_dVmt_dVmf,   dQt_dVmt_dthetat,   dQt_dVmt_dVmt;
    
    /* dQt_dthetat = Vmf*Vmt*(Btf*sin(thetatf) + Gtf*cos(thetatf)); */
    dQt_dthetat_dthetat = Vmf*Vmt*(Btf*PetscCosScalar(thetatf)  - Gtf*PetscSinScalar(thetatf));
    dQt_dthetat_dVmt    =     Vmf*(Btf*PetscSinScalar(thetatf)  + Gtf*PetscCosScalar(thetatf));
    dQt_dthetat_dthetaf = Vmf*Vmt*(-Btf*PetscCosScalar(thetatf) + Gtf*PetscSinScalar(thetatf));
    dQt_dthetat_dVmf    =     Vmf*(Btf*PetscSinScalar(thetatf) + Gtf*PetscCosScalar(thetatf));
    
    /* dQt_dVmt = -2*Btt*Vmt + Vmf*(-Btf*cos(thetatf) + Gtf*sin(thetatf)); */
    dQt_dVmt_dthetat    =  Vmf*(Btf*PetscSinScalar(thetatf) + Gtf*PetscCosScalar(thetatf)); 
    dQt_dVmt_dVmt       = -2*Btt;
    dQt_dVmt_dthetaf    =  Vmf*(-Btf*PetscSinScalar(thetatf) + Gtf*PetscCosScalar(thetatf));
    dQt_dVmt_dVmf       =      (-Btf*PetscCosScalar(thetatf) + Gtf*PetscSinScalar(thetatf));
    
    /* dQt_dthetaf = Vmf*Vmt*(-Btf*sin(thetatf) - Gtf*cos(thetatf)); */
    dQt_dthetaf_dthetat = Vmf*Vmt*(-Btf*PetscCosScalar(thetatf) + Gtf*PetscSinScalar(thetatf));
    dQt_dthetaf_dVmt    =     Vmf*(-Btf*PetscSinScalar(thetatf) - Gtf*PetscCosScalar(thetatf));
    dQt_dthetaf_dthetaf = Vmf*Vmt*( Btf*PetscCosScalar(thetatf) - Gtf*PetscSinScalar(thetatf));
    dQt_dthetaf_dVmf    =     Vmt*(-Btf*PetscSinScalar(thetatf) - Gtf*PetscCosScalar(thetatf));
    
    /* dQt_dVmf = Vmt*(-Btf*cos(thetatf) + Gtf*sin(thetatf)); */
    dQt_dVmf_dthetat    = Vmt*(Btf*PetscSinScalar(thetatf) + Gtf*PetscCosScalar(thetatf)); 
    dQt_dVmf_dVmt       =    (-Btf*PetscCosScalar(thetatf) + Gtf*PetscSinScalar(thetatf));
    dQt_dVmf_dthetaf    = Vmt*(-Btf*PetscSinScalar(thetatf) - Gtf*PetscCosScalar(thetatf));
    dQt_dVmf_dVmf       = 0.0;
    
    row[0] = lineparams->xidxt[i];
    row[1] = row[0]+1;
    col[0] = lineparams->xidxt[i]; 
    col[1] = col[0]+1; 
    col[2] = lineparams->xidxf[i]; 
    col[3] = lineparams->xidxf[i]+1;

    gloc   = lineparams->geqidxt[i];

    val[0] = val[1] = val[2] = val[3] = val[4] = val[5] = val[6] = val[7] = 0.0;

    val[0] = lambda[gloc]*dPt_dthetat_dthetat + lambda[gloc+1]*dQt_dthetat_dthetat;
    val[1] = lambda[gloc]*dPt_dthetat_dVmt    + lambda[gloc+1]*dQt_dthetat_dVmt;
    val[2] = lambda[gloc]*dPt_dthetat_dthetaf + lambda[gloc+1]*dQt_dthetat_dthetaf;
    val[3] = lambda[gloc]*dPt_dthetat_dVmf    + lambda[gloc+1]*dQt_dthetat_dVmf;

    val[4] = lambda[gloc]*dPt_dVmt_dthetat + lambda[gloc+1]*dQt_dVmt_dthetat;
    val[5] = lambda[gloc]*dPt_dVmt_dVmt    + lambda[gloc+1]*dQt_dVmt_dVmt;
    val[6] = lambda[gloc]*dPt_dVmt_dthetaf + lambda[gloc+1]*dQt_dVmt_dthetaf;
    val[7] = lambda[gloc]*dPt_dVmt_dVmf    + lambda[gloc+1]*dQt_dVmt_dVmf;
    
    ierr = MatSetValues(H,2,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);
    
    row[0] = lineparams->xidxf[i];
    row[1] = row[0]+1;
    col[0] = lineparams->xidxt[i]; 
    col[1] = col[0]+1; 
    col[2] = lineparams->xidxf[i]; 
    col[3] = lineparams->xidxf[i]+1;
    
    val[0] = val[1] = val[2] = val[3] = val[4] = val[5] = val[6] = val[7] = 0.0;
    
    val[0] = lambda[gloc]*dPt_dthetaf_dthetat + lambda[gloc+1]*dQt_dthetaf_dthetat;
    val[1] = lambda[gloc]*dPt_dthetaf_dVmt    + lambda[gloc+1]*dQt_dthetaf_dVmt;
    val[2] = lambda[gloc]*dPt_dthetaf_dthetaf + lambda[gloc+1]*dQt_dthetaf_dthetaf;
    val[3] = lambda[gloc]*dPt_dthetaf_dVmf    + lambda[gloc+1]*dQt_dthetaf_dVmf;
    
    val[4] = lambda[gloc]*dPt_dVmf_dthetat + lambda[gloc+1]*dQt_dVmf_dthetat;
    val[5] = lambda[gloc]*dPt_dVmf_dVmt    + lambda[gloc+1]*dQt_dVmf_dVmt;
    val[6] = lambda[gloc]*dPt_dVmf_dthetaf + lambda[gloc+1]*dQt_dVmf_dthetaf;
    val[7] = lambda[gloc]*dPt_dVmf_dVmf    + lambda[gloc+1]*dQt_dVmf_dVmf;
    
    ierr = MatSetValues(H,2,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);
  }

  ierr = VecRestoreArray(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(Lambda,&lambda);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}


/*
  OPFLOWComputeInequalityConstraintsHessian - Computes the Inequality Constraints Hessian

  Input Parameters:
+ opflow   - the OPFLOW object
. X        - the solution vector
- Lambda   - Lagrangian multipler vector

  Output Parameters:
+ H   - the Hessian matrix

*/
PetscErrorCode OPFLOWComputeInequalityConstraintsHessian_PBPOLRAJA(OPFLOW opflow, Vec X, Vec Lambda,Mat H)
{
  PetscErrorCode ierr;
  PetscInt       i;
//  PBPOLRAJA         pbpolraja=(PBPOLRAJA)opflow->model;
  PbpolModelRaja* pbpolraja = reinterpret_cast<PbpolModelRaja*>(opflow->model);
  LINEParams     *lineparams=&pbpolraja->lineparams;
  PetscScalar    *x;
  PetscScalar    *lambda;
  PetscInt       row[12],col[12];
  PetscScalar    val[12];

  PetscFunctionBegin;
      
  ierr = VecGetArray(X,&x);CHKERRQ(ierr);
  ierr = VecGetArray(Lambda,&lambda);CHKERRQ(ierr);

  // Hessian from line contributions 
  for(i=0; i < lineparams->nlinelim; i++) {
    int j = lineparams->linelimidx[i];
    int gloc;
      
    double Pf,Qf,Pt,Qt;
    double thetaf=x[lineparams->xidxf[j]], Vmf=x[lineparams->xidxf[j]+1];
    double thetat=x[lineparams->xidxt[j]], Vmt=x[lineparams->xidxt[j]+1];
    double thetaft=thetaf-thetat;
    double thetatf=thetat-thetaf;
    double Gff = lineparams->Gff[j], Bff = lineparams->Bff[j];
    double Gft = lineparams->Gft[j], Bft = lineparams->Bft[j];
    double Gtf = lineparams->Gtf[j], Btf = lineparams->Btf[j];
    double Gtt = lineparams->Gtt[j], Btt = lineparams->Btt[j];

    Pf =  Gff*Vmf*Vmf + Vmf*Vmt*( Gft*PetscCosScalar(thetaft) + Bft*PetscSinScalar(thetaft));
    Qf = -Bff*Vmf*Vmf + Vmf*Vmt*(-Bft*PetscCosScalar(thetaft) + Gft*PetscSinScalar(thetaft));
    
    Pt =  Gtt*Vmt*Vmt + Vmt*Vmf*( Gtf*PetscCosScalar(thetatf) + Btf*PetscSinScalar(thetatf));
    Qt = -Btt*Vmt*Vmt + Vmt*Vmf*(-Btf*PetscCosScalar(thetatf) + Gtf*PetscSinScalar(thetatf));
      
    double dSf2_dPf, dSf2_dQf, dSt2_dPt, dSt2_dQt;
    
    dSf2_dPf = 2.*Pf;
    dSf2_dQf = 2.*Qf;
    dSt2_dPt = 2.*Pt;
    dSt2_dQt = 2.*Qt;
    
    double dPf_dthetaf,dPf_dVmf,dPf_dthetat,dPf_dVmt;
    double dQf_dthetaf,dQf_dVmf,dQf_dthetat,dQf_dVmt;
    double dPt_dthetaf,dPt_dVmf,dPt_dthetat,dPt_dVmt;
    double dQt_dthetaf,dQt_dVmf,dQt_dthetat,dQt_dVmt;
    
    dPf_dthetaf = 			Vmf*Vmt*(-Gft*PetscSinScalar(thetaft) + Bft*PetscCosScalar(thetaft));
    dPf_dVmf    = 2.*Gff*Vmf + 	Vmt*( Gft*PetscCosScalar(thetaft) + Bft*PetscSinScalar(thetaft));
    dPf_dthetat = 			Vmf*Vmt*( Gft*PetscSinScalar(thetaft) - Bft*PetscCosScalar(thetaft));
    dPf_dVmt    = 				Vmf*( Gft*PetscCosScalar(thetaft) + Bft*PetscSinScalar(thetaft));
    
    dQf_dthetaf = 			Vmf*Vmt*( Bft*PetscSinScalar(thetaft) + Gft*PetscCosScalar(thetaft));
    dQf_dVmf    = -2.*Bff*Vmf + 	Vmt*(-Bft*PetscCosScalar(thetaft) + Gft*PetscSinScalar(thetaft));
    dQf_dthetat = 			Vmf*Vmt*(-Bft*PetscSinScalar(thetaft) - Gft*PetscCosScalar(thetaft));
    dQf_dVmt    = 				Vmf*(-Bft*PetscCosScalar(thetaft) + Gft*PetscSinScalar(thetaft));
    
    dPt_dthetat = 			Vmt*Vmf*(-Gtf*PetscSinScalar(thetatf) + Btf*PetscCosScalar(thetatf));
    dPt_dVmt    = 2.*Gtt*Vmt + 	Vmf*( Gtf*PetscCosScalar(thetatf) + Btf*PetscSinScalar(thetatf));
    dPt_dthetaf = 			Vmt*Vmf*( Gtf*PetscSinScalar(thetatf) - Btf*PetscCosScalar(thetatf));
    dPt_dVmf    = 				Vmt*( Gtf*PetscCosScalar(thetatf) + Btf*PetscSinScalar(thetatf));
    
    dQt_dthetat = 			Vmt*Vmf*( Btf*PetscSinScalar(thetatf) + Gtf*PetscCosScalar(thetatf));
    dQt_dVmt    = -2.*Btt*Vmt + 	Vmf*(-Btf*PetscCosScalar(thetatf) + Gtf*PetscSinScalar(thetatf));
    dQt_dthetaf = 			Vmt*Vmf*(-Btf*PetscSinScalar(thetatf) - Gtf*PetscCosScalar(thetatf));
    dQt_dVmf    = 				Vmt*(-Btf*PetscCosScalar(thetatf) + Gtf*PetscSinScalar(thetatf));
    
    double d2Pf_dthetaf_dthetaf,d2Pf_dthetaf_dVmf,d2Pf_dthetaf_dthetat,d2Pf_dthetaf_dVmt;
    double d2Pf_dVmf_dthetaf,   d2Pf_dVmf_dVmf,   d2Pf_dVmf_dthetat,   d2Pf_dVmf_dVmt;
    double d2Pf_dthetat_dthetaf,d2Pf_dthetat_dVmf,d2Pf_dthetat_dthetat,d2Pf_dthetat_dVmt;
    double d2Pf_dVmt_dthetaf,   d2Pf_dVmt_dVmf,   d2Pf_dVmt_dthetat,   d2Pf_dVmt_dVmt;
    
    /* dPf_dthetaf = Vmf*Vmt*(-Gft*sin(thetaft) + Bft*cos(thetaft)); */
    d2Pf_dthetaf_dthetaf = -Vmf*Vmt*(Gft*PetscCosScalar(thetaft) + Bft*PetscSinScalar(thetaft));
    d2Pf_dthetaf_dVmf    =     Vmt*(-Gft*PetscSinScalar(thetaft) + Bft*PetscCosScalar(thetaft));
    d2Pf_dthetaf_dthetat =  Vmf*Vmt*(Gft*PetscCosScalar(thetaft) + Bft*PetscSinScalar(thetaft));
    d2Pf_dthetaf_dVmt    =     Vmf*(-Gft*PetscSinScalar(thetaft) + Bft*PetscCosScalar(thetaft));
    
    /* dPf_Vmf  = 2*Gff*Vmf + Vmt*(Gft*cos(thetaft) + Bft*sin(thetaft)); */
    d2Pf_dVmf_dthetaf    =  Vmt*(-Gft*PetscSinScalar(thetaft) + Bft*PetscCosScalar(thetaft)); 
    d2Pf_dVmf_dVmf       =  2*Gff;
    d2Pf_dVmf_dthetat    =  Vmt*(Gft*PetscSinScalar(thetaft) - Bft*PetscCosScalar(thetaft));
    d2Pf_dVmf_dVmt       =      (Gft*PetscCosScalar(thetaft) + Bft*PetscSinScalar(thetaft));
    
    /* dPf_dthetat = Vmf*Vmt*(Gft*sin(thetaft) - Bft*cos(thetaft)); */
    d2Pf_dthetat_dthetaf = Vmf*Vmt*(Gft*PetscCosScalar(thetaft)  + Bft*PetscSinScalar(thetaft));
    d2Pf_dthetat_dVmf    =     Vmt*(Gft*PetscSinScalar(thetaft)  - Bft*PetscCosScalar(thetaft));
    d2Pf_dthetat_dthetat = Vmf*Vmt*(-Gft*PetscCosScalar(thetaft) - Bft*PetscSinScalar(thetaft));
    d2Pf_dthetat_dVmt    =     Vmf*(Gft*PetscSinScalar(thetaft)  - Bft*PetscCosScalar(thetaft));
    
    /* dPf_dVmt = Vmf*(Gft*cos(thetaft) + Bft*sin(thetaft)); */
    d2Pf_dVmt_dthetaf    = Vmf*(-Gft*PetscSinScalar(thetaft) + Bft*PetscCosScalar(thetaft)); 
    d2Pf_dVmt_dVmf       =      (Gft*PetscCosScalar(thetaft) + Bft*PetscSinScalar(thetaft));
    d2Pf_dVmt_dthetat    = Vmf*(Gft*PetscSinScalar(thetaft) - Bft*PetscCosScalar(thetaft));
    d2Pf_dVmt_dVmt       = 0.0;
    
    double d2Qf_dthetaf_dthetaf,d2Qf_dthetaf_dVmf,d2Qf_dthetaf_dthetat,d2Qf_dthetaf_dVmt;
    double d2Qf_dVmf_dthetaf,   d2Qf_dVmf_dVmf,   d2Qf_dVmf_dthetat,   d2Qf_dVmf_dVmt;
    double d2Qf_dthetat_dthetaf,d2Qf_dthetat_dVmf,d2Qf_dthetat_dthetat,d2Qf_dthetat_dVmt;
    double d2Qf_dVmt_dthetaf,   d2Qf_dVmt_dVmf,   d2Qf_dVmt_dthetat,   d2Qf_dVmt_dVmt;
    
    /* dQf_dthetaf = Vmf*Vmt*(Bft*sin(thetaft) + Gft*cos(thetaft)); */
    d2Qf_dthetaf_dthetaf = Vmf*Vmt*(Bft*PetscCosScalar(thetaft)  - Gft*PetscSinScalar(thetaft));
    d2Qf_dthetaf_dVmf    =     Vmt*(Bft*PetscSinScalar(thetaft)  + Gft*PetscCosScalar(thetaft));
    d2Qf_dthetaf_dthetat = Vmf*Vmt*(-Bft*PetscCosScalar(thetaft) + Gft*PetscSinScalar(thetaft));
    d2Qf_dthetaf_dVmt    =     Vmf*(Bft*PetscSinScalar(thetaft) + Gft*PetscCosScalar(thetaft));
    
    /* dQf_dVmf = -2*Bff*Vmf + Vmt*(-Bft*cos(thetaft) + Gft*sin(thetaft)); */
    d2Qf_dVmf_dthetaf    =  Vmt*(Bft*PetscSinScalar(thetaft) + Gft*PetscCosScalar(thetaft)); 
    d2Qf_dVmf_dVmf       = -2*Bff;
    d2Qf_dVmf_dthetat    =  Vmt*(-Bft*PetscSinScalar(thetaft) - Gft*PetscCosScalar(thetaft));
    d2Qf_dVmf_dVmt       =      (-Bft*PetscCosScalar(thetaft) + Gft*PetscSinScalar(thetaft));
    
    /* dQf_dthetat = Vmf*Vmt*(-Bft*sin(thetaft) - Gft*cos(thetaft)); */
    d2Qf_dthetat_dthetaf = Vmf*Vmt*(-Bft*PetscCosScalar(thetaft) + Gft*PetscSinScalar(thetaft));
    d2Qf_dthetat_dVmf    =     Vmt*(-Bft*PetscSinScalar(thetaft) - Gft*PetscCosScalar(thetaft));
    d2Qf_dthetat_dthetat = Vmf*Vmt*( Bft*PetscCosScalar(thetaft) - Gft*PetscSinScalar(thetaft));
    d2Qf_dthetat_dVmt    =     Vmf*(-Bft*PetscSinScalar(thetaft) - Gft*PetscCosScalar(thetaft));
    
    /* dQf_dVmt = Vmf*(-Bft*cos(thetaft) + Gft*sin(thetaft)); */
    d2Qf_dVmt_dthetaf    = Vmf*(Bft*PetscSinScalar(thetaft) + Gft*PetscCosScalar(thetaft)); 
    d2Qf_dVmt_dVmf       =    (-Bft*PetscCosScalar(thetaft) + Gft*PetscSinScalar(thetaft));
    d2Qf_dVmt_dthetat    = Vmf*(-Bft*PetscSinScalar(thetaft) - Gft*PetscCosScalar(thetaft));
    d2Qf_dVmt_dVmt       = 0.0;
    
    double d2Pt_dthetat_dthetat,d2Pt_dthetat_dVmt,d2Pt_dthetat_dthetaf,d2Pt_dthetat_dVmf;
    double d2Pt_dVmt_dthetat,   d2Pt_dVmt_dVmt,   d2Pt_dVmt_dthetaf,   d2Pt_dVmt_dVmf;
    double d2Pt_dthetaf_dthetat,d2Pt_dthetaf_dVmt,d2Pt_dthetaf_dthetaf,d2Pt_dthetaf_dVmf;
    double d2Pt_dVmf_dthetat,   d2Pt_dVmf_dVmt,   d2Pt_dVmf_dthetaf,   d2Pt_dVmf_dVmf;
    
    /* dPt_dthetat = Vmf*Vmt*(-Gtf*sin(thetatf) + Btf*cos(thetatf)); */
    d2Pt_dthetat_dthetat = Vmf*Vmt*(-Gtf*PetscCosScalar(thetatf) - Btf*PetscSinScalar(thetatf));
    d2Pt_dthetat_dVmt    =     Vmf*(-Gtf*PetscSinScalar(thetatf) + Btf*PetscCosScalar(thetatf));
    d2Pt_dthetat_dthetaf =  Vmf*Vmt*(Gtf*PetscCosScalar(thetatf) + Btf*PetscSinScalar(thetatf));
    d2Pt_dthetat_dVmf    =     Vmt*(-Gtf*PetscSinScalar(thetatf) + Btf*PetscCosScalar(thetatf));
    
    /* dPt_Vmt  = 2*Gtt*Vmt + Vmf*(Gtf*cos(thetatf) + Btf*sin(thetatf)); */
    d2Pt_dVmt_dthetat    =  Vmf*(-Gtf*PetscSinScalar(thetatf) + Bft*PetscCosScalar(thetatf)); 
    d2Pt_dVmt_dVmt       =  2*Gtt;
    d2Pt_dVmt_dthetaf    =  Vmf*(Gtf*PetscSinScalar(thetatf) - Btf*PetscCosScalar(thetatf));
    d2Pt_dVmt_dVmf       =      (Gtf*PetscCosScalar(thetatf) + Btf*PetscSinScalar(thetatf));
    
    /* dPt_dthetaf = Vmf*Vmt*(Gtf*sin(thetatf) - Btf*cos(thetatf)); */
    d2Pt_dthetaf_dthetat = Vmf*Vmt*(Gtf*PetscCosScalar(thetatf)  + Btf*PetscSinScalar(thetatf));
    d2Pt_dthetaf_dVmt    =     Vmf*(Gtf*PetscSinScalar(thetatf)  - Btf*PetscCosScalar(thetatf));
    d2Pt_dthetaf_dthetaf = Vmf*Vmt*(-Gtf*PetscCosScalar(thetatf) - Btf*PetscSinScalar(thetatf));
    d2Pt_dthetaf_dVmf    =     Vmt*(Gtf*PetscSinScalar(thetatf)  - Btf*PetscCosScalar(thetatf));
    
    /* dPt_dVmf = Vmt*(Gtf*cos(thetatf) + Btf*sin(thetatf)); */
    d2Pt_dVmf_dthetat    = Vmt*(-Gtf*PetscSinScalar(thetatf) + Btf*PetscCosScalar(thetatf)); 
    d2Pt_dVmf_dVmt       =     (Gtf*PetscCosScalar(thetatf) + Btf*PetscSinScalar(thetatf));
    d2Pt_dVmf_dthetaf    = Vmt*(Gtf*PetscSinScalar(thetatf) - Btf*PetscCosScalar(thetatf));
    d2Pt_dVmf_dVmf       = 0.0;
    
    double d2Qt_dthetaf_dthetaf,d2Qt_dthetaf_dVmf,d2Qt_dthetaf_dthetat,d2Qt_dthetaf_dVmt;
    double d2Qt_dVmf_dthetaf,   d2Qt_dVmf_dVmf,   d2Qt_dVmf_dthetat,   d2Qt_dVmf_dVmt;
    double d2Qt_dthetat_dthetaf,d2Qt_dthetat_dVmf,d2Qt_dthetat_dthetat,d2Qt_dthetat_dVmt;
    double d2Qt_dVmt_dthetaf,   d2Qt_dVmt_dVmf,   d2Qt_dVmt_dthetat,   d2Qt_dVmt_dVmt;
    
    /* dQt_dthetat = Vmf*Vmt*(Btf*sin(thetatf) + Gtf*cos(thetatf)); */
    d2Qt_dthetat_dthetat = Vmf*Vmt*(Btf*PetscCosScalar(thetatf)  - Gtf*PetscSinScalar(thetatf));
    d2Qt_dthetat_dVmt    =     Vmf*(Btf*PetscSinScalar(thetatf)  + Gtf*PetscCosScalar(thetatf));
    d2Qt_dthetat_dthetaf = Vmf*Vmt*(-Btf*PetscCosScalar(thetatf) + Gtf*PetscSinScalar(thetatf));
    d2Qt_dthetat_dVmf    =     Vmf*(Btf*PetscSinScalar(thetatf) + Gtf*PetscCosScalar(thetatf));
    
    /* dQt_dVmt = -2*Btt*Vmt + Vmf*(-Btf*cos(thetatf) + Gtf*sin(thetatf)); */
    d2Qt_dVmt_dthetat    =  Vmf*(Btf*PetscSinScalar(thetatf) + Gtf*PetscCosScalar(thetatf)); 
    d2Qt_dVmt_dVmt       = -2*Btt;
    d2Qt_dVmt_dthetaf    =  Vmf*(-Btf*PetscSinScalar(thetatf) + Gtf*PetscCosScalar(thetatf));
    d2Qt_dVmt_dVmf       =      (-Btf*PetscCosScalar(thetatf) + Gtf*PetscSinScalar(thetatf));
    
    /* dQt_dthetaf = Vmf*Vmt*(-Btf*sin(thetatf) - Gtf*cos(thetatf)); */
    d2Qt_dthetaf_dthetat = Vmf*Vmt*(-Btf*PetscCosScalar(thetatf) + Gtf*PetscSinScalar(thetatf));
    d2Qt_dthetaf_dVmt    =     Vmf*(-Btf*PetscSinScalar(thetatf) - Gtf*PetscCosScalar(thetatf));
    d2Qt_dthetaf_dthetaf = Vmf*Vmt*( Btf*PetscCosScalar(thetatf) - Gtf*PetscSinScalar(thetatf));
    d2Qt_dthetaf_dVmf    =     Vmt*(-Btf*PetscSinScalar(thetatf) - Gtf*PetscCosScalar(thetatf));
    
    /* dQt_dVmf = Vmt*(-Btf*cos(thetatf) + Gtf*sin(thetatf)); */
    d2Qt_dVmf_dthetat    = Vmt*(Btf*PetscSinScalar(thetatf) + Gtf*PetscCosScalar(thetatf)); 
    d2Qt_dVmf_dVmt       =    (-Btf*PetscCosScalar(thetatf) + Gtf*PetscSinScalar(thetatf));
    d2Qt_dVmf_dthetaf    = Vmt*(-Btf*PetscSinScalar(thetatf) - Gtf*PetscCosScalar(thetatf));
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

    row[0] = lineparams->xidxf[j];
    col[0] = lineparams->xidxf[j];
    col[1] = col[0]+1;
    col[2] = lineparams->xidxt[j];
    col[3] = col[2]+1;

    gloc = lineparams->gineqidx[i];

    val[0] = lambda[gloc]*d2Sf2_dthetaf_dthetaf + lambda[gloc+1]*d2St2_dthetaf_dthetaf;
    val[1] = lambda[gloc]*d2Sf2_dthetaf_dVmf + lambda[gloc+1]*d2St2_dthetaf_dVmf;
    val[2] = lambda[gloc]*d2Sf2_dthetaf_dthetat + lambda[gloc+1]*d2St2_dthetaf_dthetat;
    val[3] = lambda[gloc]*d2Sf2_dthetaf_dVmt + lambda[gloc+1]*d2St2_dthetaf_dVmt;
      
    ierr = MatSetValues(H,1,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);
      
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
    col[0] = lineparams->xidxf[j];
    col[1] = col[0]+1;
    col[2] = lineparams->xidxt[j];
    col[3] = col[2]+1;

    row[0] = lineparams->xidxf[j]+1;

    val[0] = lambda[gloc]*d2Sf2_dVmf_dthetaf + lambda[gloc+1]*d2St2_dVmf_dthetaf;
    val[1] = lambda[gloc]*d2Sf2_dVmf_dVmf + lambda[gloc+1]*d2St2_dVmf_dVmf;
    val[2] = lambda[gloc]*d2Sf2_dVmf_dthetat + lambda[gloc+1]*d2St2_dVmf_dthetat;
    val[3] = lambda[gloc]*d2Sf2_dVmf_dVmt + lambda[gloc+1]*d2St2_dVmf_dVmt;

    ierr = MatSetValues(H,1,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);
      
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

    col[0] = lineparams->xidxf[j];
    col[1] = col[0]+1;
    col[2] = lineparams->xidxt[j];
    col[3] = col[2]+1;

    row[0] = lineparams->xidxt[j];
      
    val[0] = lambda[gloc]*d2Sf2_dthetat_dthetaf + lambda[gloc+1]*d2St2_dthetat_dthetaf;
    val[1] = lambda[gloc]*d2Sf2_dthetat_dVmf + lambda[gloc+1]*d2St2_dthetat_dVmf;
    val[2] = lambda[gloc]*d2Sf2_dthetat_dthetat + lambda[gloc+1]*d2St2_dthetat_dthetat;
    val[3] = lambda[gloc]*d2Sf2_dthetat_dVmt + lambda[gloc+1]*d2St2_dthetat_dVmt;
    
    ierr = MatSetValues(H,1,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);
    
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
    row[0] = lineparams->xidxt[j];
    col[0] = lineparams->xidxf[j];
    col[1] = col[0]+1;
    col[2] = lineparams->xidxt[j];
    col[3] = col[2]+1;
      
    val[0] = lambda[gloc]*d2Sf2_dVmt_dthetaf + lambda[gloc+1]*d2St2_dVmt_dthetaf;
    val[1] = lambda[gloc]*d2Sf2_dVmt_dVmf + lambda[gloc+1]*d2St2_dVmt_dVmf;
    val[2] = lambda[gloc]*d2Sf2_dVmt_dthetat + lambda[gloc+1]*d2St2_dVmt_dthetat;
    val[3] = lambda[gloc]*d2Sf2_dVmt_dVmt + lambda[gloc+1]*d2St2_dVmt_dVmt;
    
    ierr = MatSetValues(H,1,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);
  }

  ierr = VecRestoreArray(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(Lambda,&lambda);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#endif
