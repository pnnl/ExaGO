static char help[] = "Prototype version of OPFLOW code for CUDA/GPU/Accelerator implementation\n\
                      In this implementation, component-class wise assembly of the different parts\n\
                      of the OPFLOW calculations\n\n";

#include <scopflow_config.h>
#include <private/psimpl.h>
#include <private/opflowimpl.h>

typedef struct {
  int    nbus;    /* Number of buses */
  int    *isref;  /* isref[i] = 1 if bus is reference bus */
  int    *isisolated; /* isisolated[i] = 1 if bus is isolated bus */
  int    *ispvpq; /* For all other buses */
  double *vmin; /* min. voltage magnitude limit */
  double *vmax; /* max. voltage magnitude limit */
  double *va;   /* bus angle (from file only used in bounds) */
  double *vm;   /* bus voltage magnitude (from file only used in bounds) */
  double *gl;  /* bus shunt (conductance) */
  double *bl;  /* bus shunt (suspectance) */
  int    *xidx; /* starting locations for bus variables in X vector */
  int    *gidx;  /* starting locations for bus balance equations in constraint vector */
}BUSParams;

typedef struct {
  int    ngenON;       /* Number of generators with STATUS ON */
  double *cost_alpha;  /* generator cost coefficients */
  double *cost_beta;   /* generator cost coefficients */
  double *cost_gamma;  /* generator cost coefficients */
  double *pt;          /* min. active power gen. limits */
  double *pb;          /* max. active power gen. limits */
  double *qt;          /* min. reactive power gen. limits */
  double *qb;          /* max. reactive power gen. limits */
  int    *xidx;        /* starting locations in X vector */
  int    *gidx;         /* starting locations in constraint vector */
}GENParams;

typedef struct{
  int    nload; /* Number of loads */
  double *pl;   /* active power demand */
  double *ql;   /* reactive power demand */
  int    *xidx; /* starting location in X vector */
  int    *gidx;  /* starting location in constraint vector */
}LOADParams;

typedef struct{
  int     nlineON; /* Number of active lines (STATUS = 1) */
  int     nlinelim; /* Active lines + limits */
  double *Gff;  /* From side self conductance */
  double *Bff;  /* From side self susceptance */
  double *Gft;  /* From-to transfer conductance */
  double *Bft;  /* From-to transfer susceptance */
  double *Gtf;  /* To-from transfer conductance */
  double *Btf;  /* To-from transfer susceptance */
  double *Gtt;  /* To side self conductance */
  double *Btt;  /* To side self susceptance */
  double *rateA; /* Line MVA rating A (normal operation) */
  int    *xidxf; /* Starting locatin of from bus voltage variables */
  int    *xidxt; /* Starting location of to bus voltage variables */
  int    *geqidxf; /* Starting location of from side to insert equality constraint contribution in constraints vector */
  int    *geqidxt; /* Starting location of to side to insert equality constraint contribution in constraints vector */
  int    *gineqidx; /* Starting location to insert contribution to inequality constraint */
  int    *linelimidx; /* Indices for subset of lines that have finite limits */
}LINEParams;

PetscErrorCode DestroyBusParams(BUSParams *busparams)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = PetscFree(busparams->isref);CHKERRQ(ierr);
  ierr = PetscFree(busparams->isisolated);CHKERRQ(ierr);
  ierr = PetscFree(busparams->ispvpq);CHKERRQ(ierr);
  ierr = PetscFree(busparams->vmin);CHKERRQ(ierr);
  ierr = PetscFree(busparams->vmax);CHKERRQ(ierr);
  ierr = PetscFree(busparams->va);CHKERRQ(ierr);
  ierr = PetscFree(busparams->vm);CHKERRQ(ierr);
  ierr = PetscFree(busparams->gl);CHKERRQ(ierr);
  ierr = PetscFree(busparams->bl);CHKERRQ(ierr);
  ierr = PetscFree(busparams->xidx);CHKERRQ(ierr);
  ierr = PetscFree(busparams->gidx);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/* Create data for buses that is used in different computations */
PetscErrorCode CreateBusParams(OPFLOW opflow,BUSParams *busparams)
{
  PS             ps=opflow->ps;
  PetscInt       loc,gloc=0;
  PSBUS          bus;
  PetscInt       i;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  
  busparams->nbus = ps->nbus;

  /* Allocate the arrays */
  ierr = PetscCalloc1(busparams->nbus,&busparams->isref);CHKERRQ(ierr);
  ierr = PetscCalloc1(busparams->nbus,&busparams->isisolated);CHKERRQ(ierr);
  ierr = PetscCalloc1(busparams->nbus,&busparams->ispvpq);CHKERRQ(ierr);
  ierr = PetscCalloc1(busparams->nbus,&busparams->vmin);CHKERRQ(ierr);
  ierr = PetscCalloc1(busparams->nbus,&busparams->vmax);CHKERRQ(ierr);
  ierr = PetscCalloc1(busparams->nbus,&busparams->va);CHKERRQ(ierr);
  ierr = PetscCalloc1(busparams->nbus,&busparams->vm);CHKERRQ(ierr);
  ierr = PetscCalloc1(busparams->nbus,&busparams->gl);CHKERRQ(ierr);
  ierr = PetscCalloc1(busparams->nbus,&busparams->bl);CHKERRQ(ierr);
  ierr = PetscCalloc1(busparams->nbus,&busparams->xidx);CHKERRQ(ierr);
  ierr = PetscCalloc1(busparams->nbus,&busparams->gidx);CHKERRQ(ierr);

  /* Populate the arrays */
  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
    
    busparams->xidx[i] = loc;
    busparams->gidx[i] = gloc;

    if(bus->ide == REF_BUS) busparams->isref[i] = 1;
    else if(bus->ide == ISOLATED_BUS) busparams->isisolated[i] = 1;
    else busparams->ispvpq[i] = 1;

    busparams->vmin[i] = bus->Vmin;
    busparams->vmax[i] = bus->Vmax;
    busparams->vm[i]   = bus->vm;
    busparams->va[i]   = bus->va;
    busparams->gl[i]   = bus->gl;
    busparams->bl[i]   = bus->bl;

    gloc += 2;
  }
    
  PetscFunctionReturn(0);
}

PetscErrorCode DestroyLineParams(LINEParams *lineparams)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = PetscFree(lineparams->Gff);CHKERRQ(ierr);
  ierr = PetscFree(lineparams->Bff);CHKERRQ(ierr);
  ierr = PetscFree(lineparams->Gft);CHKERRQ(ierr);
  ierr = PetscFree(lineparams->Bft);CHKERRQ(ierr);
  ierr = PetscFree(lineparams->Gtf);CHKERRQ(ierr);
  ierr = PetscFree(lineparams->Btf);CHKERRQ(ierr);
  ierr = PetscFree(lineparams->Gtt);CHKERRQ(ierr);
  ierr = PetscFree(lineparams->Btt);CHKERRQ(ierr);
  ierr = PetscFree(lineparams->rateA);CHKERRQ(ierr);

  ierr = PetscFree(lineparams->xidxf);CHKERRQ(ierr);
  ierr = PetscFree(lineparams->xidxt);CHKERRQ(ierr);

  ierr = PetscFree(lineparams->geqidxf);CHKERRQ(ierr);
  ierr = PetscFree(lineparams->geqidxt);CHKERRQ(ierr);
  ierr = PetscFree(lineparams->gineqidx);CHKERRQ(ierr);

  ierr = PetscFree(lineparams->linelimidx);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/* Create data for lines that is used in different computations */
PetscErrorCode CreateLineParams(OPFLOW opflow,LINEParams *lineparams)
{
  PS             ps=opflow->ps;
  PetscInt       loc,gloc=opflow->nconeq,linei=0,linelimi=0;
  PetscInt       locf,loct;
  PSLINE         line;
  PSBUS          bus;
  PetscInt       i,j;
  PetscErrorCode ierr;
  const PSBUS    *connbuses;
  PSBUS          busf,bust;

  PetscFunctionBegin;

  ierr = PSGetNumActiveLines(ps,&lineparams->nlineON,NULL);CHKERRQ(ierr);

  lineparams->nlinelim = 0;
  /* Get the number of lines that are active and have finite limits. These lines
     will be only considered in inequality constraints */
  for(i=0; i < ps->nline; i++) {
    line = &ps->line[i];

    if(!line->status || line->rateA > 1e5) continue;
    lineparams->nlinelim++;
  }

  /* Allocate arrays */
  ierr = PetscCalloc1(lineparams->nlineON,&lineparams->Gff);CHKERRQ(ierr);
  ierr = PetscCalloc1(lineparams->nlineON,&lineparams->Bff);CHKERRQ(ierr);
  ierr = PetscCalloc1(lineparams->nlineON,&lineparams->Gft);CHKERRQ(ierr);
  ierr = PetscCalloc1(lineparams->nlineON,&lineparams->Bft);CHKERRQ(ierr);
  ierr = PetscCalloc1(lineparams->nlineON,&lineparams->Gtf);CHKERRQ(ierr);
  ierr = PetscCalloc1(lineparams->nlineON,&lineparams->Btf);CHKERRQ(ierr);
  ierr = PetscCalloc1(lineparams->nlineON,&lineparams->Gtt);CHKERRQ(ierr);
  ierr = PetscCalloc1(lineparams->nlineON,&lineparams->Btt);CHKERRQ(ierr);
  ierr = PetscCalloc1(lineparams->nlineON,&lineparams->rateA);CHKERRQ(ierr);

  ierr = PetscCalloc1(lineparams->nlineON,&lineparams->xidxf);CHKERRQ(ierr);
  ierr = PetscCalloc1(lineparams->nlineON,&lineparams->xidxt);CHKERRQ(ierr);

  ierr = PetscCalloc1(lineparams->nlineON,&lineparams->geqidxf);CHKERRQ(ierr);
  ierr = PetscCalloc1(lineparams->nlineON,&lineparams->geqidxt);CHKERRQ(ierr);
  ierr = PetscCalloc1(lineparams->nlineON,&lineparams->gineqidx);CHKERRQ(ierr);

  ierr = PetscCalloc1(lineparams->nlinelim,&lineparams->linelimidx);CHKERRQ(ierr);

  /* Populate arrays */
  for(i=0; i < ps->nline; i++) {
    line = &ps->line[i];

    if(!line->status) continue;

    lineparams->Gff[linei] = line->yff[0];
    lineparams->Bff[linei] = line->yff[1];
    lineparams->Gft[linei] = line->yft[0];
    lineparams->Bft[linei] = line->yft[1];
    lineparams->Gtf[linei] = line->ytf[0];
    lineparams->Btf[linei] = line->ytf[1];
    lineparams->Gtt[linei] = line->ytt[0];
    lineparams->Btt[linei] = line->ytt[1];
    lineparams->rateA[linei] = line->rateA;

    ierr = PSLINEGetConnectedBuses(line,&connbuses);CHKERRQ(ierr);
    busf = connbuses[0];
    bust = connbuses[1];

    ierr = PSBUSGetVariableLocation(busf,&lineparams->xidxf[linei]);CHKERRQ(ierr);
    ierr = PSBUSGetVariableLocation(bust,&lineparams->xidxt[linei]);CHKERRQ(ierr);

    /* 
       Each bus has two equality (balance) constraints, hence the use of coefficient 2
       to map the location of the equality constraint for the bus
    */
    lineparams->geqidxf[linei] = 2*busf->internal_i;
    lineparams->geqidxt[linei] = 2*bust->internal_i;

    if(line->rateA < 1e5) {
      lineparams->gineqidx[linelimi] = gloc;
      lineparams->linelimidx[linelimi] = linei;
      linelimi++;
      gloc += 2;
    }
    linei++;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode DestroyLoadParams(LOADParams *loadparams)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscFree(loadparams->pl);CHKERRQ(ierr);
  ierr = PetscFree(loadparams->ql);CHKERRQ(ierr);
  ierr = PetscFree(loadparams->xidx);CHKERRQ(ierr);
  ierr = PetscFree(loadparams->gidx);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* Create data for loads that is used in different computations */
PetscErrorCode CreateLoadParams(OPFLOW opflow,LOADParams *loadparams)
{
  PS             ps=opflow->ps;
  PetscInt       loc,gloc=0,loadi=0;
  PSLOAD         load;
  PSBUS          bus;
  PetscInt       i,j;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  
  /* Get the number of active generators (STATUS ON) */
  ierr = PSGetNumLoads(ps,&loadparams->nload,NULL);CHKERRQ(ierr);

  /* Allocate arrays */
  ierr = PetscCalloc1(loadparams->nload,&loadparams->pl);CHKERRQ(ierr);
  ierr = PetscCalloc1(loadparams->nload,&loadparams->ql);CHKERRQ(ierr);
  ierr = PetscCalloc1(loadparams->nload,&loadparams->xidx);CHKERRQ(ierr);
  ierr = PetscCalloc1(loadparams->nload,&loadparams->gidx);CHKERRQ(ierr);

  /* Insert data in loadparams */
  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
    for(j=0; j < bus->nload; j++) {
      ierr = PSBUSGetLoad(bus,j,&load);CHKERRQ(ierr);
      
      loc += 2;

      loadparams->pl[loadi] = load->pl;
      loadparams->ql[loadi] = load->ql;

      loadparams->xidx[loadi] = loc;
      loadparams->gidx[loadi] = gloc;

      loadi++;
    }
    gloc += 2;
  }

  PetscFunctionReturn(0);
}

PetscErrorCode DestroyGenParams(GENParams *genparams)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = PetscFree(genparams->cost_alpha);CHKERRQ(ierr);
  ierr = PetscFree(genparams->cost_beta);CHKERRQ(ierr);
  ierr = PetscFree(genparams->cost_gamma);CHKERRQ(ierr);

  ierr = PetscFree(genparams->pt);CHKERRQ(ierr);
  ierr = PetscFree(genparams->pb);CHKERRQ(ierr);
  ierr = PetscFree(genparams->qt);CHKERRQ(ierr);
  ierr = PetscFree(genparams->qb);CHKERRQ(ierr);

  ierr = PetscFree(genparams->xidx);CHKERRQ(ierr);
  ierr = PetscFree(genparams->gidx);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
/* Create data for generators that is used in different computations */
PetscErrorCode CreateGenParams(OPFLOW opflow,GENParams *genparams)
{
  PS             ps=opflow->ps;
  PetscInt       loc,gloc=0,geni=0;
  PSGEN          gen;
  PSBUS          bus;
  PetscInt       i,j;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  
  /* Get the number of active generators (STATUS ON) */
  ierr = PSGetNumActiveGenerators(ps,&genparams->ngenON,NULL);CHKERRQ(ierr);

  /* Allocate the arrays */
  ierr = PetscCalloc1(genparams->ngenON,&genparams->cost_alpha);CHKERRQ(ierr);
  ierr = PetscCalloc1(genparams->ngenON,&genparams->cost_beta);CHKERRQ(ierr);
  ierr = PetscCalloc1(genparams->ngenON,&genparams->cost_gamma);CHKERRQ(ierr);

  ierr = PetscCalloc1(genparams->ngenON,&genparams->pt);CHKERRQ(ierr);
  ierr = PetscCalloc1(genparams->ngenON,&genparams->pb);CHKERRQ(ierr);
  ierr = PetscCalloc1(genparams->ngenON,&genparams->qt);CHKERRQ(ierr);
  ierr = PetscCalloc1(genparams->ngenON,&genparams->qb);CHKERRQ(ierr);

  ierr = PetscCalloc1(genparams->ngenON,&genparams->xidx);CHKERRQ(ierr);
  ierr = PetscCalloc1(genparams->ngenON,&genparams->gidx);CHKERRQ(ierr);

  /* Insert data in genparams */
  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
    for(j=0; j < bus->ngen; j++) {
      ierr = PSBUSGetGen(bus,j,&gen);CHKERRQ(ierr);
      if(!gen->status) continue;
      
      loc += 2;

      genparams->cost_alpha[geni] = gen->cost_alpha;
      genparams->cost_beta[geni]  = gen->cost_beta;
      genparams->cost_gamma[geni] = gen->cost_gamma;
      genparams->pt[geni]         = gen->pt;
      genparams->pb[geni]         = gen->pb;
      genparams->qt[geni]         = gen->qt;
      genparams->qb[geni]         = gen->qb;

      genparams->xidx[geni]       = loc;
      genparams->gidx[geni]       = gloc;

      geni++;
    }
    gloc += 2;
  }

  PetscFunctionReturn(0);
}

void ComputeGenEqualityConstraintContributionKernel(int *xidx,int *gidx, int i, double *x, double *g)
{
  /* atomic operation */
  g[gidx[i]]   -= x[xidx[i]];
  g[gidx[i]+1] -= x[xidx[i]+1];
}

void ComputeLoadEqualityConstraintContributionKernel(double *pl,double *ql,int *xidx,int *gidx, int i, double *x, double *g)
{
  /* Atomic operation */
  g[gidx[i]]   += pl[i];
  g[gidx[i]+1] += ql[i];
}

void ComputeBusEqualityConstraintContributionKernel(double *gl,double *bl,double *vm,double *va,int *isisolated, int* ispvpq,int *xidx, int* gidx, int i, double *x, double *g)
{
  double theta= x[xidx[i]];
  double Vm   = x[xidx[i]+1];

  g[gidx[i]]   += isisolated[i]*(theta - va[i]*PETSC_PI/180.0) + ispvpq[i]*Vm*Vm*gl[i];
  g[gidx[i]+1] += isisolated[i]*(Vm    - vm[i]) - ispvpq[i]*Vm*Vm*bl[i];
}

void ComputeLineEqualityConstraintContributionKernel(double *Gff,double *Bff,double *Gft,double *Bft,double *Gtf,double *Btf,double *Gtt,double *Btt,int *xidxf,int *xidxt,int *geqidxf,int *geqidxt,int i,double *x,double *g)
{
  double Pf,Qf,Pt,Qt;
  double thetaf=x[xidxf[i]], Vmf=x[xidxf[i]+1];
  double thetat=x[xidxt[i]], Vmt=x[xidxt[i]+1];
  double thetaft=thetaf-thetat;
  double thetatf=thetat-thetaf;

  Pf = Gff[i]*Vmf*Vmf  + Vmf*Vmt*(Gft[i]*cos(thetaft) + Bft[i]*sin(thetaft));
  Qf = -Bff[i]*Vmf*Vmf + Vmf*Vmt*(-Bft[i]*cos(thetaft) + Gft[i]*sin(thetaft));
  Pt = Gtt[i]*Vmt*Vmt  + Vmt*Vmf*(Gtf[i]*cos(thetatf) + Btf[i]*sin(thetatf));
  Qt = -Btt[i]*Vmt*Vmt + Vmt*Vmf*(-Btf[i]*cos(thetatf) + Gtf[i]*sin(thetatf));

  /* Atomic operation */
  g[geqidxf[i]]   += Pf;
  g[geqidxf[i]+1] += Qf;
  g[geqidxt[i]]   += Pt;
  g[geqidxt[i]+1] += Qt;
}

void ComputeLineInequalityConstraintKernel(double *Gff,double *Bff,double *Gft,double *Bft,double *Gtf,double *Btf,double *Gtt,double *Btt,int *linelimidx,int *xidxf,int *xidxt,int *gineqidx,int i,double *x, double *g)
{
  int    j=linelimidx[i];
  double Pf,Qf,Pt,Qt,Sf2,St2;
  double thetaf=x[xidxf[j]], Vmf=x[xidxf[j]+1];
  double thetat=x[xidxt[j]], Vmt=x[xidxt[j]+1];
  double thetaft=thetaf-thetat;
  double thetatf=thetat-thetaf;

  Pf = Gff[j]*Vmf*Vmf  + Vmf*Vmt*(Gft[j]*cos(thetaft) + Bft[j]*sin(thetaft));
  Qf = -Bff[j]*Vmf*Vmf + Vmf*Vmt*(-Bft[j]*cos(thetaft) + Gft[j]*sin(thetaft));
  Pt = Gtt[j]*Vmt*Vmt  + Vmt*Vmf*(Gtf[j]*cos(thetatf) + Btf[j]*sin(thetatf));
  Qt = -Btt[j]*Vmt*Vmt + Vmt*Vmf*(-Btf[j]*cos(thetatf) + Gtf[j]*sin(thetatf));

  Sf2 = Pf*Pf + Qf*Qf;
  St2 = Pt*Pt + Qt*Qt;

  g[gineqidx[i]]   = Sf2;
  g[gineqidx[i]+1] = St2;
}

void ComputeInequalityConstraintBoundKernel(double *rateA,double MVAbase,int *linelimidx,int *gineqidx,int i,double *gl,double *gu)
{
  int    j=linelimidx[i];
  
  gl[gineqidx[i]]   = 0.0;
  gu[gineqidx[i]]   = (rateA[j]/MVAbase)*(rateA[j]/MVAbase);
  gl[gineqidx[i]+1] = 0.0;
  gu[gineqidx[i]+1] = (rateA[j]/MVAbase)*(rateA[j]/MVAbase);
}


PetscErrorCode ComputeConstraintBounds(OPFLOW opflow,LINEParams *lineparams,Vec X, Vec Gl, Vec Gu)
{
  PetscErrorCode ierr;
  PetscScalar    *x,*gl,*gu;
  PetscInt       i;
  PS             ps=opflow->ps;
  PetscReal      gldiff,gudiff;

  PetscFunctionBegin;

  ierr = VecSet(Gl,0.0);
  ierr = VecSet(Gu,0.0);

  ierr = VecGetArray(Gl,&gl);CHKERRQ(ierr);
  ierr = VecGetArray(Gu,&gu);CHKERRQ(ierr);

  /* Inequality constraints */
  for(i=0; i < lineparams->nlinelim; i++) {
    ComputeInequalityConstraintBoundKernel(lineparams->rateA,ps->MVAbase,lineparams->linelimidx,lineparams->gineqidx,i,gl,gu);
  }

  ierr = VecGetArray(Gl,&gl);CHKERRQ(ierr);
  ierr = VecGetArray(Gu,&gu);CHKERRQ(ierr);

  ierr = VecAXPY(Gl,-1.0,opflow->Gl);CHKERRQ(ierr);
  ierr = VecNorm(Gl,NORM_INFINITY,&gldiff);CHKERRQ(ierr); 
  if(gldiff < 1e-6) {
    ierr = PetscPrintf(PETSC_COMM_SELF,"Constraint lower bound correct\n");CHKERRQ(ierr);
  } else {
    ierr = PetscPrintf(PETSC_COMM_SELF,"Constraint lower bound incorrect\n");CHKERRQ(ierr);
  }
  ierr = VecAXPY(Gu,-1.0,opflow->Gu);CHKERRQ(ierr);
  ierr = VecNorm(Gu,NORM_INFINITY,&gudiff);CHKERRQ(ierr); 
  if(gudiff < 1e-6) {
    ierr = PetscPrintf(PETSC_COMM_SELF,"Constraint upper bound correct\n");CHKERRQ(ierr);
  } else {
    ierr = PetscPrintf(PETSC_COMM_SELF,"Constraint upper bound incorrect\n");CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}
  
PetscErrorCode ComputeConstraints(OPFLOW opflow,BUSParams *busparams,LINEParams *lineparams,GENParams *genparams,LOADParams *loadparams,Vec X, Vec G)
{
  PetscErrorCode ierr;
  PetscScalar     *x,*g;
  PS              ps=opflow->ps;
  PetscInt        i;
  Vec             opflow_G;
  PetscReal       gdiff;

  PetscFunctionBegin;
  ierr = VecSet(G,0.0);CHKERRQ(ierr);
  ierr = VecGetArray(X,&x);CHKERRQ(ierr);
  ierr = VecGetArray(G,&g);CHKERRQ(ierr);

  /* Equality constraints */
  /* Generator contributions */
  for(i=0; i < genparams->ngenON; i++) {
    ComputeGenEqualityConstraintContributionKernel(genparams->xidx,genparams->gidx,i,x,g);
  }

  /* Load contributions */
  for(i=0; i < loadparams->nload; i++) {
    ComputeLoadEqualityConstraintContributionKernel(loadparams->pl,loadparams->ql,loadparams->xidx,loadparams->gidx,i,x,g);
  }

  /* Bus contributions */
  for(i=0; i < busparams->nbus; i++) {
    ComputeBusEqualityConstraintContributionKernel(busparams->gl,busparams->bl,busparams->vm,busparams->va,busparams->isisolated,busparams->ispvpq,busparams->xidx,busparams->gidx,i,x,g);
  }

  for(i=0; i < lineparams->nlineON; i++) {
    ComputeLineEqualityConstraintContributionKernel(lineparams->Gff,lineparams->Bff,lineparams->Gft,lineparams->Bft,lineparams->Gtf,lineparams->Btf,lineparams->Gtt,lineparams->Btt,lineparams->xidxf,lineparams->xidxt,lineparams->geqidxf,lineparams->geqidxt,i,x,g);
  }

  /* Inequality constraints */
  for(i=0; i < lineparams->nlinelim; i++) {
    ComputeLineInequalityConstraintKernel(lineparams->Gff,lineparams->Bff,lineparams->Gft,lineparams->Bft,lineparams->Gtf,lineparams->Btf,lineparams->Gtt,lineparams->Btt,lineparams->linelimidx,lineparams->xidxf,lineparams->xidxt,lineparams->gineqidx,i,x,g);
  }

  ierr = VecRestoreArray(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(G,&g);CHKERRQ(ierr);

  ierr = OPFLOWGetConstraints(opflow,&opflow_G);CHKERRQ(ierr);

  ierr = VecAXPY(G,-1.0,opflow_G);CHKERRQ(ierr);
  ierr = VecNorm(G,NORM_INFINITY,&gdiff);CHKERRQ(ierr);
  if(gdiff < 1e-6) {
    ierr = PetscPrintf(PETSC_COMM_SELF,"Constraint calculation correct\n");CHKERRQ(ierr);
  } else {
    ierr = PetscPrintf(PETSC_COMM_SELF,"Constraint calculation incorrect\n");CHKERRQ(ierr);
    ierr = VecView(G,0);CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

void ComputeGenObjectiveKernel(double *alpha,double *beta,double *gamma,int *xidx,double MVAbase,int i,double *x,double *obj)
{
  double Pg = x[xidx[i]]*MVAbase;

  obj[0] += alpha[i]*Pg*Pg + beta[i]*Pg + gamma[i];
}

PetscErrorCode ComputeObjective(OPFLOW opflow,GENParams *genparams,Vec X,double *obj)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscScalar    *x;
  PS             ps=opflow->ps;
  PetscScalar    opt_obj; /* Optimal objective value calculated by OPFLOW */

  PetscFunctionBegin;
  ierr = VecGetArray(X,&x);CHKERRQ(ierr);
  /* Generator objective function contributions */
  for(i=0; i < genparams->ngenON; i++) {
    ComputeGenObjectiveKernel(genparams->cost_alpha,genparams->cost_beta,genparams->cost_gamma,genparams->xidx,ps->MVAbase,i,x,obj);
  }

  ierr = OPFLOWGetObjective(opflow,&opt_obj);CHKERRQ(ierr);

  if(PetscAbsScalar(opt_obj - *obj)/opt_obj < 1e-6) {
    ierr = PetscPrintf(PETSC_COMM_SELF,"Objective value calculation correct: objective value = %f\n",opt_obj);CHKERRQ(ierr);
  } else {
    ierr = PetscPrintf(PETSC_COMM_SELF,"Optimal value calculation incorrrect: objective value = %f, pre-calculated value = %f\n",*obj,opt_obj);CHKERRQ(ierr);
  }
  ierr = VecRestoreArray(X,&x);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

void ComputeGenGradientKernel(double *alpha,double *beta,double *gamma,int *xidx,double MVAbase,int i,double *x,double *grad)
{
  double Pg = x[xidx[i]]*MVAbase;

  grad[xidx[i]] = MVAbase*(2.0*alpha[i]*Pg + beta[i]);
}

PetscErrorCode ComputeGradient(OPFLOW opflow,GENParams *genparams,Vec X,Vec Grad)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscScalar    *x,*grad;
  PS             ps=opflow->ps;

  PetscFunctionBegin;
  ierr = VecGetArray(X,&x);CHKERRQ(ierr);
  ierr = VecGetArray(Grad,&grad);CHKERRQ(ierr);

  /* Generator gradient contributions */
  for(i=0; i < genparams->ngenON; i++) {
    ComputeGenGradientKernel(genparams->cost_alpha,genparams->cost_beta,genparams->cost_gamma,genparams->xidx,ps->MVAbase,i,x,grad);
  }

  ierr = VecRestoreArray(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(Grad,&grad);CHKERRQ(ierr);

  /*  ierr = PetscPrintf(PETSC_COMM_SELF,"Gradient vector follows:\n");CHKERRQ(ierr);
  ierr = VecView(Grad,0);CHKERRQ(ierr);
  */
  PetscFunctionReturn(0);
}

void ComputeBusVariableBoundsKernel(int* isref,int* isisolated,int* ispvpq,double* vmin,double* vmax,double* vm,double* va,int *xidx,int i,double *xl,double *xu)
{
  xl[xidx[i]] = ispvpq[i]*PETSC_NINFINITY + isisolated[i]*va[i] + isref[i]*va[i]*PETSC_PI/180.0;
  xu[xidx[i]] = ispvpq[i]*PETSC_INFINITY  + isisolated[i]*va[i] + isref[i]*va[i]*PETSC_PI/180.0;

  xl[xidx[i]+1] = isref[i]*vmin[i]  + ispvpq[i]*vmin[i] + isisolated[i]*vm[i];
  xu[xidx[i]+1] = isref[i]*vmax[i]  + ispvpq[i]*vmax[i] + isisolated[i]*vm[i];

}

void ComputeGenVariableBoundsKernel(double *pb,double *pt,double *qb,double *qt,int *xidx,int i,double *xl, double *xu)
{
  xl[xidx[i]]   = pb[i];
  xu[xidx[i]]   = pt[i];
  xl[xidx[i]+1] = qb[i];
  xu[xidx[i]+1] = qt[i];
}

PetscErrorCode ComputeVariableBounds(OPFLOW opflow,BUSParams *busparams,GENParams *genparams,Vec Xl,Vec Xu)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscScalar    *x,*xl,*xu;
  PS             ps=opflow->ps;
  PetscReal      xldiff,xudiff;

  PetscFunctionBegin;
  ierr = VecGetArray(Xl,&xl);CHKERRQ(ierr);
  ierr = VecGetArray(Xu,&xu);CHKERRQ(ierr);

  /* Bounds for bus voltages */
  for(i=0; i < busparams->nbus; i++) {
    ComputeBusVariableBoundsKernel(busparams->isref,busparams->isisolated,busparams->ispvpq,busparams->vmin,busparams->vmax,busparams->vm,busparams->va,busparams->xidx,i,xl,xu);
  }

  /* Generator lower and upper bounds on variables */
  for(i=0; i < genparams->ngenON; i++) {
    ComputeGenVariableBoundsKernel(genparams->pb,genparams->pt,genparams->qb,genparams->qt,genparams->xidx,i,xl,xu);
  }

  ierr = VecRestoreArray(Xl,&xl);CHKERRQ(ierr);
  ierr = VecRestoreArray(Xu,&xu);CHKERRQ(ierr);

  ierr = VecAXPY(Xl,-1.0,opflow->Xl);CHKERRQ(ierr);
  ierr = VecNorm(Xl,NORM_INFINITY,&xldiff);CHKERRQ(ierr); 
  if(xldiff < 1e-6) {
    ierr = PetscPrintf(PETSC_COMM_SELF,"Variable lower bound correct\n");CHKERRQ(ierr);
  } else {
    ierr = PetscPrintf(PETSC_COMM_SELF,"Variable lower bound incorrect\n");CHKERRQ(ierr);
  }
  ierr = VecAXPY(Xu,-1.0,opflow->Xu);CHKERRQ(ierr);
  ierr = VecNorm(Xu,NORM_INFINITY,&xudiff);CHKERRQ(ierr); 
  if(xudiff < 1e-6) {
    ierr = PetscPrintf(PETSC_COMM_SELF,"Variable upper bound correct\n");CHKERRQ(ierr);
  } else {
    ierr = PetscPrintf(PETSC_COMM_SELF,"Variable upper bound incorrect\n");CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}
    
int main(int argc,char **argv)
{
  PetscErrorCode    ierr;
  OPFLOW            opflow;
  PS                ps;
  char              file[PETSC_MAX_PATH_LEN];
  PetscBool         flg;
  double            obj = 0.0; /* Objective function */
  BUSParams         busparams;
  GENParams         genparams;
  LOADParams        loadparams;
  LINEParams        lineparams;
  Vec               X,G;
  Vec               Grad;
  Vec               Xl,Xu;
  Vec               Gl,Gu;

  char options_pathname[200] = SCOPFLOW_OPTIONS_DIR;
  char* filename = "/opflowoptions";
  printf("%s\n", options_pathname);
  printf("%s\n", filename);
  strcat(options_pathname, filename);
  printf("%s\n", options_pathname);

  PetscInitialize(&argc,&argv,options_pathname,help);

  /* Create OPFLOW object */
  ierr = OPFLOWCreate(PETSC_COMM_WORLD,&opflow);CHKERRQ(ierr);

  /* Get network data file from command line */
  ierr = PetscOptionsGetString(NULL,NULL,"-netfile",file,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
  /* Read Network Data file */
  if(flg) {
    ierr = OPFLOWReadMatPowerData(opflow,file);CHKERRQ(ierr);
  } else {
    ierr = OPFLOWReadMatPowerData(opflow,"datafiles/case9mod.m");CHKERRQ(ierr);
  }

  /* Solve */
  ierr = OPFLOWSolve(opflow);CHKERRQ(ierr);

  /* Get solution */
  ierr = OPFLOWGetSolution(opflow,&X);CHKERRQ(ierr);
  
  /* Duplicate vector to comute gradients */
  ierr = VecDuplicate(X,&Grad);CHKERRQ(ierr);

  /* Duplicate vectors for lower and upper bounds */
  ierr = VecDuplicate(X,&Xl);CHKERRQ(ierr);
  ierr = VecDuplicate(X,&Xu);CHKERRQ(ierr);

  /* Create params for buses */
  ierr = CreateBusParams(opflow,&busparams);CHKERRQ(ierr);

  /* Create params for gens */
  ierr = CreateGenParams(opflow,&genparams);CHKERRQ(ierr);

  /* Create params for loads */
  ierr = CreateLoadParams(opflow,&loadparams);CHKERRQ(ierr);

  /* Create params for lines */
  ierr = CreateLineParams(opflow,&lineparams);CHKERRQ(ierr);

  /* Compute global objective */
  ierr = ComputeObjective(opflow,&genparams,X,&obj);CHKERRQ(ierr);

  /* Compute global gradient */
  ierr = ComputeGradient(opflow,&genparams,X,Grad);CHKERRQ(ierr);

  /* Compute variable bounds */
  ierr = ComputeVariableBounds(opflow,&busparams,&genparams,Xl,Xu);CHKERRQ(ierr);

  /* Duplicate constraint vector */
  ierr = VecDuplicate(opflow->G,&G);CHKERRQ(ierr);
  /* Duplicate vectors for constraint lower and upper bounds */
  ierr = VecDuplicate(G,&Gl);CHKERRQ(ierr);
  ierr = VecDuplicate(G,&Gu);CHKERRQ(ierr);

  /* Compute Constraint bounds */
  ierr = ComputeConstraintBounds(opflow,&lineparams,X,Gl,Gu);CHKERRQ(ierr);
  /* Compute constraints */
  ierr = ComputeConstraints(opflow,&busparams,&lineparams,&genparams,&loadparams,X,G);CHKERRQ(ierr);

  ierr = VecDestroy(&Grad);CHKERRQ(ierr);
  ierr = VecDestroy(&Xl);CHKERRQ(ierr);
  ierr = VecDestroy(&Xu);CHKERRQ(ierr);
  ierr = VecDestroy(&G);CHKERRQ(ierr);
  ierr = VecDestroy(&Gl);CHKERRQ(ierr);
  ierr = VecDestroy(&Gu);CHKERRQ(ierr);

  ierr = DestroyBusParams(&busparams);CHKERRQ(ierr);
  ierr = DestroyLineParams(&lineparams);CHKERRQ(ierr);
  ierr = DestroyGenParams(&genparams);CHKERRQ(ierr);
  ierr = DestroyLoadParams(&loadparams);CHKERRQ(ierr);

  ierr = OPFLOWDestroy(&opflow);CHKERRQ(ierr);
  PetscFinalize();

  return 0;
}
