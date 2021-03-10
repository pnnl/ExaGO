#include <exago_config.h>

#if defined(EXAGO_ENABLE_HIOP)

#include <private/psimpl.h>
#include <private/opflowimpl.h>
#include "pbpolhiop.h"

/************* NOTE ***********************/
/* No Load loss or power imbalance variables considered yet */
/********************************************/

/* Functions to create and destroy data arrays for different
   component classes
*/
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
    
    busparams->xidx[i] = opflow->idxn2sd_map[loc];
    busparams->gidx[i] = gloc;

    if(bus->ide == REF_BUS) busparams->isref[i] = 1;
    else if(bus->ide == ISOLATED_BUS) busparams->isisolated[i] = 1;
    else busparams->ispvpq[i] = 1;

    if(opflow->genbusvoltagetype == FIXED_AT_SETPOINT) {
      if(bus->ide == REF_BUS || bus->ide == PV_BUS) {
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
  ierr = PetscFree(lineparams->gbineqidx);CHKERRQ(ierr);

  ierr = PetscFree(lineparams->linelimidx);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/* Create data for lines that is used in different computations */
PetscErrorCode CreateLineParams(OPFLOW opflow,LINEParams *lineparams)
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

  ierr = PSGetNumActiveLines(ps,&lineparams->nlineON,NULL);CHKERRQ(ierr);

  lineparams->nlinelim = 0;
  /* Get the number of lines that are active and have finite limits. These lines
     will be only considered in inequality constraints */
  if(opflow->nconineq) {
    for(i=0; i < ps->nline; i++) {
      line = &ps->line[i];
      
      if(!line->status || line->rateA > 1e5) continue;
      lineparams->nlinelim++;
    }
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
  ierr = PetscCalloc1(lineparams->nlineON,&lineparams->gbineqidx);CHKERRQ(ierr);

  if(opflow->nconineq) {
    ierr = PetscCalloc1(lineparams->nlinelim,&lineparams->linelimidx);CHKERRQ(ierr);
  }

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

    int xidxf,xidxt;
    ierr = PSBUSGetVariableLocation(busf,&xidxf);CHKERRQ(ierr);
    ierr = PSBUSGetVariableLocation(bust,&xidxt);CHKERRQ(ierr);

    lineparams->xidxf[linei] = opflow->idxn2sd_map[xidxf];
    lineparams->xidxt[linei] = opflow->idxn2sd_map[xidxt];

    /* 
       Each bus has two equality (balance) constraints, hence the use of coefficient 2
       to map the location of the equality constraint for the bus
    */
    lineparams->geqidxf[linei] = 2*busf->internal_i;
    lineparams->geqidxt[linei] = 2*bust->internal_i;
    
    if(opflow->nconineq) {
      if(line->rateA < 1e5) {
	lineparams->gbineqidx[linelimi] = gbloc;
	lineparams->gineqidx[linelimi] = gloc;
	lineparams->linelimidx[linelimi] = linei;
	linelimi++;
	gbloc += 2;
	gloc += 2;
      }
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

      loadparams->xidx[loadi] = opflow->idxn2sd_map[loc];
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

  ierr = PetscFree(genparams->jacsp_idx);CHKERRQ(ierr);
  ierr = PetscFree(genparams->jacsq_idx);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
/* Create data for generators that is used in different computations */
PetscErrorCode CreateGenParams(OPFLOW opflow,GENParams *genparams)
{
  PS             ps=opflow->ps;
  PetscInt       loc,gloc=0,geni=0,nnzs=0,gi;
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

  ierr = PetscCalloc1(genparams->ngenON,&genparams->jacsp_idx);CHKERRQ(ierr);
  ierr = PetscCalloc1(genparams->ngenON,&genparams->jacsq_idx);CHKERRQ(ierr);


  /* Insert data in genparams */
  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
    gi = 0;
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

      genparams->xidx[geni]       = opflow->idxn2sd_map[loc];
      genparams->gidx[geni]       = gloc;
      genparams->jacsp_idx[geni]  = nnzs + gi;
      genparams->jacsq_idx[geni]  = nnzs + bus->ngenON + gi;

      geni++;
      gi++;
    }
    nnzs += 2*bus->ngenON;
    gloc += 2;
  }

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSetVariableandConstraintBounds_PBPOLHIOP(OPFLOW opflow,Vec Xl,Vec Xu, Vec Gl, Vec Gu)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = OPFLOWSetVariableBounds_PBPOLHIOP(opflow,Xl,Xu);CHKERRQ(ierr);
  ierr = OPFLOWSetConstraintBounds_PBPOLHIOP(opflow,Gl,Gu);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSetInitialGuess_PBPOLHIOP(OPFLOW opflow,Vec X)
{
  PetscErrorCode ierr;
  PS             ps=opflow->ps;
  const PetscScalar    *xl,*xu;
  PetscScalar    *x;
  PetscInt       i;
  PSBUS          bus;
  PetscInt       loc,loc_nat;
  int            *idxn2sd_map = opflow->idxn2sd_map;

  PetscFunctionBegin;

  /* Get array pointers */
  ierr = VecGetArray(X,&x);CHKERRQ(ierr);
  ierr = VecGetArrayRead(opflow->Xl,&xl);CHKERRQ(ierr);
  ierr = VecGetArrayRead(opflow->Xu,&xu);CHKERRQ(ierr);
  
  for(i=0; i < ps->nbus; i++) {
    PetscInt k;

    bus = &ps->bus[i];

    ierr = PSBUSGetVariableLocation(bus,&loc_nat);CHKERRQ(ierr);

    loc = idxn2sd_map[loc_nat];

    if(bus->ide == ISOLATED_BUS) {
      x[loc] = bus->va*PETSC_PI/180.0;
      x[loc+1] = bus->vm;
    } else {
      if(opflow->initializationtype == OPFLOWINIT_MIDPOINT) {
	/* Initial guess for voltage angles and bounds on voltage magnitudes */
	x[loc]   = (xl[loc] + xu[loc])/2.0;
	x[loc+1] = (xl[loc+1] + xu[loc+1])/2.0;
      } else if(opflow->initializationtype == OPFLOWINIT_FROMFILE || opflow->initializationtype == OPFLOWINIT_ACPF) {
	x[loc] = bus->va*PETSC_PI/180.0;
	x[loc+1]   = PetscMax(bus->Vmin,PetscMin(bus->vm,bus->Vmax));
      } else if(opflow->initializationtype == OPFLOWINIT_FLATSTART) {
	x[loc] = 0.0;
	x[loc+1] = 1.0;
      }
    }

    for(k=0; k < bus->ngen; k++) {
      PSGEN gen;

      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      if(!gen->status) continue;

      loc_nat = loc_nat+2;

      loc = idxn2sd_map[loc_nat];

      if(opflow->initializationtype == OPFLOWINIT_MIDPOINT || opflow->initializationtype == OPFLOWINIT_FLATSTART) {
	x[loc]   = 0.5*(xl[loc] + xu[loc]);
	x[loc+1] = 0.5*(xl[loc+1] + xu[loc+1]);
      } else if(opflow->initializationtype == OPFLOWINIT_FROMFILE || opflow->initializationtype == OPFLOWINIT_ACPF) {
	x[loc] = PetscMax(gen->pb,PetscMin(gen->pg,gen->pt));
	x[loc+1] = PetscMax(gen->qb,PetscMin(gen->qg,gen->qt));
      }
    }

    if(opflow->include_loadloss_variables) {
      for(k=0; k < bus->nload; k++) {
	PSLOAD load;
	ierr = PSBUSGetLoad(bus,k,&load);CHKERRQ(ierr);
	loc_nat += 2;

	loc = idxn2sd_map[loc_nat];

	/* Initial value for real and reactive power load loss */
	x[loc] = 0.0;
	x[loc+1] = 0.0;
      }
    } 

    if(opflow->include_powerimbalance_variables) {
      loc_nat += 2;

      loc = idxn2sd_map[loc_nat];

      x[loc] = x[loc+1] = 0.0;
    }
  }

  ierr = VecRestoreArray(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(opflow->Xl,&xl);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(opflow->Xu,&xu);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSetInitialGuessArray_PBPOLHIOP(OPFLOW opflow,double *x0)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = VecPlaceArray(opflow->X,x0);CHKERRQ(ierr);
  ierr = (*opflow->modelops.setinitialguess)(opflow,opflow->X);CHKERRQ(ierr);
  ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeConstraints_PBPOLHIOP(OPFLOW opflow,Vec X,Vec G)
{
  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeObjandGradient_PBPOLHIOP(OPFLOW opflow,Vec X,PetscScalar *obj,Vec Grad)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = PetscLogEventBegin(opflow->objlogger,0,0,0,0);CHKERRQ(ierr);
  ierr = OPFLOWComputeObjective_PBPOLHIOP(opflow,X,obj);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(opflow->objlogger,0,0,0,0);CHKERRQ(ierr);

  ierr = PetscLogEventBegin(opflow->gradlogger,0,0,0,0);CHKERRQ(ierr);
  ierr = OPFLOWComputeGradient_PBPOLHIOP(opflow,X,Grad);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(opflow->gradlogger,0,0,0,0);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWModelSetNumVariables_PBPOLHIOP(OPFLOW opflow,PetscInt *busnvar,PetscInt *branchnvar,PetscInt *nx)
{
  PetscInt i,ngen,nload,k;
  PS       ps=opflow->ps;
  PSBUS    bus;
  PSGEN    gen;
  PetscErrorCode ierr;
  PetscBool isghost;

  PetscFunctionBegin;
  
  *nx = 0;
  /* No variables for the branches */
  for(i=0; i < ps->nline; i++) {
    branchnvar[i] = 0;
    *nx += branchnvar[i];
  }

  /* Variables for the buses */
  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSIsGhosted(bus,&isghost);CHKERRQ(ierr);
    //    if(isghost) continue;
    busnvar[i] = 2; /* 2 variables for the bus */
    ierr = PSBUSGetNGen(bus,&ngen);CHKERRQ(ierr);
    for(k=0; k < ngen; k++) {
      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      if(!gen->status) continue;
      busnvar[i] += 2; /* (2 variables for voltage + Pg, Qg for each gen) */
    }

    if(opflow->include_loadloss_variables) {
      ierr = PSBUSGetNLoad(bus,&nload);CHKERRQ(ierr);
      /* Load loss variables..Real and imaginary part of the load loss */
      busnvar[i] += 2*nload;
    }

    if(opflow->include_powerimbalance_variables) busnvar[i] += 2;
    if(!isghost) *nx += busnvar[i];
  }

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWModelSetNumConstraints_PBPOLHIOP(OPFLOW opflow,PetscInt *branchnconeq,PetscInt *busnconeq,PetscInt *nconeq,PetscInt *nconineq)
{
  PetscInt i;
  PS       ps=opflow->ps;
  PSBUS    bus;
  PSLINE   line;
  PetscErrorCode ierr;
  PetscBool isghost;

  PetscFunctionBegin;
  *nconeq = *nconineq = 0;

  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSIsGhosted(bus,&isghost);CHKERRQ(ierr);
    if(isghost) continue;
    *nconeq += 2;
  }

  if(!opflow->ignore_lineflow_constraints) {
    for(i=0; i < ps->nline; i++) {
      line = &ps->line[i];
      if(line->status && line->rateA < 1e5) *nconineq += 2; /* Line flow constraints */
    }
  }

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeHessian_PBPOLHIOP(OPFLOW opflow,Vec X,Vec Lambdae,Vec Lambdai,Mat H)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = MatZeroEntries(H);CHKERRQ(ierr);

  /* Objective function Hessian */
  ierr = OPFLOWComputeObjectiveHessian_PBPOLHIOP(opflow,X,H);CHKERRQ(ierr);

  /* Equality constraints Hessian */
  ierr = OPFLOWComputeEqualityConstraintsHessian_PBPOLHIOP(opflow,X,Lambdae,H);CHKERRQ(ierr);
  
  /* Inequality constraints Hessian */
  if(opflow->nconineq) {
    ierr = OPFLOWComputeInequalityConstraintsHessian_PBPOLHIOP(opflow,X,Lambdai,H);CHKERRQ(ierr);
  }

  ierr = MatAssemblyBegin(H,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(H,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolutionToPS_PBPOLHIOP(OPFLOW opflow)
{
  PetscErrorCode ierr;
  PS             ps=(PS)opflow->ps;
  PetscInt       i,k;
  Vec            X,Lambda;
  PSBUS          bus;
  PSGEN          gen;
  PSLOAD         load;
  PSLINE         line;
  const PetscScalar *x,*lambda,*lambdae,*lambdai;
  PetscInt       loc,gloc=0;
  PetscScalar    Gff,Bff,Gft,Bft,Gtf,Btf,Gtt,Btt;
  PetscScalar    Vmf,Vmt,thetaf,thetat,thetaft,thetatf;
  PetscScalar    Pf,Qf,Pt,Qt;
  PSBUS          busf,bust;
  const PSBUS    *connbuses;
  PetscInt       xlocf,xloct;

  PetscFunctionBegin;

  ierr = OPFLOWGetSolution(opflow,&X);CHKERRQ(ierr);
  ierr = OPFLOWGetConstraintMultipliers(opflow,&Lambda);CHKERRQ(ierr);

  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecGetArrayRead(Lambda,&lambda);CHKERRQ(ierr);
  lambdae = lambda;
  if(opflow->Nconineq) {
    lambdai = lambdae + opflow->nconeq;
  }

  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    ierr = PSBUSGetVariableLocation(bus,&loc);

    bus->va = x[loc];
    bus->vm = x[loc+1];

    bus->mult_pmis = lambdae[gloc];
    bus->mult_qmis = lambdae[gloc+1];
    gloc += 2;

    for(k=0; k < bus->ngen; k++) {
      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      if(!gen->status) continue;

      loc += 2;

      gen->pg = x[loc];
      gen->qg = x[loc+1];
    }

    if(opflow->include_loadloss_variables) {
      for(k=0; k < bus->nload; k++) {
	ierr = PSBUSGetLoad(bus,k,&load);CHKERRQ(ierr);
	loc += 2;
	load->pl = load->pl - x[loc];
	load->ql = load->ql - x[loc+1];
      }
    }

    if(opflow->include_powerimbalance_variables) {
      loc += 2;
      bus->pimb = x[loc];
      bus->qimb = x[loc+1];
    }
  }

  gloc = 0;

  if(!opflow->ignore_lineflow_constraints) {
    for(i=0; i<ps->nline; i++) {
      line = &ps->line[i];
      if(!line->status) {
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
      
      ierr = PSLINEGetConnectedBuses(line,&connbuses);CHKERRQ(ierr);
      busf = connbuses[0];
      bust = connbuses[1];
      
      ierr = PSBUSGetVariableLocation(busf,&xlocf);CHKERRQ(ierr);
      ierr = PSBUSGetVariableLocation(bust,&xloct);CHKERRQ(ierr);
      
      thetaf  = x[xlocf];
      Vmf     = x[xlocf+1];
      thetat  = x[xloct];
      Vmt     = x[xloct+1];
      thetaft = thetaf - thetat;
      thetatf = thetat - thetaf;
      
      Pf = Gff*Vmf*Vmf  + Vmf*Vmt*(Gft*cos(thetaft) + Bft*sin(thetaft));
      Qf = -Bff*Vmf*Vmf + Vmf*Vmt*(-Bft*cos(thetaft) + Gft*sin(thetaft));
      
      Pt = Gtt*Vmt*Vmt  + Vmt*Vmf*(Gtf*cos(thetatf) + Btf*sin(thetatf));
      Qt = -Btt*Vmt*Vmt + Vmt*Vmf*(-Btf*cos(thetatf) + Gtf*sin(thetatf));

      line->pf = Pf;
      line->qf = Qf;
      line->pt = Pt;
      line->qt = Qt;
      line->sf = PetscSqrtScalar(Pf*Pf + Qf*Qf);
      line->st = PetscSqrtScalar(Pt*Pt + Qt*Qt);

      if(line->rateA > 1e5) {
	line->mult_sf = line->mult_st = 0.0;
      } else {
	line->mult_sf = lambdai[gloc];
	line->mult_st = lambdai[gloc+1];
	gloc += 2;
      }
    }
  }

  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(Lambda,&lambda);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWModelSetUp_PBPOLHIOP(OPFLOW opflow)
{
  PetscErrorCode ierr;
  PBPOLHIOP      pbpolhiop=(PBPOLHIOP)opflow->model;
  PetscInt       *idxn2sd_map;
  PS             ps;
  PSBUS          bus;
  int            ngen,nxsparse,nxdense;
  PSGEN          gen;
  
  PetscFunctionBegin;

  /* Create natural to sparse dense variable mapping */
  
  idxn2sd_map = opflow->idxn2sd_map;
  
  ps = opflow->ps;
  nxsparse = 2*ps->ngenON;
  nxdense  = 2*ps->nbus;
  
  int i,k;
  int spct=0,dnct=0;
  int loc;
  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    PSBUSGetVariableLocation(bus,&loc);

    idxn2sd_map[loc] = nxsparse + dnct;
    idxn2sd_map[loc+1] = nxsparse + dnct+1;

    dnct += 2;
    loc += 2;
    PSBUSGetNGen(bus,&ngen);
    for(k=0; k < ngen; k++) {
      PSBUSGetGen(bus,k,&gen);
      if(!gen->status) continue;

      idxn2sd_map[loc] = spct;
      idxn2sd_map[loc+1] = spct + 1;

      spct += 2;
      loc += 2;
    }
  }

  ierr = CreateBusParams(opflow,&pbpolhiop->busparams);
  ierr = CreateGenParams(opflow,&pbpolhiop->genparams);
  ierr = CreateLineParams(opflow,&pbpolhiop->lineparams);
  ierr = CreateLoadParams(opflow,&pbpolhiop->loadparams);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWModelDestroy_PBPOLHIOP(OPFLOW opflow)
{
  PetscErrorCode ierr;
  PBPOLHIOP pbpolhiop=(PBPOLHIOP)opflow->model;

  PetscFunctionBegin;

  ierr = DestroyBusParams(&pbpolhiop->busparams);
  ierr = DestroyGenParams(&pbpolhiop->genparams);
  ierr = DestroyLineParams(&pbpolhiop->lineparams);
  ierr = DestroyLoadParams(&pbpolhiop->loadparams);

  ierr = PetscFree(opflow->model);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWModelCreate_PBPOLHIOP(OPFLOW opflow)
{
  PBPOLHIOP pbpol;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  
  ierr = PetscCalloc1(1,&pbpol);CHKERRQ(ierr);

  opflow->model = pbpol;

  opflow->spdnordering = PETSC_TRUE;

  /* Inherit Ops */
  opflow->modelops.destroy                              = OPFLOWModelDestroy_PBPOLHIOP;
  opflow->modelops.setnumvariables                      = OPFLOWModelSetNumVariables_PBPOLHIOP;
  opflow->modelops.setnumconstraints                    = OPFLOWModelSetNumConstraints_PBPOLHIOP;
  opflow->modelops.setvariablebounds                    = OPFLOWSetVariableBounds_PBPOLHIOP;
  opflow->modelops.setvariableboundsarray               = OPFLOWSetVariableBoundsArray_PBPOLHIOP;
  opflow->modelops.setconstraintbounds                  = OPFLOWSetConstraintBounds_PBPOLHIOP;
  opflow->modelops.setconstraintboundsarray             = OPFLOWSetConstraintBoundsArray_PBPOLHIOP;
  opflow->modelops.setvariableandconstraintbounds       = OPFLOWSetVariableandConstraintBounds_PBPOLHIOP;
  opflow->modelops.setinitialguess                      = OPFLOWSetInitialGuess_PBPOLHIOP;
  opflow->modelops.setinitialguessarray                 = OPFLOWSetInitialGuessArray_PBPOLHIOP;
  opflow->modelops.computeequalityconstraints           = OPFLOWComputeEqualityConstraints_PBPOLHIOP;
  opflow->modelops.computeinequalityconstraints         = OPFLOWComputeInequalityConstraints_PBPOLHIOP;
  opflow->modelops.computeequalityconstraintsarray      = OPFLOWComputeEqualityConstraintsArray_PBPOLHIOP;
  opflow->modelops.computeinequalityconstraintsarray    = OPFLOWComputeInequalityConstraintsArray_PBPOLHIOP;
  opflow->modelops.computeequalityconstraintjacobian    = OPFLOWComputeEqualityConstraintJacobian_PBPOLHIOP;
  opflow->modelops.computeinequalityconstraintjacobian  = OPFLOWComputeInequalityConstraintJacobian_PBPOLHIOP;
  opflow->modelops.computehessian                       = OPFLOWComputeHessian_PBPOLHIOP;
  opflow->modelops.computeobjandgradient                = OPFLOWComputeObjandGradient_PBPOLHIOP;
  opflow->modelops.computeobjective                     = OPFLOWComputeObjective_PBPOLHIOP;
  opflow->modelops.computegradient                      = OPFLOWComputeGradient_PBPOLHIOP;
  opflow->modelops.computeobjectivearray                = OPFLOWComputeObjectiveArray_PBPOLHIOP;
  opflow->modelops.computegradientarray                 = OPFLOWComputeGradientArray_PBPOLHIOP;
  opflow->modelops.solutiontops                         = OPFLOWSolutionToPS_PBPOLHIOP;
  opflow->modelops.setup                                = OPFLOWModelSetUp_PBPOLHIOP;
  opflow->modelops.computesparsejacobianhiop            = OPFLOWComputeSparseJacobian_PBPOLHIOP;
  opflow->modelops.computesparsehessianhiop             = OPFLOWComputeSparseHessian_PBPOLHIOP;
  opflow->modelops.computedenseequalityconstraintjacobianhiop = OPFLOWComputeDenseEqualityConstraintJacobian_PBPOLHIOP;
  opflow->modelops.computedenseinequalityconstraintjacobianhiop = OPFLOWComputeDenseInequalityConstraintJacobian_PBPOLHIOP;
  opflow->modelops.computedensehessianhiop    = OPFLOWComputeDenseHessian_PBPOLHIOP;
  
  PetscFunctionReturn(0);
}

#endif
