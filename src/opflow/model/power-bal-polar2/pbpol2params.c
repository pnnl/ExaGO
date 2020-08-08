
#include <private/psimpl.h>
#include <private/opflowimpl.h>
#include "pbpol2.h"

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

    if(opflow->genbusVmfixed) {
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
  ierr = PetscCalloc1(lineparams->nlineON,&lineparams->gbineqidx);CHKERRQ(ierr);

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

    if(line->rateA < 1e5) {
      lineparams->gbineqidx[linelimi] = gbloc;
      lineparams->gineqidx[linelimi] = gloc;
      lineparams->linelimidx[linelimi] = linei;
      linelimi++;
      gbloc += 2;
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
