#include <scopflow_config.h>

#if defined(EXAGO_HAVE_HIOP) 

#include <private/psimpl.h>
#include <private/opflowimpl.h>
#include "pbpolhiop.h"

/************* NOTE ***********************/
/* No Load loss or power imbalance variables considered yet */
/********************************************/

PetscErrorCode OPFLOWSetConstraintBoundsArray_PBPOLHIOP(OPFLOW opflow,double *gl,double *gu)
{
  PBPOLHIOP      pbpolhiop=(PBPOLHIOP)opflow->model;
  BUSParams      *busparams=&pbpolhiop->busparams;
  LINEParams     *lineparams=&pbpolhiop->lineparams;
  PetscInt       i;
  PS             ps=opflow->ps;

  PetscFunctionBegin;

  /* Equallity constraints (all zeros) */
  for(i=0; i < busparams->nbus; i++) {
    gl[busparams->gidx[i]] = 0.0;
    gu[busparams->gidx[i]] = 0.0;

    gl[busparams->gidx[i]+1] = 0.0;
    gu[busparams->gidx[i]+1] = 0.0;
  }

  /* Inequality constraints */
  for(i=0; i < lineparams->nlinelim; i++) {
    int    j=lineparams->linelimidx[i];
    /* Not atomic */
    gl[lineparams->gbineqidx[i]]   = 0.0;
    gu[lineparams->gbineqidx[i]]   = (lineparams->rateA[j]/ps->MVAbase)*(lineparams->rateA[j]/ps->MVAbase);
    gl[lineparams->gbineqidx[i]+1] = 0.0;
    gu[lineparams->gbineqidx[i]+1] = (lineparams->rateA[j]/ps->MVAbase)*(lineparams->rateA[j]/ps->MVAbase);
  }

  PetscFunctionReturn(0);
}

/* The calculations for different routines start from here */

/** CONSTRAINT BOUNDS  **/
PetscErrorCode OPFLOWSetConstraintBounds_PBPOLHIOP(OPFLOW opflow,Vec Gl, Vec Gu)
{
  PetscErrorCode ierr;
  PetscScalar    *gl,*gu;

  PetscFunctionBegin;

  ierr = VecSet(Gl,0.0);
  ierr = VecSet(Gu,0.0);

  ierr = VecGetArray(Gl,&gl);CHKERRQ(ierr);
  ierr = VecGetArray(Gu,&gu);CHKERRQ(ierr);

  ierr = OPFLOWSetConstraintBoundsArray_PBPOLHIOP(opflow,gl,gu);CHKERRQ(ierr);

  ierr = VecRestoreArray(Gl,&gl);CHKERRQ(ierr);
  ierr = VecRestoreArray(Gu,&gu);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/** EQUALITY CONSTRAINTS */
PetscErrorCode OPFLOWComputeEqualityConstraintsArray_PBPOLHIOP(OPFLOW opflow,const double *x, double *ge)
{
  PBPOLHIOP       pbpolhiop=(PBPOLHIOP)opflow->model;
  BUSParams      *busparams=&pbpolhiop->busparams;
  GENParams      *genparams=&pbpolhiop->genparams;
  LOADParams     *loadparams=&pbpolhiop->loadparams;
  LINEParams     *lineparams=&pbpolhiop->lineparams;
  PetscInt       i;

  PetscFunctionBegin;

  for(i=0; i < opflow->nconeq; i++) ge[i] = 0.0;

  /* Generator contributions */
  for(i=0; i < genparams->ngenON; i++) {
    /* atomic operation */
    ge[genparams->gidx[i]]   -= x[genparams->xidx[i]];
    ge[genparams->gidx[i]+1] -= x[genparams->xidx[i]+1];
  }

  /* Load contributions */
  for(i=0; i < loadparams->nload; i++) {
    /* Atomic operation */
    ge[loadparams->gidx[i]]   += loadparams->pl[i];
    ge[loadparams->gidx[i]+1] += loadparams->ql[i];
  }

  /* Bus contributions */
  for(i=0; i < busparams->nbus; i++) {
    double theta= x[busparams->xidx[i]];
    double Vm   = x[busparams->xidx[i]+1];
    /* Atomic operation if launched in parallel with other component class kernels */
    ge[busparams->gidx[i]]   += busparams->isisolated[i]*(theta - busparams->va[i]*PETSC_PI/180.0) + busparams->ispvpq[i]*Vm*Vm*busparams->gl[i];
    ge[busparams->gidx[i]+1] += busparams->isisolated[i]*(Vm    - busparams->vm[i]) - busparams->ispvpq[i]*Vm*Vm*busparams->bl[i];
  }

  /* Line contributions */
  for(i=0; i < lineparams->nlineON; i++) {
    double Pf,Qf,Pt,Qt;
    double thetaf=x[lineparams->xidxf[i]], Vmf=x[lineparams->xidxf[i]+1];
    double thetat=x[lineparams->xidxt[i]], Vmt=x[lineparams->xidxt[i]+1];
    double thetaft=thetaf-thetat;
    double thetatf=thetat-thetaf;

    Pf = lineparams->Gff[i]*Vmf*Vmf  + Vmf*Vmt*(lineparams->Gft[i]*cos(thetaft) + lineparams->Bft[i]*sin(thetaft));
    Qf = -lineparams->Bff[i]*Vmf*Vmf + Vmf*Vmt*(-lineparams->Bft[i]*cos(thetaft) + lineparams->Gft[i]*sin(thetaft));
    Pt = lineparams->Gtt[i]*Vmt*Vmt  + Vmt*Vmf*(lineparams->Gtf[i]*cos(thetatf) + lineparams->Btf[i]*sin(thetatf));
    Qt = -lineparams->Btt[i]*Vmt*Vmt + Vmt*Vmf*(-lineparams->Btf[i]*cos(thetatf) + lineparams->Gtf[i]*sin(thetatf));
    
    /* Atomic operation */
    ge[lineparams->geqidxf[i]]   += Pf;
    ge[lineparams->geqidxf[i]+1] += Qf;
    ge[lineparams->geqidxt[i]]   += Pt;
    ge[lineparams->geqidxt[i]+1] += Qt;
  }

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeEqualityConstraints_PBPOLHIOP(OPFLOW opflow, Vec X, Vec Ge)
{
  PetscErrorCode ierr;
  PetscScalar    *ge;
  const PetscScalar *x;

  PetscFunctionBegin;

  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecGetArray(Ge,&ge);CHKERRQ(ierr);

  ierr = OPFLOWComputeEqualityConstraintsArray_PBPOLHIOP(opflow,x,ge);CHKERRQ(ierr);

  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(Ge,&ge);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/** INEQUALITY CONSTRAINTS **/
PetscErrorCode OPFLOWComputeInequalityConstraintsArray_PBPOLHIOP(OPFLOW opflow, const double *x, double *gi)
{
  PBPOLHIOP         pbpolhiop=(PBPOLHIOP)opflow->model;
  LINEParams     *lineparams=&pbpolhiop->lineparams;
  PetscInt       i;

  PetscFunctionBegin;

  for(i=0; i < opflow->nconineq; i++) gi[i] = 0.0;

  /* Line contributions */
  for(i=0; i < lineparams->nlinelim; i++) {
    int    j=lineparams->linelimidx[i];
    double Pf,Qf,Pt,Qt,Sf2,St2;
    double thetaf=x[lineparams->xidxf[j]], Vmf=x[lineparams->xidxf[j]+1];
    double thetat=x[lineparams->xidxt[j]], Vmt=x[lineparams->xidxt[j]+1];
    double thetaft=thetaf-thetat;
    double thetatf=thetat-thetaf;
    
    Pf = lineparams->Gff[j]*Vmf*Vmf  + Vmf*Vmt*(lineparams->Gft[j]*cos(thetaft) + lineparams->Bft[j]*sin(thetaft));
    Qf = -lineparams->Bff[j]*Vmf*Vmf + Vmf*Vmt*(-lineparams->Bft[j]*cos(thetaft) + lineparams->Gft[j]*sin(thetaft));
    Pt = lineparams->Gtt[j]*Vmt*Vmt  + Vmt*Vmf*(lineparams->Gtf[j]*cos(thetatf) + lineparams->Btf[j]*sin(thetatf));
    Qt = -lineparams->Btt[j]*Vmt*Vmt + Vmt*Vmf*(-lineparams->Btf[j]*cos(thetatf) + lineparams->Gtf[j]*sin(thetatf));
    
    Sf2 = Pf*Pf + Qf*Qf;
    St2 = Pt*Pt + Qt*Qt;

    /* Not atomic */
    gi[lineparams->gineqidx[i]]   = Sf2;
    gi[lineparams->gineqidx[i]+1] = St2;
  }

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeInequalityConstraints_PBPOLHIOP(OPFLOW opflow, Vec X, Vec Gi)
{
  PetscErrorCode ierr;
  PetscScalar    *gi;
  const PetscScalar *x;

  PetscFunctionBegin;
  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecGetArray(Gi,&gi);CHKERRQ(ierr);

  ierr = OPFLOWComputeInequalityConstraintsArray_PBPOLHIOP(opflow,x,gi);CHKERRQ(ierr);

  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(Gi,&gi);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


/** OBJECTIVE FUNCTION **/
PetscErrorCode OPFLOWComputeObjectiveArray_PBPOLHIOP(OPFLOW opflow,const double *x,double *obj)
{
  PBPOLHIOP      pbpolhiop=(PBPOLHIOP)opflow->model;
  GENParams     *genparams=&pbpolhiop->genparams;
  PetscInt       i;
  PS             ps=opflow->ps;
  double         obj_val=0.0;
  int            isobj_gencost=opflow->obj_gencost;
  double         MVAbase=ps->MVAbase;

  PetscFunctionBegin;

  /* Generator objective function contributions */
  for(i=0; i < genparams->ngenON; i++) {
    double Pg = x[genparams->xidx[i]]*MVAbase;
    obj_val += isobj_gencost*(genparams->cost_alpha[i]*Pg*Pg + genparams->cost_beta[i]*Pg + genparams->cost_gamma[i]);
  }

  *obj = obj_val;

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeObjective_PBPOLHIOP(OPFLOW opflow,Vec X,double *obj)
{
  PetscErrorCode ierr;
  const PetscScalar *x;

  PetscFunctionBegin;
  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);

  ierr = OPFLOWComputeObjectiveArray_PBPOLHIOP(opflow,x,obj);CHKERRQ(ierr);

  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/** GRADIENT **/
PetscErrorCode OPFLOWComputeGradientArray_PBPOLHIOP(OPFLOW opflow,const double *x, double* grad)
{
  PBPOLHIOP         pbpolhiop=(PBPOLHIOP)opflow->model;
  GENParams     *genparams=&pbpolhiop->genparams;
  PetscInt       i;
  PS             ps=opflow->ps;
  int            isobj_gencost=opflow->obj_gencost;
  double         MVAbase=ps->MVAbase;

  PetscFunctionBegin;

  for(i=0; i < opflow->nx; i++) grad[i] = 0.0;

  /* Generator gradient contributions */
  for(i=0; i < genparams->ngenON; i++) {
    double Pg = x[genparams->xidx[i]]*MVAbase;
    grad[genparams->xidx[i]] = isobj_gencost*MVAbase*(2.0*genparams->cost_alpha[i]*Pg + genparams->cost_beta[i]);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeGradient_PBPOLHIOP(OPFLOW opflow,Vec X,Vec Grad)
{
  PetscErrorCode ierr;
  PetscScalar    *grad;
  const PetscScalar *x;

  PetscFunctionBegin;
  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecGetArray(Grad,&grad);CHKERRQ(ierr);

  ierr = OPFLOWComputeGradientArray_PBPOLHIOP(opflow,x,grad);CHKERRQ(ierr);

  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(Grad,&grad);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSetVariableBoundsArray_PBPOLHIOP(OPFLOW opflow,double *xl,double *xu)
{
  PBPOLHIOP      pbpolhiop=(PBPOLHIOP)opflow->model;
  BUSParams      *busparams=&pbpolhiop->busparams;
  GENParams      *genparams=&pbpolhiop->genparams;
  PetscInt       i;

  PetscFunctionBegin;

  /* Bounds for bus voltages */
  for(i=0; i < busparams->nbus; i++) {
    xl[busparams->xidx[i]] = busparams->ispvpq[i]*PETSC_NINFINITY + busparams->isisolated[i]*busparams->va[i] + busparams->isref[i]*busparams->va[i]*PETSC_PI/180.0;
    xu[busparams->xidx[i]] = busparams->ispvpq[i]*PETSC_INFINITY  + busparams->isisolated[i]*busparams->va[i] + busparams->isref[i]*busparams->va[i]*PETSC_PI/180.0;
    
    xl[busparams->xidx[i]+1] = busparams->isref[i]*busparams->vmin[i]  + busparams->ispvpq[i]*busparams->vmin[i] + busparams->isisolated[i]*busparams->vm[i];
    xu[busparams->xidx[i]+1] = busparams->isref[i]*busparams->vmax[i]  + busparams->ispvpq[i]*busparams->vmax[i] + busparams->isisolated[i]*busparams->vm[i];
  }

  /* Generator lower and upper bounds on variables */
  for(i=0; i < genparams->ngenON; i++) {
    xl[genparams->xidx[i]]   = genparams->pb[i];
    xu[genparams->xidx[i]]   = genparams->pt[i];
    xl[genparams->xidx[i]+1] = genparams->qb[i];
    xu[genparams->xidx[i]+1] = genparams->qt[i];
  }

  PetscFunctionReturn(0);
}


/** VARIABLE BOUNDS **/
PetscErrorCode OPFLOWSetVariableBounds_PBPOLHIOP(OPFLOW opflow,Vec Xl,Vec Xu)
{
  PetscErrorCode ierr;
  PetscScalar    *xl,*xu;

  PetscFunctionBegin;
  ierr = VecGetArray(Xl,&xl);CHKERRQ(ierr);
  ierr = VecGetArray(Xu,&xu);CHKERRQ(ierr);

  ierr = OPFLOWSetVariableBoundsArray_PBPOLHIOP(opflow,xl,xu);CHKERRQ(ierr);

  ierr = VecRestoreArray(Xl,&xl);CHKERRQ(ierr);
  ierr = VecRestoreArray(Xu,&xu);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeEqualityConstraintJacobian_PBPOLHIOP(OPFLOW opflow,Vec X,Mat Je)
{
  PetscErrorCode ierr;
  PetscInt       i,row[2],col[4];
  PBPOLHIOP         pbpolhiop=(PBPOLHIOP)opflow->model;
  BUSParams      *busparams=&pbpolhiop->busparams;
  GENParams      *genparams=&pbpolhiop->genparams;
  LINEParams     *lineparams=&pbpolhiop->lineparams;
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

PetscErrorCode OPFLOWComputeInequalityConstraintJacobian_PBPOLHIOP(OPFLOW opflow,Vec X,Mat Ji)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PBPOLHIOP         pbpolhiop=(PBPOLHIOP)opflow->model;
  LINEParams     *lineparams=&pbpolhiop->lineparams;
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
PetscErrorCode OPFLOWComputeObjectiveHessian_PBPOLHIOP(OPFLOW opflow,Vec X,Mat H) 
{
  PetscErrorCode ierr;
  PBPOLHIOP         pbpolhiop=(PBPOLHIOP)opflow->model;
  GENParams     *genparams=&pbpolhiop->genparams;
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
PetscErrorCode OPFLOWComputeEqualityConstraintsHessian_PBPOLHIOP(OPFLOW opflow,Vec X,Vec Lambda,Mat H) 
{
  PetscErrorCode ierr;
  PBPOLHIOP         pbpolhiop=(PBPOLHIOP)opflow->model;
  BUSParams      *busparams=&pbpolhiop->busparams;
  LINEParams     *lineparams=&pbpolhiop->lineparams;
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
PetscErrorCode OPFLOWComputeInequalityConstraintsHessian_PBPOLHIOP(OPFLOW opflow, Vec X, Vec Lambda,Mat H)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PBPOLHIOP         pbpolhiop=(PBPOLHIOP)opflow->model;
  LINEParams     *lineparams=&pbpolhiop->lineparams;
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
    row[0] = lineparams->xidxt[j] + 1;
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

/** Custom routines that work with HIOP interface only */
PetscErrorCode OPFLOWSetSparseJacobianLocations_PBPOLHIOP(OPFLOW opflow,int *iJacS, int *jJacS,int *nnz)
{
  PBPOLHIOP      pbpolhiop=(PBPOLHIOP)opflow->model;
  GENParams     *genparams=&pbpolhiop->genparams;
  int            i;

  /* Generator contributions */
  for(i=0; i < genparams->ngenON; i++) {
    iJacS[genparams->jacsp_idx[i]] = genparams->gidx[i];
    jJacS[genparams->jacsp_idx[i]] = genparams->xidx[i];

    iJacS[genparams->jacsq_idx[i]] = genparams->gidx[i]+1;
    jJacS[genparams->jacsq_idx[i]] = genparams->xidx[i]+1;

    *nnz = genparams->jacsq_idx[i]+1;
  }

  PetscFunctionReturn(0);
}
    


#endif
