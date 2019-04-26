#include <private/psimpl.h>
#include <private/opflowimpl.h>
#include "IpStdCInterface.h"

/* Function Declarations */
Bool eval_opflow_f(PetscInt n, PetscScalar* x, Bool new_x,
            PetscScalar* obj_value, UserDataPtr user_data);

Bool eval_opflow_grad_f(PetscInt n, PetscScalar* x, Bool new_x,
                 PetscScalar* grad_f, UserDataPtr user_data);

Bool eval_opflow_g(PetscInt n, PetscScalar* x, Bool new_x,
            PetscInt m, PetscScalar* g, UserDataPtr user_data);

Bool eval_opflow_jac_g(PetscInt n, PetscScalar *x, Bool new_x,
                PetscInt m, PetscInt nele_jac,
                PetscInt *iRow, PetscInt *jCol, PetscScalar *values,
                UserDataPtr user_data);

Bool eval_opflow_h(PetscInt n, PetscScalar *x, Bool new_x, PetscScalar obj_factor,
            PetscInt m, PetscScalar *lambda, Bool new_lambda,
            PetscInt nele_hess, PetscInt *iRow, PetscInt *jCol,
            PetscScalar *values, UserDataPtr user_data);

/*
  OPFLOWGetConstraintJacobianNonzeros - Gets the number of nonzeros in the lagrangian hessian matrix

  Input Paramereters:
. opflow - the optimal power flow application object

  Output Parameters:
. nnz - number of nonzeros in the lagrangian hessian

  Notes:
  OPFLOWSetUp() must be called before calling this routine.
*/
PetscErrorCode OPFLOWGetLagrangianHessianNonzeros(OPFLOW opflow,PetscInt *nnz)
{
  PetscErrorCode ierr;
  PS             ps=opflow->ps;
  PetscInt       i,k;
  PSBUS          bus;
  PSLINE         line;
  PSGEN          gen;  
  PetscInt       nconnlines;
  const PSLINE   *connlines;

  PetscFunctionBegin;

  *nnz = 0;

  // for objective
  for(i=0; i < ps->nbus; i++){
    bus = &ps->bus[i];
    for(k=0; k < bus->ngen; k++) {
	  ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      if(!gen->status) continue;
	  *nnz += 1;
    }
  }

  for(i=0; i < ps->Nbranch; i++) *nnz += 10;

  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    *nnz += 1;
    ierr = PSBUSGetSupportingLines(bus,&nconnlines,&connlines);CHKERRQ(ierr);
    for(k=0; k < nconnlines;k++) {
      line = connlines[k];
      if(!line->status) continue;
      *nnz += 9;
    }
  }

  PetscFunctionReturn(0);
}



/*
  OPFLOWSetConstraintJacobianLocations - Sets the row and column nonzero locations for the
              lagrangian hessian

  Input Parameters:
+ opflow - the OPFLOW object

  Output Parameters:
+ row - the row locations
- col - the col locations

  Notes:
   The row and column locations should be such that H(row[i],col[i]) = val
*/
PetscErrorCode OPFLOWSetLagrangianHessianLocations(OPFLOW opflow, PetscInt *row, PetscInt *col)
{
  PetscErrorCode ierr;
  PS             ps=opflow->ps;
  PetscInt       ctr=0;
  PetscInt       i,k;
  PSLINE         line;
  PSBUS          bus;
  PSGEN          gen;
  PetscInt       nconnlines;
  const PSLINE   *connlines;
  const PSBUS    *connbuses;
  PetscInt       gloc=0,xloc,xlocf,xloct;
  PSBUS          busf,bust;

  PetscFunctionBegin;

  // for the part of objective
  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    ierr = PSBUSGetVariableLocation(bus,&xloc);CHKERRQ(ierr);
   
    for(k=0; k < bus->ngen; k++) {
      xloc = xloc+2;
      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      if(!gen->status) continue;
      row[ctr] = xloc;
      col[ctr] = xloc;
      ctr += 1;
    }
  }

  // for the part of line constraints
  for(i=0; i < ps->Nbranch; i++) {
    line = &ps->line[i];

    ierr = PSLINEGetConnectedBuses(line,&connbuses);CHKERRQ(ierr);
    busf = connbuses[0];
    bust = connbuses[1];

    ierr = PSBUSGetVariableLocation(busf,&xlocf);CHKERRQ(ierr);
    ierr = PSBUSGetVariableLocation(bust,&xloct);CHKERRQ(ierr);

    row[ctr] 	= xlocf; 	col[ctr] 	= xlocf;
    row[ctr+1] 	= xlocf+1; 	col[ctr+1] 	= xlocf;
    row[ctr+2] 	= xlocf+1; 	col[ctr+2] 	= xlocf+1;
    row[ctr+3] 	= xloct; 	col[ctr+3] 	= xlocf;
    row[ctr+4] 	= xloct; 	col[ctr+4] 	= xlocf+1;
    row[ctr+5] 	= xloct; 	col[ctr+5] 	= xloct;
    row[ctr+6] 	= xloct+1; 	col[ctr+6] 	= xlocf;
    row[ctr+7] 	= xloct+1; 	col[ctr+7] 	= xlocf+1;
    row[ctr+8] 	= xloct+1; 	col[ctr+8] 	= xloct;
    row[ctr+9] 	= xloct+1; 	col[ctr+9] 	= xloct+1;

    ctr += 10;
  }

  // for the part of bus constraints
  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    ierr = PSBUSGetVariableLocation(bus,&xloc);CHKERRQ(ierr);

    row[ctr] 	= xloc+1;   col[ctr]   = xloc+1;
    ctr += 1;

    ierr = PSBUSGetSupportingLines(bus,&nconnlines,&connlines);CHKERRQ(ierr);

    for(k=0; k < nconnlines; k++) {
      line = connlines[k];
     
      if(!line->status) continue;

      ierr = PSLINEGetConnectedBuses(line,&connbuses);CHKERRQ(ierr);
      busf = connbuses[0];
      bust = connbuses[1];
      
      ierr = PSBUSGetVariableLocation(busf,&xlocf);CHKERRQ(ierr);
      ierr = PSBUSGetVariableLocation(bust,&xloct);CHKERRQ(ierr);

      if(bus == busf) {
	   	row[ctr] 	= xlocf; 	col[ctr] 	= xlocf;
		row[ctr+1] 	= xlocf+1; 	col[ctr+1] 	= xlocf;
	    row[ctr+2] 	= xlocf+1; 	col[ctr+2] 	= xlocf+1;
	    row[ctr+3] 	= xloct; 	col[ctr+3] 	= xlocf;
	    row[ctr+4] 	= xloct; 	col[ctr+4] 	= xlocf+1;
	    row[ctr+5] 	= xloct; 	col[ctr+5] 	= xloct;
	    row[ctr+6] 	= xloct+1; 	col[ctr+6] 	= xlocf;
	    row[ctr+7] 	= xloct+1; 	col[ctr+7] 	= xlocf+1;
	    row[ctr+8] 	= xloct+1; 	col[ctr+8] 	= xloct;
		
	 	ctr += 9;
      } else {
		row[ctr]	= xlocf;	col[ctr]	= xlocf;
		row[ctr+1]  = xlocf+1;  col[ctr+1]  = xlocf;
		row[ctr+2]  = xloct;	col[ctr+2]  = xlocf;
		row[ctr+3]  = xloct;	col[ctr+3]  = xlocf+1;
		row[ctr+4]  = xloct;	col[ctr+4]  = xloct;
		row[ctr+5]  = xloct+1;  col[ctr+5]  = xlocf;
		row[ctr+6]  = xloct+1;  col[ctr+6]  = xlocf+1;
		row[ctr+7]  = xloct+1;  col[ctr+7]  = xloct;
		row[ctr+8]  = xloct+1;  col[ctr+8]  = xloct+1;
	
		ctr += 9;
      }
	}
  }

  PetscFunctionReturn(0);
}



/*
  OPFLOWSetHessianValues - Sets the nonzero values for the
              Lagrangian Hessian

  Input Parameters:
+ opflow - the OPFLOW object
- X      - the current iterate


  Output Parameters:
. values - the nonzero values in the Lagrangian Hessian

  Notes:
   The values should be in the same row and col order as set in OPFLOWSetLagrangianHessianLocations
*/
PetscErrorCode OPFLOWSetLagrangianHessianValues(OPFLOW opflow, PetscScalar obj_factor, Vec X, Vec Lambda, PetscScalar *values)
{
  PetscErrorCode ierr;
  const PetscScalar *x;
  const PetscScalar *lambda;
  
  PS             ps=opflow->ps;
  PetscInt       ctr=0;
  PetscInt       i,k;
  PSLINE         line;
  PSBUS          bus;
  PSGEN          gen;
  PetscInt       nconnlines;
  const PSLINE   *connlines;
  const PSBUS    *connbuses;
  PetscInt       gloc=0,xloc,xlocf,xloct;
  PSBUS          busf,bust;

  PetscFunctionBegin;

  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecGetArrayRead(Lambda,&lambda);CHKERRQ(ierr);

  // for the part of objective
  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    ierr = PSBUSGetVariableLocation(bus,&xloc);CHKERRQ(ierr);
   
    for(k=0; k < bus->ngen; k++) {
      xloc = xloc+2;
      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      if(!gen->status) continue;

	  values[ctr]	= obj_factor*2.0*gen->cost_alpha*ps->MVAbase*ps->MVAbase;
	  ctr += 1;	  
    }
  }

  // for the part of line constraints
  for(i=0; i < ps->Nbranch; i++) {
    line = &ps->line[i];

    PetscScalar Gff,Bff,Gft,Bft,Gtf,Btf,Gtt,Btt;
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

    PetscScalar Vmf,Vmt,thetaf,thetat,thetaft,thetatf;
      
    thetaf  = x[xlocf];
    Vmf     = x[xlocf+1];
    thetat  = x[xloct];
    Vmt     = x[xloct+1];
    thetaft = thetaf - thetat;
    thetatf = thetat - thetaf;

	// Sf2 and St2 are the constraints	
    PetscScalar Pf,Qf,Pt,Qt,Sf2,St2;

    Pf =  Gff*Vmf*Vmf + Vmf*Vmt*( Gft*cos(thetaft) + Bft*sin(thetaft));
    Qf = -Bff*Vmf*Vmf + Vmf*Vmt*(-Bft*cos(thetaft) + Gft*sin(thetaft));
	
    Pt =  Gtt*Vmt*Vmt + Vmt*Vmf*( Gtf*cos(thetatf) + Btf*sin(thetatf));
    Qt = -Btt*Vmt*Vmt + Vmt*Vmf*(-Btf*cos(thetatf) + Gtf*sin(thetatf));

    Sf2 = Pf*Pf + Qf*Qf;
    St2 = Pt*Pt + Qt*Qt;

    PetscScalar dSf2_dPf, dSf2_dQf, dSt2_dPt, dSt2_dQt;

    dSf2_dPf = 2.*Pf;
    dSf2_dQf = 2.*Qf;
    dSt2_dPt = 2.*Pt;
    dSt2_dQt = 2.*Qt;

    PetscScalar dSf2_dPf_dPf, dSt2_dPt_dPt;
    PetscScalar dSf2_dQf_dQf, dSt2_dQt_dQt;

    dSf2_dPf_dPf = 2.;
    dSf2_dQf_dQf = 2.;	
    dSt2_dPt_dPt = 2.;
    dSt2_dQt_dQt = 2.;
	

    PetscScalar dPf_dthetaf,dPf_dVmf,dPf_dthetat,dPf_dVmt;
    PetscScalar dQf_dthetaf,dQf_dVmf,dQf_dthetat,dQf_dVmt;
    PetscScalar dPt_dthetaf,dPt_dVmf,dPt_dthetat,dPt_dVmt;
    PetscScalar dQt_dthetaf,dQt_dVmf,dQt_dthetat,dQt_dVmt;

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

    PetscScalar dPf_dthetaf_dthetaf,dPf_dVmf_dthetaf,dPf_dthetat_dthetaf,dPf_dVmt_dthetaf;
    PetscScalar dPf_dVmf_dVmf,dPf_dthetat_dVmf,dPf_dVmt_dVmf;
    PetscScalar dPf_dthetat_dthetat,dPf_dVmt_dthetat;
    PetscScalar dPf_dVmt_dVmt;

    dPf_dthetaf_dthetaf = 			Vmf*Vmt*(-Gft*cos(thetaft) - Bft*sin(thetaft));
    dPf_dVmf_dthetaf    =				Vmt*(-Gft*sin(thetaft) + Bft*cos(thetaft));
    dPf_dthetat_dthetaf = 			Vmf*Vmt*( Gft*cos(thetaft) + Bft*sin(thetaft));
	dPf_dVmt_dthetaf    = 				Vmf*(-Gft*sin(thetaft) + Bft*cos(thetaft));

    dPf_dVmf_dVmf    	= 2.*Gff;
    dPf_dthetat_dVmf 	= 				Vmt*( Gft*sin(thetaft) - Bft*cos(thetaft));
    dPf_dVmt_dVmf    	= 					( Gft*cos(thetaft) + Bft*sin(thetaft));

    dPf_dthetat_dthetat = 			Vmf*Vmt*(-Gft*cos(thetaft) - Bft*sin(thetaft));
    dPf_dVmt_dthetat    = 				Vmf*( Gft*sin(thetaft) - Bft*cos(thetaft));

    dPf_dVmt_dVmt    	= 				0.;

    PetscScalar dQf_dthetaf_dthetaf,dQf_dVmf_dthetaf,dQf_dthetat_dthetaf,dQf_dVmt_dthetaf;
    PetscScalar dQf_dVmf_dVmf,dQf_dthetat_dVmf,dQf_dVmt_dVmf;
    PetscScalar dQf_dthetat_dthetat,dQf_dVmt_dthetat;
    PetscScalar dQf_dVmt_dVmt;

    dQf_dthetaf_dthetaf = 			Vmf*Vmt*( Bft*cos(thetaft) - Gft*sin(thetaft));
    dQf_dVmf_dthetaf    = 				Vmt*( Bft*sin(thetaft) + Gft*cos(thetaft));
    dQf_dthetat_dthetaf = 			Vmf*Vmt*(-Bft*cos(thetaft) + Gft*sin(thetaft));
    dQf_dVmt_dthetaf    = 				Vmf*( Bft*sin(thetaft) + Gft*cos(thetaft));

    dQf_dVmf_dVmf    	= -2.*Bff;
    dQf_dthetat_dVmf 	= 				Vmt*(-Bft*sin(thetaft) - Gft*cos(thetaft));
    dQf_dVmt_dVmf    	= 					(-Bft*cos(thetaft) + Gft*sin(thetaft));

    dQf_dthetat_dthetat = 			Vmf*Vmt*( Bft*cos(thetaft) - Gft*sin(thetaft));
    dQf_dVmt_dthetat    = 				Vmf*(-Bft*sin(thetaft) - Gft*cos(thetaft));

    dQf_dVmt_dVmt    	= 				0.;
	
    PetscScalar dPt_dthetat_dthetat,dPt_dVmt_dthetat,dPt_dthetaf_dthetat,dPt_dVmf_dthetat;
    PetscScalar dPt_dVmt_dVmt,dPt_dthetaf_dVmt,dPt_dVmf_dVmt;
    PetscScalar dPt_dthetaf_dthetaf,dPt_dVmf_dthetaf;
    PetscScalar dPt_dVmf_dVmf;
	
    dPt_dthetat_dthetat = 			Vmt*Vmf*(-Gtf*cos(thetatf) - Btf*sin(thetatf));
    dPt_dVmt_dthetat    = 				Vmf*(-Gtf*sin(thetatf) + Btf*cos(thetatf));
    dPt_dthetaf_dthetat = 			Vmt*Vmf*( Gtf*cos(thetatf) + Btf*sin(thetatf));
    dPt_dVmf_dthetat    = 				Vmt*(-Gtf*sin(thetatf) + Btf*cos(thetatf));

    dPt_dVmt_dVmt   	= 2.*Gtt;
    dPt_dthetaf_dVmt 	= 				Vmf*( Gtf*sin(thetatf) - Btf*cos(thetatf));
    dPt_dVmf_dVmt    	= 					( Gtf*cos(thetatf) + Btf*sin(thetatf));

    dPt_dthetaf_dthetaf = 			Vmt*Vmf*(-Gtf*cos(thetatf) - Btf*sin(thetatf));
    dPt_dVmf_dthetaf    = 				Vmt*( Gtf*sin(thetatf) - Btf*cos(thetatf));

    dPt_dVmf_dVmf    	= 				0.;

    PetscScalar dQt_dthetat_dthetat,dQt_dVmt_dthetat,dQt_dthetaf_dthetat,dQt_dVmf_dthetat;
    PetscScalar dQt_dVmt_dVmt,dQt_dthetaf_dVmt,dQt_dVmf_dVmt;
    PetscScalar dQt_dthetaf_dthetaf,dQt_dVmf_dthetaf;
    PetscScalar dQt_dVmf_dVmf;	

    dQt_dthetat_dthetat = 			Vmt*Vmf*( Btf*cos(thetatf) - Gtf*sin(thetatf));
    dQt_dVmt_dthetat    = 				Vmf*( Btf*sin(thetatf) + Gtf*cos(thetatf));
    dQt_dthetaf_dthetat = 			Vmt*Vmf*(-Btf*cos(thetatf) + Gtf*sin(thetatf));
    dQt_dVmf_dthetat    = 				Vmt*( Btf*sin(thetatf) + Gtf*cos(thetatf));

    dQt_dVmt_dVmt    	= -2.*Btt;
    dQt_dthetaf_dVmt 	= 				Vmf*(-Btf*sin(thetatf) - Gtf*cos(thetatf));
    dQt_dVmf_dVmt    	= 					(-Btf*cos(thetatf) + Gtf*sin(thetatf));

    dQt_dthetaf_dthetaf = 			Vmt*Vmf*( Btf*cos(thetatf) - Gtf*sin(thetatf));
    dQt_dVmf_dthetaf    = 				Vmt*(-Btf*sin(thetatf) - Gtf*cos(thetatf));

    dQt_dVmf_dVmf    	= 				0.;


    PetscScalar dSf2_dPf_dthetaf, dSf2_dQf_dthetaf, dSt2_dPt_dthetaf, dSt2_dQt_dthetaf;
    dSf2_dPf_dthetaf = dSf2_dPf_dPf*dPf_dthetaf;
    dSf2_dQf_dthetaf = dSf2_dQf_dQf*dQf_dthetaf;
    dSt2_dPt_dthetaf = dSt2_dPt_dPt*dPt_dthetaf;
    dSt2_dQt_dthetaf = dSt2_dQt_dQt*dQt_dthetaf;

    PetscScalar dSf2_dPf_dVmf, dSf2_dQf_dVmf, dSt2_dPt_dVmf, dSt2_dQt_dVmf;
    dSf2_dPf_dVmf = dSf2_dPf_dPf*dPf_dVmf;
    dSf2_dQf_dVmf = dSf2_dQf_dQf*dQf_dVmf;
    dSt2_dPt_dVmf = dSt2_dPt_dPt*dPt_dVmf;
    dSt2_dQt_dVmf = dSt2_dQt_dQt*dQt_dVmf;

    PetscScalar dSf2_dPf_dthetat, dSf2_dQf_dthetat, dSt2_dPt_dthetat, dSt2_dQt_dthetat;
    dSf2_dPf_dthetat = dSf2_dPf_dPf*dPf_dthetat;
    dSf2_dQf_dthetat = dSf2_dQf_dQf*dQf_dthetat;
    dSt2_dPt_dthetat = dSt2_dPt_dPt*dPt_dthetat;
    dSt2_dQt_dthetat = dSt2_dQt_dQt*dQt_dthetat;

    PetscScalar dSf2_dPf_dVmt, dSf2_dQf_dVmt, dSt2_dPt_dVmt, dSt2_dQt_dVmt;
    dSf2_dPf_dVmt = dSf2_dPf_dPf*dPf_dVmt;
    dSf2_dQf_dVmt = dSf2_dQf_dQf*dQf_dVmt;
    dSt2_dPt_dVmt = dSt2_dPt_dPt*dPt_dVmt;
    dSt2_dQt_dVmt = dSt2_dQt_dQt*dQt_dVmt;


    PetscScalar dSf2_dthetaf_dthetaf,dSf2_dVmf_dthetaf,dSf2_dthetat_dthetaf,dSf2_dVmt_dthetaf;
    dSf2_dthetaf_dthetaf =  	dSf2_dPf_dthetaf*dPf_dthetaf 	+ dSf2_dPf*dPf_dthetaf_dthetaf
						  	+ 	dSf2_dQf_dthetaf*dQf_dthetaf 	+ dSf2_dQf*dQf_dthetaf_dthetaf;
    dSf2_dVmf_dthetaf    =  	dSf2_dPf_dthetaf*dPf_dVmf 		+ dSf2_dPf*dPf_dVmf_dthetaf    
						  	+ 	dSf2_dQf_dthetaf*dQf_dVmf 		+ dSf2_dQf*dQf_dVmf_dthetaf;	
    dSf2_dthetat_dthetaf =  	dSf2_dPf_dthetaf*dPf_dthetat 	+ dSf2_dPf*dPf_dthetat_dthetaf 
						  	+ 	dSf2_dQf_dthetaf*dQf_dthetat 	+ dSf2_dQf*dQf_dthetat_dthetaf;
    dSf2_dVmt_dthetaf    = 		dSf2_dPf_dthetaf*dPf_dVmt 		+ dSf2_dPf*dPf_dVmt_dthetaf    
							+ 	dSf2_dQf_dthetaf*dQf_dVmt    	+ dSf2_dQf*dQf_dVmt_dthetaf;

    PetscScalar dSf2_dVmf_dVmf,dSf2_dthetat_dVmf,dSf2_dVmt_dVmf;
    dSf2_dVmf_dVmf    	= 		dSf2_dPf_dVmf*dPf_dVmf    + dSf2_dPf*dPf_dVmf_dVmf    
							+   dSf2_dQf_dVmf*dQf_dVmf	  + dSf2_dQf*dQf_dVmf_dVmf;
    dSf2_dthetat_dVmf 	= 		dSf2_dPf_dVmf*dPf_dthetat + dSf2_dPf*dPf_dthetat_dVmf
							+	dSf2_dQf_dVmf*dQf_dthetat + dSf2_dQf*dQf_dthetat_dVmf;	
    dSf2_dVmt_dVmf    	= 		dSf2_dPf_dVmf*dPf_dVmt    + dSf2_dPf*dPf_dVmt_dVmf
							+	dSf2_dQf_dVmf*dQf_dVmt    + dSf2_dQf*dQf_dVmt_dVmf;

    PetscScalar dSf2_dthetat_dthetat,dSf2_dVmt_dthetat;
    dSf2_dthetat_dthetat = 		dSf2_dPf_dthetat*dPf_dthetat + dSf2_dPf*dPf_dthetat_dthetat
							+	dSf2_dQf_dthetat*dQf_dthetat + dSf2_dQf*dQf_dthetat_dthetat;	
    dSf2_dVmt_dthetat    = 		dSf2_dPf_dthetat*dPf_dVmt    + dSf2_dPf*dPf_dVmt_dthetat
							+	dSf2_dQf_dthetat*dQf_dVmt    + dSf2_dQf*dQf_dVmt_dthetat;
	
    PetscScalar dSf2_dVmt_dVmt;
    dSf2_dVmt_dVmt   	 = 		dSf2_dPf_dVmt*dPf_dVmt    + dSf2_dPf*dPf_dVmt_dVmt
							+	dSf2_dQf_dVmt*dQf_dVmt    + dSf2_dQf*dQf_dVmt_dVmt;

    PetscScalar dSt2_dthetaf_dthetaf,dSt2_dVmf_dthetaf,dSt2_dthetat_dthetaf,dSt2_dVmt_dthetaf;
    dSt2_dthetaf_dthetaf =  	dSt2_dPt_dthetaf*dPt_dthetaf 	+ dSt2_dPt*dPt_dthetaf_dthetaf
						  	+ 	dSt2_dQt_dthetaf*dQt_dthetaf 	+ dSt2_dQt*dQt_dthetaf_dthetaf;
    dSt2_dVmf_dthetaf    =  	dSt2_dPt_dthetaf*dPt_dVmf 		+ dSt2_dPt*dPt_dVmf_dthetaf    
						  	+ 	dSt2_dQt_dthetaf*dQt_dVmf 		+ dSt2_dQt*dQt_dVmf_dthetaf;	
    dSt2_dthetat_dthetaf =  	dSt2_dPt_dthetaf*dPt_dthetat 	+ dSt2_dPt*dPt_dthetaf_dthetat 
						  	+ 	dSt2_dQt_dthetaf*dQt_dthetat 	+ dSt2_dQt*dQt_dthetaf_dthetat;
    dSt2_dVmt_dthetaf    = 		dSt2_dPt_dthetaf*dPt_dVmt 		+ dSt2_dPt*dPt_dthetaf_dVmt    
							+ 	dSt2_dQt_dthetaf*dQt_dVmt    	+ dSt2_dQt*dQt_dthetaf_dVmt;

    PetscScalar dSt2_dVmf_dVmf,dSt2_dthetat_dVmf,dSt2_dVmt_dVmf;
    dSt2_dVmf_dVmf    	= 		dSt2_dPt_dVmf*dPt_dVmf    + dSt2_dPt*dPt_dVmf_dVmf 
							+	dSt2_dQt_dVmf*dQt_dVmf    + dSt2_dQt*dQt_dVmf_dVmf;
    dSt2_dthetat_dVmf 	= 		dSt2_dPt_dVmf*dPt_dthetat + dSt2_dPt*dPt_dVmf_dthetat
							+	dSt2_dQt_dVmf*dQt_dthetat + dSt2_dQt*dQt_dVmf_dthetat;	
    dSt2_dVmt_dVmf    	= 		dSt2_dPt_dVmf*dPt_dVmt    + dSt2_dPt*dPt_dVmf_dVmt
							+	dSt2_dQt_dVmf*dQt_dVmt    + dSt2_dQt*dQt_dVmf_dVmt;

    PetscScalar dSt2_dthetat_dthetat,dSt2_dVmt_dthetat;
    dSt2_dthetat_dthetat = 		dSt2_dPt_dthetat*dPt_dthetat + dSt2_dPt*dPt_dthetat_dthetat
							+	dSt2_dQt_dthetat*dQt_dthetat + dSt2_dQt*dQt_dthetat_dthetat;	
    dSt2_dVmt_dthetat    = 		dSt2_dPt_dthetat*dPt_dVmt    + dSt2_dPt*dPt_dVmt_dthetat
							+	dSt2_dQt_dthetat*dQt_dVmt    + dSt2_dQt*dQt_dVmt_dthetat;
	
    PetscScalar dSt2_dVmt_dVmt;
    dSt2_dVmt_dVmt   	 = 		dSt2_dPt_dVmt*dPt_dVmt    + dSt2_dPt*dPt_dVmt_dVmt
							+	dSt2_dQt_dVmt*dQt_dVmt    + dSt2_dQt*dQt_dVmt_dVmt;


    values[ctr]   = lambda[gloc]*dSf2_dthetaf_dthetaf 	+ lambda[gloc+1]*dSt2_dthetaf_dthetaf;
    values[ctr+1] = lambda[gloc]*dSf2_dVmf_dthetaf 		+ lambda[gloc+1]*dSt2_dVmf_dthetaf;
    values[ctr+2] = lambda[gloc]*dSf2_dVmf_dVmf 		+ lambda[gloc+1]*dSt2_dVmf_dVmf;
    values[ctr+3] = lambda[gloc]*dSf2_dthetat_dthetaf 	+ lambda[gloc+1]*dSt2_dthetat_dthetaf;
    values[ctr+4] = lambda[gloc]*dSf2_dthetat_dVmf		+ lambda[gloc+1]*dSt2_dthetat_dVmf;
    values[ctr+5] = lambda[gloc]*dSf2_dthetat_dthetat 	+ lambda[gloc+1]*dSt2_dthetat_dthetat;
    values[ctr+6] = lambda[gloc]*dSf2_dVmt_dthetaf 		+ lambda[gloc+1]*dSt2_dVmt_dthetaf;
    values[ctr+7] = lambda[gloc]*dSf2_dVmt_dVmf 		+ lambda[gloc+1]*dSt2_dVmt_dVmf;
    values[ctr+8] = lambda[gloc]*dSf2_dVmt_dthetat 		+ lambda[gloc+1]*dSt2_dVmt_dthetat;
    values[ctr+9] = lambda[gloc]*dSf2_dVmt_dVmt 		+ lambda[gloc+1]*dSt2_dVmt_dVmt;

    ctr  += 10;
	gloc +=  2;
  }

  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    PetscScalar theta,Vm;

    ierr = PSBUSGetVariableLocation(bus,&xloc);CHKERRQ(ierr);

    theta = x[xloc];
    Vm    = x[xloc+1];

    values[ctr] = lambda[gloc]*2*bus->gl + lambda[gloc+1]*(-2*bus->bl);
    ctr += 1;

    ierr = PSBUSGetSupportingLines(bus,&nconnlines,&connlines);CHKERRQ(ierr);

	for(k=0; k < nconnlines; k++) {
      line = connlines[k];
     
      if(!line->status) continue;

      PetscScalar Gff,Bff,Gft,Bft,Gtf,Btf,Gtt,Btt;
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
     
      PetscScalar Vmf,Vmt,thetaf,thetat,thetaft,thetatf;
      
      thetaf  = x[xlocf];
      Vmf     = x[xlocf+1];
      thetat  = x[xloct];
      Vmt     = x[xloct+1];
      thetaft = thetaf - thetat;
      thetatf = thetat - thetaf;
    
     if(bus == busf) {
	 
		PetscScalar dPf_dthetaf_dthetaf,dPf_dVmf_dthetaf,dPf_dthetat_dthetaf,dPf_dVmt_dthetaf;
		PetscScalar dPf_dVmf_dVmf,dPf_dthetat_dVmf,dPf_dVmt_dVmf;
		PetscScalar dPf_dthetat_dthetat,dPf_dVmt_dthetat;
		PetscScalar dPf_dVmt_dVmt;

		dPf_dthetaf_dthetaf =			Vmf*Vmt*(-Gft*cos(thetaft) - Bft*sin(thetaft));
		dPf_dVmf_dthetaf	=				Vmt*(-Gft*sin(thetaft) + Bft*cos(thetaft));
		dPf_dthetat_dthetaf =			Vmf*Vmt*( Gft*cos(thetaft) + Bft*sin(thetaft));
		dPf_dVmt_dthetaf	=				Vmf*(-Gft*sin(thetaft) + Bft*cos(thetaft));
		
		dPf_dVmf_dVmf		= 2*Gff;
		dPf_dthetat_dVmf	=				Vmt*( Gft*sin(thetaft) - Bft*cos(thetaft));
		dPf_dVmt_dVmf		=					( Gft*cos(thetaft) + Bft*sin(thetaft));
		
		dPf_dthetat_dthetat =			Vmf*Vmt*(-Gft*cos(thetaft) - Bft*sin(thetaft));
		dPf_dVmt_dthetat	=				Vmf*( Gft*sin(thetaft) - Bft*cos(thetaft));
		
		dPf_dVmt_dVmt		=				0.;
		
		PetscScalar dQf_dthetaf_dthetaf,dQf_dVmf_dthetaf,dQf_dthetat_dthetaf,dQf_dVmt_dthetaf;
		PetscScalar dQf_dVmf_dVmf,dQf_dthetat_dVmf,dQf_dVmt_dVmf;
		PetscScalar dQf_dthetat_dthetat,dQf_dVmt_dthetat;
		PetscScalar dQf_dVmt_dVmt;
		
		dQf_dthetaf_dthetaf =			Vmf*Vmt*( Bft*cos(thetaft) - Gft*sin(thetaft));
		dQf_dVmf_dthetaf	=				Vmt*( Bft*sin(thetaft) + Gft*cos(thetaft));
		dQf_dthetat_dthetaf =			Vmf*Vmt*(-Bft*cos(thetaft) + Gft*sin(thetaft));
		dQf_dVmt_dthetaf	=				Vmf*( Bft*sin(thetaft) + Gft*cos(thetaft));
		
		dQf_dVmf_dVmf		= -2*Bff;
		dQf_dthetat_dVmf	=				Vmt*(-Bft*sin(thetaft) - Gft*cos(thetaft));
		dQf_dVmt_dVmf		=					(-Bft*cos(thetaft) + Gft*sin(thetaft));
		
		dQf_dthetat_dthetat =			Vmf*Vmt*( Bft*cos(thetaft) - Gft*sin(thetaft));
		dQf_dVmt_dthetat	=				Vmf*(-Bft*sin(thetaft) - Gft*cos(thetaft));
		
		dQf_dVmt_dVmt		=				0.;
		
		values[ctr]   = lambda[gloc]*dPf_dthetaf_dthetaf 	+ lambda[gloc+1]*dQf_dthetaf_dthetaf;
		values[ctr+1] = lambda[gloc]*dPf_dVmf_dthetaf 		+ lambda[gloc+1]*dQf_dVmf_dthetaf;
		values[ctr+2] = lambda[gloc]*dPf_dVmf_dVmf 			+ lambda[gloc+1]*dQf_dVmf_dVmf;
		values[ctr+3] = lambda[gloc]*dPf_dthetat_dthetaf 	+ lambda[gloc+1]*dQf_dthetat_dthetaf;
		values[ctr+4] = lambda[gloc]*dPf_dthetat_dVmf		+ lambda[gloc+1]*dQf_dthetat_dVmf;
		values[ctr+5] = lambda[gloc]*dPf_dthetat_dthetat 	+ lambda[gloc+1]*dQf_dthetat_dthetat;
		values[ctr+6] = lambda[gloc]*dPf_dVmt_dthetaf 		+ lambda[gloc+1]*dQf_dVmt_dthetaf;
		values[ctr+7] = lambda[gloc]*dPf_dVmt_dVmf 			+ lambda[gloc+1]*dQf_dVmt_dVmf;
		values[ctr+8] = lambda[gloc]*dPf_dVmt_dthetat 		+ lambda[gloc+1]*dQf_dVmt_dthetat;
		
		ctr  += 9;
     } else {
  
	  	PetscScalar dPt_dthetat_dthetat,dPt_dVmt_dthetat,dPt_dthetaf_dthetat,dPt_dVmf_dthetat;
	  	PetscScalar dPt_dVmt_dVmt,dPt_dthetaf_dVmt,dPt_dVmf_dVmt;
	  	PetscScalar dPt_dthetaf_dthetaf,dPt_dVmf_dthetaf;
	  	PetscScalar dPt_dVmf_dVmf;
	  
	  	dPt_dthetat_dthetat = 		  Vmt*Vmf*(-Gtf*cos(thetatf) - Btf*sin(thetatf));
	  	dPt_dVmt_dthetat	= 			  Vmf*(-Gtf*sin(thetatf) + Btf*cos(thetatf));
	  	dPt_dthetaf_dthetat = 		  Vmt*Vmf*( Gtf*cos(thetatf) + Btf*sin(thetatf));
	  	dPt_dVmf_dthetat	= 			  Vmt*(-Gtf*sin(thetatf) + Btf*cos(thetatf));
	  
	  	dPt_dVmt_dVmt 	  	= 2*Gtt;
	  	dPt_dthetaf_dVmt	= 			  Vmf*( Gtf*sin(thetatf) - Btf*cos(thetatf));
	  	dPt_dVmf_dVmt 	  	= 				  ( Gtf*cos(thetatf) + Btf*sin(thetatf));
	  
	  	dPt_dthetaf_dthetaf = 		  Vmt*Vmf*(-Gtf*cos(thetatf) - Btf*sin(thetatf));
	  	dPt_dVmf_dthetaf	= 			  Vmt*( Gtf*sin(thetatf) - Btf*cos(thetatf));
	  
	  	dPt_dVmf_dVmf 	  	= 			  0.;
	  
	  	PetscScalar dQt_dthetat_dthetat,dQt_dVmt_dthetat,dQt_dthetaf_dthetat,dQt_dVmf_dthetat;
	  	PetscScalar dQt_dVmt_dVmt,dQt_dthetaf_dVmt,dQt_dVmf_dVmt;
	  	PetscScalar dQt_dthetaf_dthetaf,dQt_dVmf_dthetaf;
	  	PetscScalar dQt_dVmf_dVmf;  
	  
	  	dQt_dthetat_dthetat = 		  Vmt*Vmf*( Btf*cos(thetatf) - Gtf*sin(thetatf));
	  	dQt_dVmt_dthetat	= 			  Vmf*( Btf*sin(thetatf) + Gtf*cos(thetatf));
	  	dQt_dthetaf_dthetat = 		  Vmt*Vmf*(-Btf*cos(thetatf) + Gtf*sin(thetatf));
	  	dQt_dVmf_dthetat	= 			  Vmt*( Btf*sin(thetatf) + Gtf*cos(thetatf));
	  
	  	dQt_dVmt_dVmt 	  	= -2*Btt;
	  	dQt_dthetaf_dVmt	= 			  Vmf*(-Btf*sin(thetatf) - Gtf*cos(thetatf));
	  	dQt_dVmf_dVmt 	  	= 				  (-Btf*cos(thetatf) + Gtf*sin(thetatf));
	  
	  	dQt_dthetaf_dthetaf = 		  Vmt*Vmf*( Btf*cos(thetatf) - Gtf*sin(thetatf));
	  	dQt_dVmf_dthetaf	= 			  Vmt*(-Btf*sin(thetatf) - Gtf*cos(thetatf));
	  	
	  	dQt_dVmf_dVmf 	  	= 			  0.;

		values[ctr]   = lambda[gloc]*dPt_dthetaf_dthetaf 	+ lambda[gloc+1]*dQt_dthetaf_dthetaf;
		values[ctr+1] = lambda[gloc]*dPt_dVmf_dthetaf 		+ lambda[gloc+1]*dQt_dVmf_dthetaf;
		values[ctr+2] = lambda[gloc]*dPt_dthetaf_dthetat 	+ lambda[gloc+1]*dQt_dthetaf_dthetat;
		values[ctr+3] = lambda[gloc]*dPt_dVmf_dthetat 		+ lambda[gloc+1]*dQt_dVmf_dthetat;
		values[ctr+4] = lambda[gloc]*dPt_dthetat_dthetat 	+ lambda[gloc+1]*dQt_dthetat_dthetat;
		values[ctr+5] = lambda[gloc]*dPt_dthetaf_dVmt 		+ lambda[gloc+1]*dQt_dthetaf_dVmt;
		values[ctr+6] = lambda[gloc]*dPt_dVmf_dVmt 			+ lambda[gloc+1]*dQt_dVmf_dVmt;
		values[ctr+7] = lambda[gloc]*dPt_dVmt_dthetat 		+ lambda[gloc+1]*dQt_dVmt_dthetat;
		values[ctr+8] = lambda[gloc]*dPt_dVmt_dVmt		 	+ lambda[gloc+1]*dQt_dVmt_dVmt;
	 
	 	ctr  += 9;
      }
    }
	gloc += 2;
  }
  
  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(Lambda,&lambda);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  PFLOWCreate - Creates an optimal power flow application object

  Input Parameters
. mpicomm - The MPI communicator

  Output Parameters
. opflowout - The optimal power flow application object
*/
PetscErrorCode OPFLOWCreate(MPI_Comm mpicomm, OPFLOW *opflowout)
{
  PetscErrorCode ierr;
  OPFLOW         opflow;

  PetscFunctionBegin;
  ierr = PetscCalloc1(1,&opflow);CHKERRQ(ierr);

  ierr = COMMCreate(mpicomm,&opflow->comm);CHKERRQ(ierr);

  if(opflow->comm->size > 1) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Optimal Power Flow not supported in parallel");

  ierr = PSCreate(mpicomm,&opflow->ps);CHKERRQ(ierr);

  /* Set the application with the PS object */
  ierr = PSSetApplication(opflow->ps,APP_ACOPF);CHKERRQ(ierr);

  opflow->Nconeq    = -1;
  opflow->Nconineq  = -1;
  opflow->Ncon      = -1;
  opflow->nlp       = NULL;

  opflow->setupcalled = PETSC_FALSE;
  
  *opflowout = opflow;
  PetscFunctionReturn(0);
}

/*
  OPFLOWDestroy - Destroys the optimal power flow application object

  Input Parameter
. opflow - The OPFLOW object to destroy
*/
PetscErrorCode OPFLOWDestroy(OPFLOW *opflow)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = COMMDestroy(&(*opflow)->comm);CHKERRQ(ierr);

  /* Solution vector */
  ierr = VecDestroy(&(*opflow)->X);CHKERRQ(ierr);

  /* Lower and upper bounds on X */
  ierr = VecDestroy(&(*opflow)->Xl);CHKERRQ(ierr);
  ierr = VecDestroy(&(*opflow)->Xu);CHKERRQ(ierr);
  
  /* Gradient of objective function */
  ierr = VecDestroy(&(*opflow)->gradobj);CHKERRQ(ierr);

  /* Constraints vector */
  ierr = VecDestroy(&(*opflow)->G);CHKERRQ(ierr);
  /* Lower and upper bounds on constraints */
  ierr = VecDestroy(&(*opflow)->Gl);CHKERRQ(ierr);
  ierr = VecDestroy(&(*opflow)->Gu);CHKERRQ(ierr);
  
  /* Jacobian of constraints */
  //  ierr = MatDestroy(&(*opflow)->Jac);CHKERRQ(ierr);

  /* Multipliers */
  ierr = VecDestroy(&(*opflow)->lambda_g);CHKERRQ(ierr);
  ierr = VecDestroy(&(*opflow)->lambda_xl);CHKERRQ(ierr);
  ierr = VecDestroy(&(*opflow)->lambda_xu);CHKERRQ(ierr);

  ierr = PSDestroy(&(*opflow)->ps);CHKERRQ(ierr);

  FreeIpoptProblem((*opflow)->nlp);CHKERRQ(ierr);
  ierr = PetscFree(*opflow);CHKERRQ(ierr);
  *opflow = 0;

  PetscFunctionReturn(0);
}

/*
  OPFLOWReadMatPowerData - Reads the network data given in MATPOWER data format 

  Input Parameter
+  opflow - The OPFLOW object
-  netfile - The name of the network file

*/
PetscErrorCode OPFLOWReadMatPowerData(OPFLOW opflow,const char netfile[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  /* Read MatPower data file and populate the PS data structure */
  ierr = PSReadMatPowerData(opflow->ps,netfile);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  OPFLOWCreateGlobalVector - Returns a global vector of the appropriate size
  and distribution conforming to the distribution of the PS object.

  Input Paramereters:
. opflow - the optimal power flow application object

  Output Parameters:
. vec - the global vector

  Notes:
  OPFLOWSetUp() must be called before calling this routine.
*/
PetscErrorCode OPFLOWCreateGlobalVector(OPFLOW opflow,Vec *vec)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if(!opflow->setupcalled) SETERRQ(opflow->comm->type,0,"PFLOWSetUp() must be called before calling PFLOWCreateGlobalVector");
  ierr = PSCreateGlobalVector(opflow->ps,vec);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  OPFLOWCreateConstraintJacobian - Returns a distributed matrix of appropriate size that can
                                   be used as the Jacobian


  Input Paramereters:
. opflow - the optimal power flow application object

  Output Parameters:
. mat - the jacobian of the constraints

  Notes:
  OPFLOWSetUp() must be called before calling this routine.
*/
PetscErrorCode OPFLOWCreateConstraintJacobian(OPFLOW opflow,Mat *mat)
{
  PetscErrorCode ierr;
  Mat            jac;
  PetscInt       m=opflow->m; /* Number of constraints */
  PetscInt       n=opflow->n; /* Number of variables */
  PS             ps=opflow->ps;
  PetscInt       i,k;
  PSLINE         line;
  PSBUS          bus;
  PSGEN          gen;
  PetscInt       mloc=0;
  const PSBUS    *connbuses;
  PetscInt       nconnlines;
  const PSLINE   *connlines;
  PSBUS          busf,bust;
  PetscInt       loc,locf,loct;
  PetscInt       *nnz;
  PetscInt       row[2],col[4];
  PetscScalar    val[4];

  PetscFunctionBegin;

  ierr = MatCreate(opflow->comm->type,&jac);CHKERRQ(ierr);
  ierr = MatSetSizes(jac,m,n,PETSC_DECIDE,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = MatSetType(jac,MATAIJ);CHKERRQ(ierr);

  /* Set up preallocation */
  ierr = PetscCalloc1(m,&nnz);CHKERRQ(ierr);

  for(i=0; i < ps->Nbranch; i++) {
    line = &ps->line[i];

    nnz[mloc] += 4; 
    nnz[mloc+1] += 4;

    mloc += 2;
  }

  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    ierr = PSBUSGetSupportingLines(bus,&nconnlines,&connlines);CHKERRQ(ierr);

    nnz[mloc]   += 2*nconnlines + 2 + 2*bus->ngen;
    nnz[mloc+1] += 2*nconnlines + 2 + 2*bus->ngen;

    mloc += 2;
  }

  ierr = MatSeqAIJSetPreallocation(jac,0,nnz);CHKERRQ(ierr);
  ierr = PetscFree(nnz);CHKERRQ(ierr);
  ierr = MatZeroEntries(jac);CHKERRQ(ierr);

  val[0] = val[1] = val[2] = val[3] = 0.0;
  /* Set up nonzero structure. 0 is inserted in the non-zero location */
  mloc = 0;
  for(i=0; i < ps->Nbranch; i++) {
    line = &ps->line[i];

    ierr = PSLINEGetConnectedBuses(line,&connbuses);CHKERRQ(ierr);
    busf = connbuses[0];
    bust = connbuses[1];

    ierr = PSBUSGetVariableLocation(busf,&locf);CHKERRQ(ierr);
    ierr = PSBUSGetVariableLocation(bust,&loct);CHKERRQ(ierr);

    /* From bus Sf */
    row[0] = mloc;
    col[0] = locf; col[1] = locf+1; col[2] = loct; col[3] = loct+1;
    ierr = MatSetValues(jac,1,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);

    /* To bus St */
    row[0] = mloc+1;
    col[0] = loct; col[1] = loct+1; col[2] = locf; col[3] = locf+1;
    ierr = MatSetValues(jac,1,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);
    
    mloc += 2;
  }

  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
    row[0] = mloc;
    row[1] = mloc+1;
    col[0] = loc; col[1] = loc+1;
    ierr = MatSetValues(jac,2,row,2,col,val,ADD_VALUES);CHKERRQ(ierr);


    for(k=0; k < bus->ngen; k++) {
      loc += 2;
      col[0] = loc; col[1] = loc+1;
      ierr = MatSetValues(jac,2,row,2,col,val,ADD_VALUES);CHKERRQ(ierr);
    }

    ierr = PSBUSGetSupportingLines(bus,&nconnlines,&connlines);CHKERRQ(ierr);
    for(k=0; k < nconnlines; k++) {
      line = connlines[k];
      ierr = PSLINEGetConnectedBuses(line,&connbuses);CHKERRQ(ierr);
      busf = connbuses[0];
      bust = connbuses[1];
      
      ierr = PSBUSGetVariableLocation(busf,&locf);CHKERRQ(ierr);
      ierr = PSBUSGetVariableLocation(bust,&loct);CHKERRQ(ierr);
      
      if(bus == busf) {
	row[0] = mloc;
	col[0] = locf; col[1] = locf+1; col[2] = loct; col[3] = loct+1;
	ierr = MatSetValues(jac,1,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);
	row[0] = mloc+1;
	col[0] = locf; col[1] = locf+1; col[2] = loct; col[3] = loct+1;
	ierr = MatSetValues(jac,1,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);
      } else {
	row[0] = mloc;
	col[0] = loct; col[1] = loct+1; col[2] = locf; col[3] = locf+1;
	ierr = MatSetValues(jac,1,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);
	row[0] = mloc+1;
	col[0] = loct; col[1] = loct+1; col[2] = locf; col[3] = locf+1;
	ierr = MatSetValues(jac,1,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);
      }
    }
    mloc += 2;
  }

  ierr = MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  *mat = jac;
  PetscFunctionReturn(0);
}

/*
  OPFLOWGetConstraintJacobianNonzeros - Gets the number of nonzeros in the constraint jacobian matrix

  Input Paramereters:
. opflow - the optimal power flow application object

  Output Parameters:
. nnz - number of nonzeros in the constraint Jacobian

  Notes:
  OPFLOWSetUp() must be called before calling this routine.
*/
PetscErrorCode OPFLOWGetConstraintJacobianNonzeros(OPFLOW opflow,PetscInt *nnz)
{
  PetscErrorCode ierr;
  PS             ps=opflow->ps;
  PetscInt       i,k;
  PSBUS          bus;
  PSLINE         line;
  PetscInt       nconnlines;
  const PSLINE   *connlines;

  PetscFunctionBegin;

  *nnz = 0;

  for(i=0; i < ps->Nbranch; i++) *nnz += 8;

  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    if(bus->ide == ISOLATED_BUS) {
	  *nnz += 2;
    }	
    *nnz += 2 + 2*bus->ngen;
    ierr = PSBUSGetSupportingLines(bus,&nconnlines,&connlines);CHKERRQ(ierr);
    for(k=0; k < nconnlines;k++) {
      line = connlines[k];
      if(!line->status) continue;
      *nnz += 8;
    }
  }

  PetscFunctionReturn(0);
}

/*
  OPFLOWSetVariableandConstraintBounds - Sets the bounds on variables and constraints

  Input Parameters:
. opflow - the OPFLOW object

  Output Parameters:
+ Xl     - vector of lower bound on variables
. Xu     - vector of upper bound on variables
. Gl     - vector of lower bound on constraints
. Gu     - vector of lower bound on constraints
*/
PetscErrorCode OPFLOWSetVariableandConstraintBounds(OPFLOW opflow, Vec Xl, Vec Xu, Vec Gl, Vec Gu)
{
  PetscErrorCode ierr;
  PS             ps=opflow->ps;
  PetscScalar    *xl,*xu,*gl,*gu;
  PetscInt       i;
  PSLINE         line;
  PSBUS          bus;
  PetscInt       loc,gloc=0;

  PetscFunctionBegin;

  /* Get array pointers */
  ierr = VecGetArray(Xl,&xl);CHKERRQ(ierr);
  ierr = VecGetArray(Xu,&xu);CHKERRQ(ierr);
  ierr = VecGetArray(Gl,&gl);CHKERRQ(ierr);
  ierr = VecGetArray(Gu,&gu);CHKERRQ(ierr);
  
  for(i=0; i < ps->Nbranch; i++) {

    line = &ps->line[i];

    /* Line flow inequality constraints */
    if(!line->status) gl[gloc] = gu[gloc] = gl[gloc+1] = gu[gloc+1] = 0.0;
    else {
      gl[gloc] = gl[gloc+1] = 0.0;
      gu[gloc] = gu[gloc+1] = (line->rateA == 0) ? 1E10:(line->rateA/ps->MVAbase)*(line->rateA/ps->MVAbase);
    }    
    gloc += 2;
  }

  for(i=0; i < ps->nbus; i++) {
    PetscInt k, ngen;

    bus = &ps->bus[i];

    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);

    /* Bounds on voltage angles and bounds on real power mismatch equality constraints */
    xl[loc] = -1e10; xu[loc] = 1e10;
    gl[gloc] = 0.0;   gu[gloc] = 0.0;

    /* Bounds on voltage magnitudes and bounds on reactive power mismatch equality constraints */
    xl[loc+1] = bus->Vmin; xu[loc+1] = bus->Vmax;
    gl[gloc+1] = 0.0;       gu[gloc+1] = 0.0;

    if(bus->ide == REF_BUS || bus->ide == ISOLATED_BUS) xl[loc] = xu[loc] = bus->va*PETSC_PI/180.0;
    if(bus->ide == ISOLATED_BUS) xl[loc+1] = xu[loc+1] = bus->vm;
    
    for(k=0; k < bus->ngen; k++) {
      PSGEN gen;

      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      loc = loc+2;
      if(!gen->status) xl[loc] = xu[loc] = xl[loc+1] = xu[loc+1] = 0.0;
      else {
	xl[loc] = gen->pb/ps->MVAbase; /* PGmin */
	xu[loc] = gen->pt/ps->MVAbase; /* PGmax */
	xl[loc+1] = gen->qb/ps->MVAbase; /* QGmin */
	xu[loc+1] = gen->qt/ps->MVAbase; /* QGmax */
      }
    }
	gloc += 2;
  }

  ierr = VecRestoreArray(Xl,&xl);CHKERRQ(ierr);
  ierr = VecRestoreArray(Xu,&xu);CHKERRQ(ierr);
  ierr = VecRestoreArray(Gl,&gl);CHKERRQ(ierr);
  ierr = VecRestoreArray(Gu,&gu);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  OPFLOWSetInitialGuess - Sets the initial guess for the optimization

  Input Parameters:
. opflow - the OPFLOW object

  Output Parameters:
+ X     - initial guess

  Notes:
   Sets X[i] = (Xl[i] + Xu[i])/2
*/
PetscErrorCode OPFLOWSetInitialGuess(OPFLOW opflow, Vec X)
{
  PetscErrorCode ierr;
  PS             ps=opflow->ps;
  const PetscScalar    *xl,*xu;
  PetscScalar    *x;
  PetscInt       i;
  PSBUS          bus;
  PetscInt       loc;

  PetscFunctionBegin;

  /* Get array pointers */
  ierr = VecGetArray(X,&x);CHKERRQ(ierr);
  ierr = VecGetArrayRead(opflow->Xl,&xl);CHKERRQ(ierr);
  ierr = VecGetArrayRead(opflow->Xu,&xu);CHKERRQ(ierr);
  
  for(i=0; i < ps->nbus; i++) {
    PetscInt k, ngen;

    bus = &ps->bus[i];

    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);

    /* Bounds on voltage angles and bounds on voltage magnitudes */
    x[loc]   = (xl[loc] + xu[loc])/2.0;
    x[loc+1] = (xl[loc+1] + xu[loc+1])/2.0;

    for(k=0; k < bus->ngen; k++) {
      PSGEN gen;

      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      loc = loc+2;

      x[loc]   = (xl[loc] + xu[loc])/2.0;   /* Initial guess for Pg */
      x[loc+1] = (xl[loc+1] + xu[loc+1])/2.0; /* Initial guess for Qg */
    }
  }

  ierr = VecRestoreArray(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(opflow->Xl,&xl);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(opflow->Xu,&xu);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  OPFLOWSolve - Solves the AC optimal power flow

  Input Parameters:
. opflow - the optimal power flow application object
*/
PetscErrorCode OPFLOWSolve(OPFLOW opflow)
{
  PetscErrorCode ierr;
  PetscScalar    *xl,*xu,*gl,*gu;
  char           hessiantype[100];
  PetscBool      flg;

  PetscFunctionBegin;
  if(!opflow->setupcalled) {
    ierr = OPFLOWSetUp(opflow);CHKERRQ(ierr);
  }

  if(opflow->nlp) {
    FreeIpoptProblem(opflow->nlp);
  }

  /* Set bounds on variables and constraints */
  ierr = OPFLOWSetVariableandConstraintBounds(opflow,opflow->Xl,opflow->Xu,opflow->Gl,opflow->Gu);CHKERRQ(ierr);

  ierr = VecGetArray(opflow->Xl,&xl);CHKERRQ(ierr);
  ierr = VecGetArray(opflow->Xu,&xu);CHKERRQ(ierr);
  ierr = VecGetArray(opflow->Gl,&gl);CHKERRQ(ierr);
  ierr = VecGetArray(opflow->Gu,&gu);CHKERRQ(ierr);

  /* Create IPOPT problem */
  opflow->nlp = CreateIpoptProblem(opflow->n,xl,xu,opflow->m,gl,gu,opflow->nnz_jac_g,opflow->nnz_hes,0,&eval_opflow_f,
				   &eval_opflow_g,&eval_opflow_grad_f,
				   &eval_opflow_jac_g,&eval_opflow_h);

  ierr = VecRestoreArray(opflow->Xl,&xl);CHKERRQ(ierr);
  ierr = VecRestoreArray(opflow->Xu,&xu);CHKERRQ(ierr);
  ierr = VecRestoreArray(opflow->Gl,&gl);CHKERRQ(ierr);
  ierr = VecRestoreArray(opflow->Gu,&gu);CHKERRQ(ierr);

  /* Options for IPOPT. This need to go through PetscOptionsBegin later */

  AddIpoptNumOption(opflow->nlp, "tol", 1e-4);
  AddIpoptNumOption(opflow->nlp, "acceptable_tol", 1e-4);

  AddIpoptNumOption(opflow->nlp, "mu_init", 1e-1);
/*
  AddIpoptNumOption(opflow->nlp, "bound_frac",1e-4);
  AddIpoptNumOption(opflow->nlp, "bound_push",1e-4);
  AddIpoptNumOption(opflow->nlp, "dual_inf_tol", 1e-2);
  AddIpoptNumOption(opflow->nlp, "compl_inf_tol", 1e-2);
  AddIpoptNumOption(opflow->nlp, "constr_viol_tol", 1e-2);
*/

  AddIpoptStrOption(opflow->nlp, "mu_strategy", "monotone");
  AddIpoptStrOption(opflow->nlp, "print_user_options", "yes");
  AddIpoptStrOption(opflow->nlp, "output_file", "ipopt.out");

  ierr = PetscOptionsGetString(NULL,NULL,"-opflow_hessian_type",hessiantype,sizeof(hessiantype),&flg);CHKERRQ(ierr);
  if(flg) {
    ierr = PetscStrcmp(hessiantype,"L-BFGS",&flg);CHKERRQ(ierr);
    if(flg) {
      AddIpoptStrOption(opflow->nlp, "hessian_approximation", "limited-memory");
    } else {
      ierr = PetscStrcmp(hessiantype,"EXACT",&flg);CHKERRQ(ierr);
      if(!flg) {
	SETERRQ(PETSC_COMM_SELF,0,"Incorrect input to OPFLOW Hessian type option -opflow_hessian_type\n\tAvailable types are L-BFGS and EXACT");
      }
    }
  }
  AddIpoptStrOption(opflow->nlp, "derivative_test", "first-order");

  /* Set Initial Guess */
  ierr = OPFLOWSetInitialGuess(opflow,opflow->X);CHKERRQ(ierr);

  PetscScalar *x;
  ierr = VecGetArray(opflow->X,&x);CHKERRQ(ierr);
  /* Solve */
  opflow->solve_status = IpoptSolve(opflow->nlp,x,NULL,&opflow->obj,NULL,NULL,NULL,opflow);

  ierr = VecRestoreArray(opflow->X,&x);CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_SELF,"Objective function value = %lf\n",opflow->obj);CHKERRQ(ierr);
  ierr = VecView(opflow->X,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  OPFLOWSetUp - Sets up an optimal power flow application object

  Input Parameters:
. opflow - the OPFLOW object

  Notes:
  This routine sets up the OPFLOW object and the underlying PS object. It
  also distributes the PS object when used in parallel.
*/
PetscErrorCode OPFLOWSetUp(OPFLOW opflow)
{
  PetscErrorCode ierr;
  PS             ps=opflow->ps;
  PetscScalar    *xl,*xu,*gl,*gu;

  PetscFunctionBegin;

  /* Set up PS object */
  ierr = PSSetUp(ps);CHKERRQ(ierr);

  opflow->Nconeq = 2*ps->Nbus;
  opflow->Nconineq = 2*ps->Nbranch;
  opflow->Ncon = opflow->Nconeq + opflow->Nconineq;

  /* Create the solution vector and upper and lower bounds */
  ierr = PSCreateGlobalVector(opflow->ps,&opflow->X);CHKERRQ(ierr);
  ierr = VecDuplicate(opflow->X,&opflow->Xl);CHKERRQ(ierr);
  ierr = VecDuplicate(opflow->X,&opflow->Xu);CHKERRQ(ierr);

  /* Create the constraint function vector and its bounds */
  ierr = VecCreate(ps->comm->type,&opflow->G);CHKERRQ(ierr);
  ierr = VecSetSizes(opflow->G,opflow->Ncon,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(opflow->G);CHKERRQ(ierr);
  ierr = VecDuplicate(opflow->G,&opflow->Gl);CHKERRQ(ierr);
  ierr = VecDuplicate(opflow->G,&opflow->Gu);CHKERRQ(ierr);

  /* Create the gradient vector */
  ierr = VecDuplicate(opflow->X,&opflow->gradobj);CHKERRQ(ierr);

  ierr = VecGetSize(opflow->X,&opflow->n);CHKERRQ(ierr);
  ierr = VecGetSize(opflow->G,&opflow->m);CHKERRQ(ierr);

  /* Create vectors for multipliers */
  ierr = VecDuplicate(opflow->G,&opflow->lambda_g);CHKERRQ(ierr);
  ierr = VecDuplicate(opflow->X,&opflow->lambda_xl);CHKERRQ(ierr);
  ierr = VecDuplicate(opflow->X,&opflow->lambda_xu);CHKERRQ(ierr);
  /* Constraint jacobian */
  ierr = OPFLOWGetConstraintJacobianNonzeros(opflow,&opflow->nnz_jac_g);CHKERRQ(ierr);
  //  ierr = OPFLOWCreateConstraintJacobian(opflow,&opflow->Jac);CHKERRQ(ierr);

  ierr = OPFLOWGetLagrangianHessianNonzeros(opflow,&opflow->nnz_hes);CHKERRQ(ierr);  

  opflow->setupcalled = PETSC_TRUE;
  PetscFunctionReturn(0);
}

/*
  OPFLOWObjectiveFunction - The objective function for the optimal power flow

  Input Parameters:
+ opflow - the OPFLOW object
. X      - the current iterate

  Output Parameters:
. obj - the objective function value (scalar)
*/
PetscErrorCode OPFLOWObjectiveFunction(OPFLOW opflow,Vec X, PetscScalar* obj)
{
  PetscErrorCode ierr;
  const PetscScalar *x;
  PS             ps=opflow->ps;
  PetscInt       i;
  PSBUS          bus;
  PetscInt       loc;

  PetscFunctionBegin;
  *obj = 0.0;
  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);

  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
    
    PetscInt k;
    PSGEN    gen;
    PetscScalar Pg;
    for(k=0; k < bus->ngen; k++) {
      loc = loc+2;
      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      if(!gen->status) continue;
      Pg = x[loc]*ps->MVAbase;
      *obj += gen->cost_alpha*Pg*Pg + gen->cost_beta*Pg + gen->cost_gamma;
    }
  }
  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

Bool eval_opflow_f(PetscInt n, PetscScalar* x, Bool new_x,
            PetscScalar* obj_value, UserDataPtr user_data)
{
  PetscErrorCode ierr;
  OPFLOW         opflow=(OPFLOW)user_data;
  
  ierr = VecPlaceArray(opflow->X,x);CHKERRQ(ierr);
  ierr = OPFLOWObjectiveFunction(opflow,opflow->X,obj_value);CHKERRQ(ierr);
  ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
  return TRUE;
}

/*
  OPFLOWObjGradientFunction - The gradient of the objective function for the optimal power flow

  Input Parameters:
+ opflow - the OPFLOW object
. X      - the current iterate

  Output Parameters:
. grad - the objective function gradient
*/
PetscErrorCode OPFLOWObjGradientFunction(OPFLOW opflow,Vec X, Vec grad)
{
  PetscErrorCode ierr;
  const PetscScalar *x;
  PetscScalar    *df;
  PS             ps=opflow->ps;
  PetscInt       i;
  PSBUS          bus;
  PetscInt       loc;

  PetscFunctionBegin;
  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecSet(grad,0.0);CHKERRQ(ierr);
  ierr = VecGetArray(grad,&df);CHKERRQ(ierr);

  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
    
    PetscInt k;
    PSGEN    gen;
    PetscScalar Pg;
    for(k=0; k < bus->ngen; k++) {
      loc = loc+2;
      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      if(!gen->status) continue;
      Pg = x[loc]*ps->MVAbase;
      df[loc] = ps->MVAbase*(2*gen->cost_alpha*Pg + gen->cost_beta);
    }
  }
  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(grad,&df);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

Bool eval_opflow_grad_f(PetscInt n, PetscScalar* x, Bool new_x,
                 PetscScalar* grad_f, UserDataPtr user_data)
{
  PetscErrorCode ierr;
  OPFLOW         opflow=(OPFLOW)user_data;
  
  ierr = VecPlaceArray(opflow->X,x);CHKERRQ(ierr);
  ierr = VecPlaceArray(opflow->gradobj,grad_f);CHKERRQ(ierr);
  ierr = OPFLOWObjGradientFunction(opflow,opflow->X,opflow->gradobj);CHKERRQ(ierr);
  ierr = VecResetArray(opflow->gradobj);CHKERRQ(ierr);
  ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
  return TRUE;
}

/*
  OPFLOWConstraintFunction - Evalulates the constraints for the optimal power flow

  Input Parameters:
+ opflow - the OPFLOW object
. X      - the current iterate

  Output Parameters:
. G  - vector of constraints
*/
PetscErrorCode OPFLOWConstraintFunction(OPFLOW opflow,Vec X,Vec G)
{
  PetscErrorCode ierr;
  const PetscScalar *x;
  PetscScalar    *g;
  PetscInt       m=opflow->m; /* Number of constraints */
  PetscInt       n=opflow->n; /* Number of variables */
  PS             ps=opflow->ps;
  PetscInt       i,k;
  PSLINE         line;
  PSBUS          bus;
  PSGEN          gen;
  PetscInt       gloc=0;
  const PSBUS    *connbuses;
  PetscInt       nconnlines;
  const PSLINE   *connlines;
  PSBUS          busf,bust;
  PetscInt       xloc,xlocf,xloct;
  PSLOAD         load;
  PetscScalar    theta,Vm;

  PetscFunctionBegin;
  ierr = VecSet(G,0.0);CHKERRQ(ierr);

  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecGetArray(G,&g);CHKERRQ(ierr);

  for(i=0; i < ps->Nbranch; i++) {
    line = &ps->line[i];

    if(!line->status) continue;

    PetscScalar Gff,Bff,Gft,Bft,Gtf,Btf,Gtt,Btt;
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

    PetscScalar Vmf,Vmt,thetaf,thetat,thetaft,thetatf;
      
    thetaf  = x[xlocf];
    Vmf     = x[xlocf+1];
    thetat  = x[xloct];
    Vmt     = x[xloct+1];
    thetaft = thetaf - thetat;
    thetatf = thetat - thetaf;
    
    PetscScalar Pf,Qf,Pt,Qt,Sf2,St2,flow_max,flow_max2;

    Pf = Gff*Vmf*Vmf  + Vmf*Vmt*(Gft*cos(thetaft) + Bft*sin(thetaft));
    Qf = -Bff*Vmf*Vmf + Vmf*Vmt*(-Bft*cos(thetaft) + Gft*sin(thetaft));
	
    Pt = Gtt*Vmt*Vmt  + Vmt*Vmf*(Gtf*cos(thetatf) + Btf*sin(thetatf));
    Qt = -Btt*Vmt*Vmt + Vmt*Vmf*(-Btf*cos(thetatf) + Gtf*sin(thetatf));

    Sf2 = Pf*Pf + Qf*Qf;
    St2 = Pt*Pt + Qt*Qt;

    g[gloc]   = Sf2;
    g[gloc+1] = St2;

    gloc = gloc + 2;

  }

  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    
    ierr = PSBUSGetVariableLocation(bus,&xloc);CHKERRQ(ierr);
    theta = x[xloc];
    Vm    = x[xloc+1];

    if(bus->ide == ISOLATED_BUS) { 
      g[gloc] =  theta - bus->va*PETSC_PI/180.0;
      g[gloc+1] = Vm - bus->vm;
      continue;
    }

    /* Shunt injections */
    g[gloc]   += Vm*Vm*bus->gl;
    g[gloc+1] += -Vm*Vm*bus->bl;

    for(k=0; k < bus->ngen; k++) {
      xloc += 2;
      PetscScalar Pg,Qg;
      Pg = x[xloc];
      Qg = x[xloc+1];
      
      g[gloc]   -= Pg;
      g[gloc+1] -= Qg;
    }

    for(k=0; k < bus->nload; k++) {
      ierr = PSBUSGetLoad(bus,k,&load);CHKERRQ(ierr);

      g[gloc]   += load->pl;
      g[gloc+1] += load->ql;
    }

    ierr = PSBUSGetSupportingLines(bus,&nconnlines,&connlines);CHKERRQ(ierr);

    for(k=0; k < nconnlines; k++) {
      line = connlines[k];
      
      if(!line->status) continue;
      
      PetscScalar Gff,Bff,Gft,Bft,Gtf,Btf,Gtt,Btt;
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
      
      PetscScalar Vmf,Vmt,thetaf,thetat,thetaft,thetatf;
      
      thetaf  = x[xlocf];
      Vmf     = x[xlocf+1];
      thetat  = x[xloct];
      Vmt     = x[xloct+1];
      thetaft = thetaf - thetat;
      thetatf = thetat - thetaf;
      
      PetscScalar Pf,Qf,Pt,Qt,Sf2,St2,flow_max,flow_max2;
      
      if(bus == busf) {
	Pf = Gff*Vmf*Vmf  + Vmf*Vmt*(Gft*cos(thetaft) + Bft*sin(thetaft));
	Qf = -Bff*Vmf*Vmf + Vmf*Vmt*(-Bft*cos(thetaft) + Gft*sin(thetaft));

	g[gloc]   += Pf; 
	g[gloc+1] += Qf;
      } else {
	Pt = Gtt*Vmt*Vmt  + Vmt*Vmf*(Gtf*cos(thetatf) + Btf*sin(thetatf));
	Qt = -Btt*Vmt*Vmt + Vmt*Vmf*(-Btf*cos(thetatf) + Gtf*sin(thetatf));

	g[gloc]   += Pt;
	g[gloc+1] += Qt;
      }
    }
    gloc += 2;

  }

  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(G,&g);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

Bool eval_opflow_g(PetscInt n, PetscScalar* x, Bool new_x,
             PetscInt m, PetscScalar* g, UserDataPtr user_data)
{
  PetscErrorCode ierr;
  OPFLOW         opflow=(OPFLOW)user_data;

  ierr = VecPlaceArray(opflow->X,x);CHKERRQ(ierr);
  ierr = VecPlaceArray(opflow->G,g);CHKERRQ(ierr);
  ierr = OPFLOWConstraintFunction(opflow,opflow->X,opflow->G);CHKERRQ(ierr);
  ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
  ierr = VecResetArray(opflow->G);CHKERRQ(ierr);
  return TRUE;
}

/*
  OPFLOWSetConstraintJacobianLocations - Sets the row and column nonzero locations for the
              constraint Jacobian

  Input Parameters:
+ opflow - the OPFLOW object

  Output Parameters:
+ row - the row locations
- col - the col locations

  Notes:
   The row and column locations should be such that J(row[i],col[i]) = val
*/
PetscErrorCode OPFLOWSetConstraintJacobianLocations(OPFLOW opflow, PetscInt *row, PetscInt *col)
{
  PetscErrorCode ierr;
  PS             ps=opflow->ps;
  PetscInt       ctr=0;
  PetscInt       i,k;
  PSLINE         line;
  PSBUS          bus;
  PSGEN          gen;
  PetscInt       nconnlines;
  const PSLINE   *connlines;
  const PSBUS    *connbuses;
  PetscInt       gloc=0,xloc,xlocf,xloct;
  PSBUS          busf,bust;

  PetscFunctionBegin;

  for(i=0; i < ps->Nbranch; i++) {
    line = &ps->line[i];

    ierr = PSLINEGetConnectedBuses(line,&connbuses);CHKERRQ(ierr);
    busf = connbuses[0];
    bust = connbuses[1];

    ierr = PSBUSGetVariableLocation(busf,&xlocf);CHKERRQ(ierr);
    ierr = PSBUSGetVariableLocation(bust,&xloct);CHKERRQ(ierr);

    /* From bus Sf */
    row[ctr] = row[ctr+1] = row[ctr+2] = row[ctr+3] = gloc;
    col[ctr] = xlocf; col[ctr+1] = xlocf+1; col[ctr+2] = xloct; col[ctr+3] = xloct+1;
    ctr += 4;
    /* To bus St */
    row[ctr] = row[ctr+1] = row[ctr+2] = row[ctr+3] = gloc+1;
    col[ctr] = xloct; col[ctr+1] = xloct+1; col[ctr+2] = xlocf; col[ctr+3] = xlocf+1;
    ctr += 4;

    gloc += 2;
  }

  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    ierr = PSBUSGetVariableLocation(bus,&xloc);CHKERRQ(ierr);

    if(bus->ide == ISOLATED_BUS) {
      row[ctr] = gloc; col[ctr] = xloc;
      row[ctr+1] = gloc+1; col[ctr+1] = xloc+1;
      ctr += 2;
    }

    row[ctr] = gloc;     col[ctr]   = xloc+1;
    row[ctr+1] = gloc+1; col[ctr+1] = xloc+1;
    ctr += 2;

    for(k=0; k < bus->ngen; k++) {
      xloc += 2;

      row[ctr] = gloc;     col[ctr] = xloc;
      row[ctr+1] = gloc+1; col[ctr+1] = xloc+1;
      ctr += 2;
    }

   ierr = PSBUSGetSupportingLines(bus,&nconnlines,&connlines);CHKERRQ(ierr);

   for(k=0; k < nconnlines; k++) {
     line = connlines[k];
     
     if(!line->status) continue;

     ierr = PSLINEGetConnectedBuses(line,&connbuses);CHKERRQ(ierr);
     busf = connbuses[0];
     bust = connbuses[1];
      
     ierr = PSBUSGetVariableLocation(busf,&xlocf);CHKERRQ(ierr);
     ierr = PSBUSGetVariableLocation(bust,&xloct);CHKERRQ(ierr);

     if(bus == busf) {
       /* Pf */
       row[ctr] = row[ctr+1] = row[ctr+2] = row[ctr+3] = gloc;
       col[ctr] = xlocf; col[ctr+1] = xlocf+1; col[ctr+2] = xloct; col[ctr+3] = xloct+1;
       ctr += 4;

       /* Qf */
       row[ctr] = row[ctr+1] = row[ctr+2] = row[ctr+3] = gloc+1;
       col[ctr] = xlocf; col[ctr+1] = xlocf+1; col[ctr+2] = xloct; col[ctr+3] = xloct+1;
       ctr += 4;
     } else {
       /* Pt */
       row[ctr] = row[ctr+1] = row[ctr+2] = row[ctr+3] = gloc;
       col[ctr] = xloct; col[ctr+1] = xloct+1; col[ctr+2] = xlocf; col[ctr+3] = xlocf+1;
       ctr += 4;

       /* Qt */
       row[ctr] = row[ctr+1] = row[ctr+2] = row[ctr+3] = gloc+1;
       col[ctr] = xloct; col[ctr+1] = xloct+1; col[ctr+2] = xlocf; col[ctr+3] = xlocf+1;
       ctr += 4;
     }
   }
   gloc += 2;
  }
  PetscFunctionReturn(0);
}

/*
  OPFLOWSetConstraintJacobianValues - Sets the nonzero values for the
              constraint Jacobian

  Input Parameters:
+ opflow - the OPFLOW object
- X      - the current iterate


  Output Parameters:
. values - the nonzero values in the constraint jacobian

  Notes:
   The values should be in the same row and col order as set in OPFLOWSetConstraintJacobianLocations
*/
PetscErrorCode OPFLOWSetConstraintJacobianValues(OPFLOW opflow, Vec X,PetscScalar *values)
{
  PetscErrorCode ierr;
  const PetscScalar *x;
  PS             ps=opflow->ps;
  PetscInt       ctr=0;
  PetscInt       i,k;
  PSLINE         line;
  PSBUS          bus;
  PSGEN          gen;
  PetscInt       nconnlines;
  const PSLINE   *connlines;
  const PSBUS    *connbuses;
  PetscInt       gloc=0,xloc,xlocf,xloct;
  PSBUS          busf,bust;

  PetscFunctionBegin;

  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);
  for(i=0; i < ps->Nbranch; i++) {
    line = &ps->line[i];

    PetscScalar Gff,Bff,Gft,Bft,Gtf,Btf,Gtt,Btt;
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

    PetscScalar Vmf,Vmt,thetaf,thetat,thetaft,thetatf;
      
    thetaf  = x[xlocf];
    Vmf     = x[xlocf+1];
    thetat  = x[xloct];
    Vmt     = x[xloct+1];
    thetaft = thetaf - thetat;
    thetatf = thetat - thetaf;
    
    PetscScalar Pf,Qf,Pt,Qt,Sf2,St2;

    Pf = Gff*Vmf*Vmf  + Vmf*Vmt*(Gft*cos(thetaft) + Bft*sin(thetaft));
    Qf = -Bff*Vmf*Vmf + Vmf*Vmt*(-Bft*cos(thetaft) + Gft*sin(thetaft));
	
    Pt = Gtt*Vmt*Vmt  + Vmt*Vmf*(Gtf*cos(thetatf) + Btf*sin(thetatf));
    Qt = -Btt*Vmt*Vmt + Vmt*Vmf*(-Btf*cos(thetatf) + Gtf*sin(thetatf));

    Sf2 = Pf*Pf + Qf*Qf;
    St2 = Pt*Pt + Qt*Qt;

    PetscScalar dSf2_dPf, dSf2_dQf, dSt2_dPt, dSt2_dQt;

    dSf2_dPf = 2*Pf;
    dSf2_dQf = 2*Qf;
    dSt2_dPt = 2*Pt;
    dSt2_dQt = 2*Qt;

    PetscScalar dPf_dthetaf,dPf_dVmf,dPf_dthetat,dPf_dVmt;
    PetscScalar dQf_dthetaf,dQf_dVmf,dQf_dthetat,dQf_dVmt;
    PetscScalar dPt_dthetaf,dPt_dVmf,dPt_dthetat,dPt_dVmt;
    PetscScalar dQt_dthetaf,dQt_dVmf,dQt_dthetat,dQt_dVmt;

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

    PetscScalar dSf2_dthetaf,dSf2_dVmf,dSf2_dthetat,dSf2_dVmt;

    dSf2_dthetaf = dSf2_dPf*dPf_dthetaf + dSf2_dQf*dQf_dthetaf;
    dSf2_dthetat = dSf2_dPf*dPf_dthetat + dSf2_dQf*dQf_dthetat;
    dSf2_dVmf    = dSf2_dPf*dPf_dVmf    + dSf2_dQf*dQf_dVmf;
    dSf2_dVmt    = dSf2_dPf*dPf_dVmt    + dSf2_dQf*dQf_dVmt;

    values[ctr]   = dSf2_dthetaf;
    values[ctr+1] = dSf2_dVmf;
    values[ctr+2] = dSf2_dthetat;
    values[ctr+3] = dSf2_dVmt;

    ctr += 4;

    PetscScalar dSt2_dthetaf,dSt2_dVmf,dSt2_dthetat,dSt2_dVmt;

    dSt2_dthetaf = dSt2_dPt*dPt_dthetaf + dSt2_dQt*dQt_dthetaf;
    dSt2_dthetat = dSt2_dPt*dPt_dthetat + dSt2_dQt*dQt_dthetat;
    dSt2_dVmf    = dSt2_dPt*dPt_dVmf    + dSt2_dQt*dQt_dVmf;
    dSt2_dVmt    = dSt2_dPt*dPt_dVmt    + dSt2_dQt*dQt_dVmt;

    values[ctr]   = dSt2_dthetat;
    values[ctr+1] = dSt2_dVmt;
    values[ctr+2] = dSt2_dthetaf;
    values[ctr+3] = dSt2_dVmf;

    ctr += 4;

  }

  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    PetscScalar theta,Vm;

    ierr = PSBUSGetVariableLocation(bus,&xloc);CHKERRQ(ierr);

    theta = x[xloc];
    Vm    = x[xloc+1];

    if(bus->ide == ISOLATED_BUS) {
      values[ctr] = 1.0;
      values[ctr+1] = 1.0;
      ctr += 2;
    }

    values[ctr] = 2*Vm*bus->gl;
    values[ctr+1] = -2*Vm*bus->bl;
    ctr += 2;

    for(k=0; k < bus->ngen; k++) {
      xloc += 2;

      values[ctr]   = -1;
      values[ctr+1] = -1;
      ctr += 2;
    }

   ierr = PSBUSGetSupportingLines(bus,&nconnlines,&connlines);CHKERRQ(ierr);

   for(k=0; k < nconnlines; k++) {
     line = connlines[k];
     
     if(!line->status) continue;

     PetscScalar Gff,Bff,Gft,Bft,Gtf,Btf,Gtt,Btt;
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
     
     PetscScalar Vmf,Vmt,thetaf,thetat,thetaft,thetatf;
      
     thetaf  = x[xlocf];
     Vmf     = x[xlocf+1];
     thetat  = x[xloct];
     Vmt     = x[xloct+1];
     thetaft = thetaf - thetat;
     thetatf = thetat - thetaf;
    
     PetscScalar Pf,Qf,Pt,Qt,Sf2,St2;

     if(bus == busf) {
       //       Pf = Gff*Vmf*Vmf  + Vmf*Vmt*(Gft*cos(thetaft) + Bft*sin(thetaft));
       //       Qf = -Bff*Vmf*Vmf + Vmf*Vmt*(-Bft*cos(thetaft) + Gft*sin(thetaft));

       PetscScalar dPf_dthetaf,dPf_dVmf,dPf_dthetat,dPf_dVmt;
       PetscScalar dQf_dthetaf,dQf_dVmf,dQf_dthetat,dQf_dVmt;

       dPf_dthetaf = Vmf*Vmt*(-Gft*sin(thetaft) + Bft*cos(thetaft));
       dPf_dVmf    = 2*Gff*Vmf + Vmt*(Gft*cos(thetaft) + Bft*sin(thetaft));
       dPf_dthetat = Vmf*Vmt*(Gft*sin(thetaft) - Bft*cos(thetaft));
       dPf_dVmt    = Vmf*(Gft*cos(thetaft) + Bft*sin(thetaft));
       
       values[ctr] = dPf_dthetaf;
       values[ctr+1] = dPf_dVmf;
       values[ctr+2] = dPf_dthetat;
       values[ctr+3] = dPf_dVmt;
       ctr += 4;

       dQf_dthetaf = Vmf*Vmt*(Bft*sin(thetaft) + Gft*cos(thetaft));
       dQf_dVmf    = -2*Bff*Vmf + Vmt*(-Bft*cos(thetaft) + Gft*sin(thetaft));
       dQf_dthetat = Vmf*Vmt*(-Bft*sin(thetaft) - Gft*cos(thetaft));
       dQf_dVmt    = Vmf*(-Bft*cos(thetaft) + Gft*sin(thetaft));

       values[ctr]   = dQf_dthetaf;
       values[ctr+1] = dQf_dVmf;
       values[ctr+2] = dQf_dthetat;
       values[ctr+3] = dQf_dVmt;
       
       ctr += 4;

     } else {
       //       Pt = Gtt*Vmt*Vmt  + Vmt*Vmf*(Gtf*cos(thetatf) + Btf*sin(thetatf));
       //       Qt = -Btt*Vmt*Vmt + Vmt*Vmf*(-Btf*cos(thetatf) + Gtf*sin(thetatf));

       PetscScalar dPt_dthetaf,dPt_dVmf,dPt_dthetat,dPt_dVmt;
       PetscScalar dQt_dthetaf,dQt_dVmf,dQt_dthetat,dQt_dVmt;


       dPt_dthetat = Vmt*Vmf*(-Gtf*sin(thetatf) + Btf*cos(thetatf));
       dPt_dVmt    = 2*Gtt*Vmt + Vmf*(Gtf*cos(thetatf) + Btf*sin(thetatf));
       dPt_dthetaf = Vmt*Vmf*(Gtf*sin(thetatf) - Btf*cos(thetatf));
       dPt_dVmf    = Vmt*(Gtf*cos(thetatf) + Btf*sin(thetatf));

       values[ctr]   = dPt_dthetat;
       values[ctr+1] = dPt_dVmt;
       values[ctr+2] = dPt_dthetaf;
       values[ctr+3] = dPt_dVmf;

       ctr +=4 ;

       dQt_dthetat = Vmt*Vmf*(Btf*sin(thetatf) + Gtf*cos(thetatf));
       dQt_dVmt    = -2*Btt*Vmt + Vmf*(-Btf*cos(thetatf) + Gtf*sin(thetatf));
       dQt_dthetaf = Vmt*Vmf*(-Btf*sin(thetatf) - Gtf*cos(thetatf));
       dQt_dVmf    = Vmt*(-Btf*cos(thetatf) + Gtf*sin(thetatf));

       values[ctr]   = dQt_dthetat;
       values[ctr+1] = dQt_dVmt;
       values[ctr+2] = dQt_dthetaf;
       values[ctr+3] = dQt_dVmf;

       ctr += 4;
     }
   }
  }
  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


Bool eval_opflow_jac_g(PetscInt n, PetscScalar *x, Bool new_x,
                PetscInt m, PetscInt nele_jac,
                PetscInt *iRow, PetscInt *jCol, PetscScalar *values,
                UserDataPtr user_data)
{
  PetscErrorCode ierr;
  OPFLOW         opflow=(OPFLOW)user_data;

  ierr = VecPlaceArray(opflow->X,x);CHKERRQ(ierr);
  if(values == NULL) {
    ierr = OPFLOWSetConstraintJacobianLocations(opflow,iRow,jCol);CHKERRQ(ierr);
  } else {
    ierr = OPFLOWSetConstraintJacobianValues(opflow,opflow->X,values);CHKERRQ(ierr);
  }
  ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
  return TRUE;
}

Bool eval_opflow_h(PetscInt n, PetscScalar *x, Bool new_x, PetscScalar obj_factor,
                PetscInt m, PetscScalar *lambda, Bool new_lambda, 
                PetscInt nele_hess, PetscInt *iRow, PetscInt *jCol, 
                PetscScalar *values, UserDataPtr user_data)
{
  PetscErrorCode ierr;
  OPFLOW         opflow=(OPFLOW)user_data;

  ierr = VecPlaceArray(opflow->X,x);CHKERRQ(ierr);
  ierr = VecPlaceArray(opflow->lambda_g,lambda);CHKERRQ(ierr);
  
  if(values == NULL) {
    ierr = OPFLOWSetLagrangianHessianLocations(opflow,iRow,jCol);CHKERRQ(ierr);
  } else {
    ierr = OPFLOWSetLagrangianHessianValues(opflow,obj_factor, opflow->X,opflow->lambda_g,values);CHKERRQ(ierr);
  }
  ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
  ierr = VecResetArray(opflow->lambda_g);CHKERRQ(ierr);
  
  return TRUE;
}

