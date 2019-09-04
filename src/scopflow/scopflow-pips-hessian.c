#include <private/psimpl.h>
#include <private/opflowimpl.h>
#include <private/scopflowpipsimpl.h>
#include <math.h>
#include <../src/mat/impls/sbaij/seq/sbaij.h>


extern PetscErrorCode SCOPFLOWComputeObjectiveHessian(SCOPFLOW,PetscInt,Vec,Mat);
/*
  SCOPFLOWComputeEqualityConstraintsHessian - Computes the Hessian for the equality constraints function part
  
  Input Parameters:
+ scopflow - the SCOPFLOW object
. scenario - the scenario number
. X        - solution vector X
- Lambda   - Lagrangian multiplier vector

  Output Parameters:
. H - the Hessian part for the equality constraints

*/
PetscErrorCode SCOPFLOWComputeEqualityConstraintsHessian(SCOPFLOW scopflow,PetscInt scenario,Vec X,Vec Lambda,Mat H) 
{
  PetscErrorCode ierr;
  OPFLOW         opflow=scopflow->opflows[scenario];
  PS             ps=opflow->ps;
  PetscInt       i,k;
  PSLINE         line;
  PSBUS          bus;
  PetscInt       nconnlines;
  const PSLINE   *connlines;
  const PSBUS    *connbuses;
  PetscInt       xloc,xlocf,xloct;
  PSBUS          busf,bust;
  const PetscScalar *x;
  const PetscScalar *lambda;
  PetscInt       gloc=0;
  PetscInt       row[12],col[12];
  PetscScalar    val[12];

  PetscFunctionBegin;

  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecGetArrayRead(Lambda,&lambda);CHKERRQ(ierr);

  // For equality constraints (power flow) */
  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    
    PetscScalar theta,Vm;
    
    ierr = PSBUSGetVariableLocation(bus,&xloc);CHKERRQ(ierr);
    
    theta = x[xloc];
    Vm    = x[xloc+1];
    
    row[0] = xloc + 1; col[0] = xloc + 1;
    val[0] = lambda[gloc]*2*bus->gl + lambda[gloc+1]*(-2*bus->bl);
    ierr = MatSetValues(H,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);
    
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
	
	dPf_dthetaf_dthetaf = Vmf*Vmt*(-Gft*PetscCosScalar(thetaft) - Bft*PetscSinScalar(thetaft));
	dPf_dVmf_dthetaf    = Vmt*(-Gft*PetscSinScalar(thetaft) + Bft*PetscCosScalar(thetaft));
	dPf_dthetat_dthetaf = Vmf*Vmt*( Gft*PetscCosScalar(thetaft) + Bft*PetscSinScalar(thetaft));
	dPf_dVmt_dthetaf    = Vmf*(-Gft*PetscSinScalar(thetaft) + Bft*PetscCosScalar(thetaft));
	
	dPf_dVmf_dVmf	    = 2*Gff;
	dPf_dthetat_dVmf    = Vmt*( Gft*PetscSinScalar(thetaft) - Bft*PetscCosScalar(thetaft));
	dPf_dVmt_dVmf	    =	  ( Gft*PetscCosScalar(thetaft) + Bft*PetscSinScalar(thetaft));
	
	dPf_dthetat_dthetat = Vmf*Vmt*(-Gft*PetscCosScalar(thetaft) - Bft*PetscSinScalar(thetaft));
	dPf_dVmt_dthetat    = Vmf*( Gft*PetscSinScalar(thetaft) - Bft*PetscCosScalar(thetaft));
	
	dPf_dVmt_dVmt	    = 0.;
	
	PetscScalar dQf_dthetaf_dthetaf,dQf_dVmf_dthetaf,dQf_dthetat_dthetaf,dQf_dVmt_dthetaf;
	PetscScalar dQf_dVmf_dVmf,dQf_dthetat_dVmf,dQf_dVmt_dVmf;
	PetscScalar dQf_dthetat_dthetat,dQf_dVmt_dthetat;
	PetscScalar dQf_dVmt_dVmt;
	
	dQf_dthetaf_dthetaf     = Vmf*Vmt*( Bft*PetscCosScalar(thetaft) - Gft*PetscSinScalar(thetaft));
	dQf_dVmf_dthetaf	=     Vmt*( Bft*PetscSinScalar(thetaft) + Gft*PetscCosScalar(thetaft));
	dQf_dthetat_dthetaf     = Vmf*Vmt*(-Bft*PetscCosScalar(thetaft) + Gft*PetscSinScalar(thetaft));
	dQf_dVmt_dthetaf	=     Vmf*( Bft*PetscSinScalar(thetaft) + Gft*PetscCosScalar(thetaft));
	
	dQf_dVmf_dVmf		= -2*Bff;
	dQf_dthetat_dVmf	=     Vmt*(-Bft*PetscSinScalar(thetaft) - Gft*PetscCosScalar(thetaft));
	dQf_dVmt_dVmf		=	  (-Bft*PetscCosScalar(thetaft) + Gft*PetscSinScalar(thetaft));
	
	dQf_dthetat_dthetat     = Vmf*Vmt*( Bft*PetscCosScalar(thetaft) - Gft*PetscSinScalar(thetaft));
	dQf_dVmt_dthetat	=				Vmf*(-Bft*PetscSinScalar(thetaft) - Gft*PetscCosScalar(thetaft));
	
	dQf_dVmt_dVmt		=				0.;

	row[0] 	= xlocf; 	
	row[1] 	= xlocf+1; 	col[0] 	= xlocf;
	val[0]   = lambda[gloc]*dPf_dthetaf_dthetaf 	+ lambda[gloc+1]*dQf_dthetaf_dthetaf;
	val[1] = lambda[gloc]*dPf_dVmf_dthetaf 		+ lambda[gloc+1]*dQf_dVmf_dthetaf;
	ierr = MatSetValues(H,2,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);

	row[0] 	= xlocf+1; 	col[0] 	= xlocf+1;
	val[0] = lambda[gloc]*dPf_dVmf_dVmf 			+ lambda[gloc+1]*dQf_dVmf_dVmf;
	ierr = MatSetValues(H,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);


	row[0] 	= xloct; 	
	col[0] 	= xlocf;
	col[1] 	= xlocf+1;
	col[2] 	= xloct;
	val[0] = lambda[gloc]*dPf_dthetat_dthetaf 	+ lambda[gloc+1]*dQf_dthetat_dthetaf;
	val[1] = lambda[gloc]*dPf_dthetat_dVmf		+ lambda[gloc+1]*dQf_dthetat_dVmf;
	val[2] = lambda[gloc]*dPf_dthetat_dthetat 	+ lambda[gloc+1]*dQf_dthetat_dthetat;
	ierr = MatSetValues(H,1,row,3,col,val,ADD_VALUES);CHKERRQ(ierr);

	row[0] 	= xloct+1;
 	col[0] 	= xlocf;
	col[1] 	= xlocf+1;
	col[2] 	= xloct;
	val[0] = lambda[gloc]*dPf_dVmt_dthetaf 		+ lambda[gloc+1]*dQf_dVmt_dthetaf;
	val[1] = lambda[gloc]*dPf_dVmt_dVmf 		+ lambda[gloc+1]*dQf_dVmt_dVmf;
	val[2] = lambda[gloc]*dPf_dVmt_dthetat 		+ lambda[gloc+1]*dQf_dVmt_dthetat;
	ierr = MatSetValues(H,1,row,3,col,val,ADD_VALUES);CHKERRQ(ierr);

      } else {
	
	PetscScalar dPt_dthetat_dthetat,dPt_dVmt_dthetat,dPt_dthetaf_dthetat,dPt_dVmf_dthetat;
	PetscScalar dPt_dVmt_dVmt,dPt_dthetaf_dVmt,dPt_dVmf_dVmt;
	PetscScalar dPt_dthetaf_dthetaf,dPt_dVmf_dthetaf;
	PetscScalar dPt_dVmf_dVmf;
	
	dPt_dthetat_dthetat = 		  Vmt*Vmf*(-Gtf*PetscCosScalar(thetatf) - Btf*PetscSinScalar(thetatf));
	dPt_dVmt_dthetat	= 			  Vmf*(-Gtf*PetscSinScalar(thetatf) + Btf*PetscCosScalar(thetatf));
	dPt_dthetaf_dthetat = 		  Vmt*Vmf*( Gtf*PetscCosScalar(thetatf) + Btf*PetscSinScalar(thetatf));
	dPt_dVmf_dthetat	= 			  Vmt*(-Gtf*PetscSinScalar(thetatf) + Btf*PetscCosScalar(thetatf));
	
	dPt_dVmt_dVmt 	  	= 2*Gtt;
	dPt_dthetaf_dVmt	= 			  Vmf*( Gtf*PetscSinScalar(thetatf) - Btf*PetscCosScalar(thetatf));
	dPt_dVmf_dVmt 	  	= 				  ( Gtf*PetscCosScalar(thetatf) + Btf*PetscSinScalar(thetatf));
	
	dPt_dthetaf_dthetaf = 		  Vmt*Vmf*(-Gtf*PetscCosScalar(thetatf) - Btf*PetscSinScalar(thetatf));
	dPt_dVmf_dthetaf	= 			  Vmt*( Gtf*PetscSinScalar(thetatf) - Btf*PetscCosScalar(thetatf));
	
	dPt_dVmf_dVmf 	  	= 			  0.;
	
	PetscScalar dQt_dthetat_dthetat,dQt_dVmt_dthetat,dQt_dthetaf_dthetat,dQt_dVmf_dthetat;
	PetscScalar dQt_dVmt_dVmt,dQt_dthetaf_dVmt,dQt_dVmf_dVmt;
	PetscScalar dQt_dthetaf_dthetaf,dQt_dVmf_dthetaf;
	PetscScalar dQt_dVmf_dVmf;  
	
	dQt_dthetat_dthetat = 		  Vmt*Vmf*( Btf*PetscCosScalar(thetatf) - Gtf*PetscSinScalar(thetatf));
	dQt_dVmt_dthetat	= 			  Vmf*( Btf*PetscSinScalar(thetatf) + Gtf*PetscCosScalar(thetatf));
	dQt_dthetaf_dthetat = 		  Vmt*Vmf*(-Btf*PetscCosScalar(thetatf) + Gtf*PetscSinScalar(thetatf));
	dQt_dVmf_dthetat	= 			  Vmt*( Btf*PetscSinScalar(thetatf) + Gtf*PetscCosScalar(thetatf));
	
	dQt_dVmt_dVmt 	  	= -2*Btt;
	dQt_dthetaf_dVmt	= 			  Vmf*(-Btf*PetscSinScalar(thetatf) - Gtf*PetscCosScalar(thetatf));
	dQt_dVmf_dVmt 	  	= 				  (-Btf*PetscCosScalar(thetatf) + Gtf*PetscSinScalar(thetatf));
	
	dQt_dthetaf_dthetaf = 		  Vmt*Vmf*( Btf*PetscCosScalar(thetatf) - Gtf*PetscSinScalar(thetatf));
	dQt_dVmf_dthetaf	= 			  Vmt*(-Btf*PetscSinScalar(thetatf) - Gtf*PetscCosScalar(thetatf));
	
	dQt_dVmf_dVmf 	  	= 			  0.;
	
	row[0]	= xlocf;
	row[1]  = xlocf+1;  col[0]  = xlocf;
	val[0]  = lambda[gloc]*dPt_dthetaf_dthetaf 	+ lambda[gloc+1]*dQt_dthetaf_dthetaf;
	val[1]  = lambda[gloc]*dPt_dVmf_dthetaf 	+ lambda[gloc+1]*dQt_dVmf_dthetaf;
	ierr = MatSetValues(H,2,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);

	row[0]  = xloct;
	col[0]  = xlocf;
	col[1]  = xlocf+1;
	col[2]  = xloct;
	val[0] = lambda[gloc]*dPt_dthetaf_dthetat 	+ lambda[gloc+1]*dQt_dthetaf_dthetat;
	val[1] = lambda[gloc]*dPt_dVmf_dthetat 		+ lambda[gloc+1]*dQt_dVmf_dthetat;
	val[2] = lambda[gloc]*dPt_dthetat_dthetat 	+ lambda[gloc+1]*dQt_dthetat_dthetat;
	ierr = MatSetValues(H,1,row,3,col,val,ADD_VALUES);CHKERRQ(ierr);

	row[0] = xloct + 1;
	col[0]  = xlocf;
	col[1]  = xlocf+1;
	col[2]  = xloct;
	col[3]  = xloct+1;
	val[0] = lambda[gloc]*dPt_dthetaf_dVmt 		+ lambda[gloc+1]*dQt_dthetaf_dVmt;
	val[1] = lambda[gloc]*dPt_dVmf_dVmt 			+ lambda[gloc+1]*dQt_dVmf_dVmt;
	val[2] = lambda[gloc]*dPt_dVmt_dthetat 		+ lambda[gloc+1]*dQt_dVmt_dthetat;
	val[3] = lambda[gloc]*dPt_dVmt_dVmt		 	+ lambda[gloc+1]*dQt_dVmt_dVmt;
	ierr = MatSetValues(H,1,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);

      }
    }
    gloc += 2;
  }

  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(Lambda,&lambda);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}


/*
  SCOPFLOWComputeInequalityConstraintsHessian - Computes the Inequaluty Constraints Hessian

  Input Parameters:
+ scopflow - the OPFLOW object
. scenario - the scenario number
. X        - the solution vector
- Lambda   - Lagrangian multipler vector

  Output Parameters:
+ H   - the Hessian matrix

*/
PetscErrorCode SCOPFLOWComputeInequalityConstraintsHessian(SCOPFLOW scopflow, PetscInt scenario,Vec X, Vec Lambda,Mat H)
{
  PetscErrorCode ierr;
  OPFLOW         opflow=scopflow->opflows[scenario];
  PS             ps=opflow->ps;
  PetscInt       i;
  PSLINE         line;
  const PSBUS    *connbuses;
  PetscInt       xlocf,xloct;
  PSBUS          busf,bust;
  const PetscScalar *x;
  const PetscScalar *lambda;
  PetscInt       gloc;
  PetscInt       row[12],col[12];
  PetscScalar    val[12];

  PetscFunctionBegin;

  gloc = opflow->nconeq; /* offset for the inequality constraints in the Lambda vector */
      
  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecGetArrayRead(Lambda,&lambda);CHKERRQ(ierr);

  // for the part of line constraints
  for(i=0; i < ps->Nbranch; i++) {
    line = &ps->line[i];


    if(!line->status) {
      gloc += 2;
      continue;
    }

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

    Pf =  Gff*Vmf*Vmf + Vmf*Vmt*( Gft*PetscCosScalar(thetaft) + Bft*PetscSinScalar(thetaft));
    Qf = -Bff*Vmf*Vmf + Vmf*Vmt*(-Bft*PetscCosScalar(thetaft) + Gft*PetscSinScalar(thetaft));
	
    Pt =  Gtt*Vmt*Vmt + Vmt*Vmf*( Gtf*PetscCosScalar(thetatf) + Btf*PetscSinScalar(thetatf));
    Qt = -Btt*Vmt*Vmt + Vmt*Vmf*(-Btf*PetscCosScalar(thetatf) + Gtf*PetscSinScalar(thetatf));

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

    PetscScalar dPf_dthetaf_dthetaf,dPf_dVmf_dthetaf,dPf_dthetat_dthetaf,dPf_dVmt_dthetaf;
    PetscScalar dPf_dVmf_dVmf,dPf_dthetat_dVmf,dPf_dVmt_dVmf;
    PetscScalar dPf_dthetat_dthetat,dPf_dVmt_dthetat;
    PetscScalar dPf_dVmt_dVmt;

    dPf_dthetaf_dthetaf = 			Vmf*Vmt*(-Gft*PetscCosScalar(thetaft) - Bft*PetscSinScalar(thetaft));
    dPf_dVmf_dthetaf    =				Vmt*(-Gft*PetscSinScalar(thetaft) + Bft*PetscCosScalar(thetaft));
    dPf_dthetat_dthetaf = 			Vmf*Vmt*( Gft*PetscCosScalar(thetaft) + Bft*PetscSinScalar(thetaft));
	dPf_dVmt_dthetaf    = 				Vmf*(-Gft*PetscSinScalar(thetaft) + Bft*PetscCosScalar(thetaft));

    dPf_dVmf_dVmf    	= 2.*Gff;
    dPf_dthetat_dVmf 	= 				Vmt*( Gft*PetscSinScalar(thetaft) - Bft*PetscCosScalar(thetaft));
    dPf_dVmt_dVmf    	= 					( Gft*PetscCosScalar(thetaft) + Bft*PetscSinScalar(thetaft));

    dPf_dthetat_dthetat = 			Vmf*Vmt*(-Gft*PetscCosScalar(thetaft) - Bft*PetscSinScalar(thetaft));
    dPf_dVmt_dthetat    = 				Vmf*( Gft*PetscSinScalar(thetaft) - Bft*PetscCosScalar(thetaft));

    dPf_dVmt_dVmt    	= 				0.;

    PetscScalar dQf_dthetaf_dthetaf,dQf_dVmf_dthetaf,dQf_dthetat_dthetaf,dQf_dVmt_dthetaf;
    PetscScalar dQf_dVmf_dVmf,dQf_dthetat_dVmf,dQf_dVmt_dVmf;
    PetscScalar dQf_dthetat_dthetat,dQf_dVmt_dthetat;
    PetscScalar dQf_dVmt_dVmt;

    dQf_dthetaf_dthetaf = 			Vmf*Vmt*( Bft*PetscCosScalar(thetaft) - Gft*PetscSinScalar(thetaft));
    dQf_dVmf_dthetaf    = 				Vmt*( Bft*PetscSinScalar(thetaft) + Gft*PetscCosScalar(thetaft));
    dQf_dthetat_dthetaf = 			Vmf*Vmt*(-Bft*PetscCosScalar(thetaft) + Gft*PetscSinScalar(thetaft));
    dQf_dVmt_dthetaf    = 				Vmf*( Bft*PetscSinScalar(thetaft) + Gft*PetscCosScalar(thetaft));

    dQf_dVmf_dVmf    	= -2.*Bff;
    dQf_dthetat_dVmf 	= 				Vmt*(-Bft*PetscSinScalar(thetaft) - Gft*PetscCosScalar(thetaft));
    dQf_dVmt_dVmf    	= 					(-Bft*PetscCosScalar(thetaft) + Gft*PetscSinScalar(thetaft));

    dQf_dthetat_dthetat = 			Vmf*Vmt*( Bft*PetscCosScalar(thetaft) - Gft*PetscSinScalar(thetaft));
    dQf_dVmt_dthetat    = 				Vmf*(-Bft*PetscSinScalar(thetaft) - Gft*PetscCosScalar(thetaft));

    dQf_dVmt_dVmt    	= 				0.;
	
    PetscScalar dPt_dthetat_dthetat,dPt_dVmt_dthetat,dPt_dthetaf_dthetat,dPt_dVmf_dthetat;
    PetscScalar dPt_dVmt_dVmt,dPt_dthetaf_dVmt,dPt_dVmf_dVmt;
    PetscScalar dPt_dthetaf_dthetaf,dPt_dVmf_dthetaf;
    PetscScalar dPt_dVmf_dVmf;
	
    dPt_dthetat_dthetat = 			Vmt*Vmf*(-Gtf*PetscCosScalar(thetatf) - Btf*PetscSinScalar(thetatf));
    dPt_dVmt_dthetat    = 				Vmf*(-Gtf*PetscSinScalar(thetatf) + Btf*PetscCosScalar(thetatf));
    dPt_dthetaf_dthetat = 			Vmt*Vmf*( Gtf*PetscCosScalar(thetatf) + Btf*PetscSinScalar(thetatf));
    dPt_dVmf_dthetat    = 				Vmt*(-Gtf*PetscSinScalar(thetatf) + Btf*PetscCosScalar(thetatf));

    dPt_dVmt_dVmt   	= 2.*Gtt;
    dPt_dthetaf_dVmt 	= 				Vmf*( Gtf*PetscSinScalar(thetatf) - Btf*PetscCosScalar(thetatf));
    dPt_dVmf_dVmt    	= 					( Gtf*PetscCosScalar(thetatf) + Btf*PetscSinScalar(thetatf));

    dPt_dthetaf_dthetaf = 			Vmt*Vmf*(-Gtf*PetscCosScalar(thetatf) - Btf*PetscSinScalar(thetatf));
    dPt_dVmf_dthetaf    = 				Vmt*( Gtf*PetscSinScalar(thetatf) - Btf*PetscCosScalar(thetatf));

    dPt_dVmf_dVmf    	= 				0.;

    PetscScalar dQt_dthetat_dthetat,dQt_dVmt_dthetat,dQt_dthetaf_dthetat,dQt_dVmf_dthetat;
    PetscScalar dQt_dVmt_dVmt,dQt_dthetaf_dVmt,dQt_dVmf_dVmt;
    PetscScalar dQt_dthetaf_dthetaf,dQt_dVmf_dthetaf;
    PetscScalar dQt_dVmf_dVmf;	

    dQt_dthetat_dthetat = 			Vmt*Vmf*( Btf*PetscCosScalar(thetatf) - Gtf*PetscSinScalar(thetatf));
    dQt_dVmt_dthetat    = 				Vmf*( Btf*PetscSinScalar(thetatf) + Gtf*PetscCosScalar(thetatf));
    dQt_dthetaf_dthetat = 			Vmt*Vmf*(-Btf*PetscCosScalar(thetatf) + Gtf*PetscSinScalar(thetatf));
    dQt_dVmf_dthetat    = 				Vmt*( Btf*PetscSinScalar(thetatf) + Gtf*PetscCosScalar(thetatf));

    dQt_dVmt_dVmt    	= -2.*Btt;
    dQt_dthetaf_dVmt 	= 				Vmf*(-Btf*PetscSinScalar(thetatf) - Gtf*PetscCosScalar(thetatf));
    dQt_dVmf_dVmt    	= 					(-Btf*PetscCosScalar(thetatf) + Gtf*PetscSinScalar(thetatf));

    dQt_dthetaf_dthetaf = 			Vmt*Vmf*( Btf*PetscCosScalar(thetatf) - Gtf*PetscSinScalar(thetatf));
    dQt_dVmf_dthetaf    = 				Vmt*(-Btf*PetscSinScalar(thetatf) - Gtf*PetscCosScalar(thetatf));

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

    row[0] 	= xlocf; 
    row[1] 	= xlocf+1; 
    col[0] 	= xlocf;
    val[0]   = lambda[gloc]*dSf2_dthetaf_dthetaf 	+ lambda[gloc+1]*dSt2_dthetaf_dthetaf;
    val[1] = lambda[gloc]*dSf2_dVmf_dthetaf 		+ lambda[gloc+1]*dSt2_dVmf_dthetaf;
    ierr = MatSetValues(H,2,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);

    row[0] 	= xlocf+1; 	col[0] 	= xlocf+1;
    val[0] = lambda[gloc]*dSf2_dVmf_dVmf 		+ lambda[gloc+1]*dSt2_dVmf_dVmf;
    ierr = MatSetValues(H,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);

    row[0] 	= xloct; 	
    col[0] 	= xlocf;
    col[1] 	= xlocf+1;
    col[2] 	= xloct;
    val[0] = lambda[gloc]*dSf2_dthetat_dthetaf 	+ lambda[gloc+1]*dSt2_dthetat_dthetaf;
    val[1] = lambda[gloc]*dSf2_dthetat_dVmf		+ lambda[gloc+1]*dSt2_dthetat_dVmf;
    val[2] = lambda[gloc]*dSf2_dthetat_dthetat 	+ lambda[gloc+1]*dSt2_dthetat_dthetat;
    ierr = MatSetValues(H,1,row,3,col,val,ADD_VALUES);CHKERRQ(ierr);


    row[0] 	= xloct+1; 	
    col[0] 	= xlocf;
    col[1] 	= xlocf+1;
    col[2] 	= xloct;
    col[3] 	= xloct+1;
    val[0] = lambda[gloc]*dSf2_dVmt_dthetaf 		+ lambda[gloc+1]*dSt2_dVmt_dthetaf;
    val[1] = lambda[gloc]*dSf2_dVmt_dVmf 		+ lambda[gloc+1]*dSt2_dVmt_dVmf;
    val[2] = lambda[gloc]*dSf2_dVmt_dthetat 		+ lambda[gloc+1]*dSt2_dVmt_dthetat;
    val[3] = lambda[gloc]*dSf2_dVmt_dVmt 		+ lambda[gloc+1]*dSt2_dVmt_dVmt;
    ierr = MatSetValues(H,1,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);
    
    gloc +=  2;
  }

  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(Lambda,&lambda);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  SCOPFLOWComputeHessian - Computes the Hessian

  Input Parameters:
+ scopflow - the OPFLOW object
. scenario - the scenario number
. X        - the solution vector
- Lambda   - Lagrangian multipler vector

  Output Parameters:
+ H   - the Hessian matrix

*/
PetscErrorCode SCOPFLOWComputeHessian(SCOPFLOW scopflow, PetscInt scenario,Vec X, Vec Lambda,Mat H)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = MatZeroEntries(H);CHKERRQ(ierr);

  /* Objective function Hessian */
  if(scenario == 0) {
    ierr = SCOPFLOWComputeObjectiveHessian(scopflow,scenario,X,H);CHKERRQ(ierr);
  } else {
    if(!scopflow->first_stage_gen_cost_only) {
      ierr = SCOPFLOWComputeObjectiveHessian(scopflow,scenario,X,H);CHKERRQ(ierr);
    }
  }
  
  /* Equality constraints Hessian */
  ierr = SCOPFLOWComputeEqualityConstraintsHessian(scopflow,scenario,X,Lambda,H);CHKERRQ(ierr);

  /* Inequality (branch flow) constraints Hessian */
  if(!scopflow->ignore_line_flow_constraints) {
    ierr = SCOPFLOWComputeInequalityConstraintsHessian(scopflow,scenario,X,Lambda,H);CHKERRQ(ierr);
  }

  ierr = MatAssemblyBegin(H,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(H,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/* 
   All functions related to SCOPFLOW hessian
   are in this file
*/

int str_eval_h(double* x0, double* x1, double* lambda, int* nz, double* elts,
	       int* rowidx, int* colptr, CallBackDataPtr cbd) {

  PetscErrorCode ierr;
  int row = cbd->row_node_id;
  int col = cbd->col_node_id;
  SCOPFLOW scopflow=(SCOPFLOW)cbd->prob;
  Mat_SeqAIJ  *sbaij;
  OPFLOW   opflow;
  PetscInt nrow,ncol;

  if(elts==NULL) {
    *nz = 0;
    opflow = scopflow->opflows[row];
    if(row == col) {
      PetscScalar *x;
      if(row == 0) x = x0;
      else x = x1;

      ierr = VecPlaceArray(opflow->X,x);CHKERRQ(ierr);
      ierr = VecPlaceArray(opflow->Lambda,lambda);CHKERRQ(ierr);
      /* Compute Hessian */
      ierr = SCOPFLOWComputeHessian(scopflow,row,opflow->X,opflow->Lambda,opflow->Hes);CHKERRQ(ierr);
      ierr = MatSetOption(opflow->Hes,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE);CHKERRQ(ierr);
      ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
      ierr = VecResetArray(opflow->Lambda);CHKERRQ(ierr);

      sbaij = (Mat_SeqAIJ*)opflow->Hes->data;
      *nz = sbaij->nz;
    }
  } else {
    opflow = scopflow->opflows[row];
    if(row == col) {
      PetscScalar *x;
      if(row == 0) x = x0;
      else x = x1;

      ierr = VecPlaceArray(opflow->X,x);CHKERRQ(ierr);
      ierr = VecPlaceArray(opflow->Lambda,lambda);CHKERRQ(ierr);
      /* Compute Hessian */
      ierr = SCOPFLOWComputeHessian(scopflow,row,opflow->X,opflow->Lambda,opflow->Hes);CHKERRQ(ierr);
      ierr = MatGetSize(opflow->Hes,&nrow,&ncol);CHKERRQ(ierr);

      ierr = VecResetArray(opflow->X);CHKERRQ(ierr);
      ierr = VecResetArray(opflow->Lambda);CHKERRQ(ierr);

      sbaij = (Mat_SeqAIJ*)opflow->Hes->data;

      ierr = PetscMemcpy(rowidx,sbaij->j,sbaij->nz*sizeof(PetscInt));CHKERRQ(ierr);
      ierr = PetscMemcpy(colptr,sbaij->i,(nrow+1)*sizeof(PetscInt));CHKERRQ(ierr);
      ierr = PetscMemcpy(elts,sbaij->a,sbaij->nz*sizeof(PetscScalar));CHKERRQ(ierr);

      /*      ierr = MatView(opflow->Hes,0);CHKERRQ(ierr);
      exit(1);
      */
    }
  }
  
  return 1;
}

