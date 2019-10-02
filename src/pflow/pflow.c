#include <private/psimpl.h>
#include <private/pflowimpl.h>

/*
  PFLOWCreate - Creates a power flow application object

  Input Parameters
. mpicomm - The MPI communicator

  Output Parameters
. pflowout - The power flow application object
*/
PetscErrorCode PFLOWCreate(MPI_Comm mpicomm, PFLOW *pflowout)
{
  PetscErrorCode ierr;
  PFLOW             pflow;

  PetscFunctionBegin;
  ierr = PetscCalloc1(1,&pflow);CHKERRQ(ierr);

  ierr = COMMCreate(mpicomm,&pflow->comm);CHKERRQ(ierr);

  ierr = PSCreate(mpicomm,&pflow->ps);CHKERRQ(ierr);

  /* Set the application with the PS object */
  ierr = PSSetApplication(pflow->ps,APP_ACPF);CHKERRQ(ierr);

  /* Create the nonlinear solver object */
  ierr = SNESCreate(pflow->comm->type,&pflow->snes);CHKERRQ(ierr);

  pflow->setupcalled = PETSC_FALSE;

  pflow->split_gen_power_within_limits = PETSC_FALSE;

  *pflowout = pflow;
  PetscFunctionReturn(0);
}

/*
  PFLOWDestroy - Destroys the power flow application object

  Input Parameter
. pflow - The PFLOW object to destroy
*/
PetscErrorCode PFLOWDestroy(PFLOW *pflow)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = COMMDestroy(&(*pflow)->comm);CHKERRQ(ierr);

  ierr = VecDestroy(&(*pflow)->X);CHKERRQ(ierr);
  ierr = MatDestroy(&(*pflow)->Jac);CHKERRQ(ierr);

  ierr = SNESDestroy(&(*pflow)->snes);CHKERRQ(ierr);

  ierr = PSDestroy(&(*pflow)->ps);CHKERRQ(ierr);

  ierr = PetscFree(*pflow);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  PFLOWReadMatPowerData - Reads the network data given in MATPOWER data format 

  Input Parameter
+  pflow - The PFLOW object
-  netfile - The name of the network file

*/
PetscErrorCode PFLOWReadMatPowerData(PFLOW pflow,const char netfile[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  /* Read MatPower data file and populate the PS data structure */
  ierr = PSReadMatPowerData(pflow->ps,netfile);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  PFLOWReadPSSERawData - Reads the network data given in PSSE raw data format 

  Input Parameter
+  pflow - The PFLOW object
-  netfile - The name of the network file 

*/
PetscErrorCode PFLOWReadPSSERawData(PFLOW pflow,const char netfile[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  /* Read MatPower data file and populate the PS data structure */
  ierr = PSReadPSSERawData(pflow->ps,netfile);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


/*
  PFLOWCreateGlobalVector - Returns a global vector of the appropriate size
  and distribution conforming to the distribution of the PS object.

  Input Paramereters:
. pflow - the power flow application object

  Output Parameters:
. vec - the global vector

  Notes:
  PFLOWSetUp() must be called before calling this routine.
*/
PetscErrorCode PFLOWCreateGlobalVector(PFLOW pflow,Vec *vec)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if(!pflow->setupcalled) SETERRQ(pflow->comm->type,0,"PFLOWSetUp() must be called before calling PFLOWCreateGlobalVector");
  ierr = PSCreateGlobalVector(pflow->ps,vec);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  PFLOWCreateMatrix - Returns a distributed matrix of appropriate size that can
   be used as the Jacobian


  Input Paramereters:
. pflow - the power flow application object

  Output Parameters:
. mat - the matrix

  Notes:
  PFLOWSetUp() must be called before calling this routine.
*/
PetscErrorCode PFLOWCreateMatrix(PFLOW pflow,Mat *mat)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if(!pflow->setupcalled) SETERRQ(pflow->comm->type,0,"PFLOWSetUp() must be called before calling PFLOWCreateMatrix");
  ierr = PSCreateMatrix(pflow->ps,mat);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  PFLOWJacobian - The Jacobian evaluation routine for the power flow
  
  Notes:
  See PETSc routine SNESSetJacobian

*/
PetscErrorCode PFLOWJacobian(SNES snes,Vec X,Mat J, Mat Jpre, void *ctx)
{
  PetscErrorCode ierr;
  PFLOW          pflow=(PFLOW)ctx;
  PS             ps=pflow->ps;
  Vec            localX;
  const PetscScalar *xarr;
  PetscBool      ghostbus;
  PetscInt       i,k;

  PetscFunctionBegin;
  ierr = MatZeroEntries(J);CHKERRQ(ierr);

  ierr = DMGetLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);

  ierr = DMGlobalToLocalBegin(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);

  ierr = VecGetArrayRead(localX,&xarr);CHKERRQ(ierr);

  for(i=0; i < ps->nbus; i++) {
    PetscScalar Vm;
    PetscInt    loc;
    PetscInt    locglob;
    PSBUS       bus;
    PetscInt    row[2],col[4];
    PetscScalar val[4];

    bus = &ps->bus[i];

    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
    ierr = PSBUSGetVariableGlobalLocation(bus,&locglob);CHKERRQ(ierr);
    ierr = PSBUSIsGhosted(bus,&ghostbus);CHKERRQ(ierr);

    Vm = xarr[loc+1];

    if(!ghostbus) {
      row[0] = locglob; row[1] = locglob+1;
      col[0] = locglob; col[1] = locglob+1;
      /* Isolated and reference bus */
      if(bus->ide == ISOLATED_BUS || bus->ide == REF_BUS) {
	val[0] = val[3] = 1.0;
	val[1] = val[2] = 0.0;
	ierr = MatSetValues(J,2,row,2,col,val,ADD_VALUES);CHKERRQ(ierr);
	continue;
      }
      
      /* Shunt injections and PV_BUS voltage magnitude constraint */
      val[0] = 0.0; val[1] = 2*Vm*bus->gl;
      val[2] = 0.0; 
      if(bus->ide != PV_BUS) val[3]= -2*Vm*bus->bl; /* Partial derivative for shunt contribution */
      else val[3] = 1.0; /* Partial derivative of voltage magnitude constraint */
      ierr = MatSetValues(J,2,row,2,col,val,ADD_VALUES);CHKERRQ(ierr);      
    }

    /* Partial derivatives of network equations */
    PetscInt nconnlines;
    const PSLINE *connlines;
    PSLINE line;

    /* Get the lines supporting the bus */
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

      const PSBUS *connbuses;
      PSBUS busf,bust;
      PetscInt locglobf,locglobt,locf,loct;
      PetscScalar thetaf,thetat,Vmf,Vmt,thetaft,thetatf;

      /* Get the connected buses to this line */
      ierr = PSLINEGetConnectedBuses(line,&connbuses);CHKERRQ(ierr);
      busf = connbuses[0];
      bust = connbuses[1];

      ierr = PSBUSGetVariableGlobalLocation(busf,&locglobf);CHKERRQ(ierr);
      ierr = PSBUSGetVariableGlobalLocation(bust,&locglobt);CHKERRQ(ierr);

      ierr = PSBUSGetVariableLocation(busf,&locf);CHKERRQ(ierr);
      ierr = PSBUSGetVariableLocation(bust,&loct);CHKERRQ(ierr);

      thetaf = xarr[locf];  Vmf = xarr[locf+1];
      thetat = xarr[loct];  Vmt = xarr[loct+1];
      thetaft = thetaf - thetat;
      thetatf = thetat - thetaf;

      if(bus == busf) {
	if(bus->ide != REF_BUS) {
	  row[0] = locglobf;
	  col[0] = locglobf; col[1] = locglobf+1; col[2] = locglobt; col[3] = locglobt+1;
	  val[0] = Vmf*Vmt*(-Gft*sin(thetaft) + Bft*cos(thetaft));
	  val[1] = 2*Gff*Vmf + Vmt*(Gft*cos(thetaft) + Bft*sin(thetaft));
	  val[2] = Vmf*Vmt*(Gft*sin(thetaft) - Bft*cos(thetaft));
	  val[3] = Vmf*(Gft*cos(thetaft) + Bft*sin(thetaft));
	  ierr = MatSetValues(J,1,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);

	  if(bus->ide != PV_BUS) {
	    row[0] = locglobf+1;
	    val[0] = Vmf*Vmt*(Bft*sin(thetaft) + Gft*cos(thetaft));
	    val[1] = -2*Bff*Vmf + Vmt*(-Bft*cos(thetaft) + Gft*sin(thetaft));
	    val[2] = Vmf*Vmt*(-Bft*sin(thetaft) - Gft*cos(thetaft));
	    val[3] = Vmf*(-Bft*cos(thetaft) + Gft*sin(thetaft));
	    ierr = MatSetValues(J,1,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);
	  }
	}
      } else {
	if(bus->ide != REF_BUS) {
	  row[0] = locglobt;
	  col[0] = locglobt; col[1] = locglobt+1; col[2] = locglobf; col[3] = locglobf+1;
	  val[0] = Vmt*Vmf*(-Gtf*sin(thetatf) + Btf*cos(thetatf));
	  val[1] = 2*Gtt*Vmt + Vmf*(Gtf*cos(thetatf) + Btf*sin(thetatf));
	  val[2] = Vmt*Vmf*(Gtf*sin(thetatf) - Btf*cos(thetatf));
	  val[3] = Vmt*(Gtf*cos(thetatf) + Btf*sin(thetatf));
	  ierr = MatSetValues(J,1,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);
	  
	  if(bus->ide != PV_BUS) {
	    row[0] = locglobt+1;
	    val[0] = Vmt*Vmf*(Btf*sin(thetatf) + Gtf*cos(thetatf));
	    val[1] = -2*Btt*Vmt + Vmf*(-Btf*cos(thetatf) + Gtf*sin(thetatf));
	    val[2] = Vmt*Vmf*(-Btf*sin(thetatf) - Gtf*cos(thetatf));
	    val[3] = Vmt*(-Btf*cos(thetatf) + Gtf*sin(thetatf));
	    ierr = MatSetValues(J,1,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);
	  }
	}
      }
    }
  }

  ierr = VecRestoreArrayRead(localX,&xarr);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);

  ierr = MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    
  //  ierr = MatView(J,0);
  //  exit(1);
  PetscFunctionReturn(0);
}

/*
  PFLOWFunction - The residual function for the power flow
  
  Notes:
  See PETSc routine SNESSetFunction()

  The equations are expressed in power balance form:
    S_network (=Y*V) + S_shunt - S_gen + S_load = 0

*/
PetscErrorCode PFLOWFunction(SNES snes,Vec X,Vec F,void *ctx)
{
  PetscErrorCode ierr;
  PFLOW          pflow=(PFLOW)ctx;
  PS             ps=pflow->ps;
  Vec            localX,localF;
  const PetscScalar    *xarr;
  PetscScalar     *farr;
  PetscBool      ghostbus;
  PetscInt       i;

  PetscFunctionBegin;
  ierr = VecSet(F,0.0);CHKERRQ(ierr);

  ierr = DMGetLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = VecGetArrayRead(localX,&xarr);CHKERRQ(ierr);

  ierr = DMGetLocalVector(ps->networkdm,&localF);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ps->networkdm,F,INSERT_VALUES,localF);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ps->networkdm,F,INSERT_VALUES,localF);CHKERRQ(ierr);
  ierr = VecGetArray(localF,&farr);CHKERRQ(ierr);

  for(i=0; i < ps->nbus; i++) {
    PetscScalar theta, Vm;
    PetscInt    loc;
    PSBUS       bus;
    PetscInt    j,k;

    bus = &ps->bus[i];

    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
    ierr = PSBUSIsGhosted(bus,&ghostbus);CHKERRQ(ierr);

    theta = xarr[loc]; Vm = xarr[loc+1];
    
    if(!ghostbus) {
      if(bus->ide == REF_BUS || bus->ide == ISOLATED_BUS) {
	farr[loc]   = theta - bus->va*PETSC_PI/180.0;
	farr[loc+1] = Vm - bus->vm;
	continue;
      }
  
      /* Shunt injections */
      farr[loc] += Vm*Vm*bus->gl;
      if(bus->ide != PV_BUS) farr[loc+1] += -Vm*Vm*bus->bl;

      /* Gen. injections */
      for(j=0; j < bus->ngen; j++) {
	PSGEN gen;

	ierr = PSBUSGetGen(bus,j,&gen);CHKERRQ(ierr);
	if(!gen->status) continue;
	farr[loc]  -= gen->pg;
	if(bus->ide == PV_BUS) farr[loc+1] = Vm - gen->vs;
	else farr[loc+1] -= gen->qg;
      }

      /* Load injections */
      for(j=0; j < bus->nload; j++) {
	PSLOAD load;

	ierr = PSBUSGetLoad(bus,j,&load);CHKERRQ(ierr);
	if(!load->status) continue;
	farr[loc]  += load->pl;
	if(bus->ide != PV_BUS) farr[loc+1] += load->ql;
      }
    }

    PetscInt nconnlines;
    const PSLINE *connlines;
    PSLINE line;
    /* Injection from the lines supporting the bus */
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
      
      const PSBUS *connbuses;
      PSBUS busf,bust;
      PetscInt locf,loct;
      PetscScalar thetaf,Vmf,thetat,Vmt,thetaft,thetatf;
      
      /* Get the connected buses to this line */
      ierr = PSLINEGetConnectedBuses(line,&connbuses);CHKERRQ(ierr);
      busf = connbuses[0];
      bust = connbuses[1];
      
      ierr = PSBUSGetVariableLocation(busf,&locf);CHKERRQ(ierr);
      ierr = PSBUSGetVariableLocation(bust,&loct);CHKERRQ(ierr);
      thetaf = xarr[locf]; Vmf = xarr[locf+1];
      thetat = xarr[loct]; Vmt = xarr[loct+1];
      thetaft = thetaf - thetat;
      thetatf = thetat - thetaf;
      
      if (bus == busf) {
	if(busf->ide != REF_BUS) farr[locf]   += Gff*Vmf*Vmf + Vmf*Vmt*(Gft*cos(thetaft) + Bft*sin(thetaft));
	if(busf->ide != PV_BUS && busf->ide != REF_BUS) farr[locf+1] += -Bff*Vmf*Vmf + Vmf*Vmt*(-Bft*cos(thetaft) + Gft*sin(thetaft));
      } else {
	if(bust->ide != REF_BUS) farr[loct]   += Gtt*Vmt*Vmt + Vmt*Vmf*(Gtf*cos(thetatf) + Btf*sin(thetatf));
	if(bust->ide != PV_BUS && bust->ide != REF_BUS) farr[loct+1] += -Btt*Vmt*Vmt + Vmt*Vmf*(-Btf*cos(thetatf) + Gtf*sin(thetatf));
      }
    } 
  }

  ierr = VecRestoreArrayRead(localX,&xarr);CHKERRQ(ierr);
  ierr = VecRestoreArray(localF,&farr);CHKERRQ(ierr);

  ierr = DMLocalToGlobalBegin(ps->networkdm,localF,ADD_VALUES,F);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(ps->networkdm,localF,ADD_VALUES,F);CHKERRQ(ierr);

  ierr = DMRestoreLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ps->networkdm,&localF);CHKERRQ(ierr);

  /* Only activate for debugging
  ierr = DMGetLocalVector(ps->networkdm,&localF);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ps->networkdm,F,INSERT_VALUES,localF);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ps->networkdm,F,INSERT_VALUES,localF);CHKERRQ(ierr);
  ierr = VecGetArray(localF,&farr);CHKERRQ(ierr);

  for(i=0; i < ps->nbus; i++) {
    PetscScalar ftheta, fVm;
    PetscInt    loc;
    PSBUS       bus;

    bus = &ps->bus[i];

    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
    ierr = PSBUSIsGhosted(bus,&ghostbus);CHKERRQ(ierr);

    ftheta = farr[loc]; fVm = farr[loc+1];
    
    if(!ghostbus) {
      ierr = PetscPrintf(PETSC_COMM_SELF,"Rank %d Bus %d ftheta %5.4f fVm %5.4f\n",pflow->comm->rank,bus->bus_i,ftheta,fVm);CHKERRQ(ierr);
    }
  }
  //  ierr = VecView(F,0);CHKERRQ(ierr);
  exit(1);
  ****/

  PetscFunctionReturn(0);
}

/*
  PFLOWSetInitialGuess - Initializes the network variables in the solution
                            vector X.

  Input Parameters:
. pflow - the power flow application object

  Output Parameters:
. X - the solution vector
*/
PetscErrorCode PFLOWSetInitialGuess(PFLOW pflow,Vec X)
{
  PetscErrorCode ierr;
  PS             ps=pflow->ps;
  Vec            localX;
  PetscScalar    *xarr;
  PetscInt       nbus=ps->nbus,i,loc;
  PetscBool      ghostbus;
  PSBUS          bus;
  PetscScalar    Vm=1.0,Va; /* Voltage magnitude and angle */

  PetscFunctionBegin;
  ierr = DMGetLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);

  ierr = VecSet(X,0.0);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);

  ierr = VecGetArray(localX,&xarr);CHKERRQ(ierr);

  for(i=0; i < nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSIsGhosted(bus,&ghostbus);CHKERRQ(ierr);
    if (ghostbus) continue;

    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);

    if (bus->ngen) {
      PetscInt k;
      PSGEN gen;
      PetscInt atleastonegenon=0;
      for(k=0; k < bus->ngen; k++) {
	ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
	atleastonegenon = atleastonegenon+gen->status;
	if(!gen->status) continue;
	Vm = gen->vs;
	bus->vm = Vm;
      }
      if(!atleastonegenon) Vm = bus->vm;
    } else Vm = bus->vm;
    Va = bus->va*PETSC_PI/180.0;

    xarr[loc] = Va;
    xarr[loc+1] = Vm;

    //    ierr = PetscPrintf(PETSC_COMM_SELF,"Rank %d Bus %d Type %d Vm %4.3f Va %4.3f\n",pflow->comm->rank,bus->bus_i,bus->ide,Vm,bus->va);CHKERRQ(ierr);


  }

  ierr = VecRestoreArray(localX,&xarr);CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(ps->networkdm,localX,ADD_VALUES,X);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(ps->networkdm,localX,ADD_VALUES,X);CHKERRQ(ierr);

  ierr = DMRestoreLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);

  //  ierr = VecView(X,0);

  PetscFunctionReturn(0);
}

/*
  PFLOWSolve - Solves the AC power flow equations

  Input Parameters:
. pflow - the power flow application object
*/
PetscErrorCode PFLOWSolve(PFLOW pflow)
{
  SNESConvergedReason reason;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if(!pflow->setupcalled) {
    ierr = PFLOWSetUp(pflow);CHKERRQ(ierr);
  }

  /* Check Topology */
  ierr = PSCheckTopology(pflow->ps);CHKERRQ(ierr);

  /* Set up initial guess */
  ierr = PFLOWSetInitialGuess(pflow,pflow->X);CHKERRQ(ierr);

  /* Call the nonlinear solver */
  ierr = SNESSolve(pflow->snes,NULL,pflow->X);CHKERRQ(ierr);
  ierr = SNESGetConvergedReason(pflow->snes,&reason);CHKERRQ(ierr);
  pflow->converged = reason< 0 ? PETSC_FALSE: PETSC_TRUE;
  PetscFunctionReturn(0);
}

/*
  PFLOWConverged - Returns the convergence status of the power flow

  Input Parameters:
+ pflow - the PFLOW object
- flg - flag indicating if the poer flow converges
*/
PetscErrorCode PFLOWConverged(PFLOW pflow,PetscBool *flg)
{
  PetscFunctionBegin;
  *flg = pflow->converged;
  PetscFunctionReturn(0);
}

/*
  PFLOWPostSolve - Updates the buses and the branches with the solution from the power flow

  Input Parameters:
. pflow - the PFLOW object
*/
PetscErrorCode PFLOWPostSolve(PFLOW pflow)
{
  PetscErrorCode     ierr;
  PS                 ps=pflow->ps;
  Vec                localX,localPQgen;
  const PetscScalar *xarr;
  PetscBool          ghostbus;
  PetscInt           i;
  Vec                PQgen;
  PetscScalar        Pnet,Qnet,*pqgenarr;
  PSLINE             line;

  PetscFunctionBegin;

  /* Create global vector PQgen of size 2*nbus X 1 for storing Pgen and Qgen at each bus */
  ierr = PFLOWCreateGlobalVector(pflow,&PQgen);CHKERRQ(ierr);
  ierr = VecSet(PQgen,0.0);CHKERRQ(ierr);
  ierr = DMGetLocalVector(ps->networkdm,&localPQgen);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ps->networkdm,PQgen,INSERT_VALUES,localPQgen);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ps->networkdm,PQgen,INSERT_VALUES,localPQgen);CHKERRQ(ierr);
  ierr = VecGetArray(localPQgen,&pqgenarr);CHKERRQ(ierr);
  
  ierr = DMGetLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ps->networkdm,pflow->X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ps->networkdm,pflow->X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = VecGetArrayRead(localX,&xarr);CHKERRQ(ierr);

  for(i=0; i < ps->nbranch; i++) {
    line = &ps->line[i];
    line->pf = line->qf = line->pt = line->qt = 0.0;
  }

  for(i=0; i < ps->nbus; i++) {
    PetscScalar theta, Vm;
    PSBUS bus;
    PetscInt loc,j;
    PetscScalar Pload=0.,Qload=0.,Pshunt=0.,Qshunt=0.;

    bus = &ps->bus[i];

    if(bus->ide == ISOLATED_BUS) continue;

    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
    ierr = PSBUSIsGhosted(bus,&ghostbus);CHKERRQ(ierr);

    theta = xarr[loc]; Vm = xarr[loc+1];

    bus->va = theta*180.0/PETSC_PI;
    bus->vm = Vm;
    if(!ghostbus) {
      //      ierr = PetscPrintf(PETSC_COMM_SELF,"Rank %d Bus %d Type %d Vm %4.3f Va %4.3f\n",pflow->comm->rank,bus->bus_i,bus->ide,bus->vm,bus->va);CHKERRQ(ierr);
    }

    /* Load injections */
    for(j=0; j < bus->nload; j++) {
      PSLOAD load;
      
      ierr = PSBUSGetLoad(bus,j,&load);CHKERRQ(ierr);
      if(!load->status) continue;
      Pload  += load->pl;
      Qload  += load->ql;
    }
    
    /* Shunt injections */
    Pshunt = Vm*Vm*bus->gl;
    Qshunt = -Vm*Vm*bus->bl;

    if(!ghostbus) {
      pqgenarr[loc]   = Pload + Pshunt;
      pqgenarr[loc+1] = Qload + Qshunt;
    }

    PetscInt nconnlines,k;
    const PSLINE *connlines;
    /* Injection from the lines supporting the bus */
    ierr = PSBUSGetSupportingLines(bus,&nconnlines,&connlines);CHKERRQ(ierr);
    
    Pnet = Qnet = 0.0;
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
      
      const PSBUS *connbuses;
      PSBUS busf,bust;
      PetscInt locf,loct;
      PetscScalar thetaf,Vmf,thetat,Vmt,thetaft,thetatf;
      
      /* Get the connected buses to this line */
      ierr = PSLINEGetConnectedBuses(line,&connbuses);CHKERRQ(ierr);
      busf = connbuses[0];
      bust = connbuses[1];
      
      ierr = PSBUSGetVariableLocation(busf,&locf);CHKERRQ(ierr);
      ierr = PSBUSGetVariableLocation(bust,&loct);CHKERRQ(ierr);
      thetaf = xarr[locf]; Vmf = xarr[locf+1];
      thetat = xarr[loct]; Vmt = xarr[loct+1];
      thetaft = thetaf - thetat;
      thetatf = thetat - thetaf;
      
      if (bus == busf) {
	line->pf   = Gff*Vmf*Vmf + Vmf*Vmt*(Gft*cos(thetaft) + Bft*sin(thetaft));
	line->qf   = -Bff*Vmf*Vmf + Vmf*Vmt*(-Bft*cos(thetaft) + Gft*sin(thetaft));
	Pnet += line->pf;
	Qnet += line->qf;
      } else {
	line->pt   = Gtt*Vmt*Vmt + Vmt*Vmf*(Gtf*cos(thetatf) + Btf*sin(thetatf));
	line->qt   = -Btt*Vmt*Vmt + Vmt*Vmf*(-Btf*cos(thetatf) + Gtf*sin(thetatf));
	Pnet += line->pt;
	Qnet += line->qt;
      }
    }

    pqgenarr[loc]   += Pnet;
    pqgenarr[loc+1] += Qnet;
  }

  ierr = VecRestoreArrayRead(localX,&xarr);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);

  ierr = VecRestoreArray(localPQgen,&pqgenarr);CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(ps->networkdm,localPQgen,ADD_VALUES,PQgen);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(ps->networkdm,localPQgen,ADD_VALUES,PQgen);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ps->networkdm,&localPQgen);CHKERRQ(ierr);

  ierr = DMGetLocalVector(ps->networkdm,&localPQgen);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ps->networkdm,PQgen,INSERT_VALUES,localPQgen);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ps->networkdm,PQgen,INSERT_VALUES,localPQgen);CHKERRQ(ierr);

  ierr = VecGetArray(localPQgen,&pqgenarr);CHKERRQ(ierr);

  for(i=0; i < ps->nbus; i++) {
    PSBUS bus;
    PetscInt loc,k;
    PetscScalar temptotal,pgtotalnonlt,qgtotalonlt,MVAbasenonlt,diffqg;
    bus = &ps->bus[i];

    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);

    if(!pflow->split_gen_power_within_limits) {
      if(bus->ide == REF_BUS || bus->ide == PV_BUS) {
	for(k=0; k < bus->ngen; k++) {
	  PSGEN gen;
	  ierr = PSBUSGetGen(bus,k,&gen);
	  if(bus->ide == REF_BUS) {
	    gen->pg = pqgenarr[loc]*gen->mbase/bus->MVAbasetot;
	  }
	  gen->qg = pqgenarr[loc+1]*gen->mbase/bus->MVAbasetot;
	}
      }
    } else {
      /* Split powers adhering to real and reactive power limits */
      if(bus->ide == REF_BUS || bus->ide == PV_BUS) {
	temptotal=0;
	pgtotalnonlt = bus->Pgtot; // Net Pg for generators not on VAR limits on this bus. This will be needed later if some of the buses are VAR limited during splitting.
	qgtotalonlt = 0; // Net Qg for generators on Var limit on the bus
	MVAbasenonlt = bus->MVAbasetot; // Net MVAbase for generators not on VAR limits on this bus.
	
	if (bus->ngen==1){     // If there is only one generator, do not split, this avoids that generator being VAR limited at the splitting stage. VAR limits will still be an issue if the net VAR at the bus is beyond its net VAR limit.
	  PSGEN gen;
	  ierr = PSBUSGetGen(bus,0,&gen);CHKERRQ(ierr);
	  if(bus->ide == REF_BUS) {
	    gen->pg = pqgenarr[loc];
	  }
	  gen->qg = pqgenarr[loc+1];
	  temptotal += gen->qg;
	} else {
	  for(k=0; k < bus->ngen; k++) {
	    PSGEN gen;
	    ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
	    if(gen->status) {
	      if(bus->ide == REF_BUS) {
		if(bus->Pgtot!=0) {
		  gen->pg = pqgenarr[loc]*PetscAbsScalar(gen->pg)/bus->Pgtot;
		}
		else{gen->pg = pqgenarr[loc]*gen->mbase/bus->MVAbasetot;}
	      }
	      if(PetscAbsScalar(gen->qt - gen->qb) < PETSC_SMALL) gen->qg = gen->qt; /* Fixed VAR generator */
	      else {
		if(bus->Pgtot!=0){gen->qg = pqgenarr[loc+1]*PetscAbsScalar(gen->pg)/bus->Pgtot; /* Proportional distribution */}
		else {
		  gen->qg = pqgenarr[loc+1]*gen->mbase/bus->MVAbasetot;
		}	 
		if (gen->qg > gen->qt) { // In PSSE while splitting the VARS on the generators, VAR limits are respected.
		  gen->qg = gen->qt;
		  pgtotalnonlt -= PetscAbsScalar(gen->pg); //Substract this PG from net PG not on limit
		  qgtotalonlt += gen->qg;    // Add this Qg to net Qg on limit
		  MVAbasenonlt -= gen->mbase; // Substract this MVAbase from net MVAbase not on limit
		}    
		else if (gen->qg<gen->qb) {
		  gen->qg = gen->qb;
		  pgtotalnonlt -= PetscAbsScalar(gen->pg);
		  qgtotalonlt += gen->qg;
		  MVAbasenonlt -= gen->mbase;
		}   
	      }
	    }
	    temptotal += gen->qg; // Calculate the total Qg on the bus, since some of the generators might be VAR limited, need to check if the total Qg is still the same as the total QG before splitting.
	  }
	}
	diffqg=pqgenarr[loc+1]-temptotal;  // Calculate the difference between net QG before and after splitting. 
	if (PetscAbsScalar(diffqg)>0.0001){ // If there is a difference between the net Qg before and after splitting, distribute the difference on generators  which are not on VAR limits proportional to their PG output (or MVAbase if net PG output of generators not on VAR limits is zero.)
	  for(k=0; k < bus->ngen; k++) {
	    PSGEN gen;
	    ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
	    if(PetscAbsScalar(gen->qg-gen->qt)>PETSC_SMALL && PetscAbsScalar(gen->qg-gen->qb)>PETSC_SMALL) {
	      if (PetscAbsScalar(pgtotalnonlt-0)>PETSC_SMALL) {
		gen->qg = (pqgenarr[loc+1]-qgtotalonlt)*PetscAbsScalar(gen->pg)/pgtotalnonlt;
	      }
	      else {
		gen->qg = (pqgenarr[loc+1]-qgtotalonlt)*gen->mbase/MVAbasenonlt;
	      }
	    }
	  }
	}
      }
    }
  }

#if defined PFLOW_DISPLAY_SNES_CONFIG
  /* Print to stdout */
  ierr = PetscPrintf(pflow->comm->type,"=============================================================\n");CHKERRQ(ierr);
  ierr = PetscPrintf(pflow->comm->type,"\tNonlinear solver configuration\n");CHKERRQ(ierr);
  ierr = PetscPrintf(pflow->comm->type,"=============================================================\n");CHKERRQ(ierr);
  ierr = SNESView(pflow->snes,0);CHKERRQ(ierr);
  MPI_Barrier(pflow->comm->type);
  ierr = PetscPrintf(pflow->comm->type,"\n\n");CHKERRQ(ierr);
#endif

#if defined PFLOW_DISPLAY_RESULTS
  ierr = PetscPrintf(pflow->comm->type,"=============================================================\n");CHKERRQ(ierr);
  ierr = PetscPrintf(pflow->comm->type,"\t\tBus data follows\n");CHKERRQ(ierr);
  ierr = PetscPrintf(pflow->comm->type,"=============================================================\n");CHKERRQ(ierr);
  ierr = PetscPrintf(pflow->comm->type,"%-6s %-10s %-4s %-7s %-7s %-7s %-6s %-7s %-7s\n","Rank","Bus #","Type","Pd","Qd","Pg","Qg","Vm","Va");CHKERRQ(ierr);
  MPI_Barrier(pflow->comm->type);

  for(i=0;i < ps->nbus; i++) {
    PSBUS bus;
    PSLOAD load;
    PetscInt loc;
    PetscScalar Pd,Qd;
    PetscBool   ghostbus;

    bus = &ps->bus[i];

    ierr = PSBUSIsGhosted(bus,&ghostbus);CHKERRQ(ierr);
    if(ghostbus) continue;

    Pd = Qd = 0.0;
    if(bus->nload) {
      ierr = PSBUSGetLoad(bus,0,&load);CHKERRQ(ierr);
      Pd = load->pl*ps->MVAbase;
      Qd = load->ql*ps->MVAbase;
    }

    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF,"%-6d %-6d %-4d %7.2f %7.2f %7.2f %7.2f %7.3f %7.3f\n",pflow->comm->rank,bus->bus_i,bus->ide,Pd,Qd,pqgenarr[loc]*ps->MVAbase,pqgenarr[loc+1]*ps->MVAbase,bus->vm,bus->va);CHKERRQ(ierr);
  }
  MPI_Barrier(pflow->comm->type);
  ierr = PetscPrintf(pflow->comm->type,"\n");CHKERRQ(ierr);
  ierr = PetscPrintf(pflow->comm->type,"=============================================================\n");CHKERRQ(ierr);
  ierr = PetscPrintf(pflow->comm->type,"\t\tBranch data follows\n");CHKERRQ(ierr);
  ierr = PetscPrintf(pflow->comm->type,"=============================================================\n");CHKERRQ(ierr);
  ierr = PetscPrintf(pflow->comm->type,"%-6s %-6s %-11s %-8s %-8s %-8s %-8s\n","Rank","From","To","Pf","Qf","Pt","Qt");CHKERRQ(ierr);
  MPI_Barrier(pflow->comm->type);
  for(i=0; i < ps->nbranch; i++) {
    PSLINE line;
    line = &ps->line[i];
    ierr = PetscPrintf(PETSC_COMM_SELF,"%-6d %-6d %-6d %8.2f %8.2f %8.2f %8.2f\n",pflow->comm->rank,line->fbus,line->tbus,line->pf*ps->MVAbase,line->qf*ps->MVAbase,line->pt*ps->MVAbase,line->qt*ps->MVAbase);CHKERRQ(ierr);
  }
  MPI_Barrier(pflow->comm->type);
  ierr = PetscPrintf(pflow->comm->type,"\n");CHKERRQ(ierr);
  ierr = PetscPrintf(pflow->comm->type,"=============================================================\n");CHKERRQ(ierr);
  ierr = PetscPrintf(pflow->comm->type,"\t\t Generator data follows\n");CHKERRQ(ierr);
  ierr = PetscPrintf(pflow->comm->type,"=============================================================\n");CHKERRQ(ierr);
  ierr = PetscPrintf(pflow->comm->type,"%-6s %-6s %-6s %-6s %-6s %-6s %-6s %-6s %-6s %-6s %-6s %-6s\n","Rank","Bus #","Pg","Qg","Qmin","Qmax","Qmintot","Qrange","Qgbus","Pgtotal","MVAbasetotal","Generator MVAbase");CHKERRQ(ierr);
  MPI_Barrier(pflow->comm->type);
  for(i=0; i < ps->nbus; i++) {
    PSBUS bus;
    PetscInt loc,k;
    
    bus = &ps->bus[i];
    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
    ierr = PSBUSIsGhosted(bus,&ghostbus);CHKERRQ(ierr);
    if(ghostbus) continue;

    for(k=0; k < bus->ngen; k++) {
      PSGEN gen;
      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_SELF,"%-6d %-6d %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n",pflow->comm->rank,bus->bus_i,gen->pg*ps->MVAbase,gen->qg*ps->MVAbase, gen->qb*ps->MVAbase, gen->qt*ps->MVAbase,bus->qmintot*ps->MVAbase,bus->qrange*ps->MVAbase,pqgenarr[loc+1]*ps->MVAbase,bus->Pgtot*ps->MVAbase,bus->MVAbasetot,gen->mbase);CHKERRQ(ierr);
    }
  }
  MPI_Barrier(pflow->comm->type);
#endif
  ierr = VecRestoreArray(PQgen,&pqgenarr);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ps->networkdm,&localPQgen);CHKERRQ(ierr);

  ierr = VecDestroy(&PQgen);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  PFLOWSetUp - Sets up a power flow application object

  Input Parameters:
. pflow - the PFLOW object

  Notes:
  This routine sets up the PFLOW object and the underlying PS object. It
  also distributes the PS object when used in parallel.
*/
PetscErrorCode PFLOWSetUp(PFLOW pflow)
{
 PetscErrorCode ierr;
  PS             ps=pflow->ps;

  PetscFunctionBegin;

  /* Set up PS object */
  ierr = PSSetUp(ps);CHKERRQ(ierr);

  /* Create the solution vector and the Jacobian matrix */
  ierr = PSCreateGlobalVector(pflow->ps,&pflow->X);CHKERRQ(ierr);
  ierr = PSCreateMatrix(pflow->ps,&pflow->Jac);CHKERRQ(ierr);

  /* Associate the DM object in PS with the SNES solver */
  ierr = SNESSetDM(pflow->snes,pflow->ps->networkdm);CHKERRQ(ierr);

  /* Set the prefix for this SNES.. All SNES runtime options will need to have the prefix "-pflow_" */
  ierr = SNESSetOptionsPrefix(pflow->snes,"pflow_");CHKERRQ(ierr);

  /* Set the function and Jacobian routines */
  ierr = SNESSetFunction(pflow->snes,NULL,PFLOWFunction,(void*)pflow);CHKERRQ(ierr);
  ierr = SNESSetJacobian(pflow->snes,pflow->Jac,pflow->Jac,PFLOWJacobian,(void*)pflow);CHKERRQ(ierr);

  ierr = SNESSetFromOptions(pflow->snes);CHKERRQ(ierr);

  pflow->setupcalled = PETSC_TRUE;
  pflow->converged = PETSC_FALSE;
  PetscFunctionReturn(0);
}

/*
  PFLOWSetLineStatus - Sets the status (ON/OFF) for a line

  Input Parameters
+ pflow - The PFLOW object
. fbus  - from bus
. tbus  - to bus
. id    - line id
- status - line status (0 = OFF, 1 = ON)

Notes: PFLOWSetUp() must be called before calling PFLOWSetLineStatus()

*/
PetscErrorCode PFLOWSetLineStatus(PFLOW pflow, PetscInt fbus, PetscInt tbus, const char* id,PetscInt status)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PSSetLineStatus(pflow->ps,fbus,tbus,id,status);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  PFLOWSetGenStatus - Sets the status (ON/OFF) for a generator

  Input Parameters
+ pflow - The PFLOW object
. gbus  - generator bus
. gid   - generator id
- status - generator status (0 = OFF, 1 = ON)

Notes: PFLOWSetUp() must be called before calling PFLOWSetGenStatus()

*/
PetscErrorCode PFLOWSetGenStatus(PFLOW pflow, PetscInt gbus, const char* gid,PetscInt status)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PSSetGenStatus(pflow->ps,gbus,gid,status);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  PFLOWGetBusVoltage - Returns the voltage magnitude and angle for the given bus

  Input Parameters
+ pflow - The PFLOW object
. busnum  - bus number
. Vm   - bus voltage magnitude (pu)
. Va   - bus voltage angle (degrees)
- found - TRUE if bus found, else FALSE

Notes: PFLOWSetUp() must be called before calling PFLOWGetBusVoltage()

*/
PetscErrorCode PFLOWGetBusVoltage(PFLOW pflow, PetscInt busnum, PetscScalar *Vm, PetscScalar *Va, PetscBool *found)
{
  PetscInt       intbusnum;
  PS             ps=pflow->ps;
  PSBUS          bus;

  PetscFunctionBegin;

  /* Convert from external to internal bus number */
  intbusnum = ps->busext2intmap[busnum];
  if(intbusnum == -1) {
    *found = PETSC_FALSE;
    PetscFunctionReturn(0);
  } else {
    *found = PETSC_TRUE;
    bus = &ps->bus[intbusnum];
    *Vm = bus->vm;
    *Va = bus->va;
  }

  PetscFunctionReturn(0);
}

/*
  PFLOWSetBusVoltage - Sets the voltage magnitude and angle for the given bus

  Input Parameters
+ pflow - The PFLOW object
. busnum  - bus number
. Vm   - bus voltage magnitude (pu)
- Va   - bus voltage angle (degrees)


Notes: PFLOWSetUp() must be called before calling PFLOWSetBusVoltage()

*/
PetscErrorCode PFLOWSetBusVoltage(PFLOW pflow, PetscInt busnum, PetscScalar Vm, PetscScalar Va)
{
  PetscErrorCode ierr;
  PetscInt       intbusnum;
  PS             ps=pflow->ps;
  PSBUS          bus;
  PSGEN          gen;
  PetscInt       k;

  PetscFunctionBegin;

  /* Convert from external to internal bus number */
  intbusnum = ps->busext2intmap[busnum];
  if(intbusnum != -1) {
    bus = &ps->bus[intbusnum];
    bus->vm = Vm;
    bus->va = Va;
    if(bus->ngen) {
      for(k=0; k < bus->ngen; k++) {
	ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
	gen->vs = Vm;
      }
    }
  }

  PetscFunctionReturn(0);
}


/*
  PFLOWGetGenDispatch - Gets the real and reactive power output for generator at a given bus

  Input Parameters
+ pflow   - The PFLOW object
. busnum  - bus number
. gennum  - generate number
. Pg   - real power output for the generator (MW)
. Qg   - reactive power output for the generator (MVAr)
- found - PETSC_TRUE if generator found, PETSC_FALSE otherwise

Notes: PFLOWSetUp() must be called before calling PFLOWSetBusVoltage()
*/
PetscErrorCode PFLOWGetGenDispatch(PFLOW pflow, PetscInt busnum, PetscInt gennum,PetscScalar* Pg, PetscScalar* Qg,PetscBool *found)
{
  PetscErrorCode ierr;
  PetscInt       intbusnum;
  PS             ps=pflow->ps;
  PSBUS          bus;
  PSGEN          gen;

  PetscFunctionBegin;

  /* Convert from external to internal bus number */
  intbusnum = ps->busext2intmap[busnum];
  if(intbusnum != -1) {
    bus = &ps->bus[intbusnum];
    if(!bus->ngen) {
      SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_SUP,"No generators on bus %d",busnum);
    }
    if(bus->ide == ISOLATED_BUS) {
      SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_SUP,"Bus %d is isolated",busnum);
    }
    ierr = PSBUSGetGen(bus,gennum,&gen);CHKERRQ(ierr);
    *Pg = gen->pg*ps->MVAbase;
    *Qg = gen->qg*ps->MVAbase;
    *found = PETSC_TRUE;
  } else {
    *found = PETSC_FALSE;
  }

  PetscFunctionReturn(0);
}

/*
  PFLOWSetGenDispatch - Sets the real and reactive power output for generator at a given bus

  Input Parameters
+ pflow   - The PFLOW object
. busnum  - bus number
. gennum  - generate number
. Pg   - real power outpu (MW)
- Qg   - reactive power output (MVAr)


Notes: PFLOWSetUp() must be called before calling PFLOWSetBusVoltage()
       If the bus is a reference bus, the real and reactive power set will
       not be adhered to on power flow solution.
       If the bus is a PV bus, the reactive power set will be not be adhered
       to on power flow solution.

*/
PetscErrorCode PFLOWSetGenDispatch(PFLOW pflow, PetscInt busnum, PetscInt gennum,PetscScalar Pg, PetscScalar Qg)
{
  PetscErrorCode ierr;
  PetscInt       intbusnum;
  PS             ps=pflow->ps;
  PSBUS          bus;
  PSGEN          gen;

  PetscFunctionBegin;

  /* Convert from external to internal bus number */
  intbusnum = ps->busext2intmap[busnum];
  if(intbusnum != -1) {
    bus = &ps->bus[intbusnum];
    if(!bus->ngen) {
      SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_SUP,"No generators on bus %d",busnum);
    }
    if(bus->ide == ISOLATED_BUS) {
      SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_SUP,"Bus %d is isolated",busnum);
    }
    ierr = PSBUSGetGen(bus,gennum,&gen);CHKERRQ(ierr);
    gen->pg = Pg/ps->MVAbase;
    gen->qg = Qg/ps->MVAbase;
  }

  PetscFunctionReturn(0);
}

/*
  PFLOWSetLoadPower - Sets the real and reactive power load for the given bus

  Input Parameters
+ pflow - The PFLOW object
. busnum  - bus number
. Pd   - real load (MW)
- Qd   - reactive load (MVAr)


Notes: PFLOWSetUp() must be called before calling PFLOWSetBusVoltage()

*/
PetscErrorCode PFLOWSetLoadPower(PFLOW pflow, PetscInt busnum, PetscScalar Pd, PetscScalar Qd)
{
  PetscErrorCode ierr;
  PetscInt       intbusnum;
  PS             ps=pflow->ps;
  PSBUS          bus;
  PSLOAD         load;

  PetscFunctionBegin;

  /* Convert from external to internal bus number */
  intbusnum = ps->busext2intmap[busnum];
  if(intbusnum != -1) {
    bus = &ps->bus[intbusnum];
    ierr = PSBUSGetLoad(bus,0,&load);CHKERRQ(ierr);
    load->pl = Pd/ps->MVAbase;
    load->ql = Qd/ps->MVAbase;
  }

  PetscFunctionReturn(0);
}

/*
  PFLOWGetLoadPower - Gets the real and reactive power load for the given bus

  Input Parameters
+ pflow - The PFLOW object
. busnum  - bus number
. Pd   - real load (MW)
. Qd   - reactive load (MVAr)
- found - TRUE if bus found, else FALSE


Notes: PFLOWSetUp() must be called before calling PFLOWGetBusVoltage()

*/
PetscErrorCode PFLOWGetLoadPower(PFLOW pflow, PetscInt busnum, PetscScalar *Pd, PetscScalar *Qd, PetscBool *found)
{
  PetscErrorCode ierr;
  PetscInt       intbusnum;
  PS             ps=pflow->ps;
  PSBUS          bus;
  PSLOAD         load;

  PetscFunctionBegin;

  /* Convert from external to internal bus number */
  intbusnum = ps->busext2intmap[busnum];
  if(intbusnum != -1) {
    bus = &ps->bus[intbusnum];
    ierr = PSBUSGetLoad(bus,0,&load);CHKERRQ(ierr);
    *Pd = load->pl*ps->MVAbase;
    *Qd = load->ql*ps->MVAbase;
    *found = PETSC_TRUE;
  }
  
  PetscFunctionReturn(0);
}

/*
  PFLOWAddBusShunt - Sets the real and reactive shunt for the given bus

  Input Parameters
+ pflow - The PFLOW object
. busnum  - bus number
. Gs   - shunt condutance
- Bs   - shunt susceptance


Notes: PFLOWSetUp() must be called before calling PFLOWAddBusShunt()

*/
PetscErrorCode PFLOWAddBusShunt(PFLOW pflow, PetscInt busnum, PetscScalar Gs, PetscScalar Bs)
{
  PetscErrorCode ierr;
  PetscInt       intbusnum;
  PS             ps=pflow->ps;
  PSBUS          bus;
  
  PetscFunctionBegin;
  
  /* Convert from external to internal bus number */
  intbusnum = ps->busext2intmap[busnum];
  if(intbusnum != -1) {
    bus = &ps->bus[intbusnum];
    ierr = PSBUSAddShunt(bus,Gs/ps->MVAbase,Bs/ps->MVAbase);CHKERRQ(ierr);
  }
  
  PetscFunctionReturn(0);
}
