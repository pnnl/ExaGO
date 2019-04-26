#include <dyn.h>
#include <private/psimpl.h>
#include <private/pflowimpl.h>
#include <private/dynimpl.h>
#include <dirent.h>
#include <petscdmnetwork.h>
#include <petsc/private/petscimpl.h>

#include "dynloadmodels/dynzip.h"
#include "dyngenmodels/dyncv.h"

PetscErrorCode DYNGetCostFunction(DYN dyn,PetscReal t,PetscReal *costfcn,PetscInt ncostfun);

static PetscErrorCode DYNSetUpActivePowerIS(DYN dyn,IS *ispg)
{
  IS             is;
  PS             ps=dyn->ps;
  PSBUS          bus;
  PSGEN          gen;
  PetscInt       i,k,paramloc=0,ctr=0,*idx;
  PetscBool      ghostbus;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = VecGetOwnershipRange(dyn->Xparam,&paramloc,NULL);CHKERRQ(ierr);
  ierr = PetscCalloc1(dyn->ncostfcns,&idx);CHKERRQ(ierr);
  for (i=0; i<ps->nbus; i++) {
    bus = &ps->bus[i];
    paramloc += 2; // VA and VM
    ierr = PSBUSIsGhosted(bus,&ghostbus);CHKERRQ(ierr);
    if (ghostbus) continue;
    for (k=0; k<bus->ngen; k++) {
      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      if (!gen->status) {
        paramloc += 2;
        continue;
      }
      idx[ctr++] = paramloc; // PG location
      paramloc += 2;
    }
  }
  ierr = ISCreateGeneral(dyn->comm->type,dyn->ncostfcns,idx,PETSC_COPY_VALUES,&is);CHKERRQ(ierr);
  ierr = PetscFree(idx);CHKERRQ(ierr);
  *ispg = is;
  PetscFunctionReturn(0);
}

PetscErrorCode DYNGetCostFunctionParametersFromOptions(DYN dyn)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscOptionsGetReal(NULL,NULL,"-dyn_cost_scale",&dyn->scal,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-dyn_cost_exp",&dyn->exp,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-dyn_cost_eta",&dyn->eta,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-dyn_freq_max",&dyn->freq_max,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-dyn_freq_min",&dyn->freq_min,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-dyn_sum_frequency_costfun",&dyn->sum_freq_costfun,NULL);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  DYNSetPostStepCallback - Sets a callback function for every step of DYN

  Input Parameters:
+ dyn - the dynamic simulation object
- poststepcallback - the callback function

  Notes: The prototype of the callback function should be PetscErrorCode PostStepFunction(DYN dyn)
*/
PetscErrorCode DYNSetPostStepCallback(DYN dyn,PetscErrorCode (*poststepcallback)(DYN dyn))
{

  PetscFunctionBegin;
  if (dyn->npoststepfcns == MAXPOSTSTEPFCNS) SETERRQ1(dyn->comm->type,0,"Exceeded max. post step functions allowed %D",MAXPOSTSTEPFCNS);
  dyn->poststepfcn[dyn->npoststepfcns++] = poststepcallback;
  PetscFunctionReturn(0);
}

/*
  DYNSetPreStepCallback - Sets a callback function called before every step

  Input Parameters:
+ dyn - the dynamic simulation object
- prestepcallback - the callback function

  Notes: The prototype of the callback function should be PetscErrorCode PreStepFunction(DYN dyn)
*/
PetscErrorCode DYNSetPreStepCallback(DYN dyn,PetscErrorCode (*prestepcallback)(DYN dyn))
{

  PetscFunctionBegin;
  if (dyn->nprestepfcns == MAXPOSTSTEPFCNS) SETERRQ1(dyn->comm->type,0,"Exceeded max. pre step functions allowed %D",MAXPOSTSTEPFCNS);
  dyn->prestepfcn[dyn->nprestepfcns++] = prestepcallback;
  PetscFunctionReturn(0);
}

/*
  DYNGetTS - Gets the PETSc time-stepping solver object set with DYN

Input parameters:
. dyn - the DYN object

Output parameters:
. ts - the TS object

  Notes: This returns a borrowed reference to ts. The object should not be destroyed.
*/
PetscErrorCode DYNGetTS(DYN dyn, TS *ts)
{
  PetscFunctionBegin;
  if(!dyn->setupcalled) SETERRQ(dyn->comm->type,0,"DYNSetUp() must be called before calling DYNGetTS");
  *ts = dyn->ts;
  PetscFunctionReturn(0);
}

/*
  DYNGetTS - Gets the power system object PS set with DYN

Input parameters:
. dyn - the DYN object

Output parameters:
. ps - the power system object

  Notes: This returns a borrowed reference to ps. The object should not be destroyed.
*/
PetscErrorCode DYNGetPS(DYN dyn, PS *ps)
{
  PetscFunctionBegin;
  if(!dyn->setupcalled) SETERRQ(dyn->comm->type,0,"DYNSetUp() must be called before calling DYNGetTS");
  *ps = dyn->ps;
  PetscFunctionReturn(0);
}

/*
  DYNViewer - Saves the bus voltage magnitudes and generator frequencies at each time-step
*/
PetscErrorCode DYNViewer(DYN dyn)
{
  PetscErrorCode ierr;
  TS             ts=dyn->ts;
  PetscReal      t;
  PS             ps;
  Vec            X,localX;
  char           volt_file[PETSC_MAX_PATH_LEN];
  char           freq_file[PETSC_MAX_PATH_LEN];
  char           log_file[PETSC_MAX_PATH_LEN];
  //  char           mvaflow_file[PETSC_MAX_PATH_LEN];
  PetscInt       stepnum;
  PetscViewer    vltview,frqview,logview;
  //  PetscViewer    mvaflowview;
  const PetscScalar *xarr,*xdyn;
  PetscInt       i,loc,k;
  PSBUS          bus;
  PSGEN          gen;
  DYNGenModel    dyngen;
  PetscBool      ghostbus;
  PetscScalar    VD,VQ,Vm;
  PetscScalar    frq;
  
  PetscFunctionBegin;

  ierr = TSGetTime(ts,&t);CHKERRQ(ierr);
  ierr = TSGetSolution(ts,&X);CHKERRQ(ierr);
  ierr = TSGetStepNumber(ts,&stepnum);CHKERRQ(ierr);

  ierr = PetscSNPrintf(volt_file,sizeof(volt_file),"%s/volt-%06d.out",dyn->viewer_dir,stepnum);CHKERRQ(ierr);
  ierr = PetscSNPrintf(freq_file,sizeof(freq_file),"%s/freq-%06d.out",dyn->viewer_dir,stepnum);CHKERRQ(ierr);
  ierr = PetscSNPrintf(log_file,sizeof(log_file),"%s/log-%06d.out",dyn->viewer_dir,stepnum);CHKERRQ(ierr);
  //  ierr = PetscSNPrintf(mvaflow_file,sizeof(mvaflow_file),"%s/mvaflow-%06d.out",dyn->viewer_dir,stepnum);CHKERRQ(ierr);

  //  ierr = PetscSNPrintf(volt_file,sizeof(volt_file),"dyn-output/volt-%06d.out",stepnum);CHKERRQ(ierr);
  //  ierr = PetscSNPrintf(freq_file,sizeof(freq_file),"dyn-output/freq-%06d.out",stepnum);CHKERRQ(ierr);

  ierr = PetscViewerASCIIOpen(dyn->comm->type,volt_file,&vltview);CHKERRQ(ierr);
  ierr = PetscViewerASCIIOpen(dyn->comm->type,freq_file,&frqview);CHKERRQ(ierr);
  ierr = PetscViewerASCIIOpen(dyn->comm->type,log_file,&logview);CHKERRQ(ierr);
  //  ierr = PetscViewerASCIIOpen(dyn->comm->type,mvaflow_file,&mvaflowview);CHKERRQ(ierr);
  
  ierr = PetscViewerASCIIPrintf(vltview,"Time   Bus    Vm\n");CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(frqview,"Time   Bus    GenId  Frequency\n");CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(logview,"Time   EventTag  Bus/Branch number  Message\n");CHKERRQ(ierr);
  //  ierr = PetscViewerASCIIPrintf(mvaflowview,"Time From To  Id Sf St\n");CHKERRQ(ierr);

  ps = dyn->ps;

  ierr = DMGetLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);

  ierr = VecGetArrayRead(localX,&xarr);CHKERRQ(ierr);

  for(i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSIsGhosted(bus,&ghostbus);CHKERRQ(ierr);
    if (ghostbus) continue;
    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);

    VD = xarr[loc]; VQ = xarr[loc+1];
    Vm = PetscSqrtScalar(VD*VD+VQ*VQ);

    /* Format => t   bus_num  Vm */
    ierr = PetscViewerASCIIPrintf(vltview,"%4.4f %-6d %3.4f\n",t,bus->bus_i,Vm);CHKERRQ(ierr);

    /* Format => t  tag bus number message */
    if(Vm > bus->Vmax) {
      ierr = PetscViewerASCIIPrintf(logview,"%4.4f,OVERVOLTAGE VIOLATION,bus,%-6d,voltage magnitude %3.4f over limit %3.4f\n",t,bus->bus_i,Vm,bus->Vmax);CHKERRQ(ierr);
    } else if(Vm < bus->Vmin) {
      ierr = PetscViewerASCIIPrintf(logview,"%4.4f,UNDERVOLTAGE VIOLATION,bus,%-6d,voltage magnitude %3.4f under limit %3.4f\n",t,bus->bus_i,Vm,bus->Vmin);CHKERRQ(ierr);
    }

    xdyn = xarr + loc;
    for(k = 0; k < bus->ngen; k++) {
      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      /*       if(!gen->status) continue; */

      ierr = PSGENGetDYNGen(gen,&dyngen);CHKERRQ(ierr);
      ierr = DYNGenModelGetFrequency(dyngen,t, xdyn, &frq);CHKERRQ(ierr);

      ierr = PetscViewerASCIIPrintf(frqview,"%4.4f %-6d %-6s %3.4f\n",t,gen->bus_i,gen->id, frq);CHKERRQ(ierr);
    }
  }

  ierr = VecRestoreArrayRead(localX,&xarr);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);

  ierr = PetscViewerDestroy(&vltview);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&frqview);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&logview);CHKERRQ(ierr);
  //  ierr = PetscViewerDestroy(&mvaflowview);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode DYNUserMonitor(TS ts,PetscInt stepnum,PetscReal t,Vec X,void *ctx)
{
  DYN               dyn=(DYN)ctx;
  PetscInt          i;
  PetscErrorCode    ierr;
  const PetscScalar *x;
  FILE              *fp1,*fp2;

  PetscFunctionBegin;

  fp1 = fopen("freq_violation.data","a");
  ierr = PetscFPrintf(PETSC_COMM_WORLD,fp1,"%g ",t);CHKERRQ(ierr);

  fp2 = fopen("voltage.data","a");
  ierr = PetscFPrintf(PETSC_COMM_WORLD,fp2,"%g ",t);CHKERRQ(ierr);

  Vec localX;
  PS       ps=dyn->ps;
  PetscBool ghostbus;

  ierr = DMGetLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = VecGetArrayRead(localX,&x);CHKERRQ(ierr);

  for(i=0; i < ps->nbus; i++) {
    PetscInt    loc;
    PSBUS       bus;
    PetscInt    k;

    bus = &ps->bus[i];
    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
    ierr = PSBUSIsGhosted(bus,&ghostbus);CHKERRQ(ierr);
    if (ghostbus) continue;

    if(bus->ide == ISOLATED_BUS) continue;

    const PetscScalar *xdyn;
    xdyn = x+loc;
    /* Generator frequency deviation */
    for(k=0; k < bus->ngen; k++) {
      PSGEN gen;
      DYNGenModel dyngen;
      DYNExcModel dynexc;

      PetscInt    excloc;
      PetscScalar frequency=0;

      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      ierr = PSGENGetDYNGen(gen,&dyngen);CHKERRQ(ierr);

      //if(!gen->status) continue;
      if(gen->status) {
        /* Get Machine frequency */
        ierr = DYNGenModelGetFrequency(dyngen,t,xdyn,&frequency);
        /* Print frequency to file */
        ierr = PetscFPrintf(PETSC_COMM_WORLD,fp1,"%g ",frequency);CHKERRQ(ierr);
        if(dyngen->dynexc) {
          ierr = PSGENGetDYNExc(gen,&dynexc);CHKERRQ(ierr);
          ierr = DYNExcModelGetFirstVariableLocation(dynexc,&excloc);
          /* Print index */
          //printf("loc=%d exloc=%d\n",loc,excloc);
          /* Print VR to file */
          ierr = PetscFPrintf(PETSC_COMM_WORLD,fp2,"%g ",xdyn[excloc+1]);CHKERRQ(ierr);
        }
      }else {
        /* Print zero if the generator if off */
        ierr = PetscFPrintf(PETSC_COMM_WORLD,fp1,"0 ");CHKERRQ(ierr);
        ierr = PetscFPrintf(PETSC_COMM_WORLD,fp2,"0 ");CHKERRQ(ierr);
      }
    }
  }

  ierr = VecRestoreArrayRead(localX,&x);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);

  PetscScalar *c;
  ierr = PetscCalloc1(dyn->ncostfcns,&c);CHKERRQ(ierr);
  ierr = DYNGetCostFunction(dyn,t,c,dyn->ncostfcns);CHKERRQ(ierr);
  for(i=0;i<dyn->ncostfcns;i++) {
    ierr = PetscFPrintf(PETSC_COMM_WORLD,fp1,"%g ",c[i]);CHKERRQ(ierr);
  }
  ierr = PetscFree(c);CHKERRQ(ierr);

  ierr = PetscFPrintf(PETSC_COMM_WORLD,fp1,"\n");CHKERRQ(ierr);
  ierr = PetscFPrintf(PETSC_COMM_WORLD,fp2,"\n");CHKERRQ(ierr);
  fclose(fp1);
  fclose(fp2);
#ifdef FWDSA
  if(dyn->useforward) {
    PetscInt j;
    fp1 = fopen("exc_fwdsensi.data", "a");
    ierr = PetscFPrintf(PETSC_COMM_WORLD,fp1,"%24.10f ",t);CHKERRQ(ierr);
    for(j=0;j<dyn->nparams;j++) {
      const PetscScalar *s;
      ierr = VecGetArrayRead(dyn->sp[j],&s);CHKERRQ(ierr);
      for(i=0; i < ps->nbus; i++) {
        PetscInt    loc;
        PSBUS       bus;
        PetscInt    k;

        bus = &ps->bus[i];
        ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
        ierr = PSBUSIsGhosted(bus,&ghostbus);CHKERRQ(ierr);
        if (ghostbus) continue;

        if(bus->ide == ISOLATED_BUS) continue;

        const PetscScalar *sdyn;
        sdyn = s+loc;
        /* Generator frequency deviation */
        for(k=0; k < bus->ngen; k++) {
          PSGEN gen;
          DYNExcModel dynexc;
          DYNGenModel dyngen;
          PetscInt    excloc;
          ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
          ierr = PSGENGetDYNGen(gen,&dyngen);CHKERRQ(ierr);
          //if(!gen->status) continue;
          if(gen->status) {
            if(dyngen->dynexc) {
              ierr = PSGENGetDYNExc(gen,&dynexc);CHKERRQ(ierr);
              ierr = DYNExcModelGetFirstVariableLocation(dynexc,&excloc);
              /* Print VR to file */
              ierr = PetscFPrintf(PETSC_COMM_WORLD,fp1,"%24.15e ",sdyn[excloc+1]);CHKERRQ(ierr);
            }
          }else {
            /* Print zero if the generator is off */
            ierr = PetscFPrintf(PETSC_COMM_WORLD,fp1,"%24.15e ",sdyn[excloc+1]);CHKERRQ(ierr);
          }
        }
      }
      ierr = VecRestoreArrayRead(dyn->sp[j],&s);CHKERRQ(ierr);
    }
    ierr = PetscFPrintf(PETSC_COMM_WORLD,fp1,"\n");CHKERRQ(ierr);
    fclose(fp1);
  }
#endif
  PetscFunctionReturn(0);
}

PetscErrorCode CheckAlgebraicPart(DYN dyn,Vec v)
{
  PetscInt xdynsize;
  Vec      subv;
  IS       isalge;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = VecGetSize(dyn->X,&xdynsize);CHKERRQ(ierr);
  ierr = ISComplement(dyn->isdiff,0,xdynsize,&isalge);CHKERRQ(ierr);
  ierr = VecGetSubVector(v,isalge,&subv);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"checking algebraic part\n");CHKERRQ(ierr);
  ierr = VecView(subv,0);CHKERRQ(ierr);
  ierr = VecRestoreSubVector(v,isalge, &subv);
  PetscFunctionReturn(0);
}

/*
  DYNEventMonitor_Machine - Event monitoring routine for machines
*/
PetscErrorCode DYNEventMonitor_Machine(DYNEvent dynevent, PetscReal t, Vec X, PetscScalar *fval)
{
  PetscErrorCode ierr;
  DYN            dyn=(DYN)dynevent->data;
  PS             ps=dyn->ps;
  PetscInt       nbus=ps->nbus,i,k;
  Vec            localX;
  const PetscScalar *xarr,*xdyn;
  PetscBool      ghostbus;
  PSBUS          bus;
  PSGEN          gen;
  PSLOAD         load;
  DYNGenModel    dyngen;
  DYNExcModel    dynexc;
  DYNTurbgovModel dynturbgov;
  DYNStabModel    dynstab;
  DYNLoadModel    dynload;
  PetscScalar    VD,VQ;
  PetscInt       loc,startloc;
  PetscScalar    *fdynval = fval;

  PetscFunctionBegin;

  ierr = DMGetLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = VecGetArrayRead(localX,&xarr);CHKERRQ(ierr);

  for(i=0; i < nbus; i++) {
    bus = &ps->bus[i];

    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
    ierr = PSBUSIsGhosted(bus,&ghostbus);CHKERRQ(ierr);
    if (ghostbus) continue;

    VD = xarr[loc];  VQ = xarr[loc+1];

    if (bus->ngen) {
      for(k=0; k < bus->ngen; k++) {
        ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);

        if(!gen->status) {
	  /* Check if the generator was initially ON and then turned OFF */
	  if(gen->initial_status) {
	    ierr = PSGENGetDYNGen(gen,&dyngen);CHKERRQ(ierr);
	    fdynval += dyngen->numMonitors;

	    if(dyngen->dynexc) {
	      ierr = PSGENGetDYNExc(gen,&dynexc);CHKERRQ(ierr);
	      fdynval += dynexc->numMonitors;
	    }
	  
	    if(dyngen->dynturbgov) {
	      ierr = PSGENGetDYNTurbgov(gen,&dynturbgov);CHKERRQ(ierr);
	      fdynval += dynturbgov->numMonitors;
	    }

	    if(dyngen->dynexc && dynexc->dynstab) {
	      ierr = PSGENGetDYNStab(gen,&dynstab);CHKERRQ(ierr);
	      fdynval += dynstab->numMonitors;
	    }
	  }
	  continue;
	}

        ierr = PSGENGetDYNGen(gen,&dyngen);CHKERRQ(ierr);

        startloc = loc;
        xdyn = xarr+startloc;

        /* Call generator event function */
        if(dyngen->numMonitors) {
          ierr = (*dyngen->eventfcn)(dyngen,t,VD,VQ,(PetscScalar*)xdyn,fdynval);CHKERRQ(ierr);
          fdynval += dyngen->numMonitors;
        }

        if(dyngen->dynexc) {
          ierr = PSGENGetDYNExc(gen,&dynexc);CHKERRQ(ierr);
          /* Call exciter event function */
          ierr = (*dynexc->eventfcn)(dynexc,t,VD,VQ,(PetscScalar*)xdyn,fdynval);CHKERRQ(ierr);
          fdynval += dynexc->numMonitors;
        }

        if(dyngen->dynturbgov) {
          ierr = PSGENGetDYNTurbgov(gen,&dynturbgov);CHKERRQ(ierr);
          /* Call turbine governor event function */
          ierr = (*dynturbgov->eventfcn)(dynturbgov,t,VD,VQ,(PetscScalar*)xdyn,fdynval);CHKERRQ(ierr);
          fdynval += dynturbgov->numMonitors;
        }

        if(dyngen->dynexc && dynexc->dynstab) {
          ierr = PSGENGetDYNStab(gen,&dynstab);CHKERRQ(ierr);
          /* Call turbine governor event function */
          ierr = (*dynstab->eventfcn)(dynstab,t,VD,VQ,(PetscScalar*)xdyn,fdynval);CHKERRQ(ierr);
          fdynval += dynstab->numMonitors;
        }
      }
    }

    if(bus->nload) {
      for(k=0; k < bus->nload; k++) {
	ierr = PSBUSGetLoad(bus,k,&load);CHKERRQ(ierr);
	
	if(!load->status) continue;
	
	ierr = PSLOADGetDYNLoad(load,&dynload);CHKERRQ(ierr);
	
	startloc = loc;
	xdyn = xarr+startloc;
	
	/* Call load event function */
	if(dynload->numMonitors) {
	  ierr = (*dynload->eventfcn)(dynload,t,VD,VQ,(PetscScalar*)xdyn,fdynval);CHKERRQ(ierr);
          fdynval += dynload->numMonitors;
	}
      }
    }
  }

  ierr = VecRestoreArrayRead(localX,&xarr);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  Note: X is a local vector */
/*
  DYNEventPostFunction_Machine - Post event handling routine for machines
*/
PetscErrorCode DYNEventPostFunction_Machine(DYNEvent dynevent, PetscInt nmonitors, PetscInt monidx[], PetscReal t, Vec localX, PetscBool forwardsolve, PetscBool *solve_alg)
{
  PetscErrorCode ierr;
  DYN            dyn=(DYN)dynevent->data;
  PS             ps=dyn->ps;
  PetscInt       nbus=ps->nbus,i,k,j;
  PetscScalar    *xarr,*xdyn;
  PetscBool      ghostbus;
  PSBUS          bus;
  PSGEN          gen;
  PSLOAD         load;
  DYNGenModel    dyngen;
  DYNExcModel    dynexc;
  DYNTurbgovModel dynturbgov;
  DYNStabModel    dynstab;
  DYNLoadModel   dynload;
  PetscScalar    VD,VQ;
  PetscInt       loc,startloc,nmon=0,ev_list[MAXDYNEVENTSLOCATEDATSTEP],ctr=0;
  PetscBool      solve_alg_mac=PETSC_FALSE;

  PetscFunctionBegin;

  ierr = VecGetArray(localX,&xarr);CHKERRQ(ierr);

  for(i=0; i < nbus; i++) {
    if(ctr == nmonitors) break;

    bus = &ps->bus[i];
    ierr = PSBUSIsGhosted(bus,&ghostbus);CHKERRQ(ierr);
    if (ghostbus) continue;

    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
    VD = xarr[loc];  VQ = xarr[loc+1];

    if (bus->ngen) {
      for(k=0; k < bus->ngen; k++) {
        ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
        if(!gen->status) {
	  /* Check if the generator was initially ON and then turned OFF */
	  if(gen->initial_status) {
	    ierr = PSGENGetDYNGen(gen,&dyngen);CHKERRQ(ierr);
	    nmon += dyngen->numMonitors;

	    if(dyngen->dynexc) {
	      ierr = PSGENGetDYNExc(gen,&dynexc);CHKERRQ(ierr);
	      nmon += dynexc->numMonitors;
	    }
	  
	    if(dyngen->dynturbgov) {
	      ierr = PSGENGetDYNTurbgov(gen,&dynturbgov);CHKERRQ(ierr);
	      nmon += dynturbgov->numMonitors;
	    }

	    if(dyngen->dynexc && dynexc->dynstab) {
	      ierr = PSGENGetDYNStab(gen,&dynstab);CHKERRQ(ierr);
	      nmon += dynstab->numMonitors;
	    }
	  }
	  continue;
	}
	
        ierr = PSGENGetDYNGen(gen,&dyngen);CHKERRQ(ierr);

        startloc = loc;
        xdyn = xarr+startloc;

        /* Post event handling for generator model */
        j = 0;
        while(monidx[ctr] >= nmon && monidx[ctr] < nmon + dyngen->numMonitors) {
          ev_list[j++] = monidx[ctr++]-nmon;
          if(ctr == nmonitors) break;
        }
        /* Call generator post event function */
        if(j && dyngen->posteventfcn) {
          ierr = (*dyngen->posteventfcn)(dyngen,j,ev_list,t,VD,VQ,(PetscScalar*)xdyn,forwardsolve,&solve_alg_mac);CHKERRQ(ierr);
          *solve_alg = *solve_alg || solve_alg_mac;
        }
        nmon += dyngen->numMonitors;

        if(dyngen->dynexc) {
          ierr = PSGENGetDYNExc(gen,&dynexc);CHKERRQ(ierr);
          j = 0;
          while(monidx[ctr] >= nmon && monidx[ctr] < nmon + dynexc->numMonitors) {
            ev_list[j++] = monidx[ctr++]-nmon;
            if(ctr == nmonitors) break;
          }
          /* Call exciter event function */
          if(j && dynexc->posteventfcn) {
            ierr = (*dynexc->posteventfcn)(dynexc,j,ev_list,t,VD,VQ,(PetscScalar*)xdyn,forwardsolve,solve_alg);CHKERRQ(ierr);
          }
          nmon += dynexc->numMonitors;
        }

        /* Turbine governor model */
        if(dyngen->dynturbgov) {
          ierr = PSGENGetDYNTurbgov(gen,&dynturbgov);CHKERRQ(ierr);
          j = 0;
          while(monidx[ctr] >= nmon && monidx[ctr] < nmon + dynturbgov->numMonitors) {
            ev_list[j++] = monidx[ctr++]-nmon;
            if(ctr == nmonitors) break;
          }
          /* Call turbine governor event function */
          if(j && dynturbgov->posteventfcn) {
            ierr = (*dynturbgov->posteventfcn)(dynturbgov,j,ev_list,t,VD,VQ,(PetscScalar*)xdyn,forwardsolve,solve_alg);CHKERRQ(ierr);
          }
          nmon += dynturbgov->numMonitors;
        }

        /* Stabilizer model */
        if(dyngen->dynexc && dynexc->dynstab) {
          ierr = PSGENGetDYNStab(gen,&dynstab);CHKERRQ(ierr);
          j = 0;
          while(monidx[ctr] >= nmon && monidx[ctr] < nmon + dynstab->numMonitors) {
            ev_list[j++] = monidx[ctr++]-nmon;
            if(ctr == nmonitors) break;
          }
          /* Call stabilizer event function */
          if(j && dynstab->posteventfcn) {
            ierr = (*dynstab->posteventfcn)(dynstab,j,ev_list,t,VD,VQ,(PetscScalar*)xdyn,forwardsolve,solve_alg);CHKERRQ(ierr);
          }
          nmon += dynstab->numMonitors;
        }
      }
    }

    if (bus->nload) {
      for(k=0; k < bus->nload; k++) {
        ierr = PSBUSGetLoad(bus,k,&load);CHKERRQ(ierr);
        if(!load->status) continue;

        ierr = PSLOADGetDYNLoad(load,&dynload);CHKERRQ(ierr);

        startloc = loc;
        xdyn = xarr+startloc;

        /* Post event handling for load model */
        j = 0;
        while(monidx[ctr] >= nmon && monidx[ctr] < nmon + dynload->numMonitors) {
          ev_list[j++] = monidx[ctr++]-nmon;
          if(ctr == nmonitors) break;
        }
        /* Call load post event function */
        if(j && dynload->posteventfcn) {
          ierr = (*dynload->posteventfcn)(dynload,j,ev_list,t,VD,VQ,(PetscScalar*)xdyn,forwardsolve,&solve_alg_mac);CHKERRQ(ierr);
          *solve_alg = *solve_alg || solve_alg_mac;
        }
        nmon += dynload->numMonitors;
      }
    }
  }

  ierr = VecRestoreArray(localX,&xarr);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  DYNEventPostDGamma_Machine - Post event routine for machines
*/
PetscErrorCode DYNEventPostDGamma_Machine(DYNEvent dynevent, PetscInt nmonitors, PetscInt monidx[], PetscReal t, Vec localX, PetscBool forwardsolve)
{
  PetscErrorCode ierr;
  DYN            dyn=(DYN)dynevent->data;
  PS             ps=dyn->ps;
  PetscInt       nbus=ps->nbus,i,k,j;
  Vec            localdGdX;
  const PetscScalar *xarr,*xdyn;
  PetscScalar    *dgdxarr;
  PetscScalar    *dgdx=NULL,*dgdp=NULL;
  PetscBool      ghostbus;
  PSBUS          bus;
  PSGEN          gen;
  DYNGenModel    dyngen;
  DYNExcModel    dynexc;
  DYNTurbgovModel dynturbgov;
  DYNStabModel    dynstab;
  PetscScalar    VD,VQ;
  PetscInt       loc,startloc,nmon=0,ev_list[MAXDYNEVENTSLOCATEDATSTEP],ctr=0;
  PetscInt       paramloc;

  PetscFunctionBegin;
  ierr = VecGetArray(dyn->dGammadP,&dgdp);CHKERRQ(ierr);

  ierr = VecGetArrayRead(localX,&xarr);CHKERRQ(ierr);

  /* get local dGammadX */
  ierr = DMGetLocalVector(ps->networkdm,&localdGdX);CHKERRQ(ierr);
  ierr = VecSet(dyn->dGammadX,0.0);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ps->networkdm,dyn->dGammadX,INSERT_VALUES,localdGdX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ps->networkdm,dyn->dGammadX,INSERT_VALUES,localdGdX);CHKERRQ(ierr);
  ierr = VecGetArray(localdGdX,&dgdxarr);CHKERRQ(ierr);

  ierr = VecGetOwnershipRange(dyn->Xparam,&paramloc,NULL);CHKERRQ(ierr);

  for(i=0; i < nbus; i++) {
    PetscInt PGloc,QGloc,VAloc,VMloc,ctr2=0;

    if(ctr == nmonitors) break;
    bus = &ps->bus[i];
    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr); /* ? */
    ierr = PSBUSIsGhosted(bus,&ghostbus);CHKERRQ(ierr);
    if (ghostbus) continue;

    VAloc = paramloc; VMloc = paramloc+1; paramloc += 2;
    VD = xarr[loc]; VQ = xarr[loc+1];
    for(k=0; k < bus->ngen; k++) {
      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);

      PGloc = paramloc + ctr2;
      QGloc = paramloc + ctr2 + 1;
      ctr2 += 2; // keep PG and QG even when the generator is off
      if(!gen->status) continue;
      ierr = PSGENGetDYNGen(gen,&dyngen);CHKERRQ(ierr);
      startloc = loc;
      xdyn = xarr+startloc;
      dgdx = dgdxarr+startloc;

      /* Generator gamma function handling */
      j = 0;
      while(nmon >= monidx[ctr] && monidx[ctr] < nmon + dyngen->numMonitors) {
        ev_list[j++] = monidx[ctr++]-nmon;
        if(ctr == nmonitors) break;
      }
      /* Call generator post event dgamma function */
      if(j && dyngen->posteventdgamma) {
        ierr = (*dyngen->posteventdgamma)(dyngen,j,ev_list,t,VD,VQ,(PetscScalar*)xdyn,dgdx,dgdp,VAloc,VMloc,PGloc,QGloc);CHKERRQ(ierr);
      }
      nmon += dyngen->numMonitors;

      if(dyngen->dynexc) {
        ierr = PSGENGetDYNExc(gen,&dynexc);CHKERRQ(ierr);
        j = 0;
        while(nmon >= monidx[ctr] && monidx[ctr] < nmon + dynexc->numMonitors) {
          ev_list[j++] = monidx[ctr++]-nmon;
          if(ctr == nmonitors) break;
        }
        /* Call exciter post event dgamma function */
        if(j && dynexc->posteventdgamma) {
          ierr = (*dynexc->posteventdgamma)(dynexc,j,ev_list,t,VD,VQ,(PetscScalar*)xdyn,dgdx,dgdp,VAloc,VMloc,PGloc,QGloc);CHKERRQ(ierr);
        }
        nmon += dynexc->numMonitors;
      }

      if(dyngen->dynturbgov) {
        ierr = PSGENGetDYNTurbgov(gen,&dynturbgov);CHKERRQ(ierr);
        j = 0;
        while(nmon >= monidx[ctr] && monidx[ctr] < nmon + dynturbgov->numMonitors) {
          ev_list[j++] = monidx[ctr++]-nmon;
          if(ctr == nmonitors) break;
        }
        /* Call turbine governor posteventdgamma function */
        if(j && dynturbgov->posteventdgamma) {
          ierr = (*dynturbgov->posteventdgamma)(dynturbgov,j,ev_list,t,VD,VQ,(PetscScalar*)xdyn,dgdx,dgdp,VAloc,VMloc,PGloc,QGloc);CHKERRQ(ierr);
        }
        nmon += dynturbgov->numMonitors;
      }

      if(dyngen->dynexc && dynexc->dynstab) {
        ierr = PSGENGetDYNStab(gen,&dynstab);CHKERRQ(ierr);
        j = 0;
        while(nmon >= monidx[ctr] && monidx[ctr] < nmon + dynstab->numMonitors) {
          ev_list[j++] = monidx[ctr++]-nmon;
          if(ctr == nmonitors) break;
        }
        /* Call stabilizer posteventdgamma function */
        if(j && dynstab->posteventdgamma) {
          ierr = (*dynstab->posteventdgamma)(dynstab,j,ev_list,t,VD,VQ,(PetscScalar*)xdyn,dgdx,dgdp,VAloc,VMloc,PGloc,QGloc);CHKERRQ(ierr);
        }
        nmon += dynstab->numMonitors;
      }
    }
    paramloc += ctr2;
  }

  ierr = VecRestoreArrayRead(localX,&xarr);CHKERRQ(ierr);
  ierr = VecRestoreArray(localdGdX,&dgdxarr);CHKERRQ(ierr);
  ierr = VecRestoreArray(dyn->dGammadP,&dgdp);CHKERRQ(ierr);

  ierr = DMLocalToGlobalBegin(ps->networkdm,localdGdX,INSERT_VALUES,dyn->dGammadX);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(ps->networkdm,localdGdX,INSERT_VALUES,dyn->dGammadX);CHKERRQ(ierr);

  ierr = DMRestoreLocalVector(ps->networkdm,&localdGdX);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  DYNCheckandSetMachineLimits - Checks if device limits have been exceeded and, if so, set its limit flags

  Input Parameters:
+ dyn - the DYN object
. t   - the current time
- X   - the solution vector at current time

*/
PetscErrorCode DYNCheckandSetMachineLimits(DYN dyn, PetscReal t, Vec X,PetscBool forwardsolve)
{
  PetscErrorCode ierr;
  DYNEvent       dynevent;
  PS             ps=dyn->ps;
  PetscInt       nbus=ps->nbus,i,k,j;
  Vec            localX;
  const PetscScalar *xarr,*xdyn;
  PetscBool      ghostbus;
  PSBUS          bus;
  PSGEN          gen;
  PSLOAD         load;
  DYNGenModel    dyngen;
  DYNExcModel    dynexc;
  DYNTurbgovModel dynturbgov;
  DYNStabModel    dynstab;
  DYNLoadModel   dynload;
  PetscScalar    VD,VQ;
  PetscInt       loc,startloc;
  PetscScalar    fdynval[MAXDYNEVENTSLOCATEDATSTEP];
  PetscInt       ev_list[MAXDYNEVENTSLOCATEDATSTEP];
  PetscScalar    *fprv=dyn->feventtmp;
  PetscInt       nlimshit=0;
  PetscBool     solve_alg;

  PetscFunctionBegin;

  for(i=0; i < dyn->nevents; i++) {
    dynevent = dyn->events[dyn->event_idx[i]];
    if(dynevent->eventmonitor == DYNEventMonitor_Machine) break;
    fprv += dynevent->numMonitors;
  }

  ierr = DMGetLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = VecGetArrayRead(localX,&xarr);CHKERRQ(ierr);

  for(i=0; i < nbus; i++) {
    bus = &ps->bus[i];

    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
    ierr = PSBUSIsGhosted(bus,&ghostbus);CHKERRQ(ierr);
    if (ghostbus) continue;

    VD = xarr[loc];  VQ = xarr[loc+1];

    if (bus->ngen) {
      for(k=0; k < bus->ngen; k++) {
        ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);

        if(!gen->status) {
	  /* Check if the generator was initially ON and then turned OFF */
	  if(gen->initial_status) {
	    ierr = PSGENGetDYNGen(gen,&dyngen);CHKERRQ(ierr);
	    fprv += dyngen->numMonitors;

	    if(dyngen->dynexc) {
	      ierr = PSGENGetDYNExc(gen,&dynexc);CHKERRQ(ierr);
	      fprv += dynexc->numMonitors;
	    }
	  
	    if(dyngen->dynturbgov) {
	      ierr = PSGENGetDYNTurbgov(gen,&dynturbgov);CHKERRQ(ierr);
	      fprv += dynturbgov->numMonitors;
	    }

	    if(dyngen->dynexc && dynexc->dynstab) {
	      ierr = PSGENGetDYNStab(gen,&dynstab);CHKERRQ(ierr);
	      fprv += dynstab->numMonitors;
	    }
	  }
	  continue;
	}
	
        ierr = PSGENGetDYNGen(gen,&dyngen);CHKERRQ(ierr);

	fprv += dyngen->numMonitors; /* No limits checked for generators */
	
        startloc = loc;
        xdyn = xarr+startloc;

	nlimshit = 0;
	if(dyngen->dynexc) {
	  ierr = PSGENGetDYNExc(gen,&dynexc);CHKERRQ(ierr);
	  if(dynexc->numMonitors) {
	    /* Call exciter event function */
	    ierr = (*dynexc->eventfcn)(dynexc,t,VD,VQ,(PetscScalar*)xdyn,fdynval);CHKERRQ(ierr);
	    for(j=0; j < dynexc->numMonitors; j++) {
	      if(PetscSign(PetscRealPart(fprv[j])) != PetscSign(PetscRealPart(fdynval[j]))) ev_list[nlimshit++] = j;
	    }
	    if(nlimshit && dynexc->posteventfcn) {
	      /* Call post event function */
	      ierr = (*dynexc->posteventfcn)(dynexc,nlimshit,ev_list,t,VD,VQ,(PetscScalar*)xdyn,forwardsolve,&solve_alg);CHKERRQ(ierr);
	    }
	       
	    /* Condition for event sign change and calling post event function */
	    fprv += dynexc->numMonitors;
	  }
	}
	nlimshit = 0;

	if(dyngen->dynturbgov) {
	  ierr = PSGENGetDYNTurbgov(gen,&dynturbgov);CHKERRQ(ierr);
	  if(dynturbgov->numMonitors) {
	    /* Call turbine governor event function */
	    ierr = (*dynturbgov->eventfcn)(dynturbgov,t,VD,VQ,(PetscScalar*)xdyn,fdynval);CHKERRQ(ierr);

	    for(j=0; j < dynturbgov->numMonitors; j++) {
	      if(PetscSign(PetscRealPart(fprv[j])) != PetscSign(PetscRealPart(fdynval[j]))) ev_list[nlimshit++] = j;
	    }
	    if(nlimshit && dynturbgov->posteventfcn) {
	      /* Call post event function */
	      ierr = (*dynturbgov->posteventfcn)(dynturbgov,nlimshit,ev_list,t,VD,VQ,(PetscScalar*)xdyn,forwardsolve,&solve_alg);CHKERRQ(ierr);
	    }
	    
	    /* Condition for event sign change and calling post event function */
	    fprv += dynturbgov->numMonitors;
	  }
	}

	nlimshit=0;

	if(dyngen->dynexc && dynexc->dynstab) {
	  ierr = PSGENGetDYNStab(gen,&dynstab);CHKERRQ(ierr);
	  if(dynstab->numMonitors) {
	    /* Call stabilizer event function */
	    ierr = (*dynstab->eventfcn)(dynstab,t,VD,VQ,(PetscScalar*)xdyn,fdynval);CHKERRQ(ierr);
	    
	    for(j=0; j < dynstab->numMonitors; j++) {
	      if(PetscSign(PetscRealPart(fprv[j])) != PetscSign(PetscRealPart(fdynval[j]))) ev_list[nlimshit++] = j;
	    }
	    if(nlimshit && dynstab->posteventfcn) {
	      /* Call post event function */
	      ierr = (*dynstab->posteventfcn)(dynstab,nlimshit,ev_list,t,VD,VQ,(PetscScalar*)xdyn,forwardsolve,&solve_alg);CHKERRQ(ierr);
	    }
	    
	    /* Condition for event sign change and calling post event function */
	    fprv += dynstab->numMonitors;
	  }
	}
      }
    }

    if (bus->nload) {
      for(k=0; k < bus->nload; k++) {
        ierr = PSBUSGetLoad(bus,k,&load);CHKERRQ(ierr);

        if(!load->status) continue;

        ierr = PSLOADGetDYNLoad(load,&dynload);CHKERRQ(ierr);

        startloc = loc;
        xdyn = xarr+startloc;

	nlimshit = 0;

	ierr = PSLOADGetDYNLoad(load,&dynload);CHKERRQ(ierr);
	if(dynload->numMonitors) {
	  /* Call load event function */
	  ierr = (*dynload->eventfcn)(dynload,t,VD,VQ,(PetscScalar*)xdyn,fdynval);CHKERRQ(ierr);
	  for(j=0; j < dynload->numMonitors; j++) {
	    if(PetscSign(PetscRealPart(fprv[j])) != PetscSign(PetscRealPart(fdynval[j]))) ev_list[nlimshit++] = j;
	  }
	  if(nlimshit && dynload->posteventfcn) {
	    /* Call post event function */
	    ierr = (*dynload->posteventfcn)(dynload,nlimshit,ev_list,t,VD,VQ,(PetscScalar*)xdyn,forwardsolve,&solve_alg);CHKERRQ(ierr);
	  }
	       
	  /* Condition for event sign change and calling post event function */
	  fprv += dynload->numMonitors;
	}
      }
    }
  }

  ierr = VecRestoreArrayRead(localX,&xarr);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  DYNSetMachineEvents - Sets the events on machines and their control circuity

  Input Parameters:
. dyn - the DYN object

  Notes: This routine adds events for gens, excs, and other controllers.
*/
PetscErrorCode DYNSetMachineEvents(DYN dyn)
{
  PetscErrorCode ierr;
  PS             ps=dyn->ps;
  PetscInt       nbus=ps->nbus,i,k,j,ctr=0;
  PetscBool      ghostbus;
  PSBUS          bus;
  PSGEN          gen;
  PSLOAD         load;
  DYNGenModel    dyngen;
  DYNExcModel    dynexc;
  DYNTurbgovModel dynturbgov;
  DYNStabModel    dynstab;
  DYNLoadModel    dynload;
  DYNEvent       dynevent;
  PetscInt       nmonitors=0,Nmonitors=0;

  PetscFunctionBegin;

  for(i=0; i < nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSIsGhosted(bus,&ghostbus);CHKERRQ(ierr);
    if (ghostbus) continue;

    if (bus->ngen) {
      for(k=0; k < bus->ngen; k++) {
        ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);

        if(!gen->initial_status) continue;

        ierr = PSGENGetDYNGen(gen,&dyngen);CHKERRQ(ierr);

        ierr = DYNGenModelSetEventMonitor(dyngen);CHKERRQ(ierr);
        nmonitors += dyngen->numMonitors;

        /* Exciter Model */
        if(dyngen->dynexc) {
          ierr = PSGENGetDYNExc(gen,&dynexc);CHKERRQ(ierr);
          /* Set up exciter events */
          ierr = DYNExcModelSetEventMonitor(dynexc);CHKERRQ(ierr);
          nmonitors += dynexc->numMonitors;
        }

        /* Turbine governor model */
        if(dyngen->dynturbgov) {
          ierr = PSGENGetDYNTurbgov(gen,&dynturbgov);CHKERRQ(ierr);
          /* Set up turbine governor events */
          ierr = DYNTurbgovModelSetEventMonitor(dynturbgov);CHKERRQ(ierr);
          nmonitors += dynturbgov->numMonitors;
        }

        /* Stabilizer model */
        if(dyngen->dynexc && dynexc->dynstab) {
          ierr = PSGENGetDYNStab(gen,&dynstab);CHKERRQ(ierr);
          /* Set up stabilizer events */
          ierr = DYNStabModelSetEventMonitor(dynstab);CHKERRQ(ierr);
          nmonitors += dynstab->numMonitors;
        }
      }
    }

    if(bus->nload) {
      for(k=0; k < bus->nload; k++) {
	ierr = PSBUSGetLoad(bus,k,&load);CHKERRQ(ierr);

	if(!load->status) continue;

	ierr = PSLOADGetDYNLoad(load,&dynload);CHKERRQ(ierr);

	ierr = DYNLoadModelSetEventMonitor(dynload);CHKERRQ(ierr);
	nmonitors += dynload->numMonitors;
      }
    }
  }

  ierr = MPI_Allreduce(&nmonitors,&Nmonitors,1,MPIU_INT,MPI_SUM,dyn->comm->type);CHKERRQ(ierr);
  if(Nmonitors) {
    ierr = DYNEventCreate(&dyn->events[dyn->Nevents]);CHKERRQ(ierr);
    dynevent = dyn->events[dyn->Nevents];
    dynevent->data = (void*)dyn;
    dynevent->eventmonitor = DYNEventMonitor_Machine;
    dynevent->posteventfunction = DYNEventPostFunction_Machine;
    dynevent->posteventdgamma = DYNEventPostDGamma_Machine;
    dynevent->destroy = NULL;

    dynevent->numMonitors = nmonitors;
    dyn->event_idx[dyn->nevents] = dyn->Nevents;
    dyn->toteventmonitors += dynevent->numMonitors;

    ierr = PetscCalloc1(nmonitors,&dynevent->direction);CHKERRQ(ierr);
    ierr = PetscCalloc1(nmonitors,&dynevent->terminate);CHKERRQ(ierr);

    /* Copy over terminate and direction */
    for(i=0; i < nbus; i++) {
      bus = &ps->bus[i];
      ierr = PSBUSIsGhosted(bus,&ghostbus);CHKERRQ(ierr);
      if (ghostbus) continue;

      if (bus->ngen) {
	for(k=0; k < bus->ngen; k++) {
	  ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
	  
	  if(!gen->status) continue;

	  ierr = PSGENGetDYNGen(gen,&dyngen);CHKERRQ(ierr);
	  for(j=0; j < dyngen->numMonitors; j++, ctr++) {
	    dynevent->direction[ctr] = dyngen->direction[j];
	    dynevent->terminate[ctr] = dyngen->terminate[j];
	  }

	  /* Exciter model */
	  if(dyngen->dynexc) {
	    ierr = PSGENGetDYNExc(gen,&dynexc);CHKERRQ(ierr);
	    
	    for(j=0; j < dynexc->numMonitors; j++,ctr++) {
	      dynevent->direction[ctr] = dynexc->direction[j];
	      dynevent->terminate[ctr] = dynexc->terminate[j];
	    }
	  } 

	  /* Turbine governor model */
	  if(dyngen->dynturbgov) {
	    ierr = PSGENGetDYNTurbgov(gen,&dynturbgov);CHKERRQ(ierr);
	    
	    for(j=0; j < dynturbgov->numMonitors; j++,ctr++) {
	      dynevent->direction[ctr] = dynturbgov->direction[j];
	      dynevent->terminate[ctr] = dynturbgov->terminate[j];
	    }
	  } 

	  /* Stabilizer model */
	  if(dyngen->dynexc && dynexc->dynstab) {
	    ierr = PSGENGetDYNStab(gen,&dynstab);CHKERRQ(ierr);
	    
	    for(j=0; j < dynstab->numMonitors; j++,ctr++) {
	      dynevent->direction[ctr] = dynstab->direction[j];
	      dynevent->terminate[ctr] = dynstab->terminate[j];
	    }
	  } 
	}
      }

      if(bus->nload) {
	for(k=0; k < bus->nload; k++) {
	  ierr = PSBUSGetLoad(bus,k,&load);CHKERRQ(ierr);
	  
	  if(!load->status) continue;
	  
	  ierr = PSLOADGetDYNLoad(load,&dynload);CHKERRQ(ierr);
	  for(j=0; j < dynload->numMonitors; j++, ctr++) {
	    dynevent->direction[ctr] = dynload->direction[j];
	    dynevent->terminate[ctr] = dynload->terminate[j];
	  }
	}
      }
    }
    dyn->nevents++;
    dyn->Nevents++;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode DYNAdjointMonitor(TS ts,PetscInt stepnum,PetscReal t,Vec X,PetscInt ncost,Vec *lambda,Vec *mu,void *ctx)
{
  PetscInt       i;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  for (i=0;i<ncost;i++) {
    ierr = VecView(lambda[i],0);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*
  DYNSetGenStatus - Sets the status (ON/OFF) for a generator

  Input Parameters
+ dyn - The DYN object
. gbus  - generator bus
. gid   - generator id
- status - generator status (0 = OFF, 1 = ON)

Notes: DYNSetUp() must be called before calling DYNSetGenStatus()

*/
PetscErrorCode DYNSetGenStatus(DYN dyn, PetscInt gbus, const char* gid,PetscInt status)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PSSetGenStatus(dyn->ps,gbus,gid,status);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


/*
  DYNSetGenDispatchandStatus - Sets the generator status and dispatch

  Input Parameters:
+ dyn - the DYN object
. busnum - the bus number
. gennum - generator number
. status - the generator status (0 for OFF, 1 for ON)
. pg - active power dispatch in MW
- qg - reactive power dispatch in MVAr
*/
PetscErrorCode DYNSetGenDispatchandStatus(DYN dyn, PetscInt busnum, PetscInt gennum, PetscInt status, PetscScalar Pg, PetscScalar Qg)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PSSetGenDispatchandStatus(dyn->ps,busnum,gennum,status,Pg,Qg);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

void swap_dm(DM *dm1, DM *dm2)
{
  DM temp = *dm1;
  *dm1 = *dm2;
  *dm2 = temp;
}

/*
  DYNSetUpInitPflow - Sets up the power flow solver object used in obtaining initial conditions.

  Note: This power flow solver object shares the underlying power system object and all the data. It
  is created at the end of DYNSetUp. So, it already has the distributed ps object. The only
  difference is that it uses a different PetscSection for storing the degrees of freedom that
  gets associated with the dmnetwork (and plex)
*/
PetscErrorCode DYNSetUpInitPflow(DYN dyn)
{
  PetscErrorCode ierr;
  PetscInt       vStart,vEnd,eStart,eEnd;
  DM             networkdm,plexdm;
  PetscInt       i;

  PetscFunctionBegin;

  networkdm = dyn->ps->networkdm;

  ierr = PFLOWCreate(dyn->comm->type,&dyn->initpflow);CHKERRQ(ierr);

  /* PFLOW creates a new PS object. Destroy it so that we can associate the
     PS from DYN with initpflow
  */
  ierr = PSDestroy(&dyn->initpflow->ps);CHKERRQ(ierr);

  dyn->initpflow->ps = dyn->ps;
  /* Increase the reference count for dyn->ps */
  ierr = PSIncreaseReferenceCount(dyn->ps);CHKERRQ(ierr);

  ierr = PetscSectionCreate(dyn->comm->type,&dyn->initpflowpsection);CHKERRQ(ierr);

  ierr = DMNetworkGetVertexRange(networkdm,&vStart,&vEnd);CHKERRQ(ierr);
  ierr = DMNetworkGetEdgeRange(networkdm,&eStart,&eEnd);CHKERRQ(ierr);

  ierr = PetscSectionSetChart(dyn->initpflowpsection,eStart,vEnd);CHKERRQ(ierr);

  for(i=vStart; i < vEnd; i++) {
    /* Two variables at each bus/vertex */
    ierr = PetscSectionSetDof(dyn->initpflowpsection,i,2);CHKERRQ(ierr);
  }
  ierr = PetscSectionSetUp(dyn->initpflowpsection);CHKERRQ(ierr);

  /* Clone DM to be used with initpflow */
  ierr = DMClone(networkdm,&dyn->initpflowdm);CHKERRQ(ierr);

  /* Set initpflowdm in dyn->ps->networkdm, the previous networkdm get stashed in dyn->initpflowdm */
  swap_dm(&dyn->ps->networkdm,&dyn->initpflowdm);
  networkdm = dyn->ps->networkdm;

  /* Get the plex dm */
  ierr = DMNetworkGetPlex(networkdm,&plexdm);CHKERRQ(ierr);

  /* Get default sections associated with this plex */
  ierr = DMGetDefaultSection(plexdm,&dyn->defaultsection);CHKERRQ(ierr);
  ierr = PetscObjectReference((PetscObject)dyn->defaultsection);CHKERRQ(ierr);

  ierr = DMGetDefaultGlobalSection(plexdm,&dyn->defaultglobalsection);CHKERRQ(ierr);
  ierr = PetscObjectReference((PetscObject)dyn->defaultglobalsection);CHKERRQ(ierr);

  /* Set the new section created for initial power flow */
  ierr = DMSetDefaultSection(plexdm,dyn->initpflowpsection);CHKERRQ(ierr);
  ierr = DMGetDefaultGlobalSection(plexdm,&dyn->initpflowpglobsection);CHKERRQ(ierr);

  /* Set up PFLOW object. Note pflow->ps will not be set up again as it has
     been already set up by dyn
  */
  ierr = PFLOWSetUp(dyn->initpflow);CHKERRQ(ierr);

  /*
     Update the ref. counts for init pflow sections so that they do not
     get destroyed when DMSetDefaultSection is called
  */
  // ierr = PetscObjectReference((PetscObject)dyn->initpflowpsection);CHKERRQ(ierr);
  ierr = PetscObjectReference((PetscObject)dyn->initpflowpglobsection);CHKERRQ(ierr);
  
  /* Reset the sections */
  ierr = DMSetDefaultSection(plexdm,dyn->defaultsection);CHKERRQ(ierr);
  ierr = DMSetDefaultGlobalSection(plexdm,dyn->defaultglobalsection);CHKERRQ(ierr);

  /* Reset dm */
  swap_dm(&dyn->ps->networkdm,&dyn->initpflowdm);

  PetscFunctionReturn(0);
}

PetscErrorCode PSSetEdgeandBusStartLoc(PS ps)
{
  PetscErrorCode ierr;
  PetscInt       i,vStart,vEnd,eStart,eEnd;

  PetscFunctionBegin;
  /* Reset the start locations for the buses and lines */
  ierr = DMNetworkGetVertexRange(ps->networkdm,&vStart,&vEnd);CHKERRQ(ierr);
  ierr = DMNetworkGetEdgeRange(ps->networkdm,&eStart,&eEnd);CHKERRQ(ierr);
  
  for(i=eStart; i < eEnd; i++) {
    ierr = DMNetworkGetVariableOffset(ps->networkdm,i,&ps->line[i].startloc);CHKERRQ(ierr);
    ierr = DMNetworkGetVariableGlobalOffset(ps->networkdm,i,&ps->line[i].startlocglob);CHKERRQ(ierr);
  }
  
  for(i=vStart; i < vEnd; i++) {
    ierr = DMNetworkGetVariableOffset(ps->networkdm,i,&ps->bus[i-vStart].startloc);CHKERRQ(ierr);
    ierr = DMNetworkGetVariableGlobalOffset(ps->networkdm,i,&ps->bus[i-vStart].startlocglob);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*
  DYNComputePrePflow - Computes the steady-state power flow for initializing the dynamics simulation

  Input Parameters:
. dyn - the DYN object
*/
PetscErrorCode DYNComputePrePflow(DYN dyn,PetscBool *converged)
{
  PetscErrorCode ierr;
  DM             plexdm;
  PFLOW          initpflow=dyn->initpflow;
  PS             ps=dyn->ps;

  PetscFunctionBegin;
  /* Set initpflowdm for solving the power flow */
  swap_dm(&dyn->ps->networkdm,&dyn->initpflowdm);

  ierr = DMNetworkGetPlex(dyn->ps->networkdm,&plexdm);CHKERRQ(ierr);

  /* Increase the ref. counts for the default sections so that they do not get
     destroyed when DMSetDefaultXXX is called with the initpflowxxx sections
  */
  ierr = PetscObjectReference((PetscObject)dyn->defaultsection);CHKERRQ(ierr);
  ierr = PetscObjectReference((PetscObject)dyn->defaultglobalsection);CHKERRQ(ierr);

  /* Set the new section created for initial power flow. */
  ierr = DMSetDefaultSection(plexdm,dyn->initpflowpsection);CHKERRQ(ierr);
  ierr = DMSetDefaultGlobalSection(plexdm,dyn->initpflowpglobsection);CHKERRQ(ierr);

  ierr = DMCreateDefaultSF(plexdm,dyn->initpflowpsection,dyn->initpflowpglobsection);CHKERRQ(ierr);

  /* Reset the edge and bus starting locations of variables */
  ierr = PSSetEdgeandBusStartLoc(ps);CHKERRQ(ierr);

  /* Solve */
  ierr = PFLOWSolve(initpflow);CHKERRQ(ierr);
  ierr = PFLOWConverged(initpflow,converged);CHKERRQ(ierr);

  /* Update bus and gen structs in pflow->ps */
  ierr = PFLOWPostSolve(initpflow);CHKERRQ(ierr);

  /*
     Update the ref. counts for init pflow sections so that they do not
     get destroyed when DMSetDefaultSection is called
  */
  ierr = PetscObjectReference((PetscObject)dyn->initpflowpsection);CHKERRQ(ierr);
  ierr = PetscObjectReference((PetscObject)dyn->initpflowpglobsection);CHKERRQ(ierr);

  /* Reset the sections */
  ierr = DMSetDefaultSection(plexdm,dyn->defaultsection);CHKERRQ(ierr);
  ierr = DMSetDefaultGlobalSection(plexdm,dyn->defaultglobalsection);CHKERRQ(ierr);

  /* Reset SF */
  ierr = DMCreateDefaultSF(plexdm,dyn->defaultsection,dyn->defaultglobalsection);CHKERRQ(ierr);

  /* Reset the bus and edge starting locations for variables */
  ierr = PSSetEdgeandBusStartLoc(ps);CHKERRQ(ierr);

  /* Reset dm */
  swap_dm(&dyn->initpflowdm,&dyn->ps->networkdm);
  PetscFunctionReturn(0);
}

PetscErrorCode DYNComputeSensiP2(DYN dyn)
{
  PetscErrorCode ierr;
  PS             ps=dyn->ps;
  PetscInt       i,k;
  PetscBool      ghostbus;
  PetscScalar    val;
  PSBUS          bus;
  PetscInt       paramloc,row;
  PetscInt       VA,VM,VD,VQ;

  PetscFunctionBegin;
  ierr = MatZeroEntries(dyn->ICp);CHKERRQ(ierr);
  ierr = VecGetOwnershipRange(dyn->Xparam,&paramloc,NULL);CHKERRQ(ierr);

  for (i=0; i<ps->nbus; i++) {
    PetscInt PGloc,QGloc,VAloc,VMloc,buslocglob; // should be global indices

    bus = &ps->bus[i];
    ierr = PSBUSIsGhosted(bus,&ghostbus);CHKERRQ(ierr);
    if (ghostbus) continue;

    VAloc = paramloc; VMloc = paramloc+1; paramloc += 2;
    ierr = PSBUSGetVariableGlobalLocation(bus,&buslocglob);CHKERRQ(ierr);
    if (bus->ngen) {
      PSGEN    gen;
      PetscInt atleastonegenon=0;

      for(k=0; k<bus->ngen; k++) {
        ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
        atleastonegenon = atleastonegenon+gen->status;
        if(!gen->status) continue;
        VM = gen->vs;
        bus->vm = VM;
      }
      if(!atleastonegenon) VM = bus->vm;
    } else VM = bus->vm;
    VA = bus->va*PETSC_PI/180.0;
    VD = VM*PetscCosScalar(VA);
    VQ = VM*PetscSinScalar(VA);

    row  = buslocglob;
    val  = -VM*PetscSinScalar(VA);
    ierr = SetMatrixValues(dyn->ICp,1,&row,1,&VAloc,&val);CHKERRQ(ierr);
    val  = PetscCosScalar(VA);
    ierr = SetMatrixValues(dyn->ICp,1,&row,1,&VMloc,&val);CHKERRQ(ierr);
    row  = buslocglob+1;
    val  = VM*PetscCosScalar(VA);
    ierr = SetMatrixValues(dyn->ICp,1,&row,1,&VAloc,&val);CHKERRQ(ierr);
    val  = PetscSinScalar(VA);
    ierr = SetMatrixValues(dyn->ICp,1,&row,1,&VMloc,&val);CHKERRQ(ierr);
    for(k=0; k<bus->ngen; k++) {
      PSGEN           gen;
      DYNGenModel     dyngen;
      DYNExcModel     dynexc;
      DYNTurbgovModel dynturbgov;
      DYNStabModel    dynstab;
      PetscInt        dynlocglob;

      PGloc = paramloc + 2*k;
      QGloc = paramloc + 2*k + 1;
      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      if(!gen->status) continue;
      ierr = PSGENGetDYNGen(gen,&dyngen);CHKERRQ(ierr);
      dynlocglob = buslocglob + dyngen->startloc;
      PetscScalar PG=gen->pg,QG=gen->qg;
      ierr = DYNGenModelSetInitialConditionsP(dyngen,PG,QG,VA,VM,PGloc,QGloc,VAloc,VMloc,dyn->ICp,dynlocglob);CHKERRQ(ierr);

      dyngen->freq_max = dyn->freq_max;
      dyngen->freq_min = dyn->freq_min;

      /* Exciter Model */
      if(dyngen->dynexc) {
        ierr = PSGENGetDYNExc(gen,&dynexc);CHKERRQ(ierr);
        dynlocglob = buslocglob + dynexc->startloc;
        ierr = DYNExcModelSetInitialConditionsP(dynexc,PG,QG,VA,VM,PGloc,QGloc,VAloc,VMloc,dyn->ICp,dynlocglob);CHKERRQ(ierr);
      }

      /* Turbine governor Model */
      if(dyngen->dynturbgov) {
        ierr = PSGENGetDYNTurbgov(gen,&dynturbgov);CHKERRQ(ierr);
        dynlocglob = buslocglob + dynturbgov->startloc;
        ierr = DYNTurbgovModelSetInitialConditionsP(dynturbgov,PG,QG,VA,VM,PGloc,QGloc,VAloc,VMloc,dyn->ICp,dynlocglob);CHKERRQ(ierr);
      }

      /* Stabilizer Model */
      if(dyngen->dynexc && dynexc->dynstab) {
        ierr = PSGENGetDYNStab(gen,&dynstab);CHKERRQ(ierr);
        dynlocglob = buslocglob + dynstab->startloc;
        //ierr = DYNStabModelSetInitialConditions(dynstab,VA,VM,dy0dparr[VAloc]+loc,dy0dparr[VMloc]+loc);CHKERRQ(ierr);
      }
    }
    paramloc += 2*bus->ngen;
  }
  ierr = MatAssemblyBegin(dyn->ICp,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(dyn->ICp,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

#ifdef DEBUGFD
  for (i=0; i<dyn->Nparams; i++) {
    Vec col;
    ierr = MatCreateVecs(dyn->ICp,NULL,&col);CHKERRQ(ierr);
    ierr = MatGetColumnVector(dyn->ICp,col,i);CHKERRQ(ierr);
    ierr = VecView(col,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = VecDestroy(&col);CHKERRQ(ierr);
  }
#endif

#ifdef FWDSA
  if(dyn->sp) {
    for(i=0;i<dyn->Nparams;i++){
      ierr = VecCopy(dy0dp[i],dyn->sp[i]);CHKERRQ(ierr);
    }
  }else { /* mu = \sum_i <dy0dp[i],lambda[i]>*/
#endif
    for(i=0; i<dyn->ncostfcns; i++) {
      ierr = MatMultTransposeAdd(dyn->ICp,dyn->lambda[i],dyn->mu[i],dyn->mu[i]);CHKERRQ(ierr);
    }
#ifdef FWDSA
  }
#endif
  PetscFunctionReturn(0);
}

/*
  DYNComputeSensiP - Computes the sensitivity of the dynamic constraints w.r.t. the parameters

  Input Parameters:
. dyn - the DYN object
*/
PetscErrorCode DYNComputeSensiP(DYN dyn)
{
  PetscErrorCode ierr;
  PS             ps=dyn->ps;
  PetscInt       nc=dyn->ncostfcns;
  PetscInt       i,j,nbus=ps->nbus;
  PetscScalar    eps=1e-7;
  PetscBool      ghostbus;
  Vec            X;
  Vec            *dy0dp=dyn->dy0dp;
  PSBUS          bus;
  PetscInt       paramloc,ctr,k;
  PSGEN          gen;
  PetscInt       jj;

  PetscFunctionBegin;
  ierr = VecDuplicate(dyn->X,&X);CHKERRQ(ierr);
  ierr = DYNSetInitialConditions(dyn,X);CHKERRQ(ierr);

  ierr = VecGetOwnershipRange(dyn->Xparam,&paramloc,NULL);CHKERRQ(ierr);

  for(i=0; i < nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSIsGhosted(bus,&ghostbus);CHKERRQ(ierr);
    if(ghostbus) continue;

    /* Perturb Va and compute sensitivity */
    bus->va += eps*180.0/PETSC_PI;
    ierr = DYNSetInitialConditions(dyn,dy0dp[paramloc]);CHKERRQ(ierr);
    bus->va -= eps*180.0/PETSC_PI;

    /* Finite differencing */
    ierr = VecAXPY(dy0dp[paramloc],-1,X);CHKERRQ(ierr);
    ierr = VecScale(dy0dp[paramloc],1./eps);CHKERRQ(ierr);

    /* Perturb Vm and compute sensitivity */
    if(bus->ngen) {
      for(jj=0; jj < bus->ngen; jj++) {
        ierr = PSBUSGetGen(bus,jj,&gen);CHKERRQ(ierr);
        if(gen->status) gen->vs += eps;
      }
    }
    bus->vm += eps;

    ierr = DYNSetInitialConditions(dyn,dy0dp[paramloc+1]);CHKERRQ(ierr);
    if(bus->ngen) {
      ierr = PSBUSGetGen(bus,0,&gen);CHKERRQ(ierr);
      for(jj=0; jj < bus->ngen; jj++) {
        ierr = PSBUSGetGen(bus,jj,&gen);CHKERRQ(ierr);
        if(gen->status) gen->vs -= eps;
      }
    }
    bus->vm -= eps;

    /* Finite differencing */
    ierr = VecAXPY(dy0dp[paramloc+1],-1,X);CHKERRQ(ierr);
    ierr = VecScale(dy0dp[paramloc+1],1./eps);CHKERRQ(ierr);

    paramloc += 2;
    ctr = 0;
    for(k=0; k < bus->ngen; k++) {
      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      if(!gen->status) {
        ierr = VecSet(dy0dp[paramloc],0.0);CHKERRQ(ierr);
        ierr = VecSet(dy0dp[paramloc+1],0.0);CHKERRQ(ierr);
        ctr += 2;
        continue;
      }
      /* Perturb PG and compute sensitivity */
      gen->pg += eps;
      ierr = DYNSetInitialConditions(dyn,dy0dp[paramloc+ctr]);CHKERRQ(ierr);
      gen->pg -= eps;

      /* Finite differencing */
      ierr = VecAXPY(dy0dp[paramloc+ctr],-1,X);CHKERRQ(ierr);
      ierr = VecScale(dy0dp[paramloc+ctr],1./eps);CHKERRQ(ierr);

      /* Perturb QG and compute sensitivity */
      gen->qg += eps;
      ierr = DYNSetInitialConditions(dyn,dy0dp[paramloc+ctr+1]);CHKERRQ(ierr);
      gen->qg -= eps;

      /* Finite differencing */
      ierr = VecAXPY(dy0dp[paramloc+ctr+1],-1,X);CHKERRQ(ierr);
      ierr = VecScale(dy0dp[paramloc+ctr+1],1./eps);CHKERRQ(ierr);

      ctr += 2;
    }
    paramloc += ctr;
  }

  ierr = VecDestroy(&X);CHKERRQ(ierr);

#ifdef DEBUGFD
  for (i=0; i<dyn->Nparams; i++) {
    ierr = VecView(dy0dp[i],PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }
#endif

#ifdef FWDSA
  if(dyn->sp) {
    for(j=0;j<dyn->nparams;j++){
      ierr = VecCopy(dy0dp[j],dyn->sp[j]);CHKERRQ(ierr);
    }
  }else { /* mu = \sum_i <dy0dp[i],lambda[i]>*/
#endif
    Vec         *dcdy=dyn->lambda;
    Vec         *dcdp=dyn->mu;
    PetscScalar s,*dcp;
    for(i=0; i < nc; i++) {
      ierr = VecGetArray(dcdp[i],&dcp);CHKERRQ(ierr);
      for(j=0; j < dyn->nparams; j++) {
        ierr = VecDot(dcdy[i],dy0dp[j],&s);CHKERRQ(ierr);
        dcp[j] += s;
      }
      ierr = VecRestoreArray(dcdp[i],&dcp);CHKERRQ(ierr);
    }
#ifdef FWDSA
  }
#endif
  PetscFunctionReturn(0);
}

/*
  DYNComputeAdjointJacobianP - Computes the Jacobian of the DYN equations w.r.t. parameters for the adjoint

  Input Parameters:
+ ts - the TS solver
. t  - the current time
. x  - the solution vector
- ctx - application context (dyn)

  Output Parameters
. jacP - Jacobian of the DAE equations w.r.t. parameters
*/
PetscErrorCode DYNComputeAdjointJacobianP(TS ts,PetscReal t,Vec X,Mat jacP,void *ctx)
{
  DYN               dyn=(DYN)ctx;
  const PetscScalar *x,*xdyn;
  Vec               localX;
  PetscInt          i;
  PS                dynps=dyn->ps;
  PetscBool         ghostbus;
  PetscInt          dynlocglob,startloc,paramloc;
  PetscErrorCode    ierr;

  PetscFunctionBegin;

  ierr = MatZeroEntries(jacP);CHKERRQ(ierr);

  /* Get DYN local vector */
  ierr = DMGetLocalVector(dynps->networkdm,&localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(dynps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(dynps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);

  ierr = VecGetOwnershipRange(dyn->Xparam,&paramloc,NULL);CHKERRQ(ierr);

  ierr = VecGetArrayRead(localX,&x);CHKERRQ(ierr);

  for (i=0; i<dynps->nbus; i++) {
    PetscScalar VD,VQ;
    PetscInt    dynloc;
    PSBUS       dynbus;
    PetscInt    k,PGloc,QGloc,VAloc,VMloc,ctr=0;

    dynbus = &dynps->bus[i];
    ierr = PSBUSGetVariableGlobalLocation(dynbus,&dynloc);CHKERRQ(ierr);
    ierr = PSBUSIsGhosted(dynbus,&ghostbus);CHKERRQ(ierr);
    if (ghostbus) continue;

    VAloc = paramloc; VMloc = paramloc+1; paramloc += 2;
    VD = x[dynloc]; VQ = x[dynloc+1];

    if(dynbus->ide == ISOLATED_BUS) {
      PetscInt    row[2],col[2];
      PetscScalar val[4];

      row[0] = dynloc; row[1] = dynloc+1;
      col[0] = VAloc; col[1] = VMloc;
      val[0] = VD + dynbus->vm*PetscSinScalar(dynbus->va*PETSC_PI/180.0); // _VA
      val[1] = VD - PetscCosScalar(dynbus->va*PETSC_PI/180.0); // _VM
      val[2] = VQ - dynbus->vm*PetscCosScalar(dynbus->va*PETSC_PI/180.0); // _VA
      val[3] = VQ - PetscSinScalar(dynbus->va*PETSC_PI/180.0); // _VM
      ierr = SetMatrixValues(jacP,2,row,2,col,val);CHKERRQ(ierr);
      continue;
    }

    /* Shunt injection takes no effect to jacp*/

    /* Load injection */
    startloc   = dynloc;
    xdyn       = x + startloc;
    dynlocglob = startloc;
    for(k=0; k < dynbus->nload; k++) {
      PSLOAD       load;
      DYNLoadModel dynload;
      ierr = PSBUSGetLoad(dynbus,k,&load);CHKERRQ(ierr);
      if(!load->status) continue;

      ierr = PSLOADGetDYNLoad(load,&dynload);CHKERRQ(ierr);
      ierr = DYNLoadModelDAERHSJacobianP(dynload,t,xdyn,jacP,dynlocglob,VAloc,VMloc);CHKERRQ(ierr);
    }

    /* Generator frequency deviation */
    for(k=0; k < dynbus->ngen; k++) {
      PSGEN gen;
      DYNGenModel dyngen;
      DYNExcModel dynexc;
      DYNTurbgovModel dynturbgov;
      DYNStabModel dynstab;
      const PetscScalar *xdyn;
      PetscInt    startloc;
      PetscInt    dynlocglob;

      ierr = PSBUSGetGen(dynbus,k,&gen);CHKERRQ(ierr);
      ierr = PSGENGetDYNGen(gen,&dyngen);CHKERRQ(ierr);

      startloc   = dynloc;
      xdyn       = x+startloc;
      PGloc      = paramloc + ctr;
      QGloc      = paramloc + ctr + 1;
      ctr       += 2; // keep PG and QG even when the generator is off
      dynlocglob = startloc + dyngen->startloc;
      if (!gen->status) continue;
      /* Set partial derivatives of machine DAE equations w.r.t parameters */
      ierr = DYNGenModelDAERHSJacobianP(dyngen,t,xdyn,jacP,dynlocglob,VAloc,VMloc,PGloc,QGloc);CHKERRQ(ierr);
      if (dyngen->dynexc) {
        ierr = PSGENGetDYNExc(gen,&dynexc);CHKERRQ(ierr);
        dynlocglob = startloc + dynexc->startloc;
        ierr = DYNExcModelDAERHSJacobianP(dynexc,t,xdyn,jacP,dynlocglob,VAloc,VMloc);CHKERRQ(ierr);
      }
      if (dyngen->dynturbgov) {
        ierr = PSGENGetDYNTurbgov(gen,&dynturbgov);CHKERRQ(ierr);
        dynlocglob = startloc + dynturbgov->startloc;
        ierr = DYNTurbgovModelDAERHSJacobianP(dynturbgov,t,xdyn,jacP,dynlocglob,PGloc);CHKERRQ(ierr);
      }
      if (dyngen->dynexc && dynexc->dynstab) {
        ierr = PSGENGetDYNStab(gen,&dynstab);CHKERRQ(ierr);
        dynlocglob = startloc + dynstab->startloc;
        ierr = DYNStabModelDAERHSJacobianP(dynstab,t,xdyn,jacP,dynlocglob,PGloc);CHKERRQ(ierr);
      }
    }
    paramloc += ctr;
  }

  ierr = VecRestoreArrayRead(localX,&x);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(dynps->networkdm,&localX);CHKERRQ(ierr);

  ierr = MatAssemblyBegin(jacP,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(jacP,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode DYNSetUpGenFreqIS(DYN dyn,IS *isgenfreq)
{
  IS             is;
  PS             ps=dyn->ps;
  PSBUS          bus;
  PetscInt       i,k,ctr=0,locglob,*idx;
  PetscBool      ghostbus;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscCalloc1(dyn->ncostfcns,&idx);CHKERRQ(ierr);
  for(i=0; i<ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSIsGhosted(bus,&ghostbus);CHKERRQ(ierr);
    if(!ghostbus) {
      ierr = PSBUSGetVariableGlobalLocation(bus,&locglob);CHKERRQ(ierr);
      for(k=0; k<bus->ngen; k++) {
        PSGEN       gen;
        DYNGenModel dyngen;
        ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
        if(!gen->status) continue;
        ierr = PSGENGetDYNGen(gen,&dyngen);CHKERRQ(ierr);
        idx[ctr++] = locglob + dyngen->startloc + 5;
      }
    }
  }
  ierr = ISCreateGeneral(dyn->comm->type,dyn->ncostfcns,idx,PETSC_COPY_VALUES,&is);CHKERRQ(ierr);
  ierr = PetscFree(idx);CHKERRQ(ierr);
  *isgenfreq = is;
  PetscFunctionReturn(0);
}

static PetscErrorCode DYNCostIntegrand(TS ts,PetscReal t,Vec X,Vec R,void *ctx)
{
  PetscErrorCode    ierr;
  PetscScalar       *r;
  const PetscScalar *x;
  DYN               dyn=(DYN)ctx;
  PetscInt          ctr=0;

  PetscFunctionBegin;

  ierr = VecGetArray(R,&r);CHKERRQ(ierr);
  if(dyn->sum_freq_costfun) r[0] = 0;

  Vec localX;
  PetscInt i;
  PS       ps=dyn->ps;
  PetscBool ghostbus;

  ierr = DMGetLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = VecGetArrayRead(localX,&x);CHKERRQ(ierr);

  for(i=0; i < ps->nbus; i++) {
    PetscInt    loc;
    PSBUS       bus;
    PetscInt    k;

    bus = &ps->bus[i];
    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
    ierr = PSBUSIsGhosted(bus,&ghostbus);CHKERRQ(ierr);
    if (ghostbus) continue;

    if(bus->ide == ISOLATED_BUS) continue;

    /* Generator frequency deviation */
    for(k=0; k < bus->ngen; k++) {
      PSGEN gen;
      DYNGenModel dyngen;
      const PetscScalar *xdyn;
      PetscInt    startloc;
      PetscScalar frequency;
      PetscScalar freq_min = dyn->freq_min;
      PetscScalar freq_max = dyn->freq_max;
      PetscScalar freq_viol;

      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      ierr = PSGENGetDYNGen(gen,&dyngen);CHKERRQ(ierr);

      startloc = loc;
      xdyn = x+startloc;

      if(!gen->status) {
        if(!dyn->sum_freq_costfun) r[ctr++] = 0.0;
        continue;
      }
      /* Get Machine frequency */
      ierr = DYNGenModelGetFrequency(dyngen,t,xdyn,&frequency);
      freq_viol = dyn->scal*PetscPowScalarInt(PetscMax(0.0,PetscMax(frequency-freq_max,freq_min-frequency)),dyn->exp);CHKERRQ(ierr);
      if(dyn->sum_freq_costfun) r[0] += freq_viol;
      else r[ctr++] = freq_viol;
    }
  }

  ierr = VecRestoreArrayRead(localX,&x);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);
  ierr = VecRestoreArray(R,&r);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode DYNDRDYFunction(TS ts,PetscReal t,Vec X,Vec *drdy,void *ctx)
{
  PetscErrorCode    ierr;
  PetscScalar       *ry;
  const PetscScalar *x;
  DYN               dyn=(DYN)ctx;
  PetscInt          veci=0;

  PetscFunctionBegin;

  Vec localX;
  PetscInt  i;
  PS        ps=dyn->ps;
  PetscBool ghostbus;

  ierr = DMGetLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = VecGetArrayRead(localX,&x);CHKERRQ(ierr);

  for(i=0; i < dyn->ncostfcns; i++) {
    ierr = VecSet(drdy[i],0.0);CHKERRQ(ierr);
  }

  if(dyn->sum_freq_costfun) {
    ierr = VecGetArray(drdy[0],&ry);CHKERRQ(ierr);
  }

  for(i=0; i < ps->nbus; i++) {
    PetscInt    loc;
    PSBUS       bus;
    PetscInt    k;

    bus = &ps->bus[i];
    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
    ierr = PSBUSIsGhosted(bus,&ghostbus);CHKERRQ(ierr);
    if (ghostbus) continue;

    if(bus->ide == ISOLATED_BUS) continue;

    /* Generator frequency deviation */
    for(k=0; k < bus->ngen; k++) {
      PSGEN gen;
      DYNGenModel dyngen;
      const PetscScalar *xdyn;
      PetscInt    startloc;
      PetscScalar frequency;     /* Machine frequency */
      PetscScalar dfreq_dstate; /* Partial derivative of machine frequency w.r.t. to its state (dw, dn,w, others,delta) that governs its frequency */
      PetscInt    stateloc;     /* Location of the start relative to its bus */
      PetscScalar freq_min = dyn->freq_min;
      PetscScalar freq_max = dyn->freq_max;

      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      ierr = PSGENGetDYNGen(gen,&dyngen);CHKERRQ(ierr);

      startloc = loc;
      xdyn = x+startloc;

      if(gen->status) {
        /* Get Machine frequency and its partial derivative w.r.t machine state*/
        ierr = DYNGenModelGetFrequency(dyngen,t,xdyn,&frequency);CHKERRQ(ierr);
        ierr = DYNGenModelGetdFreqdState(dyngen,t,xdyn,&dfreq_dstate,&stateloc);

        if(!dyn->sum_freq_costfun) {
          ierr = VecGetArray(drdy[veci],&ry);CHKERRQ(ierr);
        }

        if(frequency > freq_max) ry[startloc+stateloc] += dyn->exp*dyn->scal*dfreq_dstate*PetscPowScalarInt(frequency-freq_max,dyn->exp-1);
        else if(frequency < freq_min) ry[startloc+stateloc] += -dyn->exp*dyn->scal*dfreq_dstate*PetscPowScalarInt(freq_min-frequency,dyn->exp-1);
        if(!dyn->sum_freq_costfun) {
          ierr = VecRestoreArray(drdy[veci++],&ry);CHKERRQ(ierr);
        }
      } else veci++;
    }
  }

  if(dyn->sum_freq_costfun) {
    ierr = VecRestoreArray(drdy[0],&ry);CHKERRQ(ierr);
  }

  ierr = VecRestoreArrayRead(localX,&x);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode DYNDRDPFunction(TS ts,PetscReal t,Vec U,Vec *drdp,void *ctx)
{
  PetscErrorCode    ierr;
  DYN               dyn=(DYN)ctx;
  PetscInt          i;

  PetscFunctionBegin;
  for(i=0; i<dyn->ncostfcns; i++) {
    ierr = VecSet(drdp[i],0.0);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*
  DYNSetComputePrePflow - Initial power flow solution for DYN

  Inputs:
+ dyn - the DYN object
- flag - compute initial power flow (0 = No, 1 = Yes)

  Notes:
   Calling this routine with PETSC_TRUE flag computes the steady-state power flow solution for the DYN object.
*/
PetscErrorCode DYNSetComputePrePflow(DYN dyn,PetscBool flag)
{
  PetscFunctionBegin;
  dyn->prepflow = flag;
  PetscFunctionReturn(0);
}

/*
  DYNCreate - Creates a dynamic simulation application object

  Input Parameters
. mpicomm - The MPI communicator

  Output Parameters
. dynout - The dynamic simulation application object
*/
PetscErrorCode DYNCreate(MPI_Comm mpicomm, DYN *dynout)
{
  PetscErrorCode ierr;
  DYN            dyn;

  PetscFunctionBegin;
  ierr = PetscCalloc1(1,&dyn);CHKERRQ(ierr);

  ierr = COMMCreate(mpicomm,&dyn->comm);CHKERRQ(ierr);

  ierr = PSCreate(mpicomm,&dyn->ps);CHKERRQ(ierr);

  /* Set the application with the PS object */
  ierr = PSSetApplication(dyn->ps,APP_DYNSIM);CHKERRQ(ierr);

  /* Create the time-stepping solver object */
  ierr = TSCreate(dyn->comm->type,&dyn->ts);CHKERRQ(ierr);

  /* Register generator models */
  ierr = DYNGenModelRegisterAll();CHKERRQ(ierr);
  /* Register exciter models */
  ierr = DYNExcModelRegisterAll();CHKERRQ(ierr);
  /* Register turbine governor models */
  ierr = DYNTurbgovModelRegisterAll();CHKERRQ(ierr);
  /* Register stabilizer models */
  ierr = DYNStabModelRegisterAll();CHKERRQ(ierr);
  /* Register load models */
  ierr = DYNLoadModelRegisterAll();CHKERRQ(ierr);

  /* Register events */
  ierr = DYNEventRegisterAll();CHKERRQ(ierr);

  dyn->ndiff  = 0;
  dyn->isdiff = 0;
  dyn->solve_alg_only = PETSC_FALSE;
  dyn->t0 = 0.0;
  dyn->tmax = 2.0;
  dyn->dt0 = 0.01;
  dyn->maxsteps = 10000;
  dyn->prepflow = PETSC_FALSE;
  dyn->use_semiexplicit = PETSC_FALSE;
  dyn->eval_integral = PETSC_FALSE;
  dyn->npoststepfcns = 0;
  dyn->monitor = PETSC_FALSE;
  dyn->useforward = PETSC_FALSE;
  dyn->useadjoint = PETSC_FALSE;
  dyn->nparams = 0;
  dyn->ncostfcns = 0;

  /* For parameter sensitivity calculation */
  dyn->scal = 1.0;
  dyn->exp  = 2;
  dyn->eta  = 1e-4;
  dyn->freq_max = 61.8; /* Instantaneous tripping (Eastern Interconnection standard) */
  dyn->freq_min = 57.8; /* Instantaneous tripping (Eastern Interconnection standard) */
  dyn->sum_freq_costfun = PETSC_FALSE;

  dyn->visualize = PETSC_FALSE;

  dyn->setupcalled = PETSC_FALSE;
  /*Initiate Event Log Count (SAM)*/
  dyn->eventlogcount = 0;

  dyn->switching_solution = PETSC_FALSE;
  *dynout = dyn;
  PetscFunctionReturn(0);
}

/*
  DYNDestroy - Destroys the DYN object

  Input Parameter
. dyn - The DYN object to destroy
*/
PetscErrorCode DYNDestroy(DYN *dyn)
{
  PetscErrorCode ierr;
  PetscInt       i;

  PetscFunctionBegin;

  if((*dyn)->prepflow) {
    /* Destroy objects created for initial power flow */
    ierr = PFLOWDestroy(&(*dyn)->initpflow);CHKERRQ(ierr);
    
    ierr = DMDestroy(&(*dyn)->initpflowdm);CHKERRQ(ierr);

    PetscObjectDereference((PetscObject)((*dyn)->initpflowpsection));
    PetscObjectDereference((PetscObject)((*dyn)->initpflowpglobsection));
    ierr = PetscSectionDestroy(&(*dyn)->initpflowpsection);CHKERRQ(ierr);
    ierr = PetscSectionDestroy(&(*dyn)->initpflowpglobsection);CHKERRQ(ierr);
    
  }

  ierr = COMMDestroy(&(*dyn)->comm);CHKERRQ(ierr);

  ierr = VecDestroy(&(*dyn)->X);CHKERRQ(ierr);
  ierr = MatDestroy(&(*dyn)->Jac);CHKERRQ(ierr);

  ierr = TSDestroy(&(*dyn)->ts);CHKERRQ(ierr);

  ierr = PSDestroy(&(*dyn)->ps);CHKERRQ(ierr);

  /* Destroy work vectors for jump conditoin */
  if ((*dyn)->Nevents && ((*dyn)->useforward || (*dyn)->useadjoint)) {
    ierr = VecDestroy(&(*dyn)->dGammadX);CHKERRQ(ierr);
    ierr = VecDestroy(&(*dyn)->dGammadP);CHKERRQ(ierr);
  }
  /* Destroy events */
  for(i=0; i < (*dyn)->Nevents;i++) {
    ierr = DYNEventDestroy(&(*dyn)->events[i]);CHKERRQ(ierr);
  }
  ierr = PetscFree((*dyn)->eventdirection);CHKERRQ(ierr);
  ierr = PetscFree((*dyn)->eventterminate);CHKERRQ(ierr);
  ierr = PetscFree((*dyn)->feventtmp);CHKERRQ(ierr);

  ierr = ISDestroy(&(*dyn)->isdiff);CHKERRQ(ierr);

  ierr = VecDestroy(&(*dyn)->vatol);CHKERRQ(ierr);
  ierr = VecDestroy(&(*dyn)->vrtol);CHKERRQ(ierr);

  ierr = PetscFree(*dyn);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  DYNAdjointDestroy - Destroys the DYN object

  Input Parameter
. dyn - The DYN object to destroy
*/
PetscErrorCode DYNAdjointDestroy(DYN *dyn)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = VecDestroyVecs((*dyn)->ncostfcns,&(*dyn)->lambda);CHKERRQ(ierr);
  ierr = VecDestroyVecs((*dyn)->ncostfcns,&(*dyn)->mu);CHKERRQ(ierr);
  ierr = VecDestroyVecs((*dyn)->nparams,&(*dyn)->dy0dp);CHKERRQ(ierr);
  ierr = MatDestroy(&(*dyn)->Jacp);CHKERRQ(ierr);
  ierr = MatDestroy(&(*dyn)->ICp);CHKERRQ(ierr);
  ierr = VecDestroy(&(*dyn)->Xparam);CHKERRQ(ierr);
  ierr = VecDestroy(&(*dyn)->Xtmp);CHKERRQ(ierr);
  ierr = VecDestroy(&(*dyn)->Xtmp2);CHKERRQ(ierr);
  ierr = VecDestroy(&(*dyn)->Xtmp3);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* 
   DYNSetUpStaticLoadModels - Creates static load models at each load bus if there
   are no other load models given

   Input Parameter
. DYN - the DYN object
*/
PetscErrorCode DYNSetUpStaticLoadModels(DYN dyn)
{
  PetscErrorCode ierr;
  PS             ps=dyn->ps;
  PetscInt       i;
  PSLOAD         load;
  PSBUS          Bus;
  DYNZIP         zip;

  PetscFunctionBegin;

  for(i=0; i < ps->nload; i++) {
    load = &ps->load[i];
    if(load->status) {
      Bus  = &ps->bus[load->internal_i];
      if(load->dynloadsetup == 0) {
	ierr = DYNLoadModelSetType(&load->dynload,"ZIP");CHKERRQ(ierr);
	ierr = DYNLoadModelGetModelData(&load->dynload,(void*)&zip);CHKERRQ(ierr);
	zip->bus_i = load->bus_i;
	ierr = PetscStrcpy(zip->id,load->id);CHKERRQ(ierr);

	load->dyniszip = PETSC_TRUE;
	load->dynloadsetup = 1;
      }
    }
  }

  PetscFunctionReturn(0);
}

/* 
   DYNSetUpConstantVoltageGenModels - Creates constant voltage source generator models for generators with no dynamic data

   Input Parameter
. DYN - the DYN object
*/
PetscErrorCode DYNSetUpConstantVoltageGenModels(DYN dyn)
{
  PetscErrorCode ierr;
  PS             ps=dyn->ps;
  PetscInt       i;
  PSGEN          gen;
  PSBUS          Bus;
  DYNCv          cv;

  PetscFunctionBegin;

  for(i=0; i < ps->ngen; i++) {
    gen = &ps->gen[i];
    if(gen->status) {
      Bus  = &ps->bus[gen->internal_i];
      if(gen->dyngensetup == 0) {
	ierr = DYNGenModelSetType(&gen->dyngen,"CV");CHKERRQ(ierr);
	ierr = DYNGenModelGetModelData(&gen->dyngen,(void*)&cv);CHKERRQ(ierr);
	cv->bus_i = gen->bus_i;
	ierr = PetscStrcpy(cv->id,gen->id);CHKERRQ(ierr);

	gen->dyniscv = PETSC_TRUE;
	gen->dyngensetup = 1;
      }
    }
  }

  PetscFunctionReturn(0);
}

/*
  DYNReadMatPowerData - Reads the network data given in MATPOWER data format

  Input Parameter
+  DYN - The DYN object
-  netfile - The name of the network file

*/
PetscErrorCode DYNReadMatPowerData(DYN dyn,const char netfile[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  /* Read MatPower data file and populate the PS data structure */
  ierr = PSReadMatPowerData(dyn->ps,netfile);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNReadPSSERawData - Reads the network data given in PSSE raw data format

  Input Parameter
+  DYN - The DYN object
-  netfile - The name of the network file

*/
PetscErrorCode DYNReadPSSERawData(DYN dyn,const char netfile[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  /* Read PSSE raw data file and populate the PS data structure */
  ierr = PSReadPSSERawData(dyn->ps,netfile);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNReadDyrData - Reads the data file with dynamic models

  Input Parameter
+  DYN - The DYN object
-  dyrfile - The name of the dyr file

*/
PetscErrorCode DYNReadDyrData(DYN dyn,const char dyrfile[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  /* Read dyr data file and populate the PS data structure */
  ierr = PSReadDyrData(dyn->ps,dyrfile);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNSetDuration - Sets the duration of the dynamic simulation

  Input Parameters
+  DYN - the DYN object
.  max_steps  - the maximum number of steps that the time-stepping solver takes
-  tend - the end time
*/
PetscErrorCode DYNSetDuration(DYN dyn,PetscInt max_steps,PetscReal tend)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = TSSetMaxTime(dyn->ts,tend);CHKERRQ(ierr);
  dyn->tmax = tend;
  dyn->maxsteps = max_steps;
  PetscFunctionReturn(0);
}

/*
  DYNSetInitialTimeAndStep - Sets the start time and the time step of the dynamic simulation

  Input Parameters
+  DYN - the DYN object
.  start_time - start time
.  time_step  - the step size (in seconds)

   Notes:
   For variable time-stepping methods, this step is used as the initial time step.
*/
PetscErrorCode DYNSetStartTimeAndStep(DYN dyn,PetscReal start_time, PetscReal time_step)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = TSSetTime(dyn->ts,start_time);CHKERRQ(ierr);
  ierr = TSSetTimeStep(dyn->ts,time_step);CHKERRQ(ierr);
  dyn->t0 = start_time;
  dyn->dt0 = time_step;
  PetscFunctionReturn(0);
}

/*
  DYNCreateGlobalVector - Returns a global vector of the appropriate size
  and distribution conforming to the distribution of the PS object.

  Input Paramereters:
. DYN - the dynamics simulation application object

  Output Parameters:
. vec - the global vector

  Notes:
  DYNSetUp() must be called before calling this routine.
*/
PetscErrorCode DYNCreateGlobalVector(DYN dyn,Vec *vec)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if(!dyn->setupcalled) SETERRQ(dyn->comm->type,0,"DYNSetUp() must be called before calling DYNCreateGlobalVector");
  ierr = PSCreateGlobalVector(dyn->ps,vec);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNCreateMatrix - Returns a distributed matrix of appropriate size that can
   be used as the Jacobian

  Input Paramereters:
. DYN - the dynamics simulation application object

  Output Parameters:
. mat - the matrix

  Notes:
  DYNSetUp() must be called before calling this routine.
*/
PetscErrorCode DYNCreateMatrix(DYN dyn,Mat *mat)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if(!dyn->setupcalled) SETERRQ(dyn->comm->type,0,"DYNSetUp() must be called before calling DYNCreateMatrix");
  ierr = PSCreateMatrix(dyn->ps,mat);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNIJacobian - The Jacobian evaluation routine for the dynamics simulation

  Notes:
  See PETSc routine TSSetIJacobian()

*/
PetscErrorCode DYNIJacobian(TS ts,PetscReal t,Vec X,Vec Xdot,PetscReal a, Mat J, Mat Jpre, void *ctx)
{
  PetscErrorCode ierr;
  DYN            dyn=(DYN)ctx;
  PS             ps=dyn->ps;
  Vec            localX;
  const PetscScalar *xarr;
  PetscBool      ghostbus;
  PetscInt       i;

  PetscFunctionBegin;
  ierr = MatZeroEntries(J);CHKERRQ(ierr);

  ierr = DMGetLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);

  ierr = DMGlobalToLocalBegin(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);

  ierr = VecGetArrayRead(localX,&xarr);CHKERRQ(ierr);

  for(i=0; i < ps->nbus; i++) {
    PetscScalar VD,VQ;
    PetscInt    loc; /* Location in the local array for accessing the local vector values*/
    PetscInt    locglob; /* location in the global array to set the matrix entries */
    PSBUS       bus;
    PetscScalar bshunt,gshunt;
    PetscInt    j,k;
    PetscInt    row[2],col[2];
    PetscScalar val[4];
    PetscScalar *xdyn;
    PetscInt    startloc;
    PetscInt V_loc[2],I_loc[2],dynlocglob;

    bus = &ps->bus[i];

    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
    ierr = PSBUSGetVariableGlobalLocation(bus,&locglob);CHKERRQ(ierr);
    ierr = PSBUSIsGhosted(bus,&ghostbus);CHKERRQ(ierr);

    VD = xarr[loc]; VQ = xarr[loc+1];

    V_loc[0] = locglob; V_loc[1] = locglob+1; /* global location for V */
    I_loc[0] = locglob+1; I_loc[1] = locglob; /* global location for I */

    row[0] = locglob; row[1] = locglob+1;
    col[0] = locglob; col[1] = locglob+1;

    /* Isolated buses */
    if(bus->ide == ISOLATED_BUS) {
      val[0] = val[3] = 1.0;
      val[1] = val[2] = 0.0;
      ierr = MatSetValues(J,2,row,2,col,val,ADD_VALUES);CHKERRQ(ierr);
      continue;
    }

    if(!ghostbus) {
      /* Shunt injection */
      bshunt = bus->bl;
      gshunt = bus->gl;

      /* dIshuntQ_dVDQ */
      val[0] = -bshunt; val[1] = -gshunt;
      /* dIshuntD_dVDQ */
      val[2] = -gshunt; val[3] = bshunt;
      ierr = MatSetValues(J,2,row,2,col,val,ADD_VALUES);CHKERRQ(ierr);

      /* Load injection */
      for(k=0; k < bus->nload; k++) {
        PSLOAD load;
	DYNLoadModel dynload;

	ierr = PSBUSGetLoad(bus,k,&load);CHKERRQ(ierr);

	if(!load->status) continue;

	ierr = PSLOADGetDYNLoad(load,&dynload);CHKERRQ(ierr);
	
	startloc = loc;
	xdyn     = (PetscScalar*)(xarr+startloc);

	dynlocglob = locglob + dynload->startloc;  /* starting global location of xdyngen */
	ierr = DYNLoadModelDAERHSJacobian(dynload,J,t,VD,VQ,xdyn,dynlocglob,V_loc,I_loc);CHKERRQ(ierr);
	for(j=0; j < dynload->nvar; j++) {
	  row[0] = col[0] = dynlocglob+j;
	  val[0] = -a*dynload->eqtypes[j];
	  ierr = MatSetValues(J,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);
	}

	/*
        ierr = PSBUSGetLoad(bus,k,&load);CHKERRQ(ierr);
        if(!load->status) continue;
        PetscScalar Pd,Qd,yp,yq;
        Pd = load->pl;
        Qd = load->ql;

        yp = Pd/(bus->vm*bus->vm);
        yq = Qd/(bus->vm*bus->vm);

        val[0] = yq; val[1] = -yp;
        val[2] = -yp; val[3] = -yq;

        ierr = MatSetValues(J,2,row,2,col,val,ADD_VALUES);CHKERRQ(ierr);
	*/
      }

      /* Generator jacobian */
      for(k=0; k < bus->ngen; k++) {
	PSGEN gen;
	DYNGenModel dyngen;
	DYNExcModel dynexc;
	DYNTurbgovModel dynturbgov;
	DYNStabModel dynstab;

	ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);

	if(!gen->status) continue;

	ierr = PSGENGetDYNGen(gen,&dyngen);CHKERRQ(ierr);
	
	startloc = loc;
	xdyn     = (PetscScalar*)(xarr+startloc);

	dynlocglob = locglob + dyngen->startloc;  /* starting global location of xdyngen */
	ierr = DYNGenModelDAERHSJacobian(dyngen,J,t,VD,VQ,xdyn,dynlocglob,V_loc,I_loc);CHKERRQ(ierr);
	for(j=0; j < dyngen->nvar; j++) {
	  row[0] = col[0] = dynlocglob+j;
	  val[0] = -a*dyngen->eqtypes[j];
	  ierr = MatSetValues(J,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);
	}

	/* Exciter Jacobian */
	if(dyngen->dynexc) {
	  ierr = PSGENGetDYNExc(gen,&dynexc);CHKERRQ(ierr);

	  dynlocglob = locglob + dynexc->startloc;
	  ierr = DYNExcModelDAERHSJacobian(dynexc,J,t,VD,VQ,xdyn,dynlocglob,V_loc);CHKERRQ(ierr);
	  for(j=0; j < dynexc->nvar; j++) {
	    row[0] = col[0] = dynlocglob+j;
	    val[0] = -a*dynexc->eqtypes[j];
	    ierr = MatSetValues(J,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);
	  }
	} 

	/* Turbine governor Jacobian */
	if(dyngen->dynturbgov) {
	  ierr = PSGENGetDYNTurbgov(gen,&dynturbgov);CHKERRQ(ierr);

	  dynlocglob = locglob + dynturbgov->startloc;
	  ierr = DYNTurbgovModelDAERHSJacobian(dynturbgov,J,t,VD,VQ,xdyn,dynlocglob,V_loc);CHKERRQ(ierr);
	  for(j=0; j < dynturbgov->nvar; j++) {
	    row[0] = col[0] = dynlocglob+j;
	    val[0] = -a*dynturbgov->eqtypes[j];
	    ierr = MatSetValues(J,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);
	  }
	} 

	/* Stabilizer Jacobian */
	if(dyngen->dynexc && dynexc->dynstab) {
	  ierr = PSGENGetDYNStab(gen,&dynstab);CHKERRQ(ierr);

	  dynlocglob = locglob + dynstab->startloc;
	  ierr = DYNStabModelDAERHSJacobian(dynstab,J,t,VD,VQ,xdyn,dynlocglob,V_loc);CHKERRQ(ierr);
	  for(j=0; j < dynstab->nvar; j++) {
	    row[0] = col[0] = dynlocglob+j;
	    val[0] = -a*dynstab->eqtypes[j];
	    ierr = MatSetValues(J,1,row,1,col,val,ADD_VALUES);CHKERRQ(ierr);
	  }
	} 
      }
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
      PetscInt locglobf,locglobt;

      /* Get the connected buses to this line */
      ierr = PSLINEGetConnectedBuses(line,&connbuses);CHKERRQ(ierr);
      busf = connbuses[0];
      bust = connbuses[1];

      ierr = PSBUSGetVariableGlobalLocation(busf,&locglobf);CHKERRQ(ierr);
      ierr = PSBUSGetVariableGlobalLocation(bust,&locglobt);CHKERRQ(ierr);

      if(bus == busf) {
        row[0] = locglobf; row[1] = locglobf+1;
        col[0] = locglobf; col[1] = locglobf+1;
        val[0] = -Bff; val[1] = -Gff; val[2] = -Gff; val[3] = Bff;
        ierr = MatSetValues(J,2,row,2,col,val,ADD_VALUES);CHKERRQ(ierr);

        row[0] = locglobf; row[1] = locglobf+1;
        col[0] = locglobt; col[1] = locglobt+1;
        val[0] = -Bft; val[1] = -Gft; val[2] = -Gft; val[3] = Bft;
        ierr = MatSetValues(J,2,row,2,col,val,ADD_VALUES);CHKERRQ(ierr);
      } else {
        row[0] = locglobt; row[1] = locglobt+1;
        col[0] = locglobt; col[1] = locglobt+1;
        val[0] = -Btt; val[1] = -Gtt; val[2] = -Gtt; val[3] = Btt;
        ierr = MatSetValues(J,2,row,2,col,val,ADD_VALUES);CHKERRQ(ierr);

        row[0] = locglobt; row[1] = locglobt+1;
        col[0] = locglobf; col[1] = locglobf+1;
        val[0] = -Btf; val[1] = -Gtf; val[2] = -Gtf; val[3] = Btf;
        ierr = MatSetValues(J,2,row,2,col,val,ADD_VALUES);CHKERRQ(ierr);
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
  DYNRHSFunction - Computes the RHS of the DAE power grid equations (only used
                   with semi-explicit time-stepping solver, i.e., when the
                   option -dyn_use_semiexplicit is used

   Notes:
   See PETSc routine TSSetRHSFunction()
*/
PetscErrorCode DYNRHSFunction(TS ts, PetscReal t, Vec X, Vec F, void *ctx)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = DYNIFunction(ts,t,X,(Vec)NULL,F,ctx);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  DYNAlgebraicSolveOnly - Solves the algebraic part of the power grid DAE equations

   Notes:
   This function gets called only when PETSc's explicit time-stepping solvers are
   used. In order to use the explicit time-stepping solver use the option
   -dyn_use_semiexplicit.
*/
PetscErrorCode DYNAlgebraicSolveOnly(DYN dyn, PetscReal t, Vec X)
{
  PetscErrorCode ierr;
  SNES tssnes;
  void *tssnesctx;
  void *tssnesjacctx;
  PetscErrorCode (*tssnesfunc)(SNES,Vec,Vec,void*);
  PetscErrorCode (*tssnesjacfun)(SNES,Vec,Mat,Mat,void*);
  PetscInt tssneslag;
  
  PetscFunctionBegin;

  /* Set solve_alg_only flag and change the function and Jacobian evaluation routine pointers */
  dyn->solve_alg_only = PETSC_TRUE;
  /* Get Inner SNES */
  ierr = TSGetSNES(dyn->ts,&tssnes);CHKERRQ(ierr);
  /* Get the function and jacobian evaluation routines. SNES associates it with the DM attached to it */
  ierr = DMSNESGetFunction(dyn->ps->networkdm,&tssnesfunc,&tssnesctx);CHKERRQ(ierr);
  ierr = DMSNESGetJacobian(dyn->ps->networkdm,&tssnesjacfun,&tssnesjacctx);CHKERRQ(ierr);
  /* Set function and jacobian evaluation routines for the algebtraic solve */
  ierr = DMSNESSetFunction(dyn->ps->networkdm,DYNEventPostSolveAlgebraic,(void*)dyn);CHKERRQ(ierr);
  ierr = DMSNESSetJacobian(dyn->ps->networkdm,DYNEventPostSolveAlgebraicJacobian,(void*)dyn);CHKERRQ(ierr);
  /* Some methods, such as ROSW, compute Jacobian only once at beginning of time-stepping, this may result in an incorrect
     Jacobian for algbraic solve. Need to reset the lag here
  */
  ierr = SNESGetLagJacobian(tssnes,&tssneslag);CHKERRQ(ierr);
  ierr = SNESSetLagJacobian(tssnes,1);CHKERRQ(ierr);
  
  /* Solve algebraic equations */
  ierr = SNESSolve(tssnes,NULL,dyn->X);CHKERRQ(ierr);
  
  /* Reset Jacobian lag */
  ierr = SNESSetLagJacobian(tssnes,tssneslag);CHKERRQ(ierr);
  /* Revert to the original function and jacobian evaluation routines */
  ierr = DMSNESSetFunction(dyn->ps->networkdm,tssnesfunc,tssnesctx);CHKERRQ(ierr);
  ierr = DMSNESSetJacobian(dyn->ps->networkdm,tssnesjacfun,tssnesjacctx);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

/*
  DYNPostStage - Post stage function for explicit scheme, solves the algebraic part only.

   Notes:
   This function gets called only when PETSc's explicit time-stepping solvers are
   used. In order to use the explicit time-stepping solver use the option
   -dyn_use_semiexplicit.
*/
PetscErrorCode DYNPostStage(TS ts, PetscReal t, PetscInt i, Vec *X)
{
  PetscErrorCode ierr;
  DYN            dyn;
  
  PetscFunctionBegin;
  ierr = TSGetApplicationContext(ts,&dyn);CHKERRQ(ierr);
  ierr = DYNAlgebraicSolveOnly(dyn,t,X[i]);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  DYNPostEvaluate - Post evaluate function for explicit scheme, solves the algebraic part only.

   Notes:
   This function gets called only when PETSc's explicit time-stepping solvers are
   used. In order to use the explicit time-stepping solver use the option
   -dyn_use_semiexplicit.
*/
PetscErrorCode DYNPostEvaluate(TS ts)
{
  PetscErrorCode ierr;
  Vec            X;
  PetscReal      t,dt;
  DYN            dyn;
  
  PetscFunctionBegin;

  ierr = TSGetApplicationContext(ts,&dyn);CHKERRQ(ierr);
  ierr = TSGetSolution(ts,&X);CHKERRQ(ierr);
  ierr = TSGetTime(ts,&t);CHKERRQ(ierr);
  ierr = TSGetTimeStep(ts,&dt);CHKERRQ(ierr);
  ierr = DYNAlgebraicSolveOnly(dyn,t+dt,X);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  DYNIFunction - The IFunction for the time-stepping solver

  Notes:
  See PETSc routine TSSetIFunction()

  The equations are expressed in current balance form:
    I_gen - I_network (=Y*V) + I_shunt - I_load = 0

    The current balance equations for each bus are ordered as [I_bus(imag);I_bus(real)]
    Ordering them in such a manner allows having the susceptance B in the diagonal location
    of the Jacobian matrix. For transmission networks, B > G (in many cases G=0) and hence
    having this ordering is better.
*/
PetscErrorCode DYNIFunction(TS ts,PetscReal t,Vec X,Vec Xdot,Vec F,void *ctx)
{
  PetscErrorCode ierr;
  DYN            dyn=(DYN)ctx;
  PS             ps=dyn->ps;
  Vec            localX,localF,localXdot;
  PetscScalar    *farr;
  const PetscScalar *xarr,*xdotarr;
  PetscBool      ghostbus;
  PetscInt       i;
  PetscScalar *fdyn;
  const PetscScalar *xdyn,*xdotdyn;
  PetscInt    startloc;

  PetscFunctionBegin;
  ierr = VecSet(F,0.0);CHKERRQ(ierr);

  ierr = DMGetLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = VecGetArrayRead(localX,&xarr);CHKERRQ(ierr);

  if(!dyn->solve_alg_only && Xdot) {
    ierr = DMGetLocalVector(ps->networkdm,&localXdot);CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(ps->networkdm,Xdot,INSERT_VALUES,localXdot);CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(ps->networkdm,Xdot,INSERT_VALUES,localXdot);CHKERRQ(ierr);
    ierr = VecGetArrayRead(localXdot,&xdotarr);CHKERRQ(ierr);
  }

  ierr = DMGetLocalVector(ps->networkdm,&localF);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ps->networkdm,F,INSERT_VALUES,localF);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ps->networkdm,F,INSERT_VALUES,localF);CHKERRQ(ierr);
  ierr = VecGetArray(localF,&farr);CHKERRQ(ierr);

  for(i=0; i < ps->nbus; i++) {
    PetscScalar VD,VQ;
    PetscInt    loc;
    PSBUS       bus;
    PetscScalar bshunt,gshunt;
    PetscInt    j,k;

    bus = &ps->bus[i];

    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
    ierr = PSBUSIsGhosted(bus,&ghostbus);CHKERRQ(ierr);
#ifdef PINDEX
    printf("loc=%d %d\n",loc,bus->ide);
#endif
    VD = xarr[loc]; VQ = xarr[loc+1];

    /* For isolated, or out-of-service, buses the residual function serves to keep VD and VQ at the bus constant */
    if(bus->ide == ISOLATED_BUS) {
      farr[loc]   = VD - bus->vm*PetscCosScalar(bus->va*PETSC_PI/180.0);
      farr[loc+1] = VQ - bus->vm*PetscSinScalar(bus->va*PETSC_PI/180.0);
      continue;
    }

    if(!ghostbus) {
      /* Shunt injection */
      bshunt = bus->bl;
      gshunt = bus->gl;
      PetscScalar IshuntD, IshuntQ;

      IshuntD = gshunt*VD - bshunt*VQ;
      IshuntQ = bshunt*VD + gshunt*VQ;

      farr[loc]   -= IshuntQ;
      farr[loc+1] -= IshuntD;

      /* Load injection */
      /* Assume for now that all loads are modeled as constant impedances */
      for(k=0; k < bus->nload; k++) {
        PSLOAD load;
	DYNLoadModel dynload;
        PetscScalar IloadD=0.0,IloadQ=0.0;

        ierr = PSBUSGetLoad(bus,k,&load);CHKERRQ(ierr);
        if(!load->status) continue;

	ierr = PSLOADGetDYNLoad(load,&dynload);CHKERRQ(ierr);

	startloc = loc;
	xdyn = xarr + startloc;

	if(!dyn->solve_alg_only && Xdot) {
	  xdotdyn = xdotarr + startloc;
	}
	fdyn    = farr + startloc;

	/* Get the RHS of the load DAE function */
	ierr = DYNLoadModelDAERHSFunction(dynload,t,VD,VQ,(PetscScalar*)xdyn,fdyn,&IloadD,&IloadQ);CHKERRQ(ierr);
	if(!dyn->solve_alg_only && Xdot) {
	  for(j=0; j < dynload->nvar;j++) {
	    fdyn[dynload->startloc+j] -= dynload->eqtypes[j]*xdotdyn[dynload->startloc+j];
	  }
	}

        farr[loc]   -= IloadQ;
        farr[loc+1] -= IloadD;
      }
#ifdef PINDEX
      printf("busngen=%d\n",bus->ngen);
#endif
      /* Generator residual and injection */
      for(k=0; k < bus->ngen; k++) {
	PSGEN gen;
	DYNGenModel dyngen;
	DYNExcModel dynexc;
	DYNTurbgovModel dynturbgov;
	DYNStabModel    dynstab;
	PetscScalar IGD=0.0,IGQ=0.0;

	ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
#ifdef PINDEX
    printf("genstatus=%d\n",gen->status);
#endif
	if(!gen->status) continue;

	ierr = PSGENGetDYNGen(gen,&dyngen);CHKERRQ(ierr);

	startloc = loc;
	xdyn    = xarr+startloc;
	if(!dyn->solve_alg_only && Xdot) {
	  xdotdyn = xdotarr + startloc;
	}
	fdyn    = farr + startloc;

	/* Get the RHS of the machine DAE function */
	ierr = DYNGenModelDAERHSFunction(dyngen,t,VD,VQ,(PetscScalar*)xdyn,fdyn,&IGD,&IGQ);CHKERRQ(ierr);
	if(!dyn->solve_alg_only && Xdot) {
	  for(j=0; j < dyngen->nvar;j++) {
	    fdyn[dyngen->startloc+j] -= dyngen->eqtypes[j]*xdotdyn[dyngen->startloc+j];
	  }
	}
	/* Get the RHS of the exciter DAE function */
	if(dyngen->dynexc) {
	  ierr = PSGENGetDYNExc(gen,&dynexc);CHKERRQ(ierr);
	  ierr = DYNExcModelDAERHSFunction(dynexc,t,VD,VQ,(PetscScalar*)xdyn,fdyn);CHKERRQ(ierr);
	  if(!dyn->solve_alg_only && Xdot) {
	    for(j=0; j < dynexc->nvar; j++) {
	      fdyn[dynexc->startloc+j] -= dynexc->eqtypes[j]*xdotdyn[dynexc->startloc+j];
	    }
	  }
	}

	/* Get the RHS of the turbine governor DAE function */
	if(dyngen->dynturbgov) {
	  ierr = PSGENGetDYNTurbgov(gen,&dynturbgov);CHKERRQ(ierr);
	  ierr = DYNTurbgovModelDAERHSFunction(dynturbgov,t,VD,VQ,(PetscScalar*)xdyn,fdyn);CHKERRQ(ierr);
	  if(!dyn->solve_alg_only && Xdot) {
	    for(j=0; j < dynturbgov->nvar; j++) {
	      fdyn[dynturbgov->startloc+j] -= dynturbgov->eqtypes[j]*xdotdyn[dynturbgov->startloc+j];
	    }
	  }
	}

	/* Get the RHS of the stabilizer DAE function */
	if(dyngen->dynexc && dynexc->dynstab) {
	  ierr = PSGENGetDYNStab(gen,&dynstab);CHKERRQ(ierr);
	  ierr = DYNStabModelDAERHSFunction(dynstab,t,VD,VQ,(PetscScalar*)xdyn,fdyn);CHKERRQ(ierr);
	  if(!dyn->solve_alg_only && Xdot) {
	    for(j=0; j < dynstab->nvar; j++) {
	      fdyn[dynstab->startloc+j] -= dynstab->eqtypes[j]*xdotdyn[dynstab->startloc+j];
	    }
	  }
	}

	farr[loc]   += IGQ;
	farr[loc+1] += IGD;
      }
    }

    /* Current injection from lines */
    PetscInt nconnlines;
    const PSLINE   *connlines;
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
      PetscInt locf,loct;
      PetscScalar VDf,VQf,VDt,VQt;

      /* Get the connected buses to this line */
      ierr = PSLINEGetConnectedBuses(line,&connbuses);CHKERRQ(ierr);
      busf = connbuses[0];
      bust = connbuses[1];

      ierr = PSBUSGetVariableLocation(busf,&locf);CHKERRQ(ierr);
      ierr = PSBUSGetVariableLocation(bust,&loct);CHKERRQ(ierr);
      VDf = xarr[locf]; VQf = xarr[locf+1];
      VDt = xarr[loct]; VQt = xarr[loct+1];

      if(bus == busf) {
        farr[locf]   -= Bff*VDf + Gff*VQf + Bft*VDt + Gft*VQt;
        farr[locf+1] -= Gff*VDf - Bff*VQf + Gft*VDt - Bft*VQt;
      } else {
        farr[loct]   -= Btt*VDt + Gtt*VQt + Btf*VDf + Gtf*VQf;
        farr[loct+1] -= Gtt*VDt - Btt*VQt + Gtf*VDf - Btf*VQf;
      }
    }
  }

  ierr = VecRestoreArrayRead(localX,&xarr);CHKERRQ(ierr);
  ierr = VecRestoreArray(localF,&farr);CHKERRQ(ierr);

  ierr = DMLocalToGlobalBegin(ps->networkdm,localF,ADD_VALUES,F);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(ps->networkdm,localF,ADD_VALUES,F);CHKERRQ(ierr);

  ierr = DMRestoreLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ps->networkdm,&localF);CHKERRQ(ierr);

  if(!dyn->solve_alg_only && Xdot) {
    ierr = VecRestoreArrayRead(localXdot,&xdotarr);CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(ps->networkdm,&localXdot);CHKERRQ(ierr);
  }

  //  ierr = VecView(F,0);CHKERRQ(ierr);
  //  PetscReal fnorm;
  //  ierr = VecNorm(F,NORM_2,&fnorm);CHKERRQ(ierr);
  //  ierr = PetscPrintf(PETSC_COMM_WORLD,"Ifunction norm = %lf\n",fnorm);CHKERRQ(ierr);
  //  exit(1);

  PetscFunctionReturn(0);
}

/*
  DYNSetInitialConditions - Initializes the network and the dynamic state variables in the solution
                            vector X.

  Input Parameters:
. DYN - the dynamics simulation application object

  Output Parameters:
. X - the solution vector
*/
PetscErrorCode DYNSetInitialConditions(DYN dyn,Vec X)
{
  PetscErrorCode ierr;
  PS             ps=dyn->ps;
  Vec            localX;
  PetscScalar    *xarr;
  PetscInt       nbus=ps->nbus,i,loc,k,n=3;
  PetscBool      ghostbus;
  PSBUS          bus;
  PetscScalar    VD,VQ; /* Real and imaginary parts of bus voltages */
  PetscScalar    Vm,Va; /* Voltage magnitude and angle */
  PetscReal      zip_p[3]={0,0,1.0},zip_q[3]={0,0,1.0},Vm_thresh=0.8;
  PetscScalar    *xbus;

  PetscFunctionBegin;

  /* Get the zip load model composition */
  ierr = PetscOptionsGetRealArray(NULL,NULL,"-dyn_zip_p",zip_p,&n,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetRealArray(NULL,NULL,"-dyn_zip_q",zip_q,&n,NULL);CHKERRQ(ierr);

  ierr = PetscOptionsGetReal(NULL,NULL,"-dyn_zip_Vm_thresh",&Vm_thresh,NULL);CHKERRQ(ierr);

  if((zip_p[0] + zip_p[1] + zip_p[2] - 1.0) > PETSC_SMALL) SETERRQ3(PETSC_COMM_SELF,0,"Inconsistent zip load model composition: cp = %f + ip = %f + yp = %f should equal 1.0",zip_p[0],zip_p[1],zip_p[2]);
  if((zip_q[0] + zip_q[1] + zip_q[2] - 1.0) > PETSC_SMALL) SETERRQ3(PETSC_COMM_SELF,0,"Inconsistent zip load model composition: cq = %f + iq = %f + yq = %f should equal 1.0",zip_q[0],zip_q[1],zip_q[2]);
  

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
	if(bus->ide != PQ_BUS) Vm = gen->vs;
	else Vm = bus->vm;
      }
      if(!atleastonegenon) Vm = bus->vm;
    } else Vm = bus->vm;
    Va = bus->va*PETSC_PI/180.0;

    VD = Vm*PetscCosScalar(Va);
    VQ = Vm*PetscSinScalar(Va);
    xarr[loc] = VD;
    xarr[loc+1] = VQ;

    if(bus->nload) {
      PSLOAD load;
      DYNLoadModel dynload;

      for(k=0; k < bus->nload; k++) {
	ierr = PSBUSGetLoad(bus,k,&load);CHKERRQ(ierr);

	if(!load->status) continue;

	ierr = PSLOADGetDYNLoad(load,&dynload);CHKERRQ(ierr);

	if(load->dyniszip) { /* ZIP load model */
	  DYNZIP zip;

	  ierr = DYNLoadModelGetModelData(dynload,(void*)&zip);CHKERRQ(ierr);
	  zip->pl = zip_p[0]*load->pl;
	  zip->ql = zip_q[0]*load->ql;
	  zip->ip = zip_p[1]*load->pl/bus->vm;
	  zip->iq = zip_q[1]*load->ql/bus->vm;
	  zip->yp = zip_p[2]*load->pl/(bus->vm*bus->vm);
	  zip->yq = zip_p[2]*load->ql/(bus->vm*bus->vm);
	  zip->Vm0 = Vm;
	  zip->Vm_thresh = Vm_thresh;
	} else {
	  xbus = xarr+loc;
	  PetscScalar Pl=load->pl,Ql=load->ql;
	  ierr = DYNLoadModelSetInitialConditions(dynload,Pl,Ql,VD,VQ,xbus);CHKERRQ(ierr);
	}
      }
    }

    /* Set initial conditions for dyngen model */
    if(bus->ngen) {
      PSGEN gen;
      DYNGenModel dyngen;
      DYNExcModel dynexc;
      DYNTurbgovModel dynturbgov;
      DYNStabModel dynstab;

      for(k=0; k < bus->ngen; k++) {
        ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);

        if(!gen->status) continue;

        ierr = PSGENGetDYNGen(gen,&dyngen);CHKERRQ(ierr);

        xbus = xarr+loc;
        PetscScalar Pg=gen->pg,Qg=gen->qg;
        ierr = DYNGenModelSetInitialConditions(dyngen,Pg,Qg,VD,VQ,xbus);CHKERRQ(ierr);

	dyngen->freq_max = dyn->freq_max;
	dyngen->freq_min = dyn->freq_min;

	/* Exciter Model */
	if(dyngen->dynexc) {
	  ierr = PSGENGetDYNExc(gen,&dynexc);CHKERRQ(ierr);
	  ierr = DYNExcModelSetInitialConditions(dynexc,VD,VQ,xbus);CHKERRQ(ierr);
	}
	
	/* Turbine governor Model */
	if(dyngen->dynturbgov) {
	  ierr = PSGENGetDYNTurbgov(gen,&dynturbgov);CHKERRQ(ierr);
	  ierr = DYNTurbgovModelSetInitialConditions(dynturbgov,VD,VQ,xbus);CHKERRQ(ierr);
	}

	/* Stabilizer Model */
	if(dyngen->dynexc && dynexc->dynstab) {
	  ierr = PSGENGetDYNStab(gen,&dynstab);CHKERRQ(ierr);
	  ierr = DYNStabModelSetInitialConditions(dynstab,VD,VQ,xbus);CHKERRQ(ierr);
	}

      }
    }
  }

  ierr = VecRestoreArray(localX,&xarr);CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(ps->networkdm,localX,ADD_VALUES,X);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(ps->networkdm,localX,ADD_VALUES,X);CHKERRQ(ierr);

  ierr = DMRestoreLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);

  //  ierr = VecView(X,0);
  PetscFunctionReturn(0);
}

/*
  DYNSetLambdaToUVectors - Initializes the sensitivity variables in the vectors lambda.

  Input Parameters:
. DYN - the dynamics simulation application object
. ncostfcns - number of cost functions

  Output Parameters:
. lambda - the sensitivity vectors
*/
PetscErrorCode DYNSetLambdaToUVectors(DYN dyn,Vec* lambda,PetscInt ncostfcns)
{
  PetscErrorCode ierr;
  PS             ps=dyn->ps;
  PSBUS          bus;
  Vec            localLam;
  PetscBool      ghostbus;
  PetscInt       i,k,loc,veci=0;

  PetscFunctionBegin;
  ierr = DMGetLocalVector(ps->networkdm,&localLam);CHKERRQ(ierr);
  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
    ierr = PSBUSIsGhosted(bus,&ghostbus);CHKERRQ(ierr);
    if (ghostbus) continue;

    if(bus->ide == ISOLATED_BUS) continue;

    for(k=0; k < bus->ngen; k++) {
      PSGEN gen;
      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      if(gen->status) {
        DYNGenModel dyngen;
        PetscScalar *lam,*lamdyn;

        ierr = PSGENGetDYNGen(gen,&dyngen);CHKERRQ(ierr);

        ierr = DMGlobalToLocalBegin(ps->networkdm,lambda[veci],INSERT_VALUES,localLam);CHKERRQ(ierr);
        ierr = DMGlobalToLocalEnd(ps->networkdm,lambda[veci],INSERT_VALUES,localLam);CHKERRQ(ierr);

        ierr = VecGetArray(localLam,&lam);CHKERRQ(ierr);
        lamdyn = lam+loc;
        ierr = DYNGenModelSetFrequency(dyngen,lamdyn,1.0);CHKERRQ(ierr);
        ierr = VecRestoreArray(localLam,&lam);CHKERRQ(ierr);

        ierr = DMLocalToGlobalBegin(ps->networkdm,localLam,INSERT_VALUES,lambda[veci]);CHKERRQ(ierr);
        ierr = DMLocalToGlobalEnd(ps->networkdm,localLam,INSERT_VALUES,lambda[veci]);CHKERRQ(ierr);
      }
      veci++;
      if(veci>ncostfcns) {
        SETERRQ(PETSC_COMM_SELF,0,"exceeds number of cost functions when initialize lambda.");
        exit(1);
      }
    }
  }
  ierr = DMRestoreLocalVector(ps->networkdm,&localLam);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNGetCostFunction - Get the values of cost functions.

  Input Parameters:
. DYN - the dynamics simulation application object
. t   - current time

  Output Parameters:
. lambda - the sensitivity vectors
*/
PetscErrorCode DYNGetCostFunction(DYN dyn,PetscReal t,PetscReal *costfcn,PetscInt ncostfun)
{
  Vec               costfcn_seq;
  VecScatter        ctx;
  PetscScalar       *x;
  PetscErrorCode    ierr;

  PetscFunctionBegin;
  if(dyn->objfun==FREQVIOL) {
    ierr = TSGetCostIntegral(dyn->ts,&dyn->costintegral);CHKERRQ(ierr);
    ierr = VecScatterCreateToAll(dyn->costintegral,&ctx,&costfcn_seq);CHKERRQ(ierr);
    ierr = VecScatterBegin(ctx,dyn->costintegral,costfcn_seq,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
    ierr = VecScatterEnd(ctx,dyn->costintegral,costfcn_seq,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
    ierr = VecScatterDestroy(&ctx);CHKERRQ(ierr);
    ierr = VecGetArray(costfcn_seq,&x);CHKERRQ(ierr);
    ierr = PetscMemcpy(costfcn,x,ncostfun*sizeof(PetscScalar));CHKERRQ(ierr);
    ierr = VecRestoreArray(costfcn_seq,&x);CHKERRQ(ierr);
    ierr = VecDestroy(&costfcn_seq);CHKERRQ(ierr);
  }else if (dyn->objfun==FREQ) {
    IS       isgenfreq;
    PetscInt Ngen;

    ierr = PSGetNumGenerators(dyn->ps,NULL,&Ngen);CHKERRQ(ierr);
    ierr = DYNSetUpGenFreqIS(dyn,&isgenfreq);CHKERRQ(ierr);
    ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,Ngen,costfcn,&costfcn_seq);CHKERRQ(ierr);
    ierr = VecScatterCreate(dyn->X,isgenfreq,costfcn_seq,NULL,&ctx);CHKERRQ(ierr);
    ierr = VecScatterBegin(ctx,dyn->X,costfcn_seq,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
    ierr = VecScatterEnd(ctx,dyn->X,costfcn_seq,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
    ierr = VecScatterDestroy(&ctx);CHKERRQ(ierr);
    ierr = VecDestroy(&costfcn_seq);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*
  DYNSetUpDifferentialEqIS - Sets up the index set for the differential equations for the DYN

  Input Parameters:
. dyn - the dyn object

  Output Parameters:
. isdiff - index set for differential equations

  Notes:
   The IS contains the global indices
*/
PetscErrorCode DYNSetUpDifferentialEqIS(DYN dyn, IS *isdiff)
{
  PetscErrorCode ierr;
  IS             is;
  PS             ps=dyn->ps;
  PSBUS          bus;
  PetscInt       *idx,ctr=0,i,locglob,k;
  PetscBool      ghostbus;

  PetscFunctionBegin;

  ierr = PetscCalloc1(dyn->ndiff,&idx);CHKERRQ(ierr);
  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    ierr = PSBUSIsGhosted(bus,&ghostbus);CHKERRQ(ierr);
    if(!ghostbus) {
      ierr = PSBUSGetVariableGlobalLocation(bus,&locglob);CHKERRQ(ierr);

      for(k=0; k < bus->ngen; k++) {
	PSGEN gen;
	DYNGenModel dyngen;
	DYNExcModel dynexc;
	DYNTurbgovModel dynturbgov;
	DYNStabModel    dynstab;
	PetscInt    j;

	ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
	if(!gen->status) continue;

	ierr = PSGENGetDYNGen(gen,&dyngen);CHKERRQ(ierr);
	for(j=0; j < dyngen->nvar; j++) {
	  if(dyngen->eqtypes[j] == DIFF_EQ) idx[ctr++] = locglob + dyngen->startloc + j;
	}
	
	/* Exciter model */
	if(dyngen->dynexc) {
	  ierr = PSGENGetDYNExc(gen,&dynexc);CHKERRQ(ierr);
	  for(j=0; j < dynexc->nvar; j++) {
	    if(dynexc->eqtypes[j] == DIFF_EQ) idx[ctr++] = locglob + dynexc->startloc + j;
	  }
	}

	/* Turbine Governor model */
	if(dyngen->dynturbgov) {
	  ierr = PSGENGetDYNTurbgov(gen,&dynturbgov);CHKERRQ(ierr);
	  for(j=0; j < dynturbgov->nvar; j++) {
	    if(dynturbgov->eqtypes[j] == DIFF_EQ) idx[ctr++] = locglob + dynturbgov->startloc + j;
	  }
	}

	/* Stabilizer model */
	if(dyngen->dynexc && dynexc->dynstab) {
	  ierr = PSGENGetDYNStab(gen,&dynstab);CHKERRQ(ierr);
	  for(j=0; j < dynstab->nvar; j++) {
	    if(dynstab->eqtypes[j] == DIFF_EQ) idx[ctr++] = locglob + dynstab->startloc + j;
	  }
	}
      }

      for(k=0; k < bus->nload; k++) {
	PSLOAD load;
	DYNLoadModel dynload;
	PetscInt    j;

	ierr = PSBUSGetLoad(bus,k,&load);CHKERRQ(ierr);
	if(!load->status) continue;

	ierr = PSLOADGetDYNLoad(load,&dynload);CHKERRQ(ierr);
	for(j=0; j < dynload->nvar; j++) {
	  if(dynload->eqtypes[j] == DIFF_EQ) idx[ctr++] = locglob + dynload->startloc + j;
	}
      }
    }
  }

  /* Create the diff. eq. IS */
  ierr = ISCreateGeneral(dyn->comm->type,dyn->ndiff,idx,PETSC_COPY_VALUES,&is);CHKERRQ(ierr);
  ierr = PetscFree(idx);CHKERRQ(ierr);

  *isdiff = is;
  PetscFunctionReturn(0);
}

/*
  DYNPostStep - Post-step callback for DYN

  Notes: All the callback function set for DYN are called via
         this function
*/ 
PetscErrorCode DYNPostStep(TS ts)
{
  DYN dyn;
  PetscErrorCode ierr;
  PetscInt       i;

  PetscFunctionBegin;
  ierr = TSGetApplicationContext(ts,&dyn);CHKERRQ(ierr);
  for(i=0; i < dyn->npoststepfcns; i++) {
    ierr = (dyn->poststepfcn[i])(dyn);
  }
  PetscFunctionReturn(0);
}

/*
  DYNPreStep - Pre-step callback for DYN

  Notes: All the callback function set for DYN are called via
         this function
*/ 
PetscErrorCode DYNPreStep(TS ts)
{
  DYN dyn;
  PetscErrorCode ierr;
  PetscInt       i;
  PetscReal      t;
  Vec            X;

  PetscFunctionBegin;
  ierr = TSGetApplicationContext(ts,&dyn);CHKERRQ(ierr);
  ierr = TSGetTime(ts,&t);CHKERRQ(ierr);
  ierr = TSGetSolution(ts,&X);CHKERRQ(ierr);

  dyn->switching_solution = PETSC_FALSE;

  /* Store the event function values before X changes */
  ierr = DYNEventMonitor(ts,t,X,dyn->feventtmp,(void*)dyn);CHKERRQ(ierr);

  for(i=0; i < dyn->nprestepfcns; i++) {
    ierr = (dyn->prestepfcn[i])(dyn);
  }
  PetscFunctionReturn(0);
}

/*
  DYNSetUp - Sets up a DYN object

  Input Parameters:
. DYN - the DYN object

  Notes:
  This routine sets up the DYN object and the underlying PS object. It
  also distributes the PS object when used in parallel.
*/
PetscErrorCode DYNSetUp(DYN dyn)
{
  PetscErrorCode ierr;
  PS             ps=dyn->ps;
  DIR            *dir;
  PetscBool      flg;
  PetscReal      atol=1e-2,rtol=1e-2;
  const          PetscInt *idx;
  PetscInt       i,lsize;
  PetscBool      use_diff_var_only=PETSC_FALSE;
  PetscInt       rstart;
  PetscScalar    *vatolarr,*vrtolarr;

  /*DYNstart File Variables (SAM)*/
  char           start_file[PETSC_MAX_PATH_LEN];
  PetscViewer    startview;

  PetscFunctionBegin;
  dyn->useforward = dyn->useadjoint = PETSC_FALSE;

  /* Set up static load models */
  ierr = DYNSetUpStaticLoadModels(dyn);CHKERRQ(ierr);

  /* Set up constant voltage sources for generators with no dynamic data */
  ierr = DYNSetUpConstantVoltageGenModels(dyn);CHKERRQ(ierr);

  /* Set up PS object */
  ierr = PSSetUp(ps);CHKERRQ(ierr);

  /* Set up differential equation IS */
  dyn->ndiff = ps->ndiff;
  ierr = DYNSetUpDifferentialEqIS(dyn,&dyn->isdiff);CHKERRQ(ierr);

  /* Set up events */
  ierr = DYNSetUpEvents(dyn);CHKERRQ(ierr);

  /* Create the solution vector and the Jacobian matrix */
  ierr = PSCreateGlobalVector(dyn->ps,&dyn->X);CHKERRQ(ierr);
  ierr = PSCreateMatrix(dyn->ps,&dyn->Jac);CHKERRQ(ierr);

  /* Associate the DM object in PS with the time-stepping solver */
  ierr = TSSetDM(dyn->ts,dyn->ps->networkdm);CHKERRQ(ierr);

  /* Set the prefix for this TS.. All TS runtime options will need to have the prefix "-dyn_" */
  ierr = TSSetOptionsPrefix(dyn->ts,"dyn_");CHKERRQ(ierr);

  ierr = PetscOptionsGetBool(NULL,NULL,"-dyn_use_semiexplicit",&dyn->use_semiexplicit,NULL);CHKERRQ(ierr);

  if(dyn->use_semiexplicit) {
    ierr = TSSetType(dyn->ts,TSRK);CHKERRQ(ierr);
    ierr = TSSetRHSFunction(dyn->ts,NULL,DYNRHSFunction,(void*)dyn);CHKERRQ(ierr);
    ierr = TSSetPostStage(dyn->ts,DYNPostStage);CHKERRQ(ierr);
    ierr = TSSetPostEvaluate(dyn->ts,DYNPostEvaluate);CHKERRQ(ierr);
    ierr = TSSetApplicationContext(dyn->ts,(void*)dyn);CHKERRQ(ierr);
  } else {
    /* Default solver - implicit trapezoidal (set in options file) */
    ierr = TSSetProblemType(dyn->ts,TS_NONLINEAR);CHKERRQ(ierr);
    ierr = TSSetEquationType(dyn->ts,TS_EQ_IMPLICIT);CHKERRQ(ierr);
    ierr = TSARKIMEXSetFullyImplicit(dyn->ts,PETSC_TRUE);CHKERRQ(ierr);
    //    ierr = TSSetType(dyn->ts,TSBEULER);CHKERRQ(ierr);

    /* Set the function and Jacobian routines */
    ierr = TSSetIFunction(dyn->ts,NULL,DYNIFunction,(void*)dyn);CHKERRQ(ierr);
    ierr = TSSetIJacobian(dyn->ts,dyn->Jac,dyn->Jac,(TSIJacobian)DYNIJacobian,(void*)dyn);CHKERRQ(ierr);
  }

  /* Create vector of tolerances for adaptivity */
  ierr = PSCreateGlobalVector(dyn->ps,&dyn->vatol);CHKERRQ(ierr);
  ierr = PSCreateGlobalVector(dyn->ps,&dyn->vrtol);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-dyn_ts_atol",&atol,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-dyn_ts_rtol",&rtol,NULL);CHKERRQ(ierr);
  ierr = VecSet(dyn->vatol,atol);CHKERRQ(ierr);
  ierr = VecSet(dyn->vrtol,rtol);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-dyn_ts_adapt_use_diff_var_only",&use_diff_var_only,NULL);CHKERRQ(ierr);
  if(use_diff_var_only) {
    /* Set large values in the entire vector and then set the given tolerances only for the differential variables */
    ierr = VecSet(dyn->vatol,99999);CHKERRQ(ierr);
    ierr = VecSet(dyn->vrtol,99999);CHKERRQ(ierr);
    ierr = VecGetOwnershipRange(dyn->vrtol,&rstart,NULL);CHKERRQ(ierr);
    ierr = VecGetArray(dyn->vatol,&vatolarr);CHKERRQ(ierr);
    ierr = VecGetArray(dyn->vrtol,&vrtolarr);CHKERRQ(ierr);

    ierr = ISGetLocalSize(dyn->isdiff,&lsize);CHKERRQ(ierr);
    ierr = ISGetIndices(dyn->isdiff,&idx);CHKERRQ(ierr);

    for(i=0; i < lsize; i++) {
      vatolarr[idx[i]-rstart] = atol;
      vrtolarr[idx[i]-rstart] = rtol;
    }

    ierr = VecRestoreArray(dyn->vatol,&vatolarr);CHKERRQ(ierr);
    ierr = VecRestoreArray(dyn->vrtol,&vrtolarr);CHKERRQ(ierr);
  }
  ierr = TSSetTolerances(dyn->ts,PETSC_DECIDE,dyn->vatol,PETSC_DECIDE,dyn->vrtol);CHKERRQ(ierr);

  /* Stop when the final time is reached */
  ierr = TSSetExactFinalTime(dyn->ts,TS_EXACTFINALTIME_MATCHSTEP);CHKERRQ(ierr);

  if(dyn->eval_integral) {
    PetscInt Ngen;
    ierr = PSGetNumGenerators(dyn->ps,NULL,&Ngen);CHKERRQ(ierr);
    dyn->ncostfcns = Ngen;
    ierr = TSSetCostIntegrand(dyn->ts,dyn->ncostfcns,dyn->costintegral,(PetscErrorCode (*)(TS,PetscReal,Vec,Vec,void*))DYNCostIntegrand,PETSC_NULL,PETSC_NULL,PETSC_TRUE,(void*)dyn);CHKERRQ(ierr);
  }
  /* Stop the code breaking if TS fails to converge */
  ierr = TSSetErrorIfStepFails(dyn->ts,PETSC_FALSE);CHKERRQ(ierr);

  ierr = TSSetFromOptions(dyn->ts);CHKERRQ(ierr);

  ierr = TSSetApplicationContext(dyn->ts,dyn);CHKERRQ(ierr);

  ierr = PetscOptionsGetBool(NULL,NULL,"-dyn_viewer",&dyn->visualize,NULL);CHKERRQ(ierr);

  if(dyn->visualize) {
    ierr = DYNSetPostStepCallback(dyn,DYNViewer);CHKERRQ(ierr);
    ierr = PetscOptionsGetString(NULL,NULL,"-dyn_viewer_dir",dyn->viewer_dir,sizeof(dyn->viewer_dir),&flg);CHKERRQ(ierr);
    if(flg) {
      dir = opendir(dyn->viewer_dir);
      if(!dir) { /* Create directory if it does not exist */
	    ierr = PetscMkdir(dyn->viewer_dir);CHKERRQ(ierr);
      }
    } else {
      ierr = PetscRMTree("dyn-output");CHKERRQ(ierr);
      ierr = PetscMkdir("dyn-output");CHKERRQ(ierr);
      ierr = PetscStrcpy(dyn->viewer_dir,"dyn-output");CHKERRQ(ierr);
    }
    /*Create DYNstart File (SAM)*/
    ierr = PetscSNPrintf(start_file,sizeof(start_file),"%s/DYNstart.out",dyn->viewer_dir);CHKERRQ(ierr);
    ierr = PetscViewerASCIIOpen(dyn->comm->type,start_file,&startview);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(startview,"DYN Start\n");CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&startview);CHKERRQ(ierr);
  }

  /* Set the post step function for TS through which all the DYNPostStep functions will be called */
  if(dyn->npoststepfcns > 0) {
    ierr = TSSetPostStep(dyn->ts,DYNPostStep);CHKERRQ(ierr);
  }

  /* Set the pre step function for TS through which all the DYNPreStep functions will be called */
  /* Note the prestep function all resets the switching solution flag */
  ierr = TSSetPreStep(dyn->ts,DYNPreStep);CHKERRQ(ierr);

  /* Set up the event monitoring function for TS */
  ierr = TSSetEventHandler(dyn->ts,dyn->toteventmonitors,dyn->eventdirection,dyn->eventterminate,DYNEventMonitor,DYNEventPostFunction,(void*)dyn);CHKERRQ(ierr);

  /* Set up parameter sensitivity calculation options*/
  ierr = DYNGetCostFunctionParametersFromOptions(dyn);CHKERRQ(ierr);

  if(dyn->prepflow) {
    ierr = DYNSetUpInitPflow(dyn);CHKERRQ(ierr);
  }

  ierr = PetscLogEventRegister("DYNSolve",TS_CLASSID,&dyn->logdynsolve);CHKERRQ(ierr);
  dyn->setupcalled = PETSC_TRUE;

  PetscFunctionReturn(0);
}

/*
  DYNSolve - Runs the dynamic simulation

  Input Parameters:
. dyn - the dynamics simulation object
*/
PetscErrorCode DYNSolve(DYN dyn)
{
  PetscErrorCode ierr;
  TSConvergedReason reason;
  PetscBool pflowconverged;

  /*DYNend File Variables (SAM)*/
  char           end_file[PETSC_MAX_PATH_LEN];
  PetscViewer    endview;
  PetscInt       stepnum;

  PetscFunctionBegin;
  if(!dyn->setupcalled) {
    ierr = DYNSetUp(dyn);CHKERRQ(ierr);
  }

  /* Check Topology */
  ierr = PSCheckTopology(dyn->ps);CHKERRQ(ierr);
  

  if(dyn->eval_integral) {
    Vec C;
    ierr = TSGetCostIntegral(dyn->ts,&C);CHKERRQ(ierr);
    ierr = VecSet(C,0);CHKERRQ(ierr);
  }
  /* Initial steady-state power flow */
  if(dyn->prepflow) {
    ierr = DYNComputePrePflow(dyn,&pflowconverged);CHKERRQ(ierr);
  }

  /* Set up initial conditions */
  ierr = DYNSetInitialConditions(dyn,dyn->X);CHKERRQ(ierr);

  /* Set solution vector with TS */
  ierr = TSSetSolution(dyn->ts,dyn->X);CHKERRQ(ierr);

  PetscLogEventBegin(dyn->logdynsolve,0,0,0,0);
  /* Call the time-stepping solver */
  ierr = TSSolve(dyn->ts,dyn->X);CHKERRQ(ierr);
  PetscLogEventEnd(dyn->logdynsolve,0,0,0,0);

  /* Get converged reason */
  ierr = TSGetConvergedReason(dyn->ts,&reason);CHKERRQ(ierr);
  if (reason<0) {
    dyn->tsconverged = PETSC_FALSE;
  } else dyn->tsconverged = PETSC_TRUE;

  /* Get the actual finish time */
  ierr = TSGetSolveTime(dyn->ts,&dyn->finaltime);CHKERRQ(ierr);

  if(dyn->visualize) {
    /*Create DYNend File (SAM)*/
    ierr = PetscSNPrintf(end_file,sizeof(end_file),"%s/DYNend.out",dyn->viewer_dir);CHKERRQ(ierr);
    ierr = PetscViewerASCIIOpen(dyn->comm->type,end_file,&endview);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(endview,"Type\tLength\n");CHKERRQ(ierr);
    ierr = TSGetStepNumber(dyn->ts,&stepnum);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(endview,"vol,%d\n", stepnum);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(endview,"freq,%d\n", stepnum);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(endview,"log,%d\n", stepnum);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(endview,"event,%d\n", dyn->eventlogcount);CHKERRQ(ierr);
    //TODO: Add Event: ierr = PetscViewerASCIIPrintf(endview,"event,%d\n", stepnum);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&endview);CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

/*
  DYNAdjointSetUp - Sets up a DYN object

  Input Parameters:
. DYN - the DYN object

  Notes:
  This routine sets up the DYN object and the underlying PS object. It
  also distributes the PS object when used in parallel.
*/
PetscErrorCode DYNAdjointSetUp(DYN dyn)
{
  PetscErrorCode ierr;
  IS             isgenfreq;
  PS             ps=dyn->ps;
  PetscInt       i,Nbus,nghostbuses,Ngen,ngen,xdynsize; // local varaibles
  PetscBool      isghosted;

  PetscFunctionBegin;
  /* Set up static load models */
  ierr = DYNSetUpStaticLoadModels(dyn);CHKERRQ(ierr);

  /* Set up PS object */
  ierr = PSSetUp(ps);CHKERRQ(ierr);

  /* Set up differential equation IS */
  dyn->ndiff = ps->ndiff;
  ierr = DYNSetUpDifferentialEqIS(dyn,&dyn->isdiff);CHKERRQ(ierr);

  /* Create the solution vector and the Jacobian matrix */
  ierr = PSCreateGlobalVector(dyn->ps,&dyn->X);CHKERRQ(ierr);
  ierr = PSCreateMatrix(dyn->ps,&dyn->Jac);CHKERRQ(ierr);

  /* Associate the DM object in PS with the time-stepping solver */
  ierr = TSSetDM(dyn->ts,dyn->ps->networkdm);CHKERRQ(ierr);
  /* Default solver - implicit trapezoidal (can be changed via options) */
  ierr = TSSetType(dyn->ts,TSBEULER);CHKERRQ(ierr);

  /* Set the prefix for this TS.. All TS runtime options will need to have the prefix "-dyn_" */
  ierr = TSSetOptionsPrefix(dyn->ts,"dyn_");CHKERRQ(ierr);

  /* Set the function and Jacobian routines */
  ierr = TSSetIFunction(dyn->ts,NULL,DYNIFunction,(void*)dyn);CHKERRQ(ierr);
  ierr = TSSetIJacobian(dyn->ts,dyn->Jac,dyn->Jac,(TSIJacobian)DYNIJacobian,(void*)dyn);CHKERRQ(ierr);

  /* Stop when the final time is reached */
  ierr = TSSetExactFinalTime(dyn->ts,TS_EXACTFINALTIME_MATCHSTEP);CHKERRQ(ierr);

  /* Set up parameter sensitivity calculation options */
  ierr = DYNGetCostFunctionParametersFromOptions(dyn);CHKERRQ(ierr);

  ierr = PSGetNumGenerators(dyn->ps,&ngen,&Ngen);CHKERRQ(ierr);
  nghostbuses = 0;
  for (i=0;i<dyn->ps->nbus;i++) {
    ierr = PSBUSIsGhosted(&dyn->ps->bus[i],&isghosted);CHKERRQ(ierr);
    if(isghosted) nghostbuses++;
  }
  dyn->nparams = 2*(dyn->ps->nbus-nghostbuses) + 2*ngen; // number of local parameters (excluding those for ghost buses)
  ierr = PSGetNumGlobalBuses(dyn->ps,&Nbus);CHKERRQ(ierr);
  dyn->Nparams = 2*Nbus+2*Ngen; // number of global parameters
  if(dyn->sum_freq_costfun) dyn->ncostfcns = 1;
  else {
    dyn->ncostfcns = Ngen;
  }

  /* Create the Jacobian of the DAE equations w.r.t. parameters */
  ierr = VecGetLocalSize(dyn->X,&xdynsize);CHKERRQ(ierr);
  ierr = MatCreate(dyn->ps->comm->type,&dyn->Jacp);CHKERRQ(ierr);
  ierr = MatSetSizes(dyn->Jacp,xdynsize,dyn->nparams,PETSC_DETERMINE,PETSC_DETERMINE);CHKERRQ(ierr);
  ierr = MatSetFromOptions(dyn->Jacp);CHKERRQ(ierr);
  ierr = MatSetUp(dyn->Jacp);CHKERRQ(ierr); /* Skipping allocation for now */

  ierr = MatCreate(dyn->ps->comm->type,&dyn->ICp);CHKERRQ(ierr);
  ierr = MatSetSizes(dyn->ICp,xdynsize,dyn->nparams,PETSC_DETERMINE,PETSC_DETERMINE);CHKERRQ(ierr);
  ierr = MatSetFromOptions(dyn->ICp);CHKERRQ(ierr);
  ierr = MatSetUp(dyn->ICp);CHKERRQ(ierr); /* Skipping allocation for now */

  ierr = MatCreateVecs(dyn->Jacp,&dyn->Xparam,NULL);CHKERRQ(ierr);

  /* Set up vectors for sensitivities */
  //ierr = VecCreateSeq(PETSC_COMM_SELF,dyn->nparams,&dyn->Xparam);CHKERRQ(ierr);

  ierr = VecDuplicateVecs(dyn->X,dyn->ncostfcns,&dyn->lambda);CHKERRQ(ierr);
  ierr = VecDuplicateVecs(dyn->Xparam,dyn->ncostfcns,&dyn->mu);CHKERRQ(ierr);
  ierr = VecDuplicateVecs(dyn->X,dyn->Nparams,&dyn->dy0dp);CHKERRQ(ierr);

  ierr = VecDuplicate(dyn->X,&dyn->Xtmp);CHKERRQ(ierr);
  ierr = VecDuplicate(dyn->X,&dyn->Xtmp2);CHKERRQ(ierr);
  ierr = VecDuplicate(dyn->X,&dyn->Xtmp3);CHKERRQ(ierr); // costintegrand workpalce

  /* Set the number of cost gradients for dyn->ts */
  ierr = TSSetCostGradients(dyn->ts,dyn->ncostfcns,dyn->lambda,dyn->mu);CHKERRQ(ierr);

  /* Set functions needed by TS for adjoint calculation */
  ierr = TSSetSaveTrajectory(dyn->ts);CHKERRQ(ierr);
  ierr = TSAdjointSetRHSJacobian(dyn->ts,dyn->Jacp,DYNComputeAdjointJacobianP,(void*)dyn);CHKERRQ(ierr);
  if(dyn->objfun == FREQVIOL) {
    ierr = DYNSetUpGenFreqIS(dyn,&isgenfreq);CHKERRQ(ierr);
    ierr = VecGetSubVector(dyn->Xtmp3,isgenfreq,&dyn->costintegral);CHKERRQ(ierr);
    ierr = TSSetCostIntegrand(dyn->ts,dyn->ncostfcns,dyn->costintegral,(PetscErrorCode (*)(TS,PetscReal,Vec,Vec,void*))DYNCostIntegrand,
                 (PetscErrorCode (*)(TS,PetscReal,Vec,Vec*,void*))DYNDRDYFunction,
                 (PetscErrorCode (*)(TS,PetscReal,Vec,Vec*,void*))DYNDRDPFunction,PETSC_TRUE,(void*)dyn);CHKERRQ(ierr);
    ierr = VecRestoreSubVector(dyn->Xtmp3,isgenfreq,&dyn->costintegral);CHKERRQ(ierr);
  }

  /* Set up events */
  ierr = DYNSetUpEvents(dyn);CHKERRQ(ierr);

  /* Set up the event monitoring function for TS */
  ierr = TSSetEventHandler(dyn->ts,dyn->toteventmonitors,dyn->eventdirection,dyn->eventterminate,DYNEventMonitor,DYNEventPostFunction,(void*)dyn);CHKERRQ(ierr);

  //ierr = TSAdjointMonitorSet(dyn->ts,DYNAdjointMonitor,(void *)dyn,NULL);CHKERRQ(ierr);
  if(dyn->monitor) {
    ierr = TSMonitorSet(dyn->ts,DYNUserMonitor,(void *)dyn,NULL);CHKERRQ(ierr);
  }

  /* Stop the code breaking if TS fails to converge */
  ierr = TSSetErrorIfStepFails(dyn->ts,PETSC_FALSE);CHKERRQ(ierr);

  ierr = TSSetFromOptions(dyn->ts);CHKERRQ(ierr);

  if(dyn->prepflow) {
    ierr = DYNSetUpInitPflow(dyn);CHKERRQ(ierr);
  }

  dyn->setupcalled = PETSC_TRUE;
  PetscFunctionReturn(0);
}

/*
  DYNAdjointSolve - Perform adjoint sensitivity analysis of the dynamic simulation

  Input Parameters:
. dyn - the dynamics simulation object

  Notes:
   The parameters considered are the bus voltage magnitudes and angles, and the generator real and reactive power
   output. For each bus, its parameters are ordered as [Va, Vm, PG1, QG1, PG2, QG2,.., PGn, QGn]
*/
PetscErrorCode DYNAdjointSolve(DYN dyn)
{
  PetscErrorCode    ierr;
  PetscInt          i;
  FILE              *f;
  TSConvergedReason reason;

  PetscFunctionBegin;
  if(!dyn->setupcalled) {
    ierr = DYNAdjointSetUp(dyn);CHKERRQ(ierr);
  }

  /* Initial steady-state power flow */
  if(dyn->prepflow) {
    ierr = DYNComputePrePflow(dyn,&dyn->prepflowconverged);CHKERRQ(ierr);
    if (!dyn->prepflowconverged) PetscFunctionReturn(0);
  }

  /* Set up initial conditions */
  ierr = DYNSetInitialConditions(dyn,dyn->X);CHKERRQ(ierr);

  /* Set solution vector with TS */
  ierr = TSSetSolution(dyn->ts,dyn->X);CHKERRQ(ierr);

  /* Call the time-stepping solver */
  ierr = TSSolve(dyn->ts,dyn->X);CHKERRQ(ierr);

  /* Get converged reason */
  ierr = TSGetConvergedReason(dyn->ts,&reason);CHKERRQ(ierr);
  if (reason<0) {
    dyn->tsconverged = PETSC_FALSE;
  } else dyn->tsconverged = PETSC_TRUE;

  /* Get the actual finish time */
  ierr = TSGetSolveTime(dyn->ts,&dyn->finaltime);CHKERRQ(ierr);

  /* Compute parameter sensitivities */
  /* Initialize parameter sensitivity work vectors */
  for(i=0; i < dyn->ncostfcns; i++) {
    ierr = VecSet(dyn->lambda[i],0.0);CHKERRQ(ierr);
    ierr = VecSet(dyn->mu[i],0.0);CHKERRQ(ierr);
  }

  if(dyn->objfun==FREQ) {
    /* initaialize lambda as unit vectors */
    ierr = DYNSetLambdaToUVectors(dyn,dyn->lambda,dyn->ncostfcns);CHKERRQ(ierr);
  }

  /*
  Vec C;
  ierr = TSGetCostIntegral(dyn->ts,&C);CHKERRQ(ierr);
  ierr = VecSet(C,0.0);CHKERRQ(ierr);
  */

  if(dyn->printsensi){
    f = fopen("adjsensi.data", "a");
    ierr = PetscFPrintf(PETSC_COMM_WORLD,f,"%24.10f ",dyn->finaltime);CHKERRQ(ierr);
  }
  /* Do Adjoint solve */
  ierr = TSAdjointSolve(dyn->ts);CHKERRQ(ierr);

  ierr = DYNComputeSensiP2(dyn);CHKERRQ(ierr);
  //ierr = DYNComputeSensiP(dyn);CHKERRQ(ierr);

#if defined DEBUGDYN
  PetscScalar *costfcns;
  ierr = PetscMalloc(dyn->ncostfcns*sizeof(PetscScalar),&costfcns);CHKERRQ(ierr);

  /* Get cost functions */
  ierr = DYNGetCostFunction(dyn,dyn->finaltime,costfcns,dyn->ncostfcns);CHKERRQ(ierr);
  for(i=0;i<dyn->ncostfcns;i++) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"%24.15e ",costfcns[i]);CHKERRQ(ierr);
  }
  ierr = PetscFree(costfcns);CHKERRQ(ierr);
  for(i=0;i<dyn->ncostfcns;i++) {
    ierr = VecView(dyn->mu[i],0);CHKERRQ(ierr);
  }
#endif

  if(dyn->printsensi){
    const PetscScalar *s;
    PetscInt j;
    for(i=0;i<dyn->ncostfcns;i++) {
      ierr = VecGetArrayRead(dyn->mu[i],&s);CHKERRQ(ierr);
      for(j=0;j<dyn->nparams;j++) {
        ierr = PetscFPrintf(PETSC_COMM_WORLD,f,"%24.15e ",s[j]);CHKERRQ(ierr);
      }
      ierr = VecRestoreArrayRead(dyn->mu[i],&s);CHKERRQ(ierr);
    }
    ierr = PetscFPrintf(PETSC_COMM_WORLD,f,"\n");CHKERRQ(ierr);
  }
  /*
  for(i=0;i<dyn->ncostfcns;i++) {
    ierr = VecView(dyn->lambda[i],0);CHKERRQ(ierr);
    //ierr = CheckAlgebraicPart(dyn,dyn->lambda[i]);CHKERRQ(ierr);
  }

  PetscViewer viewer;
  Vec         *DCDP;
  PetscReal   distance = 0,norm1,norm2;
  char        filename[PETSC_MAX_PATH_LEN] = "finitediff_results.bin";
  ierr = VecDuplicateVecs(dyn->mu[0],dyn->ncostfcns,&DCDP);CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
  for(i=0;i<dyn->ncostfcns;i++) {
    ierr = VecLoad(DCDP[i],viewer);CHKERRQ(ierr);
    ierr = VecNorm(DCDP[i],NORM_2,&norm2);CHKERRQ(ierr);
    ierr = VecAXPY(DCDP[i],-1,dyn->mu[i]);CHKERRQ(ierr);
    ierr = VecNorm(DCDP[i],NORM_2,&norm1);CHKERRQ(ierr);
    if (norm2>0) distance += PetscSqrtScalar(norm1/norm2);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"norm1=%lf norm2=%lf\n",norm1,norm2);CHKERRQ(ierr);
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,"distance to finitediff results: %lf\n",distance);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  ierr = VecDestroyVecs(dyn->ncostfcns,&DCDP);CHKERRQ(ierr);
  */
  PetscFunctionReturn(0);
}

#ifdef FWDSA

/*
 DYNComputeForwardJacobianP - Computes the Jacobian of the DYN equations w.r.t. parameters for the forward model

 Input Parameters:
 + ts - the TS solver
 . t  - the current time
 . x  - the solution vector
 - ctx - application context (dyn)

 Output Parameters
 . jacP - Jacobian of the DAE equations w.r.t. parameters
 */
PetscErrorCode DYNComputeForwardJacobianP(TS ts,PetscReal t,Vec X,Vec *jacP,void *ctx)
{
  PetscErrorCode ierr;
  DYN            dyn=(DYN)ctx;
  const          PetscScalar *x;
  Vec            localX;
  PetscInt       i;
  PS             dynps=dyn->ps;
  PetscBool      ghostbus;
  PetscInt       paramloc;

  PetscFunctionBegin;
  for (i=0;i<dyn->nparams;i++) {
    ierr = VecSet(jacP[i],0.0);CHKERRQ(ierr);
  }
  /* Get DYN local vector */
  ierr = DMGetLocalVector(dynps->networkdm,&localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(dynps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(dynps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);

  ierr = VecGetOwnershipRange(dyn->Xparam,&paramloc,NULL);CHKERRQ(ierr);

  ierr = VecGetArrayRead(localX,&x);CHKERRQ(ierr);

  for(i=0; i < dynps->nbus; i++) {
    PetscScalar VD,VQ;
    PetscInt    dynloc;
    PSBUS       dynbus;
    PetscInt    k,PGloc,QGloc,VAloc,VMloc,ctr=0;
    //PetscScalar Vm,Vm0,Vm2;
    PetscInt    row[2];
    //PetscInt    col[1];
    PetscScalar val[2];

    dynbus = &dynps->bus[i];
    ierr = PSBUSGetVariableGlobalLocation(dynbus,&dynloc);CHKERRQ(ierr);
    ierr = PSBUSIsGhosted(dynbus,&ghostbus);CHKERRQ(ierr);
    if (ghostbus) continue;

    if (dynbus->ide == ISOLATED_BUS) continue;

    VAloc = paramloc; VMloc = paramloc+1; paramloc += 2;

    VD = x[dynloc]; VQ = x[dynloc+1];
    //Vm = PetscSqrtScalar(VD*VD + VQ*VQ);
    //Vm2 = Vm*Vm;

    for(k=0; k < dynbus->nload; k++) {
      //Vm0 = dynbus->vm;
      PSLOAD load;
      ierr = PSBUSGetLoad(dynbus,k,&load);CHKERRQ(ierr);
      if(!load->status) continue;
      PetscScalar Pd,Qd,dyp_dVm0,dyq_dVm0;
      Pd = load->pl;
      Qd = load->ql;

      dyp_dVm0 = 2*Pd/(dynbus->vm*dynbus->vm*dynbus->vm);
      dyq_dVm0 = 2*Qd/(dynbus->vm*dynbus->vm*dynbus->vm);

      row[0] = dynloc; row[1] = dynloc+1;
      //col[0] = VMloc;
      val[0] = dyq_dVm0*VD - dyp_dVm0*VQ;
      val[1] = -dyp_dVm0*VD - dyq_dVm0*VQ;

      ierr = VecSetValues(jacP[VMloc],2,row,val,INSERT_VALUES);CHKERRQ(ierr);
    }

    /* Generator frequency deviation */
    for(k=0; k < dynbus->ngen; k++) {
      PSGEN             gen;
      DYNGenModel       dyngen;
      const PetscScalar *xdyn;
      PetscInt          startloc;
      PetscInt          dynlocglob;

      ierr = PSBUSGetGen(dynbus,k,&gen);CHKERRQ(ierr);
      ierr = PSGENGetDYNGen(gen,&dyngen);CHKERRQ(ierr);

      startloc   = dynloc;
      xdyn       = x+startloc;
      PGloc      = paramloc + ctr;
      QGloc      = paramloc + ctr + 1;
      ctr       += 2; // keep PG and QG even when the generator is off
      dynlocglob = startloc + dyngen->startloc;
      if(!gen->status) continue;
        /* Set partial derivatives of machine DAE equations w.r.t parameters */
      ierr = DYNGenModelDAEFWDRHSJacobianP(dyngen,t,xdyn,jacP,dynlocglob,VAloc,VMloc,PGloc,QGloc);CHKERRQ(ierr);

      if(dyngen->dynexc) {
        DYNExcModel dynexc;
        ierr = PSGENGetDYNExc(gen,&dynexc);CHKERRQ(ierr);
        dynlocglob = startloc + dynexc->startloc;
        ierr = DYNExcModelDAEFWDRHSJacobianP(dynexc,t,xdyn,jacP,dynlocglob,VAloc,VMloc);CHKERRQ(ierr);
      }
    }
    paramloc += ctr;
  }

  ierr = VecRestoreArrayRead(localX,&x);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(dynps->networkdm,&localX);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
 DYNForwardAlgebraicSolve - Solve the algebraic part to obtain consistent sensitivity variables
 
 Input Parameters:
 . dyn - the dynamics simulation object
 */
PetscErrorCode DYNForwardAlgebraicSolve(SNES tssnes,DYN dyn,PetscReal t,Vec X)
{
    PetscErrorCode ierr;
    const PetscInt *idx;
    PetscInt       i,lsize;
    const PetscScalar    *v;
    KSP            ksp;
    Mat            J,Jp;
    Vec            subv;
    
    PetscFunctionBegin;
    ierr = ISGetLocalSize(dyn->isdiff,&lsize);CHKERRQ(ierr);
    ierr = ISGetIndices(dyn->isdiff,&idx);CHKERRQ(ierr);
    
    ierr = SNESGetJacobian(tssnes,&J,&Jp,NULL,NULL);CHKERRQ(ierr);
    /* Solve J Y = F */
    ierr = SNESComputeJacobian(tssnes,X,J,Jp);CHKERRQ(ierr);
    ierr = SNESGetKSP(tssnes,&ksp);CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp,J,Jp);CHKERRQ(ierr);
    
    //for (i=0;i<dyn->ndiff;i++) {
    //  ierr = KSPSolve(ksp,dyn->s[idx[i]],dyn->s[idx[i]]);CHKERRQ(ierr);
    //}
    ierr = DYNComputeForwardJacobianP(dyn->ts,t,X,dyn->fwdJacp,(void*)dyn);CHKERRQ(ierr);
    
    for (i=0;i<dyn->nparams;i++) {
        ierr = VecGetSubVector(dyn->sp[i],dyn->isdiff,&subv);CHKERRQ(ierr);
        ierr = VecGetArrayRead(subv,&v);CHKERRQ(ierr);
        ierr = VecSetValues(dyn->fwdJacp[i],dyn->ndiff,idx,v,INSERT_VALUES);CHKERRQ(ierr);
        
        ierr = VecAssemblyBegin(dyn->fwdJacp[i]);CHKERRQ(ierr);
        ierr = VecAssemblyEnd(dyn->fwdJacp[i]);CHKERRQ(ierr);
        
        ierr = VecRestoreArrayRead(subv,&v);CHKERRQ(ierr);
        ierr = VecRestoreSubVector(dyn->sp[i],dyn->isdiff,&subv);CHKERRQ(ierr);
        
        ierr = KSPSolve(ksp,dyn->fwdJacp[i],dyn->sp[i]);CHKERRQ(ierr);
    }
    
    PetscFunctionReturn(0);
}

/*
  DYNForwardSetUp - Sets up a DYN object

  Input Parameters:
. DYN - the DYN object

  Notes:
  This routine sets up the DYN object and the underlying PS object. It
  also distributes the PS object when used in parallel.
*/
PetscErrorCode DYNForwardSetUp(DYN dyn)
{
  PetscErrorCode ierr;
  PS             ps=dyn->ps;
  PetscInt       ngen=0,nghostbuses=0;
  PetscInt       xdynsize=0;

  PetscFunctionBegin;

  /* Set up PS object */
  ierr = PSSetUp(ps);CHKERRQ(ierr);

  /* Set up differential equation IS */
  dyn->ndiff = ps->ndiff;
  ierr = DYNSetUpDifferentialEqIS(dyn,&dyn->isdiff);CHKERRQ(ierr);

  /* Create the solution vector and the Jacobian matrix */
  ierr = PSCreateGlobalVector(dyn->ps,&dyn->X);CHKERRQ(ierr);
  ierr = PSCreateMatrix(dyn->ps,&dyn->Jac);CHKERRQ(ierr);

  /* Associate the DM object in PS with the time-stepping solver */
  ierr = TSSetDM(dyn->ts,dyn->ps->networkdm);CHKERRQ(ierr);
  /* Default solver - implicit trapezoidal (can be changed via options) */
  ierr = TSSetType(dyn->ts,TSBEULER);CHKERRQ(ierr);

  /* Set the prefix for this TS.. All TS runtime options will need to have the prefix "-dyn_" */
  ierr = TSSetOptionsPrefix(dyn->ts,"dyn_");CHKERRQ(ierr);

  /* Set the function and Jacobian routines */
  ierr = TSSetIFunction(dyn->ts,NULL,DYNIFunction,(void*)dyn);CHKERRQ(ierr);
  ierr = TSSetIJacobian(dyn->ts,dyn->Jac,dyn->Jac,(TSIJacobian)DYNIJacobian,(void*)dyn);CHKERRQ(ierr);

  /* Stop when the final time is reached */
  ierr = TSSetExactFinalTime(dyn->ts,TS_EXACTFINALTIME_MATCHSTEP);CHKERRQ(ierr);

  /* Set up parameter sensitivity calculation */
  ierr = VecGetSize(dyn->X,&xdynsize);CHKERRQ(ierr);

  ierr = DYNGetCostFunctionParametersFromOptions(dyn);CHKERRQ(ierr);

  ierr = PSGetNumGenerators(dyn->ps,NULL,Ngen);CHKERRQ(ierr);
  for (i=0;i<dyn->ps->nbus;i++) {
    ierr = PSBUSIsGhosted(&dyn->ps->bus[i],&isghosted);CHKERRQ(ierr);
    if(isghosted) nghostbuses++;
  }
  dyn->nparams = 2*(dyn->ps->nbus-nghostbuses) + 2*Ngen;
  dyn->ncostfcns = Ngen;

  /* Set up vectors for sensitivities */
  ierr = VecCreate(dyn->comm->type,&dyn->Xparam);CHKERRQ(ierr);
  ierr = VecSetSizes(dyn->Xparam,dyn->nparams,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(dyn->Xparam);CHKERRQ(ierr);

  ierr = VecDuplicateVecs(dyn->X,dyn->nparams,&dyn->sp);CHKERRQ(ierr); /* forward sensitivities of final solution wrt parameters */
  /* Comment out the following line to disable computing sensitivity to initial values */
  // ierr = VecDuplicateVecs(dyn->X,xdynsize,&dyn->s);CHKERRQ(ierr); /* forward sensitivities of final solution wrt initial values */
  ierr = VecDuplicateVecs(dyn->Xparam,dyn->ncostfcns,&dyn->intgradp);CHKERRQ(ierr); /* gradients (forward sensitivities) of integrals wrt parameters */
  /* Comment out the following line to disable computing sensitivity to initial values */
  // ierr = VecDuplicateVecs(dyn->X,dyn->ncostfcns,&dyn->intgrad);CHKERRQ(ierr); /* gradients (forward sensitivities) of integrals wrt initial values */
  ierr = VecDuplicateVecs(dyn->X,dyn->Nparams,&dyn->dy0dp);CHKERRQ(ierr);

  ierr = VecDuplicate(dyn->X,&dyn->X0);CHKERRQ(ierr); /* temporary work vector */
  ierr = VecSet(dyn->X0,0.0);CHKERRQ(ierr);
  ierr = VecDuplicate(dyn->X,&dyn->Xtmp);CHKERRQ(ierr); /* temporary work vector */
  ierr = VecDuplicate(dyn->X,&dyn->Xtmp2);CHKERRQ(ierr); /* temporary work vector */
  dyn->s = NULL;
  dyn->intgrad = NULL;
  /* Set trajectory sensitivities for dyn->ts */
  ierr = TSSetTrajectorySensitivities(dyn->ts,dyn->nparams,dyn->sp,xdynsize,dyn->s);CHKERRQ(ierr);

  /* Create the Jacobian of the DAE equations w.r.t. parameters */
  ierr = VecDuplicateVecs(dyn->X,dyn->nparams,&dyn->fwdJacp);CHKERRQ(ierr);

  ierr = TSForwardSetRHSJacobianP(dyn->ts,dyn->fwdJacp,DYNComputeForwardJacobianP,(void*)dyn);CHKERRQ(ierr);
  ierr = TSSetForwardIntegralGradients(dyn->ts,dyn->ncostfcns,dyn->intgradp,dyn->intgrad);CHKERRQ(ierr);
  ierr = TSSetCostIntegrand(dyn->ts,dyn->ncostfcns,(PetscErrorCode (*)(TS,PetscReal,Vec,Vec,void*))DYNCostIntegrand,
                 (PetscErrorCode (*)(TS,PetscReal,Vec,Vec*,void*))DYNDRDYFunction,
                 (PetscErrorCode (*)(TS,PetscReal,Vec,Vec*,void*))DYNDRDPFunction,PETSC_TRUE,(void*)dyn);CHKERRQ(ierr);

  /* Set up events */
  ierr = DYNSetUpEvents(dyn);CHKERRQ(ierr);

  /* Set up the event monitoring function for TS */
  ierr = TSSetEventHandler(dyn->ts,dyn->toteventmonitors,dyn->eventdirection,dyn->eventterminate,DYNEventMonitor,DYNEventPostFunction,(void*)dyn);CHKERRQ(ierr);

  if(dyn->monitor) {
    /* Set up user monitor print out information for visualization */
    ierr = TSMonitorSet(dyn->ts,DYNUserMonitor,(void *)dyn,NULL);CHKERRQ(ierr);
  }

  ierr = TSSetFromOptions(dyn->ts);CHKERRQ(ierr);

  if(dyn->prepflow) {
    ierr = DYNSetUpInitPflow(dyn);CHKERRQ(ierr);
  }

  dyn->setupcalled = PETSC_TRUE;
  PetscFunctionReturn(0);
}

/*
  DYNForwardSolve - Perform forward sensitivity analysis of the dynamic simulation

  Input Parameters:
. dyn - the dynamics simulation object

  Notes:
   The parameters considered are the bus voltage magnitudes and angles, and the generator real and reactive power
   output. For each bus, its parameters are ordered as [Va, Vm, PG1, QG1, PG2, QG2,.., PGn, QGn]
*/
PetscErrorCode DYNForwardSolve(DYN dyn)
{
  PetscErrorCode ierr;
  PetscInt       i;
  const PetscInt *idx;
  PetscBool      pflowconverged;

  PetscFunctionBegin;
  if(!dyn->setupcalled) {
    ierr = DYNForwardSetUp(dyn);CHKERRQ(ierr);
  }

  /* Initial steady-state power flow */
  if(dyn->prepflow) {
    ierr = DYNComputePrePflow(dyn,&pflowconverged);CHKERRQ(ierr);
    if (!pflowconverged) SETERRQ(PETSC_COMM_SELF,0,"Pflow did not converge");
  }

  /* Set up initial conditions */
  ierr = DYNSetInitialConditions(dyn,dyn->X);CHKERRQ(ierr);

  /* Set solution vector with TS */
  ierr = TSSetSolution(dyn->ts,dyn->X);CHKERRQ(ierr);

  SNES     tssnes;
  PetscInt xdynsize=0;
  Mat      tssnesjac;
  void     *tssnesjacctx;
  //void     *tssnesctx;
  //Vec      tssnesres;
  //PetscErrorCode (*tssnesfunc)(SNES,Vec,Vec,void*);
  PetscErrorCode (*tssnesjacfun)(SNES,Vec,Mat,Mat,void*);

  /* Comment out the following line to disable computing sensitivity to initial values */
  // ierr = VecGetSize(dyn->X,&xdynsize);CHKERRQ(ierr);

  /* Set solve_alg_only flag and change the function and Jacobian evaluation routine pointers */
  dyn->solve_alg_only = PETSC_TRUE;
  ierr = TSGetSNES(dyn->ts,&tssnes);CHKERRQ(ierr);
  ierr = SNESGetJacobian(tssnes,&tssnesjac,NULL,&tssnesjacfun,&tssnesjacctx);CHKERRQ(ierr);
  //ierr = SNESGetFunction(tssnes,&tssnesres,&tssnesfunc,&tssnesctx);CHKERRQ(ierr);
  //ierr = SNESSetFunction(tssnes,tssnesres,DYNEventPostSolveAlgebraicSensi,(void*)dyn);CHKERRQ(ierr);
  ierr = SNESSetJacobian(tssnes,dyn->Jac,dyn->Jac,DYNEventPostSolveAlgebraicJacobian,(void*)dyn);CHKERRQ(ierr);
  /* Initialize parameter sensitivity work vectors */
  for(i=0; i < dyn->ncostfcns; i++) {
    ierr = VecSet(dyn->intgradp[i],0.0);CHKERRQ(ierr);
    /* Comment the following line to disable sensitivity to initial values */
    // ierr = VecSet(dyn->intgrad[i],0.0);CHKERRQ(ierr);
  }
  ierr = DYNComputeSensiP(dyn);CHKERRQ(ierr);
  for(i=0; i<xdynsize; i++) {
    ierr = VecSet(dyn->s[i],0.0);CHKERRQ(ierr);
  }
  ierr = ISGetIndices(dyn->isdiff,&idx);CHKERRQ(ierr);

  // for(i=0; i<dyn->ndiff; i++) {
  //  ierr = VecSetValue(dyn->s[idx[i]],idx[i],1.0,INSERT_VALUES);CHKERRQ(ierr);
  // }

  ierr = DYNForwardAlgebraicSolve(tssnes,dyn,dyn->t0,dyn->X);

  /* Reset flag and function, jacobian evaluation routine pointers */
  dyn->solve_alg_only = PETSC_FALSE;
  //ierr = SNESSetFunction(tssnes,tssnesres,tssnesfunc,tssnesctx);CHKERRQ(ierr);
  ierr = SNESSetJacobian(tssnes,tssnesjac,tssnesjac,tssnesjacfun,tssnesjacctx);CHKERRQ(ierr);

  /*
  Vec C;
  ierr = TSGetCostIntegral(dyn->ts,&C);CHKERRQ(ierr);
  ierr = VecView(C,0);CHKERRQ(ierr);
  */

  /* Call the time-stepping solver */
  ierr = TSSolve(dyn->ts,dyn->X);CHKERRQ(ierr);

  /* Get converged reason */
  ierr = TSGetConvergedReason(dyn->ts,&reason);CHKERRQ(ierr);
  if (reason<0) {
    dyn->tsconverged = PETSC_FALSE;
  } else dyn->tsconverged = PETSC_TRUE;

  /* Get the actual finish time */
  ierr = TSGetSolveTime(dyn->ts,dyn->finaltime);CHKERRQ(ierr);

  /*
  Vec C;
  ierr = TSGetCostIntegral(dyn->ts,&C);CHKERRQ(ierr);
  ierr = VecView(C,0);CHKERRQ(ierr);

  for(i=0;i<dyn->nparams;i++) {
    ierr = VecView(dyn->s[i],0);CHKERRQ(ierr);
  }

  for(i=0;i<dyn->ncostfcns;i++) {
    ierr = VecView(dyn->intgradp[i],0);CHKERRQ(ierr);
  }
  */
  if(dyn->printsensi){
    PetscInt          j;
    const PetscScalar *s;
    FILE              *f;
    f = fopen("fwdsensi.data", "a");
    ierr = PetscFPrintf(PETSC_COMM_WORLD,f,"%24.10f ",dyn->finaltime);CHKERRQ(ierr);
    for(i=0;i<dyn->ncostfcns;i++) {
      ierr = VecGetArrayRead(dyn->intgradp[i],&s);CHKERRQ(ierr);
      for(j=0;j<dyn->nparams;j++) {
        ierr = PetscFPrintf(PETSC_COMM_WORLD,f,"%24.15e ",s[j]);CHKERRQ(ierr);
      }
      ierr = VecRestoreArrayRead(dyn->intgradp[i],&s);CHKERRQ(ierr);
    }
    ierr = PetscFPrintf(PETSC_COMM_WORLD,f,"\n");CHKERRQ(ierr);
    fclose(f);
  }
  /*
  for(i=0;i<dyn->ncostfcns;i++) {
    ierr = VecView(dyn->intgrad[i],0);CHKERRQ(ierr);
    //ierr = CheckAlgebraicPart(dyn,dyn->intgrad[i]);CHKERRQ(ierr);
  }
  PetscViewer viewer;
  Vec         *DCDP,*DCDI;
  PetscReal   distance = 0,norm1,norm2;
  char        filename[PETSC_MAX_PATH_LEN] = "finitediff_results.bin";
  ierr = VecDuplicateVecs(dyn->intgradp[0],dyn->ncostfcns,&DCDP);CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
  for(i=0;i<dyn->ncostfcns;i++) {
    ierr = VecLoad(DCDP[i],viewer);CHKERRQ(ierr);
    ierr = VecNorm(DCDP[i],NORM_2,&norm2);CHKERRQ(ierr);
    ierr = VecAXPY(DCDP[i],-1,dyn->intgradp[i]);CHKERRQ(ierr);
    ierr = VecNorm(DCDP[i],NORM_2,&norm1);CHKERRQ(ierr);
    if (norm2>0) distance += PetscSqrtScalar(norm1/norm2);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"norm1=%lf norm2=%lf\n",norm1,norm2);CHKERRQ(ierr);
  }
  ierr = VecDestroyVecs(dyn->ncostfcns,&DCDP);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"distance to finitediff results: %lf\n",distance);CHKERRQ(ierr);

  ierr = VecDuplicateVecs(dyn->intgrad[0],dyn->ncostfcns,&DCDI);CHKERRQ(ierr);
  for(i=0;i<dyn->ncostfcns;i++) {
    ierr = VecLoad(DCDI[i],viewer);CHKERRQ(ierr);
    ierr = VecNorm(DCDI[i],NORM_2,&norm2);CHKERRQ(ierr);
    ierr = VecAXPY(DCDI[i],-1,dyn->intgrad[i]);CHKERRQ(ierr);
    ierr = VecNorm(DCDI[i],NORM_2,&norm1);CHKERRQ(ierr);
    if (norm2>0) distance += PetscSqrtScalar(norm1/norm2);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"norm1=%lf norm2=%lf\n",norm1,norm2);CHKERRQ(ierr);
  }
  ierr = VecDestroyVecs(dyn->ncostfcns,&DCDI);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  */
  PetscFunctionReturn(0);
}

/*
 DYNForwardDestroy - Destroys the DYN object
 
 Input Parameter
 . dyn - The DYN object to destroy
 */
PetscErrorCode DYNForwardDestroy(DYN *dyn)
{
    PetscErrorCode ierr;
    PetscInt       xdynsize=0;
    
    PetscFunctionBegin;
    /* Comment out the following line to disable computing sensitivity to initial values */
    // ierr = VecGetSize((*dyn)->X,&xdynsize);CHKERRQ(ierr);
    ierr = VecDestroyVecs((*dyn)->ncostfcns,&(*dyn)->intgradp);CHKERRQ(ierr);
    /* Comment out the following line to disable computing sensitivity to initial values */
    // ierr = VecDestroyVecs((*dyn)->ncostfcns,&(*dyn)->intgrad);CHKERRQ(ierr);
    ierr = VecDestroyVecs((*dyn)->Nparams,&(*dyn)->dy0dp);CHKERRQ(ierr);
    ierr = VecDestroyVecs((*dyn)->Nparams,&(*dyn)->fwdJacp);CHKERRQ(ierr);
    ierr = VecDestroyVecs(xdynsize,&(*dyn)->s);CHKERRQ(ierr);
    ierr = VecDestroyVecs((*dyn)->nparams,&(*dyn)->sp);CHKERRQ(ierr);
    ierr = VecDestroy(&(*dyn)->Xparam);CHKERRQ(ierr);
    ierr = VecDestroy(&(*dyn)->X0);CHKERRQ(ierr);
    ierr = VecDestroy(&(*dyn)->Xtmp);CHKERRQ(ierr);
    ierr = VecDestroy(&(*dyn)->Xtmp2);CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

#endif

PetscErrorCode DYNReInit(DYN dyn)
{
  Vec            C;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = TSMonitorCancel(dyn->ts);CHKERRQ(ierr);
  //if(dyn->prepflow) {
  //  ierr = DYNComputePrePflow(dyn,pflowconverged);CHKERRQ(ierr);
  //}
  ierr = TSGetCostIntegral(dyn->ts,&C);CHKERRQ(ierr);
  ierr = VecSet(C,0);CHKERRQ(ierr);
  /* Set up initial conditions */
  ierr = DYNSetInitialConditions(dyn,dyn->X);CHKERRQ(ierr);
  /* Call the time-stepping solver */
  ierr = DYNSetDuration(dyn,dyn->maxsteps,dyn->tmax);CHKERRQ(ierr);
  ierr = DYNSetStartTimeAndStep(dyn,dyn->t0,dyn->dt0);CHKERRQ(ierr);
  ierr = TSSetFromOptions(dyn->ts);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DYNFDApprox(DYN dyn,Vec *gradients,const PetscScalar eps,PetscInt index)
{
  Vec               C;
  PetscInt          i;
  PetscScalar       *g;
  const PetscScalar *c;
  PetscErrorCode    ierr;

  PetscFunctionBegin;
  ierr = TSGetCostIntegral(dyn->ts,&C);CHKERRQ(ierr);
  ierr = VecGetArrayRead(C,&c);CHKERRQ(ierr);
  for(i=0;i<dyn->ncostfcns;i++) {
    ierr = VecGetArray(gradients[i],&g);CHKERRQ(ierr);
    g[index] = (c[i]-g[index])/eps;
    ierr = VecRestoreArray(gradients[i],&g);CHKERRQ(ierr);
  }
  ierr = VecRestoreArrayRead(C,&c);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNFiniteDiffSolve - Perform finite difference of the dynamic simulation

  Input Parameters:
. dyn - the dynamics simulation object

  Notes:
   The parameters considered are the bus voltage magnitudes and angles, and the generator real and reactive power
   output. For each bus, its parameters are ordered as [Va, Vm, PG1, QG1, PG2, QG2,.., PGn, QGn]
*/
PetscErrorCode DYNFiniteDiffSolve(DYN dyn)
{
  Vec               X0,C,P,*DCDP,*DCDI;
  const PetscScalar *c,eps=1e-7;
  const PetscInt    *idx;
  PetscScalar       *dcdp;
  Vec               subv;
  PetscBool         ghostbus;
  PS                ps=dyn->ps;
  PSBUS             bus;
  PSGEN             gen;
  PetscInt          i,j,k,ctr,ngen=ps->ngen,nbus=ps->nbus,index=0;
  PetscBool         pflowconverged;
  PetscErrorCode    ierr;

  PetscFunctionBegin;
  dyn->eval_integral = PETSC_TRUE;

  /* base for finite difference */
  if(!dyn->setupcalled) {
    ierr = DYNSetUp(dyn);CHKERRQ(ierr);
  }
  /* Initial steady-state power flow */
  if(dyn->prepflow) {
    ierr = DYNComputePrePflow(dyn,&pflowconverged);CHKERRQ(ierr);
    if (!pflowconverged) SETERRQ(PETSC_COMM_SELF,0,"Pflow did not converge");
  }
  ierr = DYNSetInitialConditions(dyn,dyn->X);CHKERRQ(ierr);
  /* Original IC saved to X0 */
  ierr = VecDuplicate(dyn->X,&X0);CHKERRQ(ierr);
  ierr = VecCopy(dyn->X,X0);CHKERRQ(ierr);
  ierr = TSSolve(dyn->ts,dyn->X);CHKERRQ(ierr);

  /* Set up vectors to hold gradients */
  //ierr = PetscCalloc1(dyn->ndiff,&idx);CHKERRQ(ierr);
  ierr = ISGetIndices(dyn->isdiff,&idx);CHKERRQ(ierr);
  ierr = ISView(dyn->isdiff,0);CHKERRQ(ierr);
  ierr = TSGetCostIntegral(dyn->ts,&C);CHKERRQ(ierr);
  ierr = VecGetArrayRead(C,&c);CHKERRQ(ierr);
  dyn->nparams = 2*nbus + 2*ngen;
  ierr = VecCreateSeq(PETSC_COMM_SELF,dyn->nparams,&P);CHKERRQ(ierr);
  ierr = VecDuplicateVecs(P,dyn->ncostfcns,&DCDP);CHKERRQ(ierr);
  ierr = VecDuplicateVecs(dyn->X,dyn->ncostfcns,&DCDI);CHKERRQ(ierr);
  ierr = VecDestroy(&P);CHKERRQ(ierr);
  for(i=0;i<dyn->ncostfcns;i++) {
    ierr = VecSet(DCDP[i],c[i]);CHKERRQ(ierr);
    ierr = VecSet(DCDI[i],0);CHKERRQ(ierr);
    ierr = VecGetSubVector(DCDI[i],dyn->isdiff,&subv);CHKERRQ(ierr);
    ierr = VecSet(subv,c[i]);CHKERRQ(ierr);
    ierr = VecRestoreSubVector(DCDI[i],dyn->isdiff, &subv);
  }
  ierr = VecRestoreArrayRead(C,&c);CHKERRQ(ierr);

  for(i=0;i<nbus;i++) {
    bus = &ps->bus[i];
    ierr = PSBUSIsGhosted(bus,&ghostbus);CHKERRQ(ierr);
    if(ghostbus) continue;
    /* Perturb Va and compute sensitivity */
    bus->va += eps;
    ierr = DYNReInit(dyn);CHKERRQ(ierr);
    ierr = TSSolve(dyn->ts,dyn->X);CHKERRQ(ierr);
    bus->va -= eps;
    ierr = DYNFDApprox(dyn,DCDP,eps,index);CHKERRQ(ierr);

    /* Perturb Vm and compute sensitivity */
    if(bus->ngen) {
      for(k=0;k<bus->ngen;k++) {
        ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
        if(gen->status) gen->vs += eps;
      }
    }
    bus->vm += eps;
    ierr = DYNReInit(dyn);CHKERRQ(ierr);
    ierr = TSSolve(dyn->ts,dyn->X);CHKERRQ(ierr);

    if(bus->ngen) {
      ierr = PSBUSGetGen(bus,0,&gen);CHKERRQ(ierr);
      for(k=0; k<bus->ngen; k++) {
        ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
        if(gen->status) gen->vs -= eps;
      }
    }
    bus->vm -= eps;
    ierr = DYNFDApprox(dyn,DCDP,eps,index+1);CHKERRQ(ierr);

    index += 2;
    ctr = 0;

    for(k=0;k<bus->ngen;k++) {
      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      if(!gen->status) {
        for(j=0;j<dyn->ncostfcns;j++) {
          ierr = VecGetArray(DCDP[j],&dcdp);CHKERRQ(ierr);
          dcdp[index]   = 0;
          dcdp[index+1] = 0;
          ierr = VecRestoreArray(DCDP[j],&dcdp);CHKERRQ(ierr);
        }
        ctr += 2;
        continue;
      }
      /* Perturb PG and compute sensitivity */
      gen->pg += eps;
      ierr = DYNReInit(dyn);CHKERRQ(ierr);
      ierr = TSSolve(dyn->ts,dyn->X);CHKERRQ(ierr);
      gen->pg -= eps;
      ierr = DYNFDApprox(dyn,DCDP,eps,index+ctr);CHKERRQ(ierr);
#ifdef PINDEX
      printf("index=%d\n",index+ctr);
#endif
      /* Perturb QG and compute sensitivity */
      gen->qg += eps;
      ierr = DYNReInit(dyn);CHKERRQ(ierr);
      ierr = TSSolve(dyn->ts,dyn->X);CHKERRQ(ierr);
      gen->qg -= eps;
      ierr = DYNFDApprox(dyn,DCDP,eps,index+ctr+1);CHKERRQ(ierr);
      ctr += 2;
    }
    index += ctr;
  }

  if(dyn->printsensi){
    PetscInt          j;
    PetscReal         finaltime;
    const PetscScalar *s;
    FILE              *f;
    f = fopen("finitediffsensi.data", "a");
    ierr = TSGetSolveTime(dyn->ts,&finaltime);CHKERRQ(ierr);
    ierr = PetscFPrintf(PETSC_COMM_WORLD,f,"%24.10f ",finaltime);CHKERRQ(ierr);
    for(i=0;i<dyn->ncostfcns;i++) {
      ierr = VecGetArrayRead(DCDP[i],&s);CHKERRQ(ierr);
      for(j=0;j<dyn->nparams;j++) {
        ierr = PetscFPrintf(PETSC_COMM_WORLD,f,"%24.15e ",s[j]);CHKERRQ(ierr);
      }
      ierr = VecRestoreArrayRead(DCDP[i],&s);CHKERRQ(ierr);
    }
    ierr = PetscFPrintf(PETSC_COMM_WORLD,f,"\n");CHKERRQ(ierr);
  }

  PetscViewer viewer;
  char        filename[PETSC_MAX_PATH_LEN] = "finitediff_results.bin";

  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
  for(i=0;i<dyn->ncostfcns;i++) {
    ierr = VecView(DCDP[i],viewer);CHKERRQ(ierr);
    ierr = VecView(DCDP[i],0);CHKERRQ(ierr);
  }

  ierr = VecDestroyVecs(dyn->ncostfcns,&DCDP);CHKERRQ(ierr);

  PetscInt       xdynsize;
  PetscScalar    *x;
  SNES           tssnes;
  void           *tssnesctx;
  Vec            tssnesres;
  Mat            tssnesjac;
  void           *tssnesjacctx;
  PetscErrorCode (*tssnesfunc)(SNES,Vec,Vec,void*);
  PetscErrorCode (*tssnesjacfun)(SNES,Vec,Mat,Mat,void*);

  ierr = VecGetSize(dyn->X,&xdynsize);CHKERRQ(ierr);
  for(i=0;i<dyn->ndiff;i++) {
    ierr = DYNReInit(dyn);CHKERRQ(ierr);
    ierr = VecGetArray(dyn->X,&x);CHKERRQ(ierr);
    x[idx[i]] += eps;
    ierr = VecRestoreArray(dyn->X,&x);CHKERRQ(ierr);
    /* Set solve_alg_only flag and change the function and Jacobian evaluation routine pointers */
    dyn->solve_alg_only = PETSC_TRUE;
    ierr = TSGetSNES(dyn->ts,&tssnes);CHKERRQ(ierr);
    ierr = SNESGetFunction(tssnes,&tssnesres,&tssnesfunc,&tssnesctx);CHKERRQ(ierr);
    ierr = SNESGetJacobian(tssnes,&tssnesjac,NULL,&tssnesjacfun,&tssnesjacctx);CHKERRQ(ierr);
    ierr = SNESSetJacobian(tssnes,dyn->Jac,dyn->Jac,DYNEventPostSolveAlgebraicJacobian,(void*)dyn);CHKERRQ(ierr);
    ierr = TSGetSNES(dyn->ts,&tssnes);CHKERRQ(ierr);
    ierr = SNESSetFunction(tssnes,tssnesres,DYNEventPostSolveAlgebraic,(void*)dyn);CHKERRQ(ierr);
    ierr = SNESSolve(tssnes,NULL,dyn->X);CHKERRQ(ierr);
    /* Reset flag and function, jacobian evaluation routine pointers */
    dyn->solve_alg_only = PETSC_FALSE;
    ierr = SNESSetFunction(tssnes,tssnesres,tssnesfunc,tssnesctx);CHKERRQ(ierr);
    ierr = SNESSetJacobian(tssnes,tssnesjac,tssnesjac,tssnesjacfun,tssnesjacctx);CHKERRQ(ierr);
    //ierr = PetscPrintf(PETSC_COMM_WORLD,"***difference***:\n");
    //ierr = VecAXPY(X0,-1,dyn->X);CHKERRQ(ierr);
    //ierr = VecView(X0,0);CHKERRQ(ierr);
    /* Run again */
    ierr = TSSolve(dyn->ts,dyn->X);CHKERRQ(ierr);
    /* Finite differencing */
    ierr = DYNFDApprox(dyn,DCDI,eps,idx[i]);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(C,&c);CHKERRQ(ierr);
  }
  /*
  for(i=0;i<dyn->ncostfcns;i++) {
    ierr = VecView(DCDI[i],viewer);CHKERRQ(ierr);
    ierr = VecView(DCDI[i],0);CHKERRQ(ierr);
  }
  */
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  ierr = VecDestroyVecs(dyn->ncostfcns,&DCDI);CHKERRQ(ierr);
  ierr = VecDestroy(&P);CHKERRQ(ierr);
  ierr = VecDestroy(&X0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNGetCostFunctionAndSensitivites - This function runs the forward solve and returns the
  cost functions along with their sensitivities

  Input Parameters:
+ netfilename       - Name of the network file
. dyrfilename       - Name of the dynamic data dyr file
. eventfile         - Name of the event file
. costfcn           - A double array having the cost functions (must be allocated)
. sensitivities     - An array holding the sensitivities (must be allocated)
. single_costfcn    - Flag to indicate whether a single
. active_power_only - Flag to indicate that sensitivities w.r.t. only active power
                      only needs to be computed.
. updatedispatch    - Flag to update dispatch
. genbus            - An array of generator bus numbers
. gennum            - An array of generator numbers (starting from 0)
. gen_status        - An array of generator commitments
. Pg                - An array of active power dispatch
. Qg                - An array of reactvie power dispatch
. endtime           - End time of the simulation
. stepsize          - Time step size
. print_vm          - Flag to print voltage
. cost_type         - Indicator of cost function type (frequency violation or frequency)
. freq_min          - Minimum frequency threshold
. freq_max          - Maximum frequency threshold
. monitor           - Indicator of whether to use a monitor
- finaltime         - Final time at which the simulation actually ends

  Return values:
  + 0        - success
  . 1        - power flow does not converge
  - 2        - error during midlle of time integration (e.g. SNES nonconvergence)
*/
PetscInt DYNGetCostFunctionAndSensitivities(
    int argc,char **argv,const char netfile[],const char dyrfile[],const char eventfile[],
    double *costfcnout,double *sensitivitiesout,
    PetscBool single_costfcn,PetscBool active_power_only,PetscBool updatedispatch,
    PetscInt *genbus,PetscInt* gennum,PetscInt *gen_status,PetscScalar *Pg,PetscScalar *Qg,
    PetscReal endtime,PetscReal stepsize,PetscInt print_vm,PetscInt cost_type,PetscReal freq_min,PetscReal freq_max,PetscBool monitor,PetscReal *finaltime)
{
  DYN            dyn;
  PetscInt       i,flag;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  /* Create DYN object */
  ierr = DYNCreate(PETSC_COMM_WORLD,&dyn);CHKERRQ(ierr);

  /* Do not print sensitivity results to file */
  dyn->printsensi = PETSC_FALSE;
  dyn->monitor = monitor;

  /* Set up frequency ranges */
  dyn->freq_min=freq_min;  //59.9
  dyn->freq_max=freq_max; //60.1

  /* Get network data file from command line */
  ierr = DYNReadMatPowerData(dyn,netfile);CHKERRQ(ierr);

  /* Read dyr file */
  ierr = DYNReadDyrData(dyn,dyrfile);CHKERRQ(ierr);

  /* Read event data file */
  ierr = DYNReadEventData(dyn,eventfile);CHKERRQ(ierr);

  /* Set cost function type  */
  if (cost_type == 0) dyn->objfun = FREQ;
  else {
    dyn->objfun = FREQVIOL;
    dyn->sum_freq_costfun = single_costfcn; // only makes sense for frequency voilations
  }

  if(updatedispatch) {
    /* Update generator commitmentments, dispatch */
    //printf("ngen=%d\n",dyn->ps->ngen);
    for(i=0; i < dyn->ps->ngen-1; i++) {
      //printf("%d %d %d %f %f\n",genbus[i],gennum[i],gen_status[i],Pg[i],Qg[i]);
      ierr = DYNSetGenDispatchandStatus(dyn,genbus[i],gennum[i],gen_status[i],Pg[i],Qg[i]);CHKERRQ(ierr);
    }
  }

  /* Set time step */
  ierr = DYNSetStartTimeAndStep(dyn,0.0,stepsize);CHKERRQ(ierr); //.001

  /* Set Duration */
  ierr = DYNSetDuration(dyn,10000,endtime);CHKERRQ(ierr);  //0.1+epsilon

  /* Flag to set computing initial power flow */
  ierr = DYNSetComputePrePflow(dyn,PETSC_TRUE);CHKERRQ(ierr);

  ierr = DYNAdjointSetUp(dyn);CHKERRQ(ierr);
  //if(set_load == 1) {ierr = setload(demand, dyn);CHKERRQ(ierr);}

  /* Solve */
  ierr = DYNAdjointSolve(dyn);CHKERRQ(ierr);

  /* Copy cost function values */
  ierr = DYNGetCostFunction(dyn,dyn->finaltime,costfcnout,dyn->ncostfcns);CHKERRQ(ierr);

  if(!active_power_only) {
    for (i=0; i<dyn->ncostfcns; i++) {
      Vec          sensi_seq;
      VecScatter   ctx;
      PetscScalar *x;

      ierr = VecScatterCreateToAll(dyn->mu[i],&ctx,&sensi_seq);CHKERRQ(ierr);
      VecScatterBegin(ctx,dyn->mu[i],sensi_seq,INSERT_VALUES,SCATTER_FORWARD);
      VecScatterEnd(ctx,dyn->mu[i],sensi_seq,INSERT_VALUES,SCATTER_FORWARD);
      VecScatterDestroy(&ctx);
      ierr = VecGetArray(sensi_seq,&x);CHKERRQ(ierr);
      ierr = PetscMemcpy(sensitivitiesout+dyn->Nparams*i,x,dyn->Nparams*sizeof(double));CHKERRQ(ierr);
      ierr = VecRestoreArray(sensi_seq,&x);CHKERRQ(ierr);
      VecDestroy(&sensi_seq);
    }
  } else {
    PetscInt   ngen,Ngen;
    IS         ispg;
    VecScatter ctx;
    Vec        sensi_seq;

    ierr = DYNSetUpActivePowerIS(dyn,&ispg);CHKERRQ(ierr);
    for (i=0; i<dyn->ncostfcns; i++) {
      ierr = PSGetNumGenerators(dyn->ps,&ngen,&Ngen);CHKERRQ(ierr);
      ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,Ngen,sensitivitiesout+Ngen*i,&sensi_seq);CHKERRQ(ierr);
      ierr = VecScatterCreate(dyn->mu[i],ispg,sensi_seq,NULL,&ctx);CHKERRQ(ierr);
      VecScatterBegin(ctx,dyn->mu[i],sensi_seq,INSERT_VALUES,SCATTER_FORWARD);
      VecScatterEnd(ctx,dyn->mu[i],sensi_seq,INSERT_VALUES,SCATTER_FORWARD);
      VecScatterDestroy(&ctx);
      VecDestroy(&sensi_seq);
      //PetscInt idx,gloc;
      //PS       ps=dyn->ps;
      //PSBUS    bus;
      //idx = gloc = 0;
      //ierr = VecGetArray(dyn->mu[i],&x);CHKERRQ(ierr);
      //for(j=0; j < ps->nbus; j++) {
      //  bus = &ps->bus[j];
      //  gloc += 2;
      //  if(bus->ngen) {
      //    for(k=0; k < bus->ngen; k++) {
      //      sensitivitiesout[dyn->nparams*i+idx++] = x[gloc];
      //      gloc += 2;
      //    }
      //  }
      //}
      //ierr = VecRestoreArray(dyn->mu[i],&x);CHKERRQ(ierr);
    }
  }

  if (!dyn->prepflowconverged) flag = 1;
  else if (dyn->tsconverged) flag = 0;
  else { flag = 2; *finaltime = dyn->finaltime; }

  ierr = DYNAdjointDestroy(&dyn);CHKERRQ(ierr);
  ierr = DYNDestroy(&dyn);CHKERRQ(ierr);
  return flag;
}

PetscErrorCode Initialize(int argc,char **argv)
{
//    PetscInitializeNoPointers(argc,&argv,"petscopt",NULL);
    PetscInitialize(&argc,&argv,"petscopt",NULL);
    PetscFunctionReturn(0);
}

PetscErrorCode Finalize()
{
    PetscFinalize();
    PetscFunctionReturn(0);
}

PetscErrorCode DYNGetInitPflow(DYN dyn ,PFLOW *initpflow)
{

  PetscFunctionBegin;
  if(!dyn->setupcalled) SETERRQ(dyn->comm->type,0,"DYNSetUp() must be called before calling DYNGetInitPflow");
  *initpflow = dyn->initpflow;
  PetscFunctionReturn(0);
}

PetscErrorCode DYNSetUserData(DYN dyn, void* userdata)
{
  PetscFunctionBegin;
  dyn->userdata = userdata;
  PetscFunctionReturn(0);
}

PetscErrorCode DYNGetUserData(DYN dyn, void** userdata)
{
  PetscFunctionBegin;
  *userdata = dyn->userdata;
  PetscFunctionReturn(0);
}

PetscErrorCode DYNGetTime(DYN dyn,PetscReal *time)
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = TSGetTime(dyn->ts,time);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DYNGetStepNumber(DYN dyn,PetscInt *stepnum)
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = TSGetStepNumber(dyn->ts,stepnum);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DYNSaveSolution(DYN dyn,PetscReal time,Vec X)
{
  PetscErrorCode ierr;
  PetscInt       stepnum;
  TSTrajectory   traj;

  PetscFunctionBegin;
  ierr = TSGetTime(dyn->ts,&time);CHKERRQ(ierr);
  ierr = TSGetStepNumber(dyn->ts,&stepnum);CHKERRQ(ierr);
  ierr = TSGetTrajectory(dyn->ts,&traj);CHKERRQ(ierr);
  ierr = TSTrajectorySet(traj,dyn->ts,stepnum,time,X);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNGetBusVoltage - Returns the voltage (rectangular coordinates) for the given bus

  Input Parameters
+ dyn - The DYN object
. busnum  - bus number
. VD   - real part of bus voltage vector (pu)
. VQ   - imaginary bus voltage vector (pu)
- found - TRUE if bus found, else FALSE

Notes: DYNSolve() must be called before calling DYNGetBusVoltage()

*/
PetscErrorCode DYNGetBusVoltage(DYN dyn, PetscInt busnum, PetscScalar *VD, PetscScalar *VQ, PetscBool *found)
{
  PetscErrorCode ierr;
  PetscInt       intbusnum;
  PS             ps=dyn->ps;
  PSBUS          bus;
  Vec            X,localX;
  const PetscScalar *xarr;
  PetscInt       loc;

  PetscFunctionBegin;

  ierr = TSGetSolution(dyn->ts,&X);CHKERRQ(ierr);

  ierr = DMGetLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = VecGetArrayRead(localX,&xarr);CHKERRQ(ierr);


  /* Convert from external to internal bus number */
  intbusnum = ps->busext2intmap[busnum];
  if(intbusnum == -1) {
    *found = PETSC_FALSE;
    PetscFunctionReturn(0);
  } else {
    *found = PETSC_TRUE;
    bus = &ps->bus[intbusnum];
    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
    *VD = xarr[loc];
    *VQ = xarr[loc+1];
  }

  PetscFunctionReturn(0);
}


/*
  DYNGetBusGenCurrent - Returns the generation current injection (in rectangular coordinates) for the given bus

  Input Parameters
+ dyn - The DYN object
. busnum  - bus number
. IGD   - real part of generation current injection vector (pu)
. IGQ   - imaginary part of generation current injection vector (pu)
- found - TRUE if bus found, else FALSE

Notes: DYNSolve() must be called before calling DYNGetBusGenCurrent()

*/
PetscErrorCode DYNGetBusGenCurrent(DYN dyn, PetscInt busnum, PetscScalar *IGD, PetscScalar *IGQ, PetscBool *found)
{
  PetscErrorCode ierr;
  PetscInt       intbusnum;
  PS             ps=dyn->ps;
  PSBUS          bus;
  Vec            X,localX;
  const PetscScalar *xarr,*xdyn;
  PetscInt       loc;
  PetscScalar    VD,VQ;
  PetscInt       k;

  PetscFunctionBegin;

  *IGD = *IGQ = 0.0;

  ierr = TSGetSolution(dyn->ts,&X);CHKERRQ(ierr);

  ierr = DMGetLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = VecGetArrayRead(localX,&xarr);CHKERRQ(ierr);


  /* Convert from external to internal bus number */
  intbusnum = ps->busext2intmap[busnum];
  if(intbusnum == -1) {
    *found = PETSC_FALSE;
    PetscFunctionReturn(0);
  } else {
    *found = PETSC_TRUE;
    bus = &ps->bus[intbusnum];
    ierr = PSBUSGetVariableLocation(bus,&loc);CHKERRQ(ierr);
    VD = xarr[loc];
    VQ = xarr[loc+1];

    for(k=0; k < bus->ngen; k++) {
      PSGEN gen;
      DYNGenModel dyngen;
      PetscScalar ID,IQ;

      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr); 

      if(!gen->status) continue;

      ierr = PSGENGetDYNGen(gen,&dyngen);CHKERRQ(ierr);
      xdyn = xarr + loc;

      ierr = DYNGenModelGetCurrent(dyngen,VD,VQ,(PetscScalar*)xdyn,&ID,&IQ);CHKERRQ(ierr);
      *IGD += ID; *IGQ += IQ;
    }
  }

  PetscFunctionReturn(0);
}

/*
  DYNSetBusVoltage - Sets the voltage for the given bus

  Input Parameters
+ dyn - The DYN object
. busnum  - bus number
. VD   - real part of generation current injection vector (pu)
- VQ   - imaginary part of generation current injection vector (pu)

Notes: DYNSolve() must be called before calling DYNGetBusGenCurrent()
   This routine updates the voltage set points for the generator models,
   particularly the constant voltage source generator model.

*/
PetscErrorCode DYNSetBusVoltage(DYN dyn, PetscInt busnum, PetscScalar VD, PetscScalar VQ)
{
  PetscErrorCode ierr;
  PS             ps=dyn->ps;
  PetscInt       intbusnum;
  PSBUS          bus;
  PetscInt       k;

  PetscFunctionBegin;

  /* Convert from external to internal bus number */
  intbusnum = ps->busext2intmap[busnum];
  if(intbusnum == -1) {
    PetscFunctionReturn(0);
  } else {
    bus = &ps->bus[intbusnum];

    for(k=0; k < bus->ngen; k++) {
      PSGEN gen;
      DYNGenModel dyngen;

      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr); 

      if(!gen->status) continue;

      ierr = PSGENGetDYNGen(gen,&dyngen);CHKERRQ(ierr);

      ierr = DYNGenModelSetVoltage(dyngen,VD,VQ);CHKERRQ(ierr);
    }
  }

  PetscFunctionReturn(0);
}

PetscErrorCode DYNSetBusConstantPowerLoad(DYN dyn, PetscInt busnum, PetscScalar PD, PetscScalar QD)
{
  PetscErrorCode ierr;
  PS             ps=dyn->ps;
  PetscInt       intbusnum;
  PSBUS          bus;
  PetscInt       k;

  PetscFunctionBegin;

  /* Convert from external to internal bus number */
  intbusnum = ps->busext2intmap[busnum];
  if(intbusnum == -1) {
    PetscFunctionReturn(0);
  } else {
    bus = &ps->bus[intbusnum];

    for(k=0; k < bus->nload; k++) {
      PSLOAD load;
      DYNLoadModel dynload;

      ierr = PSBUSGetLoad(bus,k,&load);CHKERRQ(ierr); 

      if(!load->status) continue;

      ierr = PSLOADGetDYNLoad(load,&dynload);CHKERRQ(ierr);

      ierr = DYNLoadModelSetConstantPowerLoad(dynload,PD,QD);CHKERRQ(ierr);
    }
  }
  
  PetscFunctionReturn(0);
}

PetscErrorCode DYNSolveAlgebraicEquationsOnly(DYN dyn)
{
  PetscErrorCode ierr;
  SNES tssnes;
  void *tssnesctx;
  void *tssnesjacctx;
  PetscErrorCode (*tssnesfunc)(SNES,Vec,Vec,void*);
  PetscErrorCode (*tssnesjacfun)(SNES,Vec,Mat,Mat,void*);
  PetscInt tssneslag;
  SNESConvergedReason reason;
  TS ts=dyn->ts;
  PetscReal t;
  Vec       X;
  
  ierr = TSGetTime(ts,&t);CHKERRQ(ierr);
  ierr = TSGetSolution(ts,&X);CHKERRQ(ierr);

  /* Check topology changes */
  //    ierr = PSCheckTopology(dyn->ps);CHKERRQ(ierr);

  /* Set solve_alg_only flag and change the function and Jacobian evaluation routine pointers */
  dyn->solve_alg_only = PETSC_TRUE;

  /* Get Inner SNES */
  ierr = TSGetSNES(dyn->ts,&tssnes);CHKERRQ(ierr);
  
  /* Get the function and jacobian evaluation routines. SNES associates it with the DM attached to it */
  ierr = DMSNESGetFunction(dyn->ps->networkdm,&tssnesfunc,&tssnesctx);CHKERRQ(ierr);
  ierr = DMSNESGetJacobian(dyn->ps->networkdm,&tssnesjacfun,&tssnesjacctx);CHKERRQ(ierr);
  
  /* Set function and jacobian evaluation routines for the algebtraic solve */
  ierr = DMSNESSetFunction(dyn->ps->networkdm,DYNEventPostSolveAlgebraic,(void*)dyn);CHKERRQ(ierr);
  ierr = DMSNESSetJacobian(dyn->ps->networkdm,DYNEventPostSolveAlgebraicJacobian,(void*)dyn);CHKERRQ(ierr);
  
  /* Some methods, such as ROSW, compute Jacobian only once at beginning of time-stepping, this may result in an incorrect
     Jacobian for algbraic solve. Need to reset the lag here
  */
  ierr = SNESGetLagJacobian(tssnes,&tssneslag);CHKERRQ(ierr);
  ierr = SNESSetLagJacobian(tssnes,1);CHKERRQ(ierr);

  ierr = SNESSetOptionsPrefix(tssnes,"dyn_");CHKERRQ(ierr);
  ierr = SNESSetFromOptions(tssnes);CHKERRQ(ierr);

  /* Solve algebraic equations */
  ierr = SNESSolve(tssnes,NULL,dyn->X);CHKERRQ(ierr);
  ierr = SNESGetConvergedReason(tssnes,&reason);CHKERRQ(ierr);
  if(reason < 1) {
    SETERRQ1(PETSC_COMM_SELF,0,"The algebraic solution for discontinuity at %3.4f sec did not converge",t);
  }

  /* Reset Jacobian lag */
  ierr = SNESSetLagJacobian(tssnes,tssneslag);CHKERRQ(ierr);
  /* Revert to the original function and jacobian evaluation routines */
  ierr = DMSNESSetFunction(dyn->ps->networkdm,tssnesfunc,tssnesctx);CHKERRQ(ierr);
  ierr = DMSNESSetJacobian(dyn->ps->networkdm,tssnesjacfun,tssnesjacctx);CHKERRQ(ierr);
    
  ierr = PetscObjectStateIncrease((PetscObject)dyn->X);CHKERRQ(ierr);

  /* Reset flag and function, jacobian evaluation routine pointers */
  dyn->solve_alg_only = PETSC_FALSE;

  PetscFunctionReturn(0);
}

/*
  DYNIsSwitchingSolutionStep - Returns if the current step has done a solution after switching, i.e.
     it has solved only the algebraic equations following an event

  Input Parameters
. dyn - the DYN object

  Output Parameters
. flg - TRUE if the current step has done a switching solution.
*/
PetscErrorCode DYNIsSwitchingSolutionStep(DYN dyn, PetscBool *flg)
{
  PetscFunctionBegin;
  *flg = dyn->switching_solution;
  PetscFunctionReturn(0);
}

extern PetscErrorCode TSEventInitialize(TS,PetscReal,Vec);

/*
  DYNInitializeEvents - Reinitializes the events. 

  Notes: An example of usage of this function is when a poststep function
  modifies the solution vector and the events need to be reinitialized
*/
PetscErrorCode DYNInitializeEvents(DYN dyn)
{
  PetscErrorCode ierr;
  TS             ts=dyn->ts;
  PetscReal      t;
  Vec            X;

  PetscFunctionBegin;

  if(!dyn->Nevents) PetscFunctionReturn(0);

  ierr = TSGetTime(ts,&t);CHKERRQ(ierr);
  ierr = TSGetSolution(ts,&X);CHKERRQ(ierr);

  ierr = TSEventInitialize(ts,t,X);CHKERRQ(ierr);

  /* Store the event function values before X changes */
  ierr = DYNEventMonitor(ts,t,X,dyn->feventtmp,(void*)dyn);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode DYNGetCurrentStepSolution(DYN dyn, Vec *X)
{
  PetscFunctionBegin;
  *X = dyn->X;
  PetscFunctionReturn(0);
}

PetscErrorCode DYNSetTimeStep(DYN dyn, PetscReal dt)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;  
  ierr = TSSetTimeStep(dyn->ts,dt);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DYNGetTimeStep(DYN dyn, PetscReal* dt)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = TSGetTimeStep(dyn->ts,dt);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DYNRollbackStep(DYN dyn)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = TSRollBack(dyn->ts);CHKERRQ(ierr);

  ierr = DYNInitializeEvents(dyn);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
