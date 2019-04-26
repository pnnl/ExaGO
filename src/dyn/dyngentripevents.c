#include <private/dynimpl.h>

/* Simple generator tripping event based on a fixed time */

typedef struct{
  PetscInt  gbus; /* Generator bus */
  char      gid[3]; /* Generator id */
  PetscReal tsw;  /* Switching time */
  PetscInt  status; /* Switching status */
}DYNEvent_GenTrip;

PetscErrorCode DYNEventDestroy_GenTrip(DYNEvent dynevent)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscFree(dynevent->data);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DYNEventReadData_GenTrip(DYNEvent event,char* line)
{
  DYNEvent_GenTrip *gentrip=(DYNEvent_GenTrip*)event->data;

  PetscFunctionBegin;
  /* Read the generator tripping event data */
  sscanf(line,"%d,'GENTRIP',%[^,],%lf,%d,",&gentrip->gbus,gentrip->gid,&gentrip->tsw,&gentrip->status);
  PetscFunctionReturn(0);
}

PetscErrorCode DYNEventGetBusNumber_GenTrip(DYNEvent dynevent,PetscInt *busnum)
{
  DYNEvent_GenTrip *gentrip=(DYNEvent_GenTrip*)dynevent->data;
  PetscFunctionBegin;
  *busnum = gentrip->gbus;
  PetscFunctionReturn(0);
}
  
/*
  DYNEventMonitor_GenTrip - Event monitoring function for the generator tripping event. 

  Input Parameters:
+ dynevent - The DYNEvent object
- t        - the current time

  Output Parameters:
. f        - event function residuals
*/
PetscErrorCode DYNEventMonitor_GenTrip(DYNEvent dynevent,PetscReal t, Vec X, PetscScalar *f)
{
  DYNEvent_GenTrip *gentrip=(DYNEvent_GenTrip*)dynevent->data;
  PetscFunctionBegin;

  *f = t - gentrip->tsw;

  PetscFunctionReturn(0);
}

/*
  DYNEventPostFunction_GenTrip - Post event function for generator tripping

  Input Parameters:
+ dynevent   - The DYNEvent object
. nmonitors  - number of monitors that have been triggered for this event
. monidx     - The index for the monitors that have been triggered
. forwarsolve - forward or adjoint solve
. t          - the current time
- X          - Solution vector at time t

  Output Parameters:
. solve_alg  - Flag to indicate whether the algebraic equations need to be resolved

  Notes:
  The forwardsolve flag comes into play only when sensitivities are being computed. The sensitivities in PETSc are computed 
  in PETSc using discrete adjoints that solve the DAE backwards in time. 
  forwardsolve = TRUE => Forward solve
  forwardsolve = FALSE => Adjoint solve
*/
PetscErrorCode DYNEventPostFunction_GenTrip(DYNEvent dynevent,PetscInt nmonitors, PetscInt monidx[], PetscReal t, Vec X, PetscBool forwardsolve,PetscBool *solve_alg)
{
  PetscErrorCode ierr;
  DYNEvent_GenTrip *gentrip=(DYNEvent_GenTrip*)dynevent->data;

  PetscFunctionBegin;
  if(forwardsolve) {
    ierr = DYNEventSetGenStatus(dynevent,gentrip->gbus,gentrip->gid,gentrip->status);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF,"[%4.3f s] Bus %d: Generator %s changed its status to %d\n",t,gentrip->gbus,gentrip->gid,gentrip->status);CHKERRQ(ierr);
    *solve_alg = PETSC_TRUE;
  } else { /* Backward solve */
    ierr = DYNEventSetGenStatus(dynevent,gentrip->gbus,gentrip->gid,1-gentrip->status);CHKERRQ(ierr);
    *solve_alg = PETSC_FALSE;
  }
  PetscFunctionReturn(0);
}


/*
  DYNEventCreate_GenTrip - Class constructor for generator trip event type

  Input Parameters:
+ dynevent - the DYNEvent object
*/
PetscErrorCode DYNEventCreate_GenTrip(DYNEvent dynevent)
{
  DYNEvent_GenTrip *gentrip;
  PetscErrorCode ierr;
  PetscInt       direction=1;
  PetscBool      terminate=PETSC_FALSE;

  PetscFunctionBegin;
  ierr = PetscCalloc1(1,&gentrip);CHKERRQ(ierr);

  dynevent->data = (void*)gentrip;
  dynevent->readdata = DYNEventReadData_GenTrip;
  dynevent->destroy  = DYNEventDestroy_GenTrip;
  dynevent->getbusnumber = DYNEventGetBusNumber_GenTrip;
  ierr = DYNEventSetMonitors(dynevent,1,&direction,&terminate,DYNEventMonitor_GenTrip,DYNEventPostFunction_GenTrip);CHKERRQ(ierr);
  ierr = DYNEventSetLocation(dynevent,DYNEVENT_ON_BUS);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
