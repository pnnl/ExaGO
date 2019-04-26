#include <private/dynimpl.h>

/* Simple line switching event based on a fixed time */

typedef struct{
  PetscInt  fbus; /* From bus */
  PetscInt  tbus; /* To bus */
  PetscReal tsw;  /* Switching time */
  PetscInt  status; /* Switching status */
}DYNEvent_LineSw;

PetscErrorCode DYNEventDestroy_LineSw(DYNEvent dynevent)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscFree(dynevent->data);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DYNEventReadData_LineSw(DYNEvent event,char* line)
{
  DYNEvent_LineSw *linesw=(DYNEvent_LineSw*)event->data;

  PetscFunctionBegin;
  /* Read the linesw data */
  sscanf(line,"%d,%d,'LINESW',%lf,%d,",&linesw->fbus,&linesw->tbus,&linesw->tsw,&linesw->status);
  PetscFunctionReturn(0);
}

PetscErrorCode DYNEventGetFromToBusNumbers_LineSw(DYNEvent dynevent,PetscInt *fbusnum,PetscInt *tbusnum)
{
  DYNEvent_LineSw *linesw=(DYNEvent_LineSw*)dynevent->data;
  PetscFunctionBegin;
  *fbusnum = linesw->fbus;
  *tbusnum = linesw->tbus;
  PetscFunctionReturn(0);
}
  
/*
  DYNEventMonitor_LineSw - Event monitoring function for the linesw event. 

  Input Parameters:
+ dynevent - The DYNEvent object
- t        - the current time

  Output Parameters:
. f        - array of event function residuals
*/
PetscErrorCode DYNEventMonitor_LineSw(DYNEvent dynevent,PetscReal t, Vec X, PetscScalar *f)
{
  DYNEvent_LineSw *linesw=(DYNEvent_LineSw*)dynevent->data;
  PetscFunctionBegin;
  *f = t - linesw->tsw;
  PetscFunctionReturn(0);
}

/*
  DYNEventPostFunction_LineSw - Post event function for linesw

  Input Parameters:
+ dynevent   - The DYNEvent object
. nmonitors  - number of monitors that have been triggered for this event
. monidx     - The index for the monitors that have been triggered
. forwardsolve - Forward or Adjoint solve?
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
PetscErrorCode DYNEventPostFunction_LineSw(DYNEvent dynevent,PetscInt nmonitors, PetscInt monidx[], PetscReal t,  Vec X, PetscBool forwardsolve,PetscBool *solve_alg)
{
  PetscErrorCode ierr;
  DYNEvent_LineSw *linesw=(DYNEvent_LineSw*)dynevent->data;

  PetscFunctionBegin;
  if(forwardsolve) {
    ierr = DYNEventSetLineStatus(dynevent,linesw->status);CHKERRQ(ierr);
    *solve_alg = PETSC_TRUE;
    if(dynevent->eventstring == NULL){
      ierr = PetscMalloc(MAXLINE,&(dynevent->eventstring));CHKERRQ(ierr);
    }
    ierr = PetscSNPrintf(dynevent->eventstring,MAXLINE,"%4.4f,OUT-OF-SERVICE,line,%-6d-%-6d,Line from %d to %d has been out-of-service\n",t,linesw->fbus,linesw->tbus,linesw->fbus,linesw->tbus);CHKERRQ(ierr);
    if(!linesw->status) {
      ierr = PetscPrintf(PETSC_COMM_SELF,"[%4.3f s] Line %d -- %d out of service\n",t,linesw->fbus,linesw->tbus);
    } else {
      ierr = PetscPrintf(PETSC_COMM_SELF,"[%4.3f s] Line %d -- %d reconnected\n",t,linesw->fbus,linesw->tbus);
    }
  } else {
    ierr = DYNEventSetLineStatus(dynevent,(1-linesw->status));CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF,"[%4.3f s] Line %d -- %d reconnected\n",t,linesw->fbus,linesw->tbus);
    *solve_alg = PETSC_FALSE;
  }

  PetscFunctionReturn(0);
}


/*
  DYNEventCreate_LineSw - Class constructor for linesw event type

  Input Parameters:
+ dynevent - the DYNEvent object
*/
PetscErrorCode DYNEventCreate_LineSw(DYNEvent dynevent)
{
  DYNEvent_LineSw *linesw;
  PetscErrorCode ierr;
  PetscInt       direction;
  PetscBool      terminate;

  PetscFunctionBegin;
  ierr = PetscCalloc1(1,&linesw);CHKERRQ(ierr);

  dynevent->data = (void*)linesw;
  dynevent->readdata = DYNEventReadData_LineSw;
  dynevent->destroy  = DYNEventDestroy_LineSw;
  dynevent->getfromtobusnumbers = DYNEventGetFromToBusNumbers_LineSw;
  direction = 1;
  terminate = PETSC_FALSE;
  ierr = DYNEventSetMonitors(dynevent,1,&direction,&terminate,DYNEventMonitor_LineSw,DYNEventPostFunction_LineSw);CHKERRQ(ierr);
  ierr = DYNEventSetLocation(dynevent,DYNEVENT_ON_LINE);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
