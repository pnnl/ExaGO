#include <private/dynimpl.h>

#define FAULT_ON 1
#define FAULT_OFF 0

/* Bus faults */
typedef struct{
  PetscInt  busnum; /* The bus number */
  PetscReal tf; /* Fault on time */
  PetscReal tcl; /* Fault clearing time */
  PetscScalar Gfault; /* Fault conductance (pu) */
  PetscScalar Bfault; /* Fault susceptance (pu) */
}DYNEvent_Fault;

PetscErrorCode DYNEventDestroy_Fault(DYNEvent dynevent)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if(dynevent->eventstring) {
    ierr = PetscFree(dynevent->eventstring);CHKERRQ(ierr);
  }
  ierr = PetscFree(dynevent->data);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DYNEventReadData_Fault(DYNEvent event,char* line)
{
  DYNEvent_Fault *fault=(DYNEvent_Fault*)event->data;

  PetscFunctionBegin;
  /* Read the fault data */
  sscanf(line,"%d,'FAULT',%lf,%lf,%lf,%lf,",&fault->busnum,&fault->tf,&fault->tcl,&fault->Gfault,&fault->Bfault);
  PetscFunctionReturn(0);
}

PetscErrorCode DYNEventGetBusNumber_Fault(DYNEvent dynevent,PetscInt *busnum)
{
  DYNEvent_Fault *fault=(DYNEvent_Fault*)dynevent->data;
  PetscFunctionBegin;
  *busnum = fault->busnum;
  PetscFunctionReturn(0);
}
  
/*
  DYNEventMonitor_Fault - Event monitoring function for the fault event. 

  Input Parameters:
+ dynevent - The DYNEvent object
- t        - the current time

  Output Parameters:
. f        - array of event function residuals
*/
PetscErrorCode DYNEventMonitor_Fault(DYNEvent dynevent,PetscReal t, Vec X, PetscScalar *f)
{
  DYNEvent_Fault *fault=(DYNEvent_Fault*)dynevent->data;
  PetscFunctionBegin;
  f[0] = t - fault->tcl;
  f[1] = t - fault->tf;

  PetscFunctionReturn(0);
}

/*
  DYNEventPostFunction_Fault - Post event function for fault

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
PetscErrorCode DYNEventPostFunction_Fault(DYNEvent dynevent,PetscInt nmonitors, PetscInt monidx[], PetscReal t, Vec X, PetscBool forwardsolve,PetscBool *solve_alg)
{
  PetscErrorCode ierr;
  DYNEvent_Fault *fault=(DYNEvent_Fault*)dynevent->data;
  PetscScalar    Gs = fault->Gfault;
  PetscScalar    Bs = fault->Bfault;
  PetscFunctionBegin;
  if(forwardsolve) {
    if(dynevent->eventstring == NULL){
      ierr = PetscMalloc(MAXLINE,&(dynevent->eventstring));CHKERRQ(ierr);
    }
    if(monidx[0] == FAULT_ON && PetscAbsScalar(t - fault->tf) < 1e-8) { /* Fault on */
      /* Add fault shunt */
      /* Note: Positive Bs means capacitance. Hence -Bs is applied to indicate fault inductance */
      ierr = DYNEventAddBusShunt(dynevent,Gs,-Bs);CHKERRQ(ierr);
      ierr = PetscSNPrintf(dynevent->eventstring,MAXLINE,"%4.4f,FAULT_ON,bus,%-6d,bus number %d has been faulty\n",t,fault->busnum,fault->busnum);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_SELF,"[%4.3f s] Bus %d: Fault ON\n",t,fault->busnum);CHKERRQ(ierr);
    } else if(monidx[0] == FAULT_OFF && PetscAbsScalar(t - fault->tcl) < 1e-8) { /* Fault off */
      /* Remove faulted shunt */
      ierr = DYNEventAddBusShunt(dynevent,Gs,Bs);CHKERRQ(ierr);
      ierr = PetscSNPrintf(dynevent->eventstring,MAXLINE,"%4.4f,FAULT_OFF,bus,%-6d,bus number %d hasn't been faulty\n",t,fault->busnum,fault->busnum);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_SELF,"[%4.3f s] Bus %d: Fault OFF\n",t,fault->busnum);CHKERRQ(ierr);
    }
    *solve_alg = PETSC_TRUE;
  } else { /* Backward solve */
    if(monidx[0] == FAULT_OFF) { 
      /* Add fault shunt */
      /* Note: Positive Bs means capacitance. Hence -Bs is applied to indicate fault inductance */
      ierr = DYNEventAddBusShunt(dynevent,Gs,-Bs);CHKERRQ(ierr);
    } else if(monidx[0] == FAULT_ON) { /* Fault off */
      /* Remove faulted shunt */
      ierr = DYNEventAddBusShunt(dynevent,Gs,Bs);CHKERRQ(ierr);
    }
    *solve_alg = PETSC_FALSE;
  }
  PetscFunctionReturn(0);
}


/*
  DYNEventCreate_Fault - Class constructor for fault event type

  Input Parameters:
+ dynevent - the DYNEvent object
*/
PetscErrorCode DYNEventCreate_Fault(DYNEvent dynevent)
{
  DYNEvent_Fault *fault;
  PetscErrorCode ierr;
  PetscInt       direction[2];
  PetscBool      terminate[2];

  PetscFunctionBegin;
  ierr = PetscCalloc1(1,&fault);CHKERRQ(ierr);

  dynevent->data = (void*)fault;
  dynevent->readdata = DYNEventReadData_Fault;
  dynevent->destroy  = DYNEventDestroy_Fault;
  dynevent->getbusnumber = DYNEventGetBusNumber_Fault;
  direction[0] = direction[1] = 1;
  terminate[0] = terminate[1] = PETSC_FALSE;
  ierr = DYNEventSetMonitors(dynevent,2,direction,terminate,DYNEventMonitor_Fault,DYNEventPostFunction_Fault);CHKERRQ(ierr);
  ierr = DYNEventSetLocation(dynevent,DYNEVENT_ON_BUS);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
