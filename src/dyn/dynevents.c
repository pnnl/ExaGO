#include <private/psimpl.h>
#include <private/dynimpl.h>
#include <petsc/private/petscimpl.h>

PetscInt ndyneventtypesregistered = 0;

struct _p_DYNEventTypesList DYNEventTypesList[MAXDYNEVENTTYPES];

/*
  DYNEventTypeRegister - Registers an event type

  Input Parameters:
+ sname     - model name (string)
- createfunction  - the class constructor
*/
PetscErrorCode DYNEventRegister(const char sname[],PetscErrorCode (*createfunction)(DYNEvent))
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscFunctionBegin;
  for(i=0; i < ndyneventtypesregistered;i++) {
    PetscBool match;
    ierr = PetscStrcmp(DYNEventTypesList[i].name,sname,&match);
    if(match) PetscFunctionReturn(0);
  }
  i = ndyneventtypesregistered;
  if(i == MAXDYNEVENTTYPES) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_SUP,"Exceeded maximum number of event types allowed. Max. allowed = %D",MAXDYNEVENTTYPES);
  ierr = PetscStrcpy(DYNEventTypesList[i].name,sname);CHKERRQ(ierr);
  DYNEventTypesList[i].create = createfunction;
  ndyneventtypesregistered++;
  PetscFunctionReturn(0);
}

extern PetscErrorCode DYNEventCreate_Fault(DYNEvent);
extern PetscErrorCode DYNEventCreate_LineSw(DYNEvent);
extern PetscErrorCode DYNEventCreate_GenTrip(DYNEvent);

/*
  DYNEventTypeRegisterAll - Registers all event types

   This routine is called only once to set up DYNEventTypesList
*/
PetscErrorCode DYNEventRegisterAll(void)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  ierr = DYNEventRegister(DYNEventFault,DYNEventCreate_Fault);CHKERRQ(ierr);
  ierr = DYNEventRegister(DYNEventLineSwitching,DYNEventCreate_LineSw);CHKERRQ(ierr);
  ierr = DYNEventRegister(DYNEventGenTrip,DYNEventCreate_GenTrip);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  DYNEventSetLineStatus - Sets the status of the line associated with this event

  Input Parameters:
+ dynevent - the dynevent
- status   - (0 or 1) the status of the line
*/
PetscErrorCode DYNEventSetLineStatus(DYNEvent dynevent,PetscInt status)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = PSLINESetStatus(dynevent->line,status);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNEventAddBusShunt - Adds shunt conductance and susceptance at the incident bus for this dynevent

  Input Parameters:
+ dynevent - The DYNEvent object
. Gs       - shunt conductance
- Bs       - shunt susceptance

  Notes: Gs and Bs should be in per unit
*/
PetscErrorCode DYNEventAddBusShunt(DYNEvent dynevent, PetscScalar Gs, PetscScalar Bs)
{
  PetscErrorCode ierr;
  PSBUS bus=dynevent->bus;

  PetscFunctionBegin;
  if(dynevent->location != DYNEVENT_ON_BUS) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Cannot add shunt for a non-bus event");

  ierr = PSBUSAddShunt(bus,Gs,Bs);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNEventSetGenStatus - Sets the status of the generator incident on the bus

  Input Parameters:
+ dynevent - The DYNEvent object
. gbus     - Generator bus
. gid      - Generator id
- status   - generator status

*/
PetscErrorCode DYNEventSetGenStatus(DYNEvent dynevent, PetscInt gbus,char gid[],PetscInt status)
{
  PetscErrorCode ierr;
  PSBUS bus=dynevent->bus;

  PetscFunctionBegin;
  if(dynevent->location != DYNEVENT_ON_BUS) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Cannot change gen. status for a non-bus event");

  ierr = PSBUSSetGenStatus(bus,gid,status);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


/*
  DYNEventReadEventType - Parses a line in the event file and reads the model type (string identifier for the
                          the registered event model)

  Input Parameters:
. line - the line to parse

. Output Parameters:
. eventtype - type of the event (string name of the event)
*/
PetscErrorCode DYNEventReadEventType(char *line, char eventtype[])
{
  PetscErrorCode ierr;
  PetscInt       i;
  char           *pos;

  PetscFunctionBegin;
  for(i=0; i < ndyneventtypesregistered; i++) {
    /* Strip new line character if any */
    if((pos = strchr(line,'\n')) != NULL) *pos = '\0';
    if(strstr(line,DYNEventTypesList[i].name) != NULL) {
      ierr = PetscStrcpy(eventtype,DYNEventTypesList[i].name);CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }
  }

  SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_SUP,"Unknown event type in event data type\n %s",line);
  PetscFunctionReturn(0);
}


PetscErrorCode DYNEventDestroy(DYNEvent *dynevent)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if((*dynevent)->destroy) {
    ierr = (*(*dynevent)->destroy)(*dynevent);CHKERRQ(ierr);
  }
  ierr = PetscFree((*dynevent)->direction);CHKERRQ(ierr);
  ierr = PetscFree((*dynevent)->terminate);CHKERRQ(ierr);
  ierr = PetscFree(*dynevent);CHKERRQ(ierr);
  *dynevent = 0;
  PetscFunctionReturn(0);
}

PetscErrorCode DYNEventSetLocation(DYNEvent dynevent,DYNEventLocation location)
{
  PetscFunctionBegin;
  dynevent->location = location;
  PetscFunctionReturn(0);
}


/*
  DYNEventCreate - Creates an instance of DYNEvent

  Output Parameters:
. dyneventout - the DYNEvent object
*/
PetscErrorCode DYNEventCreate(DYNEvent *dyneventout)
{
  PetscErrorCode ierr;
  DYNEvent       dynevent;

  PetscFunctionBegin;
  ierr = PetscCalloc1(1,&dynevent);CHKERRQ(ierr);
  dynevent->numMonitors = 0;
  dynevent->location = DYNEVENT_LOCATION_NOTSET;
  dynevent->eventstring = NULL;
  *dyneventout = dynevent;
  PetscFunctionReturn(0);
}

/*
  DYNEventSetType - Sets the type of the event and calls the implementation constructor

  Input Parameters:
+ dynevent - the event
- eventtype - the event type
*/
PetscErrorCode DYNEventSetType(DYNEvent dynevent,const char eventtype[])
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscBool      match;

  PetscFunctionBegin;
  for(i=0; i < ndyneventtypesregistered; i++) {
    ierr = PetscStrcmp(DYNEventTypesList[i].name,eventtype,&match);CHKERRQ(ierr);
    if(match) {
      ierr = (*DYNEventTypesList[i].create)(dynevent);CHKERRQ(ierr);

      ierr = PetscStrcpy(dynevent->name,DYNEventTypesList[i].name);CHKERRQ(ierr);
      break;
    }
  }
  PetscFunctionReturn(0);
}

/*
  DYNEventGetFromToBusNumbers - Gets From and To bus numbers for the line on which the event is located
*/
PetscErrorCode DYNEventGetFromToBusNumbers(DYNEvent dynevent,PetscInt *fbus,PetscInt *tbus)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dynevent->getfromtobusnumbers)(dynevent,fbus,tbus);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DYNEventGetBusNumber(DYNEvent dynevent,PetscInt *busnum)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dynevent->getbusnumber)(dynevent,busnum);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNEventMonitor - Called by TS after every step

  Input Parameters:
+ ts - the TS object
. t  - the current time
. f  - array of event function residuals
. ctx - application specific context

  Note: DYNEventMonitor() is also called by the post event function to set the feventtmp array. This array
        is used in checking if any machine device limits have been violated.
*/
PetscErrorCode DYNEventMonitor(TS ts, PetscReal t, Vec X, PetscScalar *f,void *ctx)
{
  PetscErrorCode ierr;
  DYN dyn=(DYN)ctx;
  DYNEvent dynevent;
  PetscScalar *fval;
  PetscInt   i;

  PetscFunctionBegin;
  fval = f;
  for(i=0; i < dyn->nevents; i++) {
    dynevent = dyn->events[dyn->event_idx[i]];
    ierr = (*dynevent->eventmonitor)(dynevent,t,X,fval);CHKERRQ(ierr);
    fval += dynevent->numMonitors;
  }

  PetscFunctionReturn(0);
}

/*
  DYNEventPostSolveAlgebraicJacobian - Jacobian routine for the post event algebraic solve
*/
PetscErrorCode DYNEventPostSolveAlgebraicJacobian(SNES snes,Vec X,Mat J, Mat Jpre,void* ctx)
{
  PetscErrorCode ierr;
  DYN            dyn=(DYN)ctx;
  PetscReal      t;

  PetscFunctionBegin;

  ierr = TSGetTime(dyn->ts,&t);CHKERRQ(ierr);

  ierr = DYNIJacobian(dyn->ts,t,X,(Vec)NULL, 0.0,J,Jpre,ctx);CHKERRQ(ierr);

  /* Zero out the rows corresponding to the different equations and set 1 on their diagonal */
  ierr = MatSetOption(J,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);CHKERRQ(ierr);
  ierr = MatZeroRowsIS(J,dyn->isdiff,1.0,NULL,NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNEventPostSolveAlgebraic - Residual function for SNES to solve the algebraic part when an event
    is triggered
*/
PetscErrorCode DYNEventPostSolveAlgebraic(SNES snes, Vec X, Vec F,void* ctx)
{
  PetscErrorCode ierr;
  DYN            dyn=(DYN)ctx;
  PetscReal      t;
  const PetscInt *idx;
  PetscInt lsize;
  PetscScalar val[5];
  PetscInt    nval=5,ctr;

  PetscFunctionBegin;

  ierr = TSGetTime(dyn->ts,&t);CHKERRQ(ierr);

  ierr = DYNIFunction(dyn->ts,t,X,(Vec)NULL,F,ctx);CHKERRQ(ierr);

  ierr = ISGetLocalSize(dyn->isdiff,&lsize);CHKERRQ(ierr);
  ierr = ISGetIndices(dyn->isdiff,&idx);CHKERRQ(ierr);

  val[0] = val[1] = val[2] = val[3] = val[4] = 0.0;
  ctr = 0;
  while(ctr < lsize) {
    if(ctr+nval >= lsize) nval = lsize - ctr;
    ierr = VecSetValues(F,nval,idx+ctr,val,INSERT_VALUES);CHKERRQ(ierr);
    ctr += nval;
  }
  ierr = VecAssemblyBegin(F);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F);CHKERRQ(ierr);

  //  ierr = VecView(F,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNEventPostSolveAlgebraic - Residual function for SNES to solve the algebraic part when an event
    is triggered
*/
PetscErrorCode DYNEventPostSolveAlgebraicSensi(SNES snes, Vec X, Vec F,void* ctx)
{
  PetscErrorCode ierr;
  DYN            dyn=(DYN)ctx;
  PetscReal      t;
  const PetscInt *idx;
  PetscInt lsize;
  PetscScalar val[5];
  PetscInt    nval=5,ctr;

  PetscFunctionBegin;

  ierr = TSGetTime(dyn->ts,&t);CHKERRQ(ierr);
  ierr = DYNIFunction(dyn->ts,t,X,(Vec)NULL,F,ctx);CHKERRQ(ierr);
  /* X0 is a vector of 0s*/
  ierr = DYNIFunction(dyn->ts,t,dyn->X0,(Vec)NULL,dyn->Xtmp,ctx);CHKERRQ(ierr);
  ierr = VecAXPY(F,-1,dyn->Xtmp);CHKERRQ(ierr); /* assume algebraic part is linear */

  ierr = ISGetLocalSize(dyn->isdiff,&lsize);CHKERRQ(ierr);
  ierr = ISGetIndices(dyn->isdiff,&idx);CHKERRQ(ierr);

  val[0] = val[1] = val[2] = val[3] = val[4] = 0.0;
  ctr = 0;
  while(ctr < lsize) {
    if(ctr+nval >= lsize) nval = lsize - ctr;
    ierr = VecSetValues(F,nval,idx+ctr,val,INSERT_VALUES);CHKERRQ(ierr);
    ctr += nval;
  }
  ierr = VecAssemblyBegin(F);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode DYNShiftGradients(DYN dyn)
{
  PetscScalar       innerprod1,innerprod2;
#ifdef FWDSA
  const PetscScalar *dgdp;
  PetscInt  xdynsize=0;
#endif
  PetscInt          i;
  PetscErrorCode    ierr;
  PetscFunctionBegin;

  ierr = VecDot(dyn->dGammadX,dyn->Xtmp,&innerprod1);
  if(PetscAbsScalar(innerprod1) < PETSC_SMALL) PetscFunctionReturn(0);

  ierr = VecAXPY(dyn->Xtmp2,-1,dyn->Xtmp);CHKERRQ(ierr); /* f^{(2)} - f^{(1)} */
#ifdef FWDSA
  if(dyn->useforward) {
    ierr = VecGetArrayRead(dyn->dGammadP,&dgdp);CHKERRQ(ierr);
    for(i=0;i<dyn->nparams;i++) {
      ierr = VecDot(dyn->dGammadX,dyn->sp[i],&innerprod2);
      ierr = VecAXPY(dyn->sp[i],innerprod2/innerprod1,dyn->Xtmp2);CHKERRQ(ierr);
      ierr = VecAXPY(dyn->sp[i],dgdp[i]/innerprod1,dyn->Xtmp2);CHKERRQ(ierr);
    }
    /* Comment out the following line to disable computing sensitivity to initial values */
    /* ierr = VecGetSize(dyn->X,&xdynsize);CHKERRQ(ierr); */

    for(i=0;i<xdynsize;i++) {
      ierr = VecDot(dyn->dGammadX,dyn->s[i],&innerprod2);
      ierr = VecAXPY(dyn->s[i],innerprod2/innerprod1,dyn->Xtmp2);CHKERRQ(ierr);
    }
    ierr = VecRestoreArrayRead(dyn->dGammadP,&dgdp);CHKERRQ(ierr);
  }
#endif
  if(dyn->useadjoint) {
    for(i=0; i<dyn->ncostfcns; i++) {
      ierr = VecDot(dyn->Xtmp2,dyn->lambda[i],&innerprod2); /* (f^{2}-f^{(1)})*lambda */
      ierr = VecAXPY(dyn->lambda[i],innerprod2/innerprod1,dyn->dGammadX);CHKERRQ(ierr);
      ierr = VecAXPY(dyn->mu[i],innerprod2/innerprod1,dyn->dGammadP);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

/*
  DYNEventPostFunction - Post event function routine
*/
PetscErrorCode DYNEventPostFunction(TS ts, PetscInt nmonitors, PetscInt monitor_list[], PetscReal t, Vec X, PetscBool forwardsolve,void *ctx)
{
  PetscErrorCode ierr;
  DYN dyn =(DYN)ctx;
  PS  ps  = dyn->ps;
  DYNEvent dynevent;
  PetscInt j=0,neventmonitors,eventmonitorctr=0;
  PetscInt eventmonitor_list[MAXDYNEVENTSLOCATEDATSTEP]={0};
  PetscInt eventctr=0;
  PetscBool solve_algebraic_eqs=PETSC_FALSE,snes_solve;
  Vec       localX,localXprv;

  /*Event File Variables(SAM)*/
  char           event_file[PETSC_MAX_PATH_LEN];
  PetscViewer    eventview;


  PetscFunctionBegin;
  if(nmonitors > MAXDYNEVENTSLOCATEDATSTEP) {
    SETERRQ3(PETSC_COMM_SELF,0,"Number of events located %D > Max. allowed %d at time %4.2f\n",nmonitors,MAXDYNEVENTSLOCATEDATSTEP,t);
  }

  if ((dyn->useforward && forwardsolve) || (dyn->useadjoint && !forwardsolve)) {
    ierr = VecSet(dyn->dGammadX,0.0);CHKERRQ(ierr);
    ierr = VecSet(dyn->dGammadP,0.0);CHKERRQ(ierr);
  }
    /* Xtmp will store the value for pre-event RHS function and Xtmp2 will store the value for post-event RHS function */
  if (dyn->useforward && forwardsolve) { /* first time calling postevent() */
    /* pre-event RHS evaluation */
    ierr = VecSet(dyn->Xtmp,0.0);CHKERRQ(ierr);
    ierr = DYNIFunction(dyn->ts,t,dyn->X,(Vec)NULL,dyn->Xtmp,ctx);CHKERRQ(ierr);
  }
  if(dyn->useadjoint && !forwardsolve) {
    /* post-event RHS evaluation */
    ierr = VecSet(dyn->Xtmp2,0.0);CHKERRQ(ierr);
    ierr = DYNIFunction(dyn->ts,t,dyn->X,(Vec)NULL,dyn->Xtmp2,ctx);CHKERRQ(ierr);
  }

  ierr = DMGetLocalVector(ps->networkdm,&localXprv);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ps->networkdm,dyn->X,INSERT_VALUES,localXprv);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ps->networkdm,dyn->X,INSERT_VALUES,localXprv);CHKERRQ(ierr);

  ierr = DMGetLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(ps->networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);


  while(j < nmonitors) { 
    dynevent = dyn->events[dyn->event_idx[eventctr]];
    neventmonitors = 0;
    while(monitor_list[j] < eventmonitorctr + dynevent->numMonitors) {
      eventmonitor_list[neventmonitors++] = monitor_list[j++] - eventmonitorctr;
      if(j == nmonitors) break;
    }
    if(neventmonitors && dynevent->posteventfunction) {
      if (dyn->useforward && forwardsolve && dynevent->posteventdgamma) {
        ierr = (*dynevent->posteventdgamma)(dynevent,neventmonitors,eventmonitor_list,t,localXprv,forwardsolve);CHKERRQ(ierr);
      }
      PetscBool solve_algebraic=PETSC_FALSE;
#if 0
      ierr = PetscPrintf(PETSC_COMM_SELF,"[P%d]:Event Detected: %s at t = %lf\n",dyn->comm->rank,dynevent->name,t);CHKERRQ(ierr);
#endif
      ierr = (*dynevent->posteventfunction)(dynevent,neventmonitors,eventmonitor_list,t,localX,forwardsolve,&solve_algebraic);CHKERRQ(ierr);
      solve_algebraic_eqs = solve_algebraic_eqs || solve_algebraic;


      /*Event File Create(SAM)*/
      if(dyn->visualize && dynevent->eventstring != NULL){
        ierr = PetscSNPrintf(event_file,sizeof(event_file),"%s/event-%06d.out",dyn->viewer_dir,++(dyn->eventlogcount));CHKERRQ(ierr);
        ierr = PetscViewerASCIIOpen(dyn->comm->type,event_file,&eventview);CHKERRQ(ierr);
        ierr = PetscViewerASCIIPrintf(eventview,"Time   EventTag  Bus/Branch number  Message\n");CHKERRQ(ierr);

        ierr = PetscViewerASCIIPrintf(eventview,"%s",dynevent->eventstring);CHKERRQ(ierr);
        
        /*Event File Destroy(SAM)*/
        ierr = PetscViewerDestroy(&eventview);CHKERRQ(ierr);
      }

    }
    eventmonitorctr += dynevent->numMonitors;
    eventctr++;
  }

  ierr = DMRestoreLocalVector(ps->networkdm,&localXprv);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);

  ierr = MPI_Allreduce(&solve_algebraic_eqs,&snes_solve,1,MPIU_BOOL,MPI_LOR,dyn->comm->type);CHKERRQ(ierr);
  if(snes_solve) {
    dyn->switching_solution = PETSC_TRUE;
    /* Solve only the algebraic part to obtain consistent initial conditions for the DAE. For this, the differential
       variables are held constant. We reuse the SNES associated with the TS to do the nonlinear solve
    */
    SNES tssnes;
    void *tssnesctx;
    void *tssnesjacctx;
    PetscErrorCode (*tssnesfunc)(SNES,Vec,Vec,void*);
    PetscErrorCode (*tssnesjacfun)(SNES,Vec,Mat,Mat,void*);
    PetscInt tssneslag;
    SNESConvergedReason reason;


    /* Store the event function values before X changes */
    ierr = DYNEventMonitor(ts,t,X,dyn->feventtmp,(void*)dyn);CHKERRQ(ierr);

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
    
#ifdef FWDSA
    /* Solve only the algebraic part to obtain consistency initial conditions for the sensitivity equations. */
    if (dyn->useforward) {
      ierr = TSGetSNES(dyn->ts,&tssnes);CHKERRQ(ierr);
      ierr = SNESSetFunction(tssnes,tssnesres,DYNEventPostSolveAlgebraicSensi,(void*)dyn);CHKERRQ(ierr);
      ierr = DYNForwardAlgebraicSolve(tssnes,dyn,t,X);
    }
#endif
    
    /* Check and set machine limits that may be violated */
    ierr = DYNCheckandSetMachineLimits(dyn,t,X,forwardsolve);CHKERRQ(ierr);

    /* Update feventtmp */
    ierr = DYNEventMonitor(ts,t,X,dyn->feventtmp,(void*)dyn);CHKERRQ(ierr);


  }

  if (dyn->useadjoint && !forwardsolve) {
    ierr = DMGetLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(ps->networkdm,dyn->X,INSERT_VALUES,localX);CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(ps->networkdm,dyn->X,INSERT_VALUES,localX);CHKERRQ(ierr);

    j               = 0;
    eventctr        = 0;
    eventmonitorctr = 0;
    while(j < nmonitors) {
      dynevent = dyn->events[dyn->event_idx[eventctr]];
      neventmonitors = 0;
      while(monitor_list[j] < eventmonitorctr + dynevent->numMonitors) {
        eventmonitor_list[neventmonitors++] = monitor_list[j++] - eventmonitorctr;
        if(j == nmonitors) break;
      }
      if(neventmonitors && dynevent->posteventdgamma) {
        ierr = (*dynevent->posteventdgamma)(dynevent,neventmonitors,eventmonitor_list,t,localX,forwardsolve);CHKERRQ(ierr);
      }
      eventmonitorctr += dynevent->numMonitors;
      eventctr++;
    }
    ierr = DMRestoreLocalVector(ps->networkdm,&localX);CHKERRQ(ierr);
  }

  if(dyn->useforward && forwardsolve) {
    /* post-event RHS evaluation */
    ierr = VecSet(dyn->Xtmp2,0.0);CHKERRQ(ierr);
    ierr = DYNIFunction(dyn->ts,t,dyn->X,(Vec)NULL,dyn->Xtmp2,ctx);CHKERRQ(ierr);
    /* shift the gradients */
    ierr = DYNShiftGradients(dyn);CHKERRQ(ierr);
  }
  if(dyn->useadjoint && !forwardsolve) {
    /* pre-event RHS evaluation */
    ierr = VecSet(dyn->Xtmp,0.0);CHKERRQ(ierr);
    ierr = DYNIFunction(dyn->ts,t,dyn->X,(Vec)NULL,dyn->Xtmp,ctx);CHKERRQ(ierr);
    /* shift the gradients */
    ierr = DYNShiftGradients(dyn);CHKERRQ(ierr);
  }

  ierr = PetscObjectStateIncrease((PetscObject)dyn->X);CHKERRQ(ierr);

  /* Reset flag and function, jacobian evaluation routine pointers */
  dyn->solve_alg_only = PETSC_FALSE;
  PetscFunctionReturn(0);
}

PetscErrorCode DYNSetUpEvents(DYN dyn)
{
  PetscErrorCode ierr;
  PetscInt       i,busnum;
  PS             ps=dyn->ps;
  DYNEvent       dynevent;

  PetscFunctionBegin;
  dyn->toteventmonitors = 0;
  for(i = 0; i < dyn->Nevents; i++) {
    if(dyn->events[i]->location == DYNEVENT_LOCATION_NOTSET) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Location for the event must be set using DYNEventSetLocation()");
    if(dyn->events[i]->location == DYNEVENT_ON_BUS) {       /* Event on bus */
      /* Get the bus number */
      ierr = DYNEventGetBusNumber(dyn->events[i],&busnum);CHKERRQ(ierr);
      if(ps->busext2intmap[busnum] != -1) { /* The bus is on this processor */
	dyn->events[i]->bus = &ps->bus[ps->busext2intmap[busnum]];
	dyn->event_idx[dyn->nevents] = i;
	dyn->toteventmonitors += dyn->events[i]->numMonitors;
	dyn->nevents++;
      }
    } else if(dyn->events[i]->location == DYNEVENT_ON_LINE) {
      PetscInt fbusnum,tbusnum,k;
      PSBUS    fbus;
      const PSLINE   *supplines;
      PetscInt nsupplines;
      /* Get the from and to bus numbers */
      ierr = DYNEventGetFromToBusNumbers(dyn->events[i],&fbusnum,&tbusnum);CHKERRQ(ierr);
      if(ps->busext2intmap[fbusnum] != -1) {
	fbus = &ps->bus[ps->busext2intmap[fbusnum]];
	/* Get connected lines at fbus */
	ierr = PSBUSGetSupportingLines(fbus,&nsupplines,&supplines);CHKERRQ(ierr);
	for(k=0; k < nsupplines; k++) {
	  if(supplines[k]->tbus == tbusnum) {
	    dyn->events[i]->line = supplines[k];
	    dyn->event_idx[dyn->nevents] = i;
	    dyn->toteventmonitors += dyn->events[i]->numMonitors;
	    dyn->nevents++;
	    break;
	  }
	}
      }
    }
  }

  /* Add events on machines (limiters, saturation, and others) */
  ierr = DYNSetMachineEvents(dyn);CHKERRQ(ierr);

  ierr = PetscCalloc1(dyn->toteventmonitors,&dyn->eventdirection);CHKERRQ(ierr);
  ierr = PetscCalloc1(dyn->toteventmonitors,&dyn->eventterminate);CHKERRQ(ierr);
  /* Copy over the directions and termination flags */
  PetscInt j,ctr=0;
  for(i=0; i < dyn->nevents; i++ ) {
    dynevent = dyn->events[dyn->event_idx[i]];
    for(j = 0; j < dynevent->numMonitors; j++,ctr++) {
      dyn->eventdirection[ctr] = dynevent->direction[j];
      dyn->eventterminate[ctr] = dynevent->terminate[j];
    }
  }

  /* feventtmp array is only used when checking limits at disturbance on/off times */
  ierr = PetscCalloc1(dyn->toteventmonitors,&dyn->feventtmp);CHKERRQ(ierr);

  /* Set up vectors to store derivertives needed by jump condition */
  if (dyn->useforward || dyn->useadjoint) {
    ierr = VecDuplicate(dyn->X,&dyn->dGammadX);CHKERRQ(ierr);
    ierr = VecSet(dyn->dGammadX,0.0);CHKERRQ(ierr);
    ierr = VecDuplicate(dyn->Xparam,&dyn->dGammadP);CHKERRQ(ierr);
    ierr = VecSet(dyn->dGammadP,0.0);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*
  DYNEventReadDataFromLine - Reads the event data from the line in the file
*/
PetscErrorCode DYNEventReadDataFromLine(DYNEvent dynevent,char *line)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*dynevent->readdata)(dynevent,line);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DYNEventSetMonitors - Sets the info for monitoring the event

  Input Parameters:
+ dynevent    - The DYNEvent object
. numMonitors - number of monitors for this event
. direction   - directions for zero crossing detection for each monitor ( 1 -> going positive, -1 -> going negative, 0 for both directions)
. terminate   - termination flags for each monitor if one wishes to terminate the time-stepping after an event is detected
. eventmonitor - the monitoring function for this event. This will be called after each time-step
. posteventfunction - optional postevent function to allow taking actions after an event has occured. This routine is called after                       an event has occured
*/
PetscErrorCode DYNEventSetMonitors(DYNEvent dynevent,PetscInt numMonitors,PetscInt direction[],PetscBool terminate[],PetscErrorCode (*eventmonitor)(DYNEvent,PetscReal,Vec,PetscScalar*), PetscErrorCode (*posteventfunction)(DYNEvent,PetscInt,PetscInt[],PetscReal,Vec,PetscBool,PetscBool*))
{
  PetscErrorCode ierr;
  PetscInt i;

  PetscFunctionBegin;

  dynevent->numMonitors = numMonitors;
  ierr = PetscCalloc1(numMonitors,&dynevent->direction);CHKERRQ(ierr);
  ierr = PetscCalloc1(numMonitors,&dynevent->terminate);CHKERRQ(ierr);
  for(i=0; i < numMonitors; i++) {
    dynevent->direction[i] = direction[i];
    dynevent->terminate[i] = terminate[i];
  }
  dynevent->eventmonitor = eventmonitor;
  dynevent->posteventfunction = posteventfunction;
  dynevent->posteventdgamma = NULL; /* no jump condition by default, but machine events come with NONNULL pointers */
  PetscFunctionReturn(0);
}

/*
  DYNReadEventData - Reads the event data file and sets up the events in the DYN object

  Input Parameters:
+ DYN - the dynamics simulation object
- eventfile - the name of the event file

  Notes:
    Each processor reads the data file.
*/
PetscErrorCode DYNReadEventData(DYN dyn,const char eventfile[])
{
  FILE           *fp;
  char           line[MAXLINE];
  char           eventtype[16];
  PetscErrorCode ierr;

  PetscFunctionBegin;

  fp = fopen(eventfile,"r");
  /* Check for valid file */
  if (fp == NULL) {
    SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Cannot open file %s",eventfile);
  }

  dyn->Nevents = dyn->nevents = 0;
  while(fgets(line,MAXLINE,fp) != NULL) {
    /* Get event type */
    ierr = DYNEventReadEventType(line,eventtype);CHKERRQ(ierr);
    ierr = DYNEventCreate(&dyn->events[dyn->Nevents]);CHKERRQ(ierr);
    ierr = DYNEventSetType(dyn->events[dyn->Nevents],eventtype);CHKERRQ(ierr);
    ierr = DYNEventReadDataFromLine(dyn->events[dyn->Nevents],line);CHKERRQ(ierr);

    dyn->Nevents++;
  }
  fclose(fp);

  PetscFunctionReturn(0);
}
