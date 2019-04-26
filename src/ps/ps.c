#include <private/psimpl.h>
#include <petscdmnetwork.h>

/*
  PSConnCompDestroy - Destroys the connected components struct
*/
PetscErrorCode PSConnCompDestroy(PS ps)
{
  PetscErrorCode ierr;
  PetscInt i;
  PetscFunctionBegin;
  for(i=0; i < ps->nconncomp; i++) {
    ps->conncomp[i].nv = 0;
    ierr = PetscFree(ps->conncomp[i].v);CHKERRQ(ierr);
  }
  ps->nconncomp = 0;
  PetscFunctionReturn(0);
}

/*
  PSCheckandSetRefBus - Checks for active ref. bus and sets one if not available

  Input Parameters:
. ps - The PS object

  This routine checks if there is an active ref. bus (declared in the data file and has active generators). If it is
  not defined, then the first PV bus is used.
*/
PetscErrorCode PSCheckandSetRefBus(PS ps)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PSBUS          bus;
  PetscBool      pshaslocalref=PETSC_FALSE, pshasref;
  PetscInt       firstlocpvbus=100000000,firstpvbus;

  PetscFunctionBegin;
  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    if(bus->ide == REF_BUS && bus->ngenON) {
      /* PS has active reference bus */
      pshaslocalref = PETSC_TRUE;
      break;
    }
  }

  ierr = MPI_Allreduce(&pshaslocalref,&pshasref,1,MPIU_BOOL,MPI_LOR,ps->comm->type);CHKERRQ(ierr);

  if(pshasref) PetscFunctionReturn(0);
#if defined DEBUGPS
  ierr = PetscPrintf(PETSC_COMM_SELF,"No active ref. bus in the system\n");CHKERRQ(ierr);
#endif
  /* No active ref. bus, set first PV bus to ref bus */
  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    if(bus->ide == PV_BUS && bus->ngenON && !bus->isghost) {
      firstlocpvbus = bus->bus_i;
      break;
    }
  }

  ierr = MPIU_Allreduce(&firstlocpvbus,&firstpvbus,1,MPIU_INT,MPI_MIN,ps->comm->type);CHKERRQ(ierr);

  if(firstpvbus == 100000000) {
#if defined DEBUGPS
      SETERRQ1(PETSC_COMM_SELF,0," ",firstpvbus);
#endif
  } else {
    if (ps->busext2intmap[firstpvbus] != -1) {
      ps->bus[ps->busext2intmap[firstpvbus]].ide = REF_BUS;
      ierr = PetscPrintf(PETSC_COMM_SELF," ",firstpvbus);CHKERRQ(ierr);
    }
  }

  PetscFunctionReturn(0);
}

/*
  PSIslandCheckandSetRefBus - Checks for active ref. bus and sets one if not available. This
    is a version of PSCheckandSetRefBus that allows setting reference buses on islands.

  Input Parameters:
+ ps - The PS object
- isnum  - island number

  This routine checks if there is an active ref. bus (declared in the data file and has active generators). If it is
  not defined, then the first PV bus is used.
*/
PetscErrorCode PSIslandCheckandSetRefBus(PS ps,PetscInt isnum)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PSBUS          bus;
  PetscBool      pshaslocalref=PETSC_FALSE, pshasref;
  PetscInt       firstlocpvbus=100000000,firstpvbus;

  PetscFunctionBegin;
  if(ps->conncomp[isnum].blackout) PetscFunctionReturn(0);

  for(i=0; i < ps->conncomp[isnum].nv; i++) {
    bus = &ps->bus[ps->busext2intmap[ps->conncomp[isnum].v[i]]];
    if(bus->ide == REF_BUS && bus->ngenON) {
      /* PS has active reference bus */
      pshaslocalref = PETSC_TRUE;
      break;
    }
  }

  ierr = MPI_Allreduce(&pshaslocalref,&pshasref,1,MPIU_BOOL,MPI_LOR,ps->comm->type);CHKERRQ(ierr);

  if(pshasref) PetscFunctionReturn(0);
#if defined DEBUGPS
  ierr = PetscPrintf(PETSC_COMM_SELF,"No active ref. bus in island %d\n",isnum+1);CHKERRQ(ierr);
#endif
  /* No active ref. bus, set first PV bus to ref bus */
  for(i=0; i < ps->conncomp[isnum].nv; i++) {
    bus = &ps->bus[ps->busext2intmap[ps->conncomp[isnum].v[i]]];
    if(bus->ide == PV_BUS && bus->ngenON && !bus->isghost) {
      firstlocpvbus = bus->bus_i;
      break;
    }
  }

  ierr = MPIU_Allreduce(&firstlocpvbus,&firstpvbus,1,MPIU_INT,MPI_MIN,ps->comm->type);CHKERRQ(ierr);

  if(firstpvbus == 100000000) {
#if defined DEBUGPS
      SETERRQ1(PETSC_COMM_SELF,0," ",firstpvbus);
#endif
  } else {
    if (ps->busext2intmap[firstpvbus] != -1) {
      ps->bus[ps->busext2intmap[firstpvbus]].ide = REF_BUS;
      ierr = PetscPrintf(PETSC_COMM_SELF,"Setting bus %d as the new ref. bus\n ",firstpvbus);CHKERRQ(ierr);
    }
  }

  PetscFunctionReturn(0);
}

/*
  PSSetGenDispatchandStatus - Sets the generator status and dispatch

  Input Parameters:
+ ps - the ps object
. busnum - the bus number 
. gennum - generator number
. status - the generator status (0 for OFF, 1 for ON)
. pg - active power dispatch in MW
- qg - reactive power dispatch in MVAr

  Notes: Must be called before PSSetUp() is called
*/
PetscErrorCode PSSetGenDispatchandStatus(PS ps, PetscInt busnum, PetscInt gennum, PetscInt status, PetscScalar Pg, PetscScalar Qg)
{
  PetscInt       intbusnum;
  PSBUS          bus;
  PSGEN          gen;
  PetscBool      statusupdate=PETSC_FALSE;

  PetscFunctionBegin;

  //  if(ps->setupcalled) {
  //    SETERRQ(PETSC_COMM_SELF,0,"Must call PSSetGenDispatchandStatus before PSSetUp() is called");
  //  }
  intbusnum = ps->busext2intmap[busnum];
  if(intbusnum == -1) PetscFunctionReturn(0); /* Not local to this processor */
  bus = &ps->bus[intbusnum];
  if(!bus->ngen) {
    SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_SUP,"No generators on bus %d",busnum);
  }
  if(bus->ide == ISOLATED_BUS) {
    SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_SUP,"Bus %d is isolated",busnum);
  }

  gen=&ps->gen[bus->gidx[gennum]];
  if(status != gen->status) statusupdate = PETSC_TRUE;
  gen->status = status;
  gen->pg = Pg/ps->MVAbase;
  gen->qg = Qg/ps->MVAbase;

  if(statusupdate) {
    if(status) {
      bus->ngenON++;
      if(bus->ide == PQ_BUS) bus->ide = PV_BUS;
    } else {
      bus->ngenON--;
      if(!bus->ngenON) bus->ide = PQ_BUS;
      gen->pg = gen->qg = 0.0;
    }
  }
    
  PetscFunctionReturn(0);
}

/*
  PSGetNumGenerators - Gets the number of local and global generators in this system

  Input Parameters
. ps - the PS object

  Output Parameters
+ ngen - number of local generators
- Ngen - number of global generators
 
  Notes:
   PSSetUp() must be called before a call to PSGetNumGenerators
*/
PetscErrorCode PSGetNumGenerators(PS ps,PetscInt *ngen, PetscInt *Ngen)
{
  PetscFunctionBegin;
  if(ngen) *ngen = ps->ngen;
  if(Ngen) *Ngen = ps->Ngen;
  PetscFunctionReturn(0);
}

/*
  PSLINESetStatus - Sets the status of the line

  Input Parameters:
+ line - the line
- status   - (0 or 1) the status of the line
*/
PetscErrorCode PSLINESetStatus(PSLINE line,PetscInt status)
{
  PetscFunctionBegin;
  line->status = status;
  PetscFunctionReturn(0);
}

/*
  PSGetLine - Gets the line object given the from bus, to bus, and id

  Input Parameters
+ ps - the PS object
. fbus - from bus
. tbus - to bus
- id   - line id

  Output Parameters
. line - the PSLine object

*/
PetscErrorCode PSGetLine(PS ps,PetscInt fbus, PetscInt tbus, const char* id,PSLINE *line)
{
  PetscErrorCode ierr;
  PetscInt       nsupplines;
  const PSLINE   *supplines;
  PSBUS          bus;
  PetscInt       i;
  PetscBool      flg=PETSC_FALSE;

  PetscFunctionBegin;
  if(ps->busext2intmap[fbus] != -1) {
    bus = &ps->bus[ps->busext2intmap[fbus]];
    ierr = PSBUSGetSupportingLines(bus,&nsupplines,&supplines);CHKERRQ(ierr);
    for(i=0; i < nsupplines; i++) {
      if(supplines[i]->fbus == fbus && supplines[i]->tbus == tbus) {
	ierr = PetscStrcmp(id,supplines[i]->ckt,&flg);CHKERRQ(ierr);
	if(flg) {
	  *line = supplines[i];
	  PetscFunctionReturn(0);
	}
      }
    }
  }
  *line = NULL;
  
  PetscFunctionReturn(0);
}

/*
  PSLINEGetVariableGlobalLocation - Gets the starting location for the variables for this line in the global vector

  Input Parameters:
. line - the line

  Output Parameter:
. locglob - the starting location for the variables in the global vector
*/
PetscErrorCode PSLINEGetVariableGlobalLocation(PSLINE line, PetscInt *locglob)
{

  PetscFunctionBegin;
  *locglob = line->startlocglob;
  PetscFunctionReturn(0);
}

/*
  PSLINEGetVariableLocation - Gets the starting location for the variables for this line in the local vector

  Input Parameters:
. line - the line object

  Output Parameter:
. loc - the starting location for the variables in the local vector
*/
PetscErrorCode PSLINEGetVariableLocation(PSLINE line, PetscInt *loc)
{

  PetscFunctionBegin;
  *loc = line->startloc;
  PetscFunctionReturn(0);
}

/*
  PSLINEGetConnectedBuses - Gets the buses connected to the line

  Input Parameters:
. line - the line

  Output Parameters:
. connbuses - connected buses

  Notes:
    connbuses is an PSBUS pointer object with connbuses[0] = From bus and
    connbuses[1] = To bus
*/
PetscErrorCode PSLINEGetConnectedBuses(PSLINE line,const PSBUS **connbuses)
{
  PetscFunctionBegin;
  *connbuses = line->connbuses;
  PetscFunctionReturn(0);
}

/*
  PSBUSGetSupportingLines - Gets the lines connected to this bus

  Input Parameters:
. bus - the bus

  Output Parameters:
+ nsupplines - number of support lines
. supplines - supporting lines

  Notes:
    supplines is a PSLINE pointer object with supplines[n] = nth connected line 
*/
PetscErrorCode PSBUSGetSupportingLines(PSBUS bus,PetscInt *nsupplines, const PSLINE **supplines)
{
  PetscFunctionBegin;
  *nsupplines = bus->nconnlines;
  *supplines = bus->connlines;
  PetscFunctionReturn(0);
}
  
/*
  PSBUSAddShunt - Adds shunt conductance and susceptance at this bus

  Input Parameters:
+ bus      - The PSBUS object
. Gs       - shunt conductance
- Bs       - shunt susceptance

  Notes: Gs and Bs should be in per unit
*/
PetscErrorCode PSBUSAddShunt(PSBUS bus,PetscScalar Gs,PetscScalar Bs)
{
  PetscFunctionBegin;
  bus->gl += Gs;
  bus->bl += Bs;
  PetscFunctionReturn(0);
}

/*
  PSBUSGetNGen - Gets the number of generators incident at the bus

  Input Parameters:
. bus  - the bus

  Output Parameters:
. ngen - number of generators incident at the bus
*/
PetscErrorCode PSBUSGetNGen(PSBUS bus,PetscInt *ngen)
{
  PetscFunctionBegin;
  *ngen = bus->ngen;
  PetscFunctionReturn(0);
}

/*
  PSBUSGetGen - Gets the generator incident at the bus

  Input Parameters:
+ bus  - the bus
- gen_num  - generator number (0, ..,ngen-1)

  Output Parameters:
. gen - the PSGEN generator object
*/
PetscErrorCode PSBUSGetGen(PSBUS bus,PetscInt gen_num,PSGEN *gen)
{

  PetscFunctionBegin;
  if(gen_num < 0 || gen_num > bus->ngen-1) SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,"Ngen at bus = %d, generator number requested %d",bus->ngen,gen_num);
  *gen = bus->gens[gen_num];
  PetscFunctionReturn(0);
}

/*
  PSBUSGetLoad - Gets the load incident at the bus

  Input Parameters:
+ bus  - the bus
- load_num  - load number (0, ..,nload-1)

  Output Parameters:
. load - the load PSLOAD object
*/
PetscErrorCode PSBUSGetLoad(PSBUS bus,PetscInt load_num,PSLOAD *load)
{

  PetscFunctionBegin;
  if(load_num < 0 || load_num > bus->nload-1) SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,"Nload at bus = %d, load number requested %d",bus->nload,load_num);
  *load = bus->loads[load_num];
  PetscFunctionReturn(0);
}

/*
  PSGetGen - Gets the generator object given the generator bus, and id

  Input Parameters
+ ps - the PS object
. gbus - from bus
- id   - generator id

  Output Parameters
. gen - the PSGEN object

*/
PetscErrorCode PSGetGen(PS ps,PetscInt gbus, const char* id,PSGEN *gen)
{
  PetscErrorCode ierr;
  PSBUS          bus;
  PetscInt       i;
  PSGEN          busgen;
  PetscBool      flg=PETSC_FALSE;

  PetscFunctionBegin;
  if(ps->busext2intmap[gbus] != -1) {
    bus = &ps->bus[ps->busext2intmap[gbus]];
    for(i=0; i < bus->ngen; i++) {
      ierr = PSBUSGetGen(bus,i,&busgen);CHKERRQ(ierr);
      ierr = PetscStrcmp(id,busgen->id,&flg);CHKERRQ(ierr);
      if(flg) {
	*gen = busgen;
	PetscFunctionReturn(0);
      }
    }
  }
  *gen = NULL;
  
  PetscFunctionReturn(0);
}

/*
  PSGENSetStatus - Sets the generator status (1 = ON, 0 = OFF)

  Inputs:
+ psgen - the PSGEN object
- status - the generator status
*/
PetscErrorCode PSGENSetStatus(PSGEN psgen,PetscInt status)
{
  PetscFunctionBegin;
  psgen->status = status;
  PetscFunctionReturn(0);
}

/*
  PSLOADSetStatus - Sets the load status (1 = ON, 0 = OFF)

  Inputs:
+ psload - the PSLOAD object
- status - the load status
*/
PetscErrorCode PSLOADSetStatus(PSLOAD psload,PetscInt status)
{
  PetscFunctionBegin;
  psload->status = status;
  PetscFunctionReturn(0);
}

/*
  PSLOADGetDynLoad - Returns the dynamic load object associated with the load

  Input Parameters:
. load - the PSLOAD object

  Output Parameters
. dynload - the dynamic load DYNLoadModel object
*/
PetscErrorCode PSLOADGetDYNLoad(PSLOAD load,DYNLoadModel *dynload)
{
  PetscFunctionBegin;
  *dynload = &load->dynload;
  PetscFunctionReturn(0);
}

/*
  PSGENGetDynGen - Returns the dynamic generator object associated with the generator

  Input Parameters:
. gen - the PSGEN object

  Output Parameters
. dyngen - the dynamic generator DYNGenModel object
*/
PetscErrorCode PSGENGetDYNGen(PSGEN gen,DYNGenModel *dyngen)
{
  PetscFunctionBegin;
  *dyngen = &gen->dyngen;
  PetscFunctionReturn(0);
}

/*
  PSGENGetDynExc - Returns the exciter model associated with the generator

  Input Parameters:
. gen - the PSGEN object

  Output Parameters
. dynexc - the exciter DYNExcModel object
*/
PetscErrorCode PSGENGetDYNExc(PSGEN gen,DYNExcModel *dynexc)
{
  PetscFunctionBegin;
  *dynexc = &gen->dynexc;
  PetscFunctionReturn(0);
}

/*
  PSGENGetDynTurbgov - Returns the turbgoviter model associated with the generator

  Input Parameters:
. gen - the PSGEN object

  Output Parameters
. dynturbgov - the turbgoviter DYNTurbgovModel object
*/
PetscErrorCode PSGENGetDYNTurbgov(PSGEN gen,DYNTurbgovModel *dynturbgov)
{
  PetscFunctionBegin;
  *dynturbgov = &gen->dynturbgov;
  PetscFunctionReturn(0);
}

/*
  PSGENGetDynStab - Returns the stabilizer model associated with the generator

  Input Parameters:
. gen - the PSGEN object

  Output Parameters
. dynstab - the DYNSTabModel object
*/
PetscErrorCode PSGENGetDYNStab(PSGEN gen,DYNStabModel *dynstab)
{
  PetscFunctionBegin;
  *dynstab = &gen->dynstab;
  PetscFunctionReturn(0);
}

/*
  PSBUSSetGenStatus - Sets the status of the generator
*/
PetscErrorCode PSBUSSetGenStatus(PSBUS bus,char gid[],PetscInt status)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PSGEN          gen;
  PetscBool      flg;
  PetscBool      statusupdate;

  PetscFunctionBegin;
  for(i=0; i < bus->ngen; i++) {
    ierr = PSBUSGetGen(bus,i,&gen);CHKERRQ(ierr);
    ierr = PetscStrcmp(gid,gen->id,&flg);CHKERRQ(ierr);
    if(flg) break;
  }
  
  if(status != gen->status) statusupdate = PETSC_TRUE;
  ierr = PSGENSetStatus(gen,status);CHKERRQ(ierr);
  if(statusupdate) {
    if(status) {
      bus->ngenON++;
    } else {
      bus->ngenON--;
      /* Change the bus type to PQ if no generators are active at this bus */
      if(!bus->ngenON) bus->ide = PQ_BUS;
      gen->pg = gen->qg = 0.0;
    }
  }
  PetscFunctionReturn(0);
}

/*
  PSBUSGetVariableGlobalLocation - Gets the starting location for the variables for this bus in the global vector

  Input Parameters:
. bus - the bus

  Output Parameter:
. locglob - the starting location for the variables in the global vector
*/
PetscErrorCode PSBUSGetVariableGlobalLocation(PSBUS bus, PetscInt *locglob)
{

  PetscFunctionBegin;
  *locglob = bus->startlocglob;
  PetscFunctionReturn(0);
}

/*
  PSBUSGetVariableLocation - Gets the starting location for the variables for this bus in the local vector

  Input Parameters:
. bus - the bus

  Output Parameter:
. loc - the starting location for the variables in the local vector
*/
PetscErrorCode PSBUSGetVariableLocation(PSBUS bus, PetscInt *loc)
{

  PetscFunctionBegin;
  *loc = bus->startloc;
  PetscFunctionReturn(0);
}

/*
  PSBUSIsGhosted - Returns true if the bus is a ghost bus

  Input Parameters
. PSBUS - The PSBUS object

  Output Parameters
. ghostbus - TRUE if the bus is a ghost bus
*/
PetscErrorCode PSBUSIsGhosted(PSBUS bus,PetscBool *ghostbus)
{

  PetscFunctionBegin;
  *ghostbus = bus->isghost;
  PetscFunctionReturn(0);
}

/*
  PSGENDestroy - Destroys the GEN data struct in PS

  Input Parameters
. ps - the PS object
*/
PetscErrorCode PSGENDestroy(PS ps)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  if(ps-> app == APP_DYNSIM) {
    PSGEN Gen;
    PetscInt i;

    for(i=0;i < ps->ngen; i++) {
      Gen = &ps->gen[i];
      if(Gen->initial_status) {
	ierr = DYNGenModelDestroy(&Gen->dyngen);CHKERRQ(ierr);
	if(Gen->hasexc) {
	  ierr = DYNExcModelDestroy(&Gen->dynexc);CHKERRQ(ierr);
	}
	if(Gen->hasturbgov) {
	  ierr = DYNTurbgovModelDestroy(&Gen->dynturbgov);CHKERRQ(ierr);
	}
	if(Gen->hasstab) {
	  ierr = DYNStabModelDestroy(&Gen->dynstab);CHKERRQ(ierr);
	}
      }
    }
  }
  ierr = PetscFree(ps->gen);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  PSLOADDestroy - Destroys the LOAD data struct in PS

  Input Parameters
. ps - the PS object
*/
PetscErrorCode PSLOADDestroy(PS ps)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  if(ps-> app == APP_DYNSIM) {
    PSLOAD Load;
    PetscInt i;

    for(i=0;i < ps->nload; i++) {
      Load = &ps->load[i];
      if(Load->status) {
	ierr = DYNLoadModelDestroy(&Load->dynload);CHKERRQ(ierr);
      }
    }
  }
  ierr = PetscFree(ps->load);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
    
/*
  PSSetApplication - Sets the type of application to be run on the PS object

  Input Parameters:
+ PS - the PS object
- psapp - the application (DYNSIM,ACPF,DCPF)
*/
PetscErrorCode PSSetApplication(PS ps,PSApp psapp)
{
  PetscFunctionBegin;
  ps->app = psapp;
  PetscFunctionReturn(0);
}

/*
  PSIncreaseReferenceCount - Increases the reference count of PS to indicate that is being shared
                             by other objects

  Input Parameter
. ps - The PS object
*/
PetscErrorCode PSIncreaseReferenceCount(PS ps)
{
  PetscFunctionBegin;
  ps->refct++;
  PetscFunctionReturn(0);
}

/*
  PSDecreaseReferenceCount - Decreases the reference count of PS

  Input Parameter
. ps - The PS object
*/
PetscErrorCode PSDecreaseReferenceCount(PS ps)
{
  PetscFunctionBegin;
  ps->refct--;
  PetscFunctionReturn(0);
}

/*
  PSCreate - Creates the PS object

  Input Parameter
. mpicomm - the MPI communicator

  Output Parameter
. psout - The PS object
*/
PetscErrorCode PSCreate(MPI_Comm mpicomm,PS *psout)
{
  PetscErrorCode ierr;
  PS             ps;

  PetscFunctionBegin;
  ierr = PetscCalloc1(1,&ps);CHKERRQ(ierr);

  /* Create communicator */
  ierr = COMMCreate(mpicomm,&ps->comm);CHKERRQ(ierr);

  ps->MVAbase = 0.0;
  ps->Nbus  = -1;
  ps->Ngen  = -1;
  ps->Nbranch = -1;
  ps->Nload   = -1;
  ps->refct   = 0;
  ps->app     = APP_NONE;
  ps->ndiff   = 0;
  ps->nconncomp = 0;
 
  ierr = PSIncreaseReferenceCount(ps);CHKERRQ(ierr);

  ps->setupcalled = PETSC_FALSE;

  *psout = ps;
  PetscFunctionReturn(0);
}

/*
  PSDestroy - Destroys the PS object created with PSCreate

  Input Parameter
  ps - pointer to PS object to be destroyed
*/
PetscErrorCode PSDestroy(PS *ps)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if(!(*ps)) PetscFunctionReturn(0);

  ierr = PSDecreaseReferenceCount(*ps);CHKERRQ(ierr);
  if((*ps)->refct > 0) {*ps = 0; PetscFunctionReturn(0);}

  /* Destroy PS communicator */
  ierr = COMMDestroy(&(*ps)->comm);CHKERRQ(ierr);

  /* Destroy PS bus,branch,gen,load objects */
  ierr = PetscFree((*ps)->bus);CHKERRQ(ierr);
  ierr = PSGENDestroy((*ps));CHKERRQ(ierr);
  ierr = PSLOADDestroy((*ps));CHKERRQ(ierr);
  ierr = PetscFree((*ps)->load);CHKERRQ(ierr);
  ierr = PetscFree((*ps)->line);CHKERRQ(ierr);

  ierr = PSConnCompDestroy(*ps);CHKERRQ(ierr);

  ierr = PetscFree((*ps)->busext2intmap);CHKERRQ(ierr);
  ierr = DMDestroy(&(*ps)->networkdm);CHKERRQ(ierr);
  ierr = PetscFree(*ps);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  PSReadPSSERawData - Reads the PSSE raw data file and populates the PS object

  Input Parameter
+ ps      - The power system object ps
- netfile - Name of the power system network data file in PSSE raw data format.

*/
PetscErrorCode PSReadPSSERawData(PS ps,const char netfile[])
{
  FILE           *fp;
  PetscErrorCode ierr;
  char           *out;

  PetscFunctionBegin;


  char line[MAXLINE];
  int dataformatflag = -1; 
  /*
  0 = BUS DATA, 1 = LOAD DATA, 2 = SHUNT DATA, 3 = GENERATOR DATA, 4 = BRANCH DATA, 5 = TRANSFORMER DATA, 6 = AREA DATA
  7 = TWO-TERMINAL DC DATA, 8 = VSC DC LINE DATA, 9 = IMPEDANCE CORRECTION DATA, 10 = MULTI-TERMINAL DC DATA, 11 = MULTI-SECTION LINE DATA
  12 = ZONE DATA, 13 = TRANSFER DATA, 14 = OWNER DATA, 15 = FACTS DEVICE DATA, 16 = SWITCHED SHUNT DATA, 17 = GNE DATA  
  */
  // Network Size variables
  PetscInt Nbus = 0;
  PetscInt Nload = 0;
  PetscInt Ngenerator = 0;
  PetscInt Nbranch = 0;
  PetscInt Ntransformer = 0;
  PetscInt Nshunt = 0;
  // Case Identification Data
  PetscInt IC = 0;
  PetscScalar SBASE = 100.0;
  PetscInt REV = 33;
  PetscInt XFRRAT = 0;
  PetscInt NXFRAT = 0;
  PetscScalar BASFRQ = 0.0;
  // Structure Pointers in PS
  PSBUS          Bus;
  PSLOAD         Load;
  PSGEN          Gen;
  PSLINE         Branch;
  //Unused Transformer data
  PetscInt K, CW, CZ, CM, NMETR, COD1, CONT1, NTP1, TAB1;
  PetscScalar MAG1, MAG2, NOMV1, RMA1, RMI1, VMA1, VMI1, CR1, CX1, CNXA1, WINDV2, NOMV2;
  char transname[20];
  // Temp Variables
  PetscInt loadi=0,geni=0,bri=0,busi=0,i=0;
  PetscInt MET = 1;
  PetscInt internalindex = 0;
  PetscScalar R,X,Bc,B,G,Zm,tap,shift,tap2,tapr,tapi;
  PetscInt tempbusi = 0, maxbusi = -1;
  PetscInt shuntbus;
  char     shuntid[10];
  PetscInt shuntstatus;
  PetscScalar gshunt,bshunt;

  if(ps->comm->type != PETSC_COMM_SELF && ps->comm->rank != 0) {
    ps->Nbranch = ps->Nbus = ps->Ngen = ps->Nload = 0;
    PetscFunctionReturn(0);
  }

  fp = fopen(netfile,"r");
  /* Check for valid file */
  if (fp == NULL) {
    SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Cannot open file %s",netfile);CHKERRQ(ierr);
  }

  /* Copy the network file name */
  ierr = PetscStrcpy(ps->net_file_name,netfile);CHKERRQ(ierr);

  //Read Case Identification Data
  out = fgets(line,MAXLINE,fp);
  if(out != NULL){
    sscanf(line,"%d, %lf, %d, %d, %d, %lf",&IC, &SBASE, &REV, &XFRRAT, &NXFRAT, &BASFRQ);
    out = fgets(line,MAXLINE,fp);
    if(out == NULL)  PetscFunctionReturn(0);// for commas
    out = fgets(line,MAXLINE,fp);
    if(out == NULL)  PetscFunctionReturn(0);// for commas
    dataformatflag++;
  }else{
    PetscFunctionReturn(0);
  }
  /* Setting ps->sbase to 100 */
  ps->MVAbase = SBASE;

  while((out=fgets(line,MAXLINE,fp)) != NULL) {
    if(strstr(line,"0 /") != NULL) {
      dataformatflag++;
      continue;
    }

    switch(dataformatflag){
      case 0:
        sscanf(line,"%d", &tempbusi);
        if(tempbusi > maxbusi) maxbusi = tempbusi;
        Nbus++;
      break;
      case 1:
        Nload++;
      break;
      case 2:
	Nshunt++;
      break;
      case 3:
        Ngenerator++;
      break;
      case 4:
        Nbranch++;
      break;
      case 5:
        out = fgets(line,MAXLINE,fp);
        out = fgets(line,MAXLINE,fp);
        out = fgets(line,MAXLINE,fp);
        Ntransformer++;
      break;
      default:
        //Todo: impliment parser for other data
      break;
    }
  }
  fclose(fp);
  
  ps->Nbus    = ps->nbus    = Nbus;
  ps->Nload   = ps->nload = Nload;
  ps->Ngen    = ps->ngen    = Ngenerator;
  ps->Nbranch = ps->nbranch = Nbranch + Ntransformer;
 
#if defined DEBUGPS
  ierr = PetscPrintf(PETSC_COMM_SELF,"System summary : Nbus = %d, Nload = %d, Ngenerator = %d, Nbranch = %d\n",ps->Nbus,ps->Ngen,ps->Nload,ps->Nbranch);CHKERRQ(ierr);
#endif
  ierr = PetscCalloc1(ps->Nbus,&ps->bus);CHKERRQ(ierr);
  ierr = PetscCalloc1(ps->Ngen,&ps->gen);CHKERRQ(ierr);
  ierr = PetscCalloc1(ps->Nload,&ps->load);CHKERRQ(ierr);
  ierr = PetscCalloc1(ps->Nbranch,&ps->line);CHKERRQ(ierr);
  Bus = ps->bus; Gen = ps->gen; Load = ps->load; Branch = ps->line;

  //Set initial and default data for bus
  for(i=0; i < ps->Nbus; i++) {
    ps->bus[i].ngen = ps->bus[i].nload = ps->bus[i].ngenON = 0;
    ps->bus[i].qrange = ps->bus[i].qmintot = ps->bus[i].Pgtot = ps->bus[i].MVAbasetot =  0.0 ;
    ps->bus[i].Vmax = 1.1;
    ps->bus[i].Vmin = 0.9;
  }
  /* Allocate external to internal bus number mapping array */

  PetscInt *busext2intmap;
  ps->maxbusnum = maxbusi;
  ierr = PetscCalloc1(ps->maxbusnum+1,&ps->busext2intmap);CHKERRQ(ierr);
  busext2intmap = ps->busext2intmap;
  for(i=0; i < ps->maxbusnum+1; i++) busext2intmap[i] = -1;

  dataformatflag = -1;
  fp = fopen(netfile,"r");
  if((out = fgets(line,MAXLINE,fp)) != NULL){//for case identification
    if((out = fgets(line,MAXLINE,fp)) == NULL)  PetscFunctionReturn(0);// for commas
    if((out = fgets(line,MAXLINE,fp)) == NULL)  PetscFunctionReturn(0);// for commas
    dataformatflag++;
  }else{
    PetscFunctionReturn(0);
  }
  while((out = fgets(line,MAXLINE,fp)) != NULL) {
    if(strstr(line,"0 /") != NULL) {
      dataformatflag++;
      continue;
    }
    switch(dataformatflag){
      case 0:
        sscanf(line,"%d, '%[^\t\']', %lf, %d, %d, %d, %d, %lf, %lf",&Bus[busi].bus_i, Bus[busi].name, &Bus[busi].basekV, &Bus[busi].ide, &Bus[busi].area, &Bus[busi].zone, &Bus[busi].owner, &Bus[busi].vm, &Bus[busi].va);
#if defined DEBUGPS
        if(busi < 2) ierr = PetscPrintf(PETSC_COMM_SELF,"BUSData[%d] : %d, '%s', %lf, %d, %d, %d, %d, %lf, %lf\n", busi, Bus[busi].bus_i, Bus[busi].name, Bus[busi].basekV, Bus[busi].ide, Bus[busi].area, Bus[busi].zone, Bus[busi].owner, Bus[busi].vm, Bus[busi].va);CHKERRQ(ierr);
#endif               
        busext2intmap[Bus[busi].bus_i] = busi;
        Bus[busi].internal_i = busi;
        Bus[busi].nload = 0;
        Bus[busi].ngen = 0;
        Bus[busi].ngenON = 0;
	Bus[busi].nshunt = 0;
        Bus[busi].nvhi = 1.1;
        Bus[busi].nvlo = 0.9;
        Bus[busi].evhi = 1.1;
        Bus[busi].evlo = 0.9;
        Bus[busi].Vmin = 1.1;
        Bus[busi].Vmax = 0.9;
        Bus[busi].gl = 0;
        Bus[busi].bl = 0;        
        busi++;
      break;
      case 1:
        sscanf(line,"%d, '%[^\t\']', %d, %d, %d, %lf, %lf, %lf, %lf, %lf, %lf, %d",&Load[loadi].bus_i, Load[loadi].id, &Load[loadi].status, &Load[loadi].area, &Load[loadi].zone, &Load[loadi].pl, &Load[loadi].ql, &Load[loadi].ip, &Load[loadi].iq, &Load[loadi].yp, &Load[loadi].yq, &Load[loadi].owner);
#if defined DEBUGPS
        if(loadi < 2) ierr = PetscPrintf(PETSC_COMM_SELF,"LOADData[%d] : %d, '%s', %d, %d, %d, %lf, %lf, %lf, %lf, %lf, %lf, %d\n",loadi, Load[loadi].bus_i, Load[loadi].id, Load[loadi].status, Load[loadi].area, Load[loadi].zone, Load[loadi].pl, Load[loadi].ql, Load[loadi].ip, Load[loadi].iq, Load[loadi].yp, Load[loadi].yq, Load[loadi].owner);CHKERRQ(ierr);
#endif
        Load[loadi].scale = 1;
        Load[loadi].intrpt = 0;
	Load[loadi].dynloadsetup = 0;
      	Load[loadi].pl /= ps->MVAbase;
      	Load[loadi].ql /= ps->MVAbase;
      	Load[loadi].ip /= ps->MVAbase;
      	Load[loadi].iq /= ps->MVAbase;
      	Load[loadi].yp /= ps->MVAbase;
      	Load[loadi].yq /= ps->MVAbase;
        internalindex = busext2intmap[Load[loadi].bus_i];
        Load[loadi].internal_i = internalindex;
        Bus[internalindex].lidx[Bus[internalindex].nload] = loadi;
        Bus[internalindex].nload++;
        loadi++;
      break;
      case 2:
	sscanf(line,"%d, '%[^\t\']', %d, %lf, %lf",&shuntbus,shuntid,&shuntstatus,&gshunt,&bshunt);
	if(shuntstatus) {
	  internalindex = busext2intmap[shuntbus];
	  if(Bus[internalindex].nshunt == 1) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Bus %d: No support for more than 1 fixed shunt at bus",Bus[internalindex].bus_i);
	  Bus[internalindex].nshunt++;
	  Bus[internalindex].gl = gshunt/ps->MVAbase;
	  Bus[internalindex].bl = bshunt/ps->MVAbase;
	}
      break;
      case 3:
        sscanf(line,"%d, '%[^\t\']', %lf, %lf, %lf, %lf, %lf, %d, %lf, %lf, %lf, %lf, %lf, %lf, %d, %lf, %lf, %lf, %d, %lf",&Gen[geni].bus_i, Gen[geni].id, &Gen[geni].pg, &Gen[geni].qg, &Gen[geni].qt, &Gen[geni].qb, &Gen[geni].vs, &Gen[geni].ireg, &Gen[geni].mbase, &Gen[geni].zr, &Gen[geni].zx, &Gen[geni].rt, &Gen[geni].xt, &Gen[geni].gtap, &Gen[geni].status, &Gen[geni].rmpct, &Gen[geni].pt, &Gen[geni].pb, &Gen[geni].o1, &Gen[geni].f1);
#if defined DEBUGPS
        if(geni < 2) ierr = PetscPrintf(PETSC_COMM_SELF,"GENERATORData[%d] : %d, '%s', %lf, %lf, %lf, %lf, %lf, %d, %lf, %lf, %lf, %lf, %lf, %lf, %d, %lf, %lf, %lf, %d, %lf\n",geni, Gen[geni].bus_i, Gen[geni].id, Gen[geni].pg, Gen[geni].qg, Gen[geni].qt, Gen[geni].qb, Gen[geni].vs, Gen[geni].ireg, Gen[geni].mbase, Gen[geni].zr, Gen[geni].zx, Gen[geni].rt, Gen[geni].xt, Gen[geni].gtap, Gen[geni].status, Gen[geni].rmpct, Gen[geni].pt, Gen[geni].pb, Gen[geni].o1, Gen[geni].f1);CHKERRQ(ierr);
#endif
	Gen[geni].initial_status = Gen[geni].status;
        internalindex = busext2intmap[Gen[geni].bus_i];
        Gen[geni].internal_i = internalindex;
        Bus[internalindex].gidx[Bus[internalindex].ngen] = geni;
        if(Gen[geni].status == 1) {
      	  Gen[geni].pg /= ps->MVAbase;
      	  Gen[geni].qg /= ps->MVAbase;
	  Gen[geni].qt = Gen[geni].qt/ps->MVAbase;
	  Gen[geni].qb = Gen[geni].qb/ps->MVAbase;
	  Gen[geni].dyngensetup = 0;
	  Bus[internalindex].qrange += (Gen[geni].qt - Gen[geni].qb);
	  Bus[internalindex].qmintot += Gen[geni].qb;
	  Bus[internalindex].Pgtot += PetscAbsScalar(Gen[geni].pg);
          Bus[internalindex].MVAbasetot += Gen[geni].mbase;
      	  Bus[internalindex].ngenON++;
      	}
        Bus[internalindex].ngen++;
        geni++;
      break;
      case 4:
        sscanf(line,"%d, %d, '%[^\t\']', %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %d, %d, %lf, %d, %lf",&Branch[bri].fbus, &Branch[bri].tbus, Branch[bri].ckt, &Branch[bri].r, &Branch[bri].x, &Branch[bri].b, &Branch[bri].rateA, &Branch[bri].rateB, &Branch[bri].rateC, &Branch[bri].gi, &Branch[bri].bi, &Branch[bri].gj, &Branch[bri].bj, &Branch[bri].status, &MET, &Branch[bri].length, &Branch[bri].o1, &Branch[bri].f1);
#if defined DEBUGPS
        if(bri < 2) ierr = PetscPrintf(PETSC_COMM_SELF,"BRANCHData[%d] : %d, %d, '%s', %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %d, %d, %lf, %d, %lf\n",bri, Branch[bri].fbus, Branch[bri].tbus, Branch[bri].ckt, Branch[bri].r, Branch[bri].x, Branch[bri].b, Branch[bri].rateA, Branch[bri].rateB, Branch[bri].rateC, Branch[bri].gi, Branch[bri].bi, Branch[bri].gj, Branch[bri].bj, Branch[bri].status, MET, Branch[bri].length, Branch[bri].o1, Branch[bri].f1);CHKERRQ(ierr);
#endif 
        Branch[bri].met = 1;
        internalindex = busext2intmap[Branch[bri].fbus];
        Branch[bri].internal_i = internalindex;        
        internalindex = busext2intmap[Branch[bri].tbus];
        Branch[bri].internal_j = internalindex;
        Branch[bri].tapratio = 1.0;
        Branch[bri].phaseshift = 0.0;         
        R = Branch[bri].r;
        X = Branch[bri].x;
        Bc = Branch[bri].b;

        Zm = R*R + X*X;
        G  = R/Zm;
        B  = -X/Zm;

        tap = Branch[bri].tapratio;
        shift = Branch[bri].phaseshift;
        tap2 = tap*tap;
        tapr = tap*cos(shift);
        tapi = tap*sin(shift);

        Branch[bri].yff[0] = G/tap2; 
        Branch[bri].yff[1] = (B+Bc/2.0)/tap2;
        
        Branch[bri].yft[0] = -(G*tapr - B*tapi)/tap2;
        Branch[bri].yft[1] = -(B*tapr + G*tapi)/tap2;

        Branch[bri].ytf[0] = -(G*tapr + B*tapi)/tap2;
        Branch[bri].ytf[1] = -(B*tapr - G*tapi)/tap2;

        Branch[bri].ytt[0] = G;
        Branch[bri].ytt[1] = B+Bc/2.0; 


        bri++;
      break;
      case 5:
        sscanf(line,"%d, %d, %d, '%[^\']', %d, %d, %d, %lf, %lf, %d, '%[^\t\']', %d, %d, %lf",&Branch[bri].fbus, &Branch[bri].tbus, &K, Branch[bri].ckt, &CW, &CZ, &CM, &MAG1, &MAG2, &NMETR, transname, &Branch[bri].status, &Branch[bri].o1, &Branch[bri].f1);
	if(K != 0) { /* Three winding transformer */
	  SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Found three winding transformer in the transformer data.\n\
                Three winding transformers are currently not supported.\n\
                Convert the raw file to MATPOWER format");
	}
        out = fgets(line,MAXLINE,fp);
        sscanf(line,"%lf, %lf, %lf",&Branch[bri].r, &Branch[bri].x, &Branch[bri].sbase12);
        
        out = fgets(line,MAXLINE,fp);
        sscanf(line,"%lf, %lf, %lf, %lf, %lf, %lf, %d, %d, %lf, %lf, %lf, %lf, %d, %d, %lf, %lf, %lf",&Branch[bri].tapratio, &NOMV1, &Branch[bri].phaseshift, &Branch[bri].rateA, &Branch[bri].rateB, &Branch[bri].rateC, &COD1, &CONT1, &RMA1, &RMI1, &VMA1, &VMI1, &NTP1, &TAB1, &CR1, &CX1, &CNXA1);
        
        out = fgets(line,MAXLINE,fp);
        sscanf(line,"%lf %lf",&WINDV2, &NOMV2);
#if defined DEBUGPS
        if(bri - Nbranch < 2) ierr = PetscPrintf(PETSC_COMM_SELF,"TRANSFORMERData[%d] : %d, %d, %s, %lf, %lf, %lf, %lf, %lf, %lf, %d, %d, %lf, %lf, %lf\n",bri-Nbranch, Branch[bri].fbus, Branch[bri].tbus, Branch[bri].ckt, Branch[bri].r, Branch[bri].x, Branch[bri].b, Branch[bri].rateA, Branch[bri].rateB, Branch[bri].rateC, Branch[bri].status, Branch[bri].o1, Branch[bri].f1, Branch[bri].tapratio, Branch[bri].phaseshift);CHKERRQ(ierr);
#endif  
        internalindex = busext2intmap[Branch[bri].fbus];
        Branch[bri].internal_i = internalindex;        
        internalindex = busext2intmap[Branch[bri].tbus];
        Branch[bri].internal_j = internalindex;

        Branch[bri].b = 0;
        
        R = Branch[bri].r;
        X = Branch[bri].x;
        Bc = Branch[bri].b;

        Zm = R*R + X*X;
        G  = R/Zm;
        B  = -X/Zm;

        tap = Branch[bri].tapratio;
        shift = Branch[bri].phaseshift;
        tap2 = tap*tap;
        tapr = tap*cos(shift);
        tapi = tap*sin(shift);

        Branch[bri].yff[0] = G/tap2; 
        Branch[bri].yff[1] = (B+Bc/2.0)/tap2;
        
        Branch[bri].yft[0] = -(G*tapr - B*tapi)/tap2;
        Branch[bri].yft[1] = -(B*tapr + G*tapi)/tap2;

        Branch[bri].ytf[0] = -(G*tapr + B*tapi)/tap2;
        Branch[bri].ytf[1] = -(B*tapr - G*tapi)/tap2;

        Branch[bri].ytt[0] = G;
        Branch[bri].ytt[1] = B+Bc/2.0; 
           
        bri++;
      break;
      default:
        //Todo: impliment parser for other data
      break;
    }
  }
  fclose(fp);



  PetscFunctionReturn(0);
}

/*
  PSReadMatPowerData - Reads the MATPOWER data file and populates the PS object

  Input Parameter
+ ps      - The power system object ps
- netfile - Name of the power system network data file in MATPOWER data format.

*/
PetscErrorCode PSReadMatPowerData(PS ps,const char netfile[])
{
  FILE           *fp;
  PetscErrorCode ierr;
  PSBUS          Bus;
  PSLOAD         Load;
  PSGEN          Gen;
  PSLINE         Branch;
  PetscInt       line_counter=0,linenum;
  PetscInt       bus_start_line=-1,bus_end_line=-1; /* xx_end_line points to the next line after the record ends */
  PetscInt       gen_start_line=-1,gen_end_line=-1;
  PetscInt       br_start_line=-1,br_end_line=-1;
  PetscInt       gencost_start_line=-1,gencost_end_line=-1;
  PetscInt       bus_nblank_lines=0, gen_nblank_lines=0, br_nblank_lines=0,gencost_nblank_lines=0; /* Number of blank lines in bus, gen, and gencost branch arrays */
  char           line[MAXLINE];
  PetscInt       loadi=0,geni=0,bri=0,busi=0,gencosti=0,i;
  PetscInt       extbusnum,bustype_i;
  PetscScalar    Pd,Qd;
  PetscInt       intbusnum;
  char           *str;
  char           *out;

  PetscFunctionBegin;

  if(ps->comm->type != PETSC_COMM_SELF && ps->comm->rank != 0) {
    ps->Nbranch = ps->Nbus = ps->Ngen = ps->Nload = 0;
    PetscFunctionReturn(0);
  }

  fp = fopen(netfile,"r");
  /* Check for valid file */
  if (fp == NULL) {
    SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Cannot open file %s",netfile);CHKERRQ(ierr);
  }

  /* Copy the network file name */
  ierr = PetscStrcpy(ps->net_file_name,netfile);CHKERRQ(ierr);

  ps->Nload=0;
  ps->maxbusnum = -1;
  while((out = fgets(line,MAXLINE,fp)) != NULL) {
    if(strstr(line,"mpc.baseMVA")) {
      /* Read base MVA */
      str = strtok(line," =;");
      str = strtok(NULL," =;");
      sscanf(str,"%lf",&ps->MVAbase);
    }
    if(strstr(line,"mpc.bus") != NULL && bus_start_line == -1)    bus_start_line = line_counter+1; /* Bus data starts from next line */
    if(strstr(line,"mpc.gen") != NULL && gen_start_line == -1)    gen_start_line = line_counter+1; /* Generator data starts from next line */
    if(strstr(line,"mpc.branch") != NULL) br_start_line = line_counter+1; /* Branch data starts from next line */
    if(strstr(line,"mpc.gencost") != NULL) gencost_start_line = line_counter+1; /* Gen cost data starts from next line */
    if(strstr(line,"];") != NULL) {
      if (bus_start_line != -1 && bus_end_line == -1) bus_end_line = line_counter;
      if (gen_start_line != -1 && gen_end_line == -1) gen_end_line = line_counter;
      if (br_start_line  != -1 && br_end_line == -1) br_end_line = line_counter;
      if (gencost_start_line != -1 && gencost_end_line == -1) gencost_end_line = line_counter;
    }

    if (bus_start_line != -1 && bus_end_line == -1) {
      if(strcmp(line,"\n")==0 || strcmp(line,"\r\n")==0) bus_nblank_lines++;
    }
    if (gen_start_line != -1 && gen_end_line == -1) {
      if(strcmp(line,"\n")==0 || strcmp(line,"\r\n")==0) gen_nblank_lines++;
    }

    if (br_start_line  != -1 && br_end_line == -1) {
      if(strcmp(line,"\n")==0 || strcmp(line,"\r\n")==0) br_nblank_lines++;
    }

    if (gencost_start_line != -1 && gencost_end_line == -1) {
      if(strcmp(line,"\n")==0 || strcmp(line,"\r\n")==0) gencost_nblank_lines++;
    }

    /* Count the number of pq loads */
    if(bus_start_line != -1 && line_counter >= bus_start_line && bus_end_line == -1) {
      sscanf(line,"%d %d %lf %lf",&extbusnum,&bustype_i,&Pd,&Qd);
      if(!((Pd == 0.0) && (Qd == 0.0))) ps->Nload++;
      if (extbusnum > ps->maxbusnum) ps->maxbusnum = extbusnum;
    }
    line_counter++;
  }
  fclose(fp);

  ps->Nbus    = ps->nbus    = bus_end_line - bus_start_line - bus_nblank_lines;
  ps->Ngen    = ps->ngen    = gen_end_line - gen_start_line - gen_nblank_lines;
  ps->Nbranch = ps->nbranch = br_end_line  - br_start_line - br_nblank_lines;
  ps->nload = ps->Nload;
#if defined DEBUGPS
  ierr = PetscPrintf(PETSC_COMM_SELF,"System summary : Nbuses = %d, Ngen = %d, Nload = %d, Nbranch = %d\n",ps->Nbus,ps->Ngen,ps->Nload,ps->Nbranch);CHKERRQ(ierr);
#endif
  ierr = PetscCalloc1(ps->Nbus,&ps->bus);CHKERRQ(ierr);
  ierr = PetscCalloc1(ps->Ngen,&ps->gen);CHKERRQ(ierr);
  ierr = PetscCalloc1(ps->Nload,&ps->load);CHKERRQ(ierr);
  ierr = PetscCalloc1(ps->Nbranch,&ps->line);CHKERRQ(ierr);
  Bus = ps->bus; Gen = ps->gen; Load = ps->load; Branch = ps->line;

  for(i=0; i < ps->Nbus; i++) {
    ps->bus[i].ngen = ps->bus[i].nload = ps->bus[i].ngenON = ps->bus[i].nshunt = 0;
    ps->bus[i].qrange = ps->bus[i].qmintot = ps->bus[i].Pgtot = ps->bus[i].MVAbasetot = 0.0;
  }

  /*  ierr = PetscPrintf(PETSC_COMM_SELF,"Rank[%d] maxbusnum %d",ps->comm->rank,ps->maxbusnum);CHKERRQ(ierr); */

  /* Allocate external to internal bus number mapping array */
  PetscInt *busext2intmap;
  ierr = PetscCalloc1(ps->maxbusnum+1,&ps->busext2intmap);CHKERRQ(ierr);
  busext2intmap = ps->busext2intmap;
  for(i=0; i < ps->maxbusnum+1; i++) busext2intmap[i] = -1;

  fp = fopen(netfile,"r");
  /* Reading data */
  for(i=0;i<line_counter;i++) {
    out = fgets(line,MAXLINE,fp);

    if((i >= bus_start_line) && (i < bus_end_line)) {
      if(strcmp(line,"\n")==0 || strcmp(line,"\r\n")==0) continue;
      /* Bus data */
      sscanf(line,"%d %d %lf %lf %lf %lf %d %lf %lf %lf %d %lf %lf",		\
	     &Bus[busi].bus_i,&Bus[busi].ide,&Pd,&Qd,&Bus[busi].gl,	\
	     &Bus[busi].bl,&Bus[busi].area,&Bus[busi].vm,&Bus[busi].va,&Bus[busi].basekV,&Bus[busi].zone,\
	     &Bus[busi].Vmax,&Bus[busi].Vmin);
      Bus[busi].internal_i = busi;
      busext2intmap[Bus[busi].bus_i] = busi;
      /* Convert bl and gl to per unit */
      Bus[busi].bl = Bus[busi].bl/ps->MVAbase;
      Bus[busi].gl = Bus[busi].gl/ps->MVAbase;
      Bus[busi].nshunt++;

      if(!((Pd == 0.0) && (Qd == 0.0))) {
	Load[loadi].bus_i = Bus[busi].bus_i;
	Load[loadi].status = 1;
	Load[loadi].pl = Pd/ps->MVAbase;
	Load[loadi].ql = Qd/ps->MVAbase;
	Load[loadi].area = Bus[busi].area;
	Load[loadi].internal_i = busi;

	intbusnum = busext2intmap[Load[loadi].bus_i];

	/* MatPower does not have ids for loads. Using Bus[i].nload as the id */
	snprintf(Load[loadi].id,3,"%-2d",1+Bus[intbusnum].nload);

	Bus[busi].lidx[Bus[busi].nload++] = loadi;
	if (Bus[busi].nload > NLOAD_AT_BUS_MAX)
	  SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Exceeded maximum number of loads allowed at bus");
	loadi++;
      }
      busi++;
    }

    /* Read generator data */
    if(i >= gen_start_line && i < gen_end_line) {
      if(strcmp(line,"\n")==0 || strcmp(line,"\r\n")==0) continue;
      sscanf(line,"%d %lf %lf %lf %lf %lf %lf %d %lf %lf",&Gen[geni].bus_i, \
	     &Gen[geni].pg,&Gen[geni].qg,&Gen[geni].qt,&Gen[geni].qb, \
	     &Gen[geni].vs,&Gen[geni].mbase,&Gen[geni].status,&Gen[geni].pt, \
	     &Gen[geni].pb);

      intbusnum = busext2intmap[Gen[geni].bus_i];

      Gen[geni].initial_status = Gen[geni].status;
      /* Convert Pg and Qg to per unit */
      if(Gen[geni].status) {
	Gen[geni].pg = Gen[geni].pg/ps->MVAbase;
	Gen[geni].qg = Gen[geni].qg/ps->MVAbase;
	Gen[geni].qt = Gen[geni].qt/ps->MVAbase;
	Gen[geni].qb = Gen[geni].qb/ps->MVAbase;
	Bus[intbusnum].qrange += (Gen[geni].qt - Gen[geni].qb);
	Bus[intbusnum].qmintot += Gen[geni].qb;
	Bus[intbusnum].Pgtot += PetscAbsScalar(Gen[geni].pg);
	Bus[intbusnum].MVAbasetot += Gen[geni].mbase;
	Bus[intbusnum].ngenON++;
      } else {
	Gen[geni].pg = Gen[geni].qg = 0.0;
      }

      Gen[geni].internal_i = intbusnum;

      /* MatPower does not have ids for generators. Using Bus[i].ngen as the id */
      snprintf(Gen[geni].id,3,"%-2d",1+Bus[intbusnum].ngen);

      Bus[intbusnum].gidx[Bus[intbusnum].ngen++] = geni;

      //      Bus[intbusnum].vm = Gen[geni].vs;

#if defined DEBUGPS
      ierr = PetscPrintf(PETSC_COMM_SELF,"%d %d %d %s\n",Gen[geni].status, intbusnum,Bus[intbusnum].ngen,line);CHKERRQ(ierr);
#endif
      if (Bus[intbusnum].ngen > NGEN_AT_BUS_MAX) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Exceeded maximum number of generators allowed at bus");
      geni++;
    }

    //    ierr = PetscPrintf(PETSC_COMM_SELF,"Came here\n");CHKERRQ(ierr);

    /* Read generator cost data */
    if(i >= gencost_start_line && i < gencost_end_line) {
      if(strcmp(line,"\n")==0 || strcmp(line,"\r\n")==0) continue;
      sscanf(line,"%d %lf %lf %d %lf %lf %lf",&Gen[gencosti].cost_model,	\
	     &Gen[gencosti].cost_startup,&Gen[gencosti].cost_shutdown,&Gen[gencosti].cost_ncoeffs,&Gen[gencosti].cost_alpha, \
	     &Gen[gencosti].cost_beta,&Gen[gencosti].cost_gamma);
      gencosti++;
      
    }
    
    if(i >= br_start_line && i < br_end_line) {
      if(strcmp(line,"\n")==0 || strcmp(line,"\r\n")==0) continue;
      sscanf(line,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %d",&Branch[bri].fbus,&Branch[bri].tbus, \
	     &Branch[bri].r,&Branch[bri].x,&Branch[bri].b,		\
	     &Branch[bri].rateA,&Branch[bri].rateB,&Branch[bri].rateC, \
	     &Branch[bri].tapratio,&Branch[bri].phaseshift,&Branch[bri].status);
      if(!Branch[bri].tapratio) Branch[bri].tapratio = 1.0;
      Branch[bri].phaseshift *= PETSC_PI/180.0;

      intbusnum = busext2intmap[Branch[bri].fbus];
      Branch[bri].internal_i = intbusnum;

      intbusnum = busext2intmap[Branch[bri].tbus];
      Branch[bri].internal_j = intbusnum;

      PetscInt lineididx = 0;
      for(linenum = 0; linenum < bri-1; linenum++) {
	if(Branch[bri].internal_i == Branch[linenum].internal_i && Branch[bri].internal_j == Branch[linenum].internal_j) lineididx += 1;
      }
	
      /* MatPower does not have ids for lines. Using bri+1 as the id */
      snprintf(Branch[bri].ckt,3,"%-2d",1+lineididx);

      /* Compute self and transfer admittances */
      PetscScalar R,X,Bc,B,G,Zm,tap,shift,tap2,tapr,tapi;
      R = Branch[bri].r;
      X = Branch[bri].x;
      Bc = Branch[bri].b;

      Zm = R*R + X*X;
      G  = R/Zm;
      B  = -X/Zm;

      tap = Branch[bri].tapratio;
      shift = Branch[bri].phaseshift;
      tap2 = tap*tap;
      tapr = tap*cos(shift);
      tapi = tap*sin(shift);

      Branch[bri].yff[0] = G/tap2; 
      Branch[bri].yff[1] = (B+Bc/2.0)/tap2;
      
      Branch[bri].yft[0] = -(G*tapr - B*tapi)/tap2;
      Branch[bri].yft[1] = -(B*tapr + G*tapi)/tap2;

      Branch[bri].ytf[0] = -(G*tapr + B*tapi)/tap2;
      Branch[bri].ytf[1] = -(B*tapr - G*tapi)/tap2;

      Branch[bri].ytt[0] = G;
      Branch[bri].ytt[1] = B+Bc/2.0;

      bri++;
    }
  }
  fclose(fp);
  
  PetscFunctionReturn(0);
}

/*
  PSReadDyrData - Reads the data file with dynamic models 

  Input Parameter
+  PS      - The PS object
-  dyrfile - The name of the dyr file

*/

PetscErrorCode PSReadDyrData(PS ps,const char dyrfile[])
{
  PetscErrorCode ierr;
  FILE           *fp;
  char           line[MAXLINE];
  PetscInt       extbusnum;
  PetscInt       *busext2intmap=ps->busext2intmap;
  PSBUS          Bus;
  PSGEN          Gen;
  PSLOAD         Load;
  PetscInt       i,k;
  char           gentype[16],exctype[16],turbgovtype[16],stabtype[16],loadtype[16];
  struct _p_DYNGenModel dyngen;
  struct _p_DYNExcModel dynexc;
  struct _p_DYNTurbgovModel dynturbgov;
  struct _p_DYNStabModel dynstab;
  struct _p_DYNLoadModel dynload;
  char           *out;

  PetscFunctionBegin;

  if(ps->comm->type != PETSC_COMM_SELF && ps->comm->rank != 0) { PetscFunctionReturn(0);}

  fp = fopen(dyrfile,"r");
  /* Check for valid file */
  if (fp == NULL) {
    SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Cannot open file %s",dyrfile);
  }

  /* Initialize hasxxx flag to PETSC_FALSE for all generators */
  for(i=0; i < ps->nbus; i++) {
    Bus = &ps->bus[i];
    for(k=0; k < Bus->ngen; k++) {
      Gen = &ps->gen[Bus->gidx[k]];
      Gen->hasexc = PETSC_FALSE;
      Gen->hasturbgov = PETSC_FALSE;
      Gen->hasstab = PETSC_FALSE;
    }
  }
      
  char *excid,*stabid;
  char genid[2],turbgovid[2];
  char loadid[2];
  char temp[20];
  PetscBool multiline;
  while((out = fgets(line,MAXLINE,fp)) != NULL) {
    PetscBool linehasdyngen=PETSC_FALSE,linehasdynexc=PETSC_FALSE,linehasdynturbgov=PETSC_FALSE,linehasdynstab=PETSC_FALSE,linehasdynload=PETSC_FALSE;
    
    if(strcmp(line,"\r\n") == 0 || strcmp(line,"\n") == 0) {
      continue; /* Skip blank lines */
    }

    multiline = PETSC_FALSE;
    /*Handling multiline format (seperated by ' ' and end flag is '/')*/
    char   totalline[MAXLINE];
    PetscInt totallinelen = 0;
    char*  pos;
    
    if(strstr(line,"/") == NULL ){
      multiline = PETSC_TRUE;
      memset(totalline , 0, MAXLINE);
      
      /*Read multiple line and remove newline*/
      if((pos = strstr(line, "\r\n")) != NULL){
        *pos = ' ';
        *(pos+1) = ' ';
      } 
      else if((pos = strstr(line, "\n")) != NULL) *pos = ' ';

      ierr = PetscStrcpy (totalline + totallinelen, line); CHKERRQ(ierr);
      totallinelen += strlen(line);

      while(strstr(line,"/") == NULL){        
        out = fgets(line,MAXLINE,fp);
        if((pos = strstr(line, "\r\n")) != NULL){
          *pos = ' ';
          *(pos+1) = ' ';
        } 
        else if((pos = strstr(line, "\n")) != NULL) *pos = ' ';
        
        ierr = PetscStrcpy (totalline + totallinelen, line); CHKERRQ(ierr);
        totallinelen += strlen(line);
      }
    } else {
      if(strstr(line,",") == NULL) {
	ierr = PetscStrcpy(totalline,line);CHKERRQ(ierr);
	totallinelen += strlen(line);
      }
    }

    if(multiline || (strstr(line,",") == NULL)) {
	/*Change ' ' -> ',' */
	PetscInt linepos = 0;
	i = 0;
	
	if((pos = strtok(totalline, " ")) != NULL){
	  ierr = PetscStrcpy(line + linepos, pos); CHKERRQ(ierr); 
	  linepos += strlen(pos);
	  line[linepos] = ',';
	  linepos++;
	  i++;
	}
	while((pos = strtok(NULL, " ")) != NULL){
	  ierr = PetscStrcpy(line + linepos, pos); CHKERRQ(ierr); 
	  linepos += strlen(pos);
	  
	  /*Make Third column has 2 byte*/
	  if(i == 2){
	    if(strlen(pos) == 1){
	      line[linepos] = ' ';
	      linepos++;
	    }
	  }
	  
	  /*Teminate loop*/
	  if(strstr(pos,"/") != NULL ){
	    line[linepos-2] = ' '; // remove last ',' before '/'
	    break;
	  }
	  
	  line[linepos] = ',';
	  linepos++;
	  i++;
	}
	
	line[linepos] = 0x00;
#if defined DEBUGPS
	ierr = PetscPrintf(PETSC_COMM_SELF,"Line = %s\n",line);CHKERRQ(ierr);      
#endif      
    }

    ierr = DYNGenModelReadModelType(line,gentype);CHKERRQ(ierr);
    ierr = PetscStrcmp(gentype,DYNGENNONE,&linehasdyngen);CHKERRQ(ierr);

    if(!linehasdyngen) {
      ierr = DYNGenModelSetType(&dyngen,gentype);CHKERRQ(ierr);

      /* ierr = DYNGenModelGetBusnumID(&dyngen,&extbusnum,&genid);CHKERRQ(ierr); */
      sscanf(line,"%d,%[^,],%[^,]",&extbusnum,temp,genid);

      Bus = &ps->bus[busext2intmap[extbusnum]];
      /* Find the gen with matching gen id */
      for(i=0; i < Bus->ngen;i++) {
        ierr = PetscStrcmp(ps->gen[Bus->gidx[i]].id,genid,&linehasdyngen);CHKERRQ(ierr);
        if(linehasdyngen) {
          Gen = &ps->gen[Bus->gidx[i]];
          if(Gen->status) {
            ierr = DYNGenModelReadData(&dyngen,line,Gen->mbase,ps->MVAbase);CHKERRQ(ierr);
	    /* Set the numnber of variables for this dynamic generator model */
	    ierr = DYNGenModelGetNvar(&dyngen,&dyngen.nvar);CHKERRQ(ierr);
            ierr = DYNGenModelCopy(&Gen->dyngen,&dyngen);CHKERRQ(ierr);
	    Gen->dyngensetup = PETSC_TRUE;
            break;
          } else {
            ierr = DYNGenModelDestroy(&dyngen);CHKERRQ(ierr);
          }
        }
      }
      continue;
    }

    ierr = DYNExcModelReadModelType(line,exctype);CHKERRQ(ierr);
    ierr = PetscStrcmp(exctype,DYNEXCNONE,&linehasdynexc);CHKERRQ(ierr);

    if(!linehasdynexc) {
      ierr = DYNExcModelSetType(&dynexc,exctype);CHKERRQ(ierr);
      ierr = DYNExcModelReadData(&dynexc,line);CHKERRQ(ierr);
      ierr = DYNExcModelGetBusnumID(&dynexc,&extbusnum,&excid);CHKERRQ(ierr);

      Bus = &ps->bus[busext2intmap[extbusnum]];
      /* Find the gen with matching exc id */
      for(i=0; i < Bus->ngen;i++) {
	ierr = PetscStrcmp(ps->gen[Bus->gidx[i]].id,excid,&linehasdynexc);CHKERRQ(ierr);
	if(linehasdynexc) {
	  Gen = &ps->gen[Bus->gidx[i]];
	  if(Gen->status) {
	    Gen->hasexc = PETSC_TRUE;
	    ierr = DYNExcModelCopy(&Gen->dynexc,&dynexc);CHKERRQ(ierr);
	    break;
	  } else {
	    ierr = DYNExcModelDestroy(&dynexc);CHKERRQ(ierr);
	  }
	}
      }
      continue;
    }

    ierr = DYNTurbgovModelReadModelType(line,turbgovtype);CHKERRQ(ierr);
    ierr = PetscStrcmp(turbgovtype,DYNTURBGOVNONE,&linehasdynturbgov);CHKERRQ(ierr);
    if(!linehasdynturbgov) {
      sscanf(line,"%d,%[^,],%[^,]",&extbusnum,temp,turbgovid);
      Bus = &ps->bus[busext2intmap[extbusnum]];
      ierr = DYNTurbgovModelSetType(&dynturbgov,turbgovtype);CHKERRQ(ierr);
      /* Find the gen with matching gen id */
      for(i=0; i < Bus->ngen;i++) {
	//	ierr = DYNTurbgovModelGetBusnumID(&dynturbgov,&extbusnum,&turbgovid);CHKERRQ(ierr);
	ierr = PetscStrcmp(ps->gen[Bus->gidx[i]].id,turbgovid,&linehasdynturbgov);CHKERRQ(ierr);
	if(linehasdynturbgov) {
	  Gen = &ps->gen[Bus->gidx[i]];
	  if(Gen->status) {
	    Gen->hasturbgov = PETSC_TRUE;
	    ierr = DYNTurbgovModelReadData(&dynturbgov,line,Gen->mbase,ps->MVAbase);CHKERRQ(ierr);
	    ierr = DYNTurbgovModelCopy(&Gen->dynturbgov,&dynturbgov);CHKERRQ(ierr);
	    break;
	  } else {
	    ierr = DYNTurbgovModelDestroy(&dynturbgov);CHKERRQ(ierr);
	  }
	}
      }
      continue;
    }

    ierr = DYNStabModelReadModelType(line,stabtype);CHKERRQ(ierr);
    ierr = PetscStrcmp(stabtype,DYNSTABNONE,&linehasdynstab);CHKERRQ(ierr);

    if(!linehasdynstab) {
      ierr = DYNStabModelSetType(&dynstab,stabtype);CHKERRQ(ierr);
      ierr = DYNStabModelReadData(&dynstab,line);CHKERRQ(ierr);
      ierr = DYNStabModelGetBusnumID(&dynstab,&extbusnum,&stabid);CHKERRQ(ierr);

      Bus = &ps->bus[busext2intmap[extbusnum]];
      /* Find the gen with matching stab id */
      for(i=0; i < Bus->ngen;i++) {
	ierr = PetscStrcmp(ps->gen[Bus->gidx[i]].id,stabid,&linehasdynstab);CHKERRQ(ierr);
	if(linehasdynstab) {
	  Gen = &ps->gen[Bus->gidx[i]];
	  if(Gen->status) {
	    Gen->hasstab = PETSC_TRUE;
	    ierr = DYNStabModelCopy(&Gen->dynstab,&dynstab);CHKERRQ(ierr);
	    break;
	  } else {
	    ierr = DYNStabModelDestroy(&dynstab);CHKERRQ(ierr);
	  }
	}
      }
      continue;
    }

    ierr = DYNLoadModelReadModelType(line,loadtype);CHKERRQ(ierr);
    ierr = PetscStrcmp(loadtype,DYNLOADNONE,&linehasdynload);CHKERRQ(ierr);

    if(!linehasdynload) {
      ierr = DYNLoadModelSetType(&dynload,loadtype);CHKERRQ(ierr);

      /* ierr = DYNLoadModelGetBusnumID(&dynload,&extbusnum,&loadid);CHKERRQ(ierr); */
      sscanf(line,"%d,%[^,],%[^,]",&extbusnum,temp,loadid);

      Bus = &ps->bus[busext2intmap[extbusnum]];
      /* Find the load with matching load id */
      for(i=0; i < Bus->nload;i++) {
        ierr = PetscStrcmp(ps->load[Bus->lidx[i]].id,loadid,&linehasdynload);CHKERRQ(ierr);
        if(linehasdynload) {
          Load = &ps->load[Bus->lidx[i]];
          if(Load->status) {
            ierr = DYNLoadModelReadData(&dynload,line,Load->pl*ps->MVAbase,ps->MVAbase);CHKERRQ(ierr); /* Load does not have a MBase in the power flow data, we use the real power load instead. The machine (induction motor) base will be given in the dyr file */
	    /* Set the numnber of variables for this dynamic load model */
	    ierr = DYNLoadModelGetNvar(&dynload,&dynload.nvar);CHKERRQ(ierr);
            ierr = DYNLoadModelCopy(&Load->dynload,&dynload);CHKERRQ(ierr);
	    Load->dynloadsetup = PETSC_TRUE;
            break;
          } else {
            ierr = DYNLoadModelDestroy(&dynload);CHKERRQ(ierr);
          }
        }
      }
      continue;
    }
  }
  fclose(fp);

  PetscFunctionReturn(0);
}

/*
  PSGetNumGlobalLines - Gets the total number of lines in the PS network

  Input Parameters:
+ PS     - the PS network object
- Nlines - the total number of lines (including branches, transformers) in the network
*/
PetscErrorCode PSGetNumGlobalLines(PS ps,PetscInt *Nlines)
{
  PetscFunctionBegin;
  *Nlines = ps->Nbranch;
  PetscFunctionReturn(0);
}

/*
  PSGetNumGlobalBuses - Gets the total number of buses in the PS network

  Input Parameters:
+ PS     - the PS network object
- Nbuses - the total number of buses in the network
*/
PetscErrorCode PSGetNumGlobalBuses(PS ps,PetscInt *Nbuses)
{
  PetscFunctionBegin;
  *Nbuses = ps->Nbus;
  PetscFunctionReturn(0);
}

/*
  PSGetLineConnectivity - Gets the connectivitiy (from bus - to bus connection) of all lines

  Input Parameters:
+ PS - the PS object
. Nlines - total number of lines
- lineconn  - the line connectivity area (of size 2*Nlines)

  Notes:
  The line connectivity for line i is stored as lines[2*i] = frombus_i, lines[2*i+1] = tobus_i
*/
PetscErrorCode PSGetLineConnectivity(PS ps,PetscInt Nlines,int lineconn[])
{
  PetscInt       i,fbus,tbus;
  PSLINE         line=ps->line;

  PetscFunctionBegin;
  for(i=0; i < Nlines; i++) {
    fbus = line[i].internal_i;
    tbus = line[i].internal_j;
    lineconn[2*i]   = fbus;
    lineconn[2*i+1] = tbus;
  }
  PetscFunctionReturn(0);
}

/*
  PSSetUp - Sets up the PS object to be ready to be used by the application

  Input Parameters:
. PS - the PS object

  Notes:
  PSSetUp 
    i) creates the underlying DMNetwork object
   ii) distributes the DMNetwork if used in parallel.
  iii) creates PSBUS,PSBRANCH,PSGEN,PSLOAD objects
*/
PetscErrorCode PSSetUp(PS ps)
{
  PetscErrorCode ierr;
  DM             networkdm;
  MPI_Comm       mpicomm=ps->comm->type;
  PetscInt       Nlines,Nbuses,Ngbuses=PETSC_DETERMINE,Nglines=PETSC_DETERMINE;
  PetscInt       numbusvariables; /* Number of variables at each bus..set by the application */
  PetscInt       i;
  PetscBool      match;
  DYNGenModel    dyngen;
  DYNExcModel    dynexc;
  DYNTurbgovModel dynturbgov;
  DYNStabModel    dynstab;
  DYNLoadModel    dynload;
  void            *component;
  PetscInt        key;
  PetscInt        numComponents;


  PetscFunctionBegin;
  if(ps->setupcalled) PetscFunctionReturn(0);

  /* Create empty DMNetwork object */
  ierr = DMNetworkCreate(mpicomm,&networkdm);CHKERRQ(ierr);
  
  /* Register bus,branch,gen,load objects */
  ierr = DMNetworkRegisterComponent(networkdm,"PSLINE",sizeof(struct _p_PSLINE),&ps->compkey[0]);CHKERRQ(ierr);
  ierr = DMNetworkRegisterComponent(networkdm,"PSBUS",sizeof(struct _p_PSBUS),&ps->compkey[1]);CHKERRQ(ierr);
  ierr = DMNetworkRegisterComponent(networkdm,"PSGEN",sizeof(struct _p_PSGEN),&ps->compkey[2]);CHKERRQ(ierr);
  ierr = DMNetworkRegisterComponent(networkdm,"PSLOAD",sizeof(struct _p_PSLOAD),&ps->compkey[3]);CHKERRQ(ierr);
  
  if(ps->app == APP_DYNSIM) {
    /* Register dynamic generator models with DMNetwork */
    for(i=0;i < ngenmodelsregistered;i++) { 
      ierr = DMNetworkRegisterComponent(networkdm,DYNGenModelList[i].name,DYNGenModelList[i].sizeofstruct,&DYNGenModelList[i].key);CHKERRQ(ierr);
    }
    /* Register exciter models with DMNetwork */
    for(i=0;i < nexcmodelsregistered;i++) { 
      ierr = DMNetworkRegisterComponent(networkdm,DYNExcModelList[i].name,DYNExcModelList[i].sizeofstruct,&DYNExcModelList[i].key);CHKERRQ(ierr);
    }

    /* Register turbine-governor models with DMNetwork */
    for(i=0;i < nturbgovmodelsregistered;i++) { 
      ierr = DMNetworkRegisterComponent(networkdm,DYNTurbgovModelList[i].name,DYNTurbgovModelList[i].sizeofstruct,&DYNTurbgovModelList[i].key);CHKERRQ(ierr);
    }

    /* Register stabilizer models with DMNetwork */
    for(i=0;i < nstabmodelsregistered;i++) { 
      ierr = DMNetworkRegisterComponent(networkdm,DYNStabModelList[i].name,DYNStabModelList[i].sizeofstruct,&DYNStabModelList[i].key);CHKERRQ(ierr);
    }

    /* Register load models with DMNetwork */
    for(i=0;i < nloadmodelsregistered;i++) { 
      ierr = DMNetworkRegisterComponent(networkdm,DYNLoadModelList[i].name,DYNLoadModelList[i].sizeofstruct,&DYNLoadModelList[i].key);CHKERRQ(ierr);
    }

  }
  /* Get the total number of buses and lines */
  /* Note that when the read is read from XXXReadMatPowerData, only P0 reads the data and has
     NumLines and NumBuses set, all the other processors don't have any bus, branch, gen, load data
     set.
  */
  ierr = PSGetNumGlobalLines(ps,&Nlines);CHKERRQ(ierr);
  ierr = PSGetNumGlobalBuses(ps,&Nbuses);CHKERRQ(ierr);

  /* Set up edge connectivity */
  int *lineconn;
  ierr = PetscCalloc1(2*Nlines,&lineconn);CHKERRQ(ierr);
  ierr = PSGetLineConnectivity(ps,Nlines,lineconn);CHKERRQ(ierr);

  /* Set sizes for the network */
  ierr = DMNetworkSetSizes(networkdm,1,0,&Nbuses,&Nlines,&Ngbuses,&Nglines);CHKERRQ(ierr);
  /* Set edge connectivity */
  ierr = DMNetworkSetEdgeList(networkdm,&lineconn,NULL);CHKERRQ(ierr);
  /* Set up network layout */
  ierr = DMNetworkLayoutSetUp(networkdm);CHKERRQ(ierr);
  ierr = PetscFree(lineconn);CHKERRQ(ierr);

  /* Add network components and variables */
  PetscInt eStart,eEnd,vStart,vEnd,j;
  PSGEN    gen;
  PSLOAD   load;
  /* Associate and copy line data to the edges in networkdm, add number of variables if applicable */
  ierr = DMNetworkGetEdgeRange(networkdm,&eStart,&eEnd);CHKERRQ(ierr);
  for(i=eStart; i < eEnd; i++) {
    ierr = DMNetworkAddComponent(networkdm,i,ps->compkey[0],&ps->line[i-eStart]);CHKERRQ(ierr);
  }
  /* Associate and copy bus, gen, branch to the vertices of the networkdm */
  ierr = DMNetworkGetVertexRange(networkdm,&vStart,&vEnd);CHKERRQ(ierr);
  for(i=vStart; i < vEnd; i++) {
    /* Set the number of variables for buses */
    if(ps->app == APP_DYNSIM || ps->app == APP_ACPF) numbusvariables = 2;
    else if(ps->app == APP_ACOPF) numbusvariables = 2 + 2*ps->bus[i-vStart].ngen;
    else SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Application not supported");

    ierr = DMNetworkAddNumVariables(networkdm,i,numbusvariables);CHKERRQ(ierr); /* Bus variables */
    ierr = DMNetworkAddComponent(networkdm,i,ps->compkey[1],&ps->bus[i-vStart]);CHKERRQ(ierr);

    for(j=0; j < ps->bus[i-vStart].ngen; j++) {
      /* Add generator */
      gen = &ps->gen[ps->bus[i-vStart].gidx[j]];
      ierr = DMNetworkAddComponent(networkdm,i,ps->compkey[2],gen);CHKERRQ(ierr);
      if(ps->app == APP_DYNSIM) {
        if(!gen->status) continue;
        /* Add dynamic generator model */
        PetscInt Nvar;
        ierr = DYNGenModelGetNvar(&gen->dyngen,&Nvar);CHKERRQ(ierr);
        ierr = DMNetworkAddNumVariables(networkdm,i,Nvar);CHKERRQ(ierr);
        /* Add the implementation type as a component to the networkdm */
        /* Get the component key by comparing the name (This is ugly and slow!!) */
        PetscInt k;
        for(k=0; k < ngenmodelsregistered;k++) {
          ierr = PetscStrcmp(DYNGenModelList[k].name,gen->dyngen.type,&match);CHKERRQ(ierr);
          if(match) {
            ierr = DMNetworkAddComponent(networkdm,i,DYNGenModelList[k].key,gen->dyngen.data);CHKERRQ(ierr);
            break;
          }
	}

	/* Add exciter model  */
	if(gen->hasexc) {
	  ierr = DYNExcModelGetNvar(&gen->dynexc,&Nvar);CHKERRQ(ierr);
	  ierr = DMNetworkAddNumVariables(networkdm,i,Nvar);CHKERRQ(ierr);
	  /* Add the implementation type as a component to the networkdm */
	  /* Get the component key by comparing the name (This is ugly and slow!!) */
	  for(k=0; k < nexcmodelsregistered;k++) {
	    ierr = PetscStrcmp(DYNExcModelList[k].name,gen->dynexc.type,&match);CHKERRQ(ierr);
	    if(match) {
	      ierr = DMNetworkAddComponent(networkdm,i,DYNExcModelList[k].key,gen->dynexc.data);CHKERRQ(ierr);
	      break;
	    }
	  }
	}
	
	/* Add turbine-governor model  */
	if(gen->hasturbgov) {
	  ierr = DYNTurbgovModelGetNvar(&gen->dynturbgov,&Nvar);CHKERRQ(ierr);
	  ierr = DMNetworkAddNumVariables(networkdm,i,Nvar);CHKERRQ(ierr);

	  /* Add the implementation type as a component to the networkdm */
	  /* Get the component key by comparing the name (This is ugly and slow!!) */
	  for(k=0; k < nturbgovmodelsregistered;k++) {
	    ierr = PetscStrcmp(DYNTurbgovModelList[k].name,gen->dynturbgov.type,&match);CHKERRQ(ierr);
	    if(match) {
	      ierr = DMNetworkAddComponent(networkdm,i,DYNTurbgovModelList[k].key,gen->dynturbgov.data);CHKERRQ(ierr);
	      break;
	    }
	  }
	}

	/* Add stabilizer model  */
	if(gen->hasstab) {
	  ierr = DYNStabModelGetNvar(&gen->dynstab,&Nvar);CHKERRQ(ierr);
	  ierr = DMNetworkAddNumVariables(networkdm,i,Nvar);CHKERRQ(ierr);

	  /* Add the implementation type as a component to the networkdm */
	  /* Get the component key by comparing the name (This is ugly and slow!!) */
	  for(k=0; k < nstabmodelsregistered;k++) {
	    ierr = PetscStrcmp(DYNStabModelList[k].name,gen->dynstab.type,&match);CHKERRQ(ierr);
	    if(match) {
	      ierr = DMNetworkAddComponent(networkdm,i,DYNStabModelList[k].key,gen->dynstab.data);CHKERRQ(ierr);
	      break;
	    }
	  }
	}
      }
    }

    /* Loads */
    for(j=0; j < ps->bus[i-vStart].nload; j++) { 
      /* Add load */
      load = &ps->load[ps->bus[i-vStart].lidx[j]];
      ierr = DMNetworkAddComponent(networkdm,i,ps->compkey[3],load);CHKERRQ(ierr);

      if(ps->app == APP_DYNSIM) {
	
	if(!load->status) continue;

	/* Add dynamic load model */
       	PetscInt Nvar;
	ierr = DYNLoadModelGetNvar(&load->dynload,&Nvar);CHKERRQ(ierr);
	ierr = DMNetworkAddNumVariables(networkdm,i,Nvar);CHKERRQ(ierr);
	
	/* Add the implementation type as a component to the networkdm */
	/* Get the component key by comparing the name (This is ugly and slow!!) */
	PetscInt k;
	for(k=0; k < nloadmodelsregistered;k++) {
	  ierr = PetscStrcmp(DYNLoadModelList[k].name,load->dynload.type,&match);CHKERRQ(ierr);
	  if(match) {
	    ierr = DMNetworkAddComponent(networkdm,i,DYNLoadModelList[k].key,load->dynload.data);CHKERRQ(ierr);
	    break;
	  }
	}
      }
    }
  }

  /* Set up DM for use */
  ierr = DMSetUp(networkdm);CHKERRQ(ierr);

  ierr = DMNetworkDistribute(&networkdm,0);CHKERRQ(ierr);
  ps->networkdm = networkdm;
  //  if(ps->comm->size > 1) {
  //    DM distnetworkdm;
    /* Network partitioning and distribution of data */
  //    ierr = DMNetworkDistribute(networkdm,0,&distnetworkdm);CHKERRQ(ierr);
  //    ierr = DMDestroy(&networkdm);
  //    ps->networkdm = distnetworkdm;
  //  } else ps->networkdm = networkdm;

  ierr = DMNetworkGetEdgeRange(ps->networkdm,&eStart,&eEnd);CHKERRQ(ierr);
  ierr = DMNetworkGetVertexRange(ps->networkdm,&vStart,&vEnd);CHKERRQ(ierr);

  /* Set local sizes of buses and branches */
  ps->nbus = vEnd - vStart;
  ps->nbranch = eEnd - eStart;

  if(ps->comm->size > 1) {
    if(ps->comm->rank == 0) {
      /* Free PSBUS,PSBRANCH,PSGEN, and PSLOAD on root */
      ierr = PetscFree(ps->bus);CHKERRQ(ierr);
      ierr = PetscFree(ps->line);CHKERRQ(ierr);
      ierr = PSGENDestroy(ps);CHKERRQ(ierr);
      ierr = PSLOADDestroy(ps);CHKERRQ(ierr);
      ierr = PetscFree(ps->load);CHKERRQ(ierr);      
      ierr = PetscFree(ps->busext2intmap);CHKERRQ(ierr);
    }

    /* Broadcast global Nbus,Ngen,Nbranch, Nload,and maxbusnum to all processors */
    PetscInt temp[5];
    /* Pack variables */
    temp[0] = ps->Nbus;
    temp[1] = ps->Ngen;
    temp[2] = ps->Nbranch;
    temp[3] = ps->Nload;
    temp[4] = ps->maxbusnum;
    ierr = MPI_Bcast(temp,5,MPI_INT,0,ps->comm->type);CHKERRQ(ierr);
    /* Unpack */
    ps->Nbus = temp[0];
    ps->Ngen = temp[1];
    ps->Nbranch = temp[2];
    ps->Nload   = temp[3];
    ps->maxbusnum = temp[4];

    /* Recreate busext2intmap..this will map the local bus numbers to external numbers */
    ierr = PetscCalloc1(ps->maxbusnum+1,&ps->busext2intmap);CHKERRQ(ierr);
    for(i=0; i < ps->maxbusnum+1; i++) ps->busext2intmap[i] = -1;

    ps->ngen = ps->nload = 0;
  /* Get the local number of gens and loads */
    for(i= vStart; i < vEnd; i++) {
      ierr = DMNetworkGetNumComponents(ps->networkdm,i,&numComponents);CHKERRQ(ierr);
      for(j=0; j < numComponents; j++) {
	ierr = DMNetworkGetComponent(ps->networkdm,i,j,&key,&component);CHKERRQ(ierr);
	if(key == ps->compkey[2]) ps->ngen++;
	else if(key == ps->compkey[3]) ps->nload++;
      }
    }

    /* Create local PSBUS, PSBRANCH, PSGEN, and PSLOAD */
    ierr = PetscCalloc1(ps->nbus,&ps->bus);CHKERRQ(ierr);
    ierr = PetscCalloc1(ps->nbranch,&ps->line);CHKERRQ(ierr);
    ierr = PetscCalloc1(ps->ngen,&ps->gen);CHKERRQ(ierr);
    ierr = PetscCalloc1(ps->nload,&ps->load);CHKERRQ(ierr);

    /* Copy line data from DMNetwork data array to PSLINE */
    for(i=eStart; i < eEnd; i++) {
      ierr = DMNetworkGetComponent(ps->networkdm,i,0,&key,&component);CHKERRQ(ierr);
      ierr = PetscMemcpy(&ps->line[i-eStart],component,sizeof(struct _p_PSLINE));CHKERRQ(ierr);
    }
    PetscInt genj=0,loadj=0;
    PetscInt genctr,loadctr;
    /* Copy bus, gen, branch data from DMNetwork data array to PSBUS, PSGEN, and PSBRANCH */
    for(i=vStart; i < vEnd; i++) {
      genctr = loadctr = 0;
      ierr = DMNetworkGetNumComponents(ps->networkdm,i,&numComponents);CHKERRQ(ierr);
      for(j=0; j < numComponents; j++) {
	ierr = DMNetworkGetComponent(ps->networkdm,i,j,&key,&component);CHKERRQ(ierr);
	if(key == ps->compkey[1]) {
	  ierr = PetscMemcpy(&ps->bus[i-vStart],component,sizeof(struct _p_PSBUS));CHKERRQ(ierr);
	  /* Set external to internal mapping */
	  ps->busext2intmap[ps->bus[i-vStart].bus_i] = ps->bus[i-vStart].internal_i = i-vStart;
	} else if(key == ps->compkey[2]) {
	  ierr = PetscMemcpy(&ps->gen[genj],component,sizeof(struct _p_PSGEN));CHKERRQ(ierr);
	  if(ps->app == APP_DYNSIM) {
	    if(!ps->gen[genj].status) {
	      ps->bus[i-vStart].gidx[genctr++] = genj++;
	      continue;
	    }

	    PetscInt k;
	    /* Generator model */
	    for(k=0; k < ngenmodelsregistered;k++) {
	      ierr = PetscStrcmp(ps->gen[genj].dyngen.type,DYNGenModelList[k].name,&match);CHKERRQ(ierr);
	      if(match) {
		ierr = DYNGenModelSetType(&ps->gen[genj].dyngen,DYNGenModelList[k].name);CHKERRQ(ierr);
		j++;
		ierr = DMNetworkGetComponent(ps->networkdm,i,j,&key,&component);CHKERRQ(ierr);
		ierr = PetscMemcpy(ps->gen[genj].dyngen.data,component,DYNGenModelList[k].sizeofstruct);CHKERRQ(ierr);
		ierr = DYNGenModelGetNvar(&ps->gen[genj].dyngen,&ps->gen[genj].dyngen.nvar);CHKERRQ(ierr);
		break;
	      }
	    }

	    /* Exciter model */
	    if(ps->gen[genj].hasexc) {
	      /* Get exciter model */
	      for(k=0; k < nexcmodelsregistered;k++) {
		ierr = PetscStrcmp(ps->gen[genj].dynexc.type,DYNExcModelList[k].name,&match);CHKERRQ(ierr);
		if(match) {
		  ierr = DYNExcModelSetType(&ps->gen[genj].dynexc,DYNExcModelList[k].name);CHKERRQ(ierr);
		  j++;
		  ierr = DMNetworkGetComponent(ps->networkdm,i,j,&key,&component);CHKERRQ(ierr);
		  ierr = PetscMemcpy(ps->gen[genj].dynexc.data,component,DYNExcModelList[k].sizeofstruct);CHKERRQ(ierr);
		  break;
		}
	      }
	    } 

	    /* Turbine-governor model */
	    if(ps->gen[genj].hasturbgov) {
	      /* Get turbine-governor model */
	      for(k=0; k < nturbgovmodelsregistered;k++) {
		ierr = PetscStrcmp(ps->gen[genj].dynturbgov.type,DYNTurbgovModelList[k].name,&match);CHKERRQ(ierr);
		if(match) {
		  ierr = DYNTurbgovModelSetType(&ps->gen[genj].dynturbgov,DYNTurbgovModelList[k].name);CHKERRQ(ierr);
		  j++;
		  ierr = DMNetworkGetComponent(ps->networkdm,i,j,&key,&component);CHKERRQ(ierr);
		  ierr = PetscMemcpy(ps->gen[genj].dynturbgov.data,component,DYNTurbgovModelList[k].sizeofstruct);CHKERRQ(ierr);
		  break;
		}
	      }
	    } 

	    /* Stabilizer model */
	    if(ps->gen[genj].hasstab) {
	      /* Get stabilizer model */
	      for(k=0; k < nstabmodelsregistered;k++) {
		ierr = PetscStrcmp(ps->gen[genj].dynstab.type,DYNStabModelList[k].name,&match);CHKERRQ(ierr);
		if(match) {
		  ierr = DYNStabModelSetType(&ps->gen[genj].dynstab,DYNStabModelList[k].name);CHKERRQ(ierr);
		  j++;
		  ierr = DMNetworkGetComponent(ps->networkdm,i,j,&key,&component);CHKERRQ(ierr);
		  ierr = PetscMemcpy(ps->gen[genj].dynstab.data,component,DYNStabModelList[k].sizeofstruct);CHKERRQ(ierr);
		  break;
		}
	      }
	    } 	    
	  }
	  ps->bus[i-vStart].gidx[genctr++] = genj++;
	} else if(key == ps->compkey[3]) {
	  ierr = PetscMemcpy(&ps->load[loadj],component,sizeof(struct _p_PSLOAD));CHKERRQ(ierr);
	  if(ps->app == APP_DYNSIM) {
	    if(!ps->load[loadj].status) {
	      ps->bus[i-vStart].lidx[loadctr++] = loadj++;
	      continue;
	    }
	    PetscInt k;
	    /* Dynamic load model */
	    for(k=0; k < nloadmodelsregistered;k++) {
	      ierr = PetscStrcmp(ps->load[loadj].dynload.type,DYNLoadModelList[k].name,&match);CHKERRQ(ierr);
	      if(match) {
		ierr = DYNLoadModelSetType(&ps->load[loadj].dynload,DYNLoadModelList[k].name);CHKERRQ(ierr);
		j++;
		ierr = DMNetworkGetComponent(ps->networkdm,i,j,&key,&component);CHKERRQ(ierr);
		ierr = PetscMemcpy(ps->load[loadj].dynload.data,component,DYNLoadModelList[k].sizeofstruct);CHKERRQ(ierr);
		/* Set the numnber of variables for this dynamic load model */
		ierr = DYNLoadModelGetNvar(&ps->load[loadj].dynload,&ps->load[loadj].dynload.nvar);CHKERRQ(ierr);


		break;
	      }
	    }
	  }
	  ps->bus[i-vStart].lidx[loadctr++] = loadj++;
	}
      }
    }
  } else {
    ps->ngen = ps->Ngen;
    ps->nload = ps->Nload;
  }

  /* Set up
     (a) connectivity information for lines and buses 
     (b) bus ghosted status
     (c) incident generators at bus
     (d) incident loads at bus
     (e) sets the starting location for the variables for this bus in the given application
         state vector
     (f) set up the number of differential and algebraic equations for APP_DYNSIM
     (g) sets up starting location for generator and exciter variables (for APP_DYNSIM) 
     (h) 
  */
  PetscInt nlines,k;
  const PetscInt *connnodes,*connlines;
  for(i=eStart; i < eEnd; i++) {
    ierr = DMNetworkGetConnectedVertices(ps->networkdm,i,&connnodes);CHKERRQ(ierr);
    ps->line[i-eStart].connbuses[0] = &ps->bus[connnodes[0]-vStart];
    ps->line[i-eStart].connbuses[1] = &ps->bus[connnodes[1]-vStart];
    ps->line[i-eStart].internal_i   = ps->bus[connnodes[0]-vStart].internal_i;
    ps->line[i-eStart].internal_j   = ps->bus[connnodes[1]-vStart].internal_i;

    /* Starting location in the local vector */
    ierr = DMNetworkGetVariableOffset(ps->networkdm,i,&ps->line[i].startloc);CHKERRQ(ierr);
    /* Starting location in the global vector */
    ierr = DMNetworkGetVariableGlobalOffset(ps->networkdm,i,&ps->line[i].startlocglob);CHKERRQ(ierr);
  }

  for(i=vStart; i < vEnd; i++) {
    /* Connected lines */
    ierr = DMNetworkGetSupportingEdges(ps->networkdm,i,&nlines,&connlines);CHKERRQ(ierr);
    if (nlines > MAXCONNLINES) SETERRQ3(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"%D lines connected to a bus %D exceeds max. allowed connected lines allowed %D",nlines,ps->bus[i-vStart].bus_i,MAXCONNLINES);
    ps->bus[i-vStart].nconnlines = nlines;
    for(k=0; k < nlines; k++) ps->bus[i-vStart].connlines[k] = &ps->line[connlines[k]-eStart];

    /* Is bus ghosted? */
    ierr = DMNetworkIsGhostVertex(ps->networkdm,i,&ps->bus[i-vStart].isghost);CHKERRQ(ierr);

    /* Starting location in the local array */
    ierr = DMNetworkGetVariableOffset(ps->networkdm,i,&ps->bus[i-vStart].startloc);CHKERRQ(ierr);

    /* Starting location in the global array */
    ierr = DMNetworkGetVariableGlobalOffset(ps->networkdm,i,&ps->bus[i-vStart].startlocglob);CHKERRQ(ierr);
    /* If the bus is a ghost bus then DMPlex (DMNetwork) has the global offset negative to indicate that it is a ghost bus.
       We need to convert it to the actual global offset so that matrix values can be set using it 
    */
    //    if(ps->bus[i-vStart].startlocglob < 0) ps->bus[i-vStart].startlocglob = -ps->bus[i-vStart].startlocglob - 1;

    /* Incident generators */
    PetscInt startloc=2*NPHASE;

    for(k=0; k < ps->bus[i-vStart].ngen; k++) {
      ps->bus[i-vStart].gens[k] = &ps->gen[ps->bus[i-vStart].gidx[k]];

      if(ps->app == APP_DYNSIM) {
	if(!ps->bus[i-vStart].gens[k]->status) continue;

	ierr = PSGENGetDYNGen(ps->bus[i-vStart].gens[k],&dyngen);CHKERRQ(ierr);
	ierr = PetscCalloc1(dyngen->nvar,&dyngen->eqtypes);CHKERRQ(ierr);
	ierr = DYNGenModelGetEquationTypes(dyngen,&dyngen->ndiff,&dyngen->nalg,dyngen->eqtypes);CHKERRQ(ierr);
	
	/* Update number of differential equations used in creating differential equations IS */
	if(!ps->bus[i-vStart].isghost) ps->ndiff += dyngen->ndiff;

	/* Starting location for the variables relative to the bus variables */
	dyngen->startloc = startloc;
	startloc += dyngen->nvar;

	/* Add a pointer to the bus */
	dyngen->bus = &ps->bus[i-vStart];

	/* Add a pointer to the gen */
	dyngen->psgen = dyngen->bus->gens[k];

	/* Exciter Model */
	if(ps->bus[i-vStart].gens[k]->hasexc) {
	  dyngen->dynexc = &ps->bus[i-vStart].gens[k]->dynexc;
	  
	  ierr = PSGENGetDYNExc(ps->bus[i-vStart].gens[k],&dynexc);CHKERRQ(ierr);
	  ierr = PetscCalloc1(dynexc->nvar,&dynexc->eqtypes);CHKERRQ(ierr);
	  ierr = DYNExcModelGetEquationTypes(dynexc,&dynexc->ndiff,&dynexc->nalg,dynexc->eqtypes);CHKERRQ(ierr);
	  /* Update number of differential equations */
	  if(!ps->bus[i-vStart].isghost) ps->ndiff += dynexc->ndiff;
	  
	  /* Starting location of the exciter variables relative to the bus variables */
	  dynexc->startloc = startloc;
	  startloc += dynexc->nvar;

	  dynexc->dyngen = dyngen;
	}
	  
	/* Turbine Governor Model */
	if(ps->bus[i-vStart].gens[k]->hasturbgov) {
	  dyngen->dynturbgov = &ps->bus[i-vStart].gens[k]->dynturbgov;
	  
	  ierr = PSGENGetDYNTurbgov(ps->bus[i-vStart].gens[k],&dynturbgov);CHKERRQ(ierr);
	  ierr = PetscCalloc1(dynturbgov->nvar,&dynturbgov->eqtypes);CHKERRQ(ierr);
	  ierr = DYNTurbgovModelGetEquationTypes(dynturbgov,&dynturbgov->ndiff,&dynturbgov->nalg,dynturbgov->eqtypes);CHKERRQ(ierr);
	  /* Update number of differential equations */
	  if(!ps->bus[i-vStart].isghost) ps->ndiff += dynturbgov->ndiff;
	  
	  /* Starting location of the turbine governor variables relative to the bus variables */
	  dynturbgov->startloc = startloc;
	  startloc += dynturbgov->nvar;

	  dynturbgov->dyngen = dyngen;
	}

	/* Stabilizer Model */
	if(ps->bus[i-vStart].gens[k]->hasstab && dyngen->dynexc) {
	  dyngen->dynexc->dynstab = &ps->bus[i-vStart].gens[k]->dynstab;
	  
	  ierr = PSGENGetDYNStab(ps->bus[i-vStart].gens[k],&dynstab);CHKERRQ(ierr);
	  ierr = PetscCalloc1(dynstab->nvar,&dynstab->eqtypes);CHKERRQ(ierr);
	  ierr = DYNStabModelGetEquationTypes(dynstab,&dynstab->ndiff,&dynstab->nalg,dynstab->eqtypes);CHKERRQ(ierr);
	  /* Update number of differential equations */
	  if(!ps->bus[i-vStart].isghost) ps->ndiff += dynstab->ndiff;
	  
	  /* Starting location of the turbine governor variables relative to the bus variables */
	  dynstab->startloc = startloc;
	  startloc += dynstab->nvar;

	  dynstab->dyngen = dyngen;
	}
      }
    }
    /* Change the bus type to PQ if no generators are incident */
    if(!ps->bus[i-vStart].ngenON && ps->bus[i-vStart].ide != ISOLATED_BUS) ps->bus[i-vStart].ide = PQ_BUS;
    
    /* Incident loads */
    for(k=0; k < ps->bus[i-vStart].nload; k++) {
      ps->bus[i-vStart].loads[k] = &ps->load[ps->bus[i-vStart].lidx[k]];

      if(ps->app == APP_DYNSIM) {
	if(!ps->bus[i-vStart].loads[k]->status) continue;

	ierr = PSLOADGetDYNLoad(ps->bus[i-vStart].loads[k],&dynload);CHKERRQ(ierr);
	ierr = PetscCalloc1(dynload->nvar,&dynload->eqtypes);CHKERRQ(ierr);
	if(dynload->nvar) {
	  ierr = DYNLoadModelGetEquationTypes(dynload,&dynload->ndiff,&dynload->nalg,dynload->eqtypes);CHKERRQ(ierr);
	}
	
	/* Update number of differential equations used in creating differential equations IS */
	if(!ps->bus[i-vStart].isghost) ps->ndiff += dynload->ndiff;

	/* Starting location for the variables relative to the bus variables */
	dynload->startloc = startloc;
	startloc += dynload->nvar;

	/* Add a pointer to the bus */
	dynload->bus = &ps->bus[i-vStart];

	/* Add a pointer to the load */
	dynload->psload = dynload->bus->loads[k];
      }
    }
  }
#if defined DEBUGPS
  ierr = PetscPrintf(PETSC_COMM_SELF,"Rank[%d]:nbuses = %d,nlines = %d,ngen = %d, nload = %d\n",ps->comm->rank,ps->nbus,ps->nbranch,ps->ngen,ps->nload);CHKERRQ(ierr);
  PSLINE line;
  for(i= 0; i < ps->nbranch; i++) {
    line = &ps->line[i];
    ierr = PetscPrintf(PETSC_COMM_SELF,"Rank[%d]:Line %d ----- %d\n",ps->comm->rank,line->fbus,line->tbus);CHKERRQ(ierr);
  }
#endif

  /* Broadcast MVAbase */
  ierr = MPI_Bcast(&ps->MVAbase,1,MPIU_SCALAR,0,ps->comm->type);CHKERRQ(ierr);

  //  ierr = PetscPrintf(PETSC_COMM_SELF,"Rank %d Came here\n",ps->comm->rank);CHKERRQ(ierr);
  ps->setupcalled = PETSC_TRUE;
  PetscFunctionReturn(0);
}

/*
  PSCreateGlobalVector - Returns a global vector of the appropriate size
  and distribution conforming to the distribution of the PS object.

  Input Paramereters:
. PS - the PS network object

  Output Parameters:
. vec - the global vector

  Notes:
  PSSetUp() must be called before calling this routine.
*/
PetscErrorCode PSCreateGlobalVector(PS ps,Vec *vec)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if(!ps->setupcalled) SETERRQ(ps->comm->type,0,"PSSetUp() must be called before calling PSCreateGlobalVector");
  ierr = DMCreateGlobalVector(ps->networkdm,vec);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  PSCreateMatrix - Returns a distributed matrix of appropriate size that can
   be used as the Jacobian


  Input Paramereters:
. PS - the PS network object

  Output Parameters:
. mat - the matrix

  Notes:
  PSSetUp() must be called before calling this routine.
*/
PetscErrorCode PSCreateMatrix(PS ps,Mat *mat)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if(!ps->setupcalled) SETERRQ(ps->comm->type,0,"PSSetUp() must be called before calling PSCreateMatrix");
  ierr = DMCreateMatrix(ps->networkdm,mat);CHKERRQ(ierr);
  ierr = MatSetOption(*mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* 
   PSGetConnectedComponents - Does a depth-first search to find the connected components starting from a given bus

   Input Parameters
.  start_bus - The starting bus/node

   Output Parameters
+  visited_nodes -- an array of visited nodes (if a node i is visited then visited_nodes[i] != 0, else it is 0)
-  nnodes_visited -- number of nodes visited (connected components)
*/
PetscErrorCode PSGetConnectedComponents(PSBUS start_bus, PetscInt *visited_nodes, PetscInt *nnodes_visited,PetscInt nrem)
{
  PetscErrorCode ierr;
  const PSLINE   *connlines;
  PSLINE         line;
  PetscInt       nconnlines;
  PSBUS          busf,bust,next_bus;
  PetscInt       i;
  const PSBUS    *connbuses;

  PetscFunctionBegin;

  if(visited_nodes[start_bus->internal_i] == 0) {
    visited_nodes[start_bus->internal_i] = 1;
    *nnodes_visited += 1;
  }
  if(*nnodes_visited == nrem) PetscFunctionReturn(0);

  ierr = PSBUSGetSupportingLines(start_bus,&nconnlines,&connlines);CHKERRQ(ierr);

  for(i=0; i < nconnlines; i++) {
    line = connlines[i];
    if(!line->status) continue;

    ierr = PSLINEGetConnectedBuses(line,&connbuses);CHKERRQ(ierr);
    busf = connbuses[0];
    bust = connbuses[1];

    if(start_bus == busf) next_bus = bust;
    else next_bus = busf;

    if(visited_nodes[next_bus->internal_i] == 0) {
      ierr =  PSGetConnectedComponents(next_bus, visited_nodes, nnodes_visited,nrem);CHKERRQ(ierr);
    }
  }

  PetscFunctionReturn(0);
}

/*
  PSCheckTopology - Checks the network topology and finds if islands exist

*/
PetscErrorCode PSCheckTopology(PS ps)
{
  PetscErrorCode ierr;
  PetscInt       *visited_nodes,nnodes_visited=0,nrem=ps->nbus,i,j;
  PSBUS          start_bus=0;
  PetscInt       Nconncomp=0;
  PetscInt       *scat_conncomp_sizes;
  PetscInt       *scatv_displs,*scatv_sendcounts;
  PetscInt       *scatv_sendbuf;
  PSConnCompgroup *cgroup;
  PetscInt *ncg_packet_glob,pctr=0,ii;
  PSConngroup grp;

  PetscFunctionBegin;
  //  if(ps->comm->size > 1) SETERRQ(PETSC_COMM_SELF,0,"No parallel support for finding islands yet\n");

  if(ps->nconncomp) {
    ierr = PSConnCompDestroy(ps);CHKERRQ(ierr);
  }

  ierr = PetscCalloc1(ps->nbus,&visited_nodes);CHKERRQ(ierr);

  while(nrem != 0) {
    nnodes_visited = 0;
    for(i=0;i < ps->nbus; i++) {
      if(visited_nodes[i] == 0) {
	start_bus = &ps->bus[i];
	break;
      }
    }
    ierr = PSGetConnectedComponents(start_bus,visited_nodes,&nnodes_visited,nrem);CHKERRQ(ierr);

    nrem -= nnodes_visited;
    ps->conncomp[ps->nconncomp].nv = 0;
    ierr = PetscCalloc1(nnodes_visited,&ps->conncomp[ps->nconncomp].v);CHKERRQ(ierr);
    for(i=0;i < ps->nbus; i++) {
      if(visited_nodes[i] == 1) {
	ps->conncomp[ps->nconncomp].v[ps->conncomp[ps->nconncomp].nv++] = ps->bus[i].bus_i;
	visited_nodes[i]++; /* Advance so that the node does not get included in the next iteration */
      }
    }
    ps->nconncomp++;
  }

  ierr = PetscFree(visited_nodes);CHKERRQ(ierr);

  /* At this stage, all processes know about the connected components for their subnetwork */
  /* Prepare to send size information to processor 0 */
  PetscInt ncg_packet_size=0; /* Size of the packet containing the size information */
  PetscInt *ncg_packet;
  PetscInt *tp;

  /* Each packet has the following layout */
  /* [Rank | # of connected componets | size of each connected component | vertices in each connected component] */

  for(i=0; i < ps->nconncomp; i++) {
    ncg_packet_size += ps->conncomp[i].nv + 1; /* Add 1 for its size info */
  }
  ncg_packet_size += 2; /* Add rank and number of connected componenets */

  ierr = PetscCalloc1(ncg_packet_size,&ncg_packet);CHKERRQ(ierr);
  tp = ncg_packet + 2 + ps->nconncomp;
  ncg_packet[0] = ps->comm->rank;
  ncg_packet[1] = ps->nconncomp;
  for(i=0; i < ps->nconncomp; i++) {
    ncg_packet[2+i] = ps->conncomp[i].nv;
    for(j=0; j < ps->conncomp[i].nv; j++) tp[j] = ps->conncomp[i].v[j];
    tp += ps->conncomp[i].nv;
  }

  /* Uncomment for debugging
  for(i=0; i < ncg_packet_size; i++) {
    if(i == ncg_packet_size - 1) {
      ierr = PetscPrintf(PETSC_COMM_SELF,"%d   \n",ncg_packet[i]);CHKERRQ(ierr);
    } else {
      ierr = PetscPrintf(PETSC_COMM_SELF,"%d   ",ncg_packet[i]);CHKERRQ(ierr);
    }
  }
  */
  /* Send packet size and number of connected components on each process to root */
  PetscInt sendbuf[2],*recvbuf;
  PetscInt *displs,*reccnt;
  if(!ps->comm->rank) {
    ierr = PetscCalloc1(2*ps->comm->size,&recvbuf);CHKERRQ(ierr);
    ierr = PetscCalloc1(ps->comm->size,&displs);CHKERRQ(ierr);
    ierr = PetscCalloc1(ps->comm->size,&reccnt);CHKERRQ(ierr);
  }

  sendbuf[0] = ncg_packet_size; /* Size of the packet sent for reduction */
  sendbuf[1] = ps->nconncomp; /* Number of connected components on each process */
  ierr = MPI_Gather(sendbuf,2,MPIU_INT,recvbuf,2,MPIU_INT,0,ps->comm->type);CHKERRQ(ierr);
  
  /* Gather packets from each process onto root */
  PetscInt ncg_packet_size_glob=0,nconncomp_glob=0;

  if(!ps->comm->rank) {
    for(i=0; i < ps->comm->size; i++) {
      displs[i] = ncg_packet_size_glob;
      ncg_packet_size_glob += recvbuf[2*i];
      nconncomp_glob += recvbuf[2*i+1];
      reccnt[i] = recvbuf[2*i];
    }
    /* Uncomment for debugging
    ierr = PetscPrintf(PETSC_COMM_SELF,"global packet size = %d, total number of connected components = %d\n",ncg_packet_size_glob,nconncomp_glob);CHKERRQ(ierr);
    */
    ierr = PetscCalloc1(ncg_packet_size_glob,&ncg_packet_glob);CHKERRQ(ierr);
  }
  
  ierr = MPI_Gatherv(ncg_packet,ncg_packet_size,MPIU_INT,ncg_packet_glob,reccnt,displs,MPIU_INT,0,ps->comm->type);CHKERRQ(ierr);

  if(!ps->comm->rank) {
    ierr = PetscFree(displs);CHKERRQ(ierr);
    ierr = PetscFree(reccnt);CHKERRQ(ierr);

    /* Uncomment for debugging
    for(i=0; i < ncg_packet_size_glob; i++) {
      ierr = PetscPrintf(PETSC_COMM_SELF,"%d   ",ncg_packet_glob[i]);CHKERRQ(ierr);
    }
    */


    /* Create conncompgroup structs. The number of elements equals to the
       total number of connected components. Each struct will have the rank,
       size of the connected component, and the connected component info. The
       idea is to think of each element of this connected component group struct
       as a vertex on a graph. A vertex is connected to another vertex if they have
       common nodes in the 'c' field of the struct.
    */

    PetscInt *ptr=ncg_packet_glob;
    PetscInt kk,nconncomp_proc,*data_ptr,rank;
    ierr = PetscCalloc1(nconncomp_glob,&cgroup);CHKERRQ(ierr);
    for(i=0; i < ps->comm->size; i++) {
      rank = *ptr;
      nconncomp_proc = *(ptr+1);
      ptr += 2;
      data_ptr = ptr+nconncomp_proc;
      for(kk=0;kk < nconncomp_proc; kk++) { 
	cgroup[pctr].rank = rank;
	cgroup[pctr].nc = *ptr++;
	ierr = PetscCalloc1(cgroup[pctr].nc,&cgroup[pctr].c);CHKERRQ(ierr);
	ierr = PetscMemcpy(cgroup[pctr].c,data_ptr,cgroup[pctr].nc*sizeof(PetscInt));CHKERRQ(ierr);
	data_ptr += cgroup[pctr].nc;
	pctr++;
      }
      ptr = data_ptr;
    }

    /*
    for(i=0; i < pctr; i++) {
      ierr = PetscPrintf(PETSC_COMM_SELF,"ctr = %d rank = %d\n",i,cgroup[i].rank);CHKERRQ(ierr);
      for(kk=0; kk < cgroup[i].nc; kk++) {
	ierr = PetscPrintf(PETSC_COMM_SELF,"%d  ",cgroup[i].c[kk]);CHKERRQ(ierr);
      }
      ierr = PetscPrintf(PETSC_COMM_SELF,"\n");CHKERRQ(ierr);
    }
    */

    PSConngroupi *gi;
    PetscInt *visited_nodegroups;
    PetscInt nrem=nconncomp_glob,k,nvis=0;
    PSConnCompgroup *start_group;

    nrem = nconncomp_glob;
    grp.n = 0;
    ierr = PetscCalloc1(nconncomp_glob,&grp.ci);CHKERRQ(ierr);
    ierr = PetscCalloc1(nconncomp_glob,&visited_nodegroups);CHKERRQ(ierr);

    while(nrem) {
      gi = &grp.ci[grp.n];
      gi->nc = 0;
      nvis = 0;
      ierr = PetscCalloc1(nconncomp_glob,&gi->cg);CHKERRQ(ierr);
      for(i=0; i < nconncomp_glob; i++) {
	if(visited_nodegroups[i] == 0) {
	  visited_nodegroups[i] += 1;
	  start_group = &cgroup[i];
	  gi->cg[gi->nc++] = start_group;
	  nvis++;
	  break;
	}
      }
      
      PetscInt new_node=1;
      PetscInt l;
      PSConnCompgroup *cgi;
      
      while(new_node) {
	new_node = 0;
	for(i=0; i < nconncomp_glob; i++) {
	  if(visited_nodegroups[i]) continue;
	  for(j=0; j < cgroup[i].nc; j++) {
	    //	    for(k=0; k < start_group->nc; k++) {
	    for(k=0; k < gi->nc; k++) {
	      cgi = gi->cg[k];
	      for(l=0; l < cgi->nc; l++) {
	      //	    ierr = PetscPrintf(PETSC_COMM_SELF,"start_group->c[%d] = %d  cgroup[%d].c[%d] = %d\n",k,start_group->c[k],i,j,cgroup[i].c[j]);CHKERRQ(ierr);
		if(cgi->c[l] == cgroup[i].c[j]) {
		//	      ierr = PetscPrintf(PETSC_COMM_SELF,"Came here\n");CHKERRQ(ierr);
		  visited_nodegroups[i] += 1;
		  gi->cg[gi->nc++] = &cgroup[i];
		  nvis++;
		  /* Break from the loop */
		  l = cgi->nc;
		  j = cgroup[i].nc;
		  k = gi->nc;
		  new_node = 1;
		}
	      }
	    }
	  }
	}
      }
      nrem -= nvis;
      grp.n++;
    }



    /* Uncomment for debugging
    for(i=0; i < grp.n; i++) {
      gi = &grp.ci[i];
      ierr = PetscPrintf(PETSC_COMM_SELF,"Group number = %d\n",i);CHKERRQ(ierr);
      for(j=0; j < gi->nc; j++) {
	for(k=0; k < gi->cg[j]->nc; k++) {
	  ierr = PetscPrintf(PETSC_COMM_SELF,"  %d  ",gi->cg[j]->c[k]);CHKERRQ(ierr);
	}
      }
      ierr = PetscPrintf(PETSC_COMM_SELF,"\n");CHKERRQ(ierr);
    }
    */

    /* Create buffer for scattering sizes */
    ierr = PetscCalloc1(ps->comm->size*(grp.n+1),&scat_conncomp_sizes);CHKERRQ(ierr);
    for(i=0; i < ps->comm->size; i++) {
      scat_conncomp_sizes[(grp.n+1)*i] = grp.n;
    }
    for(i=0; i < grp.n; i++) {
      gi = &grp.ci[i];
      for(j=0; j < gi->nc; j++) {
       	scat_conncomp_sizes[(gi->cg[j]->rank)*(grp.n+1)+1+i] += gi->cg[j]->nc;
      }
    }

    ierr = PetscCalloc1(ps->comm->size,&scatv_displs);CHKERRQ(ierr);
    ierr = PetscCalloc1(ps->comm->size,&scatv_sendcounts);CHKERRQ(ierr);
    scatv_displs[0] = 0;
    PetscInt scatv_sendbuf_size = 0;
    for(i=0; i < ps->comm->size; i++) {
      if(i >= 1) scatv_displs[i] = scatv_displs[i-1] + scatv_sendcounts[i-1];
      for(j=0; j < grp.n+1; j++) {
	if(j >= 1) scatv_sendcounts[i] += scat_conncomp_sizes[(grp.n+1)*i+j];
	//ierr = PetscPrintf(PETSC_COMM_SELF,"%d  ",scat_conncomp_sizes[(grp.n+1)*i+j]);CHKERRQ(ierr);
      }
      //ierr = PetscPrintf(PETSC_COMM_SELF,"\n");CHKERRQ(ierr);
      scatv_sendbuf_size += scatv_sendcounts[i];
    }
    //ierr = PetscPrintf(PETSC_COMM_SELF,"scatv_sendbuf_size = %d\n",scatv_sendbuf_size);CHKERRQ(ierr);
    
    ierr = PetscCalloc1(scatv_sendbuf_size,&scatv_sendbuf);CHKERRQ(ierr);
    PetscInt **displs_ptr;
    ierr = PetscCalloc1(ps->comm->size,&displs_ptr);CHKERRQ(ierr);
    for(i=0; i < ps->comm->size;i++) {
      /* Initialize displacement pointers to mark to the beginning of the array for each process */
      displs_ptr[i] = scatv_sendbuf + scatv_displs[i];
    }
    for(i=0; i < grp.n; i++) {
      gi = &grp.ci[i];
      for(j=0; j < gi->nc; j++) {
	ierr = PetscMemcpy(displs_ptr[gi->cg[j]->rank],gi->cg[j]->c,gi->cg[j]->nc*sizeof(PetscInt));CHKERRQ(ierr);
	displs_ptr[gi->cg[j]->rank] += gi->cg[j]->nc;
      }
    } 
    ierr = PetscFree(displs_ptr);CHKERRQ(ierr);

    /* Uncomment for debugging
    for(i=0; i < scatv_sendbuf_size; i++) {
      ierr = PetscPrintf(PETSC_COMM_SELF,"  %d  ",scatv_sendbuf[i]);CHKERRQ(ierr);
    }
    ierr = PetscPrintf(PETSC_COMM_SELF,"\n");CHKERRQ(ierr);
    */
    /*
    for(i=0; i < nconncomp_glob; i++) {
      ierr = PetscPrintf(PETSC_COMM_SELF,"Rank = %d, Ncomp = %d  ",cgroup[i].rank,cgroup[i].nc);CHKERRQ(ierr);
      for(kk=0; kk < cgroup[i].nc;kk++) {
	ierr = PetscPrintf(PETSC_COMM_SELF,"%d  ",cgroup[i].c[kk]);CHKERRQ(ierr);
      }
      ierr = PetscPrintf(PETSC_COMM_SELF,"\n");CHKERRQ(ierr);
    }
    */
    /*
    ncg_packet[0] = ps->comm->rank;
    ncg_packet[1] = ps->nconncomp;
    PetscInt *tp=ncg_packet + 2 + ps->nconncomp;
    for(i=0; i < ps->nconncomp; i++) {
      ncg_packet[2+i] = ps->conncomp[i].nv;
      for(j=0; j < ps->conncomp[i].nv; j++) tp[j] = ps->conncomp[i].v[j];
      tp += ps->conncomp[i].nv;
    }
    */

    Nconncomp = grp.n;

    ierr = PetscFree(visited_nodegroups);CHKERRQ(ierr);
  }
  
  MPI_Barrier(ps->comm->type);

  ierr = PSConnCompDestroy(ps);CHKERRQ(ierr);
  ps->nconncomp = Nconncomp;
  ierr = MPI_Bcast(&ps->nconncomp,1,MPIU_INT,0,ps->comm->type);CHKERRQ(ierr);
  //ierr = PetscPrintf(PETSC_COMM_SELF,"Rank[%d] nconncomp=%d\n",ps->comm->rank,ps->nconncomp);CHKERRQ(ierr);
  
  PetscInt *conncomp_sizes;
  ierr = PetscCalloc1(ps->nconncomp+1,&conncomp_sizes);CHKERRQ(ierr);
  ierr = MPI_Scatter(scat_conncomp_sizes,ps->nconncomp+1,MPIU_INT,conncomp_sizes,ps->nconncomp+1,MPIU_INT,0,ps->comm->type);CHKERRQ(ierr);
  for(i=0; i < ps->nconncomp; i++) {
    ps->conncomp[i].nv = conncomp_sizes[i+1];
    //ierr = PetscPrintf(PETSC_COMM_SELF,"Rank[%d] nv[%d] = %d\n",ps->comm->rank,i,ps->conncomp[i].nv);CHKERRQ(ierr);
    ierr = PetscCalloc1(ps->conncomp[i].nv,&ps->conncomp[i].v);CHKERRQ(ierr);    
  }  
  ierr = PetscFree(conncomp_sizes);CHKERRQ(ierr);

  PetscInt *scatv_recvbuf;
  ierr = PetscCalloc1(ps->nbus,&scatv_recvbuf);CHKERRQ(ierr);

  ierr = MPI_Scatterv(scatv_sendbuf,scatv_sendcounts,scatv_displs,MPIU_INT,scatv_recvbuf,ps->nbus,MPIU_INT,0,ps->comm->type);CHKERRQ(ierr);

  PetscInt *ptr=scatv_recvbuf;
  for(i=0; i < ps->nconncomp; i++) {
    ierr = PetscMemcpy(ps->conncomp[i].v,ptr,ps->conncomp[i].nv*sizeof(PetscInt));CHKERRQ(ierr);
    for(j=0; j < ps->conncomp[i].nv; j++) {
      //ierr = PetscPrintf(PETSC_COMM_SELF,"Rank[%d] Group[%d]  Node = %d\n",ps->comm->rank,i,ps->conncomp[i].v[j]);CHKERRQ(ierr);
    }
    ptr += ps->conncomp[i].nv;
  }

  if(!ps->comm->rank) {
    ierr = PetscFree(scatv_displs);CHKERRQ(ierr);
    ierr = PetscFree(scatv_sendcounts);CHKERRQ(ierr);
    ierr = PetscFree(scatv_sendbuf);CHKERRQ(ierr);
    ierr = PetscFree(scat_conncomp_sizes);CHKERRQ(ierr);
    ierr = PetscFree(ncg_packet_glob);CHKERRQ(ierr);
    ierr = PetscFree(recvbuf);CHKERRQ(ierr);

    for(ii=0; ii < pctr; ii++) {
      ierr = PetscFree(cgroup[ii].c);CHKERRQ(ierr);
    }
    ierr = PetscFree(cgroup);CHKERRQ(ierr);

    for(ii=0; ii < grp.n; ii++) {
      ierr = PetscFree(grp.ci[ii].cg);CHKERRQ(ierr);
    }
    ierr = PetscFree(grp.ci);CHKERRQ(ierr);
  }
  ierr = PetscFree(ncg_packet);CHKERRQ(ierr);
  ierr = PetscFree(scatv_recvbuf);CHKERRQ(ierr);

    
  /* Check if the island has an active generator */
  for(i=0; i < ps->nconncomp; i++) {
    PetscInt islandhasgen = 0,k,islandhasgen_glob;
    PSBUS    bus;
    for(k=0; k < ps->conncomp[i].nv; k++) {
      bus = &ps->bus[ps->busext2intmap[ps->conncomp[i].v[k]]];
      if(bus->ngenON) {
	islandhasgen = 1;
	ps->conncomp[i].blackout = 0;
	break;
      }
    }
    
    ierr = MPI_Allreduce(&islandhasgen,&islandhasgen_glob,1,MPIU_INT,MPI_SUM,ps->comm->type);CHKERRQ(ierr);

    if(!islandhasgen_glob) { /* No generators on this island, set all buses isolated on this island */
      PSLOAD load;
      PSLINE line;
      PetscInt nconnlines;
      const PSLINE *connlines;
      PetscInt n;

      //ierr = PetscPrintf(PETSC_COMM_WORLD,"Island %d black out\n",i+1);CHKERRQ(ierr);
      ps->conncomp[i].blackout = 1;
      for(k=0; k < ps->conncomp[i].nv; k++) {
	bus = &ps->bus[ps->busext2intmap[ps->conncomp[i].v[k]]];
	bus->ide = ISOLATED_BUS;
	/* Turn OFF connected lines */
	ierr = PSBUSGetSupportingLines(bus,&nconnlines,&connlines);CHKERRQ(ierr);
	for(n=0; n < nconnlines; n++) {
	  line = connlines[n];
	  line->status = 0;
	}
	
	/* Turn OFF connected loads */
	for(n=0; n < bus->nload; n++) {
	  ierr = PSBUSGetLoad(bus,n,&load);CHKERRQ(ierr);
	  load->status = 0;
	}
      }
    }
  }

  /*
  if(ps->nconncomp > 1) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Found %d islands in the network\n",ps->nconncomp);CHKERRQ(ierr);
    for(i=0; i < ps->nconncomp; i++) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Buses in island %d:\n",i+1);CHKERRQ(ierr);
      for(j=0; j < ps->conncomp[i].nv; j++) {
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Bus %d\n",ps->conncomp[i].v[j]);CHKERRQ(ierr);
      }
    }
  }
  */

  for(i=0; i < ps->nconncomp; i++) {
    ierr = PSIslandCheckandSetRefBus(ps,i);CHKERRQ(ierr);
  }

  ierr = PetscFree(visited_nodes);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  PSSetGenStatus - Sets generator status given the bus name and id
*/
PetscErrorCode PSSetGenStatus(PS ps,PetscInt gbus,const char* gid,PetscInt status)
{
  PetscErrorCode ierr;
  PSGEN          gen;
  PSBUS          bus;
  PetscBool       statusupdate=PETSC_FALSE;

  PetscFunctionBegin;
  if(!ps->setupcalled) SETERRQ(PETSC_COMM_SELF,0,"PSSetUp() must be called before calling PFLOWGetGen()\n");

  ierr = PSGetGen(ps,gbus,gid,&gen);CHKERRQ(ierr);
  
  if(gen) {
    if(status != gen->status) statusupdate = PETSC_TRUE;
    ierr = PSGENSetStatus(gen,status);CHKERRQ(ierr);
    if(statusupdate) {
      bus = &ps->bus[ps->busext2intmap[gbus]];
      if(status) {
	bus->ngenON++;
      } else {
	bus->ngenON--;
	/* Change the bus type to PQ if no generators are active at this bus */
	if(!bus->ngenON) bus->ide = PQ_BUS;
	gen->pg = gen->qg = 0.0;
      }
    }
  }

  PetscFunctionReturn(0);
}

/*
  PSSetLineStatus - Sets the line status given from bus, to bus, and line id
*/
PetscErrorCode PSSetLineStatus(PS ps,PetscInt fbus, PetscInt tbus, const char* id,PetscInt status)
{
  PetscErrorCode ierr;
  PSLINE         line;

  PetscFunctionBegin;
  if(!ps->setupcalled) SETERRQ(PETSC_COMM_SELF,0,"PSSetUp() must be called before calling PFLOWGetLine()\n");

  ierr = PSGetLine(ps,fbus,tbus,id,&line);CHKERRQ(ierr);
  
  if(line) {
    ierr = PSLINESetStatus(line,status);CHKERRQ(ierr);
  }
  
  PetscFunctionReturn(0);
}
