#include <private/psimpl.h>

const char *const PSGENFuelTypes[] = {"COAL",      "WIND",         "SOLAR",
                                      "NG",        "NUCLEAR",      "HYDRO",
                                      "UNDEFINED", "GENFUELTYPE_", NULL};

/*
  PSGetKVLevels - Gets the different KV levels in the system

  Input Parameters
. ps - PS struct

  Output Parameters
nkv - number of kv levels found
kv  - kv levels
*/
PetscErrorCode PSGetKVLevels(PS ps, PetscInt *nkv, const PetscScalar **kv) {
  PetscInt i, j;
  PSBUS bus;
  PetscScalar buskV;
  PetscBool kvlevelfound;

  PetscFunctionBegin;

  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    buskV = bus->basekV;
    kvlevelfound = PETSC_FALSE;
    for (j = 0; j < ps->nkvlevels; j++) {
      if (PetscAbsScalar(buskV - ps->kvlevels[j]) < 1e-6)
        kvlevelfound = PETSC_TRUE;
    }
    if (!kvlevelfound) {
      ps->kvlevels[ps->nkvlevels++] = buskV;
    }
  }

  *nkv = ps->nkvlevels;
  *kv = ps->kvlevels;
  PetscFunctionReturn(0);
}
/*
     PSCheckandSetRefBus - Checks for active ref. bus and sets one if not
available

  Input Parameters:
. ps - The PS object

  This routine checks if there is an active ref. bus (declared in the data file
and has active generators). If it is not defined, then the first PV bus is used.
*/
PetscErrorCode PSCheckandSetRefBus(PS ps) {
  PetscErrorCode ierr;
  PetscInt i;
  PSBUS bus;
  PetscBool pshaslocalref = PETSC_FALSE, pshasref;
  PetscInt firstlocpvbus = 100000000, firstpvbus;

  PetscFunctionBegin;
  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    if (bus->ide == REF_BUS && bus->ngenON) {
      /* PS has active reference bus */
      pshaslocalref = PETSC_TRUE;
      break;
    }
  }

  ierr = MPI_Allreduce(&pshaslocalref, &pshasref, 1, MPIU_BOOL, MPI_LOR,
                       ps->comm->type);
  CHKERRQ(ierr);

  if (pshasref)
    PetscFunctionReturn(0);
#if defined DEBUGPS
  ierr = PetscPrintf(PETSC_COMM_SELF, "No active ref. bus in the system\n");
  CHKERRQ(ierr);
#endif
  /* No active ref. bus, set first PV bus to ref bus */
  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    if (bus->ide == PV_BUS && bus->ngenON && !bus->isghost) {
      firstlocpvbus = bus->bus_i;
      break;
    }
  }

  ierr = MPIU_Allreduce(&firstlocpvbus, &firstpvbus, 1, MPIU_INT, MPI_MIN,
                        ps->comm->type);
  CHKERRQ(ierr);

  if (firstpvbus == 100000000) {
#if defined DEBUGPS
    SETERRQ(PETSC_COMM_SELF, 0, "Using first PV bus = %d", firstpvbus);
#endif
  } else {
    if (ps->busext2intmap[firstpvbus] != -1) {
      ps->bus[ps->busext2intmap[firstpvbus]].ide = REF_BUS;
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
PetscErrorCode PSSetGenDispatchandStatus(PS ps, PetscInt busnum,
                                         PetscInt gennum, PetscInt status,
                                         PetscScalar Pg, PetscScalar Qg) {
  PetscInt intbusnum;
  PSBUS bus;
  PSGEN gen;
  PetscBool statusupdate = PETSC_FALSE;

  PetscFunctionBegin;

  //  if(ps->setupcalled) {
  //    SETERRQ(PETSC_COMM_SELF,0,"Must call PSSetGenDispatchandStatus before
  //    PSSetUp() is called");
  //  }
  intbusnum = ps->busext2intmap[busnum];
  if (intbusnum == -1)
    PetscFunctionReturn(0); /* Not local to this processor */
  bus = &ps->bus[intbusnum];
  if (!bus->ngen) {
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "No generators on bus %d", busnum);
  }
  if (bus->ide == ISOLATED_BUS) {
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "Bus %d is isolated", busnum);
  }

  gen = &ps->gen[bus->gidx[gennum]];
  if (status != gen->status)
    statusupdate = PETSC_TRUE;
  gen->status = status;
  gen->pg = Pg / ps->MVAbase;
  gen->qg = Qg / ps->MVAbase;

  if (statusupdate) {
    if (status) {
      bus->ngenON++;
      ps->ngenON++;
      ps->NgenON++;
      if (bus->ide == PQ_BUS)
        bus->ide = PV_BUS;
    } else {
      bus->ngenON--;
      ps->ngenON--;
      ps->NgenON--;
      if (!bus->ngenON)
        bus->ide = PQ_BUS;
      gen->pg = gen->qg = 0.0;
    }
  }

  PetscFunctionReturn(0);
}

/*
  PSGetNumGenerators - Gets the number of local and global generators in this
system

  Input Parameters
. ps - the PS object

  Output Parameters
+ ngen - number of local generators
- Ngen - number of global generators

  Notes:
   PSSetUp() must be called before a call to PSGetNumGenerators
*/
PetscErrorCode PSGetNumGenerators(PS ps, PetscInt *ngen, PetscInt *Ngen) {
  PetscFunctionBegin;
  if (ngen)
    *ngen = ps->ngen;
  if (Ngen)
    *Ngen = ps->Ngen;
  PetscFunctionReturn(0);
}

/*
  PSGetNumActiveGenerators - Gets the number of local and global active
generators (status ON) in this system

  Input Parameters
. ps - the PS object

  Output Parameters
+ ngenON - number of local generators
- NgenON - number of global generators

  Notes:
   PSSetUp() must be called before a call to PSGetNumActiveGenerators
*/
PetscErrorCode PSGetNumActiveGenerators(PS ps, PetscInt *ngenON,
                                        PetscInt *NgenON) {
  PetscFunctionBegin;
  if (ngenON)
    *ngenON = ps->ngenON;
  if (NgenON)
    *NgenON = ps->NgenON;
  PetscFunctionReturn(0);
}

/*
  PSGetNumActiveLines - Gets the number of local and global active lines (status
ON) in this system

  Input Parameters
. ps - the PS object

  Output Parameters
+ nlineON - number of local lines with STATUS ON
- NlineON - number of global lines with STATUS ON

  Notes:
   PSSetUp() must be called before a call to PSGetNumActiveLines
*/
PetscErrorCode PSGetNumActiveLines(PS ps, PetscInt *nlineON,
                                   PetscInt *NlineON) {
  PetscFunctionBegin;
  if (nlineON)
    *nlineON = ps->nlineON;
  if (NlineON)
    *NlineON = ps->NlineON;
  PetscFunctionReturn(0);
}

/*
  PSLINESetStatus - Sets the status of the line

  Input Parameters:
+ line - the line
- status   - (0 or 1) the status of the line
*/
PetscErrorCode PSLINESetStatus(PSLINE line, PetscInt status) {
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
PetscErrorCode PSGetLine(PS ps, PetscInt fbus, PetscInt tbus, const char *id,
                         PSLINE *line) {
  PetscErrorCode ierr;
  PetscInt nsupplines;
  const PSLINE *supplines;
  PSBUS bus;
  PetscInt i;
  PetscBool flg = PETSC_FALSE;

  PetscFunctionBegin;
  if (ps->busext2intmap[fbus] != -1) {
    bus = &ps->bus[ps->busext2intmap[fbus]];
    ierr = PSBUSGetSupportingLines(bus, &nsupplines, &supplines);
    CHKERRQ(ierr);
    for (i = 0; i < nsupplines; i++) {
      if ((supplines[i]->fbus == fbus && supplines[i]->tbus == tbus) ||  (supplines[i]->reversed_ends & (supplines[i]->fbus == tbus && supplines[i]->tbus == fbus))) {
	ierr = PetscStrcmp(id, supplines[i]->ckt, &flg);
        CHKERRQ(ierr);
        if (flg) {
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
  PSLINEGetVariableGlobalLocation - Gets the starting location for the variables
for this line in the global vector

  Input Parameters:
. line - the line

  Output Parameter:
. locglob - the starting location for the variables in the global vector
*/
PetscErrorCode PSLINEGetVariableGlobalLocation(PSLINE line, PetscInt *locglob) {

  PetscFunctionBegin;
  *locglob = line->startlocglob;
  PetscFunctionReturn(0);
}

/*
  PSLINEGetVariableLocation - Gets the starting location for the variables for
this line in the local vector

  Input Parameters:
. line - the line object

  Output Parameter:
. loc - the starting location for the variables in the local vector
*/
PetscErrorCode PSLINEGetVariableLocation(PSLINE line, PetscInt *loc) {

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
PetscErrorCode PSLINEGetConnectedBuses(PSLINE line, const PSBUS **connbuses) {
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
PetscErrorCode PSBUSGetSupportingLines(PSBUS bus, PetscInt *nsupplines,
                                       const PSLINE **supplines) {
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
PetscErrorCode PSBUSAddShunt(PSBUS bus, PetscScalar Gs, PetscScalar Bs) {
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
PetscErrorCode PSBUSGetNGen(PSBUS bus, PetscInt *ngen) {
  PetscFunctionBegin;
  *ngen = bus->ngen;
  PetscFunctionReturn(0);
}

/*
  PSBUSGetNLoad - Gets the number of loads incident at the bus

  Input Parameters:
. bus  - the bus

  Output Parameters:
. nload - number of loads incident at the bus
*/
PetscErrorCode PSBUSGetNLoad(PSBUS bus, PetscInt *nload) {
  PetscFunctionBegin;
  *nload = bus->nload;
  PetscFunctionReturn(0);
}

/*
  PSGetNumLoads - Gets the number of local and global loads in this system

  Input Parameters
. ps - the PS object

  Output Parameters
+ nload - number of local loads
- Nload - number of global loads

  Notes:
   PSSetUp() must be called before a call to PSGetNumLoads
*/
PetscErrorCode PSGetNumLoads(PS ps, PetscInt *nload, PetscInt *Nload) {
  PetscFunctionBegin;
  if (nload)
    *nload = ps->nload;
  if (Nload)
    *Nload = ps->Nload;
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
PetscErrorCode PSBUSGetGen(PSBUS bus, PetscInt gen_num, PSGEN *gen) {

  PetscFunctionBegin;
  if (gen_num < 0 || gen_num > bus->ngen - 1)
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP,
            "Ngen at bus = %d, generator number requested %d", bus->ngen,
            gen_num);
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
PetscErrorCode PSBUSGetLoad(PSBUS bus, PetscInt load_num, PSLOAD *load) {

  PetscFunctionBegin;
  if (load_num < 0 || load_num > bus->nload - 1)
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP,
            "Nload at bus = %d, load number requested %d", bus->nload,
            load_num);
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
. gen - the PSGEN object (or NULL if not found)

*/
PetscErrorCode PSGetGen(PS ps, PetscInt gbus, const char *id, PSGEN *gen) {
  PetscErrorCode ierr;
  PSBUS bus;
  PetscInt i;
  PSGEN busgen;
  PetscBool flg = PETSC_FALSE;

  PetscFunctionBegin;
  if (ps->busext2intmap[gbus] != -1) {
    bus = &ps->bus[ps->busext2intmap[gbus]];
    for (i = 0; i < bus->ngen; i++) {
      ierr = PSBUSGetGen(bus, i, &busgen);
      CHKERRQ(ierr);
      ierr = PetscStrcmp(id, busgen->id, &flg);
      CHKERRQ(ierr);
      if (flg) {
        *gen = busgen;
        PetscFunctionReturn(0);
      }
    }
  }
  *gen = NULL;

  PetscFunctionReturn(0);
}

/*
  PSGetLoad - Gets the load object given the load bus, and id

  Input Parameters
+ ps - the PS object
. lbus - load bus
- id   - load id

  Output Parameters
. load - the PSLOAD object (or NULL if not found)

*/
PetscErrorCode PSGetLoad(PS ps, PetscInt lbus, const char *id, PSLOAD *load) {
  PetscErrorCode ierr;
  PSBUS bus;
  PetscInt i;
  PSLOAD busload;
  PetscBool flg = PETSC_FALSE;

  PetscFunctionBegin;
  if (ps->busext2intmap[lbus] != -1) {
    bus = &ps->bus[ps->busext2intmap[lbus]];
    for (i = 0; i < bus->nload; i++) {
      ierr = PSBUSGetLoad(bus, i, &busload);
      CHKERRQ(ierr);
      ierr = PetscStrcmp(id, busload->id, &flg);
      CHKERRQ(ierr);
      if (flg) {
        *load = busload;
        PetscFunctionReturn(0);
      }
    }
  }
  *load = NULL;

  PetscFunctionReturn(0);
}

/*
  PSGENSetStatus - Sets the generator status (1 = ON, 0 = OFF)

  Inputs:
+ psgen - the PSGEN object
- status - the generator status
*/
PetscErrorCode PSGENSetStatus(PSGEN psgen, PetscInt status) {
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
PetscErrorCode PSLOADSetStatus(PSLOAD psload, PetscInt status) {
  PetscFunctionBegin;
  psload->status = status;
  PetscFunctionReturn(0);
}

/*
  PSBUSSetGenStatus - Sets the status of the generator
*/
PetscErrorCode PSBUSSetGenStatus(PSBUS bus, char gid[], PetscInt status) {
  PetscErrorCode ierr;
  PetscInt i;
  PSGEN gen = NULL;
  PetscBool flg;
  PetscBool statusupdate;

  PetscFunctionBegin;
  if (!bus->ngen)
    PetscFunctionReturn(0);

  for (i = 0; i < bus->ngen; i++) {
    ierr = PSBUSGetGen(bus, i, &gen);
    CHKERRQ(ierr);
    ierr = PetscStrcmp(gid, gen->id, &flg);
    CHKERRQ(ierr);
    if (flg)
      break;
  }

  if (status != gen->status)
    statusupdate = PETSC_TRUE;
  ierr = PSGENSetStatus(gen, status);
  CHKERRQ(ierr);
  if (statusupdate) {
    if (status) {
      bus->ngenON++;
    } else {
      bus->ngenON--;
      /* Change the bus type to PQ if no generators are active at this bus */
      if (!bus->ngenON)
        bus->ide = PQ_BUS;
      gen->pg = gen->qg = 0.0;
    }
  }
  PetscFunctionReturn(0);
}

/*
  PSBUSGetVariableGlobalLocation - Gets the starting location for the variables
for this bus in the global vector

  Input Parameters:
. bus - the bus

  Output Parameter:
. locglob - the starting location for the variables in the global vector
*/
PetscErrorCode PSBUSGetVariableGlobalLocation(PSBUS bus, PetscInt *locglob) {

  PetscFunctionBegin;
  *locglob = bus->startlocglob;
  PetscFunctionReturn(0);
}

/*
  PSBUSGetVariableLocation - Gets the starting location for the variables for
this bus in the local vector

  Input Parameters:
. bus - the bus

  Output Parameter:
. loc - the starting location for the variables in the local vector
*/
PetscErrorCode PSBUSGetVariableLocation(PSBUS bus, PetscInt *loc) {

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
PetscErrorCode PSBUSIsGhosted(PSBUS bus, PetscBool *ghostbus) {

  PetscFunctionBegin;
  *ghostbus = bus->isghost;
  PetscFunctionReturn(0);
}

/*
  PSGENDestroy - Destroys the GEN data struct in PS

  Input Parameters
. ps - the PS object
*/
PetscErrorCode PSGENDestroy(PS ps) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = PetscFree(ps->gen);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  PSLOADDestroy - Destroys the LOAD data struct in PS

  Input Parameters
. ps - the PS object
*/
PetscErrorCode PSLOADDestroy(PS ps) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = PetscFree(ps->load);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  PSSetApplication - Sets the type of application to be run on the PS object

  Input Parameters:
+ PS - the PS object
. psapp - the application object
- psappname - the application name (PFLOW,OPFLOW)
*/
PetscErrorCode PSSetApplication(PS ps, void *psapp, PSApp psappname) {
  PetscBool app_found;
  PetscFunctionBegin;
  ps->app = psapp;
  ps->appname = psappname;
  app_found = (PetscBool)(ps->appname == APP_ACPF || ps->appname == APP_ACOPF);
  if (!app_found)
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "Application not supported");
  PetscFunctionReturn(0);
}

/*
  PSIncreaseReferenceCount - Increases the reference count of PS to indicate
that is being shared by other objects

  Input Parameter
. ps - The PS object
*/
PetscErrorCode PSIncreaseReferenceCount(PS ps) {
  PetscFunctionBegin;
  ps->refct++;
  PetscFunctionReturn(0);
}

/*
  PSDecreaseReferenceCount - Decreases the reference count of PS

  Input Parameter
. ps - The PS object
*/
PetscErrorCode PSDecreaseReferenceCount(PS ps) {
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
PetscErrorCode PSCreate(MPI_Comm mpicomm, PS *psout) {
  PetscErrorCode ierr;
  PS ps;

  PetscFunctionBegin;
  ierr = PetscCalloc1(1, &ps);
  CHKERRQ(ierr);

  /* Create communicator */
  ierr = COMMCreate(mpicomm, &ps->comm);
  CHKERRQ(ierr);

  ps->MVAbase = 0.0;
  ps->Nbus = -1;
  ps->Ngen = -1;
  ps->NgenON = -1;
  ps->Nline = -1;
  ps->nlineON = -1;
  ps->NlineON = -1;
  ps->Nload = -1;
  ps->refct = 0;
  ps->app = NULL;
  ps->appname = APP_NONE;
  ps->ndiff = 0;
  ps->nconncomp = 0;
  ps->nref = ps->Nref = 0;
  ps->ngencoal = ps->ngenwind = ps->ngensolar = 0;
  ps->ngenng = ps->ngennuclear = ps->ngenundefined = 0;

  ps->nkvlevels = 0;
  ierr = PetscMalloc1(MAX_KV_LEVELS, &ps->kvlevels);
  CHKERRQ(ierr);

  ierr = PSIncreaseReferenceCount(ps);
  CHKERRQ(ierr);

  ps->gic_file_set = PETSC_FALSE;
  ps->setupcalled = PETSC_FALSE;

  ps->solve_real_time = 0.0;
  ps->solve_cpu_time = 0.0;

  *psout = ps;
  PetscFunctionReturn(0);
}

/*
  PSDestroy - Destroys the PS object created with PSCreate

  Input Parameter
  ps - pointer to PS object to be destroyed
*/
PetscErrorCode PSDestroy(PS *ps) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (!(*ps) || !(*ps)->setupcalled)
    PetscFunctionReturn(0);

  ierr = PSDecreaseReferenceCount(*ps);
  CHKERRQ(ierr);
  if ((*ps)->refct > 0) {
    *ps = 0;
    PetscFunctionReturn(0);
  }

  /* Destroy PS communicator */
  ierr = COMMDestroy(&(*ps)->comm);
  CHKERRQ(ierr);

  /* Destroy PS bus,branch,gen,load objects */
  ierr = PetscFree((*ps)->bus);
  CHKERRQ(ierr);
  ierr = PSGENDestroy((*ps));
  CHKERRQ(ierr);
  ierr = PSLOADDestroy((*ps));
  CHKERRQ(ierr);
  ierr = PetscFree((*ps)->load);
  CHKERRQ(ierr);
  ierr = PetscFree((*ps)->line);
  CHKERRQ(ierr);

  ierr = PSConnCompDestroy(*ps);
  CHKERRQ(ierr);

  if ((*ps)->gic_file_set) {
    ierr = PetscFree((*ps)->substations);
    CHKERRQ(ierr);
  }

  ierr = PetscFree((*ps)->kvlevels);
  CHKERRQ(ierr);
  ierr = PetscFree((*ps)->busext2intmap);
  CHKERRQ(ierr);
  ierr = DMDestroy(&(*ps)->networkdm);
  CHKERRQ(ierr);
  ierr = PetscFree(*ps);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  PSGetNumGlobalLines - Gets the total number of lines in the PS network

  Input Parameters:
+ PS     - the PS network object
- Nlines - the total number of lines (including branches, transformers) in the
network
*/
PetscErrorCode PSGetNumGlobalLines(PS ps, PetscInt *Nlines) {
  PetscFunctionBegin;
  *Nlines = ps->Nline;
  PetscFunctionReturn(0);
}

/*
  PSGetNumGlobalBuses - Gets the total number of buses in the PS network

  Input Parameters:
+ PS     - the PS network object
- Nbuses - the total number of buses in the network
*/
PetscErrorCode PSGetNumGlobalBuses(PS ps, PetscInt *Nbuses) {
  PetscFunctionBegin;
  *Nbuses = ps->Nbus;
  PetscFunctionReturn(0);
}

/*
  PSGetLineConnectivity - Gets the connectivitiy (from bus - to bus connection)
of all lines

  Input Parameters:
+ PS - the PS object
. Nlines - total number of lines
- lineconn  - the line connectivity area (of size 2*Nlines)

  Notes:
  The line connectivity for line i is stored as lines[2*i] = frombus_i,
lines[2*i+1] = tobus_i
*/
PetscErrorCode PSGetLineConnectivity(PS ps, PetscInt Nlines,
                                     PetscInt lineconn[]) {
  PetscInt i, fbus, tbus;
  PSLINE line = ps->line;

  PetscFunctionBegin;
  for (i = 0; i < Nlines; i++) {
    fbus = line[i].internal_i;
    tbus = line[i].internal_j;
    lineconn[2 * i] = fbus;
    lineconn[2 * i + 1] = tbus;
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
PetscErrorCode PSSetUp(PS ps) {
  PetscErrorCode ierr;
  DM networkdm;
  MPI_Comm mpicomm = ps->comm->type;
  PetscInt Nlines = -1, Nbuses = -1;
  PetscInt numbusvariables =
      0; /* Number of variables at each bus..set by the application */
  void *component;
  PetscInt key = -1;
  PetscInt numComponents = 0;
  PetscInt i = 0, j = 0;
  PSGEN gen;
  PSLOAD load;
  PetscInt nv = 0, ne = 0;
  const PetscInt *vtx, *edge;
  PetscInt *lineconn[1];

  PetscFunctionBegin;
  if (ps->setupcalled)
    PetscFunctionReturn(0);

  /* Create empty DMNetwork object */
  ierr = DMNetworkCreate(mpicomm, &networkdm);
  CHKERRQ(ierr);

  /* Register bus,branch,gen,load objects */
  ierr = DMNetworkRegisterComponent(networkdm, "PSLINE",
                                    sizeof(struct _p_PSLINE), &ps->compkey[0]);
  CHKERRQ(ierr);
  ierr = DMNetworkRegisterComponent(networkdm, "PSBUS", sizeof(struct _p_PSBUS),
                                    &ps->compkey[1]);
  CHKERRQ(ierr);
  ierr = DMNetworkRegisterComponent(networkdm, "PSGEN", sizeof(struct _p_PSGEN),
                                    &ps->compkey[2]);
  CHKERRQ(ierr);
  ierr = DMNetworkRegisterComponent(networkdm, "PSLOAD",
                                    sizeof(struct _p_PSLOAD), &ps->compkey[3]);
  CHKERRQ(ierr);

  /* Get the total number of buses and lines */
  /* Note that when the read is read from XXXReadMatPowerData, only P0 reads the
     data and has NumLines and NumBuses set, all the other processors don't have
     any bus, branch, gen, load data set.
  */
  ierr = PSGetNumGlobalLines(ps, &Nlines);
  CHKERRQ(ierr);
  ierr = PSGetNumGlobalBuses(ps, &Nbuses);
  CHKERRQ(ierr);

  /* Set up edge connectivity */

  ierr = PetscCalloc1(2 * Nlines, &lineconn[0]);
  CHKERRQ(ierr);
  ierr = PSGetLineConnectivity(ps, Nlines, lineconn[0]);
  CHKERRQ(ierr);

  /* Set sizes for the network and provide edge connectivity information */
  ierr = DMNetworkSetNumSubNetworks(networkdm, PETSC_DECIDE, 1);
  CHKERRQ(ierr);
  ierr = DMNetworkAddSubnetwork(networkdm, "", Nlines, lineconn[0], NULL);
  CHKERRQ(ierr);

  /* Set up network layout */
  ierr = DMNetworkLayoutSetUp(networkdm);
  CHKERRQ(ierr);
  ierr = PetscFree(lineconn[0]);
  CHKERRQ(ierr);

  /* Add network components and variables */
  /* Associate and copy line data to the edges in networkdm, add number of
   * variables if applicable */
  ierr = DMNetworkGetSubnetwork(networkdm, 0, &nv, &ne, &vtx, &edge);
  CHKERRQ(ierr);
  for (i = 0; i < ne; i++) {
    ierr = DMNetworkAddComponent(networkdm, edge[i], ps->compkey[0],
                                 &ps->line[i], 0);
    CHKERRQ(ierr);
  }
  /* Associate and copy bus, gen, branch to the vertices of the networkdm */
  for (i = 0; i < nv; i++) {
    /* Set the number of variables for buses */
    if (ps->appname == APP_ACPF)
      numbusvariables = 2;
    if (ps->appname == APP_ACOPF)
      numbusvariables = 2; /* The variables are set later by the application */

    ierr = DMNetworkAddComponent(networkdm, vtx[i], ps->compkey[1], &ps->bus[i],
                                 numbusvariables);
    CHKERRQ(ierr);

    for (j = 0; j < ps->bus[i].ngen; j++) {
      /* Add generator */
      gen = &ps->gen[ps->bus[i].gidx[j]];
      ierr = DMNetworkAddComponent(networkdm, vtx[i], ps->compkey[2], gen, 0);
      CHKERRQ(ierr);
    }

    /* Loads */
    for (j = 0; j < ps->bus[i].nload; j++) {
      /* Add load */
      load = &ps->load[ps->bus[i].lidx[j]];
      ierr = DMNetworkAddComponent(networkdm, vtx[i], ps->compkey[3], load, 0);
      CHKERRQ(ierr);
    }
  }

  /* Set up DM for use */
  ierr = DMSetUp(networkdm);
  CHKERRQ(ierr);

  ierr = DMNetworkDistribute(&networkdm, 0);
  CHKERRQ(ierr);
  ps->networkdm = networkdm;

  PetscBool pdm_view = PETSC_FALSE;
  ierr = PetscOptionsHasName(NULL, NULL, "-pdm_view", &pdm_view);
  CHKERRQ(ierr);
  if (pdm_view) {
    ierr = DMView(networkdm, PETSC_VIEWER_STDOUT_WORLD);
    CHKERRQ(ierr);
  }

  ierr = DMNetworkGetSubnetwork(networkdm, 0, &nv, &ne, &vtx, &edge);
  CHKERRQ(ierr);

  /* Set local sizes of buses and branches */
  ps->nbus = nv;
  ps->nline = ne;

  if (ps->comm->size > 1) {
    if (ps->comm->rank == 0) {
      /* Free PSBUS,PSBRANCH,PSGEN, and PSLOAD on root */
      ierr = PetscFree(ps->bus);
      CHKERRQ(ierr);
      ierr = PetscFree(ps->line);
      CHKERRQ(ierr);
      ierr = PSGENDestroy(ps);
      CHKERRQ(ierr);
      ierr = PSLOADDestroy(ps);
      CHKERRQ(ierr);
      ierr = PetscFree(ps->load);
      CHKERRQ(ierr);
      ierr = PetscFree(ps->busext2intmap);
      CHKERRQ(ierr);
    }

    /* Broadcast global Nbus,Ngen,Nbranch, Nload,and maxbusnum to all processors
     */
    PetscInt temp[8];
    /* Pack variables */
    temp[0] = ps->Nbus;
    temp[1] = ps->Ngen;
    temp[2] = ps->Nline;
    temp[3] = ps->Nload;
    temp[4] = ps->maxbusnum;
    temp[5] = ps->NgenON;
    temp[6] = ps->NlineON;
    temp[7] = ps->Nref;
    ierr = MPI_Bcast(temp, 8, MPI_INT, 0, ps->comm->type);
    CHKERRQ(ierr);
    /* Unpack */
    ps->Nbus = temp[0];
    ps->Ngen = temp[1];
    ps->Nline = temp[2];
    ps->Nload = temp[3];
    ps->maxbusnum = temp[4];
    ps->NgenON = temp[5];
    ps->NlineON = temp[6];
    ps->Nref = temp[7];

    /* Recreate busext2intmap..this will map the local bus numbers to external
     * numbers */
    ierr = PetscCalloc1(ps->maxbusnum + 1, &ps->busext2intmap);
    CHKERRQ(ierr);
    for (i = 0; i < ps->maxbusnum + 1; i++)
      ps->busext2intmap[i] = -1;

    ps->ngen = ps->nload = ps->ngenON = ps->nlineON = 0;
    /* Get the local number of gens and loads */
    for (i = 0; i < nv; i++) {
      ierr = DMNetworkGetNumComponents(ps->networkdm, vtx[i], &numComponents);
      CHKERRQ(ierr);
      for (j = 0; j < numComponents; j++) {
        ierr = DMNetworkGetComponent(ps->networkdm, vtx[i], j, &key, &component,
                                     NULL);
        CHKERRQ(ierr);
        if (key == ps->compkey[2])
          ps->ngen++;
        else if (key == ps->compkey[3])
          ps->nload++;
      }
    }

    /* Create local PSBUS, PSBRANCH, PSGEN, and PSLOAD */
    ierr = PetscCalloc1(ps->nbus, &ps->bus);
    CHKERRQ(ierr);
    ierr = PetscCalloc1(ps->nline, &ps->line);
    CHKERRQ(ierr);
    ierr = PetscCalloc1(ps->ngen, &ps->gen);
    CHKERRQ(ierr);
    ierr = PetscCalloc1(ps->nload, &ps->load);
    CHKERRQ(ierr);

    /* Copy line data from DMNetwork data array to PSLINE */
    for (i = 0; i < ne; i++) {
      ierr = DMNetworkGetComponent(ps->networkdm, edge[i], 0, &key, &component,
                                   NULL);
      CHKERRQ(ierr);
      ierr = PetscMemcpy(&ps->line[i], component, sizeof(struct _p_PSLINE));
      CHKERRQ(ierr);
      ps->nlineON += ps->line[i].status;
    }
    PetscInt genj = 0, loadj = 0;
    PetscInt genctr, loadctr;
    /* Copy bus, gen, branch data from DMNetwork data array to PSBUS, PSGEN, and
     * PSBRANCH */
    for (i = 0; i < nv; i++) {
      genctr = loadctr = 0;
      ierr = DMNetworkGetNumComponents(ps->networkdm, vtx[i], &numComponents);
      CHKERRQ(ierr);
      for (j = 0; j < numComponents; j++) {
        ierr = DMNetworkGetComponent(ps->networkdm, vtx[i], j, &key, &component,
                                     NULL);
        CHKERRQ(ierr);
        if (key == ps->compkey[1]) {
          ierr = PetscMemcpy(&ps->bus[i], component, sizeof(struct _p_PSBUS));
          CHKERRQ(ierr);
          /* Set external to internal mapping */
          ps->busext2intmap[ps->bus[i].bus_i] = ps->bus[i].internal_i = i;
        } else if (key == ps->compkey[2]) {
          ierr =
              PetscMemcpy(&ps->gen[genj], component, sizeof(struct _p_PSGEN));
          CHKERRQ(ierr);
          ps->ngenON += ps->gen[genj].status;
          ps->bus[i].gidx[genctr++] = genj++;
        } else if (key == ps->compkey[3]) {
          ierr = PetscMemcpy(&ps->load[loadj], component,
                             sizeof(struct _p_PSLOAD));
          CHKERRQ(ierr);
          ps->bus[i].lidx[loadctr++] = loadj++;
        }
      }
    }
  } else {
    ps->ngen = ps->Ngen;
    ps->nload = ps->Nload;
    ps->ngenON = ps->NgenON;
    ps->nlineON = ps->NlineON;
  }

  /* Set up
     (a) connectivity information for lines and buses
     (b) bus ghosted status
     (c) incident generators at bus
     (d) incident loads at bus
     (e) kv levels for lines
     (e) sets the starting location for the variables for this bus in the given
     application state vector
  */
  PetscInt eStart, eEnd, vStart, vEnd;
  PetscInt nlines, k;
  const PetscInt *connnodes, *connlines;
  ierr = DMNetworkGetVertexRange(ps->networkdm, &vStart, &vEnd);
  CHKERRQ(ierr);
  ierr = DMNetworkGetEdgeRange(ps->networkdm, &eStart, &eEnd);
  CHKERRQ(ierr);
  for (i = 0; i < ne; i++) {
    ierr = DMNetworkGetConnectedVertices(ps->networkdm, edge[i], &connnodes);
    CHKERRQ(ierr);
    ps->line[i].connbuses[0] = &ps->bus[connnodes[0] - vStart];
    ps->line[i].connbuses[1] = &ps->bus[connnodes[1] - vStart];
    ps->line[i].internal_i = ps->bus[connnodes[0] - vStart].internal_i;
    ps->line[i].internal_j = ps->bus[connnodes[1] - vStart].internal_i;
    ps->line[i].kvlevel = PetscMax(ps->bus[connnodes[0] - vStart].basekV,
                                   ps->bus[connnodes[1] - vStart].basekV);

    /* Starting location in the local vector */
    ierr = DMNetworkGetLocalVecOffset(ps->networkdm, edge[i], ALL_COMPONENTS,
                                      &ps->line[i].startloc);
    CHKERRQ(ierr);
    /* Starting location in the global vector */
    ierr = DMNetworkGetGlobalVecOffset(ps->networkdm, edge[i], ALL_COMPONENTS,
                                       &ps->line[i].startlocglob);
    CHKERRQ(ierr);
  }

  for (i = 0; i < nv; i++) {
    /* Connected lines */
    ierr =
        DMNetworkGetSupportingEdges(ps->networkdm, vtx[i], &nlines, &connlines);
    CHKERRQ(ierr);
    if (nlines > MAXCONNLINES)
      SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
              "%d lines connected to a bus %d exceeds max. allowed connected "
              "lines allowed %d",
              nlines, ps->bus[i].bus_i, MAXCONNLINES);
    ps->bus[i].nconnlines = nlines;
    for (k = 0; k < nlines; k++)
      ps->bus[i].connlines[k] = &ps->line[connlines[k] - eStart];

    /* Is bus ghosted? */
    ierr = DMNetworkIsGhostVertex(ps->networkdm, vtx[i], &ps->bus[i].isghost);
    CHKERRQ(ierr);

    /* Starting location in the local array */
    ierr = DMNetworkGetLocalVecOffset(ps->networkdm, vtx[i], ALL_COMPONENTS,
                                      &ps->bus[i].startloc);
    CHKERRQ(ierr);

    /* Starting location in the global array */
    ierr = DMNetworkGetGlobalVecOffset(ps->networkdm, vtx[i], ALL_COMPONENTS,
                                       &ps->bus[i].startlocglob);
    CHKERRQ(ierr);

    /* Incident generators */

    for (k = 0; k < ps->bus[i].ngen; k++) {
      ps->bus[i].gens[k] = &ps->gen[ps->bus[i].gidx[k]];
    }
    /* Change the bus type to PQ if no generators are incident */
    if (!ps->bus[i].ngenON && ps->bus[i].ide != ISOLATED_BUS)
      ps->bus[i].ide = PQ_BUS;

    /* Incident loads */
    for (k = 0; k < ps->bus[i].nload; k++) {
      ps->bus[i].loads[k] = &ps->load[ps->bus[i].lidx[k]];
    }

    /* Update the number of local reference buses */
    if (ps->bus[i].ide == REF_BUS)
      ps->nref++;
  }
#if defined DEBUGPS
  ierr = PetscPrintf(PETSC_COMM_SELF,
                     "Rank[%d]:nbuses = %d,nlines = %d,ngen = %d, nload = %d\n",
                     ps->comm->rank, ps->nbus, ps->nline, ps->ngen, ps->nload);
  CHKERRQ(ierr);
  PSLINE line;
  for (i = 0; i < ps->nline; i++) {
    line = &ps->line[i];
    ierr = PetscPrintf(PETSC_COMM_SELF, "Rank[%d]:Line %d ----- %d\n",
                       ps->comm->rank, line->fbus, line->tbus);
    CHKERRQ(ierr);
  }
#endif

  /* Broadcast MVAbase */
  ierr = MPI_Bcast(&ps->MVAbase, 1, MPIU_SCALAR, 0, ps->comm->type);
  CHKERRQ(ierr);

  /* Update reference bus if needed */
  ierr = PSCheckandSetRefBus(ps);
  CHKERRQ(ierr);

  /* Get KV levels */
  ierr = PSGetKVLevels(ps, &ps->nkvlevels, (const PetscScalar **)&ps->kvlevels);
  CHKERRQ(ierr);
  //  ierr = PetscPrintf(PETSC_COMM_SELF,"Rank %d Came
  //  here\n",ps->comm->rank);CHKERRQ(ierr);
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
PetscErrorCode PSCreateGlobalVector(PS ps, Vec *vec) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (!ps->setupcalled)
    SETERRQ(ps->comm->type, 0,
            "PSSetUp() must be called before calling PSCreateGlobalVector");
  ierr = DMCreateGlobalVector(ps->networkdm, vec);
  CHKERRQ(ierr);
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
PetscErrorCode PSCreateMatrix(PS ps, Mat *mat) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (!ps->setupcalled)
    SETERRQ(ps->comm->type, 0,
            "PSSetUp() must be called before calling PSCreateMatrix");
  ierr = DMCreateMatrix(ps->networkdm, mat);
  CHKERRQ(ierr);
  ierr = MatSetOption(*mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  PSSetGenStatus - Sets generator status given the bus name and id

  Input Parameters
+ ps     - the PS objec
. gbus   - generator bus
. gid    - generator id
- status - Generator status (1 = ON, 0 = OFF)

  Notes: PSSetUp() must be called before calling PSSetGenStatus()
*/
PetscErrorCode PSSetGenStatus(PS ps, PetscInt gbus, const char *gid,
                              PetscInt status) {
  PetscErrorCode ierr;
  PSGEN gen;
  PSBUS bus;
  PetscBool statusupdate = PETSC_FALSE;

  PetscFunctionBegin;
  if (!ps->setupcalled)
    SETERRQ(PETSC_COMM_SELF, 0,
            "PSSetUp() must be called before calling PSSetGenStatus()\n");

  ierr = PSGetGen(ps, gbus, gid, &gen);
  CHKERRQ(ierr);

  if (gen) {
    if (status != gen->status)
      statusupdate = PETSC_TRUE;
    ierr = PSGENSetStatus(gen, status);
    CHKERRQ(ierr);
    if (statusupdate) {
      bus = &ps->bus[ps->busext2intmap[gbus]];
      if (status) {
        bus->ngenON++;
        ps->ngenON++;
        ps->NgenON++;
      } else {
        bus->ngenON--;
        ps->ngenON--;
        ps->NgenON--;
        /* Change the bus type to PQ if no generators are active at this bus */
        if (!bus->ngenON) {
          if (bus->ide == REF_BUS)
            ps->nref--;
          bus->ide = PQ_BUS;
        }
        gen->pg = gen->qg = 0.0;
      }
    }
  } else {
    SETERRQ(PETSC_COMM_SELF, 0, "No generator found on bus %d with %s id\n",
            gbus, gid);
  }

  /* Update reference bus if needed */
  ierr = PSCheckandSetRefBus(ps);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  PSSetGenStatus - Gets generator status given the bus name and id

  Input Parameters
+ ps     - the PS object
. gbus   - generator bus
. gid    - generator id
- status - Generator status

  Notes: PSSetUp() must be called before calling PSSetGenStatus()
*/
PetscErrorCode PSGetGenStatus(PS ps, PetscInt gbus, const char *gid,
                              PetscInt *status) {
  PetscErrorCode ierr;
  PSGEN gen;

  PetscFunctionBegin;
  if (!ps->setupcalled)
    SETERRQ(PETSC_COMM_SELF, 0,
            "PSSetUp() must be called before calling PSSetGenStatus()\n");

  ierr = PSGetGen(ps, gbus, gid, &gen);
  CHKERRQ(ierr);

  if (gen) {
    if (status)
      *status = gen->status;
  } else {
    SETERRQ(PETSC_COMM_SELF, 0, "No generator found on bus %d with %s id\n",
            gbus, gid);
  }

  PetscFunctionReturn(0);
}

/*
  PSSetLineStatus - Sets the line status given from bus, to bus, and line id
*/
PetscErrorCode PSSetLineStatus(PS ps, PetscInt fbus, PetscInt tbus,
                               const char *id, PetscInt status) {
  PetscErrorCode ierr;
  PSLINE line = NULL;
  PetscBool statusupdate = PETSC_FALSE;

  PetscFunctionBegin;
  if (!ps->setupcalled)
    SETERRQ(PETSC_COMM_SELF, 0,
            "PSSetUp() must be called before calling PFLOWGetLine()\n");

  ierr = PSGetLine(ps, fbus, tbus, id, &line);
  CHKERRQ(ierr);

  if (line) {
    if (line->status != status)
      statusupdate = PETSC_TRUE;
    ierr = PSLINESetStatus(line, status);
    CHKERRQ(ierr);
    if (statusupdate) {
      if (!status)
        ps->nlineON--; /* Line switching to OFF status */
      else
        ps->nlineON++; /* Line switching to ON status */
    }
  }
  ierr = MPI_Allreduce(&ps->nlineON, &ps->NlineON, 1, MPIU_INT, MPI_SUM,
                       ps->comm->type);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode PSSetEdgeandBusStartLoc(PS ps) {
  PetscErrorCode ierr;
  PetscInt i, vStart, vEnd, eStart, eEnd;

  PetscFunctionBegin;
  /* Reset the start locations for the buses and lines */
  ierr = DMNetworkGetVertexRange(ps->networkdm, &vStart, &vEnd);
  CHKERRQ(ierr);
  ierr = DMNetworkGetEdgeRange(ps->networkdm, &eStart, &eEnd);
  CHKERRQ(ierr);

  for (i = eStart; i < eEnd; i++) {
    ierr = DMNetworkGetLocalVecOffset(ps->networkdm, i, ALL_COMPONENTS,
                                      &ps->line[i].startloc);
    CHKERRQ(ierr);
    ierr = DMNetworkGetGlobalVecOffset(ps->networkdm, i, ALL_COMPONENTS,
                                       &ps->line[i].startlocglob);
    CHKERRQ(ierr);
  }

  for (i = vStart; i < vEnd; i++) {
    ierr = DMNetworkGetLocalVecOffset(ps->networkdm, i, ALL_COMPONENTS,
                                      &ps->bus[i - vStart].startloc);
    CHKERRQ(ierr);
    ierr = DMNetworkGetGlobalVecOffset(ps->networkdm, i, ALL_COMPONENTS,
                                       &ps->bus[i - vStart].startlocglob);
    CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode PSApplyContingency(PS ps, Contingency ctgc) {
  PetscInt j;
  PetscErrorCode ierr;
  PetscFunctionBegin;

  for (j = 0; j < ctgc.noutages; j++) {
    if (ctgc.outagelist[j].type == GEN_OUTAGE) {
      PetscInt gbus = ctgc.outagelist[j].bus;
      char *gid = ctgc.outagelist[j].id;
      PetscInt status = ctgc.outagelist[j].status;
      ierr = PSSetGenStatus(ps, gbus, gid, status);
      CHKERRQ(ierr);
    }
    if (ctgc.outagelist[j].type == BR_OUTAGE) {
      PetscInt fbus = ctgc.outagelist[j].fbus;
      PetscInt tbus = ctgc.outagelist[j].tbus;
      char *brid = ctgc.outagelist[j].id;
      PetscInt status = ctgc.outagelist[j].status;
      ierr = PSSetLineStatus(ps, fbus, tbus, brid, status);
      CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode PSApplyScenario(PS ps, Scenario scenario) {
  PetscInt i, j;
  Forecast *forecast;
  PSGEN gen;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  for (i = 0; i < scenario.nforecast; i++) {
    forecast = &scenario.forecastlist[i];
    if (forecast->type == FORECAST_WIND) {
      /* Wind forecast */
      for (j = 0; j < forecast->nele; j++) {
        ierr = PSGetGen(ps, forecast->buses[j], forecast->id[j], &gen);
        CHKERRQ(ierr);
        gen->pg = gen->pt =
            forecast->val[j] /
            ps->MVAbase;    /* Set real power generation. Note that
                 Pg upper limit is also set to the forecast value. This allows the
                 wind generator to be dispatched at its limit */
        gen->pgs = gen->pb; /* Set-point set to lower bound. This is good for
                               optimization */
        if (PetscAbsScalar(gen->pg) < 1e-6)
          gen->status =
              0; /* Generation value is zero, so switch off the generator */
      }
    }
  }

  PetscFunctionReturn(0);
}

/*
   PSGetTotalGeneration - Gets the total online generation and capacity

  Input Parameters:
. ps - power system object

  Output Parameters:
+ pgentot - total online generation in pu on system MVAbase
- pmaxtot - total online capacity in pu on system MVAbase
*/
PetscErrorCode PSGetTotalGeneration(PS ps, PetscScalar *pgtot,
                                    PetscScalar *pmaxtot) {
  PetscErrorCode ierr;
  PSBUS bus;
  PSGEN gen;
  PetscInt i, j;
  PetscScalar pgtot_loc = 0., pmaxtot_loc = 0.;
  PetscBool isghosted;
  PetscScalar in[2], out[2];

  PetscFunctionBegin;

  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSIsGhosted(bus, &isghosted);
    if (isghosted)
      continue;
    for (j = 0; j < bus->ngen; j++) {
      ierr = PSBUSGetGen(bus, j, &gen);
      if (gen->status) {
        pgtot_loc += gen->pg;
        pmaxtot_loc += gen->pt;
      }
    }
  }

  in[0] = pgtot_loc;
  in[1] = pmaxtot_loc;
  out[0] = out[1] = 0.0;
  ierr = MPI_Allreduce(in, out, 2, MPIU_SCALAR, MPI_SUM, ps->comm->type);
  CHKERRQ(ierr);

  if (pgtot)
    *pgtot = out[0];
  if (pmaxtot)
    *pmaxtot = out[1];

  PetscFunctionReturn(0);
}

/*
  PSComputeParticipationFactors - Computes generation participation factor

  Input Parameters:
+ ps - power system object
- pmaxtot - total generation capacity

  Notes:
  - Currenty, this function only computes the participation factors for the
entire system, not by areas
  - The participation factor for gen. j is p^{max}_j/pmaxtot
  - pmaxtot can be obtained by PSGetTotalGeneration
*/
PetscErrorCode PSComputeParticipationFactors(PS ps) {
  PetscErrorCode ierr;
  PSBUS bus;
  PSGEN gen;
  PetscInt i, j;
  PetscScalar pmaxtot;

  PetscFunctionBegin;

  ierr = PSGetTotalGeneration(ps, NULL, &pmaxtot);
  CHKERRQ(ierr);

  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    for (j = 0; j < bus->ngen; j++) {
      ierr = PSBUSGetGen(bus, j, &gen);
      CHKERRQ(ierr);
      if (gen->status)
        gen->apf = gen->pt / pmaxtot;
    }
  }

  PetscFunctionReturn(0);
}

/*
  PSSetGenPowerLimits - Sets generator real and reactive power limits given a
bus number and generator id

  Input Parameters
+ ps   - the PS object
. gbus - bus number for the generaor
. id   - generator id
. pt   - Max. generator real power capacity (MW)
. pb   - Min. generator real power capacity (MW)
. qt   - Max. generator reactive power capacity (MVAr)
- qb   - Min. generator reactive power capacity (MVAr)

  Notes: Use EXAGO_IGNORE to not set the value
*/
PetscErrorCode PSSetGenPowerLimits(PS ps, PetscInt gbus, const char *gid,
                                   PetscScalar pt, PetscScalar pb,
                                   PetscScalar qt, PetscScalar qb) {
  PetscErrorCode ierr;
  PSGEN gen;

  PetscFunctionBegin;
  if (!ps->setupcalled)
    SETERRQ(PETSC_COMM_SELF, 0,
            "PSSetUp() or OPFLOWSetUpPS() must be called before calling "
            "PSSetGenStatus()\n");

  ierr = PSGetGen(ps, gbus, gid, &gen);
  CHKERRQ(ierr);

  if (gen) {
    if (pt != EXAGO_IGNORE)
      gen->pt = pt / ps->MVAbase;
    if (pb != EXAGO_IGNORE)
      gen->pb = pb / ps->MVAbase;
    if (qt != EXAGO_IGNORE)
      gen->qt = qt / ps->MVAbase;
    if (qb != EXAGO_IGNORE)
      gen->qb = qb / ps->MVAbase;
  } else {
    SETERRQ(PETSC_COMM_SELF, 0, "No generator found on bus %d with %s id\n",
            gbus, gid);
  }

  PetscFunctionReturn(0);
}

/*
  PSGetGenPowerLimits - Returns generator real and reactive power limits given a
bus number and generator id

  Input Parameters
+ ps   - the PS object
. gbus - bus number for the generaor
. gid  - generator id
. pt   - Max. generator real power capacity
. pb   - Min. generator real power capacity
. qt   - Max. generator reactive power capacity
- qb   - Min. generator reactive power capacity

  Notes: Use NULL to ignore retrieving value
*/
PetscErrorCode PSSetGenPowerLimits(PS ps, PetscInt gbus, const char *gid,
                                   PetscScalar *pt, PetscScalar *pb,
                                   PetscScalar *qt, PetscScalar *qb) {
  PetscErrorCode ierr;
  PSGEN gen;

  PetscFunctionBegin;
  if (!ps->setupcalled)
    SETERRQ(PETSC_COMM_SELF, 0,
            "PSSetUp() or OPFLOWSetUpPS() must be called before calling "
            "PSSetGenStatus()\n");

  ierr = PSGetGen(ps, gbus, gid, &gen);
  CHKERRQ(ierr);

  if (gen) {
    if (pt)
      *pt = gen->pt;
    if (pb)
      *pb = gen->pb;
    if (qt)
      *qt = gen->qt;
    if (qb)
      *qb = gen->qb;
  } else {
    SETERRQ(PETSC_COMM_SELF, 0, "No generator found on bus %d with %s id\n",
            gbus, gid);
  }

  PetscFunctionReturn(0);
}

/*
  PSGetGenDispatch - Returns generator real (MW) and reactive power (MVAR)
dispatch given a bus number and generator id

  Input Parameters
+ ps   - the PS object
. gbus - bus number for the generator
- id   - generator id

  Output Paramters
+ pg   - Real power dispatch (MW)
- qg   - Reactive power dispatch (MVAr)

  Notes: This function should be called after OPFLOWSolutionToPS()
*/
PetscErrorCode PSGetGenDispatch(PS ps, PetscInt gbus, const char *gid,
                                PetscScalar *pg, PetscScalar *qg) {
  PetscErrorCode ierr;
  PSGEN gen;

  PetscFunctionBegin;

  ierr = PSGetGen(ps, gbus, gid, &gen);
  CHKERRQ(ierr);

  if (gen) {
    if (pg)
      *pg = gen->pg * ps->MVAbase;
    if (qg)
      *qg = gen->qg * ps->MVAbase;
  } else {
    SETERRQ(PETSC_COMM_SELF, 0, "No generator found on bus %d with %s id\n",
            gbus, gid);
  }

  PetscFunctionReturn(0);
}

/*
  PSCopy - Copies PS data

  Inputs:
+ psin  - PS object from which data is copied
- psout - PS object where data is copied

  Notes: This routine copies data for the bus, branch, gen, and load
         structs.

         Must be called after PSSetUp() is called.

*/
PetscErrorCode PSCopy(PS psin, PS psout) {
  PetscErrorCode ierr;
  PetscInt i, k;
  PSBUS busin, busout;
  PSGEN genin, genout;
  PSLINE linein, lineout;
  PSLOAD loadin, loadout;

  PetscFunctionBegin;

  /* Copy bus data */
  for (i = 0; i < psin->nbus; i++) {
    busin = &psin->bus[i];
    busout = &psout->bus[i];

    /* voltage magnitude and angle */
    busout->vm = busin->vm;
    busout->va = busin->va;

    /* mismatch lagrange multipliers */
    busout->mult_pmis = busin->mult_pmis;
    busout->mult_qmis = busin->mult_qmis;

    /* Power imbalance */
    busout->pimb = busin->pimb;
    busout->qimb = busin->qimb;

    /* Copy generator data */
    for (k = 0; k < busin->ngen; k++) {
      ierr = PSBUSGetGen(busin, k, &genin);
      CHKERRQ(ierr);
      ierr = PSBUSGetGen(busout, k, &genout);
      CHKERRQ(ierr);

      genout->pg = genin->pg;
      genout->qg = genin->qg;
    }

    /* Copy load data */
    for (k = 0; k < busin->nload; k++) {
      ierr = PSBUSGetLoad(busin, k, &loadin);
      CHKERRQ(ierr);
      ierr = PSBUSGetLoad(busout, k, &loadout);
      CHKERRQ(ierr);
      loadout->pl = loadin->pl;
      loadout->ql = loadin->ql;
      loadout->pl_loss = loadin->pl_loss;
      loadout->ql_loss = loadin->ql_loss;
    }
  }

  /* Copy branch data */
  for (i = 0; i < psin->nline; i++) {
    linein = &psin->line[i];
    lineout = &psout->line[i];

    lineout->pf = linein->pf;
    lineout->qf = linein->qf;
    lineout->pt = linein->pt;
    lineout->qt = linein->qt;
    lineout->sf = linein->sf;
    lineout->st = linein->st;

    lineout->mult_sf = linein->mult_sf;
    lineout->mult_st = linein->mult_st;
  }

  PetscFunctionReturn(0);
}

/*
  PSSetGICData - Sets the GIC data file

  Input Parameter
+  ps - The PS object
-  gicfile - The name of the GIC data file

  Notes: The GIC data file is only used for visualization. It contains the
substation geospatial coordinates and mapping of buses to substations. See the
Electric Grid Data Repository files for examples of GIC data files (given for
CTIVSg cases)
*/
PetscErrorCode PSSetGICData(PS ps, const char gicfile[]) {
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = PetscStrcpy(ps->gic_file_name, gicfile);
  CHKERRQ(ierr);

  ps->gic_file_set = PETSC_TRUE;
  PetscFunctionReturn(0);
}
