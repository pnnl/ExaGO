#include <private/psimpl.h>

/*
  PSMakeYbusred - Creates the reduced Ybus matrix by eliminating the rows/columns corresponding to buses having no injectors (gen, load, branches)

  Input parameters
. ps - The PS application object

  Output parameters
. Ybusred - the reduced Ybus matrix
*/
PetscErrorCode PSMakeYbusred(PS ps,Mat *Ybusred)
{
  PetscInt       i,ctr_act=0,ctr_pas=0;
  PetscInt       nact=0,npas=0; /* Number of active and passive buses */
  PetscInt       *idx_act,*idx_pas; /* Ybus row numbers for active and passive bus contributions */
  PetscErrorCode ierr;
  PSBUS          bus;
  PetscBool      isghost;
  PetscInt       locxg;
  IS             is_act,is_pas;

  PetscFunctionBegin;

  /* Count the number of active buses */
  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSIsGhosted(bus,&isghost);CHKERRQ(ierr);
    if(isghost) continue;
    if(bus->hasinj) nact++;
  }
      
  /* Number of passive buses */
  npas = ps->nbusowned - nact;
  
  /* Create arrays for storing row numbers for active and passive buses in Ybus */
  ierr = PetscCalloc1(2*nact,&idx_act);CHKERRQ(ierr);
  ierr = PetscCalloc1(2*npas,&idx_pas);CHKERRQ(ierr);

  /* Set active and passive bus row numbers */
  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    ierr = PSBUSIsGhosted(bus,&isghost);CHKERRQ(ierr);
    if(isghost) continue;

    ierr = PSBUSGetVariableGlobalLocation(bus,&locxg);CHKERRQ(ierr);
    if(bus->hasinj) {
      idx_act[ctr_act]   = locxg;
      idx_act[ctr_act+1] = locxg+1;
      ctr_act += 2;
    } else {
      idx_pas[ctr_pas]   = locxg;
      idx_pas[ctr_pas+1] = locxg+1;
      ctr_pas += 2;
    }
  }

  /* Create index sets for active and passive bus row numbers */
  ierr = ISCreateGeneral(ps->comm->type,2*nact,idx_act,PETSC_OWN_POINTER,&is_act);CHKERRQ(ierr);
  ierr = ISCreateGeneral(ps->comm->type,2*npas,idx_pas,PETSC_OWN_POINTER,&is_pas);CHKERRQ(ierr);
  
  ierr = ISView(is_act,0);CHKERRQ(ierr);
  ierr = ISView(is_pas,0);CHKERRQ(ierr);
  exit(-1);
  PetscFunctionReturn(0);
}

  

/*
  PSMakeYbus - Creates the Ybus matrix

  Input parameters
. ps - The PS application object

  Output parameters
. Ybus - the Ybus matrix

  Note: The Ybus matrix is a 2*Nbus X 2*Nbus matrix describing network connectivity. While, in most applications Ybus is a matrix with complex data type, here Ybus is a matrix of real data type (hence the expansion from Nbus to 2*Nbus). Each non-zero element of Ybus is a 2X2 matrix of the form [G -B; B G]

  PSSetUp() must be called before calling PSMakeYbus
*/
PetscErrorCode PSMakeYbus(PS ps,Mat *Ybus)
{
  PetscErrorCode ierr;
  PetscInt i;
  Mat      Yb;
  PSBUS    busf,bust;
  PetscInt locxfg,locxtg;
  PetscScalar val[4];
  PetscInt    row[2],col[2];
  PetscScalar    Gff,Bff,Gft,Bft,Gtf,Btf,Gtt,Btt;
  PSLINE   line;
  const PSBUS *connbuses;

  PetscFunctionBegin;
  if(!ps->setupcalled) SETERRQ(PETSC_COMM_SELF,0,"PSSetUp() must be called before calling PSMakeYbus()\n");
  ierr = DMCreateMatrix(ps->networkdm,&Yb);CHKERRQ(ierr);

  for(i=0; i < ps->nline; i++) {
    line = &ps->line[i];
    if(!line->status) continue;
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

    ierr = PSBUSGetVariableGlobalLocation(busf,&locxfg);CHKERRQ(ierr);
    ierr = PSBUSGetVariableGlobalLocation(bust,&locxtg);CHKERRQ(ierr);
    /* Add values at from bus diagonal */
    row[0] = locxfg; row[1] = locxfg+1;
    col[0] = row[0]; col[1] = row[1];
    val[0] = Gff; val[1] = -Bff; val[2] = Bff; val[3] = Gff;
    ierr = MatSetValues(Yb,2,row,2,col,val,ADD_VALUES);CHKERRQ(ierr);

    /* Add values at to bus diagonal */
    row[0] = locxtg; row[1] = locxtg+1;
    col[0] = row[0]; col[1] = row[1];
    val[0] = Gtt; val[1] = -Btt; val[2] = Btt; val[3] = Gtt;
    ierr = MatSetValues(Yb,2,row,2,col,val,ADD_VALUES);CHKERRQ(ierr);

    /* Add values at from-to location */
    row[0] = locxfg; row[1] = locxfg+1;
    col[0] = locxtg; col[1] = locxtg+1;
    val[0] = Gft; val[1] = -Bft; val[2] = Bft; val[3] = Gft;
    ierr = MatSetValues(Yb,2,row,2,col,val,ADD_VALUES);CHKERRQ(ierr);

    /* Add values at to-from location */
    row[0] = locxtg; row[1] = locxtg+1;
    col[0] = locxfg; col[1] = locxfg+1;
    val[0] = Gtf; val[1] = -Btf; val[2] = Btf; val[3] = Gtf;
    ierr = MatSetValues(Yb,2,row,2,col,val,ADD_VALUES);CHKERRQ(ierr);
  }

  ierr = MatAssemblyBegin(Yb,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Yb,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  //  ierr = MatView(Yb,0);
  //  exit(1);
  *Ybus = Yb;

  Mat Ybusred;
  ierr = PSMakeYbusred(ps,&Ybusred);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
