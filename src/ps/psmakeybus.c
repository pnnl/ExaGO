#include <private/psimpl.h>

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

  ierr = MatView(Yb,0);
  exit(1);
  *Ybus = Yb;
  PetscFunctionReturn(0);
}
