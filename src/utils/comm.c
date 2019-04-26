#include <common.h>

/*
  COMMCreate - Creates the communicator object COMM

  Input Parameter
. mpicomm - The MPI_Comm

  Output Parameter
. outcomm - The COMM object
*/
PetscErrorCode COMMCreate(MPI_Comm mpicomm,COMM *commout)
{
  PetscErrorCode ierr;
  COMM           comm;
  
  PetscFunctionBegin;
  ierr = PetscCalloc1(1,&comm);CHKERRQ(ierr);
  /* Set up the communicator information for later use (if need be) */
  comm->type = mpicomm;
  ierr = MPI_Comm_rank(mpicomm,&comm->rank);CHKERRQ(ierr);
  ierr = MPI_Comm_size(mpicomm,&comm->size);CHKERRQ(ierr);

  *commout = comm;
  PetscFunctionReturn(0);
}

/*
  COMMDestroy - Destroys the communicator object COMM created with COMMCreate

  Input Parameter
. comm - The COMM object
*/
PetscErrorCode COMMDestroy(COMM *comm)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if(!(*comm)) PetscFunctionReturn(0);
  ierr = PetscFree((*comm));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
