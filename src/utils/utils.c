#include <common.h>

/*
  SetMatrixValues - Sets the values in the matrix

  Input Parameters:
+ J    - the matrix
. nrow - number of rows
. row  - row locations
. ncol - number of columns
. col  - column locations
- val  - values to set

  Notes:
    Currently, only ADD_VALUES insert mode is used.
*/
PetscErrorCode SetMatrixValues(Mat J, PetscInt nrow, PetscInt row[],PetscInt ncol,PetscInt col[],PetscScalar val[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = MatSetValues(J,nrow,row,ncol,col,val,ADD_VALUES);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
