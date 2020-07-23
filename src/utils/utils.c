#include <common.h>
#include <utils.h>
#include <sys/utsname.h>
#include <sys/stat.h>
#include <dirent.h>

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

/**
 * @brief Checks whether or not a directory entry can be opened at the given
 * path using only POSIX standard C system calls
 *
 * @param[in] pth filepath to stat
 * @see stat
 **/
int doesFileExist(char* pth)
{
  struct stat path_stat;
  stat(pth, &path_stat);
  return S_ISREG(path_stat.st_mode);
}

/**
 * @brief Checks whether or not a directory entry can be opened at the given
 * path using only POSIX standard C system calls.
 *
 * @param[in] path path to verify is statable
 * @return 0 if pth==NULL or pth cannot be opened as directory with syscall
 *         1 else
 *
 * @see opendir
 **/
int doesDirExist(char* pth)
{
  if (pth==NULL) return 0;
  DIR *dp;
  dp = opendir(pth);
  if (dp == NULL) return 0;
  return 1;
}

/**
 * @brief Verifies that all paths passed are statable as regular files
 *
 * @param[in] pths  array of paths to check exist
 * @param[in] npths number of paths to check
 * @return first i in [0, npths) such that pths[i] is statable as a regular file
 *         -2 if pths==NULL or pths[i]==NULL for any i in [0, npths),
 *         -1 if !doesFileExist(pths[i]) for i in [0, npths)
 * 
 * @see doesFileExist
 **/
int doFilesExist(char** pths, int npths)
{
  if (pths==NULL)
    return -2;
  for (int i=0; i<npths; i++)
  {
    if (pths[i]==NULL)
      return -2;
    printf("-- Checking %-70s exists: ", pths[i]);
    if (doesFileExist(pths[i]))
    {
      puts("yes");
      return i;
    }
    else
    {
      puts("no");
    }
  }
  return -1;
}
