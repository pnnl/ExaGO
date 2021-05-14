static char help[] = "User example calling OPFLOW.\n\n";

#include <opflow.h>
#include <exago_config.h>
#include <utils.hpp>

int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  OPFLOW         opflow;
  OutputFormat   fmt=MATPOWER;
  char           file[PETSC_MAX_PATH_LEN];
  PetscBool      flg=PETSC_FALSE,print_output=PETSC_FALSE,save_output=PETSC_FALSE;
  PetscLogStage  stages[3];
  char           appname[]="opflow";
  MPI_Comm       comm=MPI_COMM_WORLD;
  PetscBool      selfcheck=PETSC_FALSE;

  /** Use `ExaGOLogSetLoggingFileName("opflow-logfile");` to log the output. */
  ierr = ExaGOInitialize(comm,&argc,&argv,appname,help);
  if (ierr)
  {
    fprintf(stderr,"Could not initialize ExaGO application %s.\n",appname);
    return ierr;
  }

  ierr = PetscOptionsGetBool(NULL,NULL,"-print_output",&print_output,NULL);ExaGOCheckError(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-save_output",&save_output,NULL);ExaGOCheckError(ierr);
  
  /* Register stages for profiling application code sections */
  ierr = PetscLogStageRegister("Reading Data",&stages[0]);ExaGOCheckError(ierr);
  ierr = PetscLogStageRegister("Set up",&stages[1]);ExaGOCheckError(ierr);
  ierr = PetscLogStageRegister("Solve",&stages[2]);ExaGOCheckError(ierr);

  /* Stage 1 - Application creation and reading data */
  ierr = PetscLogStagePush(stages[0]);ExaGOCheckError(ierr);

  ExaGOLog(EXAGO_LOG_INFO,"%s","Creating OPFlow\n");

  /* Create OPFLOW object */
  ierr = OPFLOWCreate(comm,&opflow);ExaGOCheckError(ierr);
  
  /* Get network data file from command line */
  ierr = PetscOptionsGetString(NULL,NULL,"-netfile",file,PETSC_MAX_PATH_LEN,&flg);ExaGOCheckError(ierr);
  if(flg) {
    ierr = OPFLOWReadMatPowerData(opflow,file);ExaGOCheckError(ierr);
  } else {
    ierr = OPFLOWReadMatPowerData(opflow,"datafiles/case9/case9mod.m");ExaGOCheckError(ierr);
  }

  ierr = PetscLogStagePop();ExaGOCheckError(ierr);

  /* Stage 2 - Set up OPFLOW application */
  ierr = PetscLogStagePush(stages[1]);ExaGOCheckError(ierr);

  /* Set up */
  ierr = OPFLOWSetUp(opflow);ExaGOCheckError(ierr);

  ierr = PetscLogStagePop();ExaGOCheckError(ierr);

  /* Stage 3 - Solve */
  ierr = PetscLogStagePush(stages[2]);ExaGOCheckError(ierr);

  /* Solve */
  ierr = OPFLOWSolve(opflow);ExaGOCheckError(ierr);

  ierr = PetscLogStagePop();ExaGOCheckError(ierr);
  
  if(print_output) {
    ierr = OPFLOWPrintSolution(opflow);ExaGOCheckError(ierr);
  }

  if(save_output) {
    ierr = OPFLOWSaveSolution(opflow,fmt,"opflowout");ExaGOCheckError(ierr);
  }

  /* Destroy OPFLOW object */
  ierr = OPFLOWDestroy(&opflow);ExaGOCheckError(ierr);

  ExaGOFinalize();
  return 0;
}
