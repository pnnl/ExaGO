static char help[] = "OPFLOW Functionality Test to test manipulation of PS object\n\
                      This example implements the following flow-chart\n\
                      1. Create OPFLOW object and read the network data\n\
                      2. Access the underlying PS object and change the wind generator dispatch\n\
                      3. Set up and solve OPFLOW, and get the objective function\n\
                      4. Retrieve the non-renewable generation dispatch\n\
                      5. Create another OPFLOW object and read the network data again\n\
                      6. Access the underlying PS object\n\
                      7. Set the non-renewable generator min. and max. capacities to those obtained in step 5\n\
                      8. Set the new wind generator dispatch\n\
                      9. Set up and solve OPFLOW\n\
                     10. Get the objective function value\n\
                     This example works with the 9-bus case example case9mod_gen3_wind.m\n\n";

#include <opflow.h>
#include <exago_config.h>
#include <opflowselfcheck.h>
#include <utils.hpp>

int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  OPFLOW         opflow,opflow2;
  PS             ps,ps2;
  char           appname[]="opflow";
  MPI_Comm       comm=MPI_COMM_SELF;
  PetscReal      obj,obj2;
  PetscReal      pg1set,pg2set;

  /** Use `ExaGOLogSetLoggingFileName("opflow-logfile");` to log the output. */
  ierr = ExaGOInitialize(comm,&argc,&argv,appname,help);
  if (ierr)
  {
    fprintf(stderr,"Could not initialize ExaGO application %s.\n",appname);
    return ierr;
  }

  ExaGOLog(EXAGO_LOG_INFO,"%s","Creating OPFlow\n");

  
  /* Create the first OPFLOW object */
  ierr = OPFLOWCreate(comm,&opflow);CHKERRQ(ierr);
  
  ierr = OPFLOWReadMatPowerData(opflow,"../../../datafiles/case9/case9mod_gen3_wind.m");CHKERRQ(ierr);

  /* Set up its PS (network) object */
  ierr = OPFLOWSetUpPS(opflow);CHKERRQ(ierr);

  /* Access the PS object */
  ierr = OPFLOWGetPS(opflow,&ps);CHKERRQ(ierr);

  /* We want to set the dispatch of the wind generator which is at bus 3 */
  /* This is done by manipulating the generator real power limits. The max.
     capacity is set to the required gen. dispatch. The min. capacity is 
     set to 0 so that the generator can be curtailed, if needed.
  */
  /* The max. capacity set to 65 MW, min. capacity set to 0, reactive power limits retained to
     already set values
  */
  ierr = PSSetGenPowerLimits(ps,3,"1 ",65.0,0.0,EXAGO_IGNORE,EXAGO_IGNORE);CHKERRQ(ierr);
 
  /* Set up */
  ierr = OPFLOWSetUp(opflow);CHKERRQ(ierr);

  /* Solve */
  ierr = OPFLOWSolve(opflow);CHKERRQ(ierr);

  /* Print solution */
  ierr = OPFLOWPrintSolution(opflow);CHKERRQ(ierr);

  /* This function needs to be called after solving OPFLOW
     so that the OPFLOW solution is set in the PS object
  */
  ierr = OPFLOWSolutionToPS(opflow);CHKERRQ(ierr);
  
  /* Get the objective function and the dispatch for
     non renewable generators 
  */
  ierr = OPFLOWGetObjective(opflow,&obj);CHKERRQ(ierr);
  ierr = PSGetGenDispatch(ps,1,"1 ",&pg1set,NULL);CHKERRQ(ierr);
  ierr = PSGetGenDispatch(ps,2,"1 ",&pg2set,NULL);CHKERRQ(ierr);

  /* Destroy OPFLOW object */
  ierr = OPFLOWDestroy(&opflow);CHKERRQ(ierr);

  /*** Now create the second
  /* Create the first OPFLOW object */
  ierr = OPFLOWCreate(comm,&opflow2);CHKERRQ(ierr);
  
  ierr = OPFLOWReadMatPowerData(opflow2,"../../../datafiles/case9/case9mod_gen3_wind.m");CHKERRQ(ierr);

  /* Set up its PS (network) object */
  ierr = OPFLOWSetUpPS(opflow2);CHKERRQ(ierr);

  /* Access the PS object */
  ierr = OPFLOWGetPS(opflow2,&ps2);CHKERRQ(ierr);

  /* We want to set the dispatch of the wind generator which is at bus 3 */
  /* This is done by manipulating the generator real power limits. The max.
     capacity is set to the required gen. dispatch. The min. capacity is 
     set to 0 so that the generator can be curtailed, if needed.
  */
  /* The max. capacity set to 105 MW, min. capacity set to 0, reactive power limits retained to
     already set values
  */
  ierr = PSSetGenPowerLimits(ps2,3,"1 ",105.0,0.0,EXAGO_IGNORE,EXAGO_IGNORE);CHKERRQ(ierr);
  /* Set dispatch for non-renewable generation to the set-points. Here, we set the same min.
     and max. capacities so that the generation is held constant
  */
  ierr = PSSetGenPowerLimits(ps2,1,"1 ",pg1set,pg1set,EXAGO_IGNORE,EXAGO_IGNORE);CHKERRQ(ierr);
  ierr = PSSetGenPowerLimits(ps2,2,"1 ",pg2set,pg2set,EXAGO_IGNORE,EXAGO_IGNORE);CHKERRQ(ierr);

  /* Set up */
  ierr = OPFLOWSetUp(opflow2);CHKERRQ(ierr);

  /* Solve */
  ierr = OPFLOWSolve(opflow2);CHKERRQ(ierr);

  /* Print solution */
  ierr = OPFLOWPrintSolution(opflow2);CHKERRQ(ierr);

  /* This function needs to be called after solving OPFLOW
     so that the OPFLOW solution is set in the PS object
  */
  ierr = OPFLOWSolutionToPS(opflow2);CHKERRQ(ierr);
  
  /* Get the objective function
  */
  ierr = OPFLOWGetObjective(opflow2,&obj2);CHKERRQ(ierr);

  /* Destroy OPFLOW object */
  ierr = OPFLOWDestroy(&opflow2);CHKERRQ(ierr);

  ExaGOFinalize();
  return 0;
}
