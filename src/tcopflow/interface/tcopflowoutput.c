#include <private/opflowimpl.h>
#include <private/tcopflowimpl.h>

/*
  TCOPFLOWPrintSolution - Prints TCOPFLOW solution to stdout for the given time-step

  Input Parameters:
+ tcopflow - the TCOPFLOW object
- t_num - the time-step index (0 is the first time-step)
*/
PetscErrorCode TCOPFLOWPrintSolution(TCOPFLOW tcopflow,PetscInt t_num)
{
  PetscErrorCode ierr;
  PetscBool      conv_status;
  PetscReal      cost;
  OPFLOW         opflow=tcopflow->opflows[t_num];

  PetscFunctionBegin;
  if(!opflow->solutiontops) {
    ierr = TCOPFLOWGetSolution(tcopflow,t_num,&opflow->X);CHKERRQ(ierr);
    ierr = TCOPFLOWGetConstraintMultipliers(tcopflow,t_num,&opflow->Lambda);CHKERRQ(ierr);
    ierr = OPFLOWSolutionToPS(opflow);CHKERRQ(ierr);
  }

  /* Print to stdout */
  ierr = PetscPrintf(tcopflow->comm->type,"=============================================================\n");CHKERRQ(ierr);
  ierr = PetscPrintf(tcopflow->comm->type,"\tMulti-Period Optimal Power Flow\n");CHKERRQ(ierr);
  ierr = PetscPrintf(tcopflow->comm->type,"=============================================================\n");CHKERRQ(ierr);
  ierr = PetscPrintf(tcopflow->comm->type,"%-35s %s\n","OPFLOW Model",opflow->modelname);CHKERRQ(ierr);
  ierr = PetscPrintf(tcopflow->comm->type,"%-35s %s\n","Solver",tcopflow->solvername);CHKERRQ(ierr);
  /*  ierr = PetscPrintf(tcopflow->comm->type,"%-35s %s\n","Initialization",OPFLOWInitializationTypes[opflow->initializationtype]);CHKERRQ(ierr); */
  ierr = PetscPrintf(tcopflow->comm->type,"%-35s %-5.2f\n","Duration (minutes)",tcopflow->duration*60.0);CHKERRQ(ierr);
  ierr = PetscPrintf(tcopflow->comm->type,"%-35s %-5.2f\n","Time-step (minutes)",tcopflow->dT);CHKERRQ(ierr);
  ierr = PetscPrintf(tcopflow->comm->type,"%-35s %d\n","Number of steps",tcopflow->Nt);CHKERRQ(ierr);
  ierr = PetscPrintf(tcopflow->comm->type,"%-35s %s\n","Active power demand profile",tcopflow->ploadprofileset?tcopflow->ploadprofile:"NOT SET");CHKERRQ(ierr);
  ierr = PetscPrintf(tcopflow->comm->type,"%-35s %s\n","Rective power demand profile",tcopflow->qloadprofileset?tcopflow->qloadprofile:"NOT SET");CHKERRQ(ierr);
  ierr = PetscPrintf(tcopflow->comm->type,"%-35s %s\n","Wind generation profile",tcopflow->windgenprofileset?tcopflow->windgenprofile:"NOT SET");CHKERRQ(ierr);

  ierr = PetscPrintf(tcopflow->comm->type,"%-35s %s\n","Load loss allowed",opflow->include_loadloss_variables?"YES":"NO");CHKERRQ(ierr);
  if(opflow->include_loadloss_variables) {
    ierr = PetscPrintf(tcopflow->comm->type,"%-35s %g\n","Load loss penalty ($)",opflow->loadloss_penalty);CHKERRQ(ierr);
  }

  ierr = PetscPrintf(tcopflow->comm->type,"%-35s %s\n","Power imbalance allowed",opflow->include_powerimbalance_variables?"YES":"NO");CHKERRQ(ierr);
  if(opflow->include_powerimbalance_variables) {
    ierr = PetscPrintf(tcopflow->comm->type,"%-35s %g\n","Power imbalance penalty ($)",opflow->powerimbalance_penalty);CHKERRQ(ierr);
  }

  ierr = PetscPrintf(tcopflow->comm->type,"%-35s %s\n","Ignore line flow constraints",opflow->ignore_lineflow_constraints?"YES":"NO");CHKERRQ(ierr);
  ierr = PetscPrintf(tcopflow->comm->type,"\n");CHKERRQ(ierr);
  MPI_Barrier(tcopflow->comm->type);

  ierr = PetscPrintf(tcopflow->comm->type,"%-35s %d\n","Number of variables",tcopflow->Nx);CHKERRQ(ierr);
  ierr = PetscPrintf(tcopflow->comm->type,"%-35s %d\n","Number of equality constraints",tcopflow->Nconeq);CHKERRQ(ierr);
  ierr = PetscPrintf(tcopflow->comm->type,"%-35s %d\n","Number of inequality constraints",tcopflow->Nconineq);CHKERRQ(ierr);
  ierr = PetscPrintf(tcopflow->comm->type,"%-35s %d\n","Number of coupling constraints",tcopflow->Nconcoup);CHKERRQ(ierr);

  ierr = PetscPrintf(tcopflow->comm->type,"\n");CHKERRQ(ierr);

  ierr = TCOPFLOWGetConvergenceStatus(tcopflow,&conv_status);CHKERRQ(ierr);
  ierr = PetscPrintf(tcopflow->comm->type,"%-35s %s\n","Convergence status",conv_status?"CONVERGED":"DID NOT CONVERGE");CHKERRQ(ierr);
  ierr = TCOPFLOWGetObjective(tcopflow,&cost);CHKERRQ(ierr);
  ierr = PetscPrintf(tcopflow->comm->type,"%-35s %-7.2f\n","Objective value",cost);CHKERRQ(ierr);
  ierr = PetscPrintf(tcopflow->comm->type,"\n");CHKERRQ(ierr);

  MPI_Barrier(tcopflow->comm->type);

  ierr = PSPrintSystemSummary(opflow->ps);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  TCOPFLOWSaveSolution - Saves the TCOPFLOW solution for the given contingency to file

  Input Parameters:
+ tcopflow - the OPFLOW object
. t_num    - time-step number (0 for first step)
. format - the output file format (csv, matpower)
- outfile  - Name of output file
*/
PetscErrorCode TCOPFLOWSaveSolution(TCOPFLOW tcopflow,PetscInt t_num,OutputFormat format,const char outfile[])
{
  PetscErrorCode ierr;
  OPFLOW         opflow=tcopflow->opflows[t_num];

  PetscFunctionBegin;

  if(!opflow->solutiontops) {
    ierr = TCOPFLOWGetSolution(tcopflow,t_num,&opflow->X);CHKERRQ(ierr);
    ierr = TCOPFLOWGetConstraintMultipliers(tcopflow,t_num,&opflow->Lambda);CHKERRQ(ierr);
    ierr = OPFLOWSolutionToPS(opflow);CHKERRQ(ierr);
  }

  ierr = OPFLOWSaveSolution(opflow,format,outfile);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  TCOPFLOWSaveSolutionAll - Saves all TCOPFLOW solutions

  Input Parameters:
+ tcopflow - the OPFLOW object
. format - the output file format (csv, matpower)
- outdir  - Name of output directory

  Notes:
   Save all TCOPFLOW solutions (one contingency per file)
*/
PetscErrorCode TCOPFLOWSaveSolutionAll(TCOPFLOW tcopflow,OutputFormat format,const char outdir[])
{
  PetscErrorCode ierr;
  OPFLOW         opflow;
  PetscInt       i;
  char           filename[64];
  char           outfile[256];

  PetscFunctionBegin;

  ierr = PetscMkdir(outdir);CHKERRQ(ierr);

  for(i=0; i < tcopflow->Nt; i++) {
    opflow = tcopflow->opflows[i];
    if(!opflow->solutiontops) {
      ierr = TCOPFLOWGetSolution(tcopflow,i,&opflow->X);CHKERRQ(ierr);
      ierr = TCOPFLOWGetConstraintMultipliers(tcopflow,i,&opflow->Lambda);CHKERRQ(ierr);
      ierr = OPFLOWSolutionToPS(opflow);CHKERRQ(ierr);
    }
    ierr = PetscSNPrintf(filename,64,"t%d",i);CHKERRQ(ierr);
    ierr = PetscStrncpy(outfile,outdir,256);CHKERRQ(ierr);
    ierr = PetscStrlcat(outfile,"/",256);CHKERRQ(ierr);
    ierr = PetscStrlcat(outfile,filename,256);CHKERRQ(ierr);

    ierr = OPFLOWSaveSolution(opflow,format,outfile);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

