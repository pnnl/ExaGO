#include <private/opflowimpl.h>
#include <private/scopflowimpl.h>

/*
  SCOPFLOWPrintSolution - Prints SCOPFLOW solution to stdout for the given contingency

  Input Parameters:
+ scopflow - the SCOPFLOW object
- cont_num - the contingency number (0 for base)
  Notes:
   Only prints out the system summary for base case

*/
PetscErrorCode SCOPFLOWPrintSolution(SCOPFLOW scopflow,PetscInt cont_num)
{
  PetscErrorCode ierr;
  PetscBool      conv_status;
  PetscReal      cost;
  OPFLOW         opflow=scopflow->opflows[cont_num];

  PetscFunctionBegin;
  if(!opflow->solutiontops) {
    ierr = SCOPFLOWGetSolution(scopflow,cont_num,&opflow->X);CHKERRQ(ierr);
    ierr = SCOPFLOWGetConstraintMultipliers(scopflow,cont_num,&opflow->Lambda);CHKERRQ(ierr);
    ierr = OPFLOWSolutionToPS(opflow);CHKERRQ(ierr);
  }

  /* Print to stdout */
  ierr = PetscPrintf(scopflow->comm->type,"=============================================================\n");CHKERRQ(ierr);
  ierr = PetscPrintf(scopflow->comm->type,"\tSecurity-Constrained Optimal Power Flow\n");CHKERRQ(ierr);
  ierr = PetscPrintf(scopflow->comm->type,"=============================================================\n");CHKERRQ(ierr);
  ierr = PetscPrintf(scopflow->comm->type,"%-35s %s\n","OPFLOW Formulation",opflow->modelname);CHKERRQ(ierr);
  ierr = PetscPrintf(scopflow->comm->type,"%-35s %s\n","Solver",scopflow->solvername);CHKERRQ(ierr);
  ierr = PetscPrintf(scopflow->comm->type,"%-35s %s\n","Initialization",OPFLOWInitializationTypes[opflow->initializationtype]);CHKERRQ(ierr);
  ierr = PetscPrintf(scopflow->comm->type,"%-35s %s\n","Operation Mode",scopflow->mode?"Corrective":"Preventive");CHKERRQ(ierr);
  ierr = PetscPrintf(scopflow->comm->type,"%-35s %d\n","Number of contingencies",scopflow->Nc-1);CHKERRQ(ierr);

  ierr = PetscPrintf(scopflow->comm->type,"%-35s %s\n","Load loss allowed",opflow->include_loadloss_variables?"YES":"NO");CHKERRQ(ierr);
  if(opflow->include_loadloss_variables) {
    ierr = PetscPrintf(scopflow->comm->type,"%-35s %g\n","Load loss penalty ($)",opflow->loadloss_penalty);CHKERRQ(ierr);
  }

  ierr = PetscPrintf(scopflow->comm->type,"%-35s %s\n","Power imbalance allowed",opflow->include_powerimbalance_variables?"YES":"NO");CHKERRQ(ierr);
  if(opflow->include_powerimbalance_variables) {
    ierr = PetscPrintf(scopflow->comm->type,"%-35s %g\n","Power imbalance penalty ($)",opflow->powerimbalance_penalty);CHKERRQ(ierr);
  }

  ierr = PetscPrintf(scopflow->comm->type,"%-35s %s\n","Ignore line flow constraints",opflow->ignore_lineflow_constraints?"YES":"NO");CHKERRQ(ierr);
  ierr = PetscPrintf(scopflow->comm->type,"\n");CHKERRQ(ierr);
  MPI_Barrier(scopflow->comm->type);

  ierr = PetscPrintf(scopflow->comm->type,"%-35s %d\n","Number of variables",scopflow->Nx);CHKERRQ(ierr);
  ierr = PetscPrintf(scopflow->comm->type,"%-35s %d\n","Number of equality constraints",scopflow->Nconeq);CHKERRQ(ierr);
  ierr = PetscPrintf(scopflow->comm->type,"%-35s %d\n","Number of inequality constraints",scopflow->Nconineq);CHKERRQ(ierr);
  ierr = PetscPrintf(scopflow->comm->type,"%-35s %d\n","Number of coupling constraints",scopflow->Nconcoup);CHKERRQ(ierr);

  ierr = PetscPrintf(scopflow->comm->type,"\n");CHKERRQ(ierr);

  ierr = SCOPFLOWGetConvergenceStatus(scopflow,&conv_status);CHKERRQ(ierr);
  ierr = PetscPrintf(scopflow->comm->type,"%-35s %s\n","Convergence status",conv_status?"CONVERGED":"DID NOT CONVERGE");CHKERRQ(ierr);
  ierr = SCOPFLOWGetObjective(scopflow,&cost);CHKERRQ(ierr);
  ierr = PetscPrintf(scopflow->comm->type,"%-35s %-7.2f\n","Objective value",cost);CHKERRQ(ierr);
  ierr = PetscPrintf(scopflow->comm->type,"\n");CHKERRQ(ierr);

  MPI_Barrier(scopflow->comm->type);

  ierr = PSPrintSystemSummary(opflow->ps);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  SCOPFLOWSaveSolution - Saves the SCOPFLOW solution for the given contingency to file

  Input Parameters:
+ scopflow - the OPFLOW object
. cont_num - the contingency number (0 for base case)
. format - the output file format (csv, matpower)
- outfile  - Name of output file
*/
PetscErrorCode SCOPFLOWSaveSolution(SCOPFLOW scopflow,PetscInt cont_num,OutputFormat format,const char outfile[])
{
  PetscErrorCode ierr;
  OPFLOW         opflow=scopflow->opflows[cont_num];

  PetscFunctionBegin;

  if(!opflow->solutiontops) {
    ierr = SCOPFLOWGetSolution(scopflow,cont_num,&opflow->X);CHKERRQ(ierr);
    ierr = SCOPFLOWGetConstraintMultipliers(scopflow,cont_num,&opflow->Lambda);CHKERRQ(ierr);
    ierr = OPFLOWSolutionToPS(opflow);CHKERRQ(ierr);
  }

  ierr = OPFLOWSaveSolution(opflow,format,outfile);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
  SCOPFLOWSaveSolutionAll - Saves all SCOPFLOW solutions

  Input Parameters:
+ scopflow - the OPFLOW object
. format - the output file format (csv, matpower)
- outdir  - Name of output directory

  Notes:
   Save all SCOPFLOW solutions (one contingency per file)
*/
PetscErrorCode SCOPFLOWSaveSolutionAll(SCOPFLOW scopflow,OutputFormat format,const char outdir[])
{
  PetscErrorCode ierr;
  OPFLOW         opflow;
  PetscInt       i;
  char           filename[64];
  char           outfile[256];

  PetscFunctionBegin;

  ierr = PetscMkdir(outdir);CHKERRQ(ierr);

  for(i=0; i < scopflow->Nc; i++) {
    opflow = scopflow->opflows[i];
    if(!opflow->solutiontops) {
      ierr = SCOPFLOWGetSolution(scopflow,i,&opflow->X);CHKERRQ(ierr);
      ierr = SCOPFLOWGetConstraintMultipliers(scopflow,i,&opflow->Lambda);CHKERRQ(ierr);
      ierr = OPFLOWSolutionToPS(opflow);CHKERRQ(ierr);
    }
    ierr = PetscSNPrintf(filename,64,"cont_%d",i);CHKERRQ(ierr);
    ierr = PetscStrncpy(outfile,outdir,256);CHKERRQ(ierr);
    ierr = PetscStrlcat(outfile,"/",256);CHKERRQ(ierr);
    ierr = PetscStrlcat(outfile,filename,256);CHKERRQ(ierr);

    ierr = OPFLOWSaveSolution(opflow,format,outfile);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

