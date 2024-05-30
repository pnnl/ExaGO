#include <exago_config.h>
#include <petsc/private/dmnetworkimpl.h>
#include <private/opflowimpl.h>
#include <private/pflowimpl.h>

/*
  OPFLOWCheckConstraints - Displays information about OPFLOW constraints

  Input Parameters:
. opflow - the OPFLOW oject

  Notes: This is called by OPFLOWSolutionToPS. If called exterenally then OPFLOWSolutionToPS() needs to be called first.
*/
PetscErrorCode OPFLOWCheckConstraints(OPFLOW opflow)
{
  PetscErrorCode ierr;
  PS             ps=opflow->ps;
  PSBUS          bus;
  PSLINE         line;
  PSGEN          gen;
  PSLOAD         load;
  PetscInt       i,j;
  
  PetscFunctionBegin;

  if(!opflow->solutiontops) {
    ierr = OPFLOWSolutionToPS(opflow);
    CHKERRQ(ierr);
  }

  if(opflow->modelops.checkconstraints) {
    ierr = (*opflow->modelops.checkconstraints)(opflow);
    CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}
