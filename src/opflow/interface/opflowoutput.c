#include <private/opflowimpl.h>

/*
  OPFLOWPrintSolution - Prints OPFLOW details to stdout.

  Input Parameters:
+ opflow - the OPFLOW object

*/
PetscErrorCode OPFLOWPrintSolution(OPFLOW opflow)
{
  PetscErrorCode ierr;
  PS    ps=(PS)opflow->ps;
  PSBUS bus;
  PSLOAD load;
  PetscBool   ghostbus;
  PetscInt    i,k;
  PetscScalar Pd,Qd;
  PetscBool   conv_status;
  PetscReal   cost;

  PetscFunctionBegin;
  if(!opflow->solutiontops) {
    ierr = OPFLOWSolutionToPS(opflow);CHKERRQ(ierr);
  }

  /* Print to stdout */
  ierr = PetscPrintf(opflow->comm->type,"=============================================================\n");CHKERRQ(ierr);
  ierr = PetscPrintf(opflow->comm->type,"\t\tOptimal Power Flow\n");CHKERRQ(ierr);
  ierr = PetscPrintf(opflow->comm->type,"=============================================================\n");CHKERRQ(ierr);
  ierr = PetscPrintf(opflow->comm->type,"%-35s %s\n","Formulation",opflow->formulationname);CHKERRQ(ierr);
  ierr = PetscPrintf(opflow->comm->type,"%-35s %s\n","Solver",opflow->solvername);CHKERRQ(ierr);
  ierr = PetscPrintf(opflow->comm->type,"%-35s %s\n","Initialization",OPFLOWInitializationTypes[opflow->initializationtype]);CHKERRQ(ierr);

  ierr = PetscPrintf(opflow->comm->type,"%-35s %s\n","Load loss allowed",opflow->include_loadloss_variables?"YES":"NO");CHKERRQ(ierr);
  if(opflow->include_loadloss_variables) {
    ierr = PetscPrintf(opflow->comm->type,"%-35s %g\n","Load loss penalty ($)",opflow->loadloss_penalty);CHKERRQ(ierr);
  }

  ierr = PetscPrintf(opflow->comm->type,"%-35s %s\n","Power imbalance allowed",opflow->include_powerimbalance_variables?"YES":"NO");CHKERRQ(ierr);
  if(opflow->include_powerimbalance_variables) {
    ierr = PetscPrintf(opflow->comm->type,"%-35s %g\n","Power imbalance penalty ($)",opflow->powerimbalance_penalty);CHKERRQ(ierr);
  }

  ierr = PetscPrintf(opflow->comm->type,"%-35s %s\n","Ignore line flow constraints",opflow->ignore_lineflow_constraints?"YES":"NO");CHKERRQ(ierr);
  ierr = PetscPrintf(opflow->comm->type,"\n");CHKERRQ(ierr);
  MPI_Barrier(opflow->comm->type);

  ierr = PetscPrintf(opflow->comm->type,"%-35s %d\n","Number of variables",opflow->Nx);CHKERRQ(ierr);
  ierr = PetscPrintf(opflow->comm->type,"%-35s %d\n","Number of equality constraints",opflow->Nconeq);CHKERRQ(ierr);
  ierr = PetscPrintf(opflow->comm->type,"%-35s %d\n","Number of inequality constraints",opflow->Nconineq);CHKERRQ(ierr);
  ierr = PetscPrintf(opflow->comm->type,"\n");CHKERRQ(ierr);

  ierr = OPFLOWGetConvergenceStatus(opflow,&conv_status);CHKERRQ(ierr);
  ierr = PetscPrintf(opflow->comm->type,"%-35s %s\n","Convergence status",conv_status?"CONVERGED":"DID NOT CONVERGE");CHKERRQ(ierr);
  ierr = OPFLOWGetObjective(opflow,&cost);CHKERRQ(ierr);
  ierr = PetscPrintf(opflow->comm->type,"%-35s %-7.2f\n","Objective value",cost);CHKERRQ(ierr);
  ierr = PetscPrintf(opflow->comm->type,"\n");CHKERRQ(ierr);

  MPI_Barrier(opflow->comm->type);

  ierr = PetscPrintf(opflow->comm->type,"----------------------------------------------------------------------\n");CHKERRQ(ierr);
  ierr = PetscPrintf(opflow->comm->type,"%-10s %-7s %-7s %-7s %-7s %-14s %-14s\n","Bus","Pd","Qd","Vm","Va","mult_Pmis","mult_Qmis");CHKERRQ(ierr);
  ierr = PetscPrintf(opflow->comm->type,"----------------------------------------------------------------------\n");CHKERRQ(ierr);

  MPI_Barrier(opflow->comm->type);

  for(i=0;i < ps->nbus; i++) {
    bus = &ps->bus[i];

    ierr = PSBUSIsGhosted(bus,&ghostbus);CHKERRQ(ierr);
    if(ghostbus) continue;

    Pd = Qd = 0.0;
    for(k=0; k < bus->nload; k++) {
      ierr = PSBUSGetLoad(bus,k,&load);CHKERRQ(ierr);
      Pd += load->pl*ps->MVAbase;
      Qd += load->ql*ps->MVAbase;
    }
    ierr = PetscPrintf(PETSC_COMM_SELF,"%-6d %7.2f %7.2f %7.2f %7.2f %12.2f %12.2f\n",bus->bus_i,Pd,Qd,bus->vm,bus->va*180.0/PETSC_PI,bus->mult_pmis,bus->mult_qmis);CHKERRQ(ierr);
  }
  ierr = PetscPrintf(opflow->comm->type,"\n");CHKERRQ(ierr);
  MPI_Barrier(opflow->comm->type);

  ierr = PetscPrintf(opflow->comm->type,"----------------------------------------------------------------------------------------\n");CHKERRQ(ierr);
  ierr = PetscPrintf(opflow->comm->type,"%-10s %-8s %-10s %-8s %-7s %-8s %-8s %-8s\n","From","To","Status","Sft","Stf","Slim","mult_Sf","mult_St");CHKERRQ(ierr);
  ierr = PetscPrintf(opflow->comm->type,"----------------------------------------------------------------------------------------\n");CHKERRQ(ierr);
  MPI_Barrier(opflow->comm->type);
  for(i=0; i < ps->nline; i++) {
    PSLINE line;
    line = &ps->line[i];
    ierr = PetscPrintf(PETSC_COMM_SELF,"%-10d %-10d %-4d %8.2f %8.2f %8.2f %8.2f %8.2f\n",line->fbus,line->tbus,line->status,line->sf*ps->MVAbase,line->st*ps->MVAbase,line->rateA,line->mult_sf,line->mult_st);CHKERRQ(ierr);
  }
  ierr = PetscPrintf(opflow->comm->type,"\n");CHKERRQ(ierr);
  MPI_Barrier(opflow->comm->type);

  ierr = PetscPrintf(opflow->comm->type,"----------------------------------------------------------------------------------------\n");CHKERRQ(ierr);
  ierr = PetscPrintf(opflow->comm->type,"%-8s %-10s %-8s %-8s %-8s %-8s %-8s %-8s %-6s\n","Gen","Status","Fuel","Pg","Qg","Pmin","Pmax","Qmin","Qmax");CHKERRQ(ierr);
  ierr = PetscPrintf(opflow->comm->type,"----------------------------------------------------------------------------------------\n");CHKERRQ(ierr);
  MPI_Barrier(opflow->comm->type);
  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    ierr = PSBUSIsGhosted(bus,&ghostbus);CHKERRQ(ierr);
    if(ghostbus) continue;

    for(k=0; k < bus->ngen; k++) {
      PSGEN gen;
      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_SELF,"%-10d %-4d %8s %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n",bus->bus_i,gen->status,PSGENFuelTypes[gen->genfuel_type],gen->pg*ps->MVAbase,gen->qg*ps->MVAbase,gen->pb*ps->MVAbase,gen->pt*ps->MVAbase,gen->qb*ps->MVAbase,gen->qt*ps->MVAbase);CHKERRQ(ierr);
    }
  }
  MPI_Barrier(opflow->comm->type);


  PetscFunctionReturn(0);
}

/* Save to CSV format */
PetscErrorCode OPFLOWSaveSolution_CSV(OPFLOW opflow,const char outfile[])
{
  PetscErrorCode ierr;
  FILE           *fp;
  PS             ps=opflow->ps;
  PSBUS          bus;
  PSLOAD         load;
  PSLINE         line;
  PetscBool      ghostbus;
  PetscInt       i,k;
  PetscScalar    Pd,Qd;
  PetscReal      cost;

  PetscFunctionBegin;

  fp = fopen(outfile,"w");
  if (fp == NULL) {
    SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Cannot open OPFLOW output file %s",outfile);CHKERRQ(ierr);
  }
  fprintf(fp,"nbus,ngen,nbranch,baseMVA\n");
  fprintf(fp,"%d,%d,%d,%3.2f\n\n",ps->nbus,ps->ngen,ps->nline,ps->MVAbase);

  fprintf(fp,"Bus,Pd,Qd,Vm,Va,mult_Pmis,mult_Qmis\n");

  for(i=0;i < ps->nbus; i++) {
    bus = &ps->bus[i];

    ierr = PSBUSIsGhosted(bus,&ghostbus);CHKERRQ(ierr);
    if(ghostbus) continue;

    Pd = Qd = 0.0;
    for(k=0; k < bus->nload; k++) {
      ierr = PSBUSGetLoad(bus,k,&load);CHKERRQ(ierr);
      Pd += load->pl*ps->MVAbase;
      Qd += load->ql*ps->MVAbase;
    }
    fprintf(fp,"%6d, %7.2f, %7.2f, %7.2f, %7.2f, %12.2f, %12.2f\n",bus->bus_i,Pd,Qd,bus->vm,bus->va*180.0/PETSC_PI,bus->mult_pmis,bus->mult_qmis);
  }
  fprintf(fp,"\n");

  fprintf(fp,"From,To,Status,Sft,Stf,Slim,mult_Sf,mult_St\n");
  for(i=0; i < ps->nline; i++) {
    line = &ps->line[i];
    fprintf(fp,"%10d, %10d, %4d, %8.2f, %8.2f, %8.2f, %8.2f, %8.2f\n",line->fbus,line->tbus,line->status,line->sf*ps->MVAbase,line->st*ps->MVAbase,line->rateA,line->mult_sf,line->mult_st);CHKERRQ(ierr);
  }

  fprintf(fp,"Gen bus,Status,Fuel,Pg,Qg,Pmin,Pmax,Qmin,Qmax\n");

  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    ierr = PSBUSIsGhosted(bus,&ghostbus);CHKERRQ(ierr);
    if(ghostbus) continue;

    for(k=0; k < bus->ngen; k++) {
      PSGEN gen;
      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      fprintf(fp,"%6d, %4d, %8s, %8.2f, %8.2f, %8.2f, %8.2f, %8.2f, %8.2f\n",bus->bus_i,gen->status,PSGENFuelTypes[gen->genfuel_type],gen->pg*ps->MVAbase,gen->qg*ps->MVAbase,gen->pb*ps->MVAbase,gen->pt*ps->MVAbase,gen->qb*ps->MVAbase,gen->qt*ps->MVAbase);CHKERRQ(ierr);
    }
  }

  fclose(fp);


  PetscFunctionReturn(0);
}

/* Save to MATPOWER format */
PetscErrorCode OPFLOWSaveSolution_MATPOWER(OPFLOW opflow,const char outfile[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

/*
  OPFLOWSaveSolution - Saves the OPFLOW solution to file

  Input Parameters:
+ opflow - the OPFLOW object
. format - the output file format (csv, matpower)
- outfile  - Name of output file
*/
PetscErrorCode OPFLOWSaveSolution(OPFLOW opflow,OutputFormat format,const char outfile[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  if(!opflow->solutiontops) {
    ierr = OPFLOWSolutionToPS(opflow);CHKERRQ(ierr);
  }

  if(format == CSV) {
    ierr = OPFLOWSaveSolution_CSV(opflow,outfile);CHKERRQ(ierr);
  } else if(format == MATPOWER) {
    ierr = OPFLOWSaveSolution_MATPOWER(opflow,outfile);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}
