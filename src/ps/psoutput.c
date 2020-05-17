#include <private/psimpl.h>

/* Save to MATPOWER format */
PetscErrorCode PSSaveSolution_MATPOWER(PS ps,const char outfile[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Saving file to MATPOWER format not supported yet\n");

  PetscFunctionReturn(0);
}

/* Save to CSV format */
PetscErrorCode PSSaveSolution_CSV(PS ps,const char outfile[])
{
  PetscErrorCode ierr;
  FILE           *fp;
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
  fprintf(fp,"%d,%d,%d,%3.2f\n",ps->nbus,ps->ngen,ps->nline,ps->MVAbase);

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

  fprintf(fp,"From,To,Status,Sft,Stf,Slim,mult_Sf,mult_St\n");
  for(i=0; i < ps->nline; i++) {
    line = &ps->line[i];
    fprintf(fp,"%10d, %10d, %4d, %8.2f, %8.2f, %8.2f, %8.2f, %8.2f\n",line->fbus,line->tbus,line->status,line->sf*ps->MVAbase,line->st*ps->MVAbase,(line->rateA>1e5)?10000:line->rateA,line->mult_sf,line->mult_st);CHKERRQ(ierr);
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

/*
  PSSaveSolution - Saves the system solution to file

  Input Parameters:
+ ps - the PS object
. format - the output file format (csv, matpower)
- outfile  - Name of output file
*/
PetscErrorCode PSSaveSolution(PS ps,OutputFormat format,const char outfile[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  if(format == CSV) {
    ierr = PSSaveSolution_CSV(ps,outfile);CHKERRQ(ierr);
  } else if(format == MATPOWER) {
    ierr = PSSaveSolution_MATPOWER(ps,outfile);CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

/* 
   PSPrintSystemSummary - Prints the bus,gen, and branch data

   Input Parameters:
+  ps - PS object
*/
PetscErrorCode PSPrintSystemSummary(PS ps)
{
  PetscErrorCode ierr;
  PSBUS       bus;
  PSLOAD      load;
  PSGEN       gen;
  PSLINE      line;
  PetscBool   ghostbus;
  PetscInt    i,k;
  PetscScalar Pd,Qd;
  PetscBool   conv_status;
  PetscReal   cost;

  PetscFunctionBegin;

  ierr = PetscPrintf(ps->comm->type,"----------------------------------------------------------------------\n");CHKERRQ(ierr);
  ierr = PetscPrintf(ps->comm->type,"%-10s %-7s %-7s %-7s %-7s %-14s %-14s\n","Bus","Pd","Qd","Vm","Va","mult_Pmis","mult_Qmis");CHKERRQ(ierr);
  ierr = PetscPrintf(ps->comm->type,"----------------------------------------------------------------------\n");CHKERRQ(ierr);

  MPI_Barrier(ps->comm->type);

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
    ierr = PetscPrintf(ps->comm->type,"%-6d %7.2f %7.2f %7.2f %7.2f %12.2f %12.2f\n",bus->bus_i,Pd,Qd,bus->vm,bus->va*180.0/PETSC_PI,bus->mult_pmis,bus->mult_qmis);CHKERRQ(ierr);
  }
  ierr = PetscPrintf(ps->comm->type,"\n");CHKERRQ(ierr);
  MPI_Barrier(ps->comm->type);

  ierr = PetscPrintf(ps->comm->type,"----------------------------------------------------------------------------------------\n");CHKERRQ(ierr);
  ierr = PetscPrintf(ps->comm->type,"%-10s %-8s %-10s %-8s %-7s %-8s %-8s %-8s\n","From","To","Status","Sft","Stf","Slim","mult_Sf","mult_St");CHKERRQ(ierr);
  ierr = PetscPrintf(ps->comm->type,"----------------------------------------------------------------------------------------\n");CHKERRQ(ierr);
  MPI_Barrier(ps->comm->type);
  for(i=0; i < ps->nline; i++) {
    line = &ps->line[i];
    ierr = PetscPrintf(ps->comm->type,"%-10d %-10d %-4d %8.2f %8.2f %8.2f %8.2f %8.2f\n",line->fbus,line->tbus,line->status,line->sf*ps->MVAbase,line->st*ps->MVAbase,(line->rateA>1e5)?10000:line->rateA,line->mult_sf,line->mult_st);CHKERRQ(ierr);
  }
  ierr = PetscPrintf(ps->comm->type,"\n");CHKERRQ(ierr);
  MPI_Barrier(ps->comm->type);

  ierr = PetscPrintf(ps->comm->type,"----------------------------------------------------------------------------------------\n");CHKERRQ(ierr);
  ierr = PetscPrintf(ps->comm->type,"%-8s %-10s %-8s %-8s %-8s %-8s %-8s %-8s %-6s\n","Gen","Status","Fuel","Pg","Qg","Pmin","Pmax","Qmin","Qmax");CHKERRQ(ierr);
  ierr = PetscPrintf(ps->comm->type,"----------------------------------------------------------------------------------------\n");CHKERRQ(ierr);
  MPI_Barrier(ps->comm->type);
  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    ierr = PSBUSIsGhosted(bus,&ghostbus);CHKERRQ(ierr);
    if(ghostbus) continue;

    for(k=0; k < bus->ngen; k++) {
      ierr = PSBUSGetGen(bus,k,&gen);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_SELF,"%-10d %-4d %8s %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n",bus->bus_i,gen->status,PSGENFuelTypes[gen->genfuel_type],gen->pg*ps->MVAbase,gen->qg*ps->MVAbase,gen->pb*ps->MVAbase,gen->pt*ps->MVAbase,gen->qb*ps->MVAbase,gen->qt*ps->MVAbase);CHKERRQ(ierr);
    }
  }
  MPI_Barrier(ps->comm->type);

  PetscFunctionReturn(0);
}
