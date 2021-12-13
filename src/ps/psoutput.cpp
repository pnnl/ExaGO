#include <private/psimpl.h>

/* Save to MATPOWER format */
PetscErrorCode PSSaveSolution_MATPOWER(PS ps, const char outfile[]) {
  PetscErrorCode ierr;
  FILE *fd;
  const char *prefix = "mpc.";
  PSBUS bus;
  PSLOAD load;
  PSGEN gen;
  PSLINE line;
  PetscScalar Pd, Qd;
  PetscInt i, k;
  PetscScalar MVAbase = ps->MVAbase;
  char filename[PETSC_MAX_PATH_LEN];
  char fcn_name[100];
  char *tok, *tok2;
  char sep[] = "/";
  char ext[] = ".m";
  char file1[PETSC_MAX_PATH_LEN];

  PetscFunctionBegin;
  strcpy(file1, outfile);
  /* Check if file has .m extension */
  tok = strtok(file1, ext);
  strcpy(filename, tok);

  strcpy(file1, filename);
  tok2 = strtok(file1, sep);
  strcpy(fcn_name, file1);
  while ((tok2 = strtok(NULL, sep)) != NULL) {
    strcpy(fcn_name, tok2);
  }

  /* Add .m extension to file name */
  ierr = PetscStrlcat(filename, ".m", 256);
  CHKERRQ(ierr);

  fd = fopen(filename, "w");
  if (fd == NULL) {
    SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN,
             "Cannot open OPFLOW output file %s", outfile);
    CHKERRQ(ierr);
  }

  /* Function header */
  fprintf(fd, "function mpc = %s\n", fcn_name);

  /* Write system MVAbase */
  fprintf(fd, "\n%%%%-----  Power Flow Data  -----%%%%\n");
  fprintf(fd, "%%%% system MVA base\n");
  fprintf(fd, "%sbaseMVA = %.9g;\n", prefix, ps->MVAbase);

  /* Write OPFLOW objective function value */
  /* Note: MATPOWER data files do not store objective function values.*/
  fprintf(fd, "\n%%%% OPF objective\n");
  fprintf(fd, "%sobj = %.9g;\n", prefix, ps->opflowobj);

  /* Write bus data */
  fprintf(fd, "\n%%%% bus data\n");
  fprintf(fd, "%%\tbus_"
              "i\ttype\tPd\tQd\tGs\tBs\tarea\tVm\tVa\tbaseKV\tzone\tVmax\tVmin"
              "\tmult_Pmis\tmult_Qmis\tPslack\tQslack");
  fprintf(fd, "\n%sbus = [\n", prefix);
  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    /* Get total load on the bus */
    Pd = Qd = 0.0;
    for (k = 0; k < bus->nload; k++) {
      ierr = PSBUSGetLoad(bus, k, &load);
      CHKERRQ(ierr);
      Pd += load->pl * ps->MVAbase;
      Qd += load->ql * ps->MVAbase;
    }
    fprintf(fd,
            "\t%d\t%d\t%.9g\t%.9g\t%.9g\t%.9g\t%d\t%.9g\t%.9g\t%.9g\t%d\t%."
            "9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g;\n",
            bus->bus_i, bus->ide, Pd, Qd, bus->gl * MVAbase, bus->bl * MVAbase,
            bus->area, bus->vm, bus->va * 180.0 / PETSC_PI, bus->basekV,
            bus->zone, bus->Vmax, bus->Vmin, bus->mult_pmis, bus->mult_qmis,
            bus->pimb * ps->MVAbase, bus->qimb * ps->MVAbase);
  }
  fprintf(fd, "];\n");

  /* Write generator data */
  fprintf(fd, "\n%%%% generator data\n");
  fprintf(fd, "%%\tbus\tPg\tQg\tQmax\tQmin\tVg\tmBase\tstatus\tPmax\tPmin");
  fprintf(fd, "\tPc1\tPc2\tQc1min\tQc1max\tQc2min\tQc2max\tramp_agc\tramp_"
              "10\tramp_30\tramp_q\tapf\tpgs\n");
  fprintf(fd, "%sgen = [\n", prefix);
  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    for (k = 0; k < bus->ngen; k++) {
      ierr = PSBUSGetGen(bus, k, &gen);
      CHKERRQ(ierr);
      fprintf(
          fd,
          "\t%d\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%d\t%.9g\t%.9g\t%.9g\t%."
          "9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g;\n",
          bus->bus_i, gen->pg * MVAbase, gen->qg * MVAbase, gen->qt * MVAbase,
          gen->qb * MVAbase, bus->vm, gen->mbase, gen->status,
          gen->pt * MVAbase, gen->pb * MVAbase, gen->pc1 * MVAbase,
          gen->pc2 * MVAbase, gen->qc1min * MVAbase, gen->qc1max * MVAbase,
          gen->qc2min * MVAbase, gen->qc2max * MVAbase,
          gen->ramp_rate_min * MVAbase, gen->ramp_rate_10min * MVAbase,
          gen->ramp_rate_30min * MVAbase, gen->ramp_rate_min_mvar * MVAbase,
          gen->apf, gen->pgs * MVAbase);
    }
  }
  fprintf(fd, "];\n");

  /* Write branch data */
  fprintf(fd, "\n%%%% branch data\n");
  fprintf(fd,
          "%%\tfbus\ttbus\tr\tx\tb\trateA\trateB\trateC\tratio\tangle\tstatus");
  fprintf(fd, "\tangmin\tangmax\tPf\tQf\tPt\tQt\n");
  fprintf(fd, "%sbranch = [\n", prefix);
  for (i = 0; i < ps->nline; i++) {
    line = &ps->line[i];
    fprintf(fd,
            "\t%d\t%d\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%d\t%."
            "9g\t%.9g\t%.4f\t%.4f\t%.4f\t%.4f;\n",
            line->fbus, line->tbus, line->r, line->x, line->b, line->rateA,
            line->rateB, line->rateC, line->tapratio, line->phaseshift,
            line->status, 0.0, 0.0, line->pf * MVAbase, line->qf * MVAbase,
            line->pt * MVAbase, line->qt * MVAbase);
  }
  fprintf(fd, "];\n");

  /* Generator cost data */
  fprintf(fd, "\n%%%% generator cost data\n");
  fprintf(fd, "%%\t2\tstartup\tshutdown\tn\tc(n-1)\t...\tc0\n");
  fprintf(fd, "%%\tUsing quadratic cost curves only\n");
  fprintf(fd, "%sgencost = [\n", prefix);
  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    for (k = 0; k < bus->ngen; k++) {
      ierr = PSBUSGetGen(bus, k, &gen);
      CHKERRQ(ierr);
      ierr = fprintf(fd, "\t%d\t%.9g\t%.9g\t%d\t%.9g\t%.9g\t%.9g;\n", 2,
                     gen->cost_startup, gen->cost_shutdown, 3, gen->cost_alpha,
                     gen->cost_beta, gen->cost_gamma);
    }
  }
  fprintf(fd, "];\n");

  fclose(fd);
  PetscFunctionReturn(0);
}

/* Save to CSV format */
PetscErrorCode PSSaveSolution_CSV(PS ps, const char outfile[]) {
  PetscErrorCode ierr;
  FILE *fp;
  PSBUS bus;
  PSLOAD load;
  PSLINE line;
  PetscBool ghostbus;
  PetscInt i, k;
  PetscScalar Pd, Qd;

  PetscFunctionBegin;

  fp = fopen(outfile, "w");
  if (fp == NULL) {
    SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN,
             "Cannot open OPFLOW output file %s", outfile);
    CHKERRQ(ierr);
  }
  fprintf(fp, "nbus,ngen,nbranch,baseMVA\n");
  fprintf(fp, "%d,%d,%d,%3.2f\n", ps->nbus, ps->ngen, ps->nline, ps->MVAbase);

  fprintf(fp, "Bus,Pd,Qd,Vm,Va,mult_Pmis,mult_Qmis\n");

  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    ierr = PSBUSIsGhosted(bus, &ghostbus);
    CHKERRQ(ierr);
    if (ghostbus)
      continue;

    Pd = Qd = 0.0;
    for (k = 0; k < bus->nload; k++) {
      ierr = PSBUSGetLoad(bus, k, &load);
      CHKERRQ(ierr);
      Pd += load->pl * ps->MVAbase;
      Qd += load->ql * ps->MVAbase;
    }
    fprintf(fp, "%6d, %7.2f, %7.2f, %7.2f, %7.2f, %12.2f, %12.2f\n", bus->bus_i,
            Pd, Qd, bus->vm, bus->va * 180.0 / PETSC_PI, bus->mult_pmis,
            bus->mult_qmis);
  }

  fprintf(fp, "From,To,Status,Sft,Stf,Slim,mult_Sf,mult_St\n");
  for (i = 0; i < ps->nline; i++) {
    line = &ps->line[i];
    fprintf(fp, "%10d, %10d, %4d, %8.2f, %8.2f, %8.2f, %8.2f, %8.2f\n",
            line->fbus, line->tbus, line->status, line->sf * ps->MVAbase,
            line->st * ps->MVAbase, (line->rateA > 1e5) ? 10000 : line->rateA,
            line->mult_sf, line->mult_st);
    CHKERRQ(ierr);
  }

  fprintf(fp, "Gen bus,Status,Fuel,Pg,Qg,Pmin,Pmax,Qmin,Qmax\n");

  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    ierr = PSBUSIsGhosted(bus, &ghostbus);
    CHKERRQ(ierr);
    if (ghostbus)
      continue;

    for (k = 0; k < bus->ngen; k++) {
      PSGEN gen;
      ierr = PSBUSGetGen(bus, k, &gen);
      CHKERRQ(ierr);
      fprintf(fp, "%6d, %4d, %8s, %8.2f, %8.2f, %8.2f, %8.2f, %8.2f, %8.2f\n",
              bus->bus_i, gen->status, PSGENFuelTypes[gen->genfuel_type],
              gen->pg * ps->MVAbase, gen->qg * ps->MVAbase,
              gen->pb * ps->MVAbase, gen->pt * ps->MVAbase,
              gen->qb * ps->MVAbase, gen->qt * ps->MVAbase);
      CHKERRQ(ierr);
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
PetscErrorCode PSSaveSolution(PS ps, OutputFormat format,
                              const char outfile[]) {
  PetscErrorCode ierr;

  PetscFunctionBegin;

  if (format == CSV) {
    ierr = PSSaveSolution_CSV(ps, outfile);
    CHKERRQ(ierr);
  } else if (format == MATPOWER) {
    ierr = PSSaveSolution_MATPOWER(ps, outfile);
    CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

/*
   PSPrintSystemSummary - Prints the bus,gen, and branch data

   Input Parameters:
+  ps - PS object
*/
PetscErrorCode PSPrintSystemSummary(PS ps) {
  PetscErrorCode ierr;
  PSBUS bus;
  PSLOAD load;
  PSGEN gen;
  PSLINE line;
  PetscBool ghostbus;
  PetscInt i, k;
  PetscScalar Pd, Qd;

  PetscFunctionBegin;

  ierr = PetscPrintf(ps->comm->type, "-----------------------------------------"
                                     "-----------------------------\n");
  CHKERRQ(ierr);
  ierr = PetscPrintf(ps->comm->type,
                     "%-10s %-7s %-7s %-7s %-7s %-14s %-14s %-14s %-14s\n",
                     "Bus", "Pd", "Qd", "Vm", "Va", "mult_Pmis", "mult_Qmis",
                     "Pslack", "Qslack");
  CHKERRQ(ierr);
  ierr = PetscPrintf(ps->comm->type, "-----------------------------------------"
                                     "-----------------------------\n");
  CHKERRQ(ierr);

  MPI_Barrier(ps->comm->type);

  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    ierr = PSBUSIsGhosted(bus, &ghostbus);
    CHKERRQ(ierr);
    if (ghostbus)
      continue;

    Pd = Qd = 0.0;
    for (k = 0; k < bus->nload; k++) {
      ierr = PSBUSGetLoad(bus, k, &load);
      CHKERRQ(ierr);
      Pd += load->pl * ps->MVAbase;
      Qd += load->ql * ps->MVAbase;
    }
    ierr = PetscPrintf(
        ps->comm->type,
        "%-6d %7.2f %7.2f %7.3f %7.3f %12.2f %12.2f %12.2f %12.2f\n",
        bus->bus_i, Pd, Qd, bus->vm, bus->va * 180.0 / PETSC_PI, bus->mult_pmis,
        bus->mult_qmis, bus->pimb * ps->MVAbase, bus->qimb * ps->MVAbase);
    CHKERRQ(ierr);
  }
  ierr = PetscPrintf(ps->comm->type, "\n");
  CHKERRQ(ierr);
  MPI_Barrier(ps->comm->type);

  ierr = PetscPrintf(ps->comm->type,
                     "---------------------------------------------------------"
                     "-------------------------------\n");
  CHKERRQ(ierr);
  ierr = PetscPrintf(
      ps->comm->type, "%-10s %-8s %-10s %-8s %-7s %-8s %-8s %-8s\n", "From",
      "To", "Status", "Sft", "Stf", "Slim", "mult_Sf", "mult_St");
  CHKERRQ(ierr);
  ierr = PetscPrintf(ps->comm->type,
                     "---------------------------------------------------------"
                     "-------------------------------\n");
  CHKERRQ(ierr);
  MPI_Barrier(ps->comm->type);
  for (i = 0; i < ps->nline; i++) {
    line = &ps->line[i];
    ierr = PetscPrintf(
        ps->comm->type, "%-10d %-10d %-4d %8.2f %8.2f %8.2f %8.2f %8.2f\n",
        line->fbus, line->tbus, line->status, line->sf * ps->MVAbase,
        line->st * ps->MVAbase, (line->rateA > 1e5) ? 10000 : line->rateA,
        line->mult_sf, line->mult_st);
    CHKERRQ(ierr);
  }
  ierr = PetscPrintf(ps->comm->type, "\n");
  CHKERRQ(ierr);
  MPI_Barrier(ps->comm->type);

  ierr = PetscPrintf(ps->comm->type,
                     "---------------------------------------------------------"
                     "-------------------------------\n");
  CHKERRQ(ierr);
  ierr = PetscPrintf(
      ps->comm->type, "%-8s %-10s %-8s %-8s %-8s %-8s %-8s %-8s %-6s\n", "Gen",
      "Status", "Fuel", "Pg", "Qg", "Pmin", "Pmax", "Qmin", "Qmax");
  CHKERRQ(ierr);
  ierr = PetscPrintf(ps->comm->type,
                     "---------------------------------------------------------"
                     "-------------------------------\n");
  CHKERRQ(ierr);
  MPI_Barrier(ps->comm->type);
  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    ierr = PSBUSIsGhosted(bus, &ghostbus);
    CHKERRQ(ierr);
    if (ghostbus)
      continue;

    for (k = 0; k < bus->ngen; k++) {
      ierr = PSBUSGetGen(bus, k, &gen);
      CHKERRQ(ierr);
      ierr = PetscPrintf(
          PETSC_COMM_SELF,
          "%-10d %-4d %8s %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n", bus->bus_i,
          gen->status, PSGENFuelTypes[gen->genfuel_type], gen->pg * ps->MVAbase,
          gen->qg * ps->MVAbase, gen->pb * ps->MVAbase, gen->pt * ps->MVAbase,
          gen->qb * ps->MVAbase, gen->qt * ps->MVAbase);
      CHKERRQ(ierr);
    }
  }
  MPI_Barrier(ps->comm->type);

  PetscFunctionReturn(0);
}
