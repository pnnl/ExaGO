#include <private/psimpl.h>

// These globals are used for JSON format
static int tablevel = 0;

static char tabstring[64] = "";

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
      Pd += (load->pl - load->pl_loss) * ps->MVAbase;
      Qd += (load->ql - load->ql_loss) * ps->MVAbase;
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

static void DecreaseTabLevel()
{
  /* Decrease tab level */
  tablevel--;
  strcpy(tabstring,"");
  for(int i = 0; i < tablevel; i++) {
    strcat(tabstring,"\t");
  }
}

static void IncreaseTabLevel()
{
  /* Increase tab level */
  tablevel++;
  strcpy(tabstring,"");
  for(int i = 0; i < tablevel; i++) {
    strcat(tabstring,"\t");
  }
}

static void PrintJSONObjectBegin(FILE *fd,const char* name)
{
  if(name) {
    fprintf(fd,"%s\"%s\": {\n",tabstring,name);
  } else {
    fprintf(fd,"%s{\n",tabstring);
  }
  IncreaseTabLevel();
}

static void PrintJSONObjectEnd(FILE *fd,bool trail_comma)
{
  std::string str = trail_comma? ",":"";

  DecreaseTabLevel();
  
  fprintf(fd,"%s}%s\n",tabstring,str.c_str());    
}

static void PrintJSONArrayBegin(FILE *fd,const char* name)
{
  if(name) {
    fprintf(fd,"%s\"%s\": [\n",tabstring,name);
  } else {
    fprintf(fd,"%s[\n",tabstring);
  }
  IncreaseTabLevel();
}

static void PrintJSONArrayEnd(FILE *fd,bool trail_comma)
{
  std::string str = trail_comma? ",":"";

  DecreaseTabLevel();
  
  fprintf(fd,"%s]%s\n",tabstring,str.c_str());    
}

static void PrintJSONInt(FILE *fd,const char* key,int value,bool trail_comma) {
  std::string str = trail_comma? ",":"";
  fprintf(fd,"%s\"%s\": %d%s\n",tabstring,key,value,str.c_str());
}

static void PrintJSONDouble(FILE *fd,const char* key,double value,bool trail_comma) {
  std::string str = trail_comma? ",":"";
  fprintf(fd,"%s\"%s\": %lf%s\n",tabstring,key,value,str.c_str());
}

static void PrintJSONArrayDouble(FILE *fd,double value,bool trail_comma) {
  std::string str = trail_comma? ",":"";
  fprintf(fd,"%s%lf%s\n",tabstring,value,str.c_str());
}


static void PrintJSONString(FILE *fd,const char* key,char *value,bool trail_comma) {
  std::string str = trail_comma? ",":"";
  fprintf(fd,"%s\"%s\": \"%s\"%s\n",tabstring,key,value,str.c_str());
}

static void PrintJSONArray(FILE *fd,const char* name,int nvals,double *values,bool trail_comma)
{
  PrintJSONArrayBegin(fd,name);

  for(int i=0; i < nvals-1; i++) {
    PrintJSONArrayDouble(fd,values[i],true);
  }
  PrintJSONArrayDouble(fd,values[nvals-1],false);

  PrintJSONArrayEnd(fd,trail_comma);
}

/* 
   PSSaveSolution_JSON - Saves the system solution to file in JSON format

  Input Parameters:
+ ps - the PS object
- outfile  - Name of output file
*/
PetscErrorCode PSSaveSolution_JSON(PS ps, const char outfile[])
{
  PetscErrorCode ierr;
  FILE *fd;

  PSBUS bus;
  PSLOAD load;
  PSGEN gen;
  PSLINE line;
  PetscScalar Pd, Qd;
  PetscInt i, k;
  PetscScalar MVAbase = ps->MVAbase;
  char filename[PETSC_MAX_PATH_LEN];
  char *tok, *tok2;
  char sep[] = "/";
  char ext[] = ".json";
  char file1[PETSC_MAX_PATH_LEN];

  PetscFunctionBegin;

  strcpy(filename, outfile);
  /* Add .json extension to file name */
  ierr = PetscStrlcat(filename, ext, 256);
  CHKERRQ(ierr);

  fd = fopen(filename, "w");
  if (fd == NULL) {
    SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN,
             "Cannot open OPFLOW output file %s", outfile);
    CHKERRQ(ierr);
  }

  if(ps->gic_file_set) {
    ierr = PSReadGICData(ps);
    CHKERRQ(ierr);
  } else {
    PSSUBST subst;
    /* Assume 1 bus per substation */
    ierr = PetscCalloc1(ps->Nbus,&ps->substations);CHKERRQ(ierr);
    ps->nsubstations = 0;

    for(int i=0; i < ps->Nbus; i++) {
      subst = &ps->substations[ps->nsubstations];
      subst->num = i+1;
      snprintf(subst->name,64,"%d",ps->bus[i].bus_i);
      subst->nbus = 1;
      /* Linear distribution of lats and long from (30.0,-80.0) with 5 degrees deviation. This is completely random
	 baseless lat/long creation */
      subst->longlat[1]  = 30.0 + i/(PetscScalar)ps->Nbus*5.0;
      subst->longlat[0] = -80.0 + i/(PetscScalar)ps->Nbus*5.0;
      subst->bus[0] = &ps->bus[i];
      ps->nsubstations++;
    }
  }

  /* Begin JSON object */
  PrintJSONObjectBegin(fd,NULL);

    /* Print input file name */
    PrintJSONString(fd,"casefile",ps->net_file_name,true);
    /* Print gic file name */
    if(ps->gic_file_set) {
      PrintJSONString(fd,"gicfile",ps->gic_file_name,true);
    } else {
      PrintJSONString(fd,"gicfile","not given",true);
    }
  
    /* Print number of lines */
    PrintJSONInt(fd,"nbranch",ps->Nline,true);

    /* Print Number of gen */
    PrintJSONInt(fd,"ngen",ps->Ngen,true);

    /* Print Number of bus */
    PrintJSONInt(fd,"nbus",ps->Nbus,true);

    /* Print KV levels */
    PrintJSONArray(fd,"KVlevels",ps->nkvlevels,ps->kvlevels,true);
  
    /* Print output file name */
    PrintJSONString(fd,"casejsonfile",filename,true);

    PrintJSONObjectBegin(fd,"geojsondata");

      PrintJSONString(fd,"type","FeatureCollection",true);

      PrintJSONArrayBegin(fd,"features");

        for(i=0; i < ps->nsubstations; i++) {
        // Features
	  PrintJSONObjectBegin(fd,NULL);                  // Feature object start

	    PrintJSONString(fd,"type","Feature",true);
    
	    PrintJSONObjectBegin(fd,"geometry");                // Geometry object start

	      PrintJSONString(fd,"type","Point",true);

	      PrintJSONArray(fd,"coordinates",2,ps->substations[i].longlat,false);

	    PrintJSONObjectEnd(fd,true);                        // Geometry object end

	    PrintJSONObjectBegin(fd,"properties");              // Properties object start

	    PrintJSONObjectEnd(fd,false);                        // Properties object end
    
          PrintJSONObjectEnd(fd,true);                    // Feature object end
	}

	// Lines
	for(i=0; i < ps->Nline; i++) {
        // Features
	  PrintJSONObjectBegin(fd,NULL);                  // Feature object start

	    PrintJSONString(fd,"type","Feature",true);
    
	    PrintJSONObjectBegin(fd,"geometry");                // Geometry object start

	      PrintJSONString(fd,"type","LineString",true);

	    PrintJSONObjectEnd(fd,true);                        // Geometry object end

	    PrintJSONObjectBegin(fd,"properties");              // Properties object start

	    PrintJSONObjectEnd(fd,false);                        // Properties object end
	  if(i == ps->Nline-1) {
	    PrintJSONObjectEnd(fd,false);                    // Feature object end
	  } else {
	    PrintJSONObjectEnd(fd,true);
	  }
	}


      PrintJSONArrayEnd(fd,false); // features array end
  
    PrintJSONObjectEnd(fd,false); // geojsondata object end

  /* End of file */
  PrintJSONObjectEnd(fd,false);

  
  fclose(fd);
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
  } else if (format == JSON) {
    ierr = PSSaveSolution_JSON(ps,outfile);
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
  PetscScalar Pd, Qd, Pdloss, Qdloss;

  PetscFunctionBegin;

  ierr = PetscPrintf(ps->comm->type, "-----------------------------------------"
                                     "-----------------------------------------"
                                     "--------------------\n");
  ierr = PetscPrintf(
      ps->comm->type,
      "%-10s %-7s %-6s %-7s %-6s %-7s %-7s %-14s %-14s %-14s %-14s\n", "Bus",
      "Pd", "Pdloss", "Qd", "Qdloss", "Vm", "Va", "mult_Pmis", "mult_Qmis",
      "Pslack", "Qslack");
  CHKERRQ(ierr);
  ierr = PetscPrintf(ps->comm->type, "-----------------------------------------"
                                     "-----------------------------------------"
                                     "--------------------\n");
  CHKERRQ(ierr);

  MPI_Barrier(ps->comm->type);

  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];

    ierr = PSBUSIsGhosted(bus, &ghostbus);
    CHKERRQ(ierr);
    if (ghostbus)
      continue;

    Pd = Qd = Pdloss = Qdloss = 0.0;
    for (k = 0; k < bus->nload; k++) {
      ierr = PSBUSGetLoad(bus, k, &load);
      CHKERRQ(ierr);
      Pd += load->pl * ps->MVAbase;
      Qd += load->ql * ps->MVAbase;
      Pdloss += load->pl_loss * ps->MVAbase;
      Qdloss += load->ql_loss * ps->MVAbase;
    }
    ierr =
        PetscPrintf(ps->comm->type,
                    "%-6d %7.2f %7.2f %7.2f %7.2f %7.3f %7.3f %12.2f %12.2f "
                    "%12.2f %12.2f\n",
                    bus->bus_i, Pd, Pdloss, Qd, Qdloss, bus->vm,
                    bus->va * 180.0 / PETSC_PI, bus->mult_pmis, bus->mult_qmis,
                    bus->pimb * ps->MVAbase, bus->qimb * ps->MVAbase);
    CHKERRQ(ierr);
  }
  ierr = PetscPrintf(ps->comm->type, "\n");
  CHKERRQ(ierr);
  MPI_Barrier(ps->comm->type);

  ierr = PetscPrintf(ps->comm->type, "-----------------------------------------"
                                     "-----------------------------------------"
                                     "--------------------\n");
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
