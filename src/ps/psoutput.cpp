#include <libgen.h>
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
  char dir[PETSC_MAX_PATH_LEN];
  bool gen_fuel_defined;

  PetscFunctionBegin;

  strncpy(filename, outfile, PETSC_MAX_PATH_LEN);
  strncpy(file1, basename(filename), PETSC_MAX_PATH_LEN);
  strncpy(filename, outfile, PETSC_MAX_PATH_LEN);
  strncpy(dir, dirname(filename), PETSC_MAX_PATH_LEN);

  /* Check if file has .m extension (find last occurance) */

  tok = &file1[0];
  tok2 = NULL;
  while (tok != NULL) {
    tok = strstr(tok, ext);
    if (tok != NULL) {
      tok2 = tok;
      tok++;
    }
  }

  /* remove extension from file1, if found */

  if (tok2 != NULL) {
    *tok2 = (char)0;
  }

  /* use this as a function name */

  strncpy(fcn_name, file1, 100);

  /* re-assemble file name */

  strncpy(filename, dir, PETSC_MAX_PATH_LEN);
  ierr = PetscStrlcat(filename, sep, PETSC_MAX_PATH_LEN);
  ierr = PetscStrlcat(filename, fcn_name, PETSC_MAX_PATH_LEN);
  ierr = PetscStrlcat(filename, ext, PETSC_MAX_PATH_LEN);

  /* save file */

  fd = fopen(filename, "w");
  if (fd == NULL) {
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN,
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

  /* generator fuel info */

  /* only output fuel type if it is defined */

  gen_fuel_defined = false;
  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    for (k = 0; k < bus->ngen; k++) {
      ierr = PSBUSGetGen(bus, k, &gen);
      CHKERRQ(ierr);
      PetscInt gfuel = gen->genfuel_type;
      if (gfuel != GENFUEL_UNDEFINED) {
        gen_fuel_defined = true;
        break;
      }
    }
  }

  if (gen_fuel_defined) {
    fprintf(fd, "\n%%%% generator fuel type\n");
    fprintf(fd, "mpc.genfuel = {\n");
    for (i = 0; i < ps->nbus; i++) {
      bus = &ps->bus[i];
      for (k = 0; k < bus->ngen; k++) {
        ierr = PSBUSGetGen(bus, k, &gen);
        CHKERRQ(ierr);
        PetscInt gfuel = gen->genfuel_type;
        switch (gfuel) {
        case GENFUEL_COAL:
          fprintf(fd, "\t\"coal\"");
          break;
        case GENFUEL_WIND:
          fprintf(fd, "\t\"wind\"");
          break;
        case GENFUEL_SOLAR:
          fprintf(fd, "\t\"solar\"");
          break;
        case GENFUEL_NG:
          fprintf(fd, "\t\"ng\"");
          break;
        case GENFUEL_NUCLEAR:
          fprintf(fd, "\t\"nuclear\"");
          break;
        case GENFUEL_HYDRO:
          fprintf(fd, "\t\"hydro\"");
          break;
        case GENFUEL_UNDEFINED:
        default:
	  fprintf(fd, "\t\"other\"");
          break;
        }
        fprintf(fd, ";\n");
      }
    }
    fprintf(fd, "\n};\n");
  }

  /* Load cost data */
  fprintf(fd, "\n%%%% load cost data\n");
  fprintf(fd, "%% %%maxallowedloadshed loadshedcost\n");
  fprintf(fd, "mpc.loadcost = [\n");
  for (i = 0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    for (k = 0; k < bus->nload; k++) {
      PSLOAD load;
      ierr = PSBUSGetLoad(bus, k, &load);
      CHKERRQ(ierr);
      fprintf(fd, "%10.5g %10.5g; \n", load->loss_frac * 100.0,
              load->loss_cost);
    }
  }
  fprintf(fd, "];\n");

  /* Solution summary info */
  fprintf(fd, "\n%%%% summary data\n");
  fprintf(fd, "%ssummary_stats = [\n", prefix);
  fprintf(fd,
          "%%\tNbus\tNgen\tNgenON\tNline\tNlineON\tNload\tGenPCap\tGenTotalP"
          "\tGenTota"
          "lQ\tGenPCapON\tLoadTotP\tLoadTotalQ\tLoadShedP\tLoadShedQ\n");
  fprintf(fd,
          "\t%d\t%d\t%d\t%d\t%d\t%d\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%."
          "9g\t%.9g\n",
          ps->sys_info.Nbus, ps->sys_info.Ngen, ps->sys_info.NgenON,
          ps->sys_info.Nline, ps->sys_info.NlineON, ps->sys_info.Nload,
          ps->sys_info.total_pgencap, ps->sys_info.total_genON[0],
          ps->sys_info.total_genON[1], ps->sys_info.total_pgencapON,
          ps->sys_info.total_load[0], ps->sys_info.total_load[1],
          ps->sys_info.total_loadshed[0], ps->sys_info.total_loadshed[1]);
  fprintf(fd, "];\n");

  fprintf(fd, "\n%%%% solve time\n");
  fprintf(fd, "%ssolve_time = %.5g;\n", prefix, ps->solve_real_time);
  fprintf(fd, "%ssolve_cpu_time = %.5g;\n", prefix, ps->solve_cpu_time);

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
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN,
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

static void DecreaseTabLevel() {
  /* Decrease tab level */
  tablevel--;
  strcpy(tabstring, "");
  for (int i = 0; i < tablevel; i++) {
    strcat(tabstring, "\t");
  }
}

static void IncreaseTabLevel() {
  /* Increase tab level */
  tablevel++;
  strcpy(tabstring, "");
  for (int i = 0; i < tablevel; i++) {
    strcat(tabstring, "\t");
  }
}

static void PrintJSONObjectBegin(FILE *fd, const char *name) {
  if (name) {
    fprintf(fd, "%s\"%s\": {\n", tabstring, name);
  } else {
    fprintf(fd, "%s{\n", tabstring);
  }
  IncreaseTabLevel();
}

static void PrintJSONObjectEnd(FILE *fd, bool trail_comma) {
  std::string str = trail_comma ? "," : "";

  DecreaseTabLevel();

  fprintf(fd, "%s}%s\n", tabstring, str.c_str());
}

static void PrintJSONArrayBegin(FILE *fd, const char *name) {
  if (name) {
    fprintf(fd, "%s\"%s\": [\n", tabstring, name);
  } else {
    fprintf(fd, "%s[\n", tabstring);
  }
  IncreaseTabLevel();
}

static void PrintJSONArrayEnd(FILE *fd, bool trail_comma) {
  std::string str = trail_comma ? "," : "";

  DecreaseTabLevel();

  fprintf(fd, "%s]%s\n", tabstring, str.c_str());
}

static void PrintJSONInt(FILE *fd, const char *key, int value,
                         bool trail_comma) {
  std::string str = trail_comma ? "," : "";
  fprintf(fd, "%s\"%s\": %d%s\n", tabstring, key, value, str.c_str());
}

static void PrintJSONDouble(FILE *fd, const char *key, double value,
                            bool trail_comma) {
  std::string str = trail_comma ? "," : "";
  fprintf(fd, "%s\"%s\": %lf%s\n", tabstring, key, value, str.c_str());
}

static void PrintJSONArrayDouble(FILE *fd, double value, bool trail_comma) {
  std::string str = trail_comma ? "," : "";
  fprintf(fd, "%s%lf%s\n", tabstring, value, str.c_str());
}

static void PrintJSONString(FILE *fd, const char *key, const char *value,
                            bool trail_comma) {
  std::string str = trail_comma ? "," : "";
  fprintf(fd, "%s\"%s\": \"%s\"%s\n", tabstring, key, value, str.c_str());
}

static void PrintJSONArray(FILE *fd, const char *name, int nvals,
                           double *values, bool trail_comma) {
  PrintJSONArrayBegin(fd, name);

  for (int i = 0; i < nvals - 1; i++) {
    PrintJSONArrayDouble(fd, values[i], true);
  }
  PrintJSONArrayDouble(fd, values[nvals - 1], false);

  PrintJSONArrayEnd(fd, trail_comma);
}

static void PrintGenData(FILE *fd, PSBUS bus, bool trail_comma,
                         PetscScalar MVAbase) {
  PSGEN gen;
  bool istrail;

  PrintJSONArrayBegin(fd, "gen");

  for (int i = 0; i < bus->ngen; i++) {
    PSBUSGetGen(bus, i, &gen);
    if (i == bus->ngen - 1)
      istrail = false;
    else
      istrail = true;

    PrintJSONObjectBegin(fd, NULL);

    // Generator bus number
    PrintJSONInt(fd, "GEN_BUS", bus->bus_i, true);

    // Generator fuel
    PrintJSONString(fd, "GEN_FUEL", PSGENFuelTypes[gen->genfuel_type], true);

    // Generator Pg
    PrintJSONDouble(fd, "PG", gen->pg * MVAbase, true);

    // Generator Qg
    PrintJSONDouble(fd, "QG", gen->qg * MVAbase, true);

    // Generator status
    PrintJSONInt(fd, "GEN_STATUS", gen->status, true);

    // Generator Pmax
    PrintJSONDouble(fd, "PMAX", gen->pt * MVAbase, true);

    // Generator Pmin
    PrintJSONDouble(fd, "PMIN", gen->pb * MVAbase, true);

    // Generator QMAX
    PrintJSONDouble(fd, "QMAX", gen->qt * MVAbase, true);

    // Generator QMIN
    PrintJSONDouble(fd, "QMIN", gen->qb * MVAbase, false);

    PrintJSONObjectEnd(fd, istrail);
  }

  PrintJSONArrayEnd(fd, trail_comma);
}

static void PrintLineData(FILE *fd, PSLINE line, bool trail_comma,
                          PetscScalar MVAbase) {
  // unused
  (void)trail_comma;
  std::string line_name;
  // Elementtype
  PrintJSONString(fd, "elementtype", "Branch", true);

  // Name
  line_name = std::string(line->subst_from->name) + " -- " +
              std::string(line->subst_to->name);

  PrintJSONString(fd, "NAME", line_name.c_str(), true);

  // From bus
  PrintJSONInt(fd, "F_BUS", line->fbus, true);

  // To bus
  PrintJSONInt(fd, "T_BUS", line->tbus, true);

  // Status
  PrintJSONInt(fd, "BR_STATUS", line->status, true);

  // KV level
  PrintJSONDouble(fd, "KV", line->kvlevel, true);

  // Rate A
  PrintJSONDouble(fd, "RATE_A", (line->rateA > 1e5) ? 10000 : line->rateA,
                  true);

  // PF,QF, PT, QT
  PrintJSONDouble(fd, "PF", line->pf * MVAbase, true);
  PrintJSONDouble(fd, "QF", line->qf * MVAbase, true);
  PrintJSONDouble(fd, "PT", line->pt * MVAbase, true);
  PrintJSONDouble(fd, "QT", line->qt * MVAbase, false);
}

static void PrintBusData(FILE *fd, PSSUBST subst, bool trail_comma,
                         PetscScalar MVAbase) {
  // unused
  (void)trail_comma;
  PSBUS bus;
  PSLOAD load;
  bool istrail;
  char bus_name[64];
  double Pd = 0.0, Qd = 0.0, Pdloss = 0.0, Qdloss = 0.0;
  double Vm_avg = 0.0;

  PrintJSONArrayBegin(fd, "bus");

  for (int i = 0; i < subst->nbus; i++) {
    bus = subst->bus[i];
    PrintJSONObjectBegin(fd, NULL);

    // elementtype
    PrintJSONString(fd, "elementtype", "bus", true);

    // bus number
    PrintJSONInt(fd, "BUS_I", bus->bus_i, true);

    // VA
    PrintJSONDouble(fd, "VA", bus->va * 180.0 / PETSC_PI, true);

    // VM
    PrintJSONDouble(fd, "VM", bus->vm, true);
    Vm_avg += bus->vm;

    // bus name
    strcpy(bus_name, subst->name);
    snprintf(bus_name, 64, " %d", i + 1);
    PrintJSONString(fd, "BUS_NAME", bus_name, true);

    // Vmin
    PrintJSONDouble(fd, "VMIN", bus->Vmin, true);

    // Vmax
    PrintJSONDouble(fd, "VMAX", bus->Vmax, true);

    // Base KV
    PrintJSONDouble(fd, "BASE_KV", bus->basekV, true);

    // PD
    if (bus->nload) {
      PSBUSGetLoad(bus, 0, &load);

      Pd = load->pl * MVAbase;
      Qd = load->ql * MVAbase;
      Pdloss = load->pl_loss * MVAbase;
      Qdloss = load->ql_loss * MVAbase;
    }

    PrintJSONDouble(fd, "PD", Pd, true);
    PrintJSONDouble(fd, "QD", Qd, true);
    PrintJSONDouble(fd, "PDloss", Pdloss, true);
    PrintJSONDouble(fd, "QDloss", Qdloss, true);

    // Lagrange multipliers
    PrintJSONDouble(fd, "LAM_P", bus->mult_pmis, true);
    PrintJSONDouble(fd, "LAM_Q", bus->mult_qmis, true);

    if (i == subst->nbus - 1)
      istrail = false;
    else
      istrail = true;

    // ngen
    PrintJSONInt(fd, "ngen", bus->ngen, true);

    PrintGenData(fd, bus, false, MVAbase);

    PrintJSONObjectEnd(fd, istrail);
  }
  Vm_avg /= subst->nbus;

  PrintJSONArrayEnd(fd, true);

  // Print Vm_avg in the substation object
  PrintJSONDouble(fd, "Vm", Vm_avg, false);
}
/*
   PSSaveSolution_JSON - Saves the system solution to file in JSON format

  Input Parameters:
+ ps - the PS object
- outfile  - Name of output file
*/
PetscErrorCode PSSaveSolution_JSON(PS ps, const char outfile[]) {
  PetscErrorCode ierr;
  FILE *fd;

  PSBUS bus;
  PSLINE line;

  PetscInt i;

  PetscScalar MVAbase = ps->MVAbase;
  char filename[PETSC_MAX_PATH_LEN];
  char ext[] = ".json";

  PetscFunctionBegin;

  strcpy(filename, outfile);
  /* Add .json extension to file name */
  ierr = PetscStrlcat(filename, ext, 256);
  CHKERRQ(ierr);

  fd = fopen(filename, "w");
  if (fd == NULL) {
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN,
            "Cannot open OPFLOW output file %s", outfile);
    CHKERRQ(ierr);
  }

  if (ps->gic_file_set) {
    ierr = PSReadGICData(ps);
    CHKERRQ(ierr);
  } else {
    PSSUBST subst;
    /* Assume 1 bus per substation */
    ierr = PetscCalloc1(ps->Nbus, &ps->substations);
    CHKERRQ(ierr);
    ps->nsubstations = 0;

    for (int i = 0; i < ps->Nbus; i++) {
      PetscInt nconnlines;
      const PSLINE *connlines;
      PSLINE branch;

      subst = &ps->substations[ps->nsubstations];
      subst->num = i + 1;
      snprintf(subst->name, 64, "%d", ps->bus[i].bus_i);
      subst->nbus = 1;
      subst->nkvlevels = 1;
      /* Circular distribution of lats and long from some location with .5
       * degrees deviation. This is completely random baseless lat/long creation
       */
      subst->longlat[1] =
          35.481918 +
          PetscCosScalar(i / (PetscScalar)ps->Nbus * 2.0 * PETSC_PI) * 0.5;
      subst->longlat[0] =
          -97.508469 +
          PetscSinScalar(i / (PetscScalar)ps->Nbus * 2.0 * PETSC_PI) * 0.5;
      bus = &ps->bus[i];
      subst->bus[0] = bus;
      subst->kvlevels[0] = bus->basekV;

      ps->nsubstations++;

      /* Get the connected lines to the bus */
      PSBUSGetSupportingLines(bus, &nconnlines, &connlines);
      for (int k = 0; k < nconnlines; k++) {
        branch = connlines[k];

        const PSBUS *connbuses;
        PSBUS busf; //, bust;

        /* Get the connected buses to this line */
        PSLINEGetConnectedBuses(branch, &connbuses);
        busf = connbuses[0];
        // bust = connbuses[1];

        if (bus == busf) { /* From bus */
          branch->subst_from = subst;
        } else {
          branch->subst_to = subst;
        }
      }
    }
  }

  /* Begin JSON object */
  PrintJSONObjectBegin(fd, NULL);

  /* Print input file name */
  PrintJSONString(fd, "casefile", ps->net_file_name, true);
  /* Print gic file name */
  if (ps->gic_file_set) {
    PrintJSONString(fd, "gicfile", ps->gic_file_name, true);
  } else {
    PrintJSONString(fd, "gicfile", "not given", true);
  }

  /* Print number of lines */
  PrintJSONInt(fd, "nbranch", ps->Nline, true);

  /* Print Number of gen */
  PrintJSONInt(fd, "ngen", ps->Ngen, true);

  /* Print Number of bus */
  PrintJSONInt(fd, "nbus", ps->Nbus, true);

  /* Print KV levels */
  PrintJSONArray(fd, "KVlevels", ps->nkvlevels, ps->kvlevels, true);

  /* Print output file name */
  PrintJSONString(fd, "casejsonfile", filename, true);

  PrintJSONObjectBegin(fd, "geojsondata");

  PrintJSONString(fd, "type", "FeatureCollection", true);

  PrintJSONArrayBegin(fd, "features");

  for (i = 0; i < ps->nsubstations; i++) {
    // Features
    PrintJSONObjectBegin(fd, NULL); // Feature object start

    PrintJSONString(fd, "type", "Feature", true);

    PrintJSONObjectBegin(fd, "geometry"); // Geometry object start

    PrintJSONString(fd, "type", "Point", true);

    PrintJSONArray(fd, "coordinates", 2, ps->substations[i].longlat, false);

    PrintJSONObjectEnd(fd, true); // Geometry object end

    PrintJSONObjectBegin(fd, "properties"); // Properties object start

    // Element type
    PrintJSONString(fd, "elementtype", "SUBSTATION", true);

    // Name
    PrintJSONString(fd, "NAME", ps->substations[i].name, true);

    // Number of buses
    PrintJSONInt(fd, "nbus", ps->substations[i].nbus, true);

    // Print KV levels
    PrintJSONArray(fd, "KVlevels", ps->substations[i].nkvlevels,
                   ps->substations[i].kvlevels, true);

    // Print the substation data
    PrintBusData(fd, &ps->substations[i], false, ps->MVAbase);

    PrintJSONObjectEnd(fd, false); // Properties object end

    PrintJSONObjectEnd(fd, true); // Feature object end
  }

  // Lines
  for (i = 0; i < ps->Nline; i++) {
    line = &ps->line[i];
    // Features
    PrintJSONObjectBegin(fd, NULL); // Feature object start

    PrintJSONString(fd, "type", "Feature", true);

    PrintJSONObjectBegin(fd, "geometry"); // Geometry object start

    PrintJSONString(fd, "type", "LineString", true);

    // Coordinates
    PrintJSONArrayBegin(fd, "coordinates");

    // Print from substation coordinates
    PrintJSONArray(fd, NULL, 2, line->subst_from->longlat, true);
    // Print to substation coordinates
    PrintJSONArray(fd, NULL, 2, line->subst_to->longlat, false);

    PrintJSONArrayEnd(fd, false);

    PrintJSONObjectEnd(fd, true); // Geometry object end

    PrintJSONObjectBegin(fd, "properties"); // Properties object start
                                            // Print line data
    PrintLineData(fd, line, false, MVAbase);

    PrintJSONObjectEnd(fd, false); // Properties object end
    if (i == ps->Nline - 1) {
      PrintJSONObjectEnd(fd, false); // Feature object end
    } else {
      PrintJSONObjectEnd(fd, true);
    }
  }

  PrintJSONArrayEnd(fd, false); // features array end

  PrintJSONObjectEnd(fd, true); // geojsondata object end

  PrintJSONObjectBegin(fd, "summary"); // System summary object start

  PrintJSONInt(fd, "NBUS", ps->sys_info.Nbus, true);
  PrintJSONInt(fd, "NGEN", ps->sys_info.Ngen, true);
  PrintJSONInt(fd, "NGENON", ps->sys_info.NgenON, true);
  PrintJSONInt(fd, "NLINE", ps->sys_info.Nline, true);
  PrintJSONInt(fd, "NLINEON", ps->sys_info.NlineON, true);
  PrintJSONInt(fd, "NLOAD", ps->sys_info.Nload, true);

  PrintJSONDouble(fd, "GENCAP", ps->sys_info.total_pgencap, true);
  PrintJSONArray(fd, "GENON", 2, &ps->sys_info.total_genON[0], true);
  PrintJSONDouble(fd, "GENCAPON", ps->sys_info.total_pgencapON, true);
  PrintJSONArray(fd, "LOAD", 2, &ps->sys_info.total_load[0], true);
  PrintJSONArray(fd, "LOADSHED", 2, &ps->sys_info.total_loadshed[0], true);

  PrintJSONDouble(fd, "SolveRealTime", ps->solve_real_time, true);
  PrintJSONDouble(fd, "SolveCPUTime", ps->solve_cpu_time, true);

  PrintJSONObjectEnd(fd, false); // System summary object start

  /* End of file */
  PrintJSONObjectEnd(fd, false);

  fclose(fd);
  PetscFunctionReturn(0);
}

/*
   PSSaveSolution_MINIMAL - Saves minimal solution information to file
  Input Parameters:
+ ps - the PS object
- outfile  - Name of output file
*/
PetscErrorCode PSSaveSolution_MINIMAL(PS ps, const char outfile[]) {
  PetscErrorCode ierr;
  FILE *fd;
  char filename[PETSC_MAX_PATH_LEN];
  char ext[] = ".txt";

  PetscFunctionBegin;

  strcpy(filename, outfile);
  /* Add extension to file name */
  ierr = PetscStrlcat(filename, ext, 256);
  CHKERRQ(ierr);

  fd = fopen(filename, "w");
  if (fd == NULL) {
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN,
            "Cannot open OPFLOW output file %s", outfile);
    CHKERRQ(ierr);
  }

  fprintf(fd, "Converged: %s\n", ps->opflow_converged ? "Yes" : "No");
  fprintf(fd, "Objective: %g\n", ps->opflowobj);

  fprintf(fd, "Summary: \n");
  fprintf(fd, "\tBuses: %d\n", ps->sys_info.Nbus);
  fprintf(fd, "\tGenerators: %d\n", ps->sys_info.Ngen);
  fprintf(fd, "\tGenerators ON: %d\n", ps->sys_info.NgenON);
  fprintf(fd, "\tLines: %d\n", ps->sys_info.Nline);
  fprintf(fd, "\tLines ON: %d\n", ps->sys_info.NlineON);
  fprintf(fd, "\tLoads: %d\n", ps->sys_info.Nload);
  fprintf(fd, "\tGeneration Capacity: %.9g\n", ps->sys_info.total_pgencap);
  fprintf(fd, "\tTotal Generation Online P, Q: %9g, %9g\n",
          ps->sys_info.total_genON[0], ps->sys_info.total_genON[1]);
  fprintf(fd, "\tGeneration Capacity ON: %.9g\n", ps->sys_info.total_pgencapON);
  fprintf(fd, "\tTotal Load P, Q: %9g, %9g\n", ps->sys_info.total_load[0],
          ps->sys_info.total_load[1]);
  fprintf(fd, "\tTotal Load Shed P, Q: %9g, %9g\n",
          ps->sys_info.total_loadshed[0], ps->sys_info.total_loadshed[1]);
  fprintf(fd, "\tSolve Time: %5g\n", ps->solve_real_time);
  fprintf(fd, "\tSolve CPU Time: %5g\n", ps->solve_cpu_time);

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
    ierr = PSSaveSolution_JSON(ps, outfile);
    CHKERRQ(ierr);
  } else if (format == MINIMAL) {
    ierr = PSSaveSolution_MINIMAL(ps, outfile);
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
