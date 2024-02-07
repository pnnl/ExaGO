#include <private/psimpl.h>

/*
  PSReadPSSERawData - Reads the PSSE raw data file and populates the PS object

  Input Parameter
+ ps      - The power system object ps
- netfile - Name of the power system network data file in PSSE raw data format.

*/
PetscErrorCode PSReadPSSERawData(PS ps, const char netfile[]) {
  FILE *fp;
  PetscErrorCode ierr;
  char *out;

  PetscFunctionBegin;

  char line[MAXLINE];
  int dataformatflag = -1;
  /*
  0 = BUS DATA, 1 = LOAD DATA, 2 = SHUNT DATA, 3 = GENERATOR DATA, 4 = BRANCH
  DATA, 5 = TRANSFORMER DATA, 6 = AREA DATA 7 = TWO-TERMINAL DC DATA, 8 = VSC DC
  LINE DATA, 9 = IMPEDANCE CORRECTION DATA, 10 = MULTI-TERMINAL DC DATA, 11 =
  MULTI-SECTION LINE DATA 12 = ZONE DATA, 13 = TRANSFER DATA, 14 = OWNER DATA,
  15 = FACTS DEVICE DATA, 16 = SWITCHED SHUNT DATA, 17 = GNE DATA
  */
  // Network Size variables
  PetscInt Nbus = 0;
  PetscInt Nload = 0;
  PetscInt Ngenerator = 0;
  PetscInt Nline = 0;
  PetscInt Ntransformer = 0;
  // Case Identification Data
  PetscInt IC = 0;
  PetscScalar SBASE = 100.0;
  PetscInt REV = 33;
  PetscInt XFRRAT = 0;
  PetscInt NXFRAT = 0;
  PetscScalar BASFRQ = 0.0;
  // Structure Pointers in PS
  PSBUS Bus;
  PSLOAD Load;
  PSGEN Gen;
  PSLINE Branch;
  // Unused Transformer data
  PetscInt K, CW, CZ, CM, NMETR, COD1, CONT1, NTP1, TAB1;
  PetscScalar MAG1, MAG2, NOMV1, RMA1, RMI1, VMA1, VMI1, CR1, CX1, CNXA1,
      WINDV2, NOMV2;
  char transname[20];
  // Temp Variables
  PetscInt loadi = 0, geni = 0, bri = 0, busi = 0, i = 0;
  PetscInt MET = 1;
  PetscInt internalindex = 0;
  PetscScalar R, X, Bc, B, G, Zm, tap, shift, tap2, tapr, tapi;
  PetscInt tempbusi = 0, maxbusi = -1;
  PetscInt shuntbus;
  char shuntid[10];
  PetscInt shuntstatus;
  PetscScalar gshunt, bshunt;

  if (ps->comm->type != PETSC_COMM_SELF && ps->comm->rank != 0) {
    ps->Nline = ps->Nbus = ps->Ngen = ps->Nload = 0;
    PetscFunctionReturn(0);
  }

  fp = fopen(netfile, "r");
  /* Check for valid file */
  if (fp == NULL) {
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Cannot open file %s",
            netfile);
    CHKERRQ(ierr);
  }

  /* Copy the network file name */
  ierr = PetscStrcpy(ps->net_file_name, netfile);
  CHKERRQ(ierr);

  // Read Case Identification Data
  out = fgets(line, MAXLINE, fp);
  if (out != NULL) {
    sscanf(line, "%d, %lf, %d, %d, %d, %lf", &IC, &SBASE, &REV, &XFRRAT,
           &NXFRAT, &BASFRQ);
    out = fgets(line, MAXLINE, fp);
    if (out == NULL)
      PetscFunctionReturn(0); // for commas
    out = fgets(line, MAXLINE, fp);
    if (out == NULL)
      PetscFunctionReturn(0); // for commas
    dataformatflag++;
  } else {
    PetscFunctionReturn(0);
  }
  /* Setting ps->sbase to 100 */
  ps->MVAbase = SBASE;

  while ((out = fgets(line, MAXLINE, fp)) != NULL) {
    if (strstr(line, "0 /") != NULL) {
      dataformatflag++;
      continue;
    }

    switch (dataformatflag) {
    case 0:
      sscanf(line, "%d", &tempbusi);
      if (tempbusi > maxbusi)
        maxbusi = tempbusi;
      Nbus++;
      break;
    case 1:
      Nload++;
      break;
    case 2:
      // Nshunt++;
      break;
    case 3:
      Ngenerator++;
      break;
    case 4:
      Nline++;
      break;
    case 5:
      out = fgets(line, MAXLINE, fp);
      out = fgets(line, MAXLINE, fp);
      out = fgets(line, MAXLINE, fp);
      Ntransformer++;
      break;
    default:
      // Todo: impliment parser for other data
      break;
    }
  }
  fclose(fp);

  ps->Nbus = ps->nbus = Nbus;
  ps->Nload = ps->nload = Nload;
  ps->Ngen = ps->ngen = Ngenerator;
  ps->Nline = ps->nline = Nline + Ntransformer;
  ps->NgenON = 0;
  ps->NlineON = 0;

#if defined DEBUGPS
  ierr = PetscPrintf(
      PETSC_COMM_SELF,
      "System summary : Nbus = %d, Nload = %d, Ngenerator = %d, Nbranch = %d\n",
      ps->Nbus, ps->Ngen, ps->Nload, ps->Nline);
  CHKERRQ(ierr);
#endif
  ierr = PetscCalloc1(ps->Nbus, &ps->bus);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(ps->Ngen, &ps->gen);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(ps->Nload, &ps->load);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(ps->Nline, &ps->line);
  CHKERRQ(ierr);
  Bus = ps->bus;
  Gen = ps->gen;
  Load = ps->load;
  Branch = ps->line;

  // Set initial and default data for bus
  for (i = 0; i < ps->Nbus; i++) {
    ps->bus[i].ngen = ps->bus[i].nload = ps->bus[i].ngenON = 0;
    ps->bus[i].qrange = ps->bus[i].qmintot = ps->bus[i].Pgtot =
        ps->bus[i].MVAbasetot = 0.0;
    ps->bus[i].Vmax = 1.1;
    ps->bus[i].Vmin = 0.9;
  }
  /* Allocate external to internal bus number mapping array */

  PetscInt *busext2intmap;
  ps->maxbusnum = maxbusi;
  ierr = PetscCalloc1(ps->maxbusnum + 1, &ps->busext2intmap);
  CHKERRQ(ierr);
  busext2intmap = ps->busext2intmap;
  for (i = 0; i < ps->maxbusnum + 1; i++)
    busext2intmap[i] = -1;

  dataformatflag = -1;
  fp = fopen(netfile, "r");
  if ((out = fgets(line, MAXLINE, fp)) != NULL) { // for case identification
    if ((out = fgets(line, MAXLINE, fp)) == NULL)
      PetscFunctionReturn(0); // for commas
    if ((out = fgets(line, MAXLINE, fp)) == NULL)
      PetscFunctionReturn(0); // for commas
    dataformatflag++;
  } else {
    PetscFunctionReturn(0);
  }
  while ((out = fgets(line, MAXLINE, fp)) != NULL) {
    if (strstr(line, "0 /") != NULL) {
      dataformatflag++;
      continue;
    }
    switch (dataformatflag) {
    case 0:
      sscanf(line, "%d, '%[^\t\']', %lf, %d, %d, %d, %d, %lf, %lf",
             &Bus[busi].bus_i, Bus[busi].name, &Bus[busi].basekV,
             &Bus[busi].ide, &Bus[busi].area, &Bus[busi].zone, &Bus[busi].owner,
             &Bus[busi].vm, &Bus[busi].va);
#if defined DEBUGPS
      if (busi < 2)
        ierr = PetscPrintf(
            PETSC_COMM_SELF,
            "BUSData[%d] : %d, '%s', %lf, %d, %d, %d, %d, %lf, %lf\n", busi,
            Bus[busi].bus_i, Bus[busi].name, Bus[busi].basekV, Bus[busi].ide,
            Bus[busi].area, Bus[busi].zone, Bus[busi].owner, Bus[busi].vm,
            Bus[busi].va);
      CHKERRQ(ierr);
#endif
      // Convert angle to radians
      Bus[busi].va *= PETSC_PI / 180.0;

      if (Bus[busi].ide == REF_BUS)
        ps->Nref++;

      busext2intmap[Bus[busi].bus_i] = busi;
      Bus[busi].internal_i = busi;
      Bus[busi].nload = 0;
      Bus[busi].ngen = 0;
      Bus[busi].ngenON = 0;
      Bus[busi].nshunt = 0;
      Bus[busi].nvhi = 1.1;
      Bus[busi].nvlo = 0.9;
      Bus[busi].evhi = 1.1;
      Bus[busi].evlo = 0.9;
      Bus[busi].Vmin = 1.1;
      Bus[busi].Vmax = 0.9;
      Bus[busi].gl = 0;
      Bus[busi].bl = 0;
      busi++;
      break;
    case 1:
      sscanf(line,
             "%d, '%[^\t\']', %d, %d, %d, %lf, %lf, %lf, %lf, %lf, %lf, %d",
             &Load[loadi].bus_i, Load[loadi].id, &Load[loadi].status,
             &Load[loadi].area, &Load[loadi].zone, &Load[loadi].pl,
             &Load[loadi].ql, &Load[loadi].ip, &Load[loadi].iq, &Load[loadi].yp,
             &Load[loadi].yq, &Load[loadi].owner);
#if defined DEBUGPS
      if (loadi < 2)
        ierr = PetscPrintf(PETSC_COMM_SELF,
                           "LOADData[%d] : %d, '%s', %d, %d, %d, %lf, %lf, "
                           "%lf, %lf, %lf, %lf, %d\n",
                           loadi, Load[loadi].bus_i, Load[loadi].id,
                           Load[loadi].status, Load[loadi].area,
                           Load[loadi].zone, Load[loadi].pl, Load[loadi].ql,
                           Load[loadi].ip, Load[loadi].iq, Load[loadi].yp,
                           Load[loadi].yq, Load[loadi].owner);
      CHKERRQ(ierr);
#endif
      Load[loadi].scale = 1;
      Load[loadi].intrpt = 0;
      Load[loadi].pl /= ps->MVAbase;
      Load[loadi].ql /= ps->MVAbase;
      Load[loadi].ip /= ps->MVAbase;
      Load[loadi].iq /= ps->MVAbase;
      Load[loadi].yp /= ps->MVAbase;
      Load[loadi].yq /= ps->MVAbase;
      internalindex = busext2intmap[Load[loadi].bus_i];
      Load[loadi].internal_i = internalindex;
      Bus[internalindex].lidx[Bus[internalindex].nload] = loadi;
      Bus[internalindex].nload++;
      loadi++;
      break;
    case 2:
      sscanf(line, "%d, '%[^\t\']', %d, %lf, %lf", &shuntbus, shuntid,
             &shuntstatus, &gshunt, &bshunt);
      if (shuntstatus) {
        internalindex = busext2intmap[shuntbus];
        if (Bus[internalindex].nshunt == 1)
          SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
                  "Bus %d: No support for more than 1 fixed shunt at bus",
                  Bus[internalindex].bus_i);
        Bus[internalindex].nshunt++;
        Bus[internalindex].gl = gshunt / ps->MVAbase;
        Bus[internalindex].bl = bshunt / ps->MVAbase;
      }
      break;
    case 3:
      sscanf(line,
             "%d, '%[^\t\']', %lf, %lf, %lf, %lf, %lf, %d, %lf, %lf, %lf, %lf, "
             "%lf, %lf, %d, %lf, %lf, %lf, %d, %lf",
             &Gen[geni].bus_i, Gen[geni].id, &Gen[geni].pg, &Gen[geni].qg,
             &Gen[geni].qt, &Gen[geni].qb, &Gen[geni].vs, &Gen[geni].ireg,
             &Gen[geni].mbase, &Gen[geni].zr, &Gen[geni].zx, &Gen[geni].rt,
             &Gen[geni].xt, &Gen[geni].gtap, &Gen[geni].status,
             &Gen[geni].rmpct, &Gen[geni].pt, &Gen[geni].pb, &Gen[geni].o1,
             &Gen[geni].f1);
#if defined DEBUGPS
      if (geni < 2)
        ierr = PetscPrintf(
            PETSC_COMM_SELF,
            "GENERATORData[%d] : %d, '%s', %lf, %lf, %lf, %lf, %lf, %d, %lf, "
            "%lf, %lf, %lf, %lf, %lf, %d, %lf, %lf, %lf, %d, %lf\n",
            geni, Gen[geni].bus_i, Gen[geni].id, Gen[geni].pg, Gen[geni].qg,
            Gen[geni].qt, Gen[geni].qb, Gen[geni].vs, Gen[geni].ireg,
            Gen[geni].mbase, Gen[geni].zr, Gen[geni].zx, Gen[geni].rt,
            Gen[geni].xt, Gen[geni].gtap, Gen[geni].status, Gen[geni].rmpct,
            Gen[geni].pt, Gen[geni].pb, Gen[geni].o1, Gen[geni].f1);
      CHKERRQ(ierr);
#endif
      Gen[geni].initial_status = Gen[geni].status;
      internalindex = busext2intmap[Gen[geni].bus_i];
      Gen[geni].internal_i = internalindex;
      Bus[internalindex].gidx[Bus[internalindex].ngen] = geni;
      if (Gen[geni].status == 1) {
        Gen[geni].pg /= ps->MVAbase;
        Gen[geni].qg /= ps->MVAbase;
        Gen[geni].pb /= ps->MVAbase;
        Gen[geni].pt /= ps->MVAbase;
        Gen[geni].qt /= ps->MVAbase;
        Gen[geni].qb /= ps->MVAbase;
        if ((Gen[geni].pb <= Gen[geni].pg) && (Gen[geni].pg <= Gen[geni].pt)) {
          Gen[geni].pgs = Gen[geni].pg;
        } else {
          Gen[geni].pgs = (Gen[geni].pb + Gen[geni].pt) / 2.0;
        }
        Bus[internalindex].qrange += (Gen[geni].qt - Gen[geni].qb);
        Bus[internalindex].qmintot += Gen[geni].qb;
        Bus[internalindex].Pgtot += PetscAbsScalar(Gen[geni].pg);
        Bus[internalindex].MVAbasetot += Gen[geni].mbase;
        Bus[internalindex].vm =
            Gen[geni].vs; /* Set bus voltage magnitude at generator set point
                             voltage. Assume all generators at the bus have the
                             same setpoint voltage in the daa file */
        Bus[internalindex].ngenON++;
        ps->NgenON++;
      }
      Bus[internalindex].ngen++;
      geni++;
      break;
    case 4:
      sscanf(line,
             "%d, %d, '%[^\t\']', %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, "
             "%lf, %d, %d, %lf, %d, %lf",
             &Branch[bri].fbus, &Branch[bri].tbus, Branch[bri].ckt,
             &Branch[bri].r, &Branch[bri].x, &Branch[bri].b, &Branch[bri].rateA,
             &Branch[bri].rateB, &Branch[bri].rateC, &Branch[bri].gi,
             &Branch[bri].bi, &Branch[bri].gj, &Branch[bri].bj,
             &Branch[bri].status, &MET, &Branch[bri].length, &Branch[bri].o1,
             &Branch[bri].f1);
#if defined DEBUGPS
      if (bri < 2)
        ierr = PetscPrintf(
            PETSC_COMM_SELF,
            "BRANCHData[%d] : %d, %d, '%s', %lf, %lf, %lf, %lf, %lf, %lf, %lf, "
            "%lf, %lf, %lf, %d, %d, %lf, %d, %lf\n",
            bri, Branch[bri].fbus, Branch[bri].tbus, Branch[bri].ckt,
            Branch[bri].r, Branch[bri].x, Branch[bri].b, Branch[bri].rateA,
            Branch[bri].rateB, Branch[bri].rateC, Branch[bri].gi,
            Branch[bri].bi, Branch[bri].gj, Branch[bri].bj, Branch[bri].status,
            MET, Branch[bri].length, Branch[bri].o1, Branch[bri].f1);
      CHKERRQ(ierr);
#endif
      if (Branch[bri].status)
        ps->NlineON++;
      Branch[bri].met = 1;
      internalindex = busext2intmap[Branch[bri].fbus];
      Branch[bri].internal_i = internalindex;
      internalindex = busext2intmap[Branch[bri].tbus];
      Branch[bri].internal_j = internalindex;
      Branch[bri].tapratio = 1.0;
      Branch[bri].phaseshift = 0.0;

      Branch[bri].rateA =
          (Branch[bri].rateA == 0 ? PETSC_INFINITY : Branch[bri].rateA);

      R = Branch[bri].r;
      X = Branch[bri].x;
      Bc = Branch[bri].b;

      Zm = R * R + X * X;
      G = R / Zm;
      B = -X / Zm;

      tap = Branch[bri].tapratio;
      shift = Branch[bri].phaseshift;
      tap2 = tap * tap;
      tapr = tap * cos(shift);
      tapi = tap * sin(shift);

      Branch[bri].yff[0] = G / tap2;
      Branch[bri].yff[1] = (B + Bc / 2.0) / tap2;

      Branch[bri].yft[0] = -(G * tapr - B * tapi) / tap2;
      Branch[bri].yft[1] = -(B * tapr + G * tapi) / tap2;

      Branch[bri].ytf[0] = -(G * tapr + B * tapi) / tap2;
      Branch[bri].ytf[1] = -(B * tapr - G * tapi) / tap2;

      Branch[bri].ytt[0] = G;
      Branch[bri].ytt[1] = B + Bc / 2.0;

      bri++;
      break;
    case 5:
      sscanf(line,
             "%d, %d, %d, '%[^\']', %d, %d, %d, %lf, %lf, %d, '%[^\t\']', %d, "
             "%d, %lf",
             &Branch[bri].fbus, &Branch[bri].tbus, &K, Branch[bri].ckt, &CW,
             &CZ, &CM, &MAG1, &MAG2, &NMETR, transname, &Branch[bri].status,
             &Branch[bri].o1, &Branch[bri].f1);
      if (K != 0) { /* Three winding transformer */
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP,
                "Found three winding transformer in the transformer data.\n\
                Three winding transformers are currently not supported.\n\
                Convert the raw file to MATPOWER format");
      }
      out = fgets(line, MAXLINE, fp);
      sscanf(line, "%lf, %lf, %lf", &Branch[bri].r, &Branch[bri].x,
             &Branch[bri].sbase12);

      out = fgets(line, MAXLINE, fp);
      sscanf(line,
             "%lf, %lf, %lf, %lf, %lf, %lf, %d, %d, %lf, %lf, %lf, %lf, %d, "
             "%d, %lf, %lf, %lf",
             &Branch[bri].tapratio, &NOMV1, &Branch[bri].phaseshift,
             &Branch[bri].rateA, &Branch[bri].rateB, &Branch[bri].rateC, &COD1,
             &CONT1, &RMA1, &RMI1, &VMA1, &VMI1, &NTP1, &TAB1, &CR1, &CX1,
             &CNXA1);

      out = fgets(line, MAXLINE, fp);
      sscanf(line, "%lf %lf", &WINDV2, &NOMV2);
#if defined DEBUGPS
      if (bri - Nline < 2)
        ierr = PetscPrintf(PETSC_COMM_SELF,
                           "TRANSFORMERData[%d] : %d, %d, %s, %lf, %lf, %lf, "
                           "%lf, %lf, %lf, %d, %d, %lf, %lf, %lf\n",
                           bri - Nline, Branch[bri].fbus, Branch[bri].tbus,
                           Branch[bri].ckt, Branch[bri].r, Branch[bri].x,
                           Branch[bri].b, Branch[bri].rateA, Branch[bri].rateB,
                           Branch[bri].rateC, Branch[bri].status,
                           Branch[bri].o1, Branch[bri].f1, Branch[bri].tapratio,
                           Branch[bri].phaseshift);
      CHKERRQ(ierr);
#endif
      internalindex = busext2intmap[Branch[bri].fbus];
      Branch[bri].internal_i = internalindex;
      internalindex = busext2intmap[Branch[bri].tbus];
      Branch[bri].internal_j = internalindex;

      Branch[bri].b = 0;

      R = Branch[bri].r;
      X = Branch[bri].x;
      Bc = Branch[bri].b;

      Zm = R * R + X * X;
      G = R / Zm;
      B = -X / Zm;

      tap = Branch[bri].tapratio;
      shift = Branch[bri].phaseshift;
      tap2 = tap * tap;
      tapr = tap * cos(shift);
      tapi = tap * sin(shift);

      Branch[bri].yff[0] = G / tap2;
      Branch[bri].yff[1] = (B + Bc / 2.0) / tap2;

      Branch[bri].yft[0] = -(G * tapr - B * tapi) / tap2;
      Branch[bri].yft[1] = -(B * tapr + G * tapi) / tap2;

      Branch[bri].ytf[0] = -(G * tapr + B * tapi) / tap2;
      Branch[bri].ytf[1] = -(B * tapr - G * tapi) / tap2;

      Branch[bri].ytt[0] = G;
      Branch[bri].ytt[1] = B + Bc / 2.0;

      bri++;
      break;
    default:
      // Todo: impliment parser for other data
      break;
    }
  }
  fclose(fp);

  PetscFunctionReturn(0);
}

/*
  PSReadMatPowerData - Reads the MATPOWER data file and populates the PS object

  Input Parameter
+ ps      - The power system object ps
- netfile - Name of the power system network data file in MATPOWER data format.

*/
PetscErrorCode PSReadMatPowerData(PS ps, const char netfile[]) {
  FILE *fp;
  PetscErrorCode ierr;
  PSBUS Bus;
  PSLOAD Load;
  PSGEN Gen;
  PSLINE Branch;
  PetscInt line_counter = 0, linenum;
  PetscInt
      bus_start_line = -1,
      bus_end_line =
          -1; /* xx_end_line points to the next line after the record ends */
  PetscInt gen_start_line = -1, gen_end_line = -1;
  PetscInt br_start_line = -1, br_end_line = -1;
  PetscInt gencost_start_line = -1, gencost_end_line = -1;
  PetscInt genfuel_start_line = -1, genfuel_end_line = -1;
  PetscInt loadcost_start_line = -1, loadcost_end_line = -1;
  /* Number of blank lines in bus, gen, br, gencost, and genfuel branch arrays
   */
  PetscInt bus_nblank_lines = 0, gen_nblank_lines = 0;
  PetscInt br_nblank_lines = 0;
  char line[MAXLINE];
  PetscInt loadi = 0, geni = 0, bri = 0, busi = 0, gencosti = 0, genfueli = 0,
           loadcosti = 0, i;
  PetscInt extbusnum, bustype_i;
  PetscScalar Pd, Qd;
  PetscInt intbusnum;
  char *str;
  char *out;

  PetscFunctionBegin;

  if (ps->comm->type != PETSC_COMM_SELF && ps->comm->rank != 0) {
    ps->Nline = ps->Nbus = ps->Ngen = ps->Nload = 0;
    PetscFunctionReturn(0);
  }

  fp = fopen(netfile, "r");
  /* Check for valid file */
  if (fp == NULL) {
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Cannot open file %s",
            netfile);
    CHKERRQ(ierr);
  }

  /* Copy the network file name */
  ierr = PetscStrcpy(ps->net_file_name, netfile);
  CHKERRQ(ierr);

  ps->Nload = 0;
  ps->maxbusnum = -1;
  while ((out = fgets(line, MAXLINE, fp)) != NULL) {
    if (strstr(line, "mpc.baseMVA")) {
      /* Read base MVA */
      str = strtok(line, " =;");
      str = strtok(NULL, " =;");
      sscanf(str, "%lf", &ps->MVAbase);
    }
    if (strstr(line, "mpc.bus") != NULL && bus_start_line == -1)
      bus_start_line = line_counter + 1; /* Bus data starts from next line */
    if (strstr(line, "mpc.gen") != NULL && gen_start_line == -1)
      gen_start_line =
          line_counter + 1; /* Generator data starts from next line */
    if (strstr(line, "mpc.branch") != NULL)
      br_start_line = line_counter + 1; /* Branch data starts from next line */
    if (strstr(line, "mpc.gencost") != NULL)
      gencost_start_line =
          line_counter + 1; /* Gen cost data starts from next line */
    if (strstr(line, "mpc.loadcost") != NULL) {
      loadcost_start_line =
          line_counter + 1; /* Load cost data starts from next line */
    }

    if (strstr(line, "mpc.genfuel") != NULL) {
      genfuel_start_line =
          line_counter + 1; /* Gen fuel data starts from next line */
    }

    if (strstr(line, "};") != NULL) {
      if (genfuel_start_line != -1 && genfuel_end_line == -1)
        genfuel_end_line = line_counter;
    }

    if (strstr(line, "];") != NULL) {
      if (bus_start_line != -1 && bus_end_line == -1)
        bus_end_line = line_counter;
      if (gen_start_line != -1 && gen_end_line == -1)
        gen_end_line = line_counter;
      if (br_start_line != -1 && br_end_line == -1)
        br_end_line = line_counter;
      if (gencost_start_line != -1 && gencost_end_line == -1)
        gencost_end_line = line_counter;
      if (loadcost_start_line != -1 && loadcost_end_line == -1)
        loadcost_end_line = line_counter;
    }

    if (bus_start_line != -1 && bus_end_line == -1) {
      if (strcmp(line, "\n") == 0 || strcmp(line, "\r\n") == 0)
        bus_nblank_lines++;
    }
    if (gen_start_line != -1 && gen_end_line == -1) {
      if (strcmp(line, "\n") == 0 || strcmp(line, "\r\n") == 0)
        gen_nblank_lines++;
    }

    if (br_start_line != -1 && br_end_line == -1) {
      if (strcmp(line, "\n") == 0 || strcmp(line, "\r\n") == 0)
        br_nblank_lines++;
    }

    /* Count the number of pq loads */
    if (bus_start_line != -1 && line_counter >= bus_start_line &&
        bus_end_line == -1) {
      sscanf(line, "%d %d %lf %lf", &extbusnum, &bustype_i, &Pd, &Qd);
      if (!((Pd == 0.0) && (Qd == 0.0)))
        ps->Nload++;
      if (extbusnum > ps->maxbusnum)
        ps->maxbusnum = extbusnum;
    }
    line_counter++;
  }
  fclose(fp);

  ps->Nbus = ps->nbus = bus_end_line - bus_start_line - bus_nblank_lines;
  ps->Ngen = ps->ngen = gen_end_line - gen_start_line - gen_nblank_lines;
  ps->Nline = ps->nline = br_end_line - br_start_line - br_nblank_lines;
  ps->nload = ps->Nload;
  ps->NgenON = 0;
  ps->NlineON = 0;

#if defined DEBUGPS
  ierr = PetscPrintf(
      PETSC_COMM_SELF,
      "System summary : Nbuses = %d, Ngen = %d, Nload = %d, Nbranch = %d\n",
      ps->Nbus, ps->Ngen, ps->Nload, ps->Nline);
  CHKERRQ(ierr);
#endif
  ierr = PetscCalloc1(ps->Nbus, &ps->bus);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(ps->Ngen, &ps->gen);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(ps->Nload, &ps->load);
  CHKERRQ(ierr);
  ierr = PetscCalloc1(ps->Nline, &ps->line);
  CHKERRQ(ierr);
  Bus = ps->bus;
  Gen = ps->gen;
  Load = ps->load;
  Branch = ps->line;

  for (i = 0; i < ps->Nbus; i++) {
    ps->bus[i].ngen = ps->bus[i].nload = ps->bus[i].ngenON = ps->bus[i].nshunt =
        0;
    ps->bus[i].qrange = ps->bus[i].qmintot = ps->bus[i].Pgtot =
        ps->bus[i].MVAbasetot = 0.0;
  }

  /*  ierr = PetscPrintf(PETSC_COMM_SELF,"Rank[%d] maxbusnum
   * %d",ps->comm->rank,ps->maxbusnum);CHKERRQ(ierr); */

  /* Allocate external to internal bus number mapping array */
  PetscInt *busext2intmap;
  ierr = PetscCalloc1(ps->maxbusnum + 1, &ps->busext2intmap);
  CHKERRQ(ierr);
  busext2intmap = ps->busext2intmap;
  for (i = 0; i < ps->maxbusnum + 1; i++)
    busext2intmap[i] = -1;

  fp = fopen(netfile, "r");
  /* Reading data */
  for (i = 0; i < line_counter; i++) {
    out = fgets(line, MAXLINE, fp);

    if ((i >= bus_start_line) && (i < bus_end_line)) {
      if (strcmp(line, "\n") == 0 || strcmp(line, "\r\n") == 0)
        continue;
      /* Bus data */
      sscanf(line, "%d %d %lf %lf %lf %lf %d %lf %lf %lf %d %lf %lf",
             &Bus[busi].bus_i, &Bus[busi].ide, &Pd, &Qd, &Bus[busi].gl,
             &Bus[busi].bl, &Bus[busi].area, &Bus[busi].vm, &Bus[busi].va,
             &Bus[busi].basekV, &Bus[busi].zone, &Bus[busi].Vmax,
             &Bus[busi].Vmin);

      Bus[busi].Vmax = Bus[busi].Vmax == 0 ? 1.1 : Bus[busi].Vmax;
      Bus[busi].Vmin = Bus[busi].Vmin == 0 ? 0.9 : Bus[busi].Vmin;

      // Convert angle to radians
      Bus[busi].va *= PETSC_PI / 180.0;

      if (Bus[busi].ide == REF_BUS)
        ps->Nref++;
      Bus[busi].internal_i = busi;
      busext2intmap[Bus[busi].bus_i] = busi;
      /* Convert bl and gl to per unit */
      Bus[busi].bl = Bus[busi].bl / ps->MVAbase;
      Bus[busi].gl = Bus[busi].gl / ps->MVAbase;
      Bus[busi].nshunt++;

      if (!((Pd == 0.0) && (Qd == 0.0))) {
        Load[loadi].bus_i = Bus[busi].bus_i;
        Load[loadi].status = 1;
        Load[loadi].pl = Pd / ps->MVAbase;
        Load[loadi].ql = Qd / ps->MVAbase;
        /* Some defaults for load shed */
        Load[loadi].loss_frac = 1.0;
        Load[loadi].loss_cost = BOGUSLOSSCOST;
        Load[loadi].area = Bus[busi].area;
        Load[loadi].internal_i = busi;

        intbusnum = busext2intmap[Load[loadi].bus_i];

        /* MatPower does not have ids for loads. Using Bus[i].nload as the id */
        snprintf(Load[loadi].id, 3, "%-2d", 1 + Bus[intbusnum].nload);

        Bus[busi].lidx[Bus[busi].nload++] = loadi;
        if (Bus[busi].nload > NLOAD_AT_BUS_MAX)
          SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP,
                  "Exceeded maximum number of loads allowed at bus");
        loadi++;
      }
      busi++;
    }

    /* Read generator data */
    if (i >= gen_start_line && i < gen_end_line) {
      if (strcmp(line, "\n") == 0 || strcmp(line, "\r\n") == 0)
        continue;
      sscanf(line,
             "%d %lf %lf %lf %lf %lf %lf %d %lf %lf %lf %lf %lf %lf %lf %lf "
             "%lf %lf %lf %lf %lf",
             &Gen[geni].bus_i, &Gen[geni].pg, &Gen[geni].qg, &Gen[geni].qt,
             &Gen[geni].qb, &Gen[geni].vs, &Gen[geni].mbase, &Gen[geni].status,
             &Gen[geni].pt, &Gen[geni].pb, &Gen[geni].pc1, &Gen[geni].pc2,
             &Gen[geni].qc1min, &Gen[geni].qc1max, &Gen[geni].qc2min,
             &Gen[geni].qc2max, &Gen[geni].ramp_rate_min,
             &Gen[geni].ramp_rate_10min, &Gen[geni].ramp_rate_30min,
             &Gen[geni].ramp_rate_min_mvar, &Gen[geni].apf);

      intbusnum = busext2intmap[Gen[geni].bus_i];

      Gen[geni].qt = Gen[geni].qt > 1e10 ? PETSC_INFINITY : Gen[geni].qt;
      Gen[geni].qb = Gen[geni].qb < -1e10 ? PETSC_NINFINITY : Gen[geni].qb;
      Gen[geni].pt = Gen[geni].pt > 1e10 ? PETSC_INFINITY : Gen[geni].pt;
      Gen[geni].pb = Gen[geni].pb < -1e10 ? PETSC_NINFINITY : Gen[geni].pb;

      Gen[geni].isrenewable = PETSC_FALSE;

      Gen[geni].initial_status = Gen[geni].status;

      Gen[geni].pg = Gen[geni].pg / ps->MVAbase;
      Gen[geni].qg = Gen[geni].qg / ps->MVAbase;
      Gen[geni].pb = Gen[geni].pb / ps->MVAbase;
      Gen[geni].pt = Gen[geni].pt / ps->MVAbase;
      Gen[geni].qt = Gen[geni].qt / ps->MVAbase;
      Gen[geni].qb = Gen[geni].qb / ps->MVAbase;
      if ((Gen[geni].pb <= Gen[geni].pg) && (Gen[geni].pg <= Gen[geni].pt)) {
        Gen[geni].pgs = Gen[geni].pg;
      } else {
        Gen[geni].pgs = (Gen[geni].pb + Gen[geni].pt) / 2.0;
      }
      Gen[geni].pc1 = Gen[geni].pc1 / ps->MVAbase;
      Gen[geni].pc2 = Gen[geni].pc2 / ps->MVAbase;
      Gen[geni].qc1min = Gen[geni].qc1min / ps->MVAbase;
      Gen[geni].qc1max = Gen[geni].qc1max / ps->MVAbase;
      Gen[geni].qc2min = Gen[geni].qc2min / ps->MVAbase;
      Gen[geni].qc2max = Gen[geni].qc2max / ps->MVAbase;
      Gen[geni].ramp_rate_min = Gen[geni].ramp_rate_min / ps->MVAbase;
      Gen[geni].ramp_rate_10min = Gen[geni].ramp_rate_10min / ps->MVAbase;
      Gen[geni].ramp_rate_30min = Gen[geni].ramp_rate_30min / ps->MVAbase;
      Gen[geni].ramp_rate_min_mvar = Gen[geni].ramp_rate_min_mvar / ps->MVAbase;

      if (genfuel_start_line == -1 && genfuel_end_line == -1) {
        /* Fuel source not defined in file so use COAL RAMP RATES FOR ALL
         * GENERATORS */
        Gen[geni].genfuel_type = GENFUEL_UNDEFINED;
        Gen[geni].ramp_rate_min =
            GENRAMPRATE_COAL / ps->MVAbase; /* Defaults to COAL ramp rate */
        Gen[geni].ramp_rate_10min = Gen[genfueli].ramp_rate_min * 10;
        Gen[geni].ramp_rate_30min = Gen[genfueli].ramp_rate_min * 30;
      }

      if (Gen[geni].status) {
        Bus[intbusnum].qrange += (Gen[geni].qt - Gen[geni].qb);
        Bus[intbusnum].qmintot += Gen[geni].qb;
        Bus[intbusnum].Pgtot += PetscAbsScalar(Gen[geni].pg);
        Bus[intbusnum].MVAbasetot += Gen[geni].mbase;
        Bus[intbusnum].ngenON++;
        ps->NgenON++;
      } else {
        Gen[geni].pg = Gen[geni].qg = 0.0;
      }

      Gen[geni].internal_i = intbusnum;

      /* MatPower does not have ids for generators. Using Bus[i].ngen as the id
       */
      snprintf(Gen[geni].id, 3, "%-2d", 1 + Bus[intbusnum].ngen);

      Bus[intbusnum].gidx[Bus[intbusnum].ngen++] = geni;

      //      Bus[intbusnum].vm = Gen[geni].vs;

#if defined DEBUGPS
      ierr = PetscPrintf(PETSC_COMM_SELF, "%d %d %d %s\n", Gen[geni].status,
                         intbusnum, Bus[intbusnum].ngen, line);
      CHKERRQ(ierr);
#endif
      if (Bus[intbusnum].ngen > NGEN_AT_BUS_MAX)
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP,
                "Exceeded maximum number of generators allowed at bus");
      geni++;
    }

    //    ierr = PetscPrintf(PETSC_COMM_SELF,"Came here\n");CHKERRQ(ierr);

    /* Read generator cost data */
    if (i >= gencost_start_line && i < gencost_end_line) {
      if (strcmp(line, "\n") == 0 || strcmp(line, "\r\n") == 0)
        continue;
      sscanf(line, "%d %lf %lf %d %lf %lf %lf", &Gen[gencosti].cost_model,
             &Gen[gencosti].cost_startup, &Gen[gencosti].cost_shutdown,
             &Gen[gencosti].cost_ncoeffs, &Gen[gencosti].cost_alpha,
             &Gen[gencosti].cost_beta, &Gen[gencosti].cost_gamma);
      gencosti++;
    }

    /* Read load cost data */
    if (i >= loadcost_start_line && i < loadcost_end_line) {
      if (strcmp(line, "\n") == 0 || strcmp(line, "\r\n") == 0)
        continue;
      sscanf(line, "%lf %lf", &Load[loadcosti].loss_frac,
             &Load[loadcosti].loss_cost);
      Load[loadcosti].loss_frac /= 100; /* Convert to fraction */
      loadcosti++;
    }

    /* Read generator fuel data */
    if (i >= genfuel_start_line && i < genfuel_end_line) {
      if (strcmp(line, "\n") == 0 || strcmp(line, "\r\n") == 0)
        continue;
      if (strstr(line, "coal") != NULL) {
        Gen[genfueli].genfuel_type = GENFUEL_COAL;
        Gen[genfueli].ramp_rate_min = GENRAMPRATE_COAL / ps->MVAbase;
        Gen[genfueli].ramp_rate_10min = Gen[genfueli].ramp_rate_min * 10;
        Gen[genfueli].ramp_rate_30min = Gen[genfueli].ramp_rate_min * 30;
        ps->ngencoal++;
      } else if (strstr(line, "wind") != NULL) {
        Gen[genfueli].genfuel_type = GENFUEL_WIND;
        Gen[genfueli].ramp_rate_min = GENRAMPRATE_WIND / ps->MVAbase;
        Gen[genfueli].ramp_rate_10min = Gen[genfueli].ramp_rate_min * 10;
        Gen[genfueli].ramp_rate_30min = Gen[genfueli].ramp_rate_min * 30;
        Gen[genfueli].pb = 0.0; /* Set lower Pg limit to 0.0 so that wind power
                                   can be curtailed if need be */
        ps->ngenwind++;
        ps->ngenrenew++;
        Gen[genfueli].isrenewable = PETSC_TRUE;
      } else if (strstr(line, "ng") != NULL) {
        Gen[genfueli].genfuel_type = GENFUEL_NG;
        Gen[genfueli].ramp_rate_min = GENRAMPRATE_NG / ps->MVAbase;
        Gen[genfueli].ramp_rate_10min = Gen[genfueli].ramp_rate_min * 10;
        Gen[genfueli].ramp_rate_30min = Gen[genfueli].ramp_rate_min * 30;
        ps->ngenng++;
      } else if (strstr(line, "solar") != NULL) {
        Gen[genfueli].genfuel_type = GENFUEL_SOLAR;
        Gen[genfueli].ramp_rate_min = GENRAMPRATE_SOLAR / ps->MVAbase;
        Gen[genfueli].ramp_rate_10min = Gen[genfueli].ramp_rate_min * 10;
        Gen[genfueli].ramp_rate_30min = Gen[genfueli].ramp_rate_min * 30;
	Gen[genfueli].pb = 0.0; /* Set lower Pg limit to 0.0 so that power
                                   can be curtailed if need be */
        ps->ngensolar++;
        ps->ngenrenew++;
        Gen[genfueli].isrenewable = PETSC_TRUE;
      } else if (strstr(line, "nuclear") != NULL) {
        Gen[genfueli].genfuel_type = GENFUEL_NUCLEAR;
        Gen[genfueli].ramp_rate_min = GENRAMPRATE_NUCLEAR / ps->MVAbase;
        Gen[genfueli].ramp_rate_10min = Gen[genfueli].ramp_rate_min * 10;
        Gen[genfueli].ramp_rate_30min = Gen[genfueli].ramp_rate_min * 30;
        ps->ngennuclear++;
      } else if (strstr(line, "hydro") != NULL) {
        Gen[genfueli].genfuel_type = GENFUEL_HYDRO;
        Gen[genfueli].ramp_rate_min = GENRAMPRATE_HYDRO / ps->MVAbase;
        Gen[genfueli].ramp_rate_10min = Gen[genfueli].ramp_rate_min * 10;
        Gen[genfueli].ramp_rate_30min = Gen[genfueli].ramp_rate_min * 30;
        ps->ngenhydro++;
	ps->ngenrenew++;
        Gen[genfueli].isrenewable = PETSC_TRUE;
      } else {
        Gen[genfueli].genfuel_type = GENFUEL_UNDEFINED;
        Gen[genfueli].ramp_rate_min =
            GENRAMPRATE_COAL / ps->MVAbase; /* Defaults to COAL ramp rate */
        Gen[genfueli].ramp_rate_10min = Gen[genfueli].ramp_rate_min * 10;
        Gen[genfueli].ramp_rate_30min = Gen[genfueli].ramp_rate_min * 30;
        ps->ngenundefined++;
      }
      genfueli++;
    }

    if (i >= br_start_line && i < br_end_line) {
      if (strcmp(line, "\n") == 0 || strcmp(line, "\r\n") == 0)
        continue;
      sscanf(line, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %d",
             &Branch[bri].fbus, &Branch[bri].tbus, &Branch[bri].r,
             &Branch[bri].x, &Branch[bri].b, &Branch[bri].rateA,
             &Branch[bri].rateB, &Branch[bri].rateC, &Branch[bri].tapratio,
             &Branch[bri].phaseshift, &Branch[bri].status);
      if (Branch[bri].status)
        ps->NlineON++;
      if (!Branch[bri].tapratio)
        Branch[bri].tapratio = 1.0;
      Branch[bri].phaseshift *= PETSC_PI / 180.0;

      Branch[bri].reversed_ends = PETSC_FALSE;
      if (Branch[bri].fbus > Branch[bri].tbus) {
        /* Swap ends.. Some solvers (like IPOPT) which
           require symmetric matrices (only lower or upper
           triangle) complain that Hessian is incorrect
           if fbus > tbus. Hence, this workaround.
        */
        PetscInt temp;
        temp = Branch[bri].fbus;
        Branch[bri].fbus = Branch[bri].tbus;
        Branch[bri].tbus = temp;
        Branch[bri].reversed_ends = PETSC_TRUE;
      }
      Branch[bri].rateA =
          (Branch[bri].rateA == 0 ? PETSC_INFINITY : Branch[bri].rateA);
      intbusnum = busext2intmap[Branch[bri].fbus];
      Branch[bri].internal_i = intbusnum;

      intbusnum = busext2intmap[Branch[bri].tbus];
      Branch[bri].internal_j = intbusnum;

      PetscInt lineididx = 0;
      for (linenum = 0; linenum < bri - 1; linenum++) {
        if (Branch[bri].internal_i == Branch[linenum].internal_i &&
            Branch[bri].internal_j == Branch[linenum].internal_j)
          lineididx += 1;
      }

      /* MatPower does not have ids for lines. Using bri+1 as the id */
      snprintf(Branch[bri].ckt, 3, "%-2d", 1 + lineididx);

      /* Compute self and transfer admittances */
      PetscScalar R, X, Bc, B, G, Zm, tap, shift, tap2, tapr, tapi;
      R = Branch[bri].r;
      X = Branch[bri].x;
      Bc = Branch[bri].b;

      Zm = R * R + X * X;
      G = R / Zm;
      B = -X / Zm;

      tap = Branch[bri].tapratio;
      shift = Branch[bri].phaseshift;
      tap2 = tap * tap;
      tapr = tap * cos(shift);
      tapi = tap * sin(shift);

      if (!Branch[bri].reversed_ends) {
        Branch[bri].yff[0] = G / tap2;
        Branch[bri].yff[1] = (B + Bc / 2.0) / tap2;

        Branch[bri].yft[0] = -(G * tapr - B * tapi) / tap2;
        Branch[bri].yft[1] = -(B * tapr + G * tapi) / tap2;

        Branch[bri].ytf[0] = -(G * tapr + B * tapi) / tap2;
        Branch[bri].ytf[1] = -(B * tapr - G * tapi) / tap2;

        Branch[bri].ytt[0] = G;
        Branch[bri].ytt[1] = B + Bc / 2.0;

        /* For DC formulation */
        Branch[bri].bdc = (1 / X) / tap;
        Branch[bri].pshift = -(1 / X) * shift;
      } else {
        Branch[bri].ytt[0] = G / tap2;
        Branch[bri].ytt[1] = (B + Bc / 2.0) / tap2;

        Branch[bri].ytf[0] = -(G * tapr - B * tapi) / tap2;
        Branch[bri].ytf[1] = -(B * tapr + G * tapi) / tap2;

        Branch[bri].yft[0] = -(G * tapr + B * tapi) / tap2;
        Branch[bri].yft[1] = -(B * tapr - G * tapi) / tap2;

        Branch[bri].yff[0] = G;
        Branch[bri].yff[1] = B + Bc / 2.0;

        /* For DC formulation */
        Branch[bri].bdc = (1 / X) / tap;
        Branch[bri].pshift = -(1 / X) * shift;
      }
      bri++;
    }
  }
  fclose(fp);

  PetscFunctionReturn(0);
}

/*
  PSReadGICData - Reads GIC data

Input Parameters:
. ps - the PS object

 Notes: The GIC data format should be the same as that
        given for the Electric Grid Data repository cases.
        The function reads the GIC data and populates the
        substation data in PS object
*/
PetscErrorCode PSReadGICData(PS ps) {
  PetscErrorCode ierr;
  FILE *fp;
  char line[MAXLINE];
  char *out;
  PSBUS bus;
  PSSUBST subst;
  int fieldsread = 0;
  int subst_num, bus_num;
  PetscInt nconnlines;
  const PSLINE *connlines;
  PSLINE branch;

  PetscFunctionBegin;

  fp = fopen(ps->gic_file_name, "r");
  /* Check for valid file */
  if (fp == NULL) {
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Cannot open gic file %s",
            ps->gic_file_name);
    CHKERRQ(ierr);
  }
  /* Skip first line */
  out = fgets(line, MAXLINE, fp);
  /* Allocate substation data,
   assume number of substations = number of buses,
   but this does not neccessarily be the case and
   hence we also keep track of the number of substations.
   Usually, number of substations < number of buses */
  ierr = PetscCalloc1(ps->Nbus, &ps->substations);
  CHKERRQ(ierr);
  ps->nsubstations = 0;

  /* Start reading substation data */
  while ((out = fgets(line, MAXLINE, fp)) != NULL) {
    if (strstr(line, "0 /") != NULL) {
      fieldsread++;
      continue;
    }

    if (fieldsread == 2)
      break;

    if (fieldsread == 0) {
      subst = &ps->substations[ps->nsubstations];
      sscanf(line, "%d '%[^\']' %*d %lf %lf", &subst->num, subst->name,
             &subst->longlat[1], &subst->longlat[0]);
      ps->nsubstations++;
      subst->nbus = subst->nkvlevels = 0;
    } else if (fieldsread == 1) {
      bool found = false;
      sscanf(line, "%d %d", &bus_num, &subst_num);
      // This is really ugly and slow, need to fix this idx business
      int idx = -1;
      for (int i = 0; i < ps->nsubstations; i++) {
        if (subst_num == ps->substations[i].num) {
          idx = i;
          break;
        }
      }
      subst = &ps->substations[idx];
      bus = &ps->bus[ps->busext2intmap[bus_num]];
      subst->bus[subst->nbus++] = bus;

      for (int j = 0; j < subst->nkvlevels; j++) {
        if (PetscAbsScalar(bus->basekV - subst->kvlevels[j]) < 1e-6) {
          found = true;
          break;
        }
      }
      if (!found) {
        subst->kvlevels[subst->nkvlevels++] = bus->basekV;
      }

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

  fclose(fp);

  PetscFunctionReturn(0);
}
