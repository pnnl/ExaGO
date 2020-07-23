static char help[] = "Prototype version of OPFLOW code for CUDA/GPU/Accelerator implementation\n\
                      This uses current-balance-cartesian-2 formulation\n\n";

#include <scopflow_config.h>
#include <private/psimpl.h>
#include <private/opflowimpl.h>

/* The parameters are packed in each busparams[i] as follows
busparams[i] = {MVAbase
                # bus type
                # of gens
                # of loads
		# of lines
		# shunt (gl,bl for each bus)
		# status,Pg,Qg,alpha,beta,gamma for each gen
		# Pd,Qd for each load
                }


The parameters are packed in each branchparams[i] as follows
branchparams[i] = {
		# status,Gff,Bff,Gft,Bft,Gtf,Btf,Gtt,Btt,rateA for each line

xidxarr is an array of arrays that contains the locations for the variables at
each bus and branch. Each xidxarr[i] has locations for its own variables AND
the locations of other variables needed to compute its equation(s).
*/

PetscLogEvent Objfun_logger;
PetscLogEvent Gradfun_logger;
PetscLogEvent Eqconsfunc_logger;
PetscLogEvent Ineqconsfunc_logger;

int ComputeBusObjective(double *busparams,int *xidx,const double *x,double *obj)
{
  double MVAbase = busparams[0];
  int ngen=(int)busparams[2];
  int gstatus;
  int i;
  int locxgen = xidx[0]+2; // +2 since first two variables are voltages */
  double alpha, beta, gamma;
  double pg;

  busparams += 7;

  for(i=0; i < ngen; i++) {
    gstatus = (int)busparams[0];
    if(!gstatus) continue;
    pg = x[locxgen]; 
    alpha = busparams[3];
    beta  = busparams[4];
    gamma = busparams[5];

    *obj += alpha*pg*pg*MVAbase*MVAbase + beta*pg*MVAbase + gamma;

    busparams += 6;
    locxgen += 2;
  }

  //  printf("Output obj %f\n",*obj);
  
  return(0);
}
/* Computes the objective (generation cost) for the optimal power flow
   obj = \sum_{i=1}^ng cost_alpha*(Pg*MVAbase)^2 + cost_beta*Pg*MVAbase + cost_gamma
   We loop over all buses and add the generation cost
*/
PetscErrorCode ComputeGlobalObjective(OPFLOW opflow,double **busparamsarr,int **xidxarr,double *obj,Vec X)
{
  PS             ps=opflow->ps;
  PetscErrorCode ierr;
  PetscInt       i;
  const double   *x;
  double         busobj;

  PetscFunctionBegin;

  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);

  for(i=0; i < ps->nbus; i++) {
    busobj = 0.0;
    /* Compute the bus objective value (cost for generation at bus i) */
    ierr = ComputeBusObjective(busparamsarr[i],xidxarr[ps->nline+i],x,&busobj);
    *obj += busobj; /* Add the bus objective value */
  }

  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);

  if ((*obj - opflow->obj)/(*obj) < 1e-6) {
    ierr = PetscPrintf(PETSC_COMM_SELF,"Objective function calculation correct\n");CHKERRQ(ierr);
  } else {
    ierr = PetscPrintf(PETSC_COMM_SELF,"Objective function %f is incorrect\n",*obj);CHKERRQ(ierr);
  }

  //  ierr = PetscPrintf(PETSC_COMM_SELF,"The objective function value is %f\n",*obj);

  PetscFunctionReturn(0);
}

int ComputeBusInequalityConstraints(double *busparams, int *xidx, int *fineqidx, const double *x, double *fineq)
{
  double Vr    = x[xidx[0]];
  double Vi    = x[xidx[0]+1];
  double Vm2   = Vr*Vr + Vi*Vi; 

  fineq[fineqidx[0]] = Vm2;

  return 0;
}

int ComputeBusEqualityConstraints(double *busparams,int *xidx,int *feqidx,const double *x, double *feq)
{
  int ngen = (int)busparams[2];
  int nload = (int)busparams[3];
  int nlines = (int)busparams[4];
  double gl  = busparams[5];
  double bl  = busparams[6];
  int    idxIr = feqidx[0],idxIi=feqidx[0]+1;
  int    idxVr = xidx[0],idxVi=xidx[0]+1;
  double Vr    = x[idxVr];
  double Vi    = x[idxVi];
  double Vm2   = Vr*Vr + Vi*Vi; 
  int    i;
  double Ir,Ii;
  PetscInt gstatus,lstatus,linestatus;
  int    locxgen = xidx[0]+2;
  double Pg,Qg,Pd,Qd;

  /* Shunt current injections */
  feq[idxIr] = gl*Vr - bl*Vi;
  feq[idxIi] = bl*Vr + gl*Vi;

  busparams += 7;
  /* Generation current injection */
  for (i=0; i < ngen; i++) {
    gstatus = (int)busparams[0];
    if(!gstatus) continue;
    Pg = x[locxgen];
    Qg = x[locxgen+1];
      
    feq[idxIr] -= (Pg*Vr + Qg*Vi)/Vm2;
    feq[idxIi] -= (-Qg*Vr + Pg*Vi)/Vm2;

    locxgen += 2;
    
    busparams += 6;
  }

  /* Load current injection */
  for (i=0; i < nload; i++) {
    lstatus = (int)busparams[0];
    Pd = busparams[1];
    Qd = busparams[2];
    
    feq[idxIr] += (Pd*Vr + Qd*Vi)/Vm2;
    feq[idxIi] += (-Qd*Vr + Pd*Vi)/Vm2;
    
    busparams += 3;
  }
  
  /* Line current injection */
  for (i=0; i < nlines; i++) {
    linestatus = (int)busparams[0];
    if (!linestatus) continue;
    
    Ir = x[xidx[1+i]]; Ii = x[xidx[1+i]+1];
    
    feq[idxIr] += Ir;
    feq[idxIi] += Ii;
    
    busparams++;
  }

  // Add equality constraint for referencce angle (TO DO)
  return 0;
}

int ComputeBranchInequalityConstraints(double *branchparams, int *xidx, int *fineqidx, const double *x, double *fineq)
{
  int status=branchparams[0];
  double rateA = branchparams[9];
  double Irf = x[xidx[0]], Iif = x[xidx[0]+1];
  double Irt = x[xidx[0]+2], Iit = x[xidx[0]+3];

  if(!status || rateA > 1e5) return 0;

  fineq[fineqidx[0]] = Irf*Irf + Iif*Iif;
  fineq[fineqidx[0]+1] = Irt*Irt + Iit*Iit;

  return 0;
}

int ComputeBranchEqualityConstraints(double *branchparams,int *xidx,int *feqidx,const double *x, double *f)
{
  int status=branchparams[0];
  double Gff,Bff,Gtt,Btt,Gft,Bft,Gtf,Btf;
  double Irf = x[xidx[0]], Iif = x[xidx[0]+1];
  double Irt = x[xidx[0]+2], Iit = x[xidx[0]+3];

  double Vrf = x[xidx[1]], Vif = x[xidx[1]+1];
  double Vrt = x[xidx[2]], Vit = x[xidx[2]+1];

  if(!status) return 0;

  Gff = branchparams[1];
  Bff = branchparams[2];
  Gft = branchparams[3];
  Bft = branchparams[4];
  Gtf = branchparams[5];
  Btf = branchparams[6];
  Gtt = branchparams[7];
  Btt = branchparams[8];

  f[feqidx[0]]   = Gff*Vrf - Bff*Vif + Gft*Vrt - Bft*Vit - Irf;
  f[feqidx[0]+1] = Bff*Vrf + Gff*Vif + Bft*Vrt + Gft*Vit - Iif;
  f[feqidx[0]+2] = Gtt*Vrt - Btt*Vit + Gtf*Vrf - Btf*Vif - Irt;
  f[feqidx[0]+3] = Btt*Vrt + Gtt*Vit + Btf*Vrf + Gtf*Vif - Iit;
  
  return 0;
}

PetscErrorCode ComputeGlobalEqualityConstraints(OPFLOW opflow,double **branchparamsarr,double **busparamsarr,int **xidxarr, int **feqidxarr,Vec X, Vec Feq)
{
  PetscErrorCode ierr;
  PS             ps=opflow->ps;
  PetscInt       i;
  const double   *x;
  double         *feq;
  PSBUS          bus;
  PSLINE         line;
  PetscInt       idxIrf,idxIif,idxIrt,idxIit;

  PetscFunctionBegin;
  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecGetArray(Feq,&feq);CHKERRQ(ierr);

  for(i=0; i < ps->nline; i++) {
    ierr = ComputeBranchEqualityConstraints(branchparamsarr[i],xidxarr[i],feqidxarr[i],x,feq);
    /* Sanity check */
    /* Display residuals for branches that are erroneous */
    line = &ps->line[i];
    if(!line->status) continue;
    idxIrf = feqidxarr[i][0]; idxIif = idxIrf + 1;
    idxIrt = idxIif + 1;      idxIit = idxIrt + 1;

    if(fabs(feq[idxIrf]) > 1e-8 || fabs(feq[idxIif]) > 1e-8 || fabs(feq[idxIrt]) > 1e-8 || fabs(feq[idxIit]) > 1e-8) {
      PetscPrintf(PETSC_COMM_SELF,"Line[%d] has equality constraint mismatch\n",i);
    }
  }

  for(i=0; i < ps->nbus; i++) {
    PetscInt idxIr = feqidxarr[ps->nline+i][0], idxIi = idxIr + 1;
    ierr = ComputeBusEqualityConstraints(busparamsarr[i],xidxarr[ps->nline+i],feqidxarr[ps->nline+i],x,feq);

    if(fabs(feq[idxIr]) > 1e-8 || fabs(feq[idxIi]) > 1e-8) {
      PetscPrintf(PETSC_COMM_SELF,"Bus[%d] has equality constraint mismatch\n",i);
    }

  }
  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(Feq,&feq);CHKERRQ(ierr);

  /* Check if equality constraints calculation is correct */
  double normF;
  ierr = VecNorm(Feq,NORM_INFINITY,&normF);CHKERRQ(ierr);
    if (normF < 1e-8) {
    ierr = PetscPrintf(PETSC_COMM_SELF,"Equality constraints calculation correct\n");CHKERRQ(ierr);
  } else {
    ierr = PetscPrintf(PETSC_COMM_SELF,"Equality constraints calculaltion incorrect\n");CHKERRQ(ierr);
  }

  //  ierr = VecView(F,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ComputeGlobalInequalityConstraints(OPFLOW opflow,double **branchparamsarr,double **busparamsarr,int **xidxarr, int **fineqidxarr,Vec X, Vec Fineq)
{
  PetscErrorCode ierr;
  PS             ps=opflow->ps;
  PetscInt       i;
  const double   *x;
  double         *fineq;

  PetscFunctionBegin;
  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecGetArray(Fineq,&fineq);CHKERRQ(ierr);

  for(i=0; i < ps->nbus; i++) {
    ierr = ComputeBusInequalityConstraints(busparamsarr[i],xidxarr[ps->nline+i],fineqidxarr[i],x,fineq);
  }

  for(i=0; i < ps->nline; i++) {
    ierr = ComputeBranchInequalityConstraints(branchparamsarr[i],xidxarr[i],fineqidxarr[ps->nbus+i],x,fineq);
  }

  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(Fineq,&fineq);CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_SELF,"Finished calculating inequality constraints\n");CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


PetscErrorCode CreateParamArray(OPFLOW opflow,double **branchparamsarr,double **busparamsarr,int **xidxarr,int **feqidxarr,int **fineqidxarr)
{
  PS ps=opflow->ps;
  PetscInt i,arr_size;
  PSBUS    bus;
  PetscInt nconnlines;
  const PSLINE *connlines;
  PSLINE   line;
  PSGEN    gen;
  PSLOAD   load;
  PetscErrorCode ierr;
  double   *busparams;
  double   *branchparams;
  int      *xidx;
  PetscInt locx,locxf,locxt;
  PSBUS    busf,bust;
  PetscInt j;
  const PSBUS *connbuses;
  PetscInt neqcons=0,nineqcons=0;
  PetscInt *feqidx;
  PetscInt *fineqidx;

  PetscFunctionBegin;

  for(i=0; i < ps->nline; i++) {
    line = &ps->line[i];
    /* size of branchparams[i] */
    arr_size = 10;

    /* Create the parameter array for branch,xidx */
    /* xidx containts the starting locations for variables for each branch in the X vector */
    /* feqidx contains the starting locations for equations for each branch in the equality constraints vector */
    ierr = PetscCalloc1(arr_size,&branchparamsarr[i]);CHKERRQ(ierr);
    ierr = PetscCalloc1(3,&xidxarr[i]);CHKERRQ(ierr);
    ierr = PetscCalloc1(1,&feqidxarr[i]);CHKERRQ(ierr);
    ierr = PetscCalloc1(1,&fineqidxarr[ps->nbus+i]);CHKERRQ(ierr);
    
    branchparams = branchparamsarr[i];

    if(!line->status) continue;

    branchparams[0] = (double)line->status;
    branchparams[1] = line->yff[0];
    branchparams[2] = line->yff[1];
    branchparams[3] = line->yft[0];
    branchparams[4] = line->yft[1];
    branchparams[5] = line->ytf[0];
    branchparams[6] = line->ytf[1];
    branchparams[7] = line->ytt[0];
    branchparams[8] = line->ytt[1];
    branchparams[9] = line->rateA;

    ierr = PSLINEGetVariableLocation(line,&locx);CHKERRQ(ierr);
    
    xidx = xidxarr[i];
    /* Starting location for variables for this branch */
    xidx[0] = locx;

    /* Starting location for the equality constraints for this branch */
    feqidx = feqidxarr[i];
    feqidx[0] = neqcons;
    neqcons += 4; /* Each branch has 4 equality constraints */

    /* Starting location for the inequaity constraints for this branch */
    if(line->rateA < 1e5) {
      /* If line rating > 1e5, then the line has infinite capacity and the
	 inequality constraints for the line can be excluded
      */
      fineqidx = fineqidxarr[ps->nbus+i];
      fineqidx[0] = ps->nbus + nineqcons; /* Branch inequality constraints are ordered after the bus inequality constraints. Hence, ps->nbus + is used here */
      nineqcons += 2; /* Each branch has 2 inequality constraints */
    }

    ierr = PSLINEGetConnectedBuses(line,&connbuses);CHKERRQ(ierr);
    busf = connbuses[0];
    bust = connbuses[1];

    ierr = PSBUSGetVariableLocation(busf,&locxf);CHKERRQ(ierr);
    ierr = PSBUSGetVariableLocation(bust,&locxt);CHKERRQ(ierr);

    /* Starting locations for from bus voltages */
    xidx[1] = locxf;
    /* Starting locations for to bus voltages */
    xidx[2] = locxt;
  }

  nineqcons = 0;
  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    arr_size = 0;
    /* Calculate size of busparams[i] */
    arr_size += 2; /* MVAbase,bus type */
    arr_size += 3; /* Number of generators,loads,lines */
    arr_size += 2; /* shunt (gl, bl) for each bus */
    arr_size += 6*bus->ngen; /* status,Pg, Qg alpha,beta,gamma for each gen */
    arr_size += 3*bus->nload; /* status,Pd, Qd for each load */

    //    ierr = PetscPrintf(PETSC_COMM_SELF,"Bus [%d]: param size = %d\n",i,arr_size);CHKERRQ(ierr);

    ierr = PSBUSGetSupportingLines(bus,&nconnlines,&connlines);CHKERRQ(ierr);
    arr_size += nconnlines;
    /* Create the parameter array for bus,xidx */
    /* xidx contains the starting locations for each bus in the X vector */
    ierr = PetscCalloc1(arr_size,&busparamsarr[i]);CHKERRQ(ierr);
    ierr = PetscCalloc1(1+nconnlines,&xidxarr[ps->nline+i]);CHKERRQ(ierr);

    ierr = PetscCalloc1(1,&feqidxarr[ps->nline+i]);CHKERRQ(ierr);
    ierr = PetscCalloc1(1,&fineqidxarr[i]);CHKERRQ(ierr);

    busparams = busparamsarr[i];
    xidx      = xidxarr[ps->nline+i];

    feqidx    = feqidxarr[ps->nline+i];
    feqidx[0] = neqcons;
    neqcons  += 2; /* 2 equality constraints at each bus */
    if(bus->ide == REF_BUS) neqcons += 1; /* Additional equality constraint for ref. bus */

    fineqidx = fineqidxarr[i];
    fineqidx[0] = nineqcons;
    nineqcons++;

    busparams[0] = ps->MVAbase;

    busparams++;

    busparams[0] = bus->ide;

    busparams++;

    /* # of gens, loads, and lines */
    busparams[0] = (double)bus->ngen;
    busparams[1] = (double)bus->nload;
    busparams[2] = (double)nconnlines;

    busparams += 3;

    /* shunt */
    busparams[0] = bus->gl;
    busparams[1] = bus->bl;

    busparams += 2;

    /* status, Pg, Qg for each gen */
    for(j=0; j < bus->ngen; j++) {
      ierr = PSBUSGetGen(bus,j,&gen);CHKERRQ(ierr);

      busparams[0] = (double)gen->status;
      busparams[1] = gen->pg;
      busparams[2] = gen->qg;
      busparams[3] = gen->cost_alpha;
      busparams[4] = gen->cost_beta;
      busparams[5] = gen->cost_gamma;

      busparams += 6;
    }

    /* status, Pd, Qd for each load */
    for(j=0; j < bus->nload; j++) {
      ierr = PSBUSGetLoad(bus,j,&load);CHKERRQ(ierr);

      busparams[0] = (double)load->status;
      busparams[1] = load->pl;
      busparams[2] = load->ql;

      busparams += 3;
    }
    
    PetscInt locx;
    ierr = PSBUSGetVariableLocation(bus,&locx);CHKERRQ(ierr);
    xidx[0] = locx;

    xidx++;
    
    for(j=0; j < nconnlines; j++) {
      line = connlines[j];
      if(!line->status) continue;

      ierr = PSLINEGetVariableLocation(line,&locx);CHKERRQ(ierr);
      ierr = PSLINEGetConnectedBuses(line,&connbuses);CHKERRQ(ierr);
      busf = connbuses[0];
      bust = connbuses[1];
      if(bus == busf) xidx[0] = locx; // From bus current injection
      else xidx[0] = locx+2; // To bus current injection

      xidx++;

      busparams[0] = line->status;

      busparams++;
    }
  }

  PetscFunctionReturn(0);
}

int main(int argc,char **argv)
{
  PetscErrorCode    ierr;
  OPFLOW            opflow;
  PS                ps;
  char              file[PETSC_MAX_PATH_LEN];
  PetscBool         flg;
  PetscLogStage     read,setup,solve,cuda_proto;
  PetscScalar       obj = 0.0; /* Objective function */
  /* array of arrays to hold parameters needed by each bus to
     compute its residual */
  double            **busparams;
  /* array of arrays to hold parameters needed by each branch to compute its
     residual */
  double            **branchparams;
  /* arrays to hold x locations needed by each bus and branch to access the variables 
     it needs to compute its residual */
  int               **xidx;
  /* arrays to hold eq. constraint locations where each bus and branch needs to insert its residuals */
  int               **feqidx;
  /* arrays to hold ineq. constraint locations where each bus and branch needs to insert its residuals */
  int               **fineqidx;
  PetscInt          i;
  char options_pathname[200] = EXAGO_OPTIONS_DIR;
  char* filename = "/opflowoptions";
  printf("%s\n", options_pathname);
  printf("%s\n", filename);
  strcat(options_pathname, filename);
  printf("%s\n", options_pathname);

  //  PetscInitialize(&argc,&argv,options_pathname,help);
  PetscInitialize(&argc,&argv,NULL,help);

  ierr = PetscLogStageRegister("ReadData",&read);CHKERRQ(ierr);
  ierr = PetscLogStageRegister("SetUp",&setup);CHKERRQ(ierr);
  ierr = PetscLogStageRegister("Solve",&solve);CHKERRQ(ierr);
  ierr = PetscLogStageRegister("CudaProto",&cuda_proto);CHKERRQ(ierr);

  /* Create OPFLOW object */
  ierr = OPFLOWCreate(PETSC_COMM_WORLD,&opflow);CHKERRQ(ierr);

  ierr = PetscLogStagePush(read);CHKERRQ(ierr);
  /* Get network data file from command line */
  ierr = PetscOptionsGetString(NULL,NULL,"-netfile",file,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
  /* Read Network Data file */
  if(flg) {
    ierr = OPFLOWReadMatPowerData(opflow,file);CHKERRQ(ierr);
  } else {
    ierr = OPFLOWReadMatPowerData(opflow,"datafiles/case9mod.m");CHKERRQ(ierr);
  }
  /* End of First Stage and Start of Second */
  PetscLogStagePop();

  ierr = PetscLogStagePush(setup);CHKERRQ(ierr);

  ierr = OPFLOWSetModel(opflow,"CURRENT_BALANCE_CARTESIAN2");CHKERRQ(ierr);

  /* Set up */
  ierr = OPFLOWSetUp(opflow);CHKERRQ(ierr);
  PetscLogStagePop();

  ierr = PetscLogStagePush(solve);CHKERRQ(ierr);
  /* Solve */
  ierr = OPFLOWSolve(opflow);CHKERRQ(ierr);
  PetscLogStagePop();

  ierr = PetscPrintf(PETSC_COMM_SELF,"\n\n**********************************\n");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF,"Begin testing GPU C prototype code\n");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF,"n**********************************\n");CHKERRQ(ierr);


  ps = opflow->ps;

  ierr = PetscLogStagePush(cuda_proto);CHKERRQ(ierr);

  /* The following array of arrays hold the data needed
     for the computations. 
  */
  /* Bus paramaters */
  ierr = PetscCalloc1(ps->nbus,&busparams);CHKERRQ(ierr);
  /* Branch parameters */
  ierr = PetscCalloc1(ps->nline,&branchparams);CHKERRQ(ierr);
  /* Indices for accessing variables */
  ierr = PetscCalloc1(ps->nbus+ps->nline,&xidx);CHKERRQ(ierr);
  /* Indices for setting equality constraints */
  ierr = PetscCalloc1(ps->nbus+ps->nline,&feqidx);CHKERRQ(ierr);
  /* Indices for setting inequaity constraints */
  ierr = PetscCalloc1(ps->nbus+ps->nline,&fineqidx);CHKERRQ(ierr);

  /* Set the data in the array of arrays */
  ierr = CreateParamArray(opflow,branchparams,busparams,xidx,feqidx,fineqidx);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF,"Finished Creating data arrays\n");CHKERRQ(ierr);

  /* Objective function */
  ierr = ComputeGlobalObjective(opflow,busparams,xidx,&obj,opflow->X);

  /* Equality constraints */
  ierr = ComputeGlobalEqualityConstraints(opflow,branchparams,busparams,xidx,feqidx,opflow->X,opflow->Ge);CHKERRQ(ierr);

  /* Inequality constraints */
  ierr = ComputeGlobalInequalityConstraints(opflow,branchparams,busparams,xidx,fineqidx,opflow->X,opflow->Gi);CHKERRQ(ierr);
  
  /* Destroy OPFLOW object */
  ierr = OPFLOWDestroy(&opflow);CHKERRQ(ierr);

  for(i=0; i < ps->nbus; i++) {
    ierr = PetscFree(busparams[i]);CHKERRQ(ierr);
  }
  ierr = PetscFree(busparams);CHKERRQ(ierr);
  for(i=0; i < ps->nline; i++) {
    ierr = PetscFree(branchparams[i]);CHKERRQ(ierr);
  }
  ierr = PetscFree(branchparams);CHKERRQ(ierr);
    for(i=0; i < ps->nline+ps->nbus; i++) {
    ierr = PetscFree(xidx[i]);CHKERRQ(ierr);
  }
  ierr = PetscFree(xidx);CHKERRQ(ierr);

  for(i=0; i < ps->nline+ps->nbus; i++) {
    ierr = PetscFree(feqidx[i]);CHKERRQ(ierr);
  }
  ierr = PetscFree(feqidx);CHKERRQ(ierr);

  for(i=0; i < ps->nline+ps->nbus; i++) {
    ierr = PetscFree(fineqidx[i]);CHKERRQ(ierr);
  }
  ierr = PetscFree(fineqidx);CHKERRQ(ierr);

  PetscFinalize();
  return 0;
}
