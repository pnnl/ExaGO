static char help[] = "Prototype version of OPFLOW code for CUDA/GPU/Accelerator implementation\n\n";

#include <scopflow_config.h>
#include <private/psimpl.h>
#include <private/opflowimpl.h>

/* The parameters are packed in each busparams[i] as follows
busparams[i] = {# of gens
                # of loads
		# of lines  
		# shunt (gl,bl for each bus)
		# status,Pg,Qg,alpha,beta,gamma for each gen
		# Pd,Qd for each load
		# status,Gff,Bff,Gft,Bft,Gtf,Btf,Gtt,Btt for each line

*/

PetscLogEvent Objfun_logger;
PetscLogEvent Gradfun_logger;
PetscLogEvent Eqconsfunc_logger;
PetscLogEvent Ineqconsfunc_logger;

/* Compute the objective value contribution (generation cost) at the bus */
  /* ComputeBusObjective should loop over the generators at the bus and add the generation cost. See OPFLOWComputeObjective_PBCAR. Note that OPFLOWComputeObjective_PBCAR loops over all buses. This is being already done in ComputeGlobalObjective so this SHOULD NOT be done in this function. This function should ONLY calculate the generation cost at the bus. All the data that is needed to compute the generation cost is in the busparams array. See CreateParamsArray to see the layout of the data. A few
  other points:
  -- Assume MVAbase = 100.0
  -- Ignore PSBUSISghost()
  -- Ignore costs for load loss and power imbalance
  */

int ComputeBusObjective(double *busparams,int *xidx,const double *x,double *obj)
{
   
   double MVAbase=100;
   int ngen = (int)busparams[0];
   int gstatus;
   int i;
   int locxgen = xidx[0]+2;
   double alpha, beta, gamma;
   double pg;

   busparams += 5;

   for(i=0; i < ngen; i++){

     gstatus = (int) busparams[0];
     if(!gstatus) continue;
	   pg = x[locxgen];
	   alpha = busparams[3];
	   beta  = busparams[4];
	   gamma = busparams[5];
     
     *obj += alpha*pg*pg*MVAbase*MVAbase + beta*pg*MVAbase + gamma;
     
     busparams +=6;
     locxgen   +=2;
   }

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
    ierr = ComputeBusObjective(busparamsarr[i],xidxarr[i],x,&busobj);
    *obj += busobj; /* Add the bus objective value */
  }

  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_SELF,"The objective function value is %f\n",*obj);

  PetscFunctionReturn(0);
}


int ComputeBusEqualityConstraints(double *busparams,int *xidx,int *feqidx,const double *x, double *feq)
{

  int ngen = (int)busparams[0];
  int nload = (int)busparams[1];
  int nlines = (int)busparams[2];
  double gl  = busparams[3];
  double bl  = busparams[4];
  int    idxp  = feqidx[0],idxq=feqidx[1];
  int    idxVr = xidx[0],idxVi=xidx[1];
  double Vr    = x[idxVr];
  double Vi    = x[idxVi];
  double Vm2   = Vr*Vr + Vi*Vi; 
  int    i;
  int gstatus,lstatus,brstatus;
  int    locxgen = xidx[0] + 2;

  double Pg,Qg,Pd,Qd;
  double Gff, Bff, Gft, Bft, Gtf, Btf, Gtt, Btt;
  double Vrf, Vif, Vrt, Vit;
  double Pline, Qline;

  busparams +=5;
  xidx      +=2;

  
  /* shunt contribution */
  feq[idxp] = gl*Vm2;
  feq[idxq] = -bl*Vm2;

  
  /* Generator Contribution */
  for(i=0; i< ngen; i++){
    gstatus = (int) busparams[0];
    Pg      = x[locxgen];
    Qg      = x[locxgen+1];

    feq[idxp] -= gstatus*Pg;
    feq[idxq] -= gstatus*Qg;
    
    busparams +=6;
    locxgen +=2;

  }

  /* Load Contribution */
  for(i=0;i<nload;i++){
    lstatus = (int) busparams[0];
    Pd      = busparams[1];
    Qd      = busparams[2];

    feq[idxp] += lstatus*Pd;
    feq[idxq] += lstatus*Qd;
    
    busparams += 3;
    
  }

  /* Line Contribution */

  for(i=0; i<nlines;i++){
   
    brstatus = busparams[0];
    Gff = busparams[1];
    Bff = busparams[2];
    Gft = busparams[3];
    Bft = busparams[4];
    Gtf = busparams[5];
    Btf = busparams[6];
    Gtt = busparams[7];
    Btt = busparams[8];
   
    busparams += 9;

      Vrf  = x[xidx[0]];
      Vif  = x[xidx[1]];
      Vrt  = x[xidx[2]];
      Vit  = x[xidx[3]];
      
      if (xidx[0] == idxVr) {
	        Pline =  Gff*(Vrf*Vrf + Vif*Vif) + Vrf*(Gft*Vrt - Bft*Vit) + Vif*(Bft*Vrt + Gft*Vit);
	        Qline = -Bff*(Vrf*Vrf + Vif*Vif) + Vif*(Gft*Vrt - Bft*Vit) - Vrf*(Bft*Vrt + Gft*Vit);

      } else {
        	Pline =  Gtt*(Vrt*Vrt + Vit*Vit) + Vrt*(Gtf*Vrf - Btf*Vif) + Vit*(Btf*Vrf + Gtf*Vif);
	        Qline = -Btt*(Vrt*Vrt + Vit*Vit) + Vit*(Gtf*Vrf - Btf*Vif) - Vrt*(Btf*Vrf + Gtf*Vif);
      }
      xidx+=4;

    feq[idxp] += brstatus*Pline;
    feq[idxq] += brstatus*Qline;
    
  }

  return 0;

}


PetscErrorCode ComputeGlobalEqualityConstraints(OPFLOW opflow,double **busparamsarr,int **xidxarr, int **feqidxarr,Vec X, Vec Feq)
{

  PetscErrorCode ierr;
  PS             ps=opflow->ps;
  PetscInt       i;
  const double   *x;
  double         *feq;
  PSBUS          bus;
  PetscInt       idxp,idxq;

  PetscFunctionBegin;
  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecGetArray(Feq,&feq);CHKERRQ(ierr);

  for(i=0; i < ps->nbus; i++) {

    bus = &ps->bus[i];

    ierr = ComputeBusEqualityConstraints(busparamsarr[i],xidxarr[i],feqidxarr[i],x,feq);
    idxp = feqidxarr[i][0]; idxq = feqidxarr[i][1];
    
    if(fabs(feq[idxp]) > 1e-8 || fabs(feq[idxq]) > 1e-8) {
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

int ComputeBusInequalityConstraints(double *busparams, int *xidx, int *fineqidx, const double *x, double *fineq)
{
  double Vr    = x[xidx[0]];
  double Vi    = x[xidx[1]];
  double Vm2   = Vr*Vr + Vi*Vi; 
  
  xidx +=2;  
  fineq[fineqidx[0]] = Vm2;
  
  return 0;
}
PetscErrorCode ComputePowerFlowInequalityConstraints(double *branchparams,int *xidx,int *fineqidx,const double *x, double *fineq)
{

  int brstatus;
  int i;

  double Gff, Bff, Gft, Bft, Gtf, Btf, Gtt, Btt;
  double Vrf, Vif, Vrt, Vit;
  double Pf, Qf,Pt,Qt;
  double Sf2,St2;

  // xidx +=2;
    
  brstatus = branchparams[0];
    
    if(!brstatus) return 0;
    
    Gff = branchparams[1];
    Bff = branchparams[2];
    Gft = branchparams[3];
    Bft = branchparams[4];
    Gtf = branchparams[5];
    Btf = branchparams[6];
    Gtt = branchparams[7];
    Btt = branchparams[8];

    
    branchparams += 9;

    Vrf  = x[xidx[0]];
    Vif  = x[xidx[1]];
    Vrt  = x[xidx[2]];
    Vit  = x[xidx[3]];

  
    Pf =  Gff*(Vrf*Vrf + Vif*Vif) + Vrf*(Gft*Vrt - Bft*Vit) + Vif*(Bft*Vrt + Gft*Vit);
    Qf = -Bff*(Vrf*Vrf + Vif*Vif) + Vif*(Gft*Vrt - Bft*Vit) - Vrf*(Bft*Vrt + Gft*Vit);
    
    Pt =  Gtt*(Vrt*Vrt + Vit*Vit) + Vrt*(Gtf*Vrf - Btf*Vif) + Vit*(Btf*Vrf + Gtf*Vif);
    Qt = -Btt*(Vrt*Vrt + Vit*Vit) + Vit*(Gtf*Vrf - Btf*Vif) - Vrt*(Btf*Vrf + Gtf*Vif);

   
    Sf2 = Pf*Pf + Qf*Qf;
    St2 = Pt*Pt + Qt*Qt;
    
    fineq[fineqidx[0]] = Sf2;
    fineq[fineqidx[0]+1] = St2;

    xidx +=4;

    
    return 0;

}


PetscErrorCode ComputeGlobalInequalityConstraints(OPFLOW opflow,double **busparamsarr,double **branchparamsarr, int **xidxarr, int **fineqidxarr,Vec X, Vec Fineq)
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
    ierr = ComputeBusInequalityConstraints(busparamsarr[i],xidxarr[i],fineqidxarr[i],x,fineq);
  }

  for(i=0; i < ps->nline; i++){
    ierr = ComputePowerFlowInequalityConstraints(branchparamsarr[i],xidxarr[ps->nbus+i],fineqidxarr[ps->nbus+i],x,fineq);
  }

  // for(i=0; i < ps->nbus+2*ps->nline; i++){
  //   printf("fineq %f \n",fineq[i]);
  // }

  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(Fineq,&fineq);CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_SELF,"Finished calculating inequality constraints\n");CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


PetscErrorCode CreateParamArray(OPFLOW opflow,double **busparamsarr, double **branchparamsarr, int **xidxarr, int **feqidxarr, int **fineqidxarr)
{
  PS ps=opflow->ps;
  PetscInt i, j, arr_size;
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
  int      *feqidx;
  int      *fineqidx;
  PetscInt neqcons=0, bineqcons=0, brineqcons=0;

  PSBUS busf,bust;
  PetscInt locx, locxf,locxt;
  const PSBUS *connbuses;

  PetscFunctionBegin;
  
  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    arr_size = 0;
    /* Calculate size of busparams[i] */
    arr_size += 3; /* Number of generators,loads,lines */
    arr_size += 2; /* shunt (gl, bl) for each bus */
    arr_size += 6*bus->ngen; /* status,Pg, Qg alpha,beta,gamma for each gen */
    arr_size += 3*bus->nload; /* status,Pd, Qd for each load */

    ierr = PSBUSGetSupportingLines(bus,&nconnlines,&connlines);CHKERRQ(ierr);

    arr_size += 9*nconnlines; /* status,Gff,Bff,Gft,Bft,Gtf,Btf,Gtt,Btt for each line */

    //    ierr = PetscPrintf(PETSC_COMM_SELF,"Bus [%d]: param size = %d\n",i,arr_size);CHKERRQ(ierr);

    /* Create the parameter array for bus,xidx, feqidx */
    /* xidx contains the starting locations for each bus in the X vector */
    ierr = PetscCalloc1(arr_size,&busparamsarr[i]);CHKERRQ(ierr);
    ierr = PetscCalloc1(2+4*nconnlines,&xidxarr[i]);CHKERRQ(ierr);
    ierr = PetscCalloc1(1,&feqidxarr[i]);CHKERRQ(ierr);
    ierr = PetscCalloc1(1,&fineqidxarr[i]);CHKERRQ(ierr);

    busparams = busparamsarr[i];
    xidx      = xidxarr[i];
    feqidx    = feqidxarr[i];
    feqidx[0] = neqcons;
    neqcons  += 2;
    fineqidx  = fineqidxarr[i];
    fineqidx[0] = bineqcons;
    bineqcons +=1;


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
    xidx[0] = feqidx[0]=locx;
    xidx[1] = feqidx[1]=locx+1;
    
    
    xidx += 2;
        
    for(j=0; j < nconnlines; j++) {
 
      line = connlines[j];

      ierr = PSLINEGetConnectedBuses(line,&connbuses);CHKERRQ(ierr);
      busf = connbuses[0];
      bust = connbuses[1];
      
      ierr = PSBUSGetVariableLocation(busf,&locxf);CHKERRQ(ierr);
      ierr = PSBUSGetVariableLocation(bust,&locxt);CHKERRQ(ierr);
      
      xidx[0] = locxf;
      xidx[1] = locxf+1;
      xidx[2] = locxt;
      xidx[3] = locxt+1;

           
      xidx += 4;

      busparams[0] = (double)line->status;
      busparams[1] = line->yff[0];
      busparams[2] = line->yff[1];
      busparams[3] = line->yft[0];
      busparams[4] = line->yft[1];
      busparams[5] = line->ytf[0];
      busparams[6] = line->ytf[1];
      busparams[7] = line->ytt[0];
      busparams[8] = line->ytt[1];

      busparams += 9;

    }
    
  }
  brineqcons = 0;

  for(i=0; i < ps->nline; i++) {
      
      line = &ps->line[i];
      /* size of branchparams[i] */
      arr_size = 10;

      /* Create the parameter array for branch,xidx */
      /* xidx containts the starting locations for variables for each branch in the X vector */
      /* feqidx contains the starting locations for equations for each branch in the equality constraints vector */
      ierr = PetscCalloc1(arr_size,&branchparamsarr[i]);CHKERRQ(ierr);
      ierr = PetscCalloc1(4,&xidxarr[ps->nbus+i]);CHKERRQ(ierr);
      ierr = PetscCalloc1(2,&fineqidxarr[ps->nbus+i]);CHKERRQ(ierr);
      
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
      
      // ierr = PSLINEGetVariableLocation(line, &locx);CHKERRQ(ierr);
      xidx = xidxarr[ps->nbus+i];
      // xidx[0] = locx;

      ierr = PSLINEGetConnectedBuses(line,&connbuses);CHKERRQ(ierr);
      busf = connbuses[0];
      bust = connbuses[1];

      ierr = PSBUSGetVariableLocation(busf,&locxf);CHKERRQ(ierr);
      ierr = PSBUSGetVariableLocation(bust,&locxt);CHKERRQ(ierr);

      /* Starting locations for from bus voltages */
      xidx[0] = locxf;
      xidx[1] = locxf+1;
      
      /* Starting locations for to bus voltages */
      xidx[2] = locxt;
      xidx[3] = locxt+1;

      xidx += 4;
            
      fineqidx = fineqidxarr[ps->nbus+i];
      fineqidx[0] = ps->nbus + brineqcons;
      brineqcons +=2;

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
  double            **branchparams;
  /* arrays to hold x locations needed by each bus to access the variables 
     it needs to compute its residual */
  int               **xidx, **feqidx, **fineqidx;
  /* arrays to hold f locations where each bus needs to insert its residuals */
  PetscInt          i;
  Vec F;
  char options_pathname[200] = EXAGO_OPTIONS_DIR;
  char* filename = "/opflowoptions";
  printf("%s\n", options_pathname);
  printf("%s\n", filename);
  strcat(options_pathname, filename);
  printf("%s\n", options_pathname);

  PetscInitialize(&argc,&argv,options_pathname,help);

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
  /* Set up */
  ierr = OPFLOWSetUp(opflow);CHKERRQ(ierr);
  PetscLogStagePop();

  ierr = PetscLogStagePush(solve);CHKERRQ(ierr);
  /* Solve */
  ierr = OPFLOWSolve(opflow);CHKERRQ(ierr);
  PetscLogStagePop();
  
  ps = opflow->ps;

  ierr = PetscLogStagePush(cuda_proto);CHKERRQ(ierr);

  ierr = PetscCalloc1(ps->nbus,&busparams);CHKERRQ(ierr);
  ierr = PetscCalloc1(ps->nline,&branchparams);CHKERRQ(ierr);
  ierr = PetscCalloc1(ps->nbus+ps->nline,&xidx);CHKERRQ(ierr);
  ierr = PetscCalloc1(ps->nbus,&feqidx);CHKERRQ(ierr);
  ierr = PetscCalloc1(ps->nbus+ps->nline,&fineqidx);CHKERRQ(ierr);

  ierr = CreateParamArray(opflow,busparams,branchparams, xidx,feqidx, fineqidx);CHKERRQ(ierr);

  ierr = ComputeGlobalObjective(opflow,busparams,xidx,&obj,opflow->X);CHKERRQ(ierr);
    /* Equality constraints */
  ierr = ComputeGlobalEqualityConstraints(opflow,busparams,xidx,feqidx,opflow->X,opflow->Ge);CHKERRQ(ierr);
  /* Inequality constraints */
  ierr = ComputeGlobalInequalityConstraints(opflow,busparams,branchparams, xidx,fineqidx,opflow->X,opflow->Gi);CHKERRQ(ierr);

  PetscLogStagePop();

  PetscFinalize();

  return 0;
}
