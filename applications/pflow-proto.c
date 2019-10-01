static char help[] = "Prototype version of C code for CUDA/GPU/Accelerator implementation\n\n";

#include <private/psimpl.h>
#include <private/pflowimpl.h>

/* The parameters are packed in each busparams[i] as follows
busparams[i] = {# of gens
                # of loads
		# of lines
		# shunt (gl,bl for each bus)
		# status,Pg,Qg for each gen
		# Pd,Qd for each load
		# status,Gff,Bff,Gft,Bft,Gtf,Btf,Gtt,Btt for each line

*/
PetscErrorCode CreateParamArray(PFLOW pflow,double **busparamsarr)
{
  PS ps=pflow->ps;
  PetscInt i,arr_size;
  PSBUS    bus;
  PetscInt nconnlines;
  const PSLINE *connlines;
  PSLINE   line;
  PSGEN    gen;
  PSLOAD   load;
  PetscErrorCode ierr;
  double   *busparams;
  PetscInt j;

  PetscFunctionBegin;

  for(i=0; i < ps->nbus; i++) {
    bus = &ps->bus[i];
    arr_size = 0;
    /* Calculate size of busparams[i] */
    arr_size += 3; /* Number of generators,loads,lines */
    arr_size += 2; /* shunt (gl, bl) for each bus */
    arr_size += 3*bus->ngen; /* status,Pg, Qg for each gen */
    arr_size += 3*bus->nload; /* status,Pd, Qd for each load */

    ierr = PSBUSGetSupportingLines(bus,&nconnlines,&connlines);CHKERRQ(ierr);

    arr_size += 9*nconnlines; /* status,Gff,Bff,Gft,Bft,Gtf,Btf,Gtt,Btt for each line */

    ierr = PetscPrintf(PETSC_COMM_SELF,"Bus [%d]: param size = %d\n",i,arr_size);CHKERRQ(ierr);

    /* Create the parameter array for bus and pack it */
    ierr = PetscCalloc1(arr_size,&busparamsarr[i]);CHKERRQ(ierr);
    busparams = busparamsarr[i];

    /* # of gens, loads, and lines */
    busparams[0]   = (double)bus->ngen;
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

      busparams += 3;
    }

    /* status, Pd, Qd for each load */
    for(j=0; j < bus->nload; j++) {
      ierr = PSBUSGetLoad(bus,j,&load);CHKERRQ(ierr);

      busparams[0] = (double)load->status;
      busparams[1] = load->pl;
      busparams[2] = load->ql;

      busparams += 3;
    }
    
    for(j=0; j < nconnlines; j++) {
      line = connlines[j];
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

  PetscFunctionReturn(0);
}

int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  PFLOW             pflow;
  char              file[PETSC_MAX_PATH_LEN];
  PetscBool         flg;
  PetscLogStage     read,setup,solve;
  Vec               X; /* Global solution vector */
  Vec               F; /* Global residual vector */
  /* array of arrays to hold parameters needed by each bus to
     compute its residual */
  double            **busparams;
  
  PetscInitialize(&argc,&argv,"options/pflowoptions",help);
  
  ierr = PetscLogStageRegister("ReadData",&read);CHKERRQ(ierr);
  ierr = PetscLogStageRegister("SetUp",&setup);CHKERRQ(ierr);
  ierr = PetscLogStageRegister("Solve",&solve);CHKERRQ(ierr);

  /* Create PFLOW object */
  ierr = PFLOWCreate(PETSC_COMM_WORLD,&pflow);CHKERRQ(ierr);

  ierr = PetscLogStagePush(read);CHKERRQ(ierr);
  /* Get network data file from command line */
  ierr = PetscOptionsGetString(NULL,NULL,"-netfile",file,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
  /* Read Network Data file */
  if(flg) {
    ierr = PFLOWReadMatPowerData(pflow,file);CHKERRQ(ierr);
  } else {
    ierr = PFLOWReadMatPowerData(pflow,"datafiles/case9mod.m");CHKERRQ(ierr);
  }
  PetscLogStagePop();

  ierr = PetscLogStagePush(setup);CHKERRQ(ierr);
  /* Set up */
  ierr = PFLOWSetUp(pflow);CHKERRQ(ierr);
  PetscLogStagePop();

  ierr = PFLOWCreateGlobalVector(pflow,&X);CHKERRQ(ierr);
  ierr = VecDuplicate(X,&F);CHKERRQ(ierr);

  /* Set Initial Guess for X */
  ierr = PFLOWSetInitialGuess(pflow,X);CHKERRQ(ierr);

  PS ps = pflow->ps;
  ierr = PetscCalloc1(ps->nbus,&busparams);CHKERRQ(ierr);

  ierr = CreateParamArray(pflow,busparams);CHKERRQ(ierr);


  //  ierr = PetscLogStagePush(solve);CHKERRQ(ierr);
  /* Solve */
  //  ierr = PFLOWSolve(pflow);CHKERRQ(ierr);
  //  PetscLogStagePop();

  /* Update line flows, Pgen, Qgen, and other parameters */
  //  ierr = PFLOWPostSolve(pflow);CHKERRQ(ierr);

  for(PetscInt i=0; i < pflow->ps->nbus; i++) {
    ierr = PetscFree(busparams[i]);CHKERRQ(ierr);
  }
  ierr = PetscFree(busparams);CHKERRQ(ierr);

  ierr = VecDestroy(&X);CHKERRQ(ierr);
  ierr = VecDestroy(&F);CHKERRQ(ierr);
  /* Destroy PFLOW object */
  ierr = PFLOWDestroy(&pflow);CHKERRQ(ierr);

  PetscFinalize();
  return 0;
}
  
