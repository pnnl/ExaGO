#include <dyn.h>
#define MAXLINE 1000

PetscErrorCode ReadParameters_Private(const char paramfile[],PetscInt n,PetscInt genbus[],PetscInt genstat[],PetscScalar Pg[],PetscScalar Qg[])
{
  char     label[20];
  FILE     *fp;
  char      line[MAXLINE],*data;
  PetscInt i,offset;

  PetscFunctionBegin;
  if ((fp = fopen(paramfile,"r"))==NULL) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Cannot open file %s",paramfile);

  while ((fgets(line,MAXLINE,fp))!=NULL) {
    sscanf(line,"%s%n",label,&offset);
    data = line+offset;
    if (!strcmp(label,"GENBUS")) {
      for (i=0;i<n;i++) {
        sscanf(data," %d%n",&genbus[i],&offset);
        data += offset;
      }
    }
    if (!strcmp(label,"GENSTAT")) {
      for (i=0;i<n;i++) {
        sscanf(data," %d%n",&genstat[i],&offset);
        data += offset;
      }
    }
    if (!strcmp(label,"PG")) {
      for (i=0;i<n;i++) {
        sscanf(data," %lf%n",&Pg[i],&offset);
        data += offset;
      }
    }
    if (!strcmp(label,"QG")) {
      for (i=0;i<n;i++) {
        sscanf(data," %lf%n",&Qg[i],&offset);
        data += offset;
      }
    }
  }
  fclose(fp);
  PetscFunctionReturn(0);
}

int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  double         *costfcn;
  double         *sensitivities;
  PetscBool      single_costfcn=PETSC_FALSE;
  PetscBool      active_power_only=PETSC_FALSE,flg;
  PetscInt       i,j,ngen=-1,nparam=-1,flag,nbus=-1,ncostfcns;
  PetscBool      updatedispatch=PETSC_FALSE,monitor=PETSC_TRUE;
  char           netfile[PETSC_MAX_PATH_LEN],dyrfile[PETSC_MAX_PATH_LEN],eventfile[PETSC_MAX_PATH_LEN];
  char           paramfile[PETSC_MAX_PATH_LEN];

  /* 9 Bus */
  PetscInt genbus[3] = {1,2,3};
  PetscInt gennum[3] = {0,0,0};
  PetscInt genstat[3] = {1,1,1};
  PetscScalar Pg[3] = {46,269,135};
  PetscScalar Qg[3]= {0.0,0.0,0.0};

  //PetscScalar demand[9] = {0,0,0,0,200,200,0,200,0};  // not used

  /* 39 bus */
  /*
  PetscInt genbus[10] = {40,41,42,43,44,45,46,47,48,49};
  PetscInt gennum[10] = {0,0,0,0,0,0,0,0,0,0};
  PetscInt genstat[10] = {1,1,1,0,0,0,1,1,1,1};
  PetscScalar Pg[10] = {500.0,400.0,1000.0,0.0, 0.0, 0.0, 537.1, 1000.0, 1660.0, 1000.0};
  PetscScalar Qg[10]= {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  */
  PetscReal endtime = 1.0;
  PetscReal finaltime;
  PetscReal stepsize = 0.01;
  PetscInt cost_type = 1; // 0 freq 1 freq violation
  PetscInt print_vm = 0; // print voltages
  //PetscInt set_load = 0; // set load to demand (assuming one load per bus - zero demand if no load is present

  ierr = Initialize(argc,argv);CHKERRQ(ierr);

  ierr = PetscOptionsGetInt(NULL,NULL,"-ngen",&ngen,&flg);CHKERRQ(ierr);
  if(ngen == -1) SETERRQ(PETSC_COMM_SELF,0,"Please set the number of generators using -ngen <ngen> run time option\n");

  ierr = PetscOptionsGetInt(NULL,NULL,"-nbus",&nbus,&flg);CHKERRQ(ierr);
  if(nbus == -1) SETERRQ(PETSC_COMM_SELF,0,"Please set the number of buses using -nbus <nbus> run time option\n");

  nparam = 2*nbus + 2*ngen;

  ierr = PetscCalloc1(ngen,&costfcn);CHKERRQ(ierr);
  if (active_power_only) {
    ierr = PetscCalloc1(ngen*ngen,&sensitivities);CHKERRQ(ierr);
  } else {
    ierr = PetscCalloc1(ngen*nparam,&sensitivities);CHKERRQ(ierr);
  }
  //ierr = PetscCalloc1(ngen,&Pg);CHKERRQ(ierr);
  //ierr = PetscCalloc1(ngen,&Qg);CHKERRQ(ierr);

  ierr = PetscStrcpy(netfile,"datafiles/case9mod.m");CHKERRQ(ierr);
  /* Get network data file from command line */
  ierr = PetscOptionsGetString(NULL,NULL,"-netfile",netfile,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);

  ierr = PetscStrcpy(dyrfile,"datafiles/case9mod.dyr");CHKERRQ(ierr);
  /* Get dyr data file from command line */
  ierr = PetscOptionsGetString(NULL,NULL,"-dyrfile",dyrfile,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);

  ierr = PetscStrcpy(eventfile,"datafiles/case9mod.event");CHKERRQ(ierr);
  /* Get event data file from command line */
  ierr = PetscOptionsGetString(NULL,NULL,"-eventfile",eventfile,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);

  ierr = PetscOptionsGetString(NULL,NULL,"-paramfile",paramfile,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = ReadParameters_Private(paramfile,ngen,genbus,genstat,Pg,Qg);CHKERRQ(ierr);
    updatedispatch = PETSC_TRUE;
  }
  flag = DYNGetCostFunctionAndSensitivities(argc,argv,netfile,dyrfile,eventfile,
                                            costfcn,sensitivities,
                                            single_costfcn,active_power_only,updatedispatch,
                                            genbus,gennum,genstat,Pg,Qg,
                                            endtime,stepsize,print_vm,cost_type,59.5,60.5,monitor,&finaltime);CHKERRQ(ierr);

  if (cost_type == 1 && single_costfcn) ncostfcns = 1;
  else ncostfcns = ngen;
  if (active_power_only) {
    printf("Printing senstivities to active power only!\n");
    for (i=0; i<ncostfcns; i++) {
      printf("costfcn[%d] = %5.4g\n",i,costfcn[i]);
      for (j=0; j<ngen; j++) {
        printf("sens[%d %d] = %5.4g\n",i,j,sensitivities[ngen*i+j]);
      }
    }
  } else {
    for (i=0; i<ncostfcns; i++) {
      printf("costfcn[%d] = %5.4g\n",i,costfcn[i]);
      for (j=0; j<nparam; j++) {
        printf("sens[%d %d] = %5.4g\n",i,j,sensitivities[nparam*i+j]);
      }
    }
  }

  Finalize();
  return 0;
}

