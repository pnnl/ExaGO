#include <private/opflowimpl.h>
#include <private/sopflowimpl.h>

extern void clean2Char(char *);
extern char** blankTokenizer(const char *str, int *numtok, int maxtokens, int maxchar);
/*
  SOPFLOWSetScenarioData - Sets the scenario data

  Input Parameter
+  sopflow - The SOPFLOW object
.  scenfileformat - the scenario file format
.  scenunctype    - type of uncertainty
-  scenfile - The name of the scenario list file

*/
PetscErrorCode SOPFLOWSetScenarioData(SOPFLOW sopflow,ScenarioFileInputFormat scenfileformat,ScenarioUncertaintyType scenunctype,const char scenfile[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = PetscMemcpy(sopflow->scenfile,scenfile,100*sizeof(char));CHKERRQ(ierr);

  sopflow->scenfileformat = scenfileformat;
  sopflow->scenfileset = PETSC_TRUE;
  sopflow->scenunctype = scenunctype;
  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWReadScenarioData_Wind(SOPFLOW sopflow,const char windgenprofile[])
{
  PetscErrorCode ierr;
  FILE           *fp;
  char           line[MAXLINE];
  char           *out;
  OPFLOW         opflow;
  PS             ps;
  PetscInt       ngen=sopflow->opflows[0]->ps->ngen,*windgenbus,nw=0;
  char           *tok,*tok2;
  char           sep[] = ",",sep2[] = "_";
  PSGEN          gen;
  PetscReal      pg;
  int            scen_num=0;
  int            genid;
  char           windgenid[100][3];

  PetscFunctionBegin;

  fp = fopen(windgenprofile,"r");
  if (fp == NULL) {
    SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Cannot open wind generation profile file %s",windgenprofile);CHKERRQ(ierr);
  }

  ierr = PetscMalloc1(ngen,&windgenbus);CHKERRQ(ierr);

  /* First line -- has the bus numbers */
  out = fgets(line,MAXLINE,fp);

  /* Parse wind generator numbers */
  tok = strtok(line,sep);
  tok = strtok(NULL,sep); /* Skip first token */
  tok = strtok(NULL,sep); /* Skip second token */
  while(tok != NULL) {
    /* Parse generator info */
    tok2 = strsep(&tok,sep2);
    sscanf(tok2,"%d",&windgenbus[nw]);
    tok2 = strsep(&tok,sep2);
    tok2 = strsep(&tok,sep2); /* Skip string "Wind" */
    sscanf(tok2,"%d",&genid);
    if(nw == 100) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Exceeded max. values set for parsing wind generators\n");
    snprintf(windgenid[nw],3,"%-2d",genid);

    nw++;
    tok = strtok(NULL,sep);
  }

  while((out = fgets(line,MAXLINE,fp)) != NULL) {
    if(strcmp(line,"\r\n") == 0 || strcmp(line,"\n") == 0) {
      continue; /* Skip blank lines */
    }

    tok = strtok(line,sep);
    tok = strtok(NULL,sep); /* Skip first token */
    sscanf(tok,"%d",&scen_num); /* Scenario number */
    opflow = sopflow->opflows[scen_num]; // Only good for serial 
    ps     = opflow->ps;
    tok = strtok(NULL,sep);
    nw = 0;
    while(tok != NULL) {
      ierr = PSGetGen(ps,windgenbus[nw],windgenid[nw],&gen);CHKERRQ(ierr);
      sscanf(tok,"%lf",&pg);
      gen->pg = gen->pt = pg/ps->MVAbase; /* Set real power generation. Note that Pg upper limit is also set to Pg */
      nw++;
      tok = strtok(NULL,sep);
    }

    if(scen_num == sopflow->Ns-1) break;
  }

  ierr = PetscFree(windgenbus);CHKERRQ(ierr);

  fclose(fp);
  PetscFunctionReturn(0);
}
  
/*
  SOPFLOWReadScenarioData - Reads the scenario data file

  Input Parameters
+ sopflow - the sopflow object
. scenfileformat - the scenario file format (NATIVE or PSSE)
- scenfile - the scenario file name

*/
PetscErrorCode SOPFLOWReadScenarioData(SOPFLOW sopflow,ScenarioFileInputFormat scenfileformat,const char scenfile[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if(sopflow->scenunctype == WIND) {
    ierr = SOPFLOWReadScenarioData_Wind(sopflow,scenfile);CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

/* SOPFLOWGetNumScenarios_Native - Gets the number of scenarios from the scenario file
 */
PetscErrorCode SOPFLOWGetNumScenarios_Native(const char scenfile[],PetscInt *Ns)
{
  PetscErrorCode ierr;
  FILE           *fp;
  char           line[MAXLINE];
  char           *out;
  char           *tok,*tok2;
  char           sep[] = ",",sep2[] = "_";
  PetscInt       t=0;
  int            scen_num,ns=0;

  PetscFunctionBegin;

  fp = fopen(scenfile,"r");
  if (fp == NULL) {
    SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Cannot open wind generation profile file %s",scenfile);CHKERRQ(ierr);
  }

  /* First line -- has the bus numbers */
  out = fgets(line,MAXLINE,fp);

  while((out = fgets(line,MAXLINE,fp)) != NULL) {
    if(strcmp(line,"\r\n") == 0 || strcmp(line,"\n") == 0) {
      continue; /* Skip blank lines */
    }

    tok = strtok(line,sep);
    tok = strtok(NULL,sep); /* Skip first token */
    sscanf(tok,"%d",&scen_num); /* Scenario number */
    if(ns < scen_num) ns = scen_num;
    else break;
  }

  fclose(fp);
  *Ns = ns;
  PetscFunctionReturn(0);
}

/*
  SOPFLOWReadScenarioData - Gets the number of scenarios from the scenario data file

  Input Parameters
+ sopflow - the sopflow object
. scenfileformat - the scenario file format (NATIVE or PSSE)
. scenfile - the scenario file name
- Ns - number of scenarios given in the scenario file

*/
PetscErrorCode SOPFLOWGetNumScenarios(SOPFLOW sopflow,ScenarioFileInputFormat scenfileformat,const char scenfile[],PetscInt *Ns)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if(scenfileformat == SOPFLOW_NATIVE) {
    ierr = SOPFLOWGetNumScenarios_Native(scenfile,Ns);CHKERRQ(ierr);
  }
  else *Ns = 1;

  PetscFunctionReturn(0);
}

