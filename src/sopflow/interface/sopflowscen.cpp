#include <private/opflowimpl.h>
#include <private/tcopflowimpl.h>
#include <private/scopflowimpl.h>
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

/*
  SOPFLOWReadScenarioData_Wind - Reads the wind data 
  Input Parameters
+ sopflow - SOPFLOW object
. windgenprofile - wind generator profile file
. c       - contingency number
. t - time step
*/
PetscErrorCode SOPFLOWReadScenarioData_Wind(SOPFLOW sopflow,const char windgenprofile[],PetscInt c,PetscInt t)
{
  PetscErrorCode ierr;
  FILE           *fp;
  char           line[MAXLINE];
  char           *out;
  OPFLOW         opflow;
  PS             ps;
  PetscInt       ngen,*windgenbus,nw=0;
  char           *tok,*tok2;
  char           sep[] = ",",sep2[] = "_";
  PSGEN          gen;
  PetscReal      pg;
  int            scen_num=0;
  int            genid;
  char           windgenid[100][3];

  PetscFunctionBegin;

  if(!sopflow->ismulticontingency) ngen = sopflow->opflows[0]->ps->ngen;
  else {
    if(!sopflow->scopflows[0]->ismultiperiod) ngen = sopflow->scopflows[0]->opflows[c]->ps->ngen;
    else ngen = sopflow->scopflows[0]->tcopflows[c]->opflows[t]->ps->ngen;
  }
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
    if(scen_num != 0 && sopflow->sstart <= scen_num && scen_num < sopflow->send) {
      scen_num -= sopflow->sstart; /* converted to local scenario number */
    } else {
      continue;
    }
    
    if(!sopflow->ismulticontingency) opflow = sopflow->opflows[scen_num];
    else {
      if(!sopflow->scopflows[scen_num]->ismultiperiod) opflow = sopflow->scopflows[scen_num]->opflows[c];
      else opflow = sopflow->scopflows[scen_num]->tcopflows[c]->opflows[t];
    }
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
  SCOPFLOW       scopflow;
  TCOPFLOW       tcopflow;
  PetscInt       c,t;

  PetscFunctionBegin;
  if(sopflow->scenunctype == WIND) {
    if(!sopflow->ismulticontingency) {
      ierr = SOPFLOWReadScenarioData_Wind(sopflow,scenfile,0,0);CHKERRQ(ierr);
    } else {
      scopflow = sopflow->scopflows[0];
      for(c=0; c < scopflow->nc; c++) {
	if(!scopflow->ismultiperiod) {
	  ierr = SOPFLOWReadScenarioData_Wind(sopflow,scenfile,c,0);CHKERRQ(ierr);
	} else {
	  tcopflow = scopflow->tcopflows[0];
	  for(t=0; t < tcopflow->Nt; t++) {
	    ierr = SOPFLOWReadScenarioData_Wind(sopflow,scenfile,c,t);CHKERRQ(ierr);
	  }
	}
      }
    } 
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

