#include <private/opflowimpl.h>
#include <private/tcopflowimpl.h>
#include <private/scopflowimpl.h>
#include <private/sopflowimpl.h>
#include <private/scenariolist.h>

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

  ierr = PetscMemcpy(sopflow->scenfile,scenfile,PETSC_MAX_PATH_LEN*sizeof(char));CHKERRQ(ierr);

  sopflow->scenfileformat = scenfileformat;
  sopflow->scenfileset = PETSC_TRUE;
  sopflow->scenunctype = scenunctype;
  PetscFunctionReturn(0);
}

#if 0
/* This version will be removed in the future release. It used to read the scenario
   data and update the PS objects. This was a cumbersome way to manipulate the data, hence its being obsoleted.
*/

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
  int            scensread=0;
  int            nscens;

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

  /* scenarios on this rank excluding base */
  //  nscens=sopflow->comm->rank?sopflow->ns:sopflow->ns-1;
  nscens = sopflow->ns;

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
    scen_num -= 1; /* Scenario numbers start with 1 in the file, convert to zero-based start */

    if(scensread == nscens) {
      ierr = PetscFree(windgenbus);CHKERRQ(ierr);
      fclose(fp);
      PetscFunctionReturn(0);
    }

    if(scen_num < sopflow->sstart) { /* Haven't yet reached scenarios for this rank */
      continue;
    } else if(sopflow->sstart <= scen_num && scen_num < sopflow->send) { /* This is my range */
      scen_num -= sopflow->sstart; /* convert to local scenario number */
      scensread++;
    } else { /* Beyond range, we can exit */
      ierr = PetscFree(windgenbus);CHKERRQ(ierr);
      fclose(fp);
      PetscFunctionReturn(0);
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
      gen->pgs = gen->pb;
      nw++;
      tok = strtok(NULL,sep);
    }
  }
  PetscFunctionReturn(0);
}
#endif

/*
  SOPFLOWReadScenarioData_Wind - Reads the wind data and populates the scenario list
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
  PetscInt       ngenwind,nw=0;
  char           *tok,*tok2;
  char           sep[] = ",",sep2[] = "_";
  PetscReal      pg;
  int            scen_num=0;
  int            genid;
  int            windgenbus[100];
  char           windgenid[100][3];
  int            i;
  ScenarioList   *scenlist=&sopflow->scenlist;
  Scenario       *scenario;
  Forecast       *forecast;

  PetscFunctionBegin;

  ngenwind = 100; // This should be increased for larger cases (?)
  fp = fopen(windgenprofile,"r");
  if (fp == NULL) {
    SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Cannot open wind generation profile file %s",windgenprofile);CHKERRQ(ierr);
  }

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
    if(nw == ngenwind) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_SUP,"Exceeded max. number of wind generators=%d\n",ngenwind);
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
    scen_num -= 1; /* Scenario numbers start with 1 in the file, convert to zero-based start */

    if(scen_num < scenlist->Nscen || scen_num == sopflow->Ns) {
      fclose(fp);
      PetscFunctionReturn(0);
    }

    scenario = &scenlist->scen[scen_num];
    forecast = &scenario->forecastlist[scenario->nforecast];
    forecast->num = scen_num;
    forecast->type = FORECAST_WIND;
    forecast->nele = nw;
    ierr = PetscCalloc1(forecast->nele,&forecast->buses);CHKERRQ(ierr);
    ierr = PetscCalloc1(forecast->nele,&forecast->id);CHKERRQ(ierr);
    ierr = PetscCalloc1(forecast->nele,&forecast->val);CHKERRQ(ierr);
    
    for(i=0; i < nw; i++) {
      forecast->buses[i] = windgenbus[i];
      ierr = PetscStrcpy(forecast->id[i],windgenid[i]);CHKERRQ(ierr);
    }

    tok = strtok(NULL,sep);
    nw = 0;
    while(tok != NULL) {
      sscanf(tok,"%lf",&pg);
      forecast->val[nw] = pg;
      nw++;
      tok = strtok(NULL,sep);
    }
    scenario->nforecast++;
    scenlist->Nscen++;
  }
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
    ierr = SOPFLOWReadScenarioData_Wind(sopflow,scenfile,0,0);CHKERRQ(ierr);
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

