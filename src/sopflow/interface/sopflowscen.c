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

/*
  SOPFLOWReadScenarioData_Wind_SinglePeriod - Reads the wind data and populates the scenario list
  Input Parameters
+ sopflow - SOPFLOW object
. windgenprofile - wind generator profile file

  Note: This function reads the wind scenario data for "single period" format files.
*/
PetscErrorCode SOPFLOWReadScenarioData_Wind_SinglePeriod(SOPFLOW sopflow,const char windgenprofile[])
{
  PetscErrorCode ierr;
  FILE           *fp;
  char           line[MAXLINE];
  char           *out;
  PetscInt       ngenwind,nw=0;
  char           *tok,*tok2;
  char           sep[] = ",",sep2[] = "_";
  PetscReal      pg,weight;
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
  while(tok != NULL) {
    if(strcmp(tok,"weight") == 0 || strcmp(tok,"weight\n") == 0) {
      tok = strtok(NULL,sep);
      continue;
    }
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
    for(i=0; i < forecast->nele; i++) {
      ierr = PetscCalloc1(3,&forecast->id[i]);
    }
    ierr = PetscCalloc1(forecast->nele,&forecast->val);CHKERRQ(ierr);
    
    tok = strtok(NULL,sep);
    for(i=0; i < nw; i++) {
      forecast->buses[i] = windgenbus[i];
      ierr = PetscStrcpy(forecast->id[i],windgenid[i]);CHKERRQ(ierr);
      sscanf(tok,"%lf",&pg);
      forecast->val[i] = pg;
      tok = strtok(NULL,sep);
    }

    /* Read scenario weight */
    sscanf(tok,"%lf",&weight);
    scenario->prob = weight;

    scenario->nforecast++;
    scenlist->Nscen++;
  }
  PetscFunctionReturn(0);
}

/*
  SOPFLOWReadScenarioData_Wind_MultiPeriod - Reads the wind data and populates the scenario list
  Input Parameters
+ sopflow - SOPFLOW object
. windgenprofile - wind generator profile file

  Note: This function reads the wind scenario data for "multi-period" format files.
*/
PetscErrorCode SOPFLOWReadScenarioData_Wind_MultiPeriod(SOPFLOW sopflow,const char windgenprofile[])
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
    for(i=0; i < forecast->nele; i++) {
      ierr = PetscCalloc1(3,&forecast->id[i]);
    }
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
    if(scenfileformat == SOPFLOW_NATIVE_SINGLEPERIOD) {
      ierr = SOPFLOWReadScenarioData_Wind_SinglePeriod(sopflow,scenfile);CHKERRQ(ierr);
    } else if(scenfileformat == SOPFLOW_NATIVE_MULTIPERIOD) {
      ierr = SOPFLOWReadScenarioData_Wind_MultiPeriod(sopflow,scenfile);CHKERRQ(ierr);
    }
  }

  PetscFunctionReturn(0);
}

/* SOPFLOWGetNumScenarios_Native_MultiPeriod - Gets the number of scenarios from the scenario file
 */
PetscErrorCode SOPFLOWGetNumScenarios_Native_MultiPeriod(const char scenfile[],PetscInt *Ns)
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

/* SOPFLOWGetNumScenarios_Native_SinglePeriod - Gets the number of scenarios from the scenario file
 */
PetscErrorCode SOPFLOWGetNumScenarios_Native_SinglePeriod(const char scenfile[],PetscInt *Ns)
{
  PetscErrorCode ierr;
  FILE           *fp;
  char           line[MAXLINE];
  char           *out;
  char           *tok,*tok2;
  char           sep[] = ",",sep2[] = "_";
  int            Nscen=0;

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
    sscanf(tok,"%d",&Nscen); /* Scenario number */
  }

  fclose(fp);
  *Ns = Nscen;
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
  if(scenfileformat == SOPFLOW_NATIVE_SINGLEPERIOD) {
    ierr = SOPFLOWGetNumScenarios_Native_SinglePeriod(scenfile,Ns);CHKERRQ(ierr);
  } else if(scenfileformat == SOPFLOW_NATIVE_MULTIPERIOD) {
    ierr = SOPFLOWGetNumScenarios_Native_MultiPeriod(scenfile,Ns);CHKERRQ(ierr);
  }
  else *Ns = 1;

  PetscFunctionReturn(0);
}

