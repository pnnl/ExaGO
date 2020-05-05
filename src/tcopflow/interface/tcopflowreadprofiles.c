#include <private/opflowimpl.h>
#include <private/tcopflowimpl.h>

/*
  TCOPFLOWReadPloadProfile - Reads the active power load profile

  Input Parameters:
+ tcopflow - the TCOPFLOW object
- ploadprofile - the file containing the real load power profiles

  Note: TCOPFLOWRReadPloadProfile parses the load profile files created by NREL
*/
PetscErrorCode TCOPFLOWReadPloadProfile(TCOPFLOW tcopflow, char ploadprofile[])
{
  PetscErrorCode ierr;
  FILE           *fp;
  char           line[MAXLINE];
  char           *out;
  PetscInt       i;
  OPFLOW         opflow;
  PS             ps;
  PetscInt       nload=tcopflow->opflows[0]->ps->nload,*lbus,nl=0;
  char           *tok;
  char           *sep = ",";
  PSLOAD         load;
  PetscInt       t=0;
  PetscReal      pl;

  PetscFunctionBegin;

  fp = fopen(ploadprofile,"r");
  if (fp == NULL) {
    SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Cannot open P load profile file %s",ploadprofile);CHKERRQ(ierr);
  }

  ierr = PetscMalloc1(nload,&lbus);CHKERRQ(ierr);
  /* First line -- has the bus numbers */
  out = fgets(line,MAXLINE,fp);

  /* Parse load bus numbers */
  tok = strtok(line,sep);
  tok = strtok(NULL,sep); /* Skip first token */
  while(tok != NULL) {
    sscanf(tok,"%d",&lbus[nl]);
    nl++;
    tok = strtok(NULL,sep);
  }

  while((out = fgets(line,MAXLINE,fp)) != NULL) {
    if(strcmp(line,"\r\n") == 0 || strcmp(line,"\n") == 0) {
      continue; /* Skip blank lines */
    }

    opflow = tcopflow->opflows[t];
    ps     = opflow->ps;
    /* Parse load bus numbers */
    tok = strtok(line,sep);
    tok = strtok(NULL,sep); /* Skip first token */
    nl = 0;
    while(tok != NULL) {
      ierr = PSGetLoad(ps,lbus[nl],"1 ",&load);CHKERRQ(ierr);
      sscanf(tok,"%lf",&pl);
      load->pl = pl/ps->MVAbase;
      nl++;
      tok = strtok(NULL,sep);
    }

    t++;
    if(t == tcopflow->Nt) break;
  }

  ierr = PetscFree(lbus);CHKERRQ(ierr);
  fclose(fp);
  PetscFunctionReturn(0);
}


/*
  TCOPFLOWReadQloadProfile - Reads the reactive power load profile

  Input Parameters:
+ tcopflow - the TCOPFLOW object
- qloadprofile - the file containing the real load power profiles

  Note: TCOPFLOWReadQloadProfile parses the load profile files created by NREL
*/
PetscErrorCode TCOPFLOWReadQloadProfile(TCOPFLOW tcopflow, char qloadprofile[])
{
  PetscErrorCode ierr;
  FILE           *fp;
  char           line[MAXLINE];
  char           *out;
  PetscInt       i;
  OPFLOW         opflow;
  PS             ps;
  PetscInt       nload=tcopflow->opflows[0]->ps->nload,*lbus,nl=0;
  char           *tok;
  char           *sep = ",";
  PSLOAD         load;
  PetscInt       t=0;
  PetscReal      ql;

  PetscFunctionBegin;

  fp = fopen(qloadprofile,"r");
  if (fp == NULL) {
    SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Cannot open Q load profile file %s",qloadprofile);CHKERRQ(ierr);
  }

  ierr = PetscMalloc1(nload,&lbus);CHKERRQ(ierr);
  /* First line -- has the bus numbers */
  out = fgets(line,MAXLINE,fp);

  /* Parse load bus numbers */
  tok = strtok(line,sep);
  tok = strtok(NULL,sep); /* Skip first token */
  while(tok != NULL) {
    sscanf(tok,"%d",&lbus[nl]);
    nl++;
    tok = strtok(NULL,sep);
  }

  while((out = fgets(line,MAXLINE,fp)) != NULL) {
    if(strcmp(line,"\r\n") == 0 || strcmp(line,"\n") == 0) {
      continue; /* Skip blank lines */
    }

    opflow = tcopflow->opflows[t];
    ps     = opflow->ps;
    /* Parse load bus numbers */
    tok = strtok(line,sep);
    tok = strtok(NULL,sep); /* Skip first token */
    nl = 0;
    while(tok != NULL) {
      ierr = PSGetLoad(ps,lbus[nl],"1 ",&load);CHKERRQ(ierr);
      sscanf(tok,"%lf",&ql);
      load->ql = ql/ps->MVAbase;
      nl++;
      tok = strtok(NULL,sep);
    }

    t++;
    if(t == tcopflow->Nt) break;
  }

  ierr = PetscFree(lbus);CHKERRQ(ierr);
  fclose(fp);
  PetscFunctionReturn(0);
}

/*
  TCOPFLOWReadWindGenProfile - Reads the wind generation

  Input Parameters:
+ tcopflow - the TCOPFLOW object
- windgenprofile - the file containing the wind generation profiles

  Note: TCOPFLOWReadWindGenProfile parses the load profile files created by NREL
*/
PetscErrorCode TCOPFLOWReadWindGenProfile(TCOPFLOW tcopflow, char windgenprofile[])
{
  PetscErrorCode ierr;
  FILE           *fp;
  char           line[MAXLINE];
  char           *out;
  PetscInt       i;
  OPFLOW         opflow;
  PS             ps;
  PetscInt       ngen=tcopflow->opflows[0]->ps->ngen,*windgenbus,nw=0;
  char           *tok,*tok2;
  char           *sep = ",",*sep2 = "_";
  PSGEN          gen;
  PetscInt       t=0;
  PetscReal      pg;
  int            scen_num;
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

    opflow = tcopflow->opflows[t];
    ps     = opflow->ps;
    /* Parse load bus numbers */
    tok = strtok(line,sep);
    tok = strtok(NULL,sep); /* Skip first token */
    sscanf(tok,"%d",&scen_num); /* Scenario number */
    if(scen_num != 1) continue;
    tok = strtok(NULL,sep); /* Scenario number */
    nw = 0;
    while(tok != NULL) {
      ierr = PSGetGen(ps,windgenbus[nw],windgenid[nw],&gen);CHKERRQ(ierr);
      sscanf(tok,"%lf",&pg);
      gen->pg = gen->pt = pg/ps->MVAbase; /* Set real power generation. Note that Pg upper limit is also set to Pg */
      nw++;
      tok = strtok(NULL,sep);
    }

    t++;
    if(t == tcopflow->Nt) break;
  }

  ierr = PetscFree(windgenbus);CHKERRQ(ierr);

  fclose(fp);
  PetscFunctionReturn(0);
}

