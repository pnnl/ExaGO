#include <private/opflowimpl.h>
#include <private/scopflowimpl.h>

extern void clean2Char(char *);
extern char** blankTokenizer(const char *str, int *numtok, int maxtokens, int maxchar);
/*
  SCOPFLOWSetContingencyData - Sets the contingency data

  Input Parameter
+  scopflow - The SCOPFLOW object
.  ctgcfileformat - the contingency file format
-  ctgcfile - The name of the contingency list file

*/
PetscErrorCode SCOPFLOWSetContingencyData(SCOPFLOW scopflow,ContingencyFileInputFormat ctgcfileformat,const char ctgcfile[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = PetscMemcpy(scopflow->ctgcfile,ctgcfile,PETSC_MAX_PATH_LEN*sizeof(char));CHKERRQ(ierr);

  scopflow->ctgcfileformat = ctgcfileformat;
  scopflow->ctgcfileset = PETSC_TRUE;
  PetscFunctionReturn(0);
}

/*
  SCOPFLOWReadContingencyData_Native - Reads the contingency list data file in the native format

  Input Parameters
+ scopflow - the scopflow object
- ctgcfile - the contingency file name

  Notes:
 * The native format for Contingency input file is as follows:
 * Notes: Each field in the contingency file has the following format
 * Num,Type,Bus,Fbus,Tbus,Id,Status,prob
 * Num - Contingency number
 * Type - Type of contingency (Generator, Branch, Transformer, Load)
 * Bus - The equipment bus number
 * Fbus - From bus number (only for branch and transformer contingencies)
 * Tbus - To bus number (only for branch and transformer contingencies)
 * Id   - The equipment ID (2-char string)
 * Status - The status to be set for the outaged equipment
 * Prob   - the probability of the outage

  Examples:
(Outage generator at bus 1 with Id 1, probability = 0.1)
    1,0,1,0,0,1 ,0,0.1  
(Outage branch connecting buses 8-9 with Id 1, probability = 0.1)
    2,1,0,8,9,1 ,0,0.1
(Multiple outages
  generator at bus 1 with Id 1, probability = 0.1
  branch connecting buses 8-9 with Id 1, probability = 0.1)
    3,0,1,0,0,1 ,0,0.1  
    3,1,0,8,9,1 ,0,0.1
*/
PetscErrorCode SCOPFLOWReadContingencyData_Native(SCOPFLOW scopflow,const char ctgcfile[])
{
  PetscErrorCode ierr;
  FILE           *fp;
  ContingencyList *ctgclist=&scopflow->ctgclist;
  Contingency    *cont;
  Outage         *outage;
  char           line[MAXLINE];
  char           *out;
  PetscInt       bus,fbus,tbus,type,num;
  char           equipid[3];
  PetscInt       status;
  PetscScalar    prob;

  PetscFunctionBegin;

  fp = fopen(ctgcfile,"r");
  if (fp == NULL) {
    SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Cannot open file %s",ctgcfile);CHKERRQ(ierr);
  }

  ctgclist->Ncont = -1;

  while((out = fgets(line,MAXLINE,fp)) != NULL) {
    if(strcmp(line,"\r\n") == 0 || strcmp(line,"\n") == 0) {
      continue; /* Skip blank lines */
    }
    sscanf(line,"%d,%d,%d,%d,%d,'%[^\t\']',%d,%lf",&num,&type,&bus,&fbus,&tbus,equipid,&status,&prob);

    if(num == scopflow->Nc) break;

    if(num == PetscMax(scopflow->Nc,MAX_CONTINGENCIES)) {
      SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_SUP,"Exceeding max. allowed contingencies = %d\n",num,MAX_CONTINGENCIES);
    }
    cont   = &ctgclist->cont[num];
    outage = &cont->outagelist[cont->noutages];
    outage->num  = num;
    outage->type = (OutageType)type;
    outage->bus  = bus;
    outage->fbus = fbus;
    outage->tbus = tbus;
    ierr = PetscMemcpy(outage->id,equipid,3*sizeof(char));CHKERRQ(ierr);
    outage->status = status;
    outage->prob   = prob;
    cont->noutages++;


    if(num > ctgclist->Ncont) ctgclist->Ncont = num;
  }
  fclose(fp);

  PetscFunctionReturn(0);
}

/*
  SCOPFLOWReadContingencyData_PSSE - Reads the contingency list data file in the PSSE format (.con file

  Input Parameters
+ scopflow - the scopflow object
- ctgcfile - the contingency file name

*/
PetscErrorCode SCOPFLOWReadContingencyData_PSSE(SCOPFLOW scopflow,const char ctgcfile[])
{
  PetscErrorCode ierr;
  FILE           *fp;
  ContingencyList *ctgclist=&scopflow->ctgclist;
  Contingency    *cont;
  Outage         *outage;
  char           *out;
  char           line[MAXLINE];
  PetscInt       i,one,two;
  PetscInt       reading=0;
  PetscInt       numCont=0;
  char           equipid[32];
  PetscInt       ntokens;
  PetscInt       maxtokens=10;
  PetscInt       maxchar=32;

  PetscFunctionBegin;

  fp = fopen(ctgcfile,"r");
  if (fp == NULL) {
    SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Cannot open file %s",ctgcfile);CHKERRQ(ierr);
  }

  ctgclist->Ncont = -1;
  one = 1;
  two = 2;
  while((out = fgets(line,MAXLINE,fp)) != NULL && numCont < scopflow->Nc-1) {
    if(strcmp(line,"\r\n") == 0 || strcmp(line,"\n") == 0) {
      continue; /* Skip blank lines */
    }
    char **tokens = blankTokenizer(line,&ntokens,maxtokens,maxchar);
    if (!reading) {
      if (!strcmp(tokens[0],"CONTINGENCY")) {
        reading = 1;
      }
    } else {
      numCont++;
      while (strcmp(tokens[0],"END")) {
        if (!strcmp(tokens[0],"REMOVE")) {
          // Generator contingency
          cont   = &ctgclist->cont[numCont];
          outage = &cont->outagelist[cont->noutages];
          outage->num  = numCont;
          outage->type = (OutageType)one;
          outage->bus  = atoi(tokens[5]);
          outage->fbus = 0;
          outage->tbus = 0;
          strcpy(equipid,tokens[2]);
          clean2Char(equipid);
          ierr = PetscMemcpy(outage->id,equipid,3*sizeof(char));CHKERRQ(ierr);
          outage->status = 0;
          outage->prob   = 0.01;
          cont->noutages++;
        } else if (!strcmp(tokens[0],"OPEN")) {
          cont   = &ctgclist->cont[numCont];
          outage = &cont->outagelist[cont->noutages];
          outage->num  = numCont;
          outage->type = (OutageType)two;
          outage->bus  = 0;
          outage->fbus = atoi(tokens[4]);
          outage->tbus = atoi(tokens[6]);
          strcpy(equipid,tokens[8]);
          clean2Char(equipid);
          ierr = PetscMemcpy(outage->id,equipid,3*sizeof(char));CHKERRQ(ierr);
          outage->status = 0;
          outage->prob   = 0.01;
          cont->noutages++;
        } else {
          free(tokens);
          SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,
             "Cannot Identify contingency\n");
        }
        for (i=0; i<maxtokens; i++) {
          free(tokens[i]);
        }
        free(tokens);
        if ((out = fgets(line,MAXLINE,fp)) == NULL) {
          SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,
              "End of file encountered before END statement\n");
        }
        tokens = blankTokenizer(line,&ntokens,maxtokens,maxchar);
      }
      reading = 0;
    }
    for (i=0; i<maxtokens; i++) {
      free(tokens[i]);
    }
    free(tokens);
    if(numCont > ctgclist->Ncont) ctgclist->Ncont = numCont;
  }
  fclose(fp);
  PetscFunctionReturn(0);
}

/*
  SCOPFLOWReadContingencyData - Reads the contingency list data file

  Input Parameters
+ scopflow - the scopflow object
. ctgcfileformat - the contingency file format (NATIVE or PSSE)
- ctgcfile - the contingency file name

*/
PetscErrorCode SCOPFLOWReadContingencyData(SCOPFLOW scopflow,ContingencyFileInputFormat ctgcfileformat,const char ctgcfile[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if(ctgcfileformat == NATIVE) {
    ierr = SCOPFLOWReadContingencyData_Native(scopflow,ctgcfile);
  } else if(ctgcfileformat == PSSE) {
    ierr = SCOPFLOWReadContingencyData_PSSE(scopflow,ctgcfile);
  } else {
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Unknown contingency input file format %s\n");
  }

  PetscFunctionReturn(0);
}
