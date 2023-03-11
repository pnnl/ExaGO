/**
 * contingencylist.h
 * Private header file for contingency list
 */

#ifndef CONTINGENCYLIST_H
#define CONTINGENCYLIST_H

#include <private/psimpl.h>
#include <ps.h>

typedef enum {
  GEN_OUTAGE = 1,
  BR_OUTAGE = 2,
  TR_OUTAGE = 3,
  LOAD_OUTAGE = 4
} OutageType;

#define MAX_SIMULTANEOUS_OUTAGES 5
#define MAX_CONTINGENCIES 20000

struct _p_Outage {
  PetscInt num;        /* Outage number */
  OutageType type;     /* Outage type */
  PetscInt bus;        /* Generator bus number */
  PetscInt fbus, tbus; /* Branch or transformer from and to buses */
  char id[3];          /* Generator, branch, or transformer id */
  PetscInt status;     /* status of the equipment */
  PetscScalar prob;    /* Probability of the outage */
};

typedef struct _p_Outage Outage;

struct _p_Contingency {
  PetscInt noutages; /* Each contingency can have one or more outages */
  Outage outagelist[MAX_SIMULTANEOUS_OUTAGES]; /* List of outages for this
                                                  contingency */
};

typedef struct _p_Contingency Contingency;

struct _p_ContingencyList {
  PetscInt Ncontinit; /* Initial size of list */
  PetscInt Ncont;     /* Number of contingencies = number of scenarios */
  Contingency *cont;  /* Contingencies */
  char inputfile[PETSC_MAX_PATH_LEN];         /* input file */
  ContingencyFileInputFormat inputfileformat; /* format of input file */

  PetscLogEvent readctgcdatalogger; /* Profile for reading contingency data */
};

typedef struct _p_ContingencyList *ContingencyList;

extern PetscErrorCode ContingencyListSetData(ContingencyList,
                                             ContingencyFileInputFormat,
                                             const char[]);
extern PetscErrorCode ContingencyListCreate(PetscInt, ContingencyList *);
extern PetscErrorCode ContingencyListDestroy(ContingencyList *);
extern PetscErrorCode ContingencyListReadData(ContingencyList, PetscInt *);
extern PetscErrorCode ContingencyWriteData(Contingency *con, int id,
                                           ContingencyFileInputFormat fmt,
                                           FILE *fd);
extern PetscErrorCode ContigencyAppendMinimal(Contingency *thecon, int cont_num,
                                              ContingencyFileInputFormat fmt,
                                              const char outfile[]);

#endif
