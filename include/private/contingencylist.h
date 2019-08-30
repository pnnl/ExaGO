/**
 * contingencylist.h
 * Private header file for contingency list
 */

#ifndef CONTINGENCYLIST_H
#define CONTINGENCYLIST_H

#include <ps.h>
#include <private/psimpl.h>

typedef enum {GEN_OUTAGE=1, BR_OUTAGE=2,TR_OUTAGE=3,LOAD_OUTAGE=4} OutageType;

#define MAX_SIMULTANEOUS_OUTAGES 5

struct _p_Outage {
  PetscInt   num;  /* Outage number */
  OutageType type; /* Outage type */
  PetscInt   bus;  /* Generator bus number */
  PetscInt   fbus,tbus; /* Branch or transformer from and to buses */
  char       id[10]; /* Generator, branch, or transformer id */
  PetscInt   status; /* status of the equipment */
  PetscScalar prob; /* Probability of the outage */
};

typedef struct _p_Outage Outage;

struct _p_Contingency {
  PetscInt noutages; /* Each contingency can have one or more outages */
  Outage   outagelist[MAX_SIMULTANEOUS_OUTAGES]; /* List of outages for this contingency */
};

typedef struct _p_Contingency Contingency;

struct _p_ContingencyList {
  PetscInt Ncont; /* Number of contingencies = number of scenarios */
  Contingency *cont; /* Contingencies */
};

typedef struct _p_ContingencyList ContingencyList;

#endif
