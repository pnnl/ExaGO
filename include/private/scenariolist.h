/**
 * scenariolist.h
 * Private header file for scenario list
 */

#ifndef SCENARIOLIST_H
#define SCENARIOLIST_H

#include <private/psimpl.h>
#include <ps.h>

#define MAX_FORECASTS_PER_SCENARIO 5
#define MAX_SCENARIOS 10

struct _p_Forecast {
  PetscInt num;      /* Scenario number */
  ForecastType type; /* Forecast type */
  PetscInt nele;     /* Number of devices/elements involved in this forecast */
  PetscInt *buses;   /* Bus numbers */
  char **id;         /* Device ids */
  PetscScalar *val;  /* forecast values */
};

typedef struct _p_Forecast Forecast;

struct _p_Scenario {
  PetscInt nforecast; /* Each scenario can have one or more forecasts */
  Forecast forecastlist[MAX_FORECASTS_PER_SCENARIO]; /* List of forecasts for
                                                        this scenario */
  PetscScalar prob; /* Probability of the scenario */
};

typedef struct _p_Scenario Scenario;

struct _p_ScenarioList {
  PetscInt Nscen; /* Number of scenarios = number of scenarios */
  Scenario *scen; /* Scenarios */
};

typedef struct _p_ScenarioList ScenarioList;

#endif
