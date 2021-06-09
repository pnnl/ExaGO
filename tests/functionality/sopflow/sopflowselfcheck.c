#include <sopflow.h>
#include <sopflowselfcheck.h>
#include <utils.hpp>
#include <libgen.h>

#define ExaGOSelfcheckSOPFLOWNumAnswers 5

/** Array of available reference solutions for PFLOW selfcheck
 * See test/functionality/utils/selfcheck.h for selcheck solution structure
 */
static const ExaGOSelfcheckSOPFLOWAnswer ExaGOSelfcheckSOPFLOWAnswers[ExaGOSelfcheckSOPFLOWNumAnswers] =
{
  /* Network, number of iterations */
  {SOPFLOWSOLVER_IPOPT   , SELFCHECK_NETWORK_CASE200   , SELFCHECK_SCENARIO_CASE200, SOPFLOW_INITIALIZATION, SOPFLOW_GENBUSVOLTAGE, 1, PETSC_FALSE,  0,                          "", PETSC_FALSE, "", "", 0.0, 0.0, 27, 2.6659883e+04},
  {SOPFLOWSOLVER_IPOPT   , SELFCHECK_NETWORK_CASE200   , SELFCHECK_SCENARIO_CASE200, SOPFLOW_INITIALIZATION, SOPFLOW_GENBUSVOLTAGE, 2, PETSC_FALSE,  0,                          "", PETSC_FALSE, "", "", 0.0, 0.0, 26, 5.3385138e+04},
  {SOPFLOWSOLVER_IPOPT   , SELFCHECK_NETWORK_CASE200   , SELFCHECK_SCENARIO_CASE200, SOPFLOW_INITIALIZATION, SOPFLOW_GENBUSVOLTAGE, 3, PETSC_FALSE,  0,                          "", PETSC_FALSE, "", "", 0.0, 0.0, 35, 8.0297715e+04},
  {SOPFLOWSOLVER_IPOPT, SELFCHECK_NETWORK_CASE9_WIND, SELFCHECK_SCENARIO_CASE9, SOPFLOW_INITIALIZATION, SOPFLOW_GENBUSVOLTAGE, -1,  PETSC_TRUE, -1, SELFCHECK_CONTINGENCY_CASE9, PETSC_FALSE, "", "", 0.0, 0.0, 17, 9.1477392e+03},
  {SOPFLOWSOLVER_IPOPT, SELFCHECK_NETWORK_CASE9_WIND, SELFCHECK_SCENARIO_CASE9, SOPFLOW_INITIALIZATION, SOPFLOW_GENBUSVOLTAGE, -1,  PETSC_TRUE, -1, SELFCHECK_CONTINGENCY_CASE9, PETSC_TRUE, SELFCHECK_PLOAD, SELFCHECK_QLOAD, 5.0, 0.16666667, 74, 3.3125311e+05}
};

/** Temporary tolerance for error for number of iterations */
static int itertol = 0;

PetscErrorCode ExaGOSelfcheckSOPFLOWFindAnswer(SOPFLOW sopflow,ExaGOSelfcheckSOPFLOWAnswer *ans)
{
  PetscErrorCode ierr;
  PetscBool flg;

  char solver[PETSC_MAX_PATH_LEN];
  char netfile[PETSC_MAX_PATH_LEN];
  char scenfile[PETSC_MAX_PATH_LEN];
  char contfile[PETSC_MAX_PATH_LEN];
  char ctgcfile[PETSC_MAX_PATH_LEN];
  char ploadprofile[PETSC_MAX_PATH_LEN];
  char qloadprofile[PETSC_MAX_PATH_LEN];

  // Verify all necessary options are being passed

  // Check that a network is being passed
  ierr = PetscOptionsGetString(NULL,NULL,"-netfile",netfile,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
  if (!flg)
  {
    ExaGOLog(EXAGO_LOG_WARN,"%s needs option -netfile.",__func__);
    return PETSC_TRUE;
  }
  strcpy(ans->networkname,basename(netfile));

  // Check that a solver model is being passed
  ierr = PetscOptionsGetString(NULL,NULL,"-sopflow_solver",ans->solver,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
  if (!flg)
  {
    ExaGOLog(EXAGO_LOG_WARN,"%s needs option -sopflow_solver.",__func__);
    return PETSC_TRUE;
  }

  // Check that a model initialization parameter is being passed
  ierr = PetscOptionsGetString(NULL,NULL,"-opflow_initialization",ans->modelinit,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
  if (!flg)
  {
    ExaGOLog(EXAGO_LOG_WARN,"%s needs option -opflow_initialization.",__func__);
    return PETSC_TRUE;
  }

  // Check that a model parameter is being passed
  ierr = PetscOptionsGetString(NULL,NULL,"-opflow_genbusvoltage",ans->genbusvoltage,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
  if (!flg)
  {
    ExaGOLog(EXAGO_LOG_WARN,"%s needs option -opflow_genbusvoltage.",__func__);
    return PETSC_TRUE;
  }

  // Check that a number of scenarios is being passed
  ierr = PetscOptionsGetInt(NULL,NULL,"-sopflow_Ns",&ans->numscenarios,&flg);CHKERRQ(ierr);
  if (!flg)
  {
    ExaGOLog(EXAGO_LOG_WARN,"%s needs option -sopflow_Ns.",__func__);
    return PETSC_TRUE;
  }

  // Check for scen-file
  ierr = PetscOptionsGetString(NULL,NULL,"-scenfile",scenfile,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
  if (!flg)
  {
    ExaGOLog(EXAGO_LOG_WARN,"%s needs option -scenfile.",__func__);
    return PETSC_TRUE;
  }
  strcpy(ans->scenarioname,basename(scenfile));

  // Check for multi-contingency
  ierr = PetscOptionsGetBool(NULL,NULL,"-sopflow_enable_multicontingency",&ans->multicontingency,&flg);CHKERRQ(ierr);
  fprintf(stderr," -- !\n");
  // Options only applicable for multicontingency...
  if(ans->multicontingency)
  {
    // Check for the number of contingencies
    ierr = PetscOptionsGetInt(NULL,NULL,"-scopflow_Nc",&ans->numcontingencies,&flg);CHKERRQ(ierr);
    if (!flg)
    {
      ExaGOLog(EXAGO_LOG_WARN, "%s needs option -scopflow_Nc with -sopflow_enable_multicontingency",__func__);
      return PETSC_TRUE;
    }

    // Check for the contingency file
    ierr = PetscOptionsGetString(NULL,NULL,"-ctgcfile",ctgcfile,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
    if (!flg)
    {
      ExaGOLog(EXAGO_LOG_WARN,"%s needs option -ctgcfile with -sopflow_enable_multicontingency.",__func__);
      return PETSC_TRUE;
    }
    strcpy(ans->contingencyname,basename(ctgcfile));

    // Check for multiperiod
    ierr = PetscOptionsGetBool(NULL,NULL,"-scopflow_enable_multiperiod",&ans->multiperiod,&flg);CHKERRQ(ierr);
    if (!flg)
    {
      // If the option isn't passed, we can default to false
      // Need this line or else option is indeterministic...
      ans->multiperiod = PETSC_FALSE;
    }
    
    if(ans->multiperiod)
    {
      // Check for the pload profile
      ierr = PetscOptionsGetString(NULL,NULL,"-scopflow_ploadprofile",ploadprofile,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
      if (!flg)
      {
        ExaGOLog(EXAGO_LOG_WARN,"%s needs option -scopflow_ploadprofile with -sopflow_enable_multiperiod.",__func__);
        return PETSC_TRUE;
      }
      strcpy(ans->ploadprofile,basename(ploadprofile));

      // Check for the qload profile
      ierr = PetscOptionsGetString(NULL,NULL,"-scopflow_qloadprofile",qloadprofile,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
      if (!flg)
      {
        ExaGOLog(EXAGO_LOG_WARN,"%s needs option -scopflow_qloadprofile with -sopflow_enable_multiperiod.",__func__);
        return PETSC_TRUE;
      }
      strcpy(ans->qloadprofile,basename(qloadprofile));

      // Check for dt
      ierr = PetscOptionsGetReal(NULL,NULL,"-scopflow_dT",&ans->dt,&flg);CHKERRQ(ierr);
      if (!flg)
      {
        ExaGOLog(EXAGO_LOG_WARN,"%s needs option -scopflow_dT with -sopflow_enable_multiperiod.",__func__);
        return PETSC_TRUE;
      }

      // Check for duration
      ierr = PetscOptionsGetReal(NULL,NULL,"-scopflow_duration",&ans->duration,&flg);CHKERRQ(ierr);
      if (!flg)
      {
        ExaGOLog(EXAGO_LOG_WARN,"%s needs option -scopflow_duration with -sopflow_enable_multiperiod.",__func__);
        return PETSC_TRUE;
      }
    }
  }

  int i;
  for (i=0;i<ExaGOSelfcheckSOPFLOWNumAnswers;i++)
  {
    const ExaGOSelfcheckSOPFLOWAnswer *ians = &ExaGOSelfcheckSOPFLOWAnswers[i];
    if (strcmp(ans->networkname,ians->networkname)==0 &&
        strcmp(ans->solver,ians->solver) == 0 &&
        strcmp(ans->modelinit,ians->modelinit) == 0 &&
        strcmp(ans->genbusvoltage,ians->genbusvoltage) == 0 &&
        ans->numscenarios == ians->numscenarios &&
        strcmp(ans->scenarioname,ians->scenarioname) == 0 &&
        ans->multicontingency == ians->multicontingency)
    {
      // Perform more checks if multicontingency
      if(ans->multicontingency)
      {
        if (ans->numcontingencies == ians->numcontingencies &&
            strcmp(ans->contingencyname,ians->contingencyname) == 0 &&
            ans->multiperiod == ians->multiperiod)
        {
          if(ans->multiperiod)
          {
            if(strcmp(ans->ploadprofile,ians->ploadprofile) == 0 &&
               strcmp(ans->qloadprofile,ians->qloadprofile) == 0 &&
               PetscEqualReal(ans->dt, ians->dt) && 
               PetscEqualReal(ans->duration, ians->duration))
            {
              ans->numiter = ians->numiter;
              ans->objective = ians->objective;
              return PETSC_FALSE;
            }
          }
          else
          {
            ans->numiter = ians->numiter;
            ans->objective = ians->objective;
            return PETSC_FALSE;
          }
        }
      }
      else
      {
        ans->numiter = ians->numiter;
        ans->objective = ians->objective;
        return PETSC_FALSE;
      }
    }
  }
  ExaGOLog(EXAGO_LOG_WARN,"%s could not find answer for option:"
      "\n\tnetwork:%s\n\tsolver:%s\n\tscenario:%s\n\tinitialization:%s\
      \n\tgenbusvoltage:%s\n\tnum_scenarios:%d\n\tmulticontingency:%d\
       \n\tnum_cont:%d\n\tcontingency:%s\n\tmultiperiod:%d\
       \n\tpload:%s\n\tqload:%s\n\tdT:%f\n\tduration:%f\n",
      __func__,ans->networkname,ans->solver,ans->scenarioname,ans->modelinit,\
               ans->genbusvoltage,ans->numscenarios,ans->multicontingency,\
               ans->numcontingencies, ans->contingencyname, ans->multiperiod,\
               ans->ploadprofile, ans->qloadprofile, ans->dt, ans->duration);
  return PETSC_TRUE;
}

PetscErrorCode ExaGOSelfcheckSOPFLOW(SOPFLOW sopflow)
{
  PetscErrorCode ierr;
  PetscInt numiter;
  double obj_value,tolerance,error;
  int fail = 0, ierror = 0;

  ExaGOSelfcheckSOPFLOWAnswer *ans;
  ans = (ExaGOSelfcheckSOPFLOWAnswer*)malloc(sizeof(ExaGOSelfcheckSOPFLOWAnswer));
  memset( ans->solver, '\0', sizeof(char)* PETSC_MAX_PATH_LEN);
  memset( ans->networkname, '\0', sizeof(char)* PETSC_MAX_PATH_LEN);
  memset( ans->scenarioname, '\0', sizeof(char)* PETSC_MAX_PATH_LEN);
  memset( ans->modelinit, '\0', sizeof(char)* PETSC_MAX_PATH_LEN);
  memset( ans->genbusvoltage, '\0', sizeof(char)* PETSC_MAX_PATH_LEN);
  memset( ans->contingencyname, '\0', sizeof(char)* PETSC_MAX_PATH_LEN);
  memset( ans->ploadprofile, '\0', sizeof(char)* PETSC_MAX_PATH_LEN);
  memset( ans->qloadprofile, '\0', sizeof(char)* PETSC_MAX_PATH_LEN);
  ans->numscenarios = 0;
  ans->multicontingency = PETSC_FALSE;
  ans->numcontingencies=0;
  ans->multiperiod=PETSC_FALSE;
  ans->dt=0.0;
  ans->duration=0.0;
  ans->numiter=0;
  ans->objective=0.0;

  ierr = ExaGOSelfcheckSOPFLOWFindAnswer(sopflow, ans);CHKERRQ(ierr);
  if (ierr)
  {
    free(ans);
    return ierr;
  }

  // Printing config makes greping for solutions/results easier
  char config[PETSC_MAX_PATH_LEN];
  sprintf(config,"%s+%s+%s+%s+%s+%d+%d+%d+%s+%d+%s+%s+%d+%d",\
          ans->networkname,ans->solver,ans->scenarioname,ans->modelinit,\
          ans->genbusvoltage,ans->numscenarios,ans->multicontingency,\
          ans->numcontingencies, ans->contingencyname, ans->multiperiod,\
          ans->ploadprofile, ans->qloadprofile, ans->dt, ans->duration);

  // Compare solution to reference

  // Check for convergence
  PetscBool conv_status=PETSC_FALSE;
  ierr = SOPFLOWGetConvergenceStatus(sopflow,&conv_status);CHKERRQ(ierr);
  if (conv_status==PETSC_FALSE)
  {
    fail++;
    ExaGOLog(EXAGO_LOG_ERROR,"ExaGO SOPFLOW %s functionality test failed to converge.\n",config);
  }

  // Compare objective value
  ierr = SOPFLOWGetObjective(sopflow, &obj_value);CHKERRQ(ierr);
  ierr = SOPFLOWGetTolerance(sopflow, &tolerance);CHKERRQ(ierr);
  if (!IsEqual(obj_value,ans->objective,tolerance,error))
  {
    fail++;
    ExaGOLog(EXAGO_LOG_ERROR,
        "ExaGO SOPFLOW %s functionality test failed to match expected objective value:\n"
        "%-30s%.12e\n"
        "%-30s%.12e\n"
        "%-30s%.12e\n"
        "%-30s%.12e\n",
        config,
        "Expected:",ans->objective,
        "Actual:",obj_value,
        "Tolerance:",tolerance,
        "Scaled Error:",error);
  }

  // Compare iterations
  ierr = SOPFLOWGetNumIterations(sopflow, &numiter);CHKERRQ(ierr);
  if (ans->numiter != numiter)
  {
    fail++;
    ExaGOLog(EXAGO_LOG_ERROR,
        "ExaGO SOPFLOW %s functionality test failed to match expected number of iterations:\n"
        "expected=%d actual=%d",
        config,ans->numiter,numiter);
  }

  ExaGOLog(EXAGO_LOG_INFO,"%s: selfcheck returning %s.",
      __func__,fail==0?"success":"failure");

  free(ans);
  return fail;
}
