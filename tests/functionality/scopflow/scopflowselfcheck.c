#include <scopflowselfcheck.h>
#include <scopflow.h>
#include <utils.h>
#include <libgen.h>

#define ExaGOSelfcheckSCOPFLOWNumAnswers 6 

/** Array of avaiable reference solutions for SCOPFLOW selfcheck */
static const ExaGOSelfcheckSCOPFLOWAnswer ExaGOSelfcheckSCOPFLOWAnswers[ExaGOSelfcheckSCOPFLOWNumAnswers] = 
{
  /** solver, network file, scenario file, num. contingencies, contingency file, 
   * mode_preventive_or_corrective, is_multiperiod, pload_profile_file, 
   * qload_profile_file, dt, duration, num. iter, unscaled objective value 
   */
  {SCOPFLOWSOLVER_IPOPT,  SELFCHECK_NETWORK_CASE200, "", SCOPFLOW_INITIALIZATION, SCOPFLOW_GENBUSVOLTAGE, 0, SELFCHECK_CONTINGENCY_CASE200, 0,PETSC_FALSE, "", "", 0.0, 0.0, 43, 2.7557570526001036e+04},

  {SCOPFLOWSOLVER_IPOPT,  SELFCHECK_NETWORK_CASE200, "", SCOPFLOW_INITIALIZATION, SCOPFLOW_GENBUSVOLTAGE, 15, SELFCHECK_CONTINGENCY_CASE200, 0,PETSC_FALSE, "", "", 0.0, 0.0, 108, 2.7557571696930871e+04},

  {SCOPFLOWSOLVER_IPOPT,  SELFCHECK_NETWORK_CASE200, "", SCOPFLOW_INITIALIZATION, SCOPFLOW_GENBUSVOLTAGE, 15, SELFCHECK_CONTINGENCY_CASE200, 1,PETSC_FALSE, "", "", 0.0, 0.0, 41, 2.7557570527594460e+04},

  {SCOPFLOWSOLVER_IPOPT,  SELFCHECK_NETWORK_CASE9_WIND, "", SCOPFLOW_INITIALIZATION, SCOPFLOW_GENBUSVOLTAGE, 9, SELFCHECK_CONTINGENCY_CASE9, 0,PETSC_FALSE, "", "", 0.0, 0.0, 14, 3.0492461566531965e+03},
 
  {SCOPFLOWSOLVER_IPOPTNEW,  SELFCHECK_NETWORK_CASE9_WIND, SELFCHECK_SCENARIO_CASE9, SCOPFLOW_INITIALIZATION, SCOPFLOW_GENBUSVOLTAGE, 9, SELFCHECK_CONTINGENCY_CASE9, 1,PETSC_TRUE, SELFCHECK_PLOAD, SELFCHECK_QLOAD, 5.0, 0.16666667, 57, 1.0256987577449478e+05},

  {SCOPFLOWSOLVER_IPOPTNEW,  SELFCHECK_NETWORK_CASE200, SELFCHECK_SCENARIO_CASE200, SCOPFLOW_INITIALIZATION, SCOPFLOW_GENBUSVOLTAGE, 10, SELFCHECK_CONTINGENCY_CASE200, 0,PETSC_TRUE, SELFCHECK_PLOAD, SELFCHECK_QLOAD, 5.0, 0.16666667, 201, 8.6231192297308554e+05}
};

/** Temporary tolerance for error for number of iterations */
static int itertol = 0;

/** 
 * @param[in] scopflow the scopflow object
 * @param[in] ans the answer obtained by running the scopflow driver 
 * @return True if runtime answer matches reference solution, else False 
*/
PetscBool ExaGOSelfcheckSCOPFLOWFindAnswer(SCOPFLOW scopflow, ExaGOSelfcheckSCOPFLOWAnswer *ans)
{
  int i; 
  PetscErrorCode ierr; 
  PetscBool flg; 

  char netfile[PETSC_MAX_PATH_LEN];
  char ctgcfile[PETSC_MAX_PATH_LEN];
  char scenariofile[PETSC_MAX_PATH_LEN];
  char ploadprofile[PETSC_MAX_PATH_LEN];
  char qloadprofile[PETSC_MAX_PATH_LEN];

  /** check for netfile argument */
  ierr = PetscOptionsGetString(NULL,NULL,"-netfile",netfile,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
  if (!flg)
  {
    ExaGOLog(EXAGO_LOG_WARN,"%s needs option -netfile.",__func__);
    return PETSC_TRUE;
  }
  strcpy(ans->networkname,basename(netfile));

  /** check for "solver model" argument */
  ierr = PetscOptionsGetString(NULL,NULL,"-scopflow_solver",&ans->solver,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
  if (!flg)
  {
    ExaGOLog(EXAGO_LOG_WARN,"%s needs option -scopflow_solver.",__func__);
    return PETSC_TRUE;
  }
  
  // Check that a model initialization parameter is being passed
  ierr = PetscOptionsGetString(NULL,NULL,"-opflow_initialization",&ans->modelinit,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
  if (!flg)
  {
    ExaGOLog(EXAGO_LOG_WARN,"%s needs option -opflow_initialization.",__func__);
    return PETSC_TRUE;
  }

  // Check that a model parameter is being passed
  ierr = PetscOptionsGetString(NULL,NULL,"-opflow_genbusvoltage",&ans->genbusvoltage,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
  if (!flg)
  {
    ExaGOLog(EXAGO_LOG_WARN,"%s needs option -opflow_genbusvoltage.",__func__);
    return PETSC_TRUE;
  }

  /** check for contingency file -- all tests has contingency file */
  ierr = PetscOptionsGetString(NULL,NULL,"-ctgcfile",ctgcfile,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
  if (!flg)
  {
    ExaGOLog(EXAGO_LOG_WARN,"%s needs option -ctgcfile with -scopflow_enable_multicontingency.",__func__);
    return PETSC_TRUE;
  }
  strcpy(ans->contingencyname,basename(ctgcfile));

  /** Option check for Multi-contingencies */
  ierr = PetscOptionsGetInt(NULL,NULL,"-scopflow_Nc",&ans->numcontingencies,&flg);CHKERRQ(ierr);
  if (!flg)
  {
    ExaGOLog(EXAGO_LOG_WARN,"%s needs option -scopflow_Nc.",__func__);
    return PETSC_TRUE;
  }

  /** Option check for Multi-contingencies */
  ierr = PetscOptionsGetInt(NULL,NULL,"-scopflow_mode",&ans->mode,&flg);CHKERRQ(ierr);
  if (!flg)
  {
    ExaGOLog(EXAGO_LOG_WARN,"%s needs option -scopflow_mode.",__func__);
    return PETSC_TRUE;
  }

  ierr = PetscOptionsGetBool(NULL,NULL,"-scopflow_enable_multiperiod",&ans->multiperiod,&flg);CHKERRQ(ierr);
  if (!flg)
  {
    ExaGOLog(EXAGO_LOG_WARN,"%s needs option -scopflow_enable_multiperiod.",__func__);
    return PETSC_TRUE;
  }

  if(ans->multiperiod)
  {
    /** Check for the pload profile */
    ierr = PetscOptionsGetString(NULL,NULL,"-scopflow_ploadprofile",ploadprofile,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
    if (!flg)
    {
      ExaGOLog(EXAGO_LOG_WARN,"%s needs option -scopflow_ploadprofile with -scopflow_enable_multiperiod.",__func__);
      return PETSC_TRUE;
    }
    strcpy(ans->ploadprofile,basename(ploadprofile));
    /** Check for the qload profile */
    ierr = PetscOptionsGetString(NULL,NULL,"-scopflow_qloadprofile",qloadprofile,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
    if (!flg)
    {
      ExaGOLog(EXAGO_LOG_WARN,"%s needs option -scopflow_qloadprofile with -scopflow_enable_multiperiod.",__func__);
      return PETSC_TRUE;
    }
    strcpy(ans->qloadprofile,basename(qloadprofile));

    /** Check for the wind scenario file -- reusing variable */
    ierr = PetscOptionsGetString(NULL,NULL,"-scopflow_windgenprofile",scenariofile,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
    if (!flg)
    {
      ExaGOLog(EXAGO_LOG_WARN,"%s needs option -scopflow_windgenprofile with -scopflow_enable_multiperiod.",__func__);
      return PETSC_TRUE;
    }
    strcpy(ans->scenarioname,basename(scenariofile));

    /** Check for dt */
    ierr = PetscOptionsGetReal(NULL,NULL,"-scopflow_dT",&ans->dt,&flg);CHKERRQ(ierr);
    if (!flg)
    {
      ExaGOLog(EXAGO_LOG_WARN,"%s needs option -scopflow_dT with -scopflow_enable_multiperiod.",__func__);
      return PETSC_TRUE;
    }

    /** Check for duration */
    ierr = PetscOptionsGetReal(NULL,NULL,"-scopflow_duration",&ans->duration,&flg);CHKERRQ(ierr);
    if (!flg)
    {
      ExaGOLog(EXAGO_LOG_WARN,"%s needs option -scopflow_duration with -scopflow_enable_multiperiod.",__func__);
      return PETSC_TRUE;
    }
  }
  for(i=0; i < ExaGOSelfcheckSCOPFLOWNumAnswers; ++i) 
  {
    const ExaGOSelfcheckSCOPFLOWAnswer *ians = &ExaGOSelfcheckSCOPFLOWAnswers[i];
    if (strcmp(ans->networkname, ians->networkname)==0 
        && strcmp(ans->solver, ians->solver) == 0 
        && strcmp(ans->modelinit,ians->modelinit) == 0
        && strcmp(ans->genbusvoltage,ians->genbusvoltage) == 0
        && strcmp(ans->contingencyname, ians->contingencyname) == 0 
        && ans->numcontingencies == ians->numcontingencies
        && ans->mode == ians->mode
        )
    {
      if(ans->multiperiod)
      {
        if(strcmp(ans->ploadprofile,ians->ploadprofile) == 0 
           && strcmp(ans->qloadprofile,ians->qloadprofile) == 0 
           && strcmp(ans->scenarioname,ians->scenarioname) == 0 
           && PetscEqualReal(ans->dt, ians->dt) 
           && PetscEqualReal(ans->duration, ians->duration)) 
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

  ExaGOLog(EXAGO_LOG_WARN, "%s could not find answer for options:"
    "\n\tnetwork:%s\n\tsolver:%s\n\tscenario:%s\n\tinitialization:%s\n\tgenbusvoltage:%s\n\tnum_cont:%d\n\tcontingency:%s\n\tmode:%d\n\tmultiperiod:%d\n\tpload:%s\n\tqload:%s\n\tdT:%f\n\tduration:%f\n",
      __func__,ans->networkname,ans->solver,ans->scenarioname,ans->modelinit,
               ans->genbusvoltage,ans->numcontingencies, ans->contingencyname,
               ans->mode, ans->multiperiod,
               ans->ploadprofile, ans->qloadprofile, ans->dt, ans->duration);
  return PETSC_TRUE;
}

PetscBool ExaGOSelfcheckSCOPFLOW(SCOPFLOW scopflow)
{
  PetscErrorCode ierr; 
  PetscInt numiter; 
  double obj_value, tolerance, error;

  int fail = 0, ierror = 0;

  ExaGOSelfcheckSCOPFLOWAnswer *ans;
  ans = (ExaGOSelfcheckSCOPFLOWAnswer*)malloc(sizeof(ExaGOSelfcheckSCOPFLOWAnswer));
  memset( ans->solver, '\0', sizeof(char)* PETSC_MAX_PATH_LEN);
  memset( ans->networkname, '\0', sizeof(char)* PETSC_MAX_PATH_LEN);
  memset( ans->scenarioname, '\0', sizeof(char)* PETSC_MAX_PATH_LEN);
  memset( ans->modelinit, '\0', sizeof(char)* PETSC_MAX_PATH_LEN);
  memset( ans->genbusvoltage, '\0', sizeof(char)* PETSC_MAX_PATH_LEN);
  memset( ans->contingencyname, '\0', sizeof(char)* PETSC_MAX_PATH_LEN);
  memset( ans->ploadprofile, '\0', sizeof(char)* PETSC_MAX_PATH_LEN);
  memset( ans->qloadprofile, '\0', sizeof(char)* PETSC_MAX_PATH_LEN);
  ans->numcontingencies=0;
  ans->mode=0;
  ans->multiperiod=PETSC_FALSE;
  ans->dt=0.0;
  ans->duration=0.0;
  ans->numiter=0;
  ans->objective=0.0;

  ierr = ExaGOSelfcheckSCOPFLOWFindAnswer(scopflow, ans);CHKERRQ(ierr);
  if (ierr)
    return ierr;

  /** Printing config for greping solution/results */
  char config[PETSC_MAX_PATH_LEN];
  sprintf(config,"%s+%s+%s+%s+%s+%d+%s+%d+%d+%s+%s+%f+%f",
          ans->networkname,ans->solver,ans->scenarioname,ans->modelinit,
          ans->genbusvoltage,ans->numcontingencies, ans->contingencyname,
          ans->mode, ans->multiperiod,
          ans->ploadprofile, ans->qloadprofile, ans->dt, ans->duration);
  
  /** Check driver convergence */
  PetscBool conv_status=PETSC_FALSE; 
  ierr = SCOPFLOWGetConvergenceStatus(scopflow,&conv_status);CHKERRQ(ierr);
  if (conv_status==PETSC_FALSE)
  {
    fail++;
    ExaGOLog(EXAGO_LOG_ERROR,"ExaGO SCOPFLOW %s functionality test failed to converge.\n",config);
  }

  /** Compare objective value with reference */
  ierr = SCOPFLOWGetObjective(scopflow, &obj_value);CHKERRQ(ierr);
  ierr = SCOPFLOWGetTolerance(scopflow, &tolerance);CHKERRQ(ierr);
  if (!isEqual(obj_value,ans->objective,tolerance,&error))
  {
    fail++;
    ExaGOLog(EXAGO_LOG_ERROR,
        "ExaGO SCOPFLOW %s functionality test failed to match expected objective value:\n"
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

  /** Compare iterations */
  ierr = SCOPFLOWGetNumIterations(scopflow, &numiter);CHKERRQ(ierr);

  if (ans->numiter != numiter)
  {
    fail++;
    ExaGOLog(EXAGO_LOG_ERROR,
        "ExaGO SCOPFLOW %s functionality test failed to match expected number of iterations:\n"
        "expected=%d actual=%d",
        config,ans->numiter,numiter);
  }
  ExaGOLog(EXAGO_LOG_INFO,"%s: selfcheck returning %s.",
      __func__,fail==0?"success":"failure");

  free(ans);
  return fail;
}