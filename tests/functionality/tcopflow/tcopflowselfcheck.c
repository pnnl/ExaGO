#include <tcopflow.h>
#include <tcopflowselfcheck.h>
#include <utils.hpp>
#include <libgen.h>

#define EXAGOSelfcheckTCOPFLOWNumAnswers 1

/** Array of available reference solutions for TCOPFLOW selfcheck
 * See test/functionality/utils/selfcheck.h for selcheck solution structure
 */
static const ExaGOSelfcheckTCOPFLOWAnswer ExaGOSelfcheckTCOPFLOWAnswers[EXAGOSelfcheckTCOPFLOWNumAnswers] = 
{
  {SELFCHECK_NETWORK_CASE9_WIND, TCOPFLOW_INITIALIZATION, TCOPFLOW_GENBUSVOLTAGE, SELFCHECK_SCENARIO_CASE9, SELFCHECK_QLOAD, SELFCHECK_PLOAD, PETSC_FALSE, 5.0, 0.5, PETSC_FALSE, 18, 2.0792822464162207e+04}
};

PetscErrorCode ExaGOSelfcheckTCOPFLOWFindAnswer(TCOPFLOW tcopflow, ExaGOSelfcheckTCOPFLOWAnswer *ans)
{
  PetscErrorCode ierr;
  PetscBool flg;

  char netfile[PETSC_MAX_PATH_LEN];
  char windfile[PETSC_MAX_PATH_LEN];
  char qloadfile[PETSC_MAX_PATH_LEN];
  char ploadfile[PETSC_MAX_PATH_LEN];

  // Verify that all the correct options are being passed

  // Check for network file
  ierr = PetscOptionsGetString(NULL,NULL,"-netfile",netfile,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
  if (!flg)
  {
    ExaGOLog(EXAGO_LOG_WARN,"%s needs option -netfile.",__func__);
    return PETSC_TRUE;
  }
  strcpy(ans->networkname,basename(netfile));

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

  // Check for windfile
  ierr = PetscOptionsGetString(NULL,NULL,"-tcopflow_windgenprofile",windfile,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
  if (!flg)
  {
    ExaGOLog(EXAGO_LOG_WARN,"%s needs option -tcopflow_windgenprofile.",__func__);
    return PETSC_TRUE;
  }
  strcpy(ans->windgenname,basename(windfile));

  // Check for the pload profile
  ierr = PetscOptionsGetString(NULL,NULL,"-tcopflow_ploadprofile",ploadfile,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
  if (!flg)
  {
    ExaGOLog(EXAGO_LOG_WARN,"%s needs option -tcopflow_ploadprofile .",__func__);
    return PETSC_TRUE;
  }
  strcpy(ans->ploadprofile,basename(ploadfile));

  // Check for the qload profile
  ierr = PetscOptionsGetString(NULL,NULL,"-tcopflow_qloadprofile",qloadfile,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
  if (!flg)
  {
    ExaGOLog(EXAGO_LOG_WARN,"%s needs option -tcopflow_qloadprofile.",__func__);
    return PETSC_TRUE;
  }
  strcpy(ans->qloadprofile,basename(qloadfile));

  // Check for iscoupling
  ierr = PetscOptionsGetBool(NULL,NULL,"-tcopflow_iscoupling",&ans->iscoupling,&flg);CHKERRQ(ierr);
  if (!flg)
  {
    ExaGOLog(EXAGO_LOG_WARN,"%s needs option -tcopflow_iscoupling.",__func__);
    return PETSC_TRUE;
  }

  // Check for dt
  ierr = PetscOptionsGetReal(NULL,NULL,"-tcopflow_dT",&ans->dt,&flg);CHKERRQ(ierr);
  if (!flg)
  {
    ExaGOLog(EXAGO_LOG_WARN,"%s needs option -tcopflow_dT with -sopflow_enable_multiperiod.",__func__);
    return PETSC_TRUE;
  }

  // Check for duration
  ierr = PetscOptionsGetReal(NULL,NULL,"-tcopflow_duration",&ans->duration,&flg);CHKERRQ(ierr);
  if (!flg)
  {
    ExaGOLog(EXAGO_LOG_WARN,"%s needs option -tcopflow_duration with -sopflow_enable_multiperiod.",__func__);
    return PETSC_TRUE;
  }

  // Check for opflow_ignore_lineflow_constraints
  ierr = PetscOptionsGetBool(NULL,NULL,"-opflow_ignore_lineflow_constraints",&ans->lineflow_constraints,&flg);CHKERRQ(ierr);
  if (!flg)
  {
    ExaGOLog(EXAGO_LOG_WARN,"%s needs option -opflow_ignore_lineflow_constraints.",__func__);
    return PETSC_TRUE;
  }

  int i;
  for (int i = 0; i < EXAGOSelfcheckTCOPFLOWNumAnswers; i++)
  {
    const ExaGOSelfcheckTCOPFLOWAnswer *ians = &ExaGOSelfcheckTCOPFLOWAnswers[i];
    if (strcmp(ans->networkname,ians->networkname)==0 &&
        strcmp(ans->modelinit,ians->modelinit) == 0 &&
        strcmp(ans->genbusvoltage,ians->genbusvoltage) == 0 &&
        strcmp(ans->windgenname,ians->windgenname)==0 &&
        strcmp(ans->qloadprofile,ians->qloadprofile)==0 &&
        strcmp(ans->ploadprofile,ians->ploadprofile)==0 &&
        ans->iscoupling == ians->iscoupling &&
        ans->lineflow_constraints == ians->lineflow_constraints &&
        PetscEqualReal(ans->dt, ians->dt) &&
        PetscEqualReal(ans->duration, ians->duration))
    {
      ans->numiter = ians->numiter;
      ans->objective = ians->objective;
      return PETSC_FALSE;
    }
  }

  ExaGOLog(EXAGO_LOG_WARN,"%s could not find answer for option:"
      "\n\tnetwork:%s\n\twindgen:%s\n\tinitialization:%s\n\tgenbusvoltage:%s\
       \n\tiscoupling:%d\n\tlineflow_constraints:%d\
       \n\tpload:%s\n\tqload:%s\n\tdT:%f\n\tduration:%f\n",
      __func__,ans->networkname,ans->modelinit,ans->genbusvoltage,ans->windgenname,\
               ans->iscoupling, ans->lineflow_constraints,\
               ans->ploadprofile, ans->qloadprofile, ans->dt, ans->duration);
  return PETSC_TRUE;

}

PetscErrorCode ExaGOSelfcheckTCOPFLOW(TCOPFLOW tcopflow)
{
  PetscErrorCode ierr;
  PetscInt numiter;
  double obj_value,tolerance,error;
  int fail = 0;

  ExaGOSelfcheckTCOPFLOWAnswer *ans;
  ans = (ExaGOSelfcheckTCOPFLOWAnswer*)malloc(sizeof(ExaGOSelfcheckTCOPFLOWAnswer));

  ierr = ExaGOSelfcheckTCOPFLOWFindAnswer(tcopflow, ans);
  if (ierr)
  {
    free(ans);
    return ierr;
  }

  // Printing config makes greping for solutions/results easier
  char config[PETSC_MAX_PATH_LEN];
  sprintf(config,"%s+%s+%s+%s+%s+%s+%d+%d+%d+%d",\
          ans->networkname,ans->modelinit,ans->genbusvoltage,ans->windgenname,\
          ans->qloadprofile,ans->ploadprofile,\
          ans->iscoupling,ans->dt,ans->duration,ans->lineflow_constraints);

  // Check for convergence
  PetscBool conv_status=PETSC_FALSE;
  ierr = TCOPFLOWGetConvergenceStatus(tcopflow,&conv_status);CHKERRQ(ierr);
  if(conv_status == PETSC_FALSE)
  {
    fail++;
    ExaGOLog(EXAGO_LOG_ERROR,"ExaGO TCOPFLOW %s functionality test failed to converge.\n",config);
  }

  // Compare objective value
  ierr = TCOPFLOWGetObjective(tcopflow,&obj_value);CHKERRQ(ierr);
  ierr = TCOPFLOWGetTolerance(tcopflow, &tolerance);CHKERRQ(ierr);
  if (!isEqual(obj_value,ans->objective,tolerance,&error))
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
  ierr = TCOPFLOWGetNumIterations(tcopflow, &numiter);CHKERRQ(ierr);
  if (ans->numiter != numiter)
  {
    fail++;
    ExaGOLog(EXAGO_LOG_ERROR,
        "ExaGO TCOPFLOW %s functionality test failed to match expected number of iterations:\n"
        "expected=%d actual=%d",
        config,ans->numiter,numiter);
  }

  
  ExaGOLog(EXAGO_LOG_INFO,"%s: selfcheck returning %s.",
      __func__,fail==0?"success":"failure");

  free(ans);
  return fail;
}
