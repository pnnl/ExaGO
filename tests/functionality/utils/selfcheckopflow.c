#include <opflow.h>
#include <selfcheck.h>
#include <utils.h>
#include <libgen.h>

#define ExaGOSelfcheckOPFLOWNumAnswers 54

/**
 * Array of available reference solutions for OPFLOW selfcheck.
 * See docs/web/opflow.md for more information on solver/model combinations.
 */
static const ExaGOSelfcheckOPFLOWAnswer ExaGOSelfcheckOPFLOWAnswers[ExaGOSelfcheckOPFLOWNumAnswers] =
{
  /* solver, model, network, number of iterations, unscaled objective value */

  // - HIOP solver
  // -- HIOP model
  {OPFLOWSOLVER_HIOP,OPFLOWMODEL_PBPOLHIOP,SELFCHECK_NETWORK_CASE9  ,24 ,4.144460534990e+03},
  {OPFLOWSOLVER_HIOP,OPFLOWMODEL_PBPOLHIOP,SELFCHECK_NETWORK_CASE118,344,1.296606940659e+05},
  {OPFLOWSOLVER_HIOP,OPFLOWMODEL_PBPOLHIOP,SELFCHECK_NETWORK_CASE200,476,2.755757097185e+04},
  // -- RAJAHIOP model
  {OPFLOWSOLVER_HIOP,OPFLOWMODEL_PBPOLRAJAHIOP,SELFCHECK_NETWORK_CASE9  ,24 ,4.144460534985e+03},
  {OPFLOWSOLVER_HIOP,OPFLOWMODEL_PBPOLRAJAHIOP,SELFCHECK_NETWORK_CASE118,344,1.296606940659e+05},
  {OPFLOWSOLVER_HIOP,OPFLOWMODEL_PBPOLRAJAHIOP,SELFCHECK_NETWORK_CASE200,476,2.755757096954e+04},
  // - end HIOP solver

  // - TAO solver
  // -- HIOP model
  {OPFLOWSOLVER_TAO,OPFLOWMODEL_PBPOLHIOP,SELFCHECK_NETWORK_CASE9  ,9  ,4.144460556432e+03},
  {OPFLOWSOLVER_TAO,OPFLOWMODEL_PBPOLHIOP,SELFCHECK_NETWORK_CASE118,151,1.296606942608e+05},
  {OPFLOWSOLVER_TAO,OPFLOWMODEL_PBPOLHIOP,SELFCHECK_NETWORK_CASE200,28 ,2.755757088778e+04},
  // -- PBCAR model
  {OPFLOWSOLVER_TAO,OPFLOWMODEL_PBCAR,SELFCHECK_NETWORK_CASE9  ,12,4.144460546094e+03},
  {OPFLOWSOLVER_TAO,OPFLOWMODEL_PBCAR,SELFCHECK_NETWORK_CASE118,25,1.296606940832e+05},
  {OPFLOWSOLVER_TAO,OPFLOWMODEL_PBCAR,SELFCHECK_NETWORK_CASE200,30,2.755757088115e+04},
  // - end TAO solver

  // - IPOPT solver
  // -- PBPOL model
  {OPFLOWSOLVER_IPOPT,OPFLOWMODEL_PBPOL,SELFCHECK_NETWORK_CASE9  ,26 ,4.144460533813e+03},
  {OPFLOWSOLVER_IPOPT,OPFLOWMODEL_PBPOL,SELFCHECK_NETWORK_CASE118,17 ,1.296606940096e+05},
  {OPFLOWSOLVER_IPOPT,OPFLOWMODEL_PBPOL,SELFCHECK_NETWORK_CASE200,214,2.755757041954e+04},
  // -- PBCAR model
  {OPFLOWSOLVER_IPOPT,OPFLOWMODEL_PBCAR,SELFCHECK_NETWORK_CASE9  ,21,4.144460534632e+03},
  {OPFLOWSOLVER_IPOPT,OPFLOWMODEL_PBCAR,SELFCHECK_NETWORK_CASE118,15,1.296606940412e+05},
  {OPFLOWSOLVER_IPOPT,OPFLOWMODEL_PBCAR,SELFCHECK_NETWORK_CASE200,80,2.755757042019e+04},
  // -- IBCAR model
  {OPFLOWSOLVER_IPOPT,OPFLOWMODEL_IBCAR,SELFCHECK_NETWORK_CASE9  ,24 ,4.144460534632e+03},
  {OPFLOWSOLVER_IPOPT,OPFLOWMODEL_IBCAR,SELFCHECK_NETWORK_CASE118,29 ,1.296606940412e+05},
  {OPFLOWSOLVER_IPOPT,OPFLOWMODEL_IBCAR,SELFCHECK_NETWORK_CASE200,256,2.755757042019e+04},
  // -- IBCAR2 model
  {OPFLOWSOLVER_IPOPT,OPFLOWMODEL_IBCAR2,SELFCHECK_NETWORK_CASE9  ,20 ,4.144460534632e+03},
  {OPFLOWSOLVER_IPOPT,OPFLOWMODEL_IBCAR2,SELFCHECK_NETWORK_CASE118,36 ,1.296606940412e+05},
  {OPFLOWSOLVER_IPOPT,OPFLOWMODEL_IBCAR2,SELFCHECK_NETWORK_CASE200,121,2.755757042019e+04},
  // - end IPOPT solver
};

/** Temporary tolerance for error for number of iterations */
static int itertol = 2;

PetscBool ExaGOSelfcheckOPFLOWFindAnswer(OPFLOW opflow,ExaGOSelfcheckOPFLOWAnswer *ans)
{
  PetscErrorCode ierr;
  PetscBool flg=PETSC_FALSE;
  char netfile[PETSC_MAX_PATH_LEN];
  char solver[PETSC_MAX_PATH_LEN];
  char model[PETSC_MAX_PATH_LEN];

  ierr = PetscOptionsGetString(NULL,NULL,"-netfile",netfile,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
  if (!flg)
  {
    ExaGOLog(EXAGO_LOG_WARN,"%s needs option -netfile.",__func__);
    return PETSC_TRUE;
  }
  strcpy(ans->networkname,basename(netfile));

  ierr = PetscOptionsGetString(NULL,NULL,"-opflow_model",model,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
  if (!flg)
  {
    ExaGOLog(EXAGO_LOG_WARN,"%s needs option -opflow_model.",__func__);
    return PETSC_TRUE;
  }
  strcpy(ans->model,model);

  ierr = PetscOptionsGetString(NULL,NULL,"-opflow_solver",solver,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
  if (!flg)
  {
    ExaGOLog(EXAGO_LOG_WARN,"%s needs option -opflow_solver.",__func__);
    return PETSC_TRUE;
  }
  strcpy(ans->solver,solver);

  int i;
  for (i=0;i<ExaGOSelfcheckOPFLOWNumAnswers;i++)
  {
    const ExaGOSelfcheckOPFLOWAnswer *ians = &ExaGOSelfcheckOPFLOWAnswers[i];
    if (strcmp(ans->model,ians->model)==0
        && strcmp(ans->solver,ians->solver)==0
        && strcmp(ans->networkname,ians->networkname)==0)
    {
      ans->numiter = ians->numiter;
      ans->objective = ians->objective;
      return PETSC_FALSE;
    }
  }
  ExaGOLog(EXAGO_LOG_WARN,"%s could not find answer for options:"
      "\n\tsolver:%s\n\tmodel:%s\n\tnetwork:%s",
      __func__,ans->solver,ans->model,ans->networkname);
  return PETSC_TRUE;
}

PetscBool ExaGOSelfcheckOPFLOW(OPFLOW opflow)
{
  PetscErrorCode ierr;
  Vec X;
  double obj_value,tolerance,error;
  int numiter = 0;
  int fail=0,ierror=0;

  ExaGOSelfcheckOPFLOWAnswer *ans;
  ans = (ExaGOSelfcheckOPFLOWAnswer*)malloc(sizeof(ExaGOSelfcheckOPFLOWAnswer));

  /** Ensure answer for given opflow configuration exists */
  ierr = ExaGOSelfcheckOPFLOWFindAnswer(opflow,ans);
  if (ierr)
    return ierr;

  /** Printing config makes greping for solutions/results easier */
  char config[PETSC_MAX_PATH_LEN];
  sprintf(config,"%s+%s+%s",ans->solver,ans->model,ans->networkname);

  /** If OPFLOW didn't converge, we know it's failed */
  PetscBool conv_status=PETSC_FALSE;
  ierr = OPFLOWGetConvergenceStatus(opflow,&conv_status);CHKERRQ(ierr);
  if (conv_status==PETSC_FALSE)
  {
    fail++;
    ExaGOLog(EXAGO_LOG_ERROR,"ExaGO OPFLOW %s functionality test failed to converge.\n",config);
  }

  /** Ensure objective value for opflow matches reference solution */
  ierr = OPFLOWGetObjective(opflow,&obj_value);CHKERRQ(ierr);
  ierr = OPFLOWGetTolerance(opflow,&tolerance);CHKERRQ(ierr);
  if (!isEqual(obj_value,ans->objective,tolerance,&error))
  {
    fail++;
    ExaGOLog(EXAGO_LOG_ERROR,
        "ExaGO OPFLOW %s functionality test failed to match expected objective value:\n"
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

  /** Ensure number of iterations for opflow matches reference solution */
  ierr = OPFLOWGetNumIterations(opflow,&numiter);CHKERRQ(ierr);
  if (numiter!=ans->numiter)
  {
    fail++;
    ExaGOLog(EXAGO_LOG_ERROR,
        "ExaGO OPFLOW %s functionality test failed to match expected number of iterations:\n"
        "expected=%d actual=%d",
        config,ans->numiter,numiter);
  }

  /** Log results of functionality test */
  ExaGOLog(EXAGO_LOG_INFO,"%s: functionality test returning %s.",
      __func__,fail==0?"success":"failure");

  free(ans);
  return fail;
}
