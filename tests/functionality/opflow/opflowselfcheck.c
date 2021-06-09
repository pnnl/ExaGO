#include <opflow.h>
#include <opflowselfcheck.h>
#include <utils.hpp>
#include <libgen.h>

#define ExaGOSelfcheckOPFLOWNumAnswers 100

/**
 * Array of available reference solutions for OPFLOW selfcheck.
 * See docs/web/opflow.md for more information on solver/model combinations.
 * See test/functionality/utils/selfcheck.h for selcheck solution structure
 */
static const ExaGOSelfcheckOPFLOWAnswer ExaGOSelfcheckOPFLOWAnswers[ExaGOSelfcheckOPFLOWNumAnswers] =
{
  /* solver, model, network, isloadloss_active,ispowerimb_active,number of iterations, unscaled objective value */

  // - HIOP solver
  // -- HIOP model
  {OPFLOWSOLVER_HIOP,OPFLOWMODEL_PBPOLHIOP,SELFCHECK_NETWORK_CASE9  ,false,false,14,4.1444511e+03},
  {OPFLOWSOLVER_HIOP,OPFLOWMODEL_PBPOLHIOP,SELFCHECK_NETWORK_CASE118,false,false,25,1.2965990e+05},
  {OPFLOWSOLVER_HIOP,OPFLOWMODEL_PBPOLHIOP,SELFCHECK_NETWORK_CASE200,false,false,130,2.7552974e+04},
  // -- RAJAHIOP model
  {OPFLOWSOLVER_HIOP,OPFLOWMODEL_PBPOLRAJAHIOP,SELFCHECK_NETWORK_CASE9  ,false,false,14,4.1444511e+03},
  {OPFLOWSOLVER_HIOP,OPFLOWMODEL_PBPOLRAJAHIOP,SELFCHECK_NETWORK_CASE118,false,false,25,1.2965990e+05},
  {OPFLOWSOLVER_HIOP,OPFLOWMODEL_PBPOLRAJAHIOP,SELFCHECK_NETWORK_CASE200,false,false,140,2.7552974e+04},

  // -- HIOP model with load loss
  {OPFLOWSOLVER_HIOP,OPFLOWMODEL_PBPOLHIOP,SELFCHECK_NETWORK_CASE9  ,true,false,14,4.1444511e+03},
  {OPFLOWSOLVER_HIOP,OPFLOWMODEL_PBPOLHIOP,SELFCHECK_NETWORK_CASE118,true,false,26,1.2965990e+05},
  {OPFLOWSOLVER_HIOP,OPFLOWMODEL_PBPOLHIOP,SELFCHECK_NETWORK_CASE200,true,false,130,2.7552974e+04},
  // -- RAJAHIOP model with load loss
  {OPFLOWSOLVER_HIOP,OPFLOWMODEL_PBPOLRAJAHIOP,SELFCHECK_NETWORK_CASE9  ,true,false,14,4.1444511e+03},
  {OPFLOWSOLVER_HIOP,OPFLOWMODEL_PBPOLRAJAHIOP,SELFCHECK_NETWORK_CASE118,true,false,25,1.2965990e+05},
  {OPFLOWSOLVER_HIOP,OPFLOWMODEL_PBPOLRAJAHIOP,SELFCHECK_NETWORK_CASE200,true,false,140,2.7552974e+04},

  // -- HIOP model with power imbalance
  {OPFLOWSOLVER_HIOP,OPFLOWMODEL_PBPOLHIOP,SELFCHECK_NETWORK_CASE9  ,false,true,14,4.1444511e+03},
  {OPFLOWSOLVER_HIOP,OPFLOWMODEL_PBPOLHIOP,SELFCHECK_NETWORK_CASE118,false,true,25,1.2965990e+05},
  {OPFLOWSOLVER_HIOP,OPFLOWMODEL_PBPOLHIOP,SELFCHECK_NETWORK_CASE200,false,true,130,2.7552974e+04},
  // -- RAJAHIOP model with power imbalance
  {OPFLOWSOLVER_HIOP,OPFLOWMODEL_PBPOLRAJAHIOP,SELFCHECK_NETWORK_CASE9  ,false,true,14,4.1444511e+03},
  {OPFLOWSOLVER_HIOP,OPFLOWMODEL_PBPOLRAJAHIOP,SELFCHECK_NETWORK_CASE118,false,true,25,1.2965990e+05},
  {OPFLOWSOLVER_HIOP,OPFLOWMODEL_PBPOLRAJAHIOP,SELFCHECK_NETWORK_CASE200,false,true,140,2.7552974e+04},

  // -- HIOP model with load loss and power imbalance
  {OPFLOWSOLVER_HIOP,OPFLOWMODEL_PBPOLHIOP,SELFCHECK_NETWORK_CASE9  ,true,true,14,4.1444511e+03},
  {OPFLOWSOLVER_HIOP,OPFLOWMODEL_PBPOLHIOP,SELFCHECK_NETWORK_CASE118,true,true,26,1.2965990e+05},
  {OPFLOWSOLVER_HIOP,OPFLOWMODEL_PBPOLHIOP,SELFCHECK_NETWORK_CASE200,true,true,130,2.7552974e+04},
  // -- RAJAHIOP model
  {OPFLOWSOLVER_HIOP,OPFLOWMODEL_PBPOLRAJAHIOP,SELFCHECK_NETWORK_CASE9  ,true,true,12,4.1444511e+03},
  {OPFLOWSOLVER_HIOP,OPFLOWMODEL_PBPOLRAJAHIOP,SELFCHECK_NETWORK_CASE118,true,true,26,1.2965990e+05},
  {OPFLOWSOLVER_HIOP,OPFLOWMODEL_PBPOLRAJAHIOP,SELFCHECK_NETWORK_CASE200,true,true,140,2.7552974e+04},

  // - end HIOP solver

  // - HIOPSPARSE solver
  // -- PBPOL model
  {OPFLOWSOLVER_HIOPSPARSE,OPFLOWMODEL_PBPOL,SELFCHECK_NETWORK_CASE9  ,false,false,12 ,4.1444511e+03},
  {OPFLOWSOLVER_HIOPSPARSE,OPFLOWMODEL_PBPOL,SELFCHECK_NETWORK_CASE118,false,false,25 ,1.2965990e+05},
  {OPFLOWSOLVER_HIOPSPARSE,OPFLOWMODEL_PBPOL,SELFCHECK_NETWORK_CASE200,false,false,130,2.7552974e+04},
  // - end HIOPSPARSE solver


  // - TAO solver
  // -- HIOP model
  {OPFLOWSOLVER_TAO,OPFLOWMODEL_PBPOLHIOP,SELFCHECK_NETWORK_CASE9  ,false,false,9  ,4.144460556432e+03},
  {OPFLOWSOLVER_TAO,OPFLOWMODEL_PBPOLHIOP,SELFCHECK_NETWORK_CASE118,false,false,151,1.296606942608e+05},
  {OPFLOWSOLVER_TAO,OPFLOWMODEL_PBPOLHIOP,SELFCHECK_NETWORK_CASE200,false,false,28 ,2.755757088778e+04},
  // -- PBCAR model
  {OPFLOWSOLVER_TAO,OPFLOWMODEL_PBCAR,SELFCHECK_NETWORK_CASE9  ,false,false,12,4.144460546094e+03},
  {OPFLOWSOLVER_TAO,OPFLOWMODEL_PBCAR,SELFCHECK_NETWORK_CASE118,false,false,25,1.296606940832e+05},
  {OPFLOWSOLVER_TAO,OPFLOWMODEL_PBCAR,SELFCHECK_NETWORK_CASE200,false,false,30,2.755757088115e+04},
  // - end TAO solver

  // - IPOPT solver
  // -- PBPOL model
  {OPFLOWSOLVER_IPOPT,OPFLOWMODEL_PBPOL,SELFCHECK_NETWORK_CASE9  ,false,false,25 ,4.1513585954144e+03},
  {OPFLOWSOLVER_IPOPT,OPFLOWMODEL_PBPOL,SELFCHECK_NETWORK_CASE118,false,false,255,1.3013187197492e+05},
  {OPFLOWSOLVER_IPOPT,OPFLOWMODEL_PBPOL,SELFCHECK_NETWORK_CASE200,false,false,191,2.7564247365138e+04},
  // -- PBCAR model
  {OPFLOWSOLVER_IPOPT,OPFLOWMODEL_PBCAR,SELFCHECK_NETWORK_CASE9  ,false,false,21,4.144460534632e+03},
  {OPFLOWSOLVER_IPOPT,OPFLOWMODEL_PBCAR,SELFCHECK_NETWORK_CASE118,false,false,15,1.296606940412e+05},
  {OPFLOWSOLVER_IPOPT,OPFLOWMODEL_PBCAR,SELFCHECK_NETWORK_CASE200,false,false,80,2.755757042019e+04},
  // -- IBCAR model
  {OPFLOWSOLVER_IPOPT,OPFLOWMODEL_IBCAR,SELFCHECK_NETWORK_CASE9  ,false,false,24 ,4.144460534632e+03},
  {OPFLOWSOLVER_IPOPT,OPFLOWMODEL_IBCAR,SELFCHECK_NETWORK_CASE118,false,false,29 ,1.296606940412e+05},
  {OPFLOWSOLVER_IPOPT,OPFLOWMODEL_IBCAR,SELFCHECK_NETWORK_CASE200,false,false,256,2.755757042019e+04},
  // -- IBCAR2 model
  {OPFLOWSOLVER_IPOPT,OPFLOWMODEL_IBCAR2,SELFCHECK_NETWORK_CASE9  ,false,false,20 ,4.144460534632e+03},
  {OPFLOWSOLVER_IPOPT,OPFLOWMODEL_IBCAR2,SELFCHECK_NETWORK_CASE118,false,false,36 ,1.296606940412e+05},
  {OPFLOWSOLVER_IPOPT,OPFLOWMODEL_IBCAR2,SELFCHECK_NETWORK_CASE200,false,false,121,2.755757042019e+04},

  // -- PBPOL model with load loss
  {OPFLOWSOLVER_IPOPT,OPFLOWMODEL_PBPOL,SELFCHECK_NETWORK_CASE9  ,true,false,29 ,4.1513585954144e+03},
  {OPFLOWSOLVER_IPOPT,OPFLOWMODEL_PBPOL,SELFCHECK_NETWORK_CASE118,true,false,716,1.3013187197492e+05},
  {OPFLOWSOLVER_IPOPT,OPFLOWMODEL_PBPOL,SELFCHECK_NETWORK_CASE200,true,false,390,2.7564247365138e+04},

    // -- PBPOL model with power imbalance
  {OPFLOWSOLVER_IPOPT,OPFLOWMODEL_PBPOL,SELFCHECK_NETWORK_CASE9  ,false,true,25 ,4.1513585954144e+03},
  {OPFLOWSOLVER_IPOPT,OPFLOWMODEL_PBPOL,SELFCHECK_NETWORK_CASE118,false,true,166,1.3013187197492e+05},
  {OPFLOWSOLVER_IPOPT,OPFLOWMODEL_PBPOL,SELFCHECK_NETWORK_CASE200,false,true,245,2.7564247365138e+04},

    // -- PBPOL model with load loss and power imbalance
  {OPFLOWSOLVER_IPOPT,OPFLOWMODEL_PBPOL,SELFCHECK_NETWORK_CASE9  ,true,true,28 ,4.1513585954144e+03},
  {OPFLOWSOLVER_IPOPT,OPFLOWMODEL_PBPOL,SELFCHECK_NETWORK_CASE118,true,true,154,1.3013187197492e+05},
  {OPFLOWSOLVER_IPOPT,OPFLOWMODEL_PBPOL,SELFCHECK_NETWORK_CASE200,true,true,493,2.7564247365138e+04},

  // - end IPOPT solver
};

/** Temporary tolerance for error for number of iterations */
static int itertol = 2;

PetscErrorCode ExaGOSelfcheckOPFLOWFindAnswer(OPFLOW opflow,ExaGOSelfcheckOPFLOWAnswer *ans)
{
  PetscErrorCode ierr;
  PetscBool flg=PETSC_FALSE;
  char netfile[PETSC_MAX_PATH_LEN];
  char solver[PETSC_MAX_PATH_LEN];
  char model[PETSC_MAX_PATH_LEN];
  PetscBool isloadloss_active;
  PetscBool ispowerimb_active;

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

  ierr = PetscOptionsGetBool(NULL,NULL,"-opflow_include_loadloss_variables",&isloadloss_active,&flg);CHKERRQ(ierr);
  if (!flg)
  {
    ExaGOLog(EXAGO_LOG_WARN,"%s needs option -opflow_include_loadloss_variables.",__func__);
    return PETSC_TRUE;
  }
  ans->isloadloss_active = isloadloss_active;

  ierr = PetscOptionsGetBool(NULL,NULL,"-opflow_include_powerimbalance_variables",&ispowerimb_active,&flg);CHKERRQ(ierr);
  if (!flg)
  {
    ExaGOLog(EXAGO_LOG_WARN,"%s needs option -opflow_include_powerimbalance_variables.",__func__);
    return PETSC_TRUE;
  }
  ans->ispowerimb_active = ispowerimb_active;


  int i;
  for (i=0;i<ExaGOSelfcheckOPFLOWNumAnswers;i++)
  {
    const ExaGOSelfcheckOPFLOWAnswer *ians = &ExaGOSelfcheckOPFLOWAnswers[i];
    if ((strcmp(ans->model,ians->model)==0
        && strcmp(ans->solver,ians->solver)==0
        && strcmp(ans->networkname,ians->networkname)==0)
        && ans->isloadloss_active == ians->isloadloss_active
	&& ans->ispowerimb_active == ians->ispowerimb_active)
    {
      ans->numiter = ians->numiter;
      ans->objective = ians->objective;
      return PETSC_FALSE;
    }
  }
  ExaGOLog(EXAGO_LOG_WARN,"%s could not find answer for options:"
	   "\n\tsolver:%s\n\tmodel:%s\n\tnetwork:%s\tisloadloss_active=%d\tispowerimbalance_active=%d\n",__func__,ans->solver,ans->model,ans->networkname,ans->isloadloss_active,ans->ispowerimb_active);
  return PETSC_TRUE;
}

PetscErrorCode ExaGOSelfcheckOPFLOW(OPFLOW opflow)
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
  {
    free(ans);
    return ierr;
  }
  
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
  if (!IsEqual(obj_value,ans->objective,tolerance,error))
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
