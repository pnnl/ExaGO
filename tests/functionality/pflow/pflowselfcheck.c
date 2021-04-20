#include<pflowselfcheck.h>
#include<pflow.h>
#include<utils.h>
#include<libgen.h>

#define ExaGOSelfcheckPFLOWNumAnswers 3

/** Array of available reference solutions for PFLOW selfcheck
 * See test/functionality/utils/selfcheck.h for selcheck solution structure
 */
static const ExaGOSelfcheckPFLOWAnswer ExaGOSelfcheckPFLOWAnswers[ExaGOSelfcheckPFLOWNumAnswers] =
{
  /* Network, number of iterations */
  {SELFCHECK_NETWORK_CASE9, 3},
  {SELFCHECK_NETWORK_CASE118, 3},
  {SELFCHECK_NETWORK_CASE200, 2}
};

PetscErrorCode ExaGOSelfcheckPFLOWFindAnswer(PFLOW pflow,ExaGOSelfcheckPFLOWAnswer *ans)
{
  PetscErrorCode ierr;
  PetscBool flg;
  char netfile[PETSC_MAX_PATH_LEN];

  // Check that a network is being passed
  ierr = PetscOptionsGetString(NULL,NULL,"-netfile",netfile,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
  if (!flg)
  {
    ExaGOLog(EXAGO_LOG_WARN,"%s needs option -netfile.",__func__);
    return PETSC_TRUE;
  }
  strcpy(ans->networkname,basename(netfile));

  int i;
  for (i=0;i<ExaGOSelfcheckPFLOWNumAnswers;i++)
  {
    const ExaGOSelfcheckPFLOWAnswer *ians = &ExaGOSelfcheckPFLOWAnswers[i];
    if (strcmp(ans->networkname,ians->networkname)==0)
    {
      ans->numiter = ians->numiter;
      return PETSC_FALSE;
    }
  }
  ExaGOLog(EXAGO_LOG_WARN,"%s could not find answer for option:"
      "\n\tnetwork:%s",
      __func__,ans->networkname);
  return PETSC_FALSE;
}

PetscErrorCode ExaGOSelfcheckPFLOW(PFLOW pflow)
{
  PetscErrorCode ierr;
  PetscInt numiter;
  int fail = 0;

  ExaGOSelfcheckPFLOWAnswer *ans;
  ans = (ExaGOSelfcheckPFLOWAnswer*)malloc(sizeof(ExaGOSelfcheckPFLOWAnswer));

  ierr = ExaGOSelfcheckPFLOWFindAnswer(pflow, ans);
  if (ierr)
  {
    free(ans);
    return ierr;
  }
    
  /* Get the current number of iterations */
  ierr = PFLOWGetNumIterations(pflow, &numiter);CHKERRQ(ierr);
  if (numiter != ans->numiter)
    fail++;

  // Printing config makes greping for solutions/results easier
  char config[PETSC_MAX_PATH_LEN];
  sprintf(config,"%s",ans->networkname);

  ExaGOLog(EXAGO_LOG_INFO,"%s numiter: expected(%d) actual(%d)",
      config,ans->numiter,numiter);
  ExaGOLog(EXAGO_LOG_INFO,"%s: selfcheck returning %s.",
      __func__,fail==0?"success":"failure");

  free(ans);
  return fail;
}
