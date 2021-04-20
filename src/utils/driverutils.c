#include <common.h>
#include <utils.h>
#include <petscsys.h>
#include <private/utilsimpl.h>
#include <exago_config.h>
#include <version.h>

const char* ExaGOVerbosityNames[] = {"INFO","WARNING","ERROR"};

static char ExaGOCurrentAppName[128];

PetscErrorCode ExagoHelpPrintf(MPI_Comm comm, const char appname[],...)
{
  PetscErrorCode ierr;

  char* versionstr;
  ierr=ExaGOVersionGetFullVersionInfo(&versionstr);CHKERRQ(ierr);
  fprintf(stderr, "===================================================================\n");
  fprintf(stderr, "ExaGO Version Info:\n\n%s", versionstr); 
  if (strcmp(appname,"opflow") == 0)
  {
    fprintf(stderr, "============== Help Options for application %s ==============\n\n", appname); 
    fprintf(stderr, " General usage: mpiexec -n <N> ./%s <options>\n", appname);
    fprintf(stderr, " Options:\n");
    fprintf(stderr, "\t -netfile <netfilename>\n");
    fprintf(stderr, "\t -opflow_model <POWER_BALANCE_POLAR|...>\n");
    fprintf(stderr, "\t -opflow_solver <IPOPT|...>\n");
    fprintf(stderr, "\t -opflow_initialization <MIDPOINT|...>\n");
    fprintf(stderr, "\t -opflow_ignore_lineflow_constraints <0|1>\n");
    fprintf(stderr, "\t -opflow_include_loadloss_variables <0|1>\n");
    fprintf(stderr, "\t -opflow_include_powerimbalance_variables <0|1>\n");
    fprintf(stderr, "\t -opflow_loadloss_penalty <Penalty ($)>\n");
    fprintf(stderr, "\t -opflow_powerimbalance_penalty <Penalty ($)>\n");
    fprintf(stderr, "\t -opflow_genbusvoltage <FIXED_WITHIN_QBOUNDS|...>\n");
    fprintf(stderr, "\t -opflow_has_gensetpoint <0|1>\n");
    fprintf(stderr, "\t -opflow_objective <MIN_GEN_COST|...>\n");
    fprintf(stderr, "\t -opflow_use_agc <0|1>\n");
    fprintf(stderr, "\t -opflow_tolerance <1e-6|...>\n");
    fprintf(stderr, "\t -hiop_compute_mode <hybrid|...>\n");
    fprintf(stderr, "\t -hiop_verbosity_level <0-10>\n");
    fprintf(stderr, "\t -hiop_tolerance <1e-6|...>\n");
    fprintf(stderr, "\t -print_output <0|1>\n");
    fprintf(stderr, "\t -save_output <0|1>\n");
    fprintf(stderr, "\n");
  } 
  else if (strcmp(appname,"pflow")==0) 
  {
    fprintf(stderr, "============== Help Options for application %s ==============\n", appname);
  } 
  else if (strcmp(appname,"sopflow") == 0) 
  {
    fprintf(stderr, "============== Help Options for application %s ==============\n", appname);
    fprintf(stderr, " General usage: mpiexec -n <N> ./%s <options>\n", appname);
    fprintf(stderr, " Options:\n");
    fprintf(stderr, "\t -netfile <netfilename>\n");
    fprintf(stderr, "\t -ctgcfile <ctgcfilename>\n");
    fprintf(stderr, "\t -scenfile <scenfilename>\n");
    fprintf(stderr, "\t -opflow_model <POWER_BALANCE_POLAR|...>\n");
    fprintf(stderr, "\t -sopflow_solver <IPOPT|...>\n");
    fprintf(stderr, "\t -sopflow_mode <0|1>\n");
    fprintf(stderr, "\t -sopflow_Ns <Ns>\n");
    fprintf(stderr, "\t -sopflow_enable_multicontingency <0|1>\n");
    fprintf(stderr, "\t -scopflow_enable_multiperiod <0|1>\n");
    fprintf(stderr, "\t -sopflow_tolerance <1e-6|...>\n");
    fprintf(stderr, "\n");
  }
  else if (strcmp(appname,"scopflow") == 0) 
  {
    fprintf(stderr, "============== Help Options for application %s ==============\n", appname);
    fprintf(stderr, " General usage: mpiexec -n <N> ./%s <options>\n", appname);
    fprintf(stderr, " Options:\n");
    fprintf(stderr, "\t -netfile <netfilename>\n");
    fprintf(stderr, "\t -ctgcfile <ctgcfilename>\n");
    fprintf(stderr, "\t -scopflow_windgenprofile <windgenproffilename>\n");
    fprintf(stderr, "\t -scopflow_ploadprofile <ploadprofile_filename>\n");
    fprintf(stderr, "\t -scopflow_qloadprofile <qloadprofile_filename>\n");
    fprintf(stderr, "\t -opflow_model <POWER_BALANCE_POLAR|...>\n");
    fprintf(stderr, "\t -scopflow_solver <IPOPT|...>\n");
    fprintf(stderr, "\t -scopflow_mode <0|1>\n");
    fprintf(stderr, "\t -sopflow_Nc <Nc>\n");
    fprintf(stderr, "\t -scopflow_enable_multiperiod <0|1>\n");
    fprintf(stderr, "\t -scopflow_duration <Hours>\n");
    fprintf(stderr, "\t -scopflow_dT <Minutes>\n");
    fprintf(stderr, "\t -scopflow_tolerance <1e-6|...>\n");
    fprintf(stderr, "\n");

  }
  else if (strcmp(appname,"tcopflow") == 0) 
  {
    fprintf(stderr, "============== Help Options for application %s ==============\n", appname);
    fprintf(stderr, " General usage: mpiexec -n <N> ./%s <options>\n", appname);
    fprintf(stderr, " Options:\n");
    fprintf(stderr, "\t -netfile <netfilename>\n");
    fprintf(stderr, "\t -options_file <tcopflowoptionsfilename>\n");
    fprintf(stderr, "\t -tcopflow_windgenprofile <windgenproffilename>\n");
    fprintf(stderr, "\t -tcopflow_ploadprofile <ploadprofile_filename>\n");
    fprintf(stderr, "\t -tcopflow_qloadprofile <qloadprofile_filename>\n");
    fprintf(stderr, "\t -tcopflow_iscoupling <0|1>\n");
    fprintf(stderr, "\t -opflow_ignore_lineflow_constraints <0|1>\n");
    fprintf(stderr, "\t -tcopflow_duration <Hours>\n");
    fprintf(stderr, "\t -tcopflow_dT <Minutes>\n");
    fprintf(stderr, "\t -tcopflow_tolerance <1e-6|...>\n");
    fprintf(stderr, "\n");
  }
  else 
  {
    fprintf(stderr, "Please enter a valid application name.\n");
  }
	exit(0);
}

PetscErrorCode ExaGOInitialize(MPI_Comm comm, int *argc, char ***argv,
    char* appname, char *help)
{
  int i;
  PetscErrorCode ierr;

  strcpy(ExaGOCurrentAppName,appname);
  for(i=1; i<*argc; i++)
  {
    // Added null terminator here as this would cause errors with
    if ((strncmp((*argv)[i],"-help\0",5)==0) || 
        (strncmp((*argv)[i],"-h\0",3)==0))
    {
      PetscHelpPrintf = &ExagoHelpPrintf;
      (*PetscHelpPrintf)(MPI_COMM_NULL, appname);
    }
  }

  int initialized;
  MPI_Initialized(&initialized);
  if(!initialized)
    MPI_Init(argc,argv);

  ExaGOLogSetComm(comm);
  ierr = ExaGOLogInitialize();
  if(ierr)
  {
    fprintf(stderr, "Could not initialize ExaGO logger.\n");
    return ierr;
  }

  FILE *fp;
  ExaGOLogGetLoggingFilePointer(&fp);
  if(fp==NULL)
    ExaGOLogSetLoggingFilePointer(stderr);

  char options_pathname[200] = EXAGO_OPTIONS_DIR;
  char filename[64];
  ierr=sprintf(filename,"%soptions",ExaGOCurrentAppName);
  if(ierr<0)
  {
    fprintf(stderr,
        "Got error setting options path for ExaGO application %s.\n",
        ExaGOCurrentAppName);
    return 1;
  }

  const int noptfiles = 4;
  char **optfiles = (char**)malloc(sizeof(char*)*noptfiles);
  for(i=0; i<noptfiles; i++)
    optfiles[i] = (char*)malloc(PETSC_MAX_PATH_LEN);


  // Get -option_file option from command line
  char options_file[PETSC_MAX_PATH_LEN];
  PetscBool flg = PETSC_FALSE;
  optfiles[0] = NULL;
  for(i=1; i<*argc; i++)
  {
    // Added null terminator here as this would cause errors with
    if (0==strncmp((*argv)[i],"-options_file\0",13))
    {
      strcpy(optfiles[0],(*argv)[i + 1]);
      flg = PETSC_TRUE;
      break;
    }
  } 
  if (!flg)
  {
    ExaGOLog(EXAGO_LOG_INFO, "%s", "-options_file not passed.");
  }

  sprintf(optfiles[1],"%s/%s",options_pathname,filename);

  ierr=sprintf(optfiles[2],"./%s",filename);
  if(ierr<0)
  {
    fprintf(stderr,
        "Got error setting options path for ExaGO application %s.\n",
        ExaGOCurrentAppName);
    return 1;
  }

  ierr=sprintf(optfiles[3],"./options/%s",filename);
  if(ierr<0)
  {
    fprintf(stderr,
        "Got error setting options path for ExaGO application %s.\n",
        ExaGOCurrentAppName);
    return 1;
  }

  for(i=1; i<*argc; i++)
  {
    // Added null terminator here as this would cause errors with
    if (0==strncmp((*argv)[i],"-v\0",3) || 
        0==strncmp((*argv)[i],"-version\0",4))
    {
      char* versionstr;
      ierr=ExaGOVersionGetFullVersionInfo(&versionstr);CHKERRQ(ierr);
      ExaGOLog(EXAGO_LOG_INFO, "ExaGO Version Info:\n\n%s", versionstr);
      ExaGOLog(EXAGO_LOG_INFO, "ExaGO Help Message:\n%s", help);
      free(versionstr);
      for(i=0; i<noptfiles; i++)
        free(optfiles[i]);
      MPI_Finalize();
      exit(0);
    }
  }

  const int idx = anyFileExist(optfiles, noptfiles);
  if (idx < 0)
  {
    ExaGOLog(EXAGO_LOG_ERROR,
        "Could not find options file for application %s.\n",
        ExaGOCurrentAppName);
    for(i=0; i<noptfiles; i++)
      free(optfiles[i]);
    free(optfiles);
    return 1;
  }
  
  PetscInitialize(argc,argv,optfiles[idx],help);
  
  for(i=0; i<noptfiles; i++)
    free(optfiles[i]);
  free(optfiles);
  return 0;
}


PetscErrorCode ExaGOFinalize()
{
  ExaGOLog(EXAGO_LOG_INFO,"Finalizing %s application.\n",ExaGOCurrentAppName);
  PetscBool flg;
  ExaGOLogIsUsingLogFile(&flg);
  if (flg)
  {
    char* filename;
    ExaGOLogGetLoggingFileName(&filename);
    ExaGOLog(EXAGO_LOG_INFO,"See logfile %s for output.\n",filename);
  }
  ExaGOLogFinalize();
  PetscFinalize();
  MPI_Finalize();
}
