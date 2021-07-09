#include <iostream>
#include <cstdio>
#include <string>

#include <private/opflowimpl.h>
#include <exago_config.h>
#include <utils.hpp>

#include "include/OpflowTests.hpp"
#include "TestAcopfUtils.hpp"


#if defined(EXAGO_ENABLE_RAJA)
#include <umpire/Allocator.hpp>
#include <umpire/ResourceManager.hpp>
#include <RAJA/RAJA.hpp>

#ifdef EXAGO_ENABLE_GPU
  using exago_raja_exec = RAJA::cuda_exec<128>;
  using exago_raja_reduce = RAJA::cuda_reduce;
  using exago_raja_atomic = RAJA::cuda_atomic;
  #define RAJA_LAMBDA [=] __device__
#else
  using exago_raja_exec = RAJA::omp_parallel_for_exec;
  using exago_raja_reduce = RAJA::omp_reduce;
  using exago_raja_atomic = RAJA::omp_atomic;
  #define RAJA_LAMBDA [=]
#endif

#endif

/**
 * @brief Converts an array xin in natural ordering to an array xout in
 * sparse-dense ordering
 */
void naturaltospdense(const double *xin,double *xout,int *idxn2sd_map,int nx)
{
  int i;

  for(i=0; i < nx; i++) {
    xout[idxn2sd_map[i]] = xin[i];
  }
}

/**
 * @brief Converts an array xin in sparse dense ordering to an array xout in
 * natural ordering
 */
void spdensetonatural(const double *xin,double *xout,int *idxn2sd_map,int nx)
{
  int i;

  for(i=0; i < nx; i++) {
    xout[i] = xin[idxn2sd_map[i]];
  }
}

/**
 * @brief Unit test driver for objective function 
 * @see opflow/OpflowTests.hpp for kernel tested by this driver
 *
 * You can pass two options to the objectiveAcopf executatable through the command line (implemented using PETSc options):
 *
 *    ~ -netfile <data_file> : Specifies the input data file to test against. Default value is `/<exago_dir>/datafiles/case9/case9mod.m`.
 *                             See directory datafiles for other potential inputs.
 *
 *    ~ -num_copies <number> : Specifies the number of replications of the network given through `-netfile`.
 *                             If this is not set properly, test may fail
 *      
 */
int main(int argc, char** argv)
{
  const bool     isTestOpflowModelPBPOL     = true;
#if not defined(EXAGO_ENABLE_RAJA)
  const bool     isTestOpflowModelPBPOLRAJAHIOP = false;
  const bool     isTestOpflowModelPBPOLHIOP = false;
#else
  const bool     isTestOpflowModelPBPOLRAJAHIOP = false;
  const bool     isTestOpflowModelPBPOLHIOP = false;
#endif

  PetscErrorCode ierr;
  PetscBool      flg;
  Vec            X, Lambda;
  int            fail=0;
  PetscLogStage  stages[0];
  double         obj_value;
  char           file_c_str[PETSC_MAX_PATH_LEN];
  std::string    file;
  char           appname[]="opflow";
  MPI_Comm       comm=MPI_COMM_WORLD;
  int            num_copies = 1;

  char help[] = "Unit tests for objective function running opflow\n";
  
  /** Use `ExaGOLogSetLoggingFileName("opflow-logfile");` to log the output. */
  ierr = ExaGOInitialize(comm,&argc,&argv,appname,help);
  if (ierr)
  {
    fprintf(stderr,"Could not initialize ExaGO application %s.\n",appname);
    return ierr;
  }

  /* Get network data file from command line */
  ierr = PetscOptionsGetString(NULL,NULL,"-netfile",file_c_str,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);

  /* Get network data file from command line */
  ierr = PetscOptionsGetInt(NULL,NULL,"-num_copies",&num_copies,&flg);CHKERRQ(ierr);


  if(!flg){
    file = "../datafiles/case9/case9mod.m";
  } 
  else{
    file.assign(file_c_str);
  }

  std::cout << file << std::endl;


  ierr = PetscLogStageRegister("Test stage",&stages[0]);CHKERRQ(ierr);

  // Set obj_value as reference solution, and run as usual
  obj_value = 10.0 * num_copies;

  ierr = PetscLogStagePush(stages[0]);CHKERRQ(ierr);

  if(isTestOpflowModelPBPOLHIOP)
  {
    OPFLOW opflowtest;
    exago::tests::TestOpflow test;

    std::cout << "\nTesting custom power balance model in polar coordinates for HIOP"
              << "(componentwise assembly) ... \n";

    // Create optimal power flow model
    ierr = OPFLOWCreate(PETSC_COMM_WORLD,&opflowtest);CHKERRQ(ierr);

    /* Read Network data */
    ierr = OPFLOWReadMatPowerData(opflowtest,file.c_str());CHKERRQ(ierr);

    /* Set opflow model type to power balance polar2 (component assembly) */
    ierr = OPFLOWSetModel(opflowtest,OPFLOWMODEL_PBPOLHIOP);CHKERRQ(ierr);

    /* Set solver to HIOP */
    ierr = OPFLOWSetSolver(opflowtest,OPFLOWSOLVER_HIOP);CHKERRQ(ierr);

    /* Set up */
    ierr = OPFLOWSetUp(opflowtest);CHKERRQ(ierr);
    ierr = OPFLOWGetSolution(opflowtest, &X);CHKERRQ(ierr);
    ierr = OPFLOWGetConstraintMultipliers(opflowtest,&Lambda);CHKERRQ(ierr);

    int nx,nconeq,nconineq,*idxn2sd_map;
    ierr = OPFLOWGetSizes(opflowtest,&nx,&nconeq,&nconineq);CHKERRQ(ierr);
    ierr = OPFLOWGetVariableOrdering(opflowtest,&idxn2sd_map);CHKERRQ(ierr);

    double *x_vec, *x_ref;

    ierr = PetscMalloc1(nx,&x_ref);CHKERRQ(ierr);

    ierr = VecGetArray(X,&x_vec);CHKERRQ(ierr);

    /* Convert from natural to sparse dense ordering */
    naturaltospdense(x_vec,x_ref,idxn2sd_map,nx);CHKERRQ(ierr);

    /* _ref pointers are now in sparse-dense ordering */

    fail += test.computeObjective(opflowtest,x_ref,obj_value);

    ierr = PetscFree(x_ref);CHKERRQ(ierr);

    ierr = VecRestoreArray(X,&x_vec);CHKERRQ(ierr);

    ierr = OPFLOWDestroy(&opflowtest);CHKERRQ(ierr);
  }

#if defined(EXAGO_ENABLE_RAJA)
  if(isTestOpflowModelPBPOLRAJAHIOP)
  {
    OPFLOW opflowtest;
    exago::tests::TestOpflow test;

    std::cout << "\nTesting custom power balance model in polar coordinates for HIOP using RAJA"
              << "(PBPOLHIOPRAJA) ... \n";

    // Create optimal power flow model
    ierr = OPFLOWCreate(PETSC_COMM_WORLD,&opflowtest);CHKERRQ(ierr);

    /* Read Network data */
    ierr = OPFLOWReadMatPowerData(opflowtest,file.c_str());CHKERRQ(ierr);

    /* Set opflow model type to custom model for hiop using RAJA */
    ierr = OPFLOWSetModel(opflowtest,OPFLOWMODEL_PBPOLRAJAHIOP);CHKERRQ(ierr);

    /* Set solver to HIOP */
    ierr = OPFLOWSetSolver(opflowtest,OPFLOWSOLVER_HIOP);CHKERRQ(ierr);

    /* Set up */
    ierr = OPFLOWSetUp(opflowtest);CHKERRQ(ierr);
    ierr = OPFLOWGetSolution(opflowtest, &X);CHKERRQ(ierr);
    ierr = OPFLOWGetConstraintMultipliers(opflowtest,&Lambda);CHKERRQ(ierr);

    int nx,nconeq,nconineq,*idxn2sd_map;
    ierr = OPFLOWGetSizes(opflowtest,&nx,&nconeq,&nconineq);CHKERRQ(ierr);
    ierr = OPFLOWGetVariableOrdering(opflowtest,&idxn2sd_map);CHKERRQ(ierr);

   double *x_vec, *x_ref;

    ierr = PetscMalloc1(nx,&x_ref);CHKERRQ(ierr);

    // Petsc Documentation mentions "You MUST call VecRestoreArray() when you no longer need access to the array."...
    ierr = VecGetArray(X,&x_vec);CHKERRQ(ierr);


    /* Convert from natural to sparse dense ordering */
    naturaltospdense(x_vec,x_ref,idxn2sd_map,nx);CHKERRQ(ierr);

    /* _ref pointers are now in sparse-dense ordering */
    
    // Get resource manager instance
    auto& resmgr = umpire::ResourceManager::getInstance();

    // Get Allocator
    umpire::Allocator h_allocator = resmgr.getAllocator("HOST");

    // Register array xref with umpire
    umpire::util::AllocationRecord record_x{x_ref,sizeof(double)*nx,h_allocator.getAllocationStrategy()};
    resmgr.registerAllocation(x_ref,record_x);

    // Allocate and copy xref to device
    double *x_ref_dev;
#ifdef EXAGO_ENABLE_GPU
    umpire::Allocator d_allocator = resmgr.getAllocator("DEVICE");
    x_ref_dev = static_cast<double*>(d_allocator.allocate(nx*sizeof(double)));
#else
    x_ref_dev = x_ref;
#endif

    resmgr.copy(x_ref_dev,x_ref);
					    
    // Tests
    fail += test.computeObjective(opflowtest,x_ref_dev,obj_value);

#ifdef EXAGO_ENABLE_GPU
    d_allocator.deallocate(x_ref_dev);
#endif

    ierr = PetscFree(x_ref);CHKERRQ(ierr);

    ierr = VecRestoreArray(X,&x_vec);CHKERRQ(ierr);

    ierr = OPFLOWDestroy(&opflowtest);CHKERRQ(ierr);
  }
#endif

  if (isTestOpflowModelPBPOL)
  {
    OPFLOW                   opflowtest;
    exago::tests::TestOpflow test;

    std::cout << "\nTesting power balance model in polar coordinates "
              << " ... \n";

    /* Set up test opflow */
    ierr = OPFLOWCreate(PETSC_COMM_WORLD,&opflowtest);CHKERRQ(ierr);
    ierr = OPFLOWReadMatPowerData(opflowtest,file.c_str());CHKERRQ(ierr);
    ierr = OPFLOWSetSolver(opflowtest,OPFLOWSOLVER_IPOPT);CHKERRQ(ierr);
    ierr = OPFLOWSetModel(opflowtest,OPFLOWMODEL_PBPOL);CHKERRQ(ierr);
    ierr = OPFLOWSetUp(opflowtest);CHKERRQ(ierr);
    ierr = OPFLOWGetSolution(opflowtest, &X);CHKERRQ(ierr);

    fail += test.computeObjective(opflowtest,X,obj_value);

    ierr = OPFLOWDestroy(&opflowtest);CHKERRQ(ierr);
  }

  ExaGOFinalize();
  return fail;
}
