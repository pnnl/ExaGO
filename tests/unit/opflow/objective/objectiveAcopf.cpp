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
 * @brief Unit test driver for ACOPF models
 * @see opflow/OpflowTests.hpp for kernels tested by this driver
 *
 * You can pass several options to the TestAcopf executatable through the command line (implemented using PETSc options):
 *
 *    ~ -netfile <data_file> : Specifies the input data file to test against. Default value is `/<exago_dir>/datafiles/case9/case9mod.m`.
 *                             See directory datafiles for other potential inputs.
 *
 *    ~ -gen_test_data       : If used, generates an answer key using IPOPT to test against.
 *                             If not used, uses existing answer keys located in `/<exago_dir>/datafiles/test_validation`.
 *
 *    ~ -write_test_data     : If used, generates test data, and then writes out the results to `/<install_dir>/tests/datafiles/test_validation/<data_file>/`.
 *                             In order to save generated results, you should copy them into `/<exago_dir>/datafiles/test_validation/<data_file>/` and add them to git.
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
  PetscBool      flg, gen_test_data, write_test_data;
  bool           ineq_present = false;
  /*
  Vec            X, Xl, Xu, G, Gl, Gu, grad, Lambda;
  */
  Vec            X, Lambda;
  // Mat            Jeq, Jineq, Hess;
  int            fail=0;
  PetscLogStage  stages[0];
  double         obj_value, obj_factor;
  char           file_c_str[PETSC_MAX_PATH_LEN];
  char           validation_dir_c_str[PETSC_MAX_PATH_LEN];
  std::string    file, validation_path;
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


  // Vec X = -
  // x_ref = -

  ierr = PetscLogStagePop();CHKERRQ(ierr);

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

    /*
    double *x_vec, *xl_vec, *xu_vec, *grad_vec, \
           *x_ref, *xl_ref, *xu_ref, *grad_ref, *g_ref, *gl_ref, *gu_ref, *lambda_ref;
           */
    double *x_vec, *x_ref, *lambda_ref;

    ierr = PetscMalloc1(nx,&x_ref);CHKERRQ(ierr);
    // ierr = PetscMalloc1(nx,&xl_ref);CHKERRQ(ierr);
    // ierr = PetscMalloc1(nx,&xu_ref);CHKERRQ(ierr);
    // ierr = PetscMalloc1(nx,&grad_ref);CHKERRQ(ierr);

    ierr = VecGetArray(X,&x_vec);CHKERRQ(ierr);
    // ierr = VecGetArray(Xl,&xl_vec);CHKERRQ(ierr);
    // ierr = VecGetArray(Xu,&xu_vec);CHKERRQ(ierr);
    // ierr = VecGetArray(grad,&grad_vec);CHKERRQ(ierr);


    /**
     * @note IPOPT passes the scaled Lagrangian multipliers for Hessian 
     * calculation. The scaling factor is the Hessian objective factor it 
     * uses in the Hessian calculation, so we need to scale Lambda. Note 
     * that obj_factor = 1.0 except when the solver is IPOPT which sets 
     * obj_factor
     */
    ierr = VecScale(Lambda,obj_factor);CHKERRQ(ierr);

    ierr = VecGetArray(Lambda, &lambda_ref);CHKERRQ(ierr);

    // ierr = VecGetArray(G,&g_ref);CHKERRQ(ierr);
    // ierr = VecGetArray(Gl,&gl_ref);CHKERRQ(ierr);
    // ierr = VecGetArray(Gu,&gu_ref);CHKERRQ(ierr);

    /* Convert from natural to sparse dense ordering */
    naturaltospdense(x_vec,x_ref,idxn2sd_map,nx);CHKERRQ(ierr);
    // naturaltospdense(xl_vec,xl_ref,idxn2sd_map,nx);CHKERRQ(ierr);
    // naturaltospdense(xu_vec,xu_ref,idxn2sd_map,nx);CHKERRQ(ierr);
    // naturaltospdense(grad_vec,grad_ref,idxn2sd_map,nx);CHKERRQ(ierr);
    /* _ref pointers are now in sparse-dense ordering */
    // fail += test.computeVariableBounds(opflowtest,xl_ref,xu_ref);
    //
    fail += test.computeObjective(opflowtest,x_ref,obj_value);
    // fail += test.computeGradient(opflowtest,x_ref,grad_ref);
    // fail += test.computeConstraints(opflowtest,x_ref,g_ref);
    // fail += test.computeConstraintBounds(opflowtest,gl_ref,gu_ref);
    // fail += test.computeConstraintJacobian(opflowtest, x_ref, Jeq, Jineq);
    // fail += test.computeHessian(opflowtest, x_ref, lambda_ref, obj_factor, Hess);

    ierr = PetscFree(x_ref);CHKERRQ(ierr);
    //  ierr = PetscFree(xl_ref);CHKERRQ(ierr);
    //  ierr = PetscFree(xu_ref);CHKERRQ(ierr);
    //  ierr = PetscFree(grad_ref);CHKERRQ(ierr);

    ierr = VecRestoreArray(X,&x_vec);CHKERRQ(ierr);
    // ierr = VecRestoreArray(Xl,&xl_vec);CHKERRQ(ierr);
    // ierr = VecRestoreArray(Xu,&xu_vec);CHKERRQ(ierr);
    // ierr = VecRestoreArray(grad,&grad_vec);CHKERRQ(ierr);
    // ierr = VecRestoreArray(G,&g_ref);CHKERRQ(ierr);
    // ierr = VecRestoreArray(Gl,&gl_ref);CHKERRQ(ierr);
    // ierr = VecRestoreArray(Gu,&gu_ref);CHKERRQ(ierr);
    ierr = VecRestoreArray(Lambda,&lambda_ref);CHKERRQ(ierr);

    ierr = VecScale(Lambda,1/obj_factor);CHKERRQ(ierr);
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

    std::cout << "nx = " << nx << std::endl  << "nconeq = " << nconeq << std::endl
              << "nconineq = " << nconineq << std::endl;

    //double *x_vec, *xl_vec, *xu_vec, *grad_vec, \
      *x_ref, *xl_ref, *xu_ref, *grad_ref, *g_ref, *gl_ref, *gu_ref,*lambda_ref;
   double *x_vec, *x_ref, *lambda_ref;

    ierr = PetscMalloc1(nx,&x_ref);CHKERRQ(ierr);
    // ierr = PetscMalloc1(nx,&xl_ref);CHKERRQ(ierr);
    // ierr = PetscMalloc1(nx,&xu_ref);CHKERRQ(ierr);
    // ierr = PetscMalloc1(nx,&grad_ref);CHKERRQ(ierr);

    // Petsc Documentation mentions "You MUST call VecRestoreArray() when you no longer need access to the array."...
    ierr = VecGetArray(X,&x_vec);CHKERRQ(ierr);
    /*
    ierr = VecGetArray(Xl,&xl_vec);CHKERRQ(ierr);
    ierr = VecGetArray(Xu,&xu_vec);CHKERRQ(ierr);
    ierr = VecGetArray(grad,&grad_vec);CHKERRQ(ierr);
    */

    /**
     * @note IPOPT passes the scaled Lagrangian multipliers for Hessian 
     * calculation. The scaling factor is the Hessian objective factor it 
     * uses in the Hessian calculation, so we need to scale Lambda. Note 
     * that obj_factor = 1.0 except when the solver is IPOPT which sets 
     * obj_factor
     */
    ierr = VecScale(Lambda,obj_factor);CHKERRQ(ierr);

    ierr = VecGetArray(Lambda,&lambda_ref);CHKERRQ(ierr);

    /*
    ierr = VecGetArray(G,&g_ref);CHKERRQ(ierr);
    ierr = VecGetArray(Gl,&gl_ref);CHKERRQ(ierr);
    ierr = VecGetArray(Gu,&gu_ref);CHKERRQ(ierr);
    */

    /* Convert from natural to sparse dense ordering */
    naturaltospdense(x_vec,x_ref,idxn2sd_map,nx);CHKERRQ(ierr);
    /*
    naturaltospdense(xl_vec,xl_ref,idxn2sd_map,nx);CHKERRQ(ierr);
    naturaltospdense(xu_vec,xu_ref,idxn2sd_map,nx);CHKERRQ(ierr);
    naturaltospdense(grad_vec,grad_ref,idxn2sd_map,nx);CHKERRQ(ierr);
    */
    /* _ref pointers are now in sparse-dense ordering */
    
    // Get resource manager instance
    auto& resmgr = umpire::ResourceManager::getInstance();

    // Get Allocator
    umpire::Allocator h_allocator = resmgr.getAllocator("HOST");

    // Register array xref and lambdaref with umpire
    umpire::util::AllocationRecord record_x{x_ref,sizeof(double)*nx,h_allocator.getAllocationStrategy()};
    resmgr.registerAllocation(x_ref,record_x);

    umpire::util::AllocationRecord record_lam{lambda_ref,sizeof(double)*(nconeq+nconineq),h_allocator.getAllocationStrategy()};
    resmgr.registerAllocation(lambda_ref,record_lam);

    // Allocate and copy xref and lambdaref to device
    double *x_ref_dev, *lambda_ref_dev;
#ifdef EXAGO_ENABLE_GPU
    umpire::Allocator d_allocator = resmgr.getAllocator("DEVICE");
    x_ref_dev = static_cast<double*>(d_allocator.allocate(nx*sizeof(double)));
    lambda_ref_dev = static_cast<double*>(d_allocator.allocate((nconeq+nconineq)*sizeof(double)));
#else
    x_ref_dev = x_ref;
    lambda_ref_dev = lambda_ref;
#endif

    resmgr.copy(x_ref_dev,x_ref);
    resmgr.copy(lambda_ref_dev,lambda_ref);
					    
    // Tests
    // fail += test.computeVariableBounds(opflowtest,xl_ref,xu_ref,resmgr);
    fail += test.computeObjective(opflowtest,x_ref_dev,obj_value);
    // fail += test.computeGradient(opflowtest,x_ref_dev,grad_ref,resmgr);
    // fail += test.computeConstraints(opflowtest,x_ref_dev,g_ref,resmgr);
    // fail += test.computeConstraintBounds(opflowtest,gl_ref,gu_ref,resmgr);
    // fail += test.computeConstraintJacobian(opflowtest,x_ref_dev,Jeq,Jineq,resmgr);

    /*
    int    nxdense = 2*opflowtest->ps->nbus;
    double *hess_dense,*hess_dense_dev;

    hess_dense     = static_cast<double*>(h_allocator.allocate(nxdense*nxdense*sizeof(double*)));
#ifdef EXAGO_ENABLE_GPU
    hess_dense_dev = static_cast<double*>(d_allocator.allocate(nxdense*nxdense*sizeof(double*)));
#else
    hess_dense_dev = hess_dense;
#endif

    fail += test.computeHessian(opflowtest,x_ref_dev,lambda_ref_dev,obj_factor,Hess,resmgr,hess_dense,hess_dense_dev);

    // Cleanup
    h_allocator.deallocate(hess_dense);
    */
#ifdef EXAGO_ENABLE_GPU
    d_allocator.deallocate(x_ref_dev);
    d_allocator.deallocate(lambda_ref_dev);
    // d_allocator.deallocate(hess_dense_dev);
#endif

    ierr = PetscFree(x_ref);CHKERRQ(ierr);
    // ierr = PetscFree(xl_ref);CHKERRQ(ierr);
    // ierr = PetscFree(xu_ref);CHKERRQ(ierr);
    // ierr = PetscFree(grad_ref);CHKERRQ(ierr);

    ierr = VecRestoreArray(X,&x_vec);CHKERRQ(ierr);
    /*
    ierr = VecRestoreArray(Xl,&xl_vec);CHKERRQ(ierr);
    ierr = VecRestoreArray(Xu,&xu_vec);CHKERRQ(ierr);
    ierr = VecRestoreArray(grad,&grad_vec);CHKERRQ(ierr);
    */
    ierr = VecRestoreArray(Lambda,&lambda_ref);CHKERRQ(ierr);

    /*
    ierr = VecRestoreArray(G,&g_ref);CHKERRQ(ierr);
    ierr = VecRestoreArray(Gl,&gl_ref);CHKERRQ(ierr);
    ierr = VecRestoreArray(Gu,&gu_ref);CHKERRQ(ierr);
    */

    ierr = VecScale(Lambda,1/obj_factor);CHKERRQ(ierr);
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
    // VecView:
    //    See the value of Lambda
    //    See the value of X
    //    Exit after
    ierr = OPFLOWGetConstraintMultipliers(opflowtest,&Lambda);CHKERRQ(ierr);

    // fail += test.computeVariableBounds(opflowtest,Xl,Xu);
    fail += test.computeObjective(opflowtest,X,obj_value);
    // fail += test.computeGradient(opflowtest,X,grad);
    // fail += test.computeConstraints(opflowtest,X,G);
    // fail += test.computeConstraintBounds(opflowtest,Gl,Gu);
    // fail += test.computeConstraintJacobian(opflowtest,X,Jeq,Jineq);
    // fail += test.computeHessian(opflowtest,X,Lambda,obj_factor,Hess);

    ierr = OPFLOWDestroy(&opflowtest);CHKERRQ(ierr);
  }

  /* Destroy OPFLOW objects */
  /*
  ierr = VecDestroy(&G);CHKERRQ(ierr);
  ierr = VecDestroy(&Gl);CHKERRQ(ierr);
  ierr = VecDestroy(&Gu);CHKERRQ(ierr);
  */
  // ierr = VecDestroy(&X);CHKERRQ(ierr);
  /*
  ierr = VecDestroy(&Xl);CHKERRQ(ierr);
  ierr = VecDestroy(&Xu);CHKERRQ(ierr);
  ierr = VecDestroy(&grad);CHKERRQ(ierr);
  */
  // ierr = VecDestroy(&Lambda);CHKERRQ(ierr);
  /*
  ierr = MatDestroy(&Jeq);CHKERRQ(ierr);
  if (ineq_present)
    ierr = MatDestroy(&Jineq);CHKERRQ(ierr);
  ierr = MatDestroy(&Hess);CHKERRQ(ierr);
  */
  ExaGOFinalize();
  // PetscFinalize();
  return fail;
}
