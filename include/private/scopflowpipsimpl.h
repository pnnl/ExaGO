/**
 * @file scopflowimpl.h
 * @brief Private header file that defines data and structures for Security constrained Optimal Power Flow application
 */
#ifndef SCOPFLOWIMPL_H
#define SCOPFLOWIMPL_H

#include <ps.h>
#include <private/psimpl.h>
#include <scopflow.h>
#if defined(PSAPPS_HAVE_PIPS)
#include <Drivers/parallelPipsNlp_C_Callback.h>
#endif

/**
  We want to solve the 2-stage optimization problem with Ns scenarios
  min \sum_{i=0}^Ns f(i)
  s.t.
     g(xi) >= 0           i = 0,Ns
  x^- <= xi <= x+         i = 0,Ns
  -|xi - x0| >= -\deltax  i = 1,Ns 

We can combine the first and third constraint to have an inequality constraint vector,
   c(x) = [g(xi) ; d(xi,x0) ]
where d = -|xi - x0| + \deltax

to have
min \sum_{i=0}^Ns f(i)                                                                                                                     
  s.t.                                                                                                                
     c(xi) >= 0      i = 0,Ns
  x^- <= xi <= x+    i = 0,Ns

-- Let each processor handle ns scenarios with \sum ns = Ns. Note that
   ns can be different on different processors.
-- Each processor handles local ns scenario(s) (\sum ns = Ns)

Obj. function routine(input Vec X, output scalar obj):
 -- Global X to Local x (VecPlaceArray)
 -- Loop k=1:ns
    -- Get vector xk for scenario k from vector x (VecPlaceArray)
    -- Compute obj. function f(xk) for scenario k. (OPFLOW)
 -- Sum objective function for ns scenarios f_rank = \sum_{k=1}^ns f(k)
 -- Reduction sum operation to get the total objective function value \sum_{rank=0}^size f_rank

Constraints function routine(input Vec X, output Vec C):.
 -- Scatter x0 from Global X to local x0 on all processors (Does PIPS do this?)
 -- Global X to local x (VecPlaceArray)
 
 -- Global C to local c (VecPlaceArray)
 -- Loop k=1:ns
    -- Get vector xk for scenario k from vector x (VecPlaceArray)
    -- Get vector ck for scenario k from vector c (VecPlaceArray)
    -- Get vector gk for scenario k from vector ck (VecPlaceArray)
    -- Compute g(xk) for scenario k (OPFLOW)
    -- Get vector dk for scenario k from vector k (Only rank > 0)
    -- Compute d(xk,x0) for scenario k (Only rank > 0)


**/

 /**
  * @brief private struct for optimal power flow application
  */
struct _p_SCOPFLOW{
  /* Sizes */
  PetscInt ns,Ns; /* Local and global number of scenarios */
  PetscInt nx,Nx; /* Local and global size of the decision variable vector */
  PetscInt ng,Ng; /* Local and global size of equality and inequality constraints */
  PetscInt nd,Nd; /* Local and global size of coupling constraints */
  PetscInt nc,Nc; /* Local and global size of equality,inequality + coupling constraints */
  
  COMM comm; /**< Communicator context */
  OPFLOW *opflow; /* Array of optimal power flow application objects.
		   Each processor creates ns objects, one for each 
		  scenario */
   /* Solution vector for the entire problem. 
     X = [X0;X1;...;Xn]*/  
  Vec  X;    /* Parallel vector X \in R^{Ns*nx} to hold the solution for the entire problem */
  Vec  Xloc; /* Serial local X \in R^{ns*nx} vector to hold the solution of ns scenarios */
  Vec  X0;   /* Serial vector X0 \in R^{nx} to hold the solution of the first stage problem. This will be 
                non-empty vector except rank 0. The values from X that correspond to the stage 1 solution, 
                are scattered to X0 */
  Vec  Xi;   /* Serial vector Xs \in R^{nx} to hold the solution for the ith scenario */
  /*Inequality and equality constraint function */
  Vec  C;    /* Parallel vector \in R^{Ns*ng + (Ns-1)*nd} to hold the inequality constraint vector */
  Vec  Cloc; /* Serial vector in R^{ns*(ng + nd)} to hold the local inequality constraint vector for ns scenarios. Note that
		size of Cloc on rank 0 (which has the base case) is R^{ns*ng + (ns-1)*nd} */
  Vec  Ci;  /* Serial vector in R^{ng+nd} to hold the local inequality constraints for ith scenario */
  Vec  Di;  /* Serial vector in R^{nd} to hold the d type inequality constraints for the ith scenario */

  Vec Cl; /* Parallel vector in R^{Ns*ng + (Ns-1)*nd} for lower bounds */
  Vec Cu; /* Parallel vector in R^{Ns*ng + (Ns-1)*nd} for upper bounds */
           

  Vec Xl; /* Parallel vector in R^{Ns*nx} for lower bounds */
  Vec Xu; /* Parallel vector in R^{Ns*nx} for upper bounds */

  Vec         F;   /* Parallel vector \in R^{Ns} to hold the objective functions for each processor.
		      Each processor contributes ns elements, i.e., the objective functions for
		      its scenarios. */
  PetscScalar obj; /* Objective function value \sum_{i=0}^Ns f(i) */

  Vec gradobj; /**< Gradient of the objective function */

  PetscBool setupcalled; /* SCOPFLOWSetUp called? */


  PetscBool converged; // Convergence status

  /* For PIPS */
  Mat  Jac;  /* Jacobian of constraints */
  Mat  Hes;  /* Lagrangian Hessian */

  PetscScalar *x; /**< Solution array - same as the array for X */

  PetscInt nnz_jac_g; /**< Number of nonzeros in the jacobian of the constraints */
  PetscInt nnz_hes; /**< Number of nonzeros in the Lagrangian Hessian */

  /* Lagrange multipliers */
  Vec lambda_g;
  Vec lambda_xl;
  Vec lambda_xu;

#if defined(PSAPPS_HAVE_PIPS)
  /* PIPS specific terms */
  PipsNlpProblemStructPtr nlp_pips; /**< PIPS solver */
#endif

};


#endif
