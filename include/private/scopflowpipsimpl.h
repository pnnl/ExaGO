/**
 * @file scopflowimpl.h
 * @brief Private header file that defines data and structures for Security constrained Optimal Power Flow application
 */
#ifndef SCOPFLOWIMPL_H
#define SCOPFLOWIMPL_H

#include <ps.h>
#include <private/psimpl.h>
#include <scopflow.h>
#if defined(SCOPFLOW_HAVE_PIPS)
#include <Drivers/parallelPipsNlp_C_Callback.h>
#endif

#include <private/contingencylist.h>

/**
  We want to solve the 2-stage optimization problem with Ns scenarios
  min \sum_{i=0}^Ns f(xi)
  s.t.
     g(xi) = 0              i = 0,Ns  (Equality constraints)
   h^- <= h(xi) <= h^+      i = 0,Ns  (Inequality constraints)
    x^- <= xi <= x+         i = 0,Ns  (Variable bounds)
 -dx <= xi - x0 <= dx       i = 1,Ns  (Coupling constraints between first stage and each scenario)

This can be decomposed for first stage and each scenario as
First stage:
  min f(x0)
  s.t.
     g(x0) = 0              (Equality constraints)
   h^- <= h(x0) <= h^+      (Inequality constraints)
    x^- <= x0 <= x+         (Variable bounds)

For each scenario i \in Ns:
  min f(xi)
  s.t.
     g(xi) = 0              (Equality constraints)
   h^- <= h(xi) <= h^+      (Inequality constraints)
    x^- <= xi <= x+         (Variable bounds)
 -dx <= xi - x0 <= dx       (Coupling constraints between first stage and each scenario)

Notes for PIPS-NLP:
 -- Each rank has the following:
   a) First stage problem
   b) One or more scenarios.
   The first stage problem is duplicated across all ranks.
 -- The first stage and the scenarios are identified by "row" and "col" variables that PIPS-NLP sets.
    0 is for the first stage, and row,col > 0 is for the scenarios.
 -- "row" is for the equations (objective, constraints, etc.) and "col" is for the variables. For example,
    with row=col=0, PIPS-NLP expects the application to set the equations for the first stage using the 
    variables for the first stage. Another example is a Jacobian which can have row \neq col. For instance
    in the Jacobian calculation with row = 1, col = 0; PIPS-NLP is expecting the application to compute the
    partial derivatives of the equations for the first scenario w.r.t. variables for the first stage, i.e. x0.
--  PIPS-NLP expects the Jacobian to be set in a compressed sparse column (CSC) format.
--  PIPS-NLP expects the Lagrangian Hessian (dL/dX) to be set in a symmetric compressed sparse column (SCSC) format
    where the values only in the lower triangle are set.
**/

 /**
  * @brief private struct for security optimal power flow application
  */

struct _p_SCOPFLOW{
  /* Sizes */
  PetscInt ns,Ns; /* Local and global number of scenarios */
  
  COMM comm; /**< Communicator context */
  OPFLOW *opflows; /* Array of optimal power flow application objects.
  		      Each processor creates ns objects, one for each 
		      scenario */

  PetscBool setupcalled; /* SCOPFLOWSetUp called? */

  PetscBool converged; /* Convergence status */

  char netfile[100]; /* Network data file */

  /* Note that with the above coupling constraints, the Jacobian
     for the coupling constraints is a constant matrix, so we 
     can use just one matrix, instead of one per scenario
  */
  Mat Jcoup;  /* Jacobian for the coupling constraints */
  Mat JcoupT; /* Transpose of the coupling Jacobian */

  PetscBool iscoupling; /* Is each scenario coupled with base scenario? */
  PetscBool first_stage_gen_cost_only; /* Only include the gen cost for first stage only */
  PetscBool ignore_line_flow_constraints; /* Ignore line flow constraints */
  PetscBool replicate_basecase;           /* Replicate base case for all scenarios */

  ContingencyList ctgclist;
  PetscBool       ctgcfileset;
  char            ctgcfile[100];

#if defined(SCOPFLOW_HAVE_PIPS)
  /* PIPS specific terms */
  PipsNlpProblemStructPtr nlp_pips; /**< PIPS solver */
#endif

};


#endif
