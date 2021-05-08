#include <exago_config.h>

#if defined(EXAGO_ENABLE_HIOP)
#if defined(EXAGO_ENABLE_HIOP_DISTRIBUTED)
/* 
   Primal decomposition based solver from HIOP
*/
#ifndef SOPFLOWHIOP_HPP
#define SOPFLOWHIOP_HPP

#include <cstdlib>
#include <string>
#include <cstring>
#include <cstdio>

#include <hiopInterfacePrimalDecomp.hpp>
#include <hiopAlgPrimalDecomp.hpp>

class SOPFLOWHIOPInterface : public hiop::hiopInterfacePriDecProblem
{
public:
  SOPFLOWHIOPInterface(SOPFLOW);

  ~SOPFLOWHIOPInterface();

  hiop::hiopSolveStatus solve_master(double* x,
                                     const bool& include_r,
                                     const double& rval, 
                                     const double* grad,
                                     const double*hess);


  bool eval_f_rterm(size_t idx, const int& n, const double* x, double& rval);
  bool eval_grad_rterm(size_t idx, const int& n, double* x, double* grad);

  size_t get_num_rterms() const;

  size_t get_num_vars() const;

  void get_solution(double* x) const;
  double get_objective();

  bool set_recourse_approx_evaluator(const int n,hiopInterfacePriDecProblem::RecourseApproxEvaluator* evaluator);

  int      nxcoup;    /* Number of coupling variables */
  std::vector<int> loc_xcoup; /* Indices for the coupling variables in the base problem */
  bool     include_r_;
  hiopInterfacePriDecProblem::RecourseApproxEvaluator* rec_evaluator;
private:
  SOPFLOW sopflow;
  OPFLOW   opflowscen;

};

typedef struct _p_SOPFLOWSolver_HIOP *SOPFLOWSolver_HIOP;

struct _p_SOPFLOWSolver_HIOP {
  SOPFLOWHIOPInterface  *pridecompprob;
  hiop::hiopAlgPrimalDecomposition *pridecsolver;
  hiop::hiopSolveStatus status;
};

#endif
#endif
#endif
