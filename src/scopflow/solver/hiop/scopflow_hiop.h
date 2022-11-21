#include <exago_config.h>

#if defined(EXAGO_ENABLE_HIOP)
/*
   Primal decomposition based solver from HIOP
*/
#ifndef SCOPFLOWHIOP_HPP
#define SCOPFLOWHIOP_HPP

#include <scopflow.h>

#include <hiopAlgPrimalDecomp.hpp>
#include <hiopInterfacePrimalDecomp.hpp>

class SCOPFLOWHIOPInterface : public hiop::hiopInterfacePriDecProblem {
public:
  SCOPFLOWHIOPInterface(SCOPFLOW);

  ~SCOPFLOWHIOPInterface();

  hiop::hiopSolveStatus solve_master(hiop::hiopVector &x, const bool &include_r,
                                     const double &rval, const double *grad,
                                     const double *hess,
                                     const char *master_options_file);

  bool eval_f_rterm(size_t idx, const int &n, const double *x, double &rval);
  bool eval_grad_rterm(size_t idx, const int &n, double *x,
                       hiop::hiopVector &grad);

  size_t get_num_rterms() const;

  size_t get_num_vars() const;

  void get_solution(double *x) const;

  double get_objective();

  bool set_recourse_approx_evaluator(
      const int n,
      hiopInterfacePriDecProblem::RecourseApproxEvaluator *evaluator);

  int nxcoup;     /* Number of coupling variables */
  int *loc_xcoup; /* Indices for the coupling variables in the base problem */
  bool include_r_;
  hiopInterfacePriDecProblem::RecourseApproxEvaluator *rec_evaluator;

private:
  SCOPFLOW scopflow;
  OPFLOW opflowctgc;
};

typedef struct _p_SCOPFLOWSolver_HIOP *SCOPFLOWSolver_HIOP;

struct _p_SCOPFLOWSolver_HIOP {
  SCOPFLOWHIOPInterface *pridecompprob;
  hiop::hiopAlgPrimalDecomposition *pridecsolver;
  hiop::hiopSolveStatus status;
};

#endif
#endif
