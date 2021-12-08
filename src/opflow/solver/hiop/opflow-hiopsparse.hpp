#include <exago_config.h>

#if defined(EXAGO_ENABLE_HIOP)
#if defined(EXAGO_ENABLE_HIOP_SPARSE)

#ifndef OPFLOWHIOPSPARSE_HPP
#define OPFFLOWHIOPSPARSE_HPP
#include <hiopAlgFilterIPM.hpp>
#include <hiopInterface.hpp>
#include <hiopLinAlgFactory.hpp>
#include <hiopMatrix.hpp>
#include <hiopNlpFormulation.hpp>
#include <opflow.h>

#if defined(EXAGO_ENABLE_IPOPT)
#include <IpIpoptApplication.hpp>
#include <IpoptAdapter.hpp> // ipopt adapter from hiop
#endif

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

#ifdef HIOP_USE_MAGMA
#include "magma_v2.h"
#endif

class OPFLOWHIOPSPARSEInterface : public hiop::hiopInterfaceSparse {
public:
  OPFLOWHIOPSPARSEInterface(OPFLOW);

  ~OPFLOWHIOPSPARSEInterface() {}

  bool get_prob_sizes(hiop::size_type &n, hiop::size_type &m);

  bool get_vars_info(const hiop::size_type &n, double *xlow, double *xupp,
                     NonlinearityType *type);

  bool get_cons_info(const hiop::size_type &m, double *clow, double *cupp,
                     NonlinearityType *type);

  bool get_sparse_blocks_info(hiop::size_type &nx,
                              hiop::size_type &nnz_sparse_Jaceq,
                              hiop::size_type &nnz_sparse_Jacineq,
                              hiop::size_type &nnz_sparse_Hess_Lagr);

  bool eval_f(const hiop::size_type &n, const double *x, bool new_x,
              double &obj_value);

  bool eval_cons(const hiop::size_type &n, const hiop::size_type &m,
                 const hiop::size_type &num_cons,
                 const hiop::size_type *idx_cons, const double *x, bool new_x,
                 double *cons);

  bool eval_cons(const hiop::size_type &n, const hiop::size_type &m,
                 const double *x, bool new_x, double *cons);

  bool eval_grad_f(const hiop::size_type &n, const double *x, bool new_x,
                   double *gradf);

  bool eval_Jac_cons(const hiop::size_type &n, const hiop::size_type &m,
                     const hiop::size_type &num_cons,
                     const hiop::index_type *idx_cons, const double *x,
                     bool new_x, const hiop::size_type &nnzJacS,
                     hiop::index_type *iJacS, hiop::index_type *jJacS,
                     double *MJacS);

  bool eval_Jac_cons(const hiop::size_type &n, const hiop::size_type &m,
                     const double *x, bool new_x,
                     const hiop::size_type &nnzJacS, hiop::index_type *iJacS,
                     hiop::index_type *jJacS, double *MJacS);

  bool get_starting_point(const hiop::size_type &n, double *x0);
  bool eval_Hess_Lagr(const hiop::size_type &n, const hiop::size_type &m,
                      const double *x, bool new_x, const double &obj_factor,
                      const double *lambda, bool new_lambda,
                      const hiop::size_type &nnzHSS, hiop::index_type *iHSS,
                      hiop::index_type *jHSS, double *MHSS);

  void solution_callback(hiop::hiopSolveStatus status, int n, const double *x,
                         const double *z_L, const double *z_U, int m,
                         const double *g, const double *lambda,
                         double obj_value);

  bool iterate_callback(int iter, double obj_value, double logbar_obj_value,
                        int n, const double *x, const double *z_L,
                        const double *z_U, int m_ineq, const double *s, int m,
                        const double *g, const double *lambda, double inf_pr,
                        double inf_du, double onenorm_pr_, double mu,
                        double alpha_du, double alpha_pr, int ls_trials);

  bool get_MPI_comm(MPI_Comm &comm_out) {
    comm_out = MPI_COMM_SELF;
    return true;
  }

private:
  OPFLOW opflow;
};

typedef struct _p_OPFLOWSolver_HIOPSPARSE *OPFLOWSolver_HIOPSPARSE;

struct _p_OPFLOWSolver_HIOPSPARSE {

  OPFLOWHIOPSPARSEInterface *nlp;
  hiop::hiopSolveStatus status;
  hiop::hiopNlpSparse *sp;
  hiop::hiopAlgFilterIPMNewton *solver;

#if defined(EXAGO_ENABLE_IPOPT)
  // Ipopt Adapter structs
  SmartPtr<hiop::hiopSparse2IpoptTNLP> ipoptTNLP;
  SmartPtr<Ipopt::IpoptApplication> ipoptApp;
  PetscBool ipopt_debug;
#endif
};

#endif
#endif
#endif
