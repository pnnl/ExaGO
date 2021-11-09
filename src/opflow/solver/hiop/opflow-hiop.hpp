#include <exago_config.h>

#if defined(EXAGO_ENABLE_HIOP)

#ifndef OPFLOWHIOP_HPP
#define OPFLOWHIOP_HPP
#include <opflow.h>
#include <hiopInterface.hpp>
#include <hiopMatrix.hpp>
#include <hiopNlpFormulation.hpp>
#include <hiopAlgFilterIPM.hpp>

#if defined(EXAGO_ENABLE_IPOPT)
#include <IpoptAdapter.hpp>
#include <IpIpoptApplication.hpp>
#endif

#include <cstdlib>
#include <string>
#include <cstring>
#include <cstdio>

#ifdef HIOP_USE_MAGMA
#include "magma_v2.h"
#endif

class OPFLOWHIOPInterface : public hiop::hiopInterfaceMDS
{
public:
  OPFLOWHIOPInterface(OPFLOW);


  ~OPFLOWHIOPInterface() {
  }

  bool get_prob_sizes(hiop::size_type& n, hiop::size_type& m);

  bool get_vars_info(const hiop::size_type& n, double *xlow, double* xupp, NonlinearityType* type);

  bool get_cons_info(const hiop::size_type& m, double* clow, double* cupp, NonlinearityType* type);

  bool get_sparse_dense_blocks_info(int& nx_sparse, int& nx_dense,
				    int& nnz_sparse_Jace, int& nnz_sparse_Jaci,
				    int& nnz_sparse_Hess_Lagr_SS, int& nnz_sparse_Hess_Lagr_SD);

  bool eval_f(const hiop::size_type& n, const double* x, bool new_x, double& obj_value);

  bool eval_cons(const hiop::size_type& n, const hiop::size_type& m, 
		 const hiop::size_type& num_cons, const hiop::size_type* idx_cons,  
		 const double* x, bool new_x, double* cons);

  bool eval_grad_f(const hiop::size_type& n, const double* x, bool new_x, double* gradf);

  bool eval_Jac_cons(const hiop::size_type& n, const hiop::size_type& m, 
		     const hiop::size_type& num_cons, const hiop::size_type* idx_cons,
		     const double* x, bool new_x,
		     const hiop::size_type& nsparse, const hiop::size_type& ndense, 
		     const int& nnzJacS, int* iJacS, int* jJacS, double* MJacS, 
		     double* JacD);

  bool eval_Hess_Lagr(const hiop::size_type& n, const hiop::size_type& m, 
		      const double* x, bool new_x, const double& obj_factor,
		      const double* lambda, bool new_lambda,
		      const hiop::size_type& nsparse, const hiop::size_type& ndense, 
		      const int& nnzHSS, int* iHSS, int* jHSS, double* MHSS, 
		      double* HDD,
		      int& nnzHSD, int* iHSD, int* jHSD, double* MHSD);

  bool get_starting_point(const hiop::size_type& global_n, double* x0);

  void solution_callback(hiop::hiopSolveStatus status,
                                 int n, const double* x,
                                 const double* z_L,
                                 const double* z_U,
                                 int m, const double* g,
                                 const double* lambda,
			         double obj_value);

  bool iterate_callback(int iter,
         double obj_value,
         double logbar_obj_value,
			int n,
         const double* x,
			const double* z_L,
			const double* z_U,
         int m_ineq,
         const double *s,
			int m,
         const double* g,
			const double* lambda,
			double inf_pr,
         double inf_du,
         double onenorm_pr_,
			double mu,
			double alpha_du,
         double alpha_pr,
			int ls_trials);

  bool get_MPI_comm(MPI_Comm& comm_out) {comm_out = MPI_COMM_SELF; return true;}

  void naturaltospdense(const double*,double*);

  void spdensetonatural(const double*,double*);
private:
  OPFLOW   opflow;
};

typedef struct _p_OPFLOWSolver_HIOP *OPFLOWSolver_HIOP;

struct _p_OPFLOWSolver_HIOP {
  
  OPFLOWHIOPInterface           *nlp;
  hiop::hiopSolveStatus         status;
  hiop::hiopNlpMDS              *mds;
  hiop::hiopAlgFilterIPMNewton  *solver;

#if defined(EXAGO_ENABLE_IPOPT)
  // Ipopt Adapter structs
  SmartPtr<hiop::hiopMDS2IpoptTNLP>       ipoptTNLP;
  SmartPtr<Ipopt::IpoptApplication>       ipoptApp;
  PetscBool                               ipopt_debug;
#endif
};

#endif
#endif
