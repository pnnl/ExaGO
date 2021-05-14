#include <exago_config.h>

#if defined(EXAGO_ENABLE_HIOP)
#if defined(EXAGO_ENABLE_HIOP_SPARSE)

#ifndef OPFLOWHIOPSPARSE_HPP
#define OPFFLOWHIOPSPARSE_HPP
#include <opflow.h>
#include <hiopLinAlgFactory.hpp>
#include <hiopInterface.hpp>
#include <hiopMatrix.hpp>
#include <hiopNlpFormulation.hpp>
#include <hiopAlgFilterIPM.hpp>

#if defined(EXAGO_ENABLE_IPOPT)
#include <IpoptAdapter.hpp> // ipopt adapter from hiop
#include <IpIpoptApplication.hpp>
#endif

#include <cstdlib>
#include <string>
#include <cstring>
#include <cstdio>

#ifdef HIOP_USE_MAGMA
#include "magma_v2.h"
#endif

class OPFLOWHIOPSPARSEInterface : public hiop::hiopInterfaceSparse
{
public:
  OPFLOWHIOPSPARSEInterface(OPFLOW);


  ~OPFLOWHIOPSPARSEInterface() {
  }

  bool get_prob_sizes(long long& n, long long& m);

  bool get_vars_info(const long long& n, double *xlow, double* xupp, NonlinearityType* type);

  bool get_cons_info(const long long& m, double* clow, double* cupp, NonlinearityType* type);

  bool get_sparse_blocks_info(int& nx,int& nnz_sparse_Jaceq, int& nnz_sparse_Jacineq,int& nnz_sparse_Hess_Lagr);

  bool eval_f(const long long& n, const double* x, bool new_x, double& obj_value);

  bool eval_cons(const long long& n, const long long& m, 
		 const long long& num_cons, const long long* idx_cons,  
		 const double* x, bool new_x, double* cons);

  bool eval_cons(const long long& n, const long long& m,
                         const double* x, bool new_x,
                         double* cons);

  bool eval_grad_f(const long long& n, const double* x, bool new_x, double* gradf);

  virtual bool eval_Jac_cons(const long long& n, const long long& m,
			     const long long& num_cons, const long long* idx_cons,
			     const double* x, bool new_x,
			     const int& nnzJacS, int* iJacS, int* jJacS, double* MJacS);
  bool eval_Jac_cons(const long long& n, const long long& m,
			     const double* x, bool new_x,
			     const int& nnzJacS, int* iJacS, int* jJacS, double* MJacS);
  bool get_starting_point(const long long&n, double* x0);
  bool eval_Hess_Lagr(const long long& n, const long long& m,
			      const double* x, bool new_x, const double& obj_factor,
			      const double* lambda, bool new_lambda,
			      const int& nnzHSS, int* iHSS, int* jHSS, double* MHSS);

  void solution_callback(hiop::hiopSolveStatus status,
			 int n, const double* x,
			 const double* z_L,
			 const double* z_U,
			 int m, const double* g,
			 const double* lambda,
			 double obj_value);

  bool get_MPI_comm(MPI_Comm& comm_out) {comm_out = MPI_COMM_SELF; return true;}

private:
  OPFLOW   opflow;
};

typedef struct _p_OPFLOWSolver_HIOPSPARSE *OPFLOWSolver_HIOPSPARSE;

struct _p_OPFLOWSolver_HIOPSPARSE {
  
  OPFLOWHIOPSPARSEInterface     *nlp;
  hiop::hiopSolveStatus         status;
  hiop::hiopNlpSparse           *sp;
  hiop::hiopAlgFilterIPMNewton  *solver;

#if defined(EXAGO_ENABLE_IPOPT)
  // Ipopt Adapter structs
  SmartPtr<hiop::hiopSparse2IpoptTNLP>    ipoptTNLP;
  SmartPtr<Ipopt::IpoptApplication>       ipoptApp;
  PetscBool                               ipopt_debug;
#endif

};

#endif
#endif
#endif
