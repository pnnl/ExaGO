#include <exago_config.h>

#if defined(EXAGO_ENABLE_RAJA)

#pragma once

#include <umpire/Allocator.hpp>
#include <umpire/ResourceManager.hpp>

#include "pbpolrajahiop.h"
#include <opflow.h>

/**

  @brief Takes a raw pointer to some data and registers it with the
  allocator and resource manager.

  TODO: Ensure that we do not have to deregister the allocation with
  the resource manager. This is done with `resmgr.deregisterAllocation(T* ptr)`,
  but must we do this with ever allocation? Will we incurr a memory leak if not?

*/
template <typename T, typename SizeType>
void registerWith(T *ptr, SizeType N, umpire::ResourceManager &resmgr,
                  umpire::Allocator &allocator) {
  umpire::util::AllocationRecord record{ptr, sizeof(T) * N,
                                        allocator.getAllocationStrategy()};
  resmgr.registerAllocation(ptr, record);
}

template <typename T, typename I> T *paramAlloc(umpire::Allocator &a, I N) {
  return static_cast<T *>(a.allocate(N * sizeof(T)));
}

/// Class storing bus parameters
struct BUSParamsRajaHiop {
  // Host data
  int nbus;                       /* Number of buses */
  int *isref;                     /* isref[i] = 1 if bus is reference bus */
  int *isisolated;                /* isisolated[i] = 1 if bus is isolated bus */
  int *ispvpq;                    /* For all other buses */
  double *powerimbalance_penalty; /* Penalty for power imbalance */
  double *vmin;                   /* min. voltage magnitude limit */
  double *vmax;                   /* max. voltage magnitude limit */
  double *va;      /* bus angle (from file only used in bounds) */
  double *vm;      /* bus voltage magnitude (from file only used in bounds) */
  double *gl;      /* bus shunt (conductance) */
  double *bl;      /* bus shunt (suspectance) */
  int *xidx;       /* starting locations for bus variables in X vector */
  int *xidxpimb;   /* starting locations for power imbalance bus variables in X
                      vector */
  int *gidx;       /* starting locations for bus balance equations in constraint
                      vector */
  int *jacsp_idx;  /* Location number in the sparse Jacobian for Pimb */
  int *jacsq_idx;  /* Location number in the sparse Jacobian for Qimb */
  int *hesssp_idx; /* Location number in the Hessian */

  // Device data
  int *isref_dev_;      /* isref[i] = 1 if bus is reference bus */
  int *isisolated_dev_; /* isisolated[i] = 1 if bus is isolated bus */
  int *ispvpq_dev_;     /* For all other buses */
  double *powerimbalance_penalty_dev_; /* Penalty for power imbalance */
  double *vmin_dev_;                   /* min. voltage magnitude limit */
  double *vmax_dev_;                   /* max. voltage magnitude limit */
  double *va_dev_; /* bus angle (from file only used in bounds) */
  double *vm_dev_; /* bus voltage magnitude (from file only used in bounds) */
  double *gl_dev_; /* bus shunt (conductance) */
  double *bl_dev_; /* bus shunt (suspectance) */
  int *xidx_dev_;  /* starting locations for bus variables in X vector */
  int *xidxpimb_dev_; /* starting locations for power imbalance bus variables in
                         X vector */
  int *gidx_dev_; /* starting locations for bus balance equations in constraint
                     vector */
  int *jacsp_idx_dev_;  /* Location number in the sparse Jacobian for Pimb */
  int *jacsq_idx_dev_;  /* Location number in the sparse Jacobian for Qimb */
  int *hesssp_idx_dev_; /* Location number in the Hessian */

  int allocate(OPFLOW);
  int destroy(OPFLOW);
  int copy(OPFLOW);

private:
  // Umpire memory allocators
  umpire::Allocator h_allocator_;
  umpire::Allocator d_allocator_;
};

struct GENParamsRajaHiop {
public:
  // Host data
  int ngenON;         /* Number of generators with STATUS ON */
  double *cost_alpha; /* generator cost coefficients */
  double *cost_beta;  /* generator cost coefficients */
  double *cost_gamma; /* generator cost coefficients */
  double *pt;         /* min. active power gen. limits */
  double *pb;         /* max. active power gen. limits */
  double *qt;         /* min. reactive power gen. limits */
  double *qb;         /* max. reactive power gen. limits */
  double *pgs;        /* real power output setpoint */
  int *isrenewable;   /* Is renewable generator? */

  int *xidx; /* starting locations in X vector */
  int *
      gidxbus; /* starting locations in constraint vector for bus constraints */
  int *geqidxgen;    /* starting locations in equality constraint vector for gen
                        constraints */
  int *gineqidxgen;  /* starting locations in inequality constraint vector for
                        gen constraints */
  int *gbineqidxgen; /* Starting location to insert contribution to inequality
                        constraint bound */

  /* The following members are only used with HIOP */
  int *eqjacspbus_idx; /* Location number in the bus equality constraints sparse
                          Jacobian for Pg */
  int *eqjacsqbus_idx; /* Location number in the bus equality constraints sparse
                          Jacobian for Qg */
  int *eqjacspgen_idx; /* Location number in the gen equality constraints sparse
                          Jacobian for Pg */
  int *ineqjacspgen_idx; /* Location number in the bus equality constraints
                            sparse Jacobian for Pg */
  int *hesssp_idx;       /* Location number in the Hessian */

  // Device data
  double *cost_alpha_dev_; /* generator cost coefficients */
  double *cost_beta_dev_;  /* generator cost coefficients */
  double *cost_gamma_dev_; /* generator cost coefficients */
  double *pt_dev_;         /* min. active power gen. limits */
  double *pb_dev_;         /* max. active power gen. limits */
  double *qt_dev_;         /* min. reactive power gen. limits */
  double *qb_dev_;         /* max. reactive power gen. limits */
  double *pgs_dev_;        /* real power output setpoint */
  int *isrenewable_dev_;   /* Is renewable generator? */

  int *xidx_dev_;        /* starting locations in X vector */
  int *gidxbus_dev_;     /* starting locations in constraint vector for bus
                            constraints */
  int *geqidxgen_dev_;   /* starting locations in equality constraint vector for
                            gen constraints */
  int *gineqidxgen_dev_; /* starting locations in inequality constraint vector
                            for gen constraints */
  int *gbineqidxgen_dev_; /* Starting location to insert contribution to
                             inequality constraint bound */

  /* The following members are only used with HIOP */
  int *eqjacspbus_idx_dev_;   /* Location number in the bus equality constraints
                                 sparse Jacobian for Pg */
  int *eqjacsqbus_idx_dev_;   /* Location number in the bus equality constraints
                                 sparse Jacobian for Qg */
  int *eqjacspgen_idx_dev_;   /* Location number in the gen equality constraints
                                 sparse Jacobian for Pg */
  int *ineqjacspgen_idx_dev_; /* Location number in the bus equality constraints
                                 sparse Jacobian for Pg */
  int *hesssp_idx_dev_;       /* Location number in the Hessian */

  int allocate(OPFLOW);
  int destroy(OPFLOW);
  int copy(OPFLOW);

private:
  // Umpire memory allocators
  umpire::Allocator h_allocator_;
  umpire::Allocator d_allocator_;
};

struct LOADParamsRajaHiop {
  // Host data
  int nload;                /* Number of loads */
  double *pl;               /* active power demand */
  double *ql;               /* reactive power demand */
  double *loadloss_penalty; /* Penalty for load loss */
  int *xidx;                /* starting location in X vector */
  int *gidx;                /* starting location in constraint vector */

  /* The following members are only used with HIOP */
  int *jacsp_idx;  /* Location number in the sparse Jacobian for delPload */
  int *jacsq_idx;  /* Location number in the sparse Jacobian for delQload */
  int *hesssp_idx; /* Location number in the Hessian */

  double *pl_dev_;               /* active power demand */
  double *ql_dev_;               /* reactive power demand */
  double *loadloss_penalty_dev_; /* Penalty for load loss */
  int *xidx_dev_;                /* starting location in X vector */
  int *gidx_dev_;                /* starting location in constraint vector */

  int *jacsp_idx_dev_; /* Location number in the sparse Jacobian for delPload */
  int *jacsq_idx_dev_; /* Location number in the sparse Jacobian for delQload */
  int *hesssp_idx_dev_; /* Location number in the Hessian */

  int allocate(OPFLOW);
  int destroy(OPFLOW);
  int copy(OPFLOW);

private:
  // Umpire memory allocators
  umpire::Allocator h_allocator_;
  umpire::Allocator d_allocator_;
};

struct LINEParamsRajaHiop {
  // Host data
  int nlineON;   /* Number of active lines (STATUS = 1) */
  int nlinelim;  /* Active lines + limits */
  double *Gff;   /* From side self conductance */
  double *Bff;   /* From side self susceptance */
  double *Gft;   /* From-to transfer conductance */
  double *Bft;   /* From-to transfer susceptance */
  double *Gtf;   /* To-from transfer conductance */
  double *Btf;   /* To-from transfer susceptance */
  double *Gtt;   /* To side self conductance */
  double *Btt;   /* To side self susceptance */
  double *rateA; /* Line MVA rating A (normal operation) */
  int *xidxf;    /* Starting locatin of from bus voltage variables */
  int *xidxt;    /* Starting location of to bus voltage variables */
  int *geqidxf;  /* Starting location of from side to insert equality constraint
                    contribution in constraints vector */
  int *geqidxt;  /* Starting location of to side to insert equality constraint
                    contribution in constraints vector */
  int *gineqidx; /* Starting location to insert contribution to inequality
                    constraint */
  int *gbineqidx;  /* Starting location to insert contribution to inequality
                      constraint bound */
  int *linelimidx; /* Indices for subset of lines that have finite limits */

  // Device data
  double *Gff_dev_;    /* From side self conductance */
  double *Bff_dev_;    /* From side self susceptance */
  double *Gft_dev_;    /* From-to transfer conductance */
  double *Bft_dev_;    /* From-to transfer susceptance */
  double *Gtf_dev_;    /* To-from transfer conductance */
  double *Btf_dev_;    /* To-from transfer susceptance */
  double *Gtt_dev_;    /* To side self conductance */
  double *Btt_dev_;    /* To side self susceptance */
  double *rateA_dev_;  /* Line MVA rating A (normal operation) */
  int *xidxf_dev_;     /* Starting locatin of from bus voltage variables */
  int *xidxt_dev_;     /* Starting location of to bus voltage variables */
  int *geqidxf_dev_;   /* Starting location of from side to insert equality
                          constraint contribution in constraints vector */
  int *geqidxt_dev_;   /* Starting location of to side to insert equality
                          constraint contribution in constraints vector */
  int *gineqidx_dev_;  /* Starting location to insert contribution to inequality
                          constraint */
  int *gbineqidx_dev_; /* Starting location to insert contribution to inequality
                          constraint bound */
  int *
      linelimidx_dev_; /* Indices for subset of lines that have finite limits */

  int allocate(OPFLOW);
  int destroy(OPFLOW);
  int copy(OPFLOW);

private:
  // Umpire memory allocators
  umpire::Allocator h_allocator_;
  umpire::Allocator d_allocator_;
};

struct PbpolModelRajaHiop : public _p_FormPBPOLRAJAHIOP {
  PbpolModelRajaHiop(OPFLOW opflow) {}

  void destroy(OPFLOW opflow) {
    loadparams.destroy(opflow);
    lineparams.destroy(opflow);
    busparams.destroy(opflow);
    genparams.destroy(opflow);
  }

  ~PbpolModelRajaHiop() {}

  GENParamsRajaHiop genparams;
  LOADParamsRajaHiop loadparams;
  LINEParamsRajaHiop lineparams;
  BUSParamsRajaHiop busparams;
};

#endif // EXAGO_ENABLE_HIOP
