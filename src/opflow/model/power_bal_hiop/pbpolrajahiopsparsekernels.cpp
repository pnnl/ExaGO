
#include <iostream>
#include <iomanip>
#include <vector>
#include <tuple>
#include <algorithm>

#include <exago_config.h>

#if defined(EXAGO_ENABLE_RAJA)
#if defined(EXAGO_ENABLE_HIOP_SPARSE)

#include <umpire/Allocator.hpp>
#include <umpire/ResourceManager.hpp>

#include <RAJA/RAJA.hpp>
#include <private/raja_exec_config.h>

#include <private/psimpl.h>
#include <private/opflowimpl.h>
#include "pbpolrajahiopsparsekernels.hpp"
#include "pbpolrajahiopsparse.hpp"

static const bool debugmsg(true);
static const bool oldhostway(false);

PetscErrorCode OPFLOWSetInitialGuessArray_PBPOLRAJAHIOPSPARSE(OPFLOW opflow,
                                                              double *x0_dev) {
  PetscErrorCode ierr;
  double *x;
  auto &resmgr = umpire::ResourceManager::getInstance();

  PetscFunctionBegin;
  //  ierr = PetscPrintf(MPI_COMM_SELF,"Entered
  //  OPFLOWInitialization\n");CHKERRQ(ierr);
  // Do initialization on Host
  ierr = (*opflow->modelops.setinitialguess)(opflow, opflow->X, opflow->Lambda);
  CHKERRQ(ierr);
  ierr = VecGetArray(opflow->X, &x);
  CHKERRQ(ierr);

  // Copy from host to device
  umpire::Allocator h_allocator_ = resmgr.getAllocator("HOST");
  registerWith(x, opflow->nx, resmgr, h_allocator_);
  resmgr.copy(x0_dev, x);

  ierr = VecRestoreArray(opflow->X, &x);
  CHKERRQ(ierr);

  //  ierr = PetscPrintf(MPI_COMM_SELF,"Exit Initialization\n");CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode
OPFLOWSetVariableBoundsArray_PBPOLRAJAHIOPSPARSE(OPFLOW opflow, double *xl_dev,
                                                 double *xu_dev) {
  PetscErrorCode ierr;
  double *xl, *xu;
  auto &resmgr = umpire::ResourceManager::getInstance();

  PetscFunctionBegin;

  // Compute variable bounds on the host
  ierr = (*opflow->modelops.setvariablebounds)(opflow, opflow->Xl, opflow->Xu);
  CHKERRQ(ierr);
  ierr = VecGetArray(opflow->Xl, &xl);
  CHKERRQ(ierr);
  ierr = VecGetArray(opflow->Xu, &xu);
  CHKERRQ(ierr);

  // Copy host to device
  umpire::Allocator h_allocator_ = resmgr.getAllocator("HOST");
  registerWith(xl, opflow->nx, resmgr, h_allocator_);
  registerWith(xu, opflow->nx, resmgr, h_allocator_);
  resmgr.copy(xl_dev, xl);
  resmgr.copy(xu_dev, xu);

  ierr = VecRestoreArray(opflow->Xl, &xl);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(opflow->Xu, &xu);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSetConstraintBoundsArray_PBPOLRAJAHIOPSPARSE(
    OPFLOW opflow, double *gl_dev, double *gu_dev) {

  PetscErrorCode ierr;
  double *gl, *gu;
  auto &resmgr = umpire::ResourceManager::getInstance();

  PetscFunctionBegin;

  // Calculate constraint bounds on host
  ierr =
      (*opflow->modelops.setconstraintbounds)(opflow, opflow->Gl, opflow->Gu);
  CHKERRQ(ierr);

  ierr = VecGetArray(opflow->Gl, &gl);
  CHKERRQ(ierr);
  ierr = VecGetArray(opflow->Gu, &gu);
  CHKERRQ(ierr);

  umpire::Allocator h_allocator_ = resmgr.getAllocator("HOST");
  registerWith(gl, opflow->ncon, resmgr, h_allocator_);
  registerWith(gu, opflow->ncon, resmgr, h_allocator_);

  // Copy from host to device
  resmgr.copy(gl_dev, gl);
  resmgr.copy(gu_dev, gu);

  ierr = VecRestoreArray(opflow->Gl, &gl);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(opflow->Gu, &gu);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/** EQUALITY CONSTRAINTS */
PetscErrorCode OPFLOWComputeEqualityConstraintsArray_PBPOLRAJAHIOPSPARSE(
    OPFLOW opflow, const double *x_dev, double *ge_dev) {
  PbpolModelRajaHiop *pbpolrajahiopsparse =
      reinterpret_cast<PbpolModelRajaHiop *>(opflow->model);
  BUSParamsRajaHiop *busparams = &pbpolrajahiopsparse->busparams;
  GENParamsRajaHiop *genparams = &pbpolrajahiopsparse->genparams;
  LOADParamsRajaHiop *loadparams = &pbpolrajahiopsparse->loadparams;
  LINEParamsRajaHiop *lineparams = &pbpolrajahiopsparse->lineparams;
  PetscErrorCode ierr;
  PetscInt flps = 0;
  int include_loadloss_variables = opflow->include_loadloss_variables;
  int include_powerimbalance_variables =
      opflow->include_powerimbalance_variables;

  PetscFunctionBegin;
  //  PetscPrintf(MPI_COMM_SELF,"Entered Equality constraints\n");

  // Zero out array
  auto &resmgr = umpire::ResourceManager::getInstance();
  resmgr.memset(ge_dev, 0, opflow->nconeq * sizeof(double));

  /* Generator contributions */
  int *g_gidxbus = genparams->gidxbus_dev_;
  int *g_xidx = genparams->xidx_dev_;
  RAJA::forall<exago_raja_exec>(
      RAJA::RangeSegment(0, genparams->ngenON),
      RAJA_LAMBDA(RAJA::Index_type i) {
        RAJA::atomicSub<exago_raja_atomic>(&ge_dev[g_gidxbus[i]],
                                           x_dev[g_xidx[i]]);
        RAJA::atomicSub<exago_raja_atomic>(&ge_dev[g_gidxbus[i] + 1],
                                           x_dev[g_xidx[i] + 1]);
      });
  flps += genparams->ngenON * 2;

  if (opflow->has_gensetpoint) {
    int *g_geqidxgen = genparams->geqidxgen_dev_;
    double *g_pgs = genparams->pgs_dev_;
    RAJA::forall<exago_raja_exec>(
        RAJA::RangeSegment(0, genparams->ngenON),
        RAJA_LAMBDA(RAJA::Index_type i) {
          double Pg, delPg, Pgset;
          Pg = x_dev[g_xidx[i]];
          delPg = x_dev[g_xidx[i] + 2];
          Pgset = x_dev[g_xidx[i] + 3];

          ge_dev[g_geqidxgen[i]] = Pgset + delPg - Pg;
          ge_dev[g_geqidxgen[i] + 1] = Pgset - g_pgs[i];
        });
  }

  /* Load contributions */
  double *pl = loadparams->pl_dev_;
  double *ql = loadparams->ql_dev_;
  int *l_gidx = loadparams->gidx_dev_;
  int *l_xidx = loadparams->xidx_dev_;
  RAJA::forall<exago_raja_exec>(
      RAJA::RangeSegment(0, loadparams->nload),
      RAJA_LAMBDA(RAJA::Index_type i) {
        if (include_loadloss_variables) {
          RAJA::atomicAdd<exago_raja_atomic>(&ge_dev[l_gidx[i]], pl[i]);
          RAJA::atomicAdd<exago_raja_atomic>(&ge_dev[l_gidx[i] + 1], ql[i]);
          RAJA::atomicSub<exago_raja_atomic>(&ge_dev[l_gidx[i]],
                                             x_dev[l_xidx[i]]);
          RAJA::atomicSub<exago_raja_atomic>(&ge_dev[l_gidx[i] + 1],
                                             x_dev[l_xidx[i] + 1]);
        } else {
          RAJA::atomicAdd<exago_raja_atomic>(&ge_dev[l_gidx[i]], pl[i]);
          RAJA::atomicAdd<exago_raja_atomic>(&ge_dev[l_gidx[i] + 1], ql[i]);
        }
      });
  flps += loadparams->nload * 2;

  /* Bus contributions */
  int *isisolated = busparams->isisolated_dev_;
  int *ispvpq = busparams->ispvpq_dev_;
  double *gl = busparams->gl_dev_;
  double *bl = busparams->bl_dev_;
  double *vm = busparams->vm_dev_;
  double *va = busparams->va_dev_;
  int *b_xidx = busparams->xidx_dev_;
  int *b_gidx = busparams->gidx_dev_;
  int *b_xidxpimb = busparams->xidxpimb_dev_;

  RAJA::forall<exago_raja_exec>(
      RAJA::RangeSegment(0, busparams->nbus), RAJA_LAMBDA(RAJA::Index_type i) {
        double theta = x_dev[b_xidx[i]];
        double Vm = x_dev[b_xidx[i] + 1];
        RAJA::atomicAdd<exago_raja_atomic>(
            &ge_dev[b_gidx[i]],
            isisolated[i] * (theta - va[i] * PETSC_PI / 180.0) +
                ispvpq[i] * Vm * Vm * gl[i]);

        RAJA::atomicAdd<exago_raja_atomic>(&ge_dev[b_gidx[i] + 1],
                                           isisolated[i] * (Vm - vm[i]) -
                                               ispvpq[i] * Vm * Vm * bl[i]);
        if (include_powerimbalance_variables) {
          double Pimb = x_dev[b_xidxpimb[i]];
          double Qimb = x_dev[b_xidxpimb[i] + 1];
          RAJA::atomicAdd<exago_raja_atomic>(&ge_dev[b_gidx[i]], Pimb);
          RAJA::atomicAdd<exago_raja_atomic>(&ge_dev[b_gidx[i] + 1], Qimb);
        }
      });
  flps += busparams->nbus * 14;

  /* Line contributions */
  double *Gff = lineparams->Gff_dev_;
  double *Gtt = lineparams->Gtt_dev_;
  double *Gft = lineparams->Gft_dev_;
  double *Gtf = lineparams->Gtf_dev_;

  double *Bff = lineparams->Bff_dev_;
  double *Btt = lineparams->Btt_dev_;
  double *Bft = lineparams->Bft_dev_;
  double *Btf = lineparams->Btf_dev_;

  int *xidxf = lineparams->xidxf_dev_;
  int *xidxt = lineparams->xidxt_dev_;
  int *geqidxf = lineparams->geqidxf_dev_;
  int *geqidxt = lineparams->geqidxt_dev_;

  RAJA::forall<exago_raja_exec>(
      RAJA::RangeSegment(0, lineparams->nlineON),
      RAJA_LAMBDA(RAJA::Index_type i) {
        double Pf, Qf, Pt, Qt;
        double thetaf = x_dev[xidxf[i]], Vmf = x_dev[xidxf[i] + 1];
        double thetat = x_dev[xidxt[i]], Vmt = x_dev[xidxt[i] + 1];
        double thetaft = thetaf - thetat;
        double thetatf = thetat - thetaf;

        Pf = Gff[i] * Vmf * Vmf +
             Vmf * Vmt * (Gft[i] * cos(thetaft) + Bft[i] * sin(thetaft));
        Qf = -Bff[i] * Vmf * Vmf +
             Vmf * Vmt * (-Bft[i] * cos(thetaft) + Gft[i] * sin(thetaft));
        Pt = Gtt[i] * Vmt * Vmt +
             Vmt * Vmf * (Gtf[i] * cos(thetatf) + Btf[i] * sin(thetatf));
        Qt = -Btt[i] * Vmt * Vmt +
             Vmt * Vmf * (-Btf[i] * cos(thetatf) + Gtf[i] * sin(thetatf));

        RAJA::atomicAdd<exago_raja_atomic>(&ge_dev[geqidxf[i]], Pf);
        RAJA::atomicAdd<exago_raja_atomic>(&ge_dev[geqidxf[i] + 1], Qf);
        RAJA::atomicAdd<exago_raja_atomic>(&ge_dev[geqidxt[i]], Pt);
        RAJA::atomicAdd<exago_raja_atomic>(&ge_dev[geqidxt[i] + 1], Qt);
      });
  flps += lineparams->nlineON *
          (46 + (4 * EXAGO_FLOPS_COSOP) + (4 * EXAGO_FLOPS_SINOP));
  //  PetscPrintf(MPI_COMM_SELF,"Exit Equality Constraints\n");

  ierr = PetscLogFlops(flps);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/** INEQUALITY CONSTRAINTS **/
PetscErrorCode OPFLOWComputeInequalityConstraintsArray_PBPOLRAJAHIOPSPARSE(
    OPFLOW opflow, const double *x_dev, double *gi_dev) {
  PbpolModelRajaHiop *pbpolrajahiopsparse =
      reinterpret_cast<PbpolModelRajaHiop *>(opflow->model);
  LINEParamsRajaHiop *lineparams = &pbpolrajahiopsparse->lineparams;
  PetscInt flps = 0;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  //  PetscPrintf(MPI_COMM_SELF,"Entered Inequality Constraints\n");

  // Zero out array
  auto &resmgr = umpire::ResourceManager::getInstance();
  resmgr.memset(gi_dev, 0, opflow->nconineq * sizeof(double));

  /* Line contributions */
  double *Gff = lineparams->Gff_dev_;
  double *Gtt = lineparams->Gtt_dev_;
  double *Gft = lineparams->Gft_dev_;
  double *Gtf = lineparams->Gtf_dev_;

  double *Bff = lineparams->Bff_dev_;
  double *Btt = lineparams->Btt_dev_;
  double *Bft = lineparams->Bft_dev_;
  double *Btf = lineparams->Btf_dev_;

  int *linelimidx = lineparams->linelimidx_dev_;
  int *xidxf = lineparams->xidxf_dev_;
  int *xidxt = lineparams->xidxt_dev_;
  int *gineqidx = lineparams->gineqidx_dev_;
  if (lineparams->nlinelim) {
    RAJA::forall<exago_raja_exec>(
        RAJA::RangeSegment(0, lineparams->nlinelim),
        RAJA_LAMBDA(RAJA::Index_type i) {
          int j = linelimidx[i];
          double Pf, Qf, Pt, Qt, Sf2, St2;
          double thetaf = x_dev[xidxf[j]], Vmf = x_dev[xidxf[j] + 1];
          double thetat = x_dev[xidxt[j]], Vmt = x_dev[xidxt[j] + 1];
          double thetaft = thetaf - thetat;
          double thetatf = thetat - thetaf;

          Pf = Gff[j] * Vmf * Vmf +
               Vmf * Vmt * (Gft[j] * cos(thetaft) + Bft[j] * sin(thetaft));
          Qf = -Bff[j] * Vmf * Vmf +
               Vmf * Vmt * (-Bft[j] * cos(thetaft) + Gft[j] * sin(thetaft));
          Pt = Gtt[j] * Vmt * Vmt +
               Vmt * Vmf * (Gtf[j] * cos(thetatf) + Btf[j] * sin(thetatf));
          Qt = -Btt[j] * Vmt * Vmt +
               Vmt * Vmf * (-Btf[j] * cos(thetatf) + Gtf[j] * sin(thetatf));

          Sf2 = Pf * Pf + Qf * Qf;
          St2 = Pt * Pt + Qt * Qt;

          gi_dev[gineqidx[i]] = Sf2;
          gi_dev[gineqidx[i] + 1] = St2;
        });
    flps += lineparams->nlinelim *
            (72 + (4 * EXAGO_FLOPS_COSOP) + (4 * EXAGO_FLOPS_SINOP));
  }
  //  PetscPrintf(MPI_COMM_SELF,"Exit Inequality Constraints\n");
  ierr = PetscLogFlops(flps);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/** OBJECTIVE FUNCTION **/
// Note: This kernel (and all the kernels for this model assume that the data
// has been already allocated on the device. x_dev is pointer to array on the
// GPU
PetscErrorCode OPFLOWComputeObjectiveArray_PBPOLRAJAHIOPSPARSE(
    OPFLOW opflow, const double *x_dev, double *obj) {
  PbpolModelRajaHiop *pbpolrajahiopsparse =
      reinterpret_cast<PbpolModelRajaHiop *>(opflow->model);
  GENParamsRajaHiop *genparams = &pbpolrajahiopsparse->genparams;
  LOADParamsRajaHiop *loadparams = &pbpolrajahiopsparse->loadparams;
  BUSParamsRajaHiop *busparams = &pbpolrajahiopsparse->busparams;
  PetscErrorCode ierr;
  PS ps = opflow->ps;
  int isobj_gencost = opflow->obj_gencost;
  double MVAbase = ps->MVAbase;

  PetscFunctionBegin;
  //  PetscPrintf(MPI_COMM_SELF,"Entered objective function\n");

  // You need local copies of device pointers so that lambda can capture them
  double *cost_alpha = genparams->cost_alpha_dev_;
  double *cost_beta = genparams->cost_beta_dev_;
  double *cost_gamma = genparams->cost_gamma_dev_;
  int *xidx = genparams->xidx_dev_;
  int *l_xidx = loadparams->xidx_dev_;
  int *b_xidxpimb = busparams->xidxpimb_dev_;

  /* Generator objective function contributions */
  // Set up reduce sum object
  RAJA::ReduceSum<exago_raja_reduce, double> obj_val_sum(0.0);
  // Compute reduction on CUDA device
  RAJA::forall<exago_raja_exec>(
      RAJA::RangeSegment(0, genparams->ngenON),
      RAJA_LAMBDA(RAJA::Index_type i) {
        double Pg = x_dev[xidx[i]] * MVAbase;
        obj_val_sum += isobj_gencost * (cost_alpha[i] * Pg * Pg +
                                        cost_beta[i] * Pg + cost_gamma[i]);
      });

  if (opflow->include_loadloss_variables) {
    RAJA::forall<exago_raja_exec>(
        RAJA::RangeSegment(0, loadparams->nload),
        RAJA_LAMBDA(RAJA::Index_type i) {
          double Pdloss = x_dev[l_xidx[i]];
          double Qdloss = x_dev[l_xidx[i] + 1];
          obj_val_sum += loadparams->loadloss_penalty_dev_[i] * ps->MVAbase *
                         ps->MVAbase * (Pdloss * Pdloss + Qdloss * Qdloss);
        });
  }

  /* Powerimbalance contributions */
  if (opflow->include_powerimbalance_variables) {
    RAJA::forall<exago_raja_exec>(
        RAJA::RangeSegment(0, busparams->nbus),
        RAJA_LAMBDA(RAJA::Index_type i) {
          double Pimb = x_dev[b_xidxpimb[i]];
          double Qimb = x_dev[b_xidxpimb[i] + 1];
          obj_val_sum += busparams->powerimbalance_penalty_dev_[i] *
                         ps->MVAbase * ps->MVAbase *
                         (Pimb * Pimb + Qimb * Qimb);
        });
  }

  *obj = static_cast<double>(obj_val_sum.get());
  ierr = PetscLogFlops(genparams->ngenON * 8.0);
  CHKERRQ(ierr);

  //  PetscPrintf(MPI_COMM_SELF,"Exit objective function\n");
  PetscFunctionReturn(0);
}

/** GRADIENT **/
PetscErrorCode OPFLOWComputeGradientArray_PBPOLRAJAHIOPSPARSE(
    OPFLOW opflow, const double *x_dev, double *grad_dev) {
  PbpolModelRajaHiop *pbpolrajahiopsparse =
      reinterpret_cast<PbpolModelRajaHiop *>(opflow->model);
  GENParamsRajaHiop *genparams = &pbpolrajahiopsparse->genparams;
  LOADParamsRajaHiop *loadparams = &pbpolrajahiopsparse->loadparams;
  BUSParamsRajaHiop *busparams = &pbpolrajahiopsparse->busparams;
  PS ps = opflow->ps;
  int isobj_gencost = opflow->obj_gencost;
  double MVAbase = ps->MVAbase;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  //  PetscPrintf(MPI_COMM_SELF,"Entered gradient function\n");
  // Zero out array
  auto &resmgr = umpire::ResourceManager::getInstance();
  resmgr.memset(grad_dev, 0);

  /* Generator gradient contributions */
  double *cost_alpha = genparams->cost_alpha_dev_;
  double *cost_beta = genparams->cost_beta_dev_;
  int *xidx = genparams->xidx_dev_;
  int *l_xidx = loadparams->xidx_dev_;
  int *b_xidxpimb = busparams->xidxpimb_dev_;

  RAJA::forall<exago_raja_exec>(
      RAJA::RangeSegment(0, genparams->ngenON),
      RAJA_LAMBDA(RAJA::Index_type i) {
        double Pg = x_dev[xidx[i]] * MVAbase;
        grad_dev[xidx[i]] =
            isobj_gencost * MVAbase * (2.0 * cost_alpha[i] * Pg + cost_beta[i]);
      });

  if (opflow->include_loadloss_variables) {
    RAJA::forall<exago_raja_exec>(
        RAJA::RangeSegment(0, loadparams->nload),
        RAJA_LAMBDA(RAJA::Index_type i) {
          double Pdloss = x_dev[l_xidx[i]];
          double Qdloss = x_dev[l_xidx[i] + 1];
          grad_dev[l_xidx[i]] = loadparams->loadloss_penalty_dev_[i] *
                                ps->MVAbase * ps->MVAbase * 2 * Pdloss;
          grad_dev[l_xidx[i] + 1] = loadparams->loadloss_penalty_dev_[i] *
                                    ps->MVAbase * ps->MVAbase * 2 * Qdloss;
        });
  }

  if (opflow->include_powerimbalance_variables) {
    RAJA::forall<exago_raja_exec>(
        RAJA::RangeSegment(0, busparams->nbus),
        RAJA_LAMBDA(RAJA::Index_type i) {
          double Pimb = x_dev[b_xidxpimb[i]];
          double Qimb = x_dev[b_xidxpimb[i] + 1];
          grad_dev[b_xidxpimb[i]] = busparams->powerimbalance_penalty_dev_[i] *
                                    ps->MVAbase * ps->MVAbase * 2 * Pimb;
          grad_dev[b_xidxpimb[i] + 1] =
              busparams->powerimbalance_penalty_dev_[i] * ps->MVAbase *
              ps->MVAbase * 2 * Qimb;
        });
  }

  ierr = PetscLogFlops(genparams->ngenON * 6.0);
  CHKERRQ(ierr);
  //  PetscPrintf(MPI_COMM_SELF,"Exit gradient function\n");

  PetscFunctionReturn(0);
}

// A routine to sort triplet matrix indexes.  The index arrays are on
// the device.  This sorts them on the host.  
static void
SortIndexes(const int& n, int *i_dev, int *j_dev, int *idx_perm_dev)
{
  auto &resmgr = umpire::ResourceManager::getInstance();
  umpire::Allocator h_allocator_ = resmgr.getAllocator("HOST");

  std::vector< std::tuple<int, int, int> > idxvect;
  idxvect.reserve(n);

  int *itemp(NULL);
  itemp = (int *)(h_allocator_.allocate(n * sizeof(int)));
  resmgr.copy(itemp, i_dev, n*sizeof(int));

  int *jtemp(NULL);
  jtemp = (int *)(h_allocator_.allocate(n * sizeof(int)));
  resmgr.copy(jtemp, j_dev, n*sizeof(int));

  for (int idx = 0; idx < n; idx++) {
    idxvect.push_back(std::make_tuple(itemp[idx], jtemp[idx], idx));
  }

  std::sort(idxvect.begin(), idxvect.end(),
            [] (std::tuple<int,int,int> const &t1,
                std::tuple<int,int,int> const &t2) {
              if (std::get<0>(t1) == std::get<0>(t2)) {
                return (std::get<1>(t1) < std::get<1>(t2));
              } 
              return (std::get<0>(t1) < std::get<0>(t2));
            });

  int *idx_perm;
  idx_perm = (int *)(h_allocator_.allocate(n * sizeof(int)));

  for (int idx = 0; idx < n; idx++) {
    itemp[idx] = std::get<0>(idxvect[idx]);
    jtemp[idx] = std::get<1>(idxvect[idx]);
    int i(std::get<2>(idxvect[idx]));
    idx_perm[i] = idx;
  }

  // std::cout << "Permuted Indexes: " << std::endl;
  // for (int idx = 0; idx < n; idx++) {
  //   std::cout << std::setw(5) << std::right << itemp[idx] << " "
  //             << std::setw(5) << std::right << jtemp[idx] << " "
  //             << std::setw(5) << std::right << idx_perm[idx] << " "
  //             << std::endl;
  // }
  
   resmgr.copy(i_dev, itemp, n*sizeof(int));
   resmgr.copy(j_dev, jtemp, n*sizeof(int));
   resmgr.copy(idx_perm_dev, idx_perm, n*sizeof(int));

   h_allocator_.deallocate(itemp);
   h_allocator_.deallocate(jtemp);
   h_allocator_.deallocate(idx_perm);
}

// A routine to get the triplet arrays from the device and print them out
static void
PrintTriplets(const std::string& title, const int& n, int *iperm,
              int *i, int *j, double *v)
{
  auto &resmgr = umpire::ResourceManager::getInstance();
  umpire::Allocator h_allocator_ = resmgr.getAllocator("HOST");

  int *ipermtemp(NULL);

  if (iperm != NULL) {
    ipermtemp = (int *)(h_allocator_.allocate(n * sizeof(int)));
    resmgr.copy(ipermtemp, iperm, n*sizeof(int));
  }

  int *itemp(NULL);

  if (i != NULL) {
    itemp = (int *)(h_allocator_.allocate(n * sizeof(int)));
    resmgr.copy(itemp, i, n*sizeof(int));
  }

  int *jtemp(NULL);
  if (j != 0) {
    jtemp = (int *)(h_allocator_.allocate(n * sizeof(int)));
    resmgr.copy(jtemp, j, n*sizeof(int));
  }

  double *vtemp(NULL);
  if (v != NULL) {
    vtemp = (double *)(h_allocator_.allocate(n * sizeof(double)));
    resmgr.copy(vtemp, v, n*sizeof(double));
  }

  std::cout << title << std::endl;
  for (int idx = 0; idx < n; ++idx) {
    std::cout << std::setw(5) << idx << " ";
    if (ipermtemp != NULL) {
      std::cout << std::setw(5) << std::right << ipermtemp[idx] << " ";
    }
    if (itemp != NULL) {
      std::cout << std::setw(5) << std::right << itemp[idx] << " ";
    }
    if (jtemp != NULL) {
      std::cout << std::setw(5) << std::right << jtemp[idx];
    }
    if (vtemp != NULL) {
      std::cout << std::setw(12) << std::right
                << std::scientific << std::setprecision(3)
                << vtemp[idx];
    }
    std::cout << std::endl;
  }
  h_allocator_.deallocate(itemp);
  h_allocator_.deallocate(jtemp);
  if (vtemp != NULL)   h_allocator_.deallocate(vtemp);
}
  

PetscErrorCode
OPFLOWComputeSparseInequalityConstraintJacobian_PBPOLRAJAHIOPSPARSE(
    OPFLOW opflow, const double *x_dev, int *iJacS_dev, int *jJacS_dev,
    double *MJacS_dev) {
  PbpolModelRajaHiop *pbpolrajahiopsparse =
      reinterpret_cast<PbpolModelRajaHiop *>(opflow->model);
  PetscErrorCode ierr;
  double *x, *values;
  PetscInt *iRowstart, *jColstart;
  PetscInt roffset, coffset, idxoffset;
  PetscInt nrow, ncol;
  PetscInt nvals;
  const PetscInt *cols;
  const PetscScalar *vals;
  PetscInt i, j;
  auto &resmgr = umpire::ResourceManager::getInstance();

  PetscFunctionBegin;

  // Create arrays on host to store i,j, and val arrays
  umpire::Allocator h_allocator_ = resmgr.getAllocator("HOST");
  umpire::Allocator d_allocator_ = resmgr.getAllocator("DEVICE");
  

  if (MJacS_dev == NULL) {

    if (debugmsg)
      std::cout << "Official Inequality Jacobian nonzero count: "
                << opflow->nnz_ineqjacsp << std::endl;
    
    /* Set locations only */

    if (opflow->Nconineq) {
      ierr = PetscLogEventBegin(opflow->ineqconsjaclogger, 0, 0, 0, 0);

      /* Inequality constraints start after equality constraints
         Hence the offset
      */
      roffset = opflow->nconeq;
      coffset = 0;
      idxoffset = opflow->nnz_eqjacsp;

      resmgr.memset(iJacS_dev + idxoffset, 0,
                    opflow->nnz_ineqjacsp*sizeof(int));
      resmgr.memset(jJacS_dev + idxoffset, 0,
                    opflow->nnz_ineqjacsp*sizeof(int));

      if (!opflow->ignore_lineflow_constraints) {
        LINEParamsRajaHiop *lineparams = &pbpolrajahiopsparse->lineparams;
        int *jac_ieq_idx = lineparams->jac_ieq_idx_dev_;
        int *linelimidx = lineparams->linelimidx_dev_;
        int *xidxf = lineparams->xidxf_dev_;
        int *xidxt = lineparams->xidxt_dev_;
        int *gbineqidx = lineparams->gbineqidx_dev_;

        RAJA::forall<exago_raja_exec>(
            RAJA::RangeSegment(0, lineparams->nlinelim),
            RAJA_LAMBDA(RAJA::Index_type i) {
              int iline(linelimidx[i]);
              int offset(0);
              
              iJacS_dev[jac_ieq_idx[i] + offset] = gbineqidx[i];
              jJacS_dev[jac_ieq_idx[i] + offset] = xidxf[iline];
              offset++;
              
              iJacS_dev[jac_ieq_idx[i] + offset] = gbineqidx[i];
              jJacS_dev[jac_ieq_idx[i] + offset] = xidxf[iline] + 1;
              offset++;
              
              iJacS_dev[jac_ieq_idx[i] + offset] = gbineqidx[i];
              jJacS_dev[jac_ieq_idx[i] + offset] = xidxt[iline];
              offset++;

              iJacS_dev[jac_ieq_idx[i] + offset] = gbineqidx[i];
              jJacS_dev[jac_ieq_idx[i] + offset] = xidxt[iline] + 1;
              offset++;
              
              iJacS_dev[jac_ieq_idx[i] + offset] = gbineqidx[i] + 1;
              jJacS_dev[jac_ieq_idx[i] + offset] = xidxf[iline];
              offset++;
              
              iJacS_dev[jac_ieq_idx[i] + offset] = gbineqidx[i] + 1;
              jJacS_dev[jac_ieq_idx[i] + offset] = xidxf[iline] + 1;
              offset++;
              
              iJacS_dev[jac_ieq_idx[i] + offset] = gbineqidx[i] + 1;
              jJacS_dev[jac_ieq_idx[i] + offset] = xidxt[iline];
              offset++;

              iJacS_dev[jac_ieq_idx[i] + offset] = gbineqidx[i] + 1;
              jJacS_dev[jac_ieq_idx[i] + offset] = xidxt[iline] + 1;
            });
        
      }

      if (pbpolrajahiopsparse->idx_jacineq_dev_ == NULL) {
        pbpolrajahiopsparse->idx_jacineq_dev_ =
          (int *) d_allocator_.allocate(opflow->nnz_ineqjacsp * sizeof(int));
      }

      SortIndexes(opflow->nnz_ineqjacsp,
                  iJacS_dev + opflow->nnz_eqjacsp,
                  jJacS_dev + opflow->nnz_eqjacsp,
                  pbpolrajahiopsparse->idx_jacineq_dev_);
      
      if (debugmsg)
        PrintTriplets("Nonzero indexes for Inequality Constraint Jacobian (GPU):",
                      opflow->nnz_ineqjacsp,
                      pbpolrajahiopsparse->idx_jacineq_dev_,
                      iJacS_dev + opflow->nnz_eqjacsp,
                      jJacS_dev + opflow->nnz_eqjacsp,
                      NULL);

      if (oldhostway) {
        
      /* Inequality constraints start after equality constraints
         Hence the offset
      */
      roffset = opflow->nconeq;
      coffset = 0;

      if (pbpolrajahiopsparse->i_jacineq == NULL) {
        pbpolrajahiopsparse->i_jacineq =
          (int *)(h_allocator_.allocate(opflow->nnz_ineqjacsp * sizeof(int)));
        pbpolrajahiopsparse->j_jacineq =
          (int *)(h_allocator_.allocate(opflow->nnz_ineqjacsp * sizeof(int)));
        pbpolrajahiopsparse->val_jacineq =
          (double *)(h_allocator_.allocate(opflow->nnz_ineqjacsp * sizeof(double)));
      }

      iRowstart = pbpolrajahiopsparse->i_jacineq;
      jColstart = pbpolrajahiopsparse->j_jacineq;

      ierr = (*opflow->modelops.computeinequalityconstraintjacobian)(
          opflow, opflow->X, opflow->Jac_Gi);
      CHKERRQ(ierr);

      ierr = MatGetSize(opflow->Jac_Gi, &nrow, &ncol);
      CHKERRQ(ierr);
      /* Copy over locations to triplet format */
      for (i = 0; i < nrow; i++) {
        ierr = MatGetRow(opflow->Jac_Gi, i, &nvals, &cols, &vals);
        CHKERRQ(ierr);
        for (j = 0; j < nvals; j++) {
          iRowstart[j] = roffset + i;
          jColstart[j] = coffset + cols[j];
        }
        /* Increment iRow,jCol pointers */
        iRowstart += nvals;
        jColstart += nvals;
        ierr = MatRestoreRow(opflow->Jac_Gi, i, &nvals, &cols, &vals);
        CHKERRQ(ierr);
      }

      // Copy over i_jacineq and j_jacineq arrays to device
      resmgr.copy(iJacS_dev + opflow->nnz_eqjacsp,
                  pbpolrajahiopsparse->i_jacineq);
      resmgr.copy(jJacS_dev + opflow->nnz_eqjacsp,
                  pbpolrajahiopsparse->j_jacineq);

      if (debugmsg)
        PrintTriplets("Nonzero indexes for Inequality Constraint Jacobian:",
                      opflow->nnz_ineqjacsp,
                      NULL,
                      iJacS_dev + opflow->nnz_eqjacsp,
                      jJacS_dev + opflow->nnz_eqjacsp,
                      NULL);

      ierr = PetscLogEventEnd(opflow->ineqconsjaclogger, 0, 0, 0, 0);
      CHKERRQ(ierr);
      }
    }
  } else {
    if (opflow->Nconineq) {
      ierr = PetscLogEventBegin(opflow->ineqconsjaclogger, 0, 0, 0, 0);
      CHKERRQ(ierr);

      int *iperm = pbpolrajahiopsparse->idx_jacineq_dev_;

      if (!opflow->ignore_lineflow_constraints) {
        LINEParamsRajaHiop *lineparams = &pbpolrajahiopsparse->lineparams;
        double *Gff_arr = lineparams->Gff_dev_;
        double *Gtt_arr = lineparams->Gtt_dev_;
        double *Gft_arr = lineparams->Gft_dev_;
        double *Gtf_arr = lineparams->Gtf_dev_;
        
        double *Bff_arr = lineparams->Bff_dev_;
        double *Btt_arr = lineparams->Btt_dev_;
        double *Bft_arr = lineparams->Bft_dev_;
        double *Btf_arr = lineparams->Btf_dev_;
      
        int *linelimidx = lineparams->linelimidx_dev_;
        int *xidxf = lineparams->xidxf_dev_;
        int *xidxt = lineparams->xidxt_dev_;
        int *jac_ieq_idx = lineparams->jac_ieq_idx_dev_;
        RAJA::forall<exago_raja_exec>(
            RAJA::RangeSegment(0, lineparams->nlinelim),
            RAJA_LAMBDA(RAJA::Index_type i) {
              int j = linelimidx[i];
              double val[4];
              double Pf, Qf, Pt, Qt;
              double thetaf = x_dev[xidxf[j]], Vmf = x_dev[xidxf[j] + 1];
              double thetat = x_dev[xidxt[j]], Vmt = x_dev[xidxt[j] + 1];
              double thetaft = thetaf - thetat;
              double thetatf = thetat - thetaf;
              double dSf2_dPf, dSf2_dQf, dSt2_dPt, dSt2_dQt;
              double dPf_dthetaf, dPf_dVmf, dPf_dthetat, dPf_dVmt;
              double dQf_dthetaf, dQf_dVmf, dQf_dthetat, dQf_dVmt;
              double dPt_dthetaf, dPt_dVmf, dPt_dthetat, dPt_dVmt;
              double dQt_dthetaf, dQt_dVmf, dQt_dthetat, dQt_dVmt;
              double dSf2_dthetaf, dSf2_dVmf, dSf2_dthetat, dSf2_dVmt;
              double dSt2_dthetaf, dSt2_dVmf, dSt2_dthetat, dSt2_dVmt;
              double Gff = Gff_arr[j], Bff = Bff_arr[j];
              double Gft = Gft_arr[j], Bft = Bft_arr[j];
              double Gtf = Gtf_arr[j], Btf = Btf_arr[j];
              double Gtt = Gtt_arr[j], Btt = Btt_arr[j];
              
              Pf = Gff * Vmf * Vmf +
                Vmf * Vmt * (Gft * cos(thetaft) + Bft * sin(thetaft));
              Qf = -Bff * Vmf * Vmf +
                Vmf * Vmt * (-Bft * cos(thetaft) + Gft * sin(thetaft));
              Pt = Gtt * Vmt * Vmt +
                Vmt * Vmf * (Gtf * cos(thetatf) + Btf * sin(thetatf));
              Qt = -Btt * Vmt * Vmt +
                Vmt * Vmf * (-Btf * cos(thetatf) + Gtf * sin(thetatf));
              
              dSf2_dPf = 2 * Pf;
              dSf2_dQf = 2 * Qf;
              dSt2_dPt = 2 * Pt;
              dSt2_dQt = 2 * Qt;
              
              dPf_dthetaf = Vmf * Vmt * (-Gft * sin(thetaft) + Bft * cos(thetaft));
              dPf_dVmf =
                2 * Gff * Vmf + Vmt * (Gft * cos(thetaft) + Bft * sin(thetaft));
              dPf_dthetat = Vmf * Vmt * (Gft * sin(thetaft) - Bft * cos(thetaft));
              dPf_dVmt = Vmf * (Gft * cos(thetaft) + Bft * sin(thetaft));
              
              dQf_dthetaf = Vmf * Vmt * (Bft * sin(thetaft) + Gft * cos(thetaft));
              dQf_dVmf =
                -2 * Bff * Vmf + Vmt * (-Bft * cos(thetaft) + Gft * sin(thetaft));
              dQf_dthetat = Vmf * Vmt * (-Bft * sin(thetaft) - Gft * cos(thetaft));
              dQf_dVmt = Vmf * (-Bft * cos(thetaft) + Gft * sin(thetaft));
              
              dPt_dthetat = Vmt * Vmf * (-Gtf * sin(thetatf) + Btf * cos(thetatf));
              dPt_dVmt =
                2 * Gtt * Vmt + Vmf * (Gtf * cos(thetatf) + Btf * sin(thetatf));
              dPt_dthetaf = Vmt * Vmf * (Gtf * sin(thetatf) - Btf * cos(thetatf));
              dPt_dVmf = Vmt * (Gtf * cos(thetatf) + Btf * sin(thetatf));
              
              dQt_dthetat = Vmt * Vmf * (Btf * sin(thetatf) + Gtf * cos(thetatf));
              dQt_dVmt =
                -2 * Btt * Vmt + Vmf * (-Btf * cos(thetatf) + Gtf * sin(thetatf));
              dQt_dthetaf = Vmt * Vmf * (-Btf * sin(thetatf) - Gtf * cos(thetatf));
              dQt_dVmf = Vmt * (-Btf * cos(thetatf) + Gtf * sin(thetatf));
              
              dSf2_dthetaf = dSf2_dPf * dPf_dthetaf + dSf2_dQf * dQf_dthetaf;
              dSf2_dthetat = dSf2_dPf * dPf_dthetat + dSf2_dQf * dQf_dthetat;
              dSf2_dVmf = dSf2_dPf * dPf_dVmf + dSf2_dQf * dQf_dVmf;
              dSf2_dVmt = dSf2_dPf * dPf_dVmt + dSf2_dQf * dQf_dVmt;

              val[0] = dSf2_dthetaf;
              val[1] = dSf2_dVmf;
              val[2] = dSf2_dthetat;
              val[3] = dSf2_dVmt;

              MJacS_dev[jac_ieq_idx[i] + 0] = val[0];
              MJacS_dev[jac_ieq_idx[i] + 1] = val[1];
              MJacS_dev[jac_ieq_idx[i] + 2] = val[2];
              MJacS_dev[jac_ieq_idx[i] + 3] = val[3];

              dSt2_dthetaf = dSt2_dPt * dPt_dthetaf + dSt2_dQt * dQt_dthetaf;
              dSt2_dthetat = dSt2_dPt * dPt_dthetat + dSt2_dQt * dQt_dthetat;
              dSt2_dVmf = dSt2_dPt * dPt_dVmf + dSt2_dQt * dQt_dVmf;
              dSt2_dVmt = dSt2_dPt * dPt_dVmt + dSt2_dQt * dQt_dVmt;
              
              val[2] = dSt2_dthetat;
              val[3] = dSt2_dVmt;
              val[0] = dSt2_dthetaf;
              val[1] = dSt2_dVmf;

              MJacS_dev[jac_ieq_idx[i] + 4] = val[0];
              MJacS_dev[jac_ieq_idx[i] + 5] = val[1];
              MJacS_dev[jac_ieq_idx[i] + 6] = val[2];
              MJacS_dev[jac_ieq_idx[i] + 7] = val[3];

            });

      }

      int *ipermout = (int *)d_allocator_.allocate(opflow->nnz_ineqjacsp*sizeof(int));
      resmgr.copy(ipermout, iperm);

      RAJA::stable_sort_pairs<exago_raja_exec>
        (RAJA::make_span(ipermout, opflow->nnz_ineqjacsp),
         RAJA::make_span(MJacS_dev + opflow->nnz_eqjacsp, opflow->nnz_ineqjacsp),
         RAJA::operators::less<int>{});

      if (debugmsg)
        PrintTriplets("Inequality Constraint Jacobian (GPU):",
                      opflow->nnz_ineqjacsp,
                      iperm, 
                      (iJacS_dev == NULL ? NULL : iJacS_dev + opflow->nnz_eqjacsp),
                      (jJacS_dev == NULL ? NULL : jJacS_dev + opflow->nnz_eqjacsp),
                      MJacS_dev);

      d_allocator_.deallocate(ipermout);
      
      if (oldhostway) {

      ierr = VecGetArray(opflow->X, &x);
      CHKERRQ(ierr);

      // Copy from device to host
      umpire::Allocator h_allocator_ = resmgr.getAllocator("HOST");
      registerWith(x, opflow->nx, resmgr, h_allocator_);
      resmgr.copy((double *)x, (double *)x_dev);

      ierr = VecRestoreArray(opflow->X, &x);
      CHKERRQ(ierr);

      /* Compute inequality constraint jacobian */
      ierr = (*opflow->modelops.computeinequalityconstraintjacobian)(
          opflow, opflow->X, opflow->Jac_Gi);
      CHKERRQ(ierr);

      ierr = MatGetSize(opflow->Jac_Gi, &nrow, &ncol);
      CHKERRQ(ierr);

      values = pbpolrajahiopsparse->val_jacineq;
      /* Copy over values */
      for (i = 0; i < nrow; i++) {
        ierr = MatGetRow(opflow->Jac_Gi, i, &nvals, &cols, &vals);
        CHKERRQ(ierr);
        for (j = 0; j < nvals; j++) {
          values[j] = vals[j];
        }
        values += nvals;
        ierr = MatRestoreRow(opflow->Jac_Gi, i, &nvals, &cols, &vals);
        CHKERRQ(ierr);
      }
      // Copy over val_jacineq to device
      resmgr.copy(MJacS_dev + opflow->nnz_eqjacsp,
                  pbpolrajahiopsparse->val_jacineq);

      if (debugmsg)
        PrintTriplets("Inequality Constraint Jacobian:",
                      opflow->nnz_ineqjacsp,
                      NULL, 
                      (iJacS_dev == NULL ? NULL : iJacS_dev + opflow->nnz_eqjacsp),
                      (jJacS_dev == NULL ? NULL : jJacS_dev + opflow->nnz_eqjacsp),
                      MJacS_dev + opflow->nnz_eqjacsp);

      ierr = PetscLogEventEnd(opflow->ineqconsjaclogger, 0, 0, 0, 0);
      CHKERRQ(ierr);
      }
    }
  }

  PetscFunctionReturn(0);
}


PetscErrorCode
OPFLOWComputeSparseEqualityConstraintJacobian_PBPOLRAJAHIOPSPARSE(
    OPFLOW opflow, const double *x_dev, int *iJacS_dev, int *jJacS_dev,
    double *MJacS_dev) {
  PbpolModelRajaHiop *pbpolrajahiopsparse =
      reinterpret_cast<PbpolModelRajaHiop *>(opflow->model);
  GENParamsRajaHiop *genparams = &pbpolrajahiopsparse->genparams;
  LOADParamsRajaHiop *loadparams = &pbpolrajahiopsparse->loadparams;
  BUSParamsRajaHiop *busparams = &pbpolrajahiopsparse->busparams;
  LINEParamsRajaHiop *lineparams = &pbpolrajahiopsparse->lineparams;

  PetscErrorCode ierr;
  PetscInt *iRowstart, *jColstart;
  PetscScalar *x, *values;
  PetscInt roffset, coffset;
  PetscInt nrow, ncol;
  PetscInt nvals;
  const PetscInt *cols;
  const PetscScalar *vals;
  PetscInt i, j;
  auto &resmgr = umpire::ResourceManager::getInstance();

  PetscFunctionBegin;

  ierr = PetscLogEventBegin(opflow->eqconsjaclogger, 0, 0, 0, 0);
  CHKERRQ(ierr);

  umpire::Allocator h_allocator_ = resmgr.getAllocator("HOST");
  umpire::Allocator d_allocator_ = resmgr.getAllocator("DEVICE");
  
  /* Using OPFLOWComputeEqualityConstraintJacobian_PBPOL() as a guide */

  if (MJacS_dev == NULL) {

    if (debugmsg)
      std::cout << "Official Equality Jacobian nonzero count: "
                << opflow->nnz_eqjacsp << std::endl;
    
    /* Set locations only */

    resmgr.memset(iJacS_dev, 0, opflow->nnz_eqjacsp*sizeof(int));
    resmgr.memset(jJacS_dev, 0, opflow->nnz_eqjacsp*sizeof(int));

    /* Bus power imbalance contribution */
    int *b_xidxpimb = busparams->xidxpimb_dev_;
    int *b_gidx = busparams->gidx_dev_;
    int *b_xidx = busparams->xidx_dev_;
    int *b_jacsp_idx = busparams->jacsp_idx_dev_;
    int *b_jacsq_idx = busparams->jacsq_idx_dev_;

    /* Bus */
    if (debugmsg) std::cout << "Begin with buses" << std::endl;
    RAJA::forall<exago_raja_exec>(
      RAJA::RangeSegment(0, busparams->nbus),
        RAJA_LAMBDA(RAJA::Index_type i) {
          iJacS_dev[b_jacsp_idx[i]] = b_gidx[i];
          jJacS_dev[b_jacsp_idx[i]] = b_xidx[i];
          iJacS_dev[b_jacsp_idx[i] + 1] = b_gidx[i];
          jJacS_dev[b_jacsp_idx[i] + 1] = b_xidx[i] + 1;
          
          iJacS_dev[b_jacsq_idx[i]] = b_gidx[i] + 1;
          jJacS_dev[b_jacsq_idx[i]] = b_xidx[i];
          iJacS_dev[b_jacsq_idx[i] + 1] = b_gidx[i] + 1;
          jJacS_dev[b_jacsq_idx[i] + 1] = b_xidx[i] + 1;
      });

    if (opflow->include_powerimbalance_variables) {
      if (debugmsg) std::cout << "Bus power imbalance variables" << std::endl;
      RAJA::forall<exago_raja_exec>(
        RAJA::RangeSegment(0, busparams->nbus),
          RAJA_LAMBDA(RAJA::Index_type i) {
            iJacS_dev[b_jacsp_idx[i]] = b_gidx[i];
            jJacS_dev[b_jacsp_idx[i]] = b_xidxpimb[i];
            iJacS_dev[b_jacsp_idx[i] + 1] = b_gidx[i];
            jJacS_dev[b_jacsp_idx[i] + 1] = b_xidxpimb[i] + 1;
            
            iJacS_dev[b_jacsq_idx[i]] = b_gidx[i] + 1;
            jJacS_dev[b_jacsq_idx[i]] = b_xidxpimb[i] + 2;
            iJacS_dev[b_jacsq_idx[i] + 1] = b_gidx[i] + 1;
            jJacS_dev[b_jacsq_idx[i] + 1] = b_xidxpimb[i] + 3;
        });
    }

    /* generation contributions */

    if (debugmsg) std::cout << "Generators " << std::endl;
    
    int *g_gidxbus = genparams->gidxbus_dev_;
    int *g_xidx = genparams->xidx_dev_;
    int *eqjacspbus_idx = genparams->eqjacspbus_idx_dev_;
    int *eqjacsqbus_idx = genparams->eqjacsqbus_idx_dev_;
    RAJA::forall<exago_raja_exec>(
        RAJA::RangeSegment(0, genparams->ngenON),
        RAJA_LAMBDA(RAJA::Index_type i) {
          iJacS_dev[eqjacspbus_idx[i]] = g_gidxbus[i];
          jJacS_dev[eqjacspbus_idx[i]] = g_xidx[i];

          iJacS_dev[eqjacsqbus_idx[i]] = g_gidxbus[i] + 1;
          jJacS_dev[eqjacsqbus_idx[i]] = g_xidx[i] + 1;
        });

    /* Loadloss contributions */
    
    if (opflow->include_loadloss_variables) {

      if (debugmsg) std::cout << "Load Loss" << std::endl;
      int *l_gidx = loadparams->gidx_dev_;
      int *l_xidx = loadparams->xidx_dev_;
      int *l_jacsp_idx = loadparams->jacsp_idx_dev_;

      RAJA::forall<exago_raja_exec>(
          RAJA::RangeSegment(0, loadparams->nload),
          RAJA_LAMBDA(RAJA::Index_type i) {
            iJacS_dev[l_jacsp_idx[i]] = l_gidx[i];
            jJacS_dev[l_jacsp_idx[i]] = l_xidx[i];
            iJacS_dev[l_jacsp_idx[i] + 1] = l_gidx[i] + 1;
            jJacS_dev[l_jacsp_idx[i] + 1] = l_xidx[i] + 1;
          });
    }

    /* Connected lines */

    if (debugmsg) std::cout << "Connected Lines" << std::endl;
    
    int *xidxf = lineparams->xidxf_dev_;
    int *xidxt = lineparams->xidxt_dev_;
    int *geqidxf = lineparams->geqidxf_dev_;
    int *geqidxt = lineparams->geqidxt_dev_;
    int *jacf_idx = lineparams->jacf_idx_dev_;
    int *jact_idx = lineparams->jact_idx_dev_;

    RAJA::forall<exago_raja_exec>(
      RAJA::RangeSegment(0, lineparams->nlineON),
      RAJA_LAMBDA(RAJA::Index_type i) {

        int offset;

        offset = 0;
        
        iJacS_dev[jacf_idx[i] + offset] = geqidxf[i];
        jJacS_dev[jacf_idx[i] + offset] = xidxt[i];
        offset++;

        iJacS_dev[jacf_idx[i] + offset] = geqidxf[i];
        jJacS_dev[jacf_idx[i] + offset] = xidxt[i] + 1;
        offset++;

        iJacS_dev[jacf_idx[i] + offset] = geqidxf[i] + 1;
        jJacS_dev[jacf_idx[i] + offset] = xidxt[i];
        offset++;

        iJacS_dev[jacf_idx[i] + offset] = geqidxf[i] + 1;
        jJacS_dev[jacf_idx[i] + offset] = xidxt[i] + 1;
        offset++;

        // to bus indexes

        offset = 0;

        iJacS_dev[jact_idx[i] + offset] = geqidxt[i];
        jJacS_dev[jact_idx[i] + offset] = xidxf[i];
        offset++;

        iJacS_dev[jact_idx[i] + offset] = geqidxt[i];
        jJacS_dev[jact_idx[i] + offset] = xidxf[i] + 1;
        offset++;

        iJacS_dev[jact_idx[i] + offset] = geqidxt[i] + 1;
        jJacS_dev[jact_idx[i] + offset] = xidxf[i];
        offset++;

        iJacS_dev[jact_idx[i] + offset] = geqidxt[i] + 1;
        jJacS_dev[jact_idx[i] + offset] = xidxf[i] + 1;
        offset++;

      });
    
    
    
    if (opflow->has_gensetpoint) {

      if (debugmsg) std::cout << "Generator set point" << std::endl;
      int *eqjacspgen_idx = genparams->eqjacspgen_idx_dev_;
      int *g_geqidxgen = genparams->geqidxgen_dev_;
      int *g_xidx = genparams->xidx_dev_;
      int *g_isrenewable = genparams->isrenewable_dev_;

      RAJA::forall<exago_raja_exec>(
          RAJA::RangeSegment(0, genparams->ngenON),
          RAJA_LAMBDA(RAJA::Index_type i) {
            if (!g_isrenewable[i]) {
              iJacS_dev[eqjacspgen_idx[i]] = g_geqidxgen[i];
              jJacS_dev[eqjacspgen_idx[i]] = g_xidx[i];

              iJacS_dev[eqjacspgen_idx[i] + 1] = g_geqidxgen[i];
              jJacS_dev[eqjacspgen_idx[i] + 1] = g_xidx[i] + 2;

              iJacS_dev[eqjacspgen_idx[i] + 2] = g_geqidxgen[i];
              jJacS_dev[eqjacspgen_idx[i] + 2] = g_xidx[i] + 3;

              iJacS_dev[eqjacspgen_idx[i] + 3] = g_geqidxgen[i] + 1;
              jJacS_dev[eqjacspgen_idx[i] + 3] = g_xidx[i] + 3;
            }
          });
    }

    if (pbpolrajahiopsparse->idx_jaceq_dev_ == NULL) {
      pbpolrajahiopsparse->idx_jaceq_dev_ =
        (int *) d_allocator_.allocate(opflow->nnz_eqjacsp * sizeof(int));
    }

    SortIndexes(opflow->nnz_eqjacsp, iJacS_dev, jJacS_dev,
                pbpolrajahiopsparse->idx_jaceq_dev_);
    
    if (debugmsg)
      PrintTriplets("Non-zero indexes for Equality Constraint Jacobian (GPU):",
                    opflow->nnz_eqjacsp,
                    pbpolrajahiopsparse->idx_jaceq_dev_,
                    iJacS_dev, jJacS_dev, NULL);

    if (oldhostway) {
    roffset = 0;
    coffset = 0;

    if (pbpolrajahiopsparse->i_jaceq == NULL) {
      pbpolrajahiopsparse->i_jaceq =
        (int *)(h_allocator_.allocate(opflow->nnz_eqjacsp * sizeof(int)));
      pbpolrajahiopsparse->j_jaceq =
        (int *)(h_allocator_.allocate(opflow->nnz_eqjacsp * sizeof(int)));
      pbpolrajahiopsparse->val_jaceq =
        (double *)(h_allocator_.allocate(opflow->nnz_eqjacsp * sizeof(double)));
    }
    
    iRowstart = pbpolrajahiopsparse->i_jaceq;
    jColstart = pbpolrajahiopsparse->j_jaceq;

    ierr = (*opflow->modelops.computeequalityconstraintjacobian)(
        opflow, opflow->X, opflow->Jac_Ge);
    CHKERRQ(ierr);

    ierr = MatGetSize(opflow->Jac_Ge, &nrow, &ncol);
    CHKERRQ(ierr);

    /* Copy over locations to triplet format */
    for (i = 0; i < nrow; i++) {
      ierr = MatGetRow(opflow->Jac_Ge, i, &nvals, &cols, &vals);
      CHKERRQ(ierr);
      for (j = 0; j < nvals; j++) {
        iRowstart[j] = roffset + i;
        jColstart[j] = coffset + cols[j];
      }
      /* Increment iRow,jCol pointers */
      iRowstart += nvals;
      jColstart += nvals;
      ierr = MatRestoreRow(opflow->Jac_Ge, i, &nvals, &cols, &vals);
      CHKERRQ(ierr);
    }

    // Copy over i_jaceq and j_jaceq arrays to device
    resmgr.copy(iJacS_dev, pbpolrajahiopsparse->i_jaceq);
    resmgr.copy(jJacS_dev, pbpolrajahiopsparse->j_jaceq);

    if (debugmsg)
      PrintTriplets("Non-zero indexes for Equality Constraint Jacobian:",
                    opflow->nnz_eqjacsp, NULL, iJacS_dev, jJacS_dev, NULL);
    }
    
  } else {

    // Bus Contribution
    int *b_jacsp_idx = busparams->jacsp_idx_dev_;
    int *b_jacsq_idx = busparams->jacsq_idx_dev_;
    int *isisolated = busparams->isisolated_dev_;
    int *ispvpq = busparams->ispvpq_dev_;
    double *gl = busparams->gl_dev_;
    double *bl = busparams->bl_dev_;
    int *b_xidx = busparams->xidx_dev_;

    int *iperm = pbpolrajahiopsparse->idx_jaceq_dev_;

    // Basic bus contribution
    RAJA::forall<exago_raja_exec>(
        RAJA::RangeSegment(0, busparams->nbus),
        RAJA_LAMBDA(RAJA::Index_type i) {
          double Vm = x_dev[b_xidx[i] + 1];
          MJacS_dev[b_jacsp_idx[i]] = isisolated[i] * 1.0 + ispvpq[i] * 0.0;
          MJacS_dev[b_jacsp_idx[i]+1] = isisolated[i] * 0.0 + ispvpq[i] * 2 * Vm * gl[i];
          MJacS_dev[b_jacsq_idx[i]] = 0.0;
          MJacS_dev[b_jacsq_idx[i]+1] = isisolated[i] * 1.0 + ispvpq[i] * -2 * Vm * bl[i];
        });
    

    // Power imbalance 
    if (opflow->include_powerimbalance_variables) {
      RAJA::forall<exago_raja_exec>(
          RAJA::RangeSegment(0, busparams->nbus),
          RAJA_LAMBDA(RAJA::Index_type i) {
            MJacS_dev[b_jacsp_idx[i]] = 1.0;
            MJacS_dev[b_jacsp_idx[i] + 1] = -1.0;
            MJacS_dev[b_jacsq_idx[i]] = 1.0;
            MJacS_dev[b_jacsq_idx[i] + 1] = -1.0;
          });
    }

    // Generator contributions 
    int *eqjacspbus_idx = genparams->eqjacspbus_idx_dev_;
    int *eqjacsqbus_idx = genparams->eqjacsqbus_idx_dev_;
    RAJA::forall<exago_raja_exec>(
        RAJA::RangeSegment(0, genparams->ngenON),
        RAJA_LAMBDA(RAJA::Index_type i) {
          MJacS_dev[eqjacspbus_idx[i]] = -1.0;
          MJacS_dev[eqjacsqbus_idx[i]] = -1.0;
        });

    if (opflow->has_gensetpoint) {
      int *eqjacspgen_idx = genparams->eqjacspgen_idx_dev_;
      int *g_isrenewable = genparams->isrenewable_dev_;

      RAJA::forall<exago_raja_exec>(
          RAJA::RangeSegment(0, genparams->ngenON),
          RAJA_LAMBDA(RAJA::Index_type i) {
            if (!g_isrenewable[i]) {
              MJacS_dev[eqjacspgen_idx[i]] = -1.0;
              MJacS_dev[eqjacspgen_idx[i] + 1] = 1.0;
              MJacS_dev[eqjacspgen_idx[i] + 2] = 1.0;
              MJacS_dev[eqjacspgen_idx[i] + 3] = 1.0;
            }
          });
    }

    /* Loadloss contributions - 2 contributions expected */
    if (opflow->include_loadloss_variables) {
      int *l_jacsp_idx = loadparams->jacsp_idx_dev_;
      RAJA::forall<exago_raja_exec>(
          RAJA::RangeSegment(0, loadparams->nload),
          RAJA_LAMBDA(RAJA::Index_type i) {
            MJacS_dev[l_jacsp_idx[i]] = -1;
            MJacS_dev[l_jacsp_idx[i] + 1] = -1;
          });
    }

    // Line contributions

    double *Gff = lineparams->Gff_dev_;
    double *Gtt = lineparams->Gtt_dev_;
    double *Gft = lineparams->Gft_dev_;
    double *Gtf = lineparams->Gtf_dev_;

    double *Bff = lineparams->Bff_dev_;
    double *Btt = lineparams->Btt_dev_;
    double *Bft = lineparams->Bft_dev_;
    double *Btf = lineparams->Btf_dev_;

    int *xidxf = lineparams->xidxf_dev_;
    int *xidxt = lineparams->xidxt_dev_;
    int *busf_idx = lineparams->busf_idx_dev_;
    int *bust_idx = lineparams->bust_idx_dev_;
    int *jacf_idx = lineparams->jacf_idx_dev_;
    int *jact_idx = lineparams->jact_idx_dev_;
    
    RAJA::forall<exago_raja_exec>(
        RAJA::RangeSegment(0, lineparams->nlineON),
        RAJA_LAMBDA(RAJA::Index_type i) {
          double thetaf = x_dev[xidxf[i]], Vmf = x_dev[xidxf[i] + 1];
          double thetat = x_dev[xidxt[i]], Vmt = x_dev[xidxt[i] + 1];
          double thetaft = thetaf - thetat;
          double thetatf = thetat - thetaf;
          int ifrom(busf_idx[i]), ito(bust_idx[i]);

          // This confusing and could probably be done in a clearer
          // way. Two indexing schemes are needed. Some line
          // contributions (from/from, to/to) need to be added to the
          // original bus contribution entries -- get those
          // from the bus index list.  A separate indexing system is for
          // the extra (to/from, from/to) entries.

          // for reference, these are computed in the same order as
          // OPFLOWComputeDenseEqualityConstraintJacobian_PBPOLRAJAHIOP()

          // from bus real entries
          
          /* dPf_dthetaf */
          RAJA::atomicAdd<exago_raja_atomic>
            (&(MJacS_dev[b_jacsp_idx[ifrom]]),
             Vmf * Vmt * (-Gft[i] * sin(thetaft) + Bft[i] * cos(thetaft)));
          /*dPf_dVmf */
          RAJA::atomicAdd<exago_raja_atomic>
            (&(MJacS_dev[b_jacsp_idx[ifrom] + 1]),
             2 * Gff[i] * Vmf + Vmt * (Gft[i] * cos(thetaft) + Bft[i] * sin(thetaft)));
          /*dPf_dthetat */
          MJacS_dev[jacf_idx[i] + 0] =
            Vmf * Vmt * (Gft[i] * sin(thetaft) - Bft[i] * cos(thetaft));
          /* dPf_dVmt */
          MJacS_dev[jacf_idx[i] + 1] =
            Vmf * (Gft[i] * cos(thetaft) + Bft[i] * sin(thetaft));

          // from bus reactive entries
          
          /* dQf_dthetaf */
          RAJA::atomicAdd<exago_raja_atomic>
            (&(MJacS_dev[b_jacsq_idx[ifrom]]),
             Vmf * Vmt * (Bft[i] * sin(thetaft) + Gft[i] * cos(thetaft)));
          /* dQf_dVmf */
          RAJA::atomicAdd<exago_raja_atomic>
            (&(MJacS_dev[b_jacsq_idx[ifrom] + 1]),
             -2 * Bff[i] * Vmf +
             Vmt * (-Bft[i] * cos(thetaft) + Gft[i] * sin(thetaft)));
          /* dQf_dthetat */
          MJacS_dev[jacf_idx[i] + 2] =
            Vmf * Vmt * (-Bft[i] * sin(thetaft) - Gft[i] * cos(thetaft));
          /* dQf_dVmt */
          MJacS_dev[jacf_idx[i] + 3] =
            Vmf * (-Bft[i] * cos(thetaft) + Gft[i] * sin(thetaft));

          // to bus real entries
          
          /* dPt_dthetat */
          RAJA::atomicAdd<exago_raja_atomic>
             (&(MJacS_dev[b_jacsp_idx[ito]]), 
              Vmt * Vmf * (-Gtf[i] * sin(thetatf) + Btf[i] * cos(thetatf)));
          /* dPt_dVmt */
          RAJA::atomicAdd<exago_raja_atomic> 
            (&(MJacS_dev[b_jacsp_idx[ito] + 1]), 2 * Gtt[i] * Vmt +
               Vmf * (Gtf[i] * cos(thetatf) + Btf[i] * sin(thetatf)));
          /* dPt_dthetaf */  
          MJacS_dev[jact_idx[i] + 0] =
            Vmt * Vmf * (Gtf[i] * sin(thetatf) - Btf[i] * cos(thetatf));
          /* dPt_dVmf */
          MJacS_dev[jact_idx[i] + 1] = 
            Vmt * (Gtf[i] * cos(thetatf) + Btf[i] * sin(thetatf));

          // to bus reactive entries
          
          /* dQt_dthetat */
          RAJA::atomicAdd<exago_raja_atomic> 
             (&(MJacS_dev[b_jacsq_idx[ito]]),
              Vmt * Vmf * (Btf[i] * sin(thetatf) + Gtf[i] * cos(thetatf)));
          /* dQt_dVmt */
          RAJA::atomicAdd<exago_raja_atomic> 
             (&(MJacS_dev[b_jacsq_idx[ito] + 1]), -2 * Btt[i] * Vmt +
              Vmf * (-Btf[i] * cos(thetatf) + Gtf[i] * sin(thetatf)));
          /* dQt_dthetaf */
          MJacS_dev[jact_idx[i] + 2] =
            Vmt * Vmf * (-Btf[i] * sin(thetatf) - Gtf[i] * cos(thetatf));
          /* dQt_dVmf */
          MJacS_dev[jact_idx[i] + 3] =
            Vmt * (-Btf[i] * cos(thetatf) + Gtf[i] * sin(thetatf));
        });
    
    int *ipermout = (int *)d_allocator_.allocate(opflow->nnz_eqjacsp*sizeof(int));
    resmgr.copy(ipermout, iperm);
    
    RAJA::stable_sort_pairs<exago_raja_exec>
      (RAJA::make_span(ipermout, opflow->nnz_eqjacsp),
       RAJA::make_span(MJacS_dev, opflow->nnz_eqjacsp),
       RAJA::operators::less<int>{});

    d_allocator_.deallocate(ipermout);
    
    if (debugmsg)
      PrintTriplets("Equality Constraint Jacobian (GPU):",
                    opflow->nnz_eqjacsp, iperm, iJacS_dev, jJacS_dev, MJacS_dev);
    
    if (oldhostway) {
    ierr = VecGetArray(opflow->X, &x);
    CHKERRQ(ierr);

    // Copy from device to host
    registerWith(x, opflow->nx, resmgr, h_allocator_);
    resmgr.copy((double *)x, (double *)x_dev);

    ierr = VecRestoreArray(opflow->X, &x);
    CHKERRQ(ierr);

    /* Compute equality constraint jacobian */
    ierr = (*opflow->modelops.computeequalityconstraintjacobian)(
        opflow, opflow->X, opflow->Jac_Ge);
    CHKERRQ(ierr);

    ierr = MatGetSize(opflow->Jac_Ge, &nrow, &ncol);
    CHKERRQ(ierr);

    values = pbpolrajahiopsparse->val_jaceq;

    /* Copy over values */
    for (i = 0; i < nrow; i++) {
      ierr = MatGetRow(opflow->Jac_Ge, i, &nvals, &cols, &vals);
      CHKERRQ(ierr);
      for (j = 0; j < nvals; j++) {
        values[j] = vals[j];
      }
      values += nvals;
      ierr = MatRestoreRow(opflow->Jac_Ge, i, &nvals, &cols, &vals);
      CHKERRQ(ierr);
    }

    // Copy over val_ineq to device
    resmgr.copy(MJacS_dev, pbpolrajahiopsparse->val_jaceq);

    if (debugmsg)
      PrintTriplets("Equality Constraint Jacobian:",
                    opflow->nnz_eqjacsp, NULL, iJacS_dev, jJacS_dev, MJacS_dev);
    }
  }

  ierr = PetscLogEventEnd(opflow->eqconsjaclogger, 0, 0, 0, 0);
  CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWComputeSparseHessian_PBPOLRAJAHIOPSPARSE(
    OPFLOW opflow, const double *x_dev, const double *lambda_dev, int *iHSS_dev,
    int *jHSS_dev, double *MHSS_dev) {
  PbpolModelRajaHiop *pbpolrajahiopsparse =
      reinterpret_cast<PbpolModelRajaHiop *>(opflow->model);
  GENParamsRajaHiop *genparams = &pbpolrajahiopsparse->genparams;
  LOADParamsRajaHiop *loadparams = &pbpolrajahiopsparse->loadparams;
  BUSParamsRajaHiop *busparams = &pbpolrajahiopsparse->busparams;
  LINEParamsRajaHiop *lineparams = &pbpolrajahiopsparse->lineparams;

  PetscErrorCode ierr;
  PetscInt *iRow, *jCol;
  PetscScalar *x, *values, *lambda;
  PetscInt nrow, ncol;
  PetscInt nvals;
  const PetscInt *cols;
  const PetscScalar *vals;
  PetscInt i, j;
  PetscInt ctr = 0;
  auto &resmgr = umpire::ResourceManager::getInstance();

  umpire::Allocator h_allocator_ = resmgr.getAllocator("HOST");
  umpire::Allocator d_allocator_ = resmgr.getAllocator("DEVICE");

  PetscFunctionBegin;

  if (iHSS_dev != NULL && jHSS_dev != NULL) {

    if (debugmsg)
      std::cout << "Official Hessian nonzero count: "
                << opflow->nnz_hesssp << std::endl;
    
    resmgr.memset(iHSS_dev, 0, opflow->nnz_hesssp*sizeof(int));
    resmgr.memset(jHSS_dev, 0, opflow->nnz_hesssp*sizeof(int));

    // Bus contributions
    
    int *b_xidx = busparams->xidx_dev_;
    int *b_hesssp_idx = busparams->hesssp_idx_dev_;
    
    RAJA::forall<exago_raja_exec>(
        RAJA::RangeSegment(0, busparams->nbus), RAJA_LAMBDA(RAJA::Index_type i) {
          int off(0);
          iHSS_dev[b_hesssp_idx[i] + off] = b_xidx[i];
          jHSS_dev[b_hesssp_idx[i] + off] = b_xidx[i];
          off++;
          
          iHSS_dev[b_hesssp_idx[i] + off] = b_xidx[i];
          jHSS_dev[b_hesssp_idx[i] + off] = b_xidx[i] + 1;
          off++;

          // upper triangular only
          // iHSS_dev[b_hesssp_idx[i] + off] = b_xidx[i] + 1;
          // jHSS_dev[b_hesssp_idx[i] + off] = b_xidx[i];
          // off++;

          iHSS_dev[b_hesssp_idx[i] + off] = b_xidx[i] + 1;
          jHSS_dev[b_hesssp_idx[i] + off] = b_xidx[i] + 1;
          off++;
        });

    if (opflow->include_powerimbalance_variables) {
      int *b_xidxpimb = busparams->xidxpimb_dev_;
      RAJA::forall<exago_raja_exec>(
          RAJA::RangeSegment(0, busparams->nbus), RAJA_LAMBDA(RAJA::Index_type i) {
            int off(2);
            
            iHSS_dev[b_hesssp_idx[i] + off] = b_xidxpimb[i];
            jHSS_dev[b_hesssp_idx[i] + off] = b_xidxpimb[i];
            off++;

            iHSS_dev[b_hesssp_idx[i] + off] = b_xidxpimb[i] + 1;
            jHSS_dev[b_hesssp_idx[i] + off] = b_xidxpimb[i] + 1;
          });
    }

    /* Generator contributions for row,col numbers */
    int *g_xidx = genparams->xidx_dev_;
    int *g_hesssp_idx = genparams->hesssp_idx_dev_;
    RAJA::forall<exago_raja_exec>(
        RAJA::RangeSegment(0, genparams->ngenON),
        RAJA_LAMBDA(RAJA::Index_type i) {
          iHSS_dev[g_hesssp_idx[i]] = g_xidx[i];
          jHSS_dev[g_hesssp_idx[i]] = g_xidx[i];
          iHSS_dev[g_hesssp_idx[i] + 1] = g_xidx[i] + 1;
          jHSS_dev[g_hesssp_idx[i] + 1] = g_xidx[i] + 1;
        });

    int *xidxf = lineparams->xidxf_dev_;
    int *xidxt = lineparams->xidxt_dev_;
    int *ln_hessp_idx = lineparams->hesssp_idx_dev_;
    int *linelimidx = lineparams->linelimidx_dev_;

    RAJA::forall<exago_raja_exec>(
        RAJA::RangeSegment(0, lineparams->nlinelim),
        RAJA_LAMBDA(RAJA::Index_type i) {
          int off(0);

          // from-bus diagonal entries already defined

          // iHSS_dev[ln_hessp_idx[i] + off] = xidxf[i]; 
          // jHSS_dev[ln_hessp_idx[i] + off] = xidxf[i];
          // off++;
          
          // iHSS_dev[ln_hessp_idx[i] + off] = xidxf[i];
          // jHSS_dev[ln_hessp_idx[i] + off] = xidxf[i] + 1;
          // off++;

          // iHSS_dev[ln_hessp_idx[i] + off] = xidxf[i] + 1;
          // jHSS_dev[ln_hessp_idx[i] + off] = xidxf[i];
          // off++;
          
          // iHSS_dev[ln_hessp_idx[i] + off] = xidxf[i] + 1;
          // jHSS_dev[ln_hessp_idx[i] + off] = xidxf[i] + 1;
          // off++;

          // from-bus off diagonal entries only there if in upper part
          if (xidxt[i] > xidxf[i]) {
            
            iHSS_dev[ln_hessp_idx[i] + off] = xidxf[i];
            jHSS_dev[ln_hessp_idx[i] + off] = xidxt[i];
            off++;
          
            iHSS_dev[ln_hessp_idx[i] + off] = xidxf[i];
            jHSS_dev[ln_hessp_idx[i] + off] = xidxt[i] + 1;
            off++;

            iHSS_dev[ln_hessp_idx[i] + off] = xidxf[i] + 1;
            jHSS_dev[ln_hessp_idx[i] + off] = xidxt[i];
            off++;
          
            iHSS_dev[ln_hessp_idx[i] + off] = xidxf[i] + 1;
            jHSS_dev[ln_hessp_idx[i] + off] = xidxt[i] + 1;
            off++;
          }

          // to-bus diagonal entries already defined

          // iHSS_dev[ln_hessp_idx[i] + off] = xidxt[i];
          // jHSS_dev[ln_hessp_idx[i] + off] = xidxt[i];
          // off++;
          
          // iHSS_dev[ln_hessp_idx[i] + off] = xidxt[i];
          // jHSS_dev[ln_hessp_idx[i] + off] = xidxt[i] + 1;
          // off++;

          // iHSS_dev[ln_hessp_idx[i] + off] = xidxt[i] + 1;
          // jHSS_dev[ln_hessp_idx[i] + off] = xidxt[i];
          // off++;
          
          // iHSS_dev[ln_hessp_idx[i] + off] = xidxt[i] + 1;
          // jHSS_dev[ln_hessp_idx[i] + off] = xidxt[i] + 1;
          // off++;

          // to-bus off diagonal entries only there if in upper part
          if (xidxf[i] > xidxt[i]) {

            iHSS_dev[ln_hessp_idx[i] + off] = xidxt[i];
            jHSS_dev[ln_hessp_idx[i] + off] = xidxf[i];
            off++;
          
            iHSS_dev[ln_hessp_idx[i] + off] = xidxt[i];
            jHSS_dev[ln_hessp_idx[i] + off] = xidxf[i] + 1;
            off++;

            iHSS_dev[ln_hessp_idx[i] + off] = xidxt[i] + 1;
            jHSS_dev[ln_hessp_idx[i] + off] = xidxf[i];
            off++;
          
            iHSS_dev[ln_hessp_idx[i] + off] = xidxt[i] + 1;
            jHSS_dev[ln_hessp_idx[i] + off] = xidxf[i] + 1;
            off++;
          }
        });
    
    /* Loadloss contributions - two contributions*/
    if (opflow->include_loadloss_variables) {
      int *l_xidx = loadparams->xidx_dev_;
      int *l_hesssp_idx = loadparams->hesssp_idx_dev_;
      RAJA::forall<exago_raja_exec>(
          RAJA::RangeSegment(0, loadparams->nload),
          RAJA_LAMBDA(RAJA::Index_type i) {
            iHSS_dev[l_hesssp_idx[i]] = l_xidx[i];
            jHSS_dev[l_hesssp_idx[i]] = l_xidx[i];
            iHSS_dev[l_hesssp_idx[i] + 1] = l_xidx[i] + 1;
            jHSS_dev[l_hesssp_idx[i] + 1] = l_xidx[i] + 1;
          });
    }

    if (pbpolrajahiopsparse->idx_hess_dev_ == NULL) {
      pbpolrajahiopsparse->idx_hess_dev_ =
        (int *) d_allocator_.allocate(opflow->nnz_hesssp * sizeof(int));
    }

    SortIndexes(opflow->nnz_hesssp, iHSS_dev, jHSS_dev,
                pbpolrajahiopsparse->idx_hess_dev_);
    
    if (debugmsg)
      PrintTriplets("Hessian Indexes (GPU):",
                    opflow->nnz_hesssp, pbpolrajahiopsparse->idx_hess_dev_,
                    iHSS_dev, jHSS_dev, NULL);
    
    // Create arrays on host to store i,j, and val arrays

    if (pbpolrajahiopsparse->i_hess == NULL) { 
      pbpolrajahiopsparse->i_hess =
        (int *)(h_allocator_.allocate(opflow->nnz_hesssp * sizeof(int)));
      pbpolrajahiopsparse->j_hess =
        (int *)(h_allocator_.allocate(opflow->nnz_hesssp * sizeof(int)));
      pbpolrajahiopsparse->val_hess =
        (double *)(h_allocator_.allocate(opflow->nnz_hesssp * sizeof(double)));
    }
    
    iRow = pbpolrajahiopsparse->i_hess;
    jCol = pbpolrajahiopsparse->j_hess;

    ierr = (*opflow->modelops.computehessian)(
        opflow, opflow->X, opflow->Lambdae, opflow->Lambdai, opflow->Hes);
    CHKERRQ(ierr);
    ierr = MatGetSize(opflow->Hes, &nrow, &ncol);
    CHKERRQ(ierr);

    if (debugmsg)
      std::cout << "Official Hessian Size: "
                << nrow << " rows x " << ncol << " cols"
                << "(should be" << opflow->Nx << " x " << opflow->Nx << ")"
                << std::endl;

    /* Copy over locations to triplet format */
    /* Note that HIOP requires a upper triangular Hessian as oppposed
       to IPOPT which requires a lower triangular Hessian
    */
    for (i = 0; i < nrow; i++) {
      ierr = MatGetRow(opflow->Hes, i, &nvals, &cols, &vals);
      CHKERRQ(ierr);
      ctr = 0;
      for (j = 0; j < nvals; j++) {
        if (cols[j] >= i) { /* upper triangle */
          /* save as upper triangle locations */
          iRow[ctr] = i;
          jCol[ctr] = cols[j];
          ctr++;
        }
      }
      iRow += ctr;
      jCol += ctr;
      ierr = MatRestoreRow(opflow->Hes, i, &nvals, &cols, &vals);
      CHKERRQ(ierr);
    }

    // Copy over i_hess and j_hess arrays to device
    resmgr.copy(iHSS_dev, pbpolrajahiopsparse->i_hess);
    resmgr.copy(jHSS_dev, pbpolrajahiopsparse->j_hess);

    if (debugmsg) {
      PrintTriplets("Hessian Indexes:",
                    opflow->nnz_hesssp, NULL, iHSS_dev, jHSS_dev, NULL);
    }
    
  } else {

    resmgr.memset(MHSS_dev, 0, opflow->nnz_hesssp*sizeof(double));


    // Bus contributions
    
    int *b_hesssp_idx = busparams->hesssp_idx_dev_;
    int *b_gidx = busparams->gidx_dev_;
    int *ispvpq = busparams->ispvpq_dev_;
    double *gl = busparams->gl_dev_;
    double *bl = busparams->bl_dev_;

    RAJA::forall<exago_raja_exec>(
        RAJA::RangeSegment(0, busparams->nbus), RAJA_LAMBDA(RAJA::Index_type i) {
          // int row, col;
          double val;
          // row = b_xidx[i] + 1 - nxsparse;
          // col = row;
          val = ispvpq[i] * (lambda_dev[b_gidx[i]] * 2 * gl[i] +
                             lambda_dev[b_gidx[i] + 1] * (-2 * bl[i]));
          RAJA::atomicAdd<exago_raja_atomic>
            (&MHSS_dev[b_hesssp_idx[i] + 2], val);
      });

    if (opflow->objectivetype == MIN_GEN_COST) {
      int *hesssp_idx = genparams->hesssp_idx_dev_;
      double *cost_alpha = genparams->cost_alpha_dev_;
      double obj_factor = opflow->obj_factor;
      int isobj_gencost = opflow->obj_gencost;
      double MVAbase = opflow->ps->MVAbase;
      double weight = opflow->weight;

      RAJA::forall<exago_raja_exec>(
          RAJA::RangeSegment(0, genparams->ngenON),
          RAJA_LAMBDA(RAJA::Index_type i) {
            MHSS_dev[hesssp_idx[i]] =
              weight * isobj_gencost * obj_factor *
              2.0 * cost_alpha[i] * MVAbase * MVAbase;
            MHSS_dev[hesssp_idx[i] + 1] = 0.0;
          });
    } else if (opflow->objectivetype == NO_OBJ) {
      int *hesssp_idx = genparams->hesssp_idx_dev_;
      RAJA::forall<exago_raja_exec>(
          RAJA::RangeSegment(0, genparams->ngenON),
          RAJA_LAMBDA(RAJA::Index_type i) {
            MHSS_dev[hesssp_idx[i]] = 0.0;
            MHSS_dev[hesssp_idx[i] + 1] = 0.0;
          });
    }

    // Line contributions

    double *Gff_arr = lineparams->Gff_dev_;
    double *Gtt_arr = lineparams->Gtt_dev_;
    double *Gft_arr = lineparams->Gft_dev_;
    double *Gtf_arr = lineparams->Gtf_dev_;

    double *Bff_arr = lineparams->Bff_dev_;
    double *Btt_arr = lineparams->Btt_dev_;
    double *Bft_arr = lineparams->Bft_dev_;
    double *Btf_arr = lineparams->Btf_dev_;

    int *busf_idx = lineparams->busf_idx_dev_;
    int *bust_idx = lineparams->bust_idx_dev_;
    int *xidxf = lineparams->xidxf_dev_;
    int *xidxt = lineparams->xidxt_dev_;
    int *geqidxf = lineparams->geqidxf_dev_;
    int *geqidxt = lineparams->geqidxt_dev_;
    int *ln_hessp_idx = lineparams->hesssp_idx_dev_;

    RAJA::forall<exago_raja_exec>(
        RAJA::RangeSegment(0, lineparams->nlineON),
        RAJA_LAMBDA(RAJA::Index_type i) {
          int gloc;
          // int row[2], col[4];
          double val[8];
          double Gff, Bff, Gft, Bft, Gtf, Btf, Gtt, Btt;
          int ibusf(busf_idx[i]), ibust(bust_idx[i]);
          int fbusidx(ibusf), tbusidx(ibust);
          Gff = Gff_arr[i];
          Bff = Bff_arr[i];
          Gft = Gft_arr[i];
          Bft = Bft_arr[i];
          Gtf = Gtf_arr[i];
          Btf = Btf_arr[i];
          Gtt = Gtt_arr[i];
          Btt = Btt_arr[i];

          double thetaf = x_dev[xidxf[i]], Vmf = x_dev[xidxf[i] + 1];
          double thetat = x_dev[xidxt[i]], Vmt = x_dev[xidxt[i] + 1];
          double thetaft = thetaf - thetat;
          double thetatf = thetat - thetaf;

          double dPf_dthetaf_dthetaf, dPf_dthetaf_dVmf, dPf_dthetaf_dthetat,
            dPf_dthetaf_dVmt;
          double dPf_dVmf_dthetaf, dPf_dVmf_dVmf, dPf_dVmf_dthetat, dPf_dVmf_dVmt;
          double dPf_dthetat_dthetaf, dPf_dthetat_dVmf, dPf_dthetat_dthetat,
            dPf_dthetat_dVmt;
          double dPf_dVmt_dthetaf, dPf_dVmt_dVmf, dPf_dVmt_dthetat, dPf_dVmt_dVmt;

          /* dPf_dthetaf = Vmf*Vmt*(-Gft*sin(thetaft) + Bft*cos(thetaft)); */
          dPf_dthetaf_dthetaf =
            -Vmf * Vmt * (Gft * cos(thetaft) + Bft * sin(thetaft));
          dPf_dthetaf_dVmf = Vmt * (-Gft * sin(thetaft) + Bft * cos(thetaft));
          dPf_dthetaf_dthetat =
            Vmf * Vmt * (Gft * cos(thetaft) + Bft * sin(thetaft));
          dPf_dthetaf_dVmt = Vmf * (-Gft * sin(thetaft) + Bft * cos(thetaft));

          /* dPf_Vmf  = 2*Gff*Vmf + Vmt*(Gft*cos(thetaft) + Bft*sin(thetaft)); */
          dPf_dVmf_dthetaf = Vmt * (-Gft * sin(thetaft) + Bft * cos(thetaft));
          dPf_dVmf_dVmf = 2 * Gff;
          dPf_dVmf_dthetat = Vmt * (Gft * sin(thetaft) - Bft * cos(thetaft));
          dPf_dVmf_dVmt = (Gft * cos(thetaft) + Bft * sin(thetaft));

          /* dPf_dthetat = Vmf*Vmt*(Gft*sin(thetaft) - Bft*cos(thetaft)); */
          dPf_dthetat_dthetaf =
            Vmf * Vmt * (Gft * cos(thetaft) + Bft * sin(thetaft));
          dPf_dthetat_dVmf = Vmt * (Gft * sin(thetaft) - Bft * cos(thetaft));
          dPf_dthetat_dthetat =
            Vmf * Vmt * (-Gft * cos(thetaft) - Bft * sin(thetaft));
          dPf_dthetat_dVmt = Vmf * (Gft * sin(thetaft) - Bft * cos(thetaft));

          /* dPf_dVmt = Vmf*(Gft*cos(thetaft) + Bft*sin(thetaft)); */
          dPf_dVmt_dthetaf = Vmf * (-Gft * sin(thetaft) + Bft * cos(thetaft));
          dPf_dVmt_dVmf = (Gft * cos(thetaft) + Bft * sin(thetaft));
          dPf_dVmt_dthetat = Vmf * (Gft * sin(thetaft) - Bft * cos(thetaft));
          dPf_dVmt_dVmt = 0.0;

          double dQf_dthetaf_dthetaf, dQf_dthetaf_dVmf, dQf_dthetaf_dthetat,
            dQf_dthetaf_dVmt;
          double dQf_dVmf_dthetaf, dQf_dVmf_dVmf, dQf_dVmf_dthetat, dQf_dVmf_dVmt;
          double dQf_dthetat_dthetaf, dQf_dthetat_dVmf, dQf_dthetat_dthetat,
            dQf_dthetat_dVmt;
          double dQf_dVmt_dthetaf, dQf_dVmt_dVmf, dQf_dVmt_dthetat, dQf_dVmt_dVmt;

          /* dQf_dthetaf = Vmf*Vmt*(Bft*sin(thetaft) + Gft*cos(thetaft)); */
          dQf_dthetaf_dthetaf =
            Vmf * Vmt * (Bft * cos(thetaft) - Gft * sin(thetaft));
          dQf_dthetaf_dVmf = Vmt * (Bft * sin(thetaft) + Gft * cos(thetaft));
          dQf_dthetaf_dthetat =
            Vmf * Vmt * (-Bft * cos(thetaft) + Gft * sin(thetaft));
          dQf_dthetaf_dVmt = Vmf * (Bft * sin(thetaft) + Gft * cos(thetaft));

          /* dQf_dVmf = -2*Bff*Vmf + Vmt*(-Bft*cos(thetaft) + Gft*sin(thetaft));
           */
          dQf_dVmf_dthetaf = Vmt * (Bft * sin(thetaft) + Gft * cos(thetaft));
          dQf_dVmf_dVmf = -2 * Bff;
          dQf_dVmf_dthetat = Vmt * (-Bft * sin(thetaft) - Gft * cos(thetaft));
          dQf_dVmf_dVmt = (-Bft * cos(thetaft) + Gft * sin(thetaft));

          /* dQf_dthetat = Vmf*Vmt*(-Bft*sin(thetaft) - Gft*cos(thetaft)); */
          dQf_dthetat_dthetaf =
            Vmf * Vmt * (-Bft * cos(thetaft) + Gft * sin(thetaft));
          dQf_dthetat_dVmf = Vmt * (-Bft * sin(thetaft) - Gft * cos(thetaft));
          dQf_dthetat_dthetat =
            Vmf * Vmt * (Bft * cos(thetaft) - Gft * sin(thetaft));
          dQf_dthetat_dVmt = Vmf * (-Bft * sin(thetaft) - Gft * cos(thetaft));

          /* dQf_dVmt = Vmf*(-Bft*cos(thetaft) + Gft*sin(thetaft)); */
          dQf_dVmt_dthetaf = Vmf * (Bft * sin(thetaft) + Gft * cos(thetaft));
          dQf_dVmt_dVmf = (-Bft * cos(thetaft) + Gft * sin(thetaft));
          dQf_dVmt_dthetat = Vmf * (-Bft * sin(thetaft) - Gft * cos(thetaft));
          dQf_dVmt_dVmt = 0.0;

          // row[0] = xidxf[i] - nxsparse;
          // row[1] = xidxf[i] + 1 - nxsparse;
          // col[0] = xidxf[i] - nxsparse;
          // col[1] = xidxf[i] + 1 - nxsparse;
          // col[2] = xidxt[i] - nxsparse;
          // col[3] = xidxt[i] + 1 - nxsparse;

          gloc = geqidxf[i];

          val[0] = val[1] = val[2] = val[3] = val[4] = val[5] = val[6] = val[7] =
            0.0;

          val[0] = lambda_dev[gloc] * dPf_dthetaf_dthetaf +
            lambda_dev[gloc + 1] * dQf_dthetaf_dthetaf;
          val[1] = lambda_dev[gloc] * dPf_dthetaf_dVmf +
            lambda_dev[gloc + 1] * dQf_dthetaf_dVmf;
          val[2] = lambda_dev[gloc] * dPf_dthetaf_dthetat +
            lambda_dev[gloc + 1] * dQf_dthetaf_dthetat;
          val[3] = lambda_dev[gloc] * dPf_dthetaf_dVmt +
            lambda_dev[gloc + 1] * dQf_dthetaf_dVmt;

          val[4] = lambda_dev[gloc] * dPf_dVmf_dthetaf +
            lambda_dev[gloc + 1] * dQf_dVmf_dthetaf;
          val[5] = lambda_dev[gloc] * dPf_dVmf_dVmf +
            lambda_dev[gloc + 1] * dQf_dVmf_dVmf;
          val[6] = lambda_dev[gloc] * dPf_dVmf_dthetat +
            lambda_dev[gloc + 1] * dQf_dVmf_dthetat;
          val[7] = lambda_dev[gloc] * dPf_dVmf_dVmt +
            lambda_dev[gloc + 1] * dQf_dVmf_dVmt;

          // Remember central bus locations were reserved and indexed
          // by bus (from-from)

          RAJA::atomicAdd<exago_raja_atomic>
            (&MHSS_dev[b_hesssp_idx[fbusidx] + 0], val[0]);
          RAJA::atomicAdd<exago_raja_atomic>
            (&MHSS_dev[b_hesssp_idx[fbusidx] + 1], val[1]);
          // not in upper triangle
          // RAJA::atomicAdd<exago_raja_atomic>
          //   (&MHSS_dev[b_hesssp_idx[fbusidx] + 2], val[4]);
          RAJA::atomicAdd<exago_raja_atomic>
            (&MHSS_dev[b_hesssp_idx[fbusidx] + 2], val[5]);

          // Off-center entries (from-to bus) were reserved and
          // indexed by line only if in upper triangle

          if (xidxt[i] > xidxf[i]) {
            RAJA::atomicAdd<exago_raja_atomic>
              (&MHSS_dev[ln_hessp_idx[i] + 0], val[2]);
            RAJA::atomicAdd<exago_raja_atomic>
              (&MHSS_dev[ln_hessp_idx[i] + 1], val[3]);
            RAJA::atomicAdd<exago_raja_atomic>
              (&MHSS_dev[ln_hessp_idx[i] + 2], val[6]);
            RAJA::atomicAdd<exago_raja_atomic>
              (&MHSS_dev[ln_hessp_idx[i] + 3], val[7]);
          }
            
          // row[0] = xidxt[i] - nxsparse;
          // row[1] = xidxt[i] + 1 - nxsparse;

          // col[0] = xidxf[i] - nxsparse;
          // col[1] = xidxf[i] + 1 - nxsparse;
          // col[2] = xidxt[i] - nxsparse;
          // col[3] = xidxt[i] + 1 - nxsparse;

          val[0] = val[1] = val[2] = val[3] = val[4] = val[5] = val[6] = val[7] =
            0.0;

          val[0] = lambda_dev[gloc] * dPf_dthetat_dthetaf +
            lambda_dev[gloc + 1] * dQf_dthetat_dthetaf;
          val[1] = lambda_dev[gloc] * dPf_dthetat_dVmf +
            lambda_dev[gloc + 1] * dQf_dthetat_dVmf;
          val[2] = lambda_dev[gloc] * dPf_dthetat_dthetat +
            lambda_dev[gloc + 1] * dQf_dthetat_dthetat;
          val[3] = lambda_dev[gloc] * dPf_dthetat_dVmt +
            lambda_dev[gloc + 1] * dQf_dthetat_dVmt;

          val[4] = lambda_dev[gloc] * dPf_dVmt_dthetaf +
            lambda_dev[gloc + 1] * dQf_dVmt_dthetaf;
          val[5] = lambda_dev[gloc] * dPf_dVmt_dVmf +
            lambda_dev[gloc + 1] * dQf_dVmt_dVmf;
          val[6] = lambda_dev[gloc] * dPf_dVmt_dthetat +
            lambda_dev[gloc + 1] * dQf_dVmt_dthetat;
          val[7] = lambda_dev[gloc] * dPf_dVmt_dVmt +
            lambda_dev[gloc + 1] * dQf_dVmt_dVmt;

          
          // Remember central bus locations were reserved and indexed
          // by bus (to-to)
          
          RAJA::atomicAdd<exago_raja_atomic>
            (&MHSS_dev[b_hesssp_idx[tbusidx] + 0], val[2]);
          RAJA::atomicAdd<exago_raja_atomic>
            (&MHSS_dev[b_hesssp_idx[tbusidx] + 1], val[3]);
          // not in upper triangle
          // RAJA::atomicAdd<exago_raja_atomic>
          //   (&MHSS_dev[b_hesssp_idx[tbusidx] + 2], val[6]);
          RAJA::atomicAdd<exago_raja_atomic>
            (&MHSS_dev[b_hesssp_idx[tbusidx] + 2], val[7]);

          // Off-center entries (to-from bus) were reserved and
          // indexed by line only if in upper triangle

          if (xidxf[i] > xidxt[i]) {
            RAJA::atomicAdd<exago_raja_atomic>
              (&MHSS_dev[ln_hessp_idx[i] + 0], val[0]);
            RAJA::atomicAdd<exago_raja_atomic>
              (&MHSS_dev[ln_hessp_idx[i] + 1], val[1]);
            RAJA::atomicAdd<exago_raja_atomic>
              (&MHSS_dev[ln_hessp_idx[i] + 2], val[4]);
            RAJA::atomicAdd<exago_raja_atomic>
              (&MHSS_dev[ln_hessp_idx[i] + 3], val[5]);
          }

          //    ierr = MatSetValues(H,2,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);

          double dPt_dthetat_dthetat, dPt_dthetat_dVmt, dPt_dthetat_dthetaf,
            dPt_dthetat_dVmf;
          double dPt_dVmt_dthetat, dPt_dVmt_dVmt, dPt_dVmt_dthetaf, dPt_dVmt_dVmf;
          double dPt_dthetaf_dthetat, dPt_dthetaf_dVmt, dPt_dthetaf_dthetaf,
            dPt_dthetaf_dVmf;
          double dPt_dVmf_dthetat, dPt_dVmf_dVmt, dPt_dVmf_dthetaf, dPt_dVmf_dVmf;

          /* dPt_dthetat = Vmf*Vmt*(-Gtf*sin(thetatf) + Btf*cos(thetatf)); */
          dPt_dthetat_dthetat =
            Vmf * Vmt * (-Gtf * cos(thetatf) - Btf * sin(thetatf));
          dPt_dthetat_dVmt = Vmf * (-Gtf * sin(thetatf) + Btf * cos(thetatf));
          dPt_dthetat_dthetaf =
            Vmf * Vmt * (Gtf * cos(thetatf) + Btf * sin(thetatf));
          dPt_dthetat_dVmf = Vmt * (-Gtf * sin(thetatf) + Btf * cos(thetatf));

          /* dPt_Vmt  = 2*Gtt*Vmt + Vmf*(Gtf*cos(thetatf) + Btf*sin(thetatf)); */
          dPt_dVmt_dthetat = Vmf * (-Gtf * sin(thetatf) + Bft * cos(thetatf));
          dPt_dVmt_dVmt = 2 * Gtt;
          dPt_dVmt_dthetaf = Vmf * (Gtf * sin(thetatf) - Btf * cos(thetatf));
          dPt_dVmt_dVmf = (Gtf * cos(thetatf) + Btf * sin(thetatf));

          /* dPt_dthetaf = Vmf*Vmt*(Gtf*sin(thetatf) - Btf*cos(thetatf)); */
          dPt_dthetaf_dthetat =
            Vmf * Vmt * (Gtf * cos(thetatf) + Btf * sin(thetatf));
          dPt_dthetaf_dVmt = Vmf * (Gtf * sin(thetatf) - Btf * cos(thetatf));
          dPt_dthetaf_dthetaf =
            Vmf * Vmt * (-Gtf * cos(thetatf) - Btf * sin(thetatf));
          dPt_dthetaf_dVmf = Vmt * (Gtf * sin(thetatf) - Btf * cos(thetatf));

          /* dPt_dVmf = Vmt*(Gtf*cos(thetatf) + Btf*sin(thetatf)); */
          dPt_dVmf_dthetat = Vmt * (-Gtf * sin(thetatf) + Btf * cos(thetatf));
          dPt_dVmf_dVmt = (Gtf * cos(thetatf) + Btf * sin(thetatf));
          dPt_dVmf_dthetaf = Vmt * (Gtf * sin(thetatf) - Btf * cos(thetatf));
          dPt_dVmf_dVmf = 0.0;

          double dQt_dthetaf_dthetaf, dQt_dthetaf_dVmf, dQt_dthetaf_dthetat,
            dQt_dthetaf_dVmt;
          double dQt_dVmf_dthetaf, dQt_dVmf_dVmf, dQt_dVmf_dthetat, dQt_dVmf_dVmt;
          double dQt_dthetat_dthetaf, dQt_dthetat_dVmf, dQt_dthetat_dthetat,
            dQt_dthetat_dVmt;
          double dQt_dVmt_dthetaf, dQt_dVmt_dVmf, dQt_dVmt_dthetat, dQt_dVmt_dVmt;

          /* dQt_dthetat = Vmf*Vmt*(Btf*sin(thetatf) + Gtf*cos(thetatf)); */
          dQt_dthetat_dthetat =
            Vmf * Vmt * (Btf * cos(thetatf) - Gtf * sin(thetatf));
          dQt_dthetat_dVmt = Vmf * (Btf * sin(thetatf) + Gtf * cos(thetatf));
          dQt_dthetat_dthetaf =
            Vmf * Vmt * (-Btf * cos(thetatf) + Gtf * sin(thetatf));
          dQt_dthetat_dVmf = Vmf * (Btf * sin(thetatf) + Gtf * cos(thetatf));

          /* dQt_dVmt = -2*Btt*Vmt + Vmf*(-Btf*cos(thetatf) + Gtf*sin(thetatf));
           */
          dQt_dVmt_dthetat = Vmf * (Btf * sin(thetatf) + Gtf * cos(thetatf));
          dQt_dVmt_dVmt = -2 * Btt;
          dQt_dVmt_dthetaf = Vmf * (-Btf * sin(thetatf) + Gtf * cos(thetatf));
          dQt_dVmt_dVmf = (-Btf * cos(thetatf) + Gtf * sin(thetatf));

          /* dQt_dthetaf = Vmf*Vmt*(-Btf*sin(thetatf) - Gtf*cos(thetatf)); */
          dQt_dthetaf_dthetat =
            Vmf * Vmt * (-Btf * cos(thetatf) + Gtf * sin(thetatf));
          dQt_dthetaf_dVmt = Vmf * (-Btf * sin(thetatf) - Gtf * cos(thetatf));
          dQt_dthetaf_dthetaf =
            Vmf * Vmt * (Btf * cos(thetatf) - Gtf * sin(thetatf));
          dQt_dthetaf_dVmf = Vmt * (-Btf * sin(thetatf) - Gtf * cos(thetatf));

          /* dQt_dVmf = Vmt*(-Btf*cos(thetatf) + Gtf*sin(thetatf)); */
          dQt_dVmf_dthetat = Vmt * (Btf * sin(thetatf) + Gtf * cos(thetatf));
          dQt_dVmf_dVmt = (-Btf * cos(thetatf) + Gtf * sin(thetatf));
          dQt_dVmf_dthetaf = Vmt * (-Btf * sin(thetatf) - Gtf * cos(thetatf));
          dQt_dVmf_dVmf = 0.0;

          // row[0] = xidxt[i] - nxsparse;
          // row[1] = xidxt[i] + 1 - nxsparse;
          // col[0] = xidxt[i] - nxsparse;
          // col[1] = xidxt[i] + 1 - nxsparse;
          // col[2] = xidxf[i] - nxsparse;
          // col[3] = xidxf[i] + 1 - nxsparse;

          gloc = geqidxt[i];

          val[0] = val[1] = val[2] = val[3] = val[4] = val[5] = val[6] = val[7] =
            0.0;

          val[0] = lambda_dev[gloc] * dPt_dthetat_dthetat +
            lambda_dev[gloc + 1] * dQt_dthetat_dthetat;
          val[1] = lambda_dev[gloc] * dPt_dthetat_dVmt +
            lambda_dev[gloc + 1] * dQt_dthetat_dVmt;
          val[2] = lambda_dev[gloc] * dPt_dthetat_dthetaf +
            lambda_dev[gloc + 1] * dQt_dthetat_dthetaf;
          val[3] = lambda_dev[gloc] * dPt_dthetat_dVmf +
            lambda_dev[gloc + 1] * dQt_dthetat_dVmf;

          val[4] = lambda_dev[gloc] * dPt_dVmt_dthetat +
            lambda_dev[gloc + 1] * dQt_dVmt_dthetat;
          val[5] = lambda_dev[gloc] * dPt_dVmt_dVmt +
            lambda_dev[gloc + 1] * dQt_dVmt_dVmt;
          val[6] = lambda_dev[gloc] * dPt_dVmt_dthetaf +
            lambda_dev[gloc + 1] * dQt_dVmt_dthetaf;
          val[7] = lambda_dev[gloc] * dPt_dVmt_dVmf +
            lambda_dev[gloc + 1] * dQt_dVmt_dVmf;

          // to-to diagonal bus entries
          RAJA::atomicAdd<exago_raja_atomic>
            (&MHSS_dev[b_hesssp_idx[tbusidx] + 0], val[0]);
          RAJA::atomicAdd<exago_raja_atomic>
            (&MHSS_dev[b_hesssp_idx[tbusidx] + 1], val[1]);
          // not in upper triangle
          // RAJA::atomicAdd<exago_raja_atomic>
          //   (&MHSS_dev[b_hesssp_idx[tbusidx] + 2], val[4]);
          RAJA::atomicAdd<exago_raja_atomic>
            (&MHSS_dev[b_hesssp_idx[tbusidx] + 2], val[5]);

          // off-center to-from entries
          if (xidxf[i] > xidxt[i]) {
            RAJA::atomicAdd<exago_raja_atomic>
              (&MHSS_dev[ln_hessp_idx[i] + 0], val[2]);
            RAJA::atomicAdd<exago_raja_atomic>
              (&MHSS_dev[ln_hessp_idx[i] + 1], val[3]);
            RAJA::atomicAdd<exago_raja_atomic>
              (&MHSS_dev[ln_hessp_idx[i] + 2], val[6]);
            RAJA::atomicAdd<exago_raja_atomic>
              (&MHSS_dev[ln_hessp_idx[i] + 3], val[7]);
          }

          // row[0] = xidxf[i] - nxsparse;
          // row[1] = xidxf[i] + 1 - nxsparse;
          // col[0] = xidxt[i] - nxsparse;
          // col[1] = xidxt[i] + 1 - nxsparse;
          // col[2] = xidxf[i] - nxsparse;
          // col[3] = xidxf[i] + 1 - nxsparse;

          val[0] = val[1] = val[2] = val[3] = val[4] = val[5] = val[6] = val[7] =
            0.0;

          val[0] = lambda_dev[gloc] * dPt_dthetaf_dthetat +
            lambda_dev[gloc + 1] * dQt_dthetaf_dthetat;
          val[1] = lambda_dev[gloc] * dPt_dthetaf_dVmt +
            lambda_dev[gloc + 1] * dQt_dthetaf_dVmt;
          val[2] = lambda_dev[gloc] * dPt_dthetaf_dthetaf +
            lambda_dev[gloc + 1] * dQt_dthetaf_dthetaf;
          val[3] = lambda_dev[gloc] * dPt_dthetaf_dVmf +
            lambda_dev[gloc + 1] * dQt_dthetaf_dVmf;

          val[4] = lambda_dev[gloc] * dPt_dVmf_dthetat +
            lambda_dev[gloc + 1] * dQt_dVmf_dthetat;
          val[5] = lambda_dev[gloc] * dPt_dVmf_dVmt +
            lambda_dev[gloc + 1] * dQt_dVmf_dVmt;
          val[6] = lambda_dev[gloc] * dPt_dVmf_dthetaf +
            lambda_dev[gloc + 1] * dQt_dVmf_dthetaf;
          val[7] = lambda_dev[gloc] * dPt_dVmf_dVmf +
            lambda_dev[gloc + 1] * dQt_dVmf_dVmf;

          // from-from bus entries
          RAJA::atomicAdd<exago_raja_atomic>
            (&MHSS_dev[b_hesssp_idx[fbusidx] + 0], val[0]);
          RAJA::atomicAdd<exago_raja_atomic>
            (&MHSS_dev[b_hesssp_idx[fbusidx] + 1], val[1]);
          // RAJA::atomicAdd<exago_raja_atomic>
          //   (&MHSS_dev[b_hesssp_idx[fbusidx] + 2], val[4]);
          RAJA::atomicAdd<exago_raja_atomic>
            (&MHSS_dev[b_hesssp_idx[fbusidx] + 2], val[7]);

          // off-center from-to entries
          if (xidxt[i] > xidxf[i]) {
            RAJA::atomicAdd<exago_raja_atomic>
              (&MHSS_dev[ln_hessp_idx[i] + 0], val[0]);
            RAJA::atomicAdd<exago_raja_atomic>
              (&MHSS_dev[ln_hessp_idx[i] + 1], val[1]);
            RAJA::atomicAdd<exago_raja_atomic>
              (&MHSS_dev[ln_hessp_idx[i] + 2], val[4]);
            RAJA::atomicAdd<exago_raja_atomic>
              (&MHSS_dev[ln_hessp_idx[i] + 3], val[5]);
          }
        });

    /* Loadloss contributions - 2 contributions expected */
    if (opflow->include_loadloss_variables) {
      int *l_hesssp_idx = loadparams->hesssp_idx_dev_;
      RAJA::forall<exago_raja_exec>(
          RAJA::RangeSegment(0, loadparams->nload),
          RAJA_LAMBDA(RAJA::Index_type i) {
            MHSS_dev[l_hesssp_idx[i]] = 0.0;
            MHSS_dev[l_hesssp_idx[i] + 1] = 0.0;
          });
    }

    if (!opflow->ignore_lineflow_constraints) {
      int *linelimidx = lineparams->linelimidx_dev_;
      int *gineqidx = lineparams->gineqidx_dev_;
      RAJA::forall<exago_raja_exec>(
          RAJA::RangeSegment(0, lineparams->nlinelim),
          RAJA_LAMBDA(RAJA::Index_type i) {
            int j = linelimidx[i];
            int gloc;
            // int row[2], col[4];
            double val[8];

            double Pf, Qf, Pt, Qt;
            double thetaf = x_dev[xidxf[j]], Vmf = x_dev[xidxf[j] + 1];
            double thetat = x_dev[xidxt[j]], Vmt = x_dev[xidxt[j] + 1];
            double thetaft = thetaf - thetat;
            double thetatf = thetat - thetaf;
            double Gff = Gff_arr[j], Bff = Bff_arr[j];
            double Gft = Gft_arr[j], Bft = Bft_arr[j];
            double Gtf = Gtf_arr[j], Btf = Btf_arr[j];
            double Gtt = Gtt_arr[j], Btt = Btt_arr[j];
            int fbusidx(busf_idx[j]), tbusidx(bust_idx[j]);

            Pf = Gff * Vmf * Vmf +
              Vmf * Vmt * (Gft * cos(thetaft) + Bft * sin(thetaft));
            Qf = -Bff * Vmf * Vmf +
              Vmf * Vmt * (-Bft * cos(thetaft) + Gft * sin(thetaft));

            Pt = Gtt * Vmt * Vmt +
              Vmt * Vmf * (Gtf * cos(thetatf) + Btf * sin(thetatf));
            Qt = -Btt * Vmt * Vmt +
              Vmt * Vmf * (-Btf * cos(thetatf) + Gtf * sin(thetatf));

            double dSf2_dPf, dSf2_dQf, dSt2_dPt, dSt2_dQt;

            dSf2_dPf = 2. * Pf;
            dSf2_dQf = 2. * Qf;
            dSt2_dPt = 2. * Pt;
            dSt2_dQt = 2. * Qt;

            double dPf_dthetaf, dPf_dVmf, dPf_dthetat, dPf_dVmt;
            double dQf_dthetaf, dQf_dVmf, dQf_dthetat, dQf_dVmt;
            double dPt_dthetaf, dPt_dVmf, dPt_dthetat, dPt_dVmt;
            double dQt_dthetaf, dQt_dVmf, dQt_dthetat, dQt_dVmt;

            dPf_dthetaf = Vmf * Vmt * (-Gft * sin(thetaft) + Bft * cos(thetaft));
            dPf_dVmf =
              2. * Gff * Vmf + Vmt * (Gft * cos(thetaft) + Bft * sin(thetaft));
            dPf_dthetat = Vmf * Vmt * (Gft * sin(thetaft) - Bft * cos(thetaft));
            dPf_dVmt = Vmf * (Gft * cos(thetaft) + Bft * sin(thetaft));

            dQf_dthetaf = Vmf * Vmt * (Bft * sin(thetaft) + Gft * cos(thetaft));
            dQf_dVmf =
              -2. * Bff * Vmf + Vmt * (-Bft * cos(thetaft) + Gft * sin(thetaft));
            dQf_dthetat = Vmf * Vmt * (-Bft * sin(thetaft) - Gft * cos(thetaft));
            dQf_dVmt = Vmf * (-Bft * cos(thetaft) + Gft * sin(thetaft));

            dPt_dthetat = Vmt * Vmf * (-Gtf * sin(thetatf) + Btf * cos(thetatf));
            dPt_dVmt =
              2. * Gtt * Vmt + Vmf * (Gtf * cos(thetatf) + Btf * sin(thetatf));
            dPt_dthetaf = Vmt * Vmf * (Gtf * sin(thetatf) - Btf * cos(thetatf));
            dPt_dVmf = Vmt * (Gtf * cos(thetatf) + Btf * sin(thetatf));

            dQt_dthetat = Vmt * Vmf * (Btf * sin(thetatf) + Gtf * cos(thetatf));
            dQt_dVmt =
              -2. * Btt * Vmt + Vmf * (-Btf * cos(thetatf) + Gtf * sin(thetatf));
            dQt_dthetaf = Vmt * Vmf * (-Btf * sin(thetatf) - Gtf * cos(thetatf));
            dQt_dVmf = Vmt * (-Btf * cos(thetatf) + Gtf * sin(thetatf));

            double d2Pf_dthetaf_dthetaf, d2Pf_dthetaf_dVmf, d2Pf_dthetaf_dthetat,
              d2Pf_dthetaf_dVmt;
            double d2Pf_dVmf_dthetaf, d2Pf_dVmf_dVmf, d2Pf_dVmf_dthetat,
              d2Pf_dVmf_dVmt;
            double d2Pf_dthetat_dthetaf, d2Pf_dthetat_dVmf, d2Pf_dthetat_dthetat,
              d2Pf_dthetat_dVmt;
            double d2Pf_dVmt_dthetaf, d2Pf_dVmt_dVmf, d2Pf_dVmt_dthetat,
              d2Pf_dVmt_dVmt;

            /* dPf_dthetaf = Vmf*Vmt*(-Gft*sin(thetaft) + Bft*cos(thetaft)); */
            d2Pf_dthetaf_dthetaf =
              -Vmf * Vmt * (Gft * cos(thetaft) + Bft * sin(thetaft));
            d2Pf_dthetaf_dVmf = Vmt * (-Gft * sin(thetaft) + Bft * cos(thetaft));
            d2Pf_dthetaf_dthetat =
              Vmf * Vmt * (Gft * cos(thetaft) + Bft * sin(thetaft));
            d2Pf_dthetaf_dVmt = Vmf * (-Gft * sin(thetaft) + Bft * cos(thetaft));

            /* dPf_Vmf  = 2*Gff*Vmf + Vmt*(Gft*cos(thetaft) + Bft*sin(thetaft)); */
            d2Pf_dVmf_dthetaf = Vmt * (-Gft * sin(thetaft) + Bft * cos(thetaft));
            d2Pf_dVmf_dVmf = 2 * Gff;
            d2Pf_dVmf_dthetat = Vmt * (Gft * sin(thetaft) - Bft * cos(thetaft));
            d2Pf_dVmf_dVmt = (Gft * cos(thetaft) + Bft * sin(thetaft));

            /* dPf_dthetat = Vmf*Vmt*(Gft*sin(thetaft) - Bft*cos(thetaft)); */
            d2Pf_dthetat_dthetaf =
              Vmf * Vmt * (Gft * cos(thetaft) + Bft * sin(thetaft));
            d2Pf_dthetat_dVmf = Vmt * (Gft * sin(thetaft) - Bft * cos(thetaft));
            d2Pf_dthetat_dthetat =
              Vmf * Vmt * (-Gft * cos(thetaft) - Bft * sin(thetaft));
            d2Pf_dthetat_dVmt = Vmf * (Gft * sin(thetaft) - Bft * cos(thetaft));

            /* dPf_dVmt = Vmf*(Gft*cos(thetaft) + Bft*sin(thetaft)); */
            d2Pf_dVmt_dthetaf = Vmf * (-Gft * sin(thetaft) + Bft * cos(thetaft));
            d2Pf_dVmt_dVmf = (Gft * cos(thetaft) + Bft * sin(thetaft));
            d2Pf_dVmt_dthetat = Vmf * (Gft * sin(thetaft) - Bft * cos(thetaft));
            d2Pf_dVmt_dVmt = 0.0;

            double d2Qf_dthetaf_dthetaf, d2Qf_dthetaf_dVmf, d2Qf_dthetaf_dthetat,
              d2Qf_dthetaf_dVmt;
            double d2Qf_dVmf_dthetaf, d2Qf_dVmf_dVmf, d2Qf_dVmf_dthetat,
              d2Qf_dVmf_dVmt;
            double d2Qf_dthetat_dthetaf, d2Qf_dthetat_dVmf, d2Qf_dthetat_dthetat,
              d2Qf_dthetat_dVmt;
            double d2Qf_dVmt_dthetaf, d2Qf_dVmt_dVmf, d2Qf_dVmt_dthetat,
              d2Qf_dVmt_dVmt;

            /* dQf_dthetaf = Vmf*Vmt*(Bft*sin(thetaft) + Gft*cos(thetaft)); */
            d2Qf_dthetaf_dthetaf =
              Vmf * Vmt * (Bft * cos(thetaft) - Gft * sin(thetaft));
            d2Qf_dthetaf_dVmf = Vmt * (Bft * sin(thetaft) + Gft * cos(thetaft));
            d2Qf_dthetaf_dthetat =
              Vmf * Vmt * (-Bft * cos(thetaft) + Gft * sin(thetaft));
            d2Qf_dthetaf_dVmt = Vmf * (Bft * sin(thetaft) + Gft * cos(thetaft));

            /* dQf_dVmf = -2*Bff*Vmf + Vmt*(-Bft*cos(thetaft) + Gft*sin(thetaft));
             */
            d2Qf_dVmf_dthetaf = Vmt * (Bft * sin(thetaft) + Gft * cos(thetaft));
            d2Qf_dVmf_dVmf = -2 * Bff;
            d2Qf_dVmf_dthetat = Vmt * (-Bft * sin(thetaft) - Gft * cos(thetaft));
            d2Qf_dVmf_dVmt = (-Bft * cos(thetaft) + Gft * sin(thetaft));

            /* dQf_dthetat = Vmf*Vmt*(-Bft*sin(thetaft) - Gft*cos(thetaft)); */
            d2Qf_dthetat_dthetaf =
              Vmf * Vmt * (-Bft * cos(thetaft) + Gft * sin(thetaft));
            d2Qf_dthetat_dVmf = Vmt * (-Bft * sin(thetaft) - Gft * cos(thetaft));
            d2Qf_dthetat_dthetat =
              Vmf * Vmt * (Bft * cos(thetaft) - Gft * sin(thetaft));
            d2Qf_dthetat_dVmt = Vmf * (-Bft * sin(thetaft) - Gft * cos(thetaft));

            /* dQf_dVmt = Vmf*(-Bft*cos(thetaft) + Gft*sin(thetaft)); */
            d2Qf_dVmt_dthetaf = Vmf * (Bft * sin(thetaft) + Gft * cos(thetaft));
            d2Qf_dVmt_dVmf = (-Bft * cos(thetaft) + Gft * sin(thetaft));
            d2Qf_dVmt_dthetat = Vmf * (-Bft * sin(thetaft) - Gft * cos(thetaft));
            d2Qf_dVmt_dVmt = 0.0;

            double d2Pt_dthetat_dthetat, d2Pt_dthetat_dVmt, d2Pt_dthetat_dthetaf,
              d2Pt_dthetat_dVmf;
            double d2Pt_dVmt_dthetat, d2Pt_dVmt_dVmt, d2Pt_dVmt_dthetaf,
              d2Pt_dVmt_dVmf;
            double d2Pt_dthetaf_dthetat, d2Pt_dthetaf_dVmt, d2Pt_dthetaf_dthetaf,
              d2Pt_dthetaf_dVmf;
            double d2Pt_dVmf_dthetat, d2Pt_dVmf_dVmt, d2Pt_dVmf_dthetaf,
              d2Pt_dVmf_dVmf;

            /* dPt_dthetat = Vmf*Vmt*(-Gtf*sin(thetatf) + Btf*cos(thetatf)); */
            d2Pt_dthetat_dthetat =
              Vmf * Vmt * (-Gtf * cos(thetatf) - Btf * sin(thetatf));
            d2Pt_dthetat_dVmt = Vmf * (-Gtf * sin(thetatf) + Btf * cos(thetatf));
            d2Pt_dthetat_dthetaf =
              Vmf * Vmt * (Gtf * cos(thetatf) + Btf * sin(thetatf));
            d2Pt_dthetat_dVmf = Vmt * (-Gtf * sin(thetatf) + Btf * cos(thetatf));

            /* dPt_Vmt  = 2*Gtt*Vmt + Vmf*(Gtf*cos(thetatf) + Btf*sin(thetatf)); */
            d2Pt_dVmt_dthetat = Vmf * (-Gtf * sin(thetatf) + Bft * cos(thetatf));
            d2Pt_dVmt_dVmt = 2 * Gtt;
            d2Pt_dVmt_dthetaf = Vmf * (Gtf * sin(thetatf) - Btf * cos(thetatf));
            d2Pt_dVmt_dVmf = (Gtf * cos(thetatf) + Btf * sin(thetatf));

            /* dPt_dthetaf = Vmf*Vmt*(Gtf*sin(thetatf) - Btf*cos(thetatf)); */
            d2Pt_dthetaf_dthetat =
              Vmf * Vmt * (Gtf * cos(thetatf) + Btf * sin(thetatf));
            d2Pt_dthetaf_dVmt = Vmf * (Gtf * sin(thetatf) - Btf * cos(thetatf));
            d2Pt_dthetaf_dthetaf =
              Vmf * Vmt * (-Gtf * cos(thetatf) - Btf * sin(thetatf));
            d2Pt_dthetaf_dVmf = Vmt * (Gtf * sin(thetatf) - Btf * cos(thetatf));

            /* dPt_dVmf = Vmt*(Gtf*cos(thetatf) + Btf*sin(thetatf)); */
            d2Pt_dVmf_dthetat = Vmt * (-Gtf * sin(thetatf) + Btf * cos(thetatf));
            d2Pt_dVmf_dVmt = (Gtf * cos(thetatf) + Btf * sin(thetatf));
            d2Pt_dVmf_dthetaf = Vmt * (Gtf * sin(thetatf) - Btf * cos(thetatf));
            d2Pt_dVmf_dVmf = 0.0;

            double d2Qt_dthetaf_dthetaf, d2Qt_dthetaf_dVmf, d2Qt_dthetaf_dthetat,
              d2Qt_dthetaf_dVmt;
            double d2Qt_dVmf_dthetaf, d2Qt_dVmf_dVmf, d2Qt_dVmf_dthetat,
              d2Qt_dVmf_dVmt;
            double d2Qt_dthetat_dthetaf, d2Qt_dthetat_dVmf, d2Qt_dthetat_dthetat,
              d2Qt_dthetat_dVmt;
            double d2Qt_dVmt_dthetaf, d2Qt_dVmt_dVmf, d2Qt_dVmt_dthetat,
              d2Qt_dVmt_dVmt;

            /* dQt_dthetat = Vmf*Vmt*(Btf*sin(thetatf) + Gtf*cos(thetatf)); */
            d2Qt_dthetat_dthetat =
              Vmf * Vmt * (Btf * cos(thetatf) - Gtf * sin(thetatf));
            d2Qt_dthetat_dVmt = Vmf * (Btf * sin(thetatf) + Gtf * cos(thetatf));
            d2Qt_dthetat_dthetaf =
              Vmf * Vmt * (-Btf * cos(thetatf) + Gtf * sin(thetatf));
            d2Qt_dthetat_dVmf = Vmf * (Btf * sin(thetatf) + Gtf * cos(thetatf));

            /* dQt_dVmt = -2*Btt*Vmt + Vmf*(-Btf*cos(thetatf) + Gtf*sin(thetatf));
             */
            d2Qt_dVmt_dthetat = Vmf * (Btf * sin(thetatf) + Gtf * cos(thetatf));
            d2Qt_dVmt_dVmt = -2 * Btt;
            d2Qt_dVmt_dthetaf = Vmf * (-Btf * sin(thetatf) + Gtf * cos(thetatf));
            d2Qt_dVmt_dVmf = (-Btf * cos(thetatf) + Gtf * sin(thetatf));

            /* dQt_dthetaf = Vmf*Vmt*(-Btf*sin(thetatf) - Gtf*cos(thetatf)); */
            d2Qt_dthetaf_dthetat =
              Vmf * Vmt * (-Btf * cos(thetatf) + Gtf * sin(thetatf));
            d2Qt_dthetaf_dVmt = Vmf * (-Btf * sin(thetatf) - Gtf * cos(thetatf));
            d2Qt_dthetaf_dthetaf =
              Vmf * Vmt * (Btf * cos(thetatf) - Gtf * sin(thetatf));
            d2Qt_dthetaf_dVmf = Vmt * (-Btf * sin(thetatf) - Gtf * cos(thetatf));

            /* dQt_dVmf = Vmt*(-Btf*cos(thetatf) + Gtf*sin(thetatf)); */
            d2Qt_dVmf_dthetat = Vmt * (Btf * sin(thetatf) + Gtf * cos(thetatf));
            d2Qt_dVmf_dVmt = (-Btf * cos(thetatf) + Gtf * sin(thetatf));
            d2Qt_dVmf_dthetaf = Vmt * (-Btf * sin(thetatf) - Gtf * cos(thetatf));
            d2Qt_dVmf_dVmf = 0.0;

            double d2Sf2_dthetaf_dthetaf = 0.0, d2Sf2_dthetaf_dVmf = 0.0,
              d2Sf2_dthetaf_dthetat = 0.0, d2Sf2_dthetaf_dVmt = 0.0;
            double d2St2_dthetaf_dthetaf = 0.0, d2St2_dthetaf_dVmf = 0.0,
              d2St2_dthetaf_dthetat = 0.0, d2St2_dthetaf_dVmt = 0.0;

            d2Sf2_dthetaf_dthetaf =
              2 * dPf_dthetaf * dPf_dthetaf + dSf2_dPf * d2Pf_dthetaf_dthetaf +
              2 * dQf_dthetaf * dQf_dthetaf + dSf2_dQf * d2Qf_dthetaf_dthetaf;
            d2Sf2_dthetaf_dVmf =
              2 * dPf_dVmf * dPf_dthetaf + dSf2_dPf * d2Pf_dthetaf_dVmf +
              2 * dQf_dVmf * dQf_dthetaf + dSf2_dQf * d2Qf_dthetaf_dVmf;
            d2Sf2_dthetaf_dthetat =
              2 * dPf_dthetat * dPf_dthetaf + dSf2_dPf * d2Pf_dthetaf_dthetat +
              2 * dQf_dthetat * dQf_dthetaf + dSf2_dQf * d2Qf_dthetaf_dthetat;
            d2Sf2_dthetaf_dVmt =
              2 * dPf_dVmt * dPf_dthetaf + dSf2_dPf * d2Pf_dthetaf_dVmt +
              2 * dQf_dVmt * dQf_dthetaf + dSf2_dQf * d2Qf_dthetaf_dVmt;

            d2St2_dthetaf_dthetaf =
              2 * dPt_dthetaf * dPt_dthetaf + dSt2_dPt * d2Pt_dthetaf_dthetaf +
              2 * dQt_dthetaf * dQt_dthetaf + dSt2_dQt * d2Qt_dthetaf_dthetaf;
            d2St2_dthetaf_dVmf =
              2 * dPt_dVmf * dPt_dthetaf + dSt2_dPt * d2Pt_dthetaf_dVmf +
              2 * dQt_dVmf * dQt_dthetaf + dSt2_dQt * d2Qt_dthetaf_dVmf;
            d2St2_dthetaf_dthetat =
              2 * dPt_dthetat * dPt_dthetaf + dSt2_dPt * d2Pt_dthetaf_dthetat +
              2 * dQt_dthetat * dQt_dthetaf + dSt2_dQt * d2Qt_dthetaf_dthetat;
            d2St2_dthetaf_dVmt =
              2 * dPt_dVmt * dPt_dthetaf + dSt2_dPt * d2Pt_dthetaf_dVmt +
              2 * dQt_dVmt * dQt_dthetaf + dSt2_dQt * d2Qt_dthetaf_dVmt;

            val[0] = val[1] = val[2] = val[3] = 0.0;

            // row[0] = xidxf[j] - nxsparse;
            // col[0] = xidxf[j] - nxsparse;
            // col[1] = xidxf[j] + 1 - nxsparse;
            // col[2] = xidxt[j] - nxsparse;
            // col[3] = xidxt[j] + 1 - nxsparse;

            gloc = gineqidx[i];

            val[0] = lambda_dev[gloc] * d2Sf2_dthetaf_dthetaf +
              lambda_dev[gloc + 1] * d2St2_dthetaf_dthetaf;
            val[1] = lambda_dev[gloc] * d2Sf2_dthetaf_dVmf +
              lambda_dev[gloc + 1] * d2St2_dthetaf_dVmf;
            val[2] = lambda_dev[gloc] * d2Sf2_dthetaf_dthetat +
              lambda_dev[gloc + 1] * d2St2_dthetaf_dthetat;
            val[3] = lambda_dev[gloc] * d2Sf2_dthetaf_dVmt +
              lambda_dev[gloc + 1] * d2St2_dthetaf_dVmt;

            RAJA::atomicAdd<exago_raja_atomic>
              (&MHSS_dev[fbusidx + 0], val[0]);
            RAJA::atomicAdd<exago_raja_atomic>
              (&MHSS_dev[fbusidx + 1], val[1]);

            if (xidxt[j] > xidxf[j]) {
              RAJA::atomicAdd<exago_raja_atomic>
                (&MHSS_dev[ln_hessp_idx[j] + 0], val[2]);
              RAJA::atomicAdd<exago_raja_atomic>
                (&MHSS_dev[ln_hessp_idx[j] + 1], val[3]);
            }

            double d2Sf2_dVmf_dthetaf, d2Sf2_dVmf_dVmf, d2Sf2_dVmf_dthetat,
              d2Sf2_dVmf_dVmt;
            double d2St2_dVmf_dthetaf, d2St2_dVmf_dVmf, d2St2_dVmf_dthetat,
              d2St2_dVmf_dVmt;

            d2Sf2_dVmf_dthetaf =
              2 * dPf_dthetaf * dPf_dVmf + dSf2_dPf * d2Pf_dVmf_dthetaf +
              2 * dQf_dthetaf * dQf_dVmf + dSf2_dQf * d2Qf_dVmf_dthetaf;
            d2Sf2_dVmf_dVmf = 2 * dPf_dVmf * dPf_dVmf + dSf2_dPf * d2Pf_dVmf_dVmf +
              2 * dQf_dVmf * dQf_dVmf + dSf2_dQf * d2Qf_dVmf_dVmf;
            d2Sf2_dVmf_dthetat =
              2 * dPf_dthetat * dPf_dVmf + dSf2_dPf * d2Pf_dVmf_dthetat +
              2 * dQf_dthetat * dQf_dVmf + dSf2_dQf * d2Qf_dVmf_dthetat;
            d2Sf2_dVmf_dVmt = 2 * dPf_dVmt * dPf_dVmf + dSf2_dPf * d2Pf_dVmf_dVmt +
              2 * dQf_dVmt * dQf_dVmf + dSf2_dQf * d2Qf_dVmf_dVmt;

            d2St2_dVmf_dthetaf =
              2 * dPt_dthetaf * dPt_dVmf + dSt2_dPt * d2Pt_dVmf_dthetaf +
              2 * dQt_dthetaf * dQt_dVmf + dSt2_dQt * d2Qt_dVmf_dthetaf;
            d2St2_dVmf_dVmf = 2 * dPt_dVmf * dPt_dVmf + dSt2_dPt * d2Pt_dVmf_dVmf +
              2 * dQt_dVmf * dQt_dVmf + dSt2_dQt * d2Qt_dVmf_dVmf;
            d2St2_dVmf_dthetat =
              2 * dPt_dthetat * dPt_dVmf + dSt2_dPt * d2Pt_dVmf_dthetat +
              2 * dQt_dthetat * dQt_dVmf + dSt2_dQt * d2Qt_dVmf_dthetat;
            d2St2_dVmf_dVmt = 2 * dPt_dVmt * dPt_dVmf + dSt2_dPt * d2Pt_dVmf_dVmt +
              2 * dQt_dVmt * dQt_dVmf + dSt2_dQt * d2Qt_dVmf_dVmt;

            val[0] = val[1] = val[2] = val[3] = 0.0;

            // row[0] = xidxf[j] + 1 - nxsparse;
            // col[0] = xidxf[j] - nxsparse;
            // col[1] = xidxf[j] + 1 - nxsparse;
            // col[2] = xidxt[j] - nxsparse;
            // col[3] = xidxt[j] + 1 - nxsparse;


            val[0] = lambda_dev[gloc] * d2Sf2_dVmf_dthetaf +
              lambda_dev[gloc + 1] * d2St2_dVmf_dthetaf;
            val[1] = lambda_dev[gloc] * d2Sf2_dVmf_dVmf +
              lambda_dev[gloc + 1] * d2St2_dVmf_dVmf;
            val[2] = lambda_dev[gloc] * d2Sf2_dVmf_dthetat +
              lambda_dev[gloc + 1] * d2St2_dVmf_dthetat;
            val[3] = lambda_dev[gloc] * d2Sf2_dVmf_dVmt +
              lambda_dev[gloc + 1] * d2St2_dVmf_dVmt;

            // not in upper triangle
            // RAJA::atomicAdd<exago_raja_atomic>
            //   (&MHSS_dev[fbusidx + 2], val[0]);
            RAJA::atomicAdd<exago_raja_atomic>
              (&MHSS_dev[fbusidx + 2], val[1]);

            if (xidxt[j] > xidxf[j]) {
              RAJA::atomicAdd<exago_raja_atomic>
                (&MHSS_dev[ln_hessp_idx[j] + 2], val[2]);
              RAJA::atomicAdd<exago_raja_atomic>
                (&MHSS_dev[ln_hessp_idx[j] + 3], val[3]);
            }

            //    ierr = MatSetValues(H,1,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);

            double d2Sf2_dthetat_dthetaf, d2Sf2_dthetat_dVmf, d2Sf2_dthetat_dthetat,
              d2Sf2_dthetat_dVmt;
            double d2St2_dthetat_dthetaf, d2St2_dthetat_dVmf, d2St2_dthetat_dthetat,
              d2St2_dthetat_dVmt;

            d2Sf2_dthetat_dthetaf =
              2 * dPf_dthetaf * dPf_dthetat + dSf2_dPf * d2Pf_dthetat_dthetaf +
              2 * dQf_dthetat * dQf_dthetaf + dSf2_dQf * d2Qf_dthetat_dthetaf;
            d2Sf2_dthetat_dVmf =
              2 * dPf_dVmf * dPf_dthetat + dSf2_dPf * d2Pf_dthetat_dVmf +
              2 * dQf_dthetat * dQf_dVmf + dSf2_dQf * d2Qf_dthetat_dVmf;
            d2Sf2_dthetat_dthetat =
              2 * dPf_dthetat * dPf_dthetat + dSf2_dPf * d2Pf_dthetat_dthetat +
              2 * dQf_dthetat * dQf_dthetat + dSf2_dQf * d2Qf_dthetat_dthetat;
            d2Sf2_dthetat_dVmt =
              2 * dPf_dVmt * dPf_dthetat + dSf2_dPf * d2Pf_dthetat_dVmt +
              2 * dQf_dthetat * dQf_dVmt + dSf2_dQf * d2Qf_dthetat_dVmt;

            d2St2_dthetat_dthetaf =
              2 * dPt_dthetaf * dPt_dthetat + dSt2_dPt * d2Pt_dthetat_dthetaf +
              2 * dQt_dthetaf * dQt_dthetat + dSt2_dQt * d2Qt_dthetat_dthetaf;
            d2St2_dthetat_dVmf =
              2 * dPt_dVmf * dPt_dthetat + dSt2_dPt * d2Pt_dthetat_dVmf +
              2 * dQt_dVmf * dQt_dthetat + dSt2_dQt * d2Qt_dthetat_dVmf;
            d2St2_dthetat_dthetat =
              2 * dPt_dthetat * dPt_dthetat + dSt2_dPt * d2Pt_dthetat_dthetat +
              2 * dQt_dthetat * dQt_dthetat + dSt2_dQt * d2Qt_dthetat_dthetat;
            d2St2_dthetat_dVmt =
              2 * dPt_dVmt * dPt_dthetat + dSt2_dPt * d2Pt_dthetat_dVmt +
              2 * dQt_dVmt * dQt_dthetat + dSt2_dQt * d2Qt_dthetat_dVmt;

            val[0] = val[1] = val[2] = val[3] = 0.0;

            // row[0] = xidxt[j] - nxsparse;
            // col[0] = xidxf[j] - nxsparse;
            // col[1] = xidxf[j] + 1 - nxsparse;
            // col[2] = xidxt[j] - nxsparse;
            // col[3] = xidxt[j] + 1 - nxsparse;


            val[0] = lambda_dev[gloc] * d2Sf2_dthetat_dthetaf +
              lambda_dev[gloc + 1] * d2St2_dthetat_dthetaf;
            val[1] = lambda_dev[gloc] * d2Sf2_dthetat_dVmf +
              lambda_dev[gloc + 1] * d2St2_dthetat_dVmf;
            val[2] = lambda_dev[gloc] * d2Sf2_dthetat_dthetat +
              lambda_dev[gloc + 1] * d2St2_dthetat_dthetat;
            val[3] = lambda_dev[gloc] * d2Sf2_dthetat_dVmt +
              lambda_dev[gloc + 1] * d2St2_dthetat_dVmt;

            RAJA::atomicAdd<exago_raja_atomic>
              (&MHSS_dev[tbusidx + 0], val[2]);
            RAJA::atomicAdd<exago_raja_atomic>
              (&MHSS_dev[tbusidx + 1], val[3]);

            if (xidxf[j] > xidxt[j]) {
              RAJA::atomicAdd<exago_raja_atomic>
                (&MHSS_dev[ln_hessp_idx[j] + 0], val[0]);
              RAJA::atomicAdd<exago_raja_atomic>
                (&MHSS_dev[ln_hessp_idx[j] + 1], val[1]);
            }
            //    ierr = MatSetValues(H,1,row,4,col,val,ADD_VALUES);CHKERRQ(ierr);

            double d2Sf2_dVmt_dthetaf, d2Sf2_dVmt_dVmf, d2Sf2_dVmt_dthetat,
              d2Sf2_dVmt_dVmt;
            double d2St2_dVmt_dthetaf, d2St2_dVmt_dVmf, d2St2_dVmt_dthetat,
              d2St2_dVmt_dVmt;

            d2Sf2_dVmt_dthetaf =
              2 * dPf_dthetaf * dPf_dVmt + dSf2_dPf * d2Pf_dVmt_dthetaf +
              2 * dQf_dthetaf * dQf_dVmt + dSf2_dQf * d2Qf_dVmt_dthetaf;
            d2Sf2_dVmt_dVmf = 2 * dPf_dVmf * dPf_dVmt + dSf2_dPf * d2Pf_dVmt_dVmf +
              2 * dQf_dVmf * dQf_dVmt + dSf2_dQf * d2Qf_dVmt_dVmf;
            d2Sf2_dVmt_dthetat =
              2 * dPf_dthetat * dPf_dVmt + dSf2_dPf * d2Pf_dVmt_dthetat +
              2 * dQf_dthetat * dQf_dVmt + dSf2_dQf * d2Qf_dVmt_dthetat;
            d2Sf2_dVmt_dVmt = 2 * dPf_dVmt * dPf_dVmt + dSf2_dPf * d2Pf_dVmt_dVmt +
              2 * dQf_dVmt * dQf_dVmt + dSf2_dQf * d2Qf_dVmt_dVmt;

            d2St2_dVmt_dthetaf =
              2 * dPt_dthetaf * dPt_dVmt + dSt2_dPt * d2Pt_dVmt_dthetaf +
              2 * dQt_dthetaf * dQt_dVmt + dSt2_dQt * d2Qt_dVmt_dthetaf;
            d2St2_dVmt_dVmf = 2 * dPt_dVmf * dPt_dVmt + dSt2_dPt * d2Pt_dVmt_dVmf +
              2 * dQt_dVmf * dQt_dVmt + dSt2_dQt * d2Qt_dVmt_dVmf;
            d2St2_dVmt_dthetat =
              2 * dPt_dthetat * dPt_dVmt + dSt2_dPt * d2Pt_dVmt_dthetat +
              2 * dQt_dthetat * dQt_dVmt + dSt2_dQt * d2Qt_dVmt_dthetat;
            d2St2_dVmt_dVmt = 2 * dPt_dVmt * dPt_dVmt + dSt2_dPt * d2Pt_dVmt_dVmt +
              2 * dQt_dVmt * dQt_dVmt + dSt2_dQt * d2Qt_dVmt_dVmt;

            val[0] = val[1] = val[2] = val[3] = 0.0;

            // row[0] = xidxt[j] + 1 - nxsparse;
            // col[0] = xidxf[j] - nxsparse;
            // col[1] = xidxf[j] + 1 - nxsparse;
            // col[2] = xidxt[j] - nxsparse;
            // col[3] = xidxt[j] + 1 - nxsparse;

            val[0] = lambda_dev[gloc] * d2Sf2_dVmt_dthetaf +
              lambda_dev[gloc + 1] * d2St2_dVmt_dthetaf;
            val[1] = lambda_dev[gloc] * d2Sf2_dVmt_dVmf +
              lambda_dev[gloc + 1] * d2St2_dVmt_dVmf;
            val[2] = lambda_dev[gloc] * d2Sf2_dVmt_dthetat +
              lambda_dev[gloc + 1] * d2St2_dVmt_dthetat;
            val[3] = lambda_dev[gloc] * d2Sf2_dVmt_dVmt +
              lambda_dev[gloc + 1] * d2St2_dVmt_dVmt;

            // not in upper triangle
            // RAJA::atomicAdd<exago_raja_atomic>
            //   (&MHSS_dev[tbusidx + 2], val[2]);
            RAJA::atomicAdd<exago_raja_atomic>
              (&MHSS_dev[tbusidx + 2], val[3]);
            
            RAJA::atomicAdd<exago_raja_atomic>
              (&MHSS_dev[ln_hessp_idx[j] + 2], val[0]);
            RAJA::atomicAdd<exago_raja_atomic>
              (&MHSS_dev[ln_hessp_idx[j] + 3], val[1]);
          });
    }
    
    int *iperm = pbpolrajahiopsparse->idx_hess_dev_;

    int *ipermout = (int *)d_allocator_.allocate(opflow->nnz_hesssp*sizeof(int));
    resmgr.copy(ipermout, iperm);
    
    RAJA::stable_sort_pairs<exago_raja_exec>
      (RAJA::make_span(ipermout, opflow->nnz_hesssp),
       RAJA::make_span(MHSS_dev, opflow->nnz_hesssp),
       RAJA::operators::less<int>{});

    if (debugmsg) {
      PrintTriplets("Hessian Values (GPU):",
                    opflow->nnz_hesssp, NULL, iHSS_dev, jHSS_dev, MHSS_dev);
    }
    
    d_allocator_.deallocate(ipermout);
    
    resmgr.memset(MHSS_dev, 0, opflow->nnz_hesssp*sizeof(double));

    ierr = VecGetArray(opflow->X, &x);
    CHKERRQ(ierr);

    // Copy from device to host
    umpire::Allocator h_allocator_ = resmgr.getAllocator("HOST");
    registerWith(x, opflow->nx, resmgr, h_allocator_);
    resmgr.copy((double *)x, (double *)x_dev);

    ierr = VecRestoreArray(opflow->X, &x);
    CHKERRQ(ierr);

    ierr = VecGetArray(opflow->Lambda, &lambda);

    registerWith(lambda, opflow->ncon, resmgr, h_allocator_);
    // copy lambda from device to host
    resmgr.copy((double *)lambda, (double *)lambda_dev);

    ierr = VecPlaceArray(opflow->Lambdae, lambda);
    CHKERRQ(ierr);
    if (opflow->Nconineq) {
      ierr = VecPlaceArray(opflow->Lambdai, lambda + opflow->nconeq);
      CHKERRQ(ierr);
    }

    /* Compute Hessian */
    ierr = (*opflow->modelops.computehessian)(
        opflow, opflow->X, opflow->Lambdae, opflow->Lambdai, opflow->Hes);
    CHKERRQ(ierr);

    ierr = VecResetArray(opflow->Lambdae);
    CHKERRQ(ierr);
    if (opflow->Nconineq) {
      ierr = VecResetArray(opflow->Lambdai);
      CHKERRQ(ierr);
    }

    ierr = VecRestoreArray(opflow->Lambda, &lambda);
    CHKERRQ(ierr);

    if (opflow->modelops.computeauxhessian) {
      ierr = VecGetArray(opflow->X, &x);
      CHKERRQ(ierr);
      ierr = (*opflow->modelops.computeauxhessian)(opflow, x, opflow->Hes,
                                                   opflow->userctx);
      CHKERRQ(ierr);
      ierr = VecRestoreArray(opflow->X, &x);
      CHKERRQ(ierr);

      ierr = MatAssemblyBegin(opflow->Hes, MAT_FINAL_ASSEMBLY);
      CHKERRQ(ierr);
      ierr = MatAssemblyEnd(opflow->Hes, MAT_FINAL_ASSEMBLY);
      CHKERRQ(ierr);
    }

    /* Copy over values */
    ierr = MatGetSize(opflow->Hes, &nrow, &nrow);
    CHKERRQ(ierr);

    values = pbpolrajahiopsparse->val_hess;

    for (i = 0; i < nrow; i++) {
      ierr = MatGetRow(opflow->Hes, i, &nvals, &cols, &vals);
      CHKERRQ(ierr);
      ctr = 0;
      for (j = 0; j < nvals; j++) {
        if (cols[j] >= i) { /* Upper triangle values (same as lower triangle) */
          values[ctr] = vals[j];
          ctr++;
        }
      }
      values += ctr;
      ierr = MatRestoreRow(opflow->Hes, i, &nvals, &cols, &vals);
      CHKERRQ(ierr);
    }

    // Copy over val_ineq to device
    resmgr.copy(MHSS_dev, pbpolrajahiopsparse->val_hess);

    if (debugmsg) {
      PrintTriplets("Hessian Values:",
                    opflow->nnz_hesssp, NULL, iHSS_dev, jHSS_dev, MHSS_dev);
    }

 }

  PetscFunctionReturn(0);
}

PetscErrorCode OPFLOWSolutionCallback_PBPOLRAJAHIOPSPARSE(
    OPFLOW opflow, const double *xsol, const double *z_L, const double *z_U,
    const double *gsol, const double *lamsol, double obj_value) {
  (void)z_L;
  (void)z_U;
  (void)obj_value;

  PetscErrorCode ierr;
  PetscScalar *x, *lam, *g;

  auto &resmgr = umpire::ResourceManager::getInstance();
  umpire::Allocator h_allocator_ = resmgr.getAllocator("HOST");

  /* FIXME: It is assumed that array arguments are on the device. If
     there are other posibilities, opflow->mem_space should be set
     properly and used to decide how to get at those arrays */

  ierr = VecGetArray(opflow->X, &x);
  CHKERRQ(ierr);
  /* Copy xsol from device to host */
  resmgr.copy(x, (double *)xsol);
  ierr = VecRestoreArray(opflow->X, &x);
  CHKERRQ(ierr);

  if (lamsol) {
    /* HIOP returns a NULL for lamsol - probably lamsol needs to be added to
     HIOP. Need to remove this condition once it is fixed
    */
    ierr = VecGetArray(opflow->Lambda, &lam);
    CHKERRQ(ierr);

    /* Create temporary vectors for copying values to HOST */
    double *lamsol_host =
        (double *)h_allocator_.allocate(opflow->ncon * sizeof(double));

    /* Copy lamsol from device to host */
    resmgr.copy(lamsol_host, (double *)lamsol);

    ierr = PetscMemcpy(lam, (double *)lamsol_host,
                       opflow->nconeq * sizeof(PetscScalar));
    CHKERRQ(ierr);
    if (opflow->Nconineq) {
      ierr = PetscMemcpy(lam + opflow->nconeq,
                         (double *)(lamsol_host + opflow->nconeq),
                         opflow->nconineq * sizeof(PetscScalar));
      CHKERRQ(ierr);
    }
    ierr = VecRestoreArray(opflow->Lambda, &lam);
    CHKERRQ(ierr);
    h_allocator_.deallocate(lamsol_host);
  } else {
    ierr = VecSet(opflow->Lambda, -9999.0);
    CHKERRQ(ierr);
  }

  if (gsol) {
    ierr = VecGetArray(opflow->G, &g);
    CHKERRQ(ierr);

    /* Create temporary vectors for copying values to HOST */
    double *gsol_host =
        (double *)h_allocator_.allocate(opflow->ncon * sizeof(double));

    /* Copy gsol from device to host */
    resmgr.copy(gsol_host, (double *)gsol);

    ierr = PetscMemcpy(g, (double *)gsol_host,
                       opflow->nconeq * sizeof(PetscScalar));
    CHKERRQ(ierr);
    if (opflow->Nconineq) {
      ierr = PetscMemcpy(g + opflow->nconeq,
                         (double *)(gsol_host + opflow->nconeq),
                         opflow->nconineq * sizeof(PetscScalar));
      CHKERRQ(ierr);
    }
    ierr = VecRestoreArray(opflow->G, &g);
    CHKERRQ(ierr);
    h_allocator_.deallocate(gsol_host);
  }
  PetscFunctionReturn(0);
}

#endif
#endif
