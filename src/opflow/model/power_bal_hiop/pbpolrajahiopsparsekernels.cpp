
#include <iostream>
#include <iomanip>

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
static const bool oldhostway(true);

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

// A routine to get the triplet arrays from the device and print them out
static void
PrintTriplets(const std::string& title, const int& n, int *i, int *j, double *v)
{
  auto &resmgr = umpire::ResourceManager::getInstance();
  umpire::Allocator h_allocator_ = resmgr.getAllocator("HOST");

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
      
      if (debugmsg)
        PrintTriplets("Nonzero indexes for Inequality Constraint Jacobian (GPU):",
                      opflow->nnz_ineqjacsp,
                      iJacS_dev + opflow->nnz_eqjacsp,
                      jJacS_dev + opflow->nnz_eqjacsp,
                      NULL);

      if (oldhostway) {
        
      /* Inequality constraints start after equality constraints
         Hence the offset
      */
      roffset = opflow->nconeq;
      coffset = 0;
      
      // Create arrays on host to store i,j, and val arrays
      umpire::Allocator h_allocator_ = resmgr.getAllocator("HOST");

      pbpolrajahiopsparse->i_jacineq =
          (int *)(h_allocator_.allocate(opflow->nnz_ineqjacsp * sizeof(int)));
      pbpolrajahiopsparse->j_jacineq =
          (int *)(h_allocator_.allocate(opflow->nnz_ineqjacsp * sizeof(int)));
      pbpolrajahiopsparse->val_jacineq = (double *)(h_allocator_.allocate(
          opflow->nnz_ineqjacsp * sizeof(double)));

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

              // RAJA::atomicAdd<exago_raja_atomic>
              //   (&(MJacS_dev[jac_ieq_idx[i] + 0]), val[0]);
              // RAJA::atomicAdd<exago_raja_atomic>
              //   (&(MJacS_dev[jac_ieq_idx[i] + 1]), val[1]);
              // RAJA::atomicAdd<exago_raja_atomic>
              //   (&(MJacS_dev[jac_ieq_idx[i] + 2]), val[2]);
              // RAJA::atomicAdd<exago_raja_atomic>
              //   (&(MJacS_dev[jac_ieq_idx[i] + 3]), val[3]);

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

              // RAJA::atomicAdd<exago_raja_atomic>
              //   (&(MJacS_dev[jac_ieq_idx[i] + 4]), val[0]);
              // RAJA::atomicAdd<exago_raja_atomic>
              //   (&(MJacS_dev[jac_ieq_idx[i] + 5]), val[1]);
              // RAJA::atomicAdd<exago_raja_atomic>
              //   (&(MJacS_dev[jac_ieq_idx[i] + 6]), val[2]);
              // RAJA::atomicAdd<exago_raja_atomic>
              //   (&(MJacS_dev[jac_ieq_idx[i] + 7]), val[3]);
            });

      }

      if (debugmsg)
        PrintTriplets("Inequality Constraint Jacobian (GPU):",
                      opflow->nnz_ineqjacsp,
                      (iJacS_dev == NULL ? NULL : iJacS_dev + opflow->nnz_eqjacsp),
                      (jJacS_dev == NULL ? NULL : jJacS_dev + opflow->nnz_eqjacsp),
                      MJacS_dev + opflow->nnz_eqjacsp);
      
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

    if (debugmsg)
      PrintTriplets("Non-zero indexes for Equality Constraint Jacobian (GPU):",
                    opflow->nnz_eqjacsp, iJacS_dev, jJacS_dev, NULL);

    if (oldhostway) {
    roffset = 0;
    coffset = 0;

    pbpolrajahiopsparse->i_jaceq =
        (int *)(h_allocator_.allocate(opflow->nnz_eqjacsp * sizeof(int)));
    pbpolrajahiopsparse->j_jaceq =
        (int *)(h_allocator_.allocate(opflow->nnz_eqjacsp * sizeof(int)));
    pbpolrajahiopsparse->val_jaceq =
        (double *)(h_allocator_.allocate(opflow->nnz_eqjacsp * sizeof(double)));

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
                    opflow->nnz_eqjacsp, iJacS_dev, jJacS_dev, NULL);
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

    if (debugmsg)
      PrintTriplets("Equality Constraint Jacobian (GPU):",
                    opflow->nnz_eqjacsp, iJacS_dev, jJacS_dev, MJacS_dev);
      
    
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
                    opflow->nnz_eqjacsp, iJacS_dev, jJacS_dev, MJacS_dev);
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
  PetscErrorCode ierr;
  PetscInt *iRow, *jCol;
  PetscScalar *x, *values, *lambda;
  PetscInt nrow;
  PetscInt nvals;
  const PetscInt *cols;
  const PetscScalar *vals;
  PetscInt i, j;
  PetscInt ctr = 0;
  auto &resmgr = umpire::ResourceManager::getInstance();

  PetscFunctionBegin;

  if (iHSS_dev != NULL && jHSS_dev != NULL) {

    // Create arrays on host to store i,j, and val arrays
    umpire::Allocator h_allocator_ = resmgr.getAllocator("HOST");

    pbpolrajahiopsparse->i_hess =
        (int *)(h_allocator_.allocate(opflow->nnz_hesssp * sizeof(int)));
    pbpolrajahiopsparse->j_hess =
        (int *)(h_allocator_.allocate(opflow->nnz_hesssp * sizeof(int)));
    pbpolrajahiopsparse->val_hess =
        (double *)(h_allocator_.allocate(opflow->nnz_hesssp * sizeof(double)));

    iRow = pbpolrajahiopsparse->i_hess;
    jCol = pbpolrajahiopsparse->j_hess;

    ierr = (*opflow->modelops.computehessian)(
        opflow, opflow->X, opflow->Lambdae, opflow->Lambdai, opflow->Hes);
    CHKERRQ(ierr);
    ierr = MatGetSize(opflow->Hes, &nrow, &nrow);
    CHKERRQ(ierr);

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
  } else {

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
