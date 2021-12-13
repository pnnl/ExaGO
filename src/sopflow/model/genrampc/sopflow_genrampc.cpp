
#include "sopflow_genrampc.h"
#include <private/opflowimpl.h>
#include <private/scopflowimpl.h>
#include <private/sopflowimpl.h>
#include <private/tcopflowimpl.h>

PetscErrorCode SOPFLOWModelDestroy_GENRAMPC(SOPFLOW sopflow) {
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = PetscFree(sopflow->model);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSetVariableBounds_GENRAMPC(SOPFLOW sopflow, Vec Xl,
                                                 Vec Xu) {
  PetscErrorCode ierr;
  SCOPFLOW scopflow;
  PetscScalar *xl, *xu, *xli, *xui;
  PetscInt i;
  PetscFunctionBegin;

  ierr = VecGetArray(Xl, &xl);
  CHKERRQ(ierr);
  ierr = VecGetArray(Xu, &xu);
  CHKERRQ(ierr);

  for (i = 0; i < sopflow->ns; i++) {
    scopflow = sopflow->scopflows[i];

    /* Set bounds on variables */
    xli = xl + sopflow->xstarti[i];
    xui = xu + sopflow->xstarti[i];

    ierr = VecPlaceArray(scopflow->Xl, xli);
    CHKERRQ(ierr);
    ierr = VecPlaceArray(scopflow->Xu, xui);
    CHKERRQ(ierr);

    /* Set bounds */
    ierr = (*scopflow->modelops.setvariablebounds)(scopflow, scopflow->Xl,
                                                   scopflow->Xu);
    CHKERRQ(ierr);

    ierr = VecResetArray(scopflow->Xl);
    CHKERRQ(ierr);
    ierr = VecResetArray(scopflow->Xu);
    CHKERRQ(ierr);
  }

  ierr = VecRestoreArray(Xl, &xl);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(Xu, &xu);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSetConstraintBounds_GENRAMPC(SOPFLOW sopflow, Vec Gl,
                                                   Vec Gu) {
  PetscInt i, j, k, ctr;
  PetscErrorCode ierr;
  PetscScalar *gl, *gu, *gli, *gui;
  SCOPFLOW scopflow, scopflow0;
  TCOPFLOW tcopflow, tcopflow0;
  OPFLOW opflow, opflow0;
  PS ps, ps0;
  PSBUS bus, bus0;
  PSGEN gen, gen0;

  PetscFunctionBegin;

  ierr = VecGetArray(Gl, &gl);
  CHKERRQ(ierr);
  ierr = VecGetArray(Gu, &gu);
  CHKERRQ(ierr);

  scopflow0 = sopflow->scopflows[0];
  if (!scopflow0->ismultiperiod)
    opflow0 = scopflow0->opflows[0];
  else {
    tcopflow0 = scopflow0->tcopflows[0];
    opflow0 = tcopflow0->opflows[0];
  }
  ps0 = opflow0->ps;
  for (i = 0; i < sopflow->ns; i++) {
    scopflow = sopflow->scopflows[i];

    /* Set bounds on constraints */
    gli = gl + sopflow->gstarti[i];
    gui = gu + sopflow->gstarti[i];

    ierr = VecPlaceArray(scopflow->Gl, gli);
    CHKERRQ(ierr);
    ierr = VecPlaceArray(scopflow->Gu, gui);
    CHKERRQ(ierr);

    ierr = (*scopflow->modelops.setconstraintbounds)(scopflow, scopflow->Gl,
                                                     scopflow->Gu);
    CHKERRQ(ierr);

    ierr = VecResetArray(scopflow->Gl);
    CHKERRQ(ierr);
    ierr = VecResetArray(scopflow->Gu);
    CHKERRQ(ierr);

    gli += scopflow->ncon;
    gui += scopflow->ncon;

    if (sopflow->nconineqcoup[i]) {
      ctr = 0;
      if (!scopflow->ismultiperiod)
        opflow = scopflow->opflows[0];
      else {
        tcopflow = scopflow->tcopflows[0];
        opflow = tcopflow->opflows[0];
      }

      ps = opflow->ps;
      /* Bounds on coupling constraints */
      for (j = 0; j < ps->nbus; j++) {
        bus = &ps->bus[j];
        bus0 = &ps0->bus[j];

        for (k = 0; k < bus->ngen; k++) {
          ierr = PSBUSGetGen(bus, k, &gen);
          CHKERRQ(ierr);
          ierr = PSBUSGetGen(bus0, k, &gen0);
          CHKERRQ(ierr);
          if (!gen->status || !gen0->status)
            continue;

          if (sopflow->mode == 0) {
            /* Only ref. bus and renewable generation can deviate from base-case
             * set points */
            if (bus->ide == REF_BUS || gen->genfuel_type == GENFUEL_WIND) {
              gli[ctr] = -10000;
              gui[ctr] = 10000;
            } else {
              gli[ctr] = 0.0;
              gui[ctr] = 0.0;
            }
          } else {
            gli[ctr] = -gen->ramp_rate_30min;
            gui[ctr] = gen->ramp_rate_30min;
          }

          ctr++;
        }
      }
    }
  }

  ierr = VecRestoreArray(Gl, &gl);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(Gu, &gu);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSetVariableandConstraintBounds_GENRAMPC(SOPFLOW sopflow,
                                                              Vec Xl, Vec Xu,
                                                              Vec Gl, Vec Gu) {
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = SOPFLOWSetVariableBounds_GENRAMPC(sopflow, Xl, Xu);
  CHKERRQ(ierr);
  ierr = SOPFLOWSetConstraintBounds_GENRAMPC(sopflow, Gl, Gu);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWSetInitialGuess_GENRAMPC(SOPFLOW sopflow, Vec X) {
  PetscErrorCode ierr;
  SCOPFLOW scopflow;
  OPFLOW opflow;
  PetscScalar *x, *xi, *xl, *xu, *xli, *xui;
  PetscInt i;
  PetscFunctionBegin;

  ierr = VecGetArray(X, &x);
  CHKERRQ(ierr);
  ierr = VecGetArray(sopflow->Xl, &xl);
  CHKERRQ(ierr);
  ierr = VecGetArray(sopflow->Xu, &xu);
  CHKERRQ(ierr);

  for (i = 0; i < sopflow->ns; i++) {
    scopflow = sopflow->scopflows[i];
    /* Set initial guess and bounds on variables */
    xi = x + sopflow->xstarti[i];
    xli = xl + sopflow->xstarti[i];
    xui = xu + sopflow->xstarti[i];

    ierr = VecPlaceArray(scopflow->X, xi);
    CHKERRQ(ierr);
    ierr = VecPlaceArray(scopflow->Xl, xli);
    CHKERRQ(ierr);
    ierr = VecPlaceArray(scopflow->Xu, xui);
    CHKERRQ(ierr);

    /* Set initial guess */
    ierr = (*scopflow->modelops.setinitialguess)(scopflow, scopflow->X);
    CHKERRQ(ierr);

    ierr = VecResetArray(scopflow->X);
    CHKERRQ(ierr);
    ierr = VecResetArray(scopflow->Xl);
    CHKERRQ(ierr);
    ierr = VecResetArray(scopflow->Xu);
    CHKERRQ(ierr);
  }

  ierr = VecRestoreArray(X, &x);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(sopflow->Xl, &xl);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(sopflow->Xu, &xu);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWComputeJacobian_GENRAMPC(SOPFLOW sopflow, Vec X, Mat J) {
  PetscErrorCode ierr;
  SCOPFLOW scopflow, scopflow0;
  TCOPFLOW tcopflow, tcopflow0;
  OPFLOW opflow, opflow0;
  PetscInt roffset, coffset;
  PetscInt nrow, ncol;
  PetscScalar *xi, *x;
  PetscInt i, j, k, loc, loc0, x0loc, xiloc;
  PS ps, ps0;
  PSBUS bus, bus0;
  PSGEN gen, gen0;
  PetscInt nvals;
  const PetscInt *cols;
  const PetscScalar *vals;
  PetscInt row, col;
  PetscScalar val;

  PetscFunctionBegin;

  ierr = VecGetArray(X, &x);
  CHKERRQ(ierr);

  scopflow0 = sopflow->scopflows[0];
  if (!scopflow0->ismultiperiod)
    opflow0 = scopflow0->opflows[0];
  else {
    tcopflow0 = scopflow0->tcopflows[0];
    opflow0 = tcopflow0->opflows[0];
  }

  for (i = 0; i < sopflow->ns; i++) {
    scopflow = sopflow->scopflows[i];
    if (!scopflow->ismultiperiod)
      opflow = scopflow->opflows[0];
    else {
      tcopflow = scopflow->tcopflows[0];
      opflow = tcopflow->opflows[0];
    }

    roffset = sopflow->gstarti[i];
    coffset = sopflow->xstarti[i];

    xi = x + sopflow->xstarti[i];

    ierr = VecPlaceArray(scopflow->X, xi);
    CHKERRQ(ierr);

    ierr = (*scopflow->modelops.computejacobian)(scopflow, scopflow->X,
                                                 scopflow->Jac);
    CHKERRQ(ierr);

    ierr = MatGetSize(scopflow->Jac, &nrow, &ncol);
    CHKERRQ(ierr);
    /* Copy over locations and values to Jac */
    for (j = 0; j < nrow; j++) {
      ierr = MatGetRow(scopflow->Jac, j, &nvals, &cols, &vals);
      CHKERRQ(ierr);
      for (k = 0; k < nvals; k++) {
        row = roffset + j;
        col = coffset + cols[k];
        val = vals[k];
        ierr = MatSetValues(J, 1, &row, 1, &col, &val, INSERT_VALUES);
      }
      ierr = MatRestoreRow(scopflow->Jac, j, &nvals, &cols, &vals);
      CHKERRQ(ierr);
    }

    roffset += scopflow->ncon;

    ierr = VecResetArray(scopflow->X);
    CHKERRQ(ierr);

    if (sopflow->nconineqcoup[i]) {
      ps = opflow->ps;
      ps0 = opflow0->ps;
      for (j = 0; j < ps->nbus; j++) {
        bus = &ps->bus[j];
        bus0 = &ps0->bus[j];
        ierr = PSBUSGetVariableLocation(bus, &loc);
        CHKERRQ(ierr);
        ierr = PSBUSGetVariableLocation(bus0, &loc0);
        CHKERRQ(ierr);
        for (k = 0; k < bus->ngen; k++) {
          ierr = PSBUSGetGen(bus, k, &gen);
          CHKERRQ(ierr);
          ierr = PSBUSGetGen(bus0, k, &gen0);
          CHKERRQ(ierr);

          if (!gen->status) {
            if (gen0->status)
              loc0 += 2;
            continue;
          } else {
            loc += 2;
            if (!gen0->status)
              continue;
            loc0 += 2;
          }

          x0loc = sopflow->xstarti[0] + loc0;
          xiloc = sopflow->xstarti[i] + loc;
          row = roffset;
          col = x0loc;
          val = -1.;
          ierr = MatSetValues(J, 1, &row, 1, &col, &val, INSERT_VALUES);
          CHKERRQ(ierr);
          col = xiloc;
          val = 1.;
          ierr = MatSetValues(J, 1, &row, 1, &col, &val, INSERT_VALUES);
          CHKERRQ(ierr);

          roffset += 1;
        }
      }
    }
  }

  ierr = VecRestoreArray(X, &x);
  CHKERRQ(ierr);

  ierr = MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWComputeConstraints_GENRAMPC(SOPFLOW sopflow, Vec X,
                                                  Vec G) {
  PetscErrorCode ierr;
  SCOPFLOW scopflow0, scopflow;
  TCOPFLOW tcopflow0, tcopflow;
  OPFLOW opflow0, opflow;
  PetscInt i, j, k, loc, loc0, ctr;
  PetscScalar *x0, *x, *xi, *g, *gi;
  PS ps, ps0;
  PSBUS bus, bus0;
  PSGEN gen, gen0;

  PetscFunctionBegin;

  ierr = VecGetArray(X, &x);
  CHKERRQ(ierr);
  ierr = VecGetArray(G, &g);
  CHKERRQ(ierr);

  scopflow0 = sopflow->scopflows[0];
  if (!scopflow0->ismultiperiod)
    opflow0 = scopflow0->opflows[0];
  else {
    tcopflow0 = scopflow0->tcopflows[0];
    opflow0 = tcopflow0->opflows[0];
  }

  x0 = x + sopflow->xstarti[0];

  for (i = 0; i < sopflow->ns; i++) {
    xi = x + sopflow->xstarti[i];
    gi = g + sopflow->gstarti[i];

    scopflow = sopflow->scopflows[i];
    if (!scopflow->ismultiperiod)
      opflow = scopflow->opflows[0];
    else {
      tcopflow = scopflow->tcopflows[0];
      opflow = tcopflow->opflows[0];
    }

    ierr = VecPlaceArray(scopflow->X, xi);
    CHKERRQ(ierr);
    ierr = VecPlaceArray(scopflow->G, gi);
    CHKERRQ(ierr);

    ierr = (*scopflow->modelops.computeconstraints)(scopflow, scopflow->X,
                                                    scopflow->G);
    CHKERRQ(ierr);

    ierr = VecResetArray(scopflow->X);
    CHKERRQ(ierr);
    ierr = VecResetArray(scopflow->G);
    CHKERRQ(ierr);

    gi += scopflow->ncon;

    if (sopflow->nconineqcoup[i]) {
      ctr = 0;
      ps = opflow->ps;
      ps0 = opflow0->ps;
      for (j = 0; j < ps->nbus; j++) {
        bus = &ps->bus[j];
        bus0 = &ps0->bus[j];
        ierr = PSBUSGetVariableLocation(bus, &loc);
        CHKERRQ(ierr);
        ierr = PSBUSGetVariableLocation(bus0, &loc0);
        CHKERRQ(ierr);

        for (k = 0; k < bus->ngen; k++) {
          ierr = PSBUSGetGen(bus, k, &gen);
          CHKERRQ(ierr);
          ierr = PSBUSGetGen(bus0, k, &gen0);
          CHKERRQ(ierr);

          if (!gen->status) {
            if (gen0->status)
              loc0 += 2;
            continue;
          } else {
            loc += 2;
            if (!gen0->status)
              continue;
            loc0 += 2;
          }

          gi[ctr] = xi[loc] - x0[loc0]; /* PG(i) - PG(0) */
          ctr++;
        }
      }
    }
  }

  ierr = VecRestoreArray(X, &x);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(G, &g);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWComputeObjective_GENRAMPC(SOPFLOW sopflow, Vec X,
                                                PetscScalar *obj) {
  PetscErrorCode ierr;
  SCOPFLOW scopflow;
  PetscInt i;
  PetscScalar *xi;
  PetscScalar scopflowobj = 0.0;
  PetscScalar *x;

  PetscFunctionBegin;
  *obj = 0.0;

  ierr = VecGetArray(X, &x);
  CHKERRQ(ierr);
  for (i = 0; i < sopflow->ns; i++) {
    scopflowobj = 0.0;
    xi = x + sopflow->xstarti[i];
    scopflow = sopflow->scopflows[i];
    ierr = VecPlaceArray(scopflow->X, xi);
    CHKERRQ(ierr);
    ierr = (*scopflow->modelops.computeobjective)(scopflow, scopflow->X,
                                                  &scopflowobj);
    CHKERRQ(ierr);
    *obj += scopflowobj;
    ierr = VecResetArray(scopflow->X);
    CHKERRQ(ierr);
  }
  ierr = VecRestoreArray(X, &x);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWComputeGradient_GENRAMPC(SOPFLOW sopflow, Vec X,
                                               Vec Grad) {
  PetscErrorCode ierr;
  SCOPFLOW scopflow;
  PetscInt i;
  PetscScalar *x, *xi, *grad, *gradi;

  ierr = VecGetArray(X, &x);
  CHKERRQ(ierr);
  ierr = VecGetArray(Grad, &grad);
  CHKERRQ(ierr);

  for (i = 0; i < sopflow->ns; i++) {
    xi = x + sopflow->xstarti[i];
    gradi = grad + sopflow->xstarti[i];

    scopflow = sopflow->scopflows[i];

    ierr = VecPlaceArray(scopflow->X, xi);
    CHKERRQ(ierr);
    ierr = VecPlaceArray(scopflow->gradobj, gradi);
    CHKERRQ(ierr);
    ierr = VecSet(scopflow->gradobj, 0.0);
    CHKERRQ(ierr);

    ierr = (*scopflow->modelops.computegradient)(scopflow, scopflow->X,
                                                 scopflow->gradobj);
    CHKERRQ(ierr);

    ierr = VecResetArray(scopflow->X);
    CHKERRQ(ierr);
    ierr = VecResetArray(scopflow->gradobj);
    CHKERRQ(ierr);
  }

  ierr = VecRestoreArray(X, &x);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(Grad, &grad);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWComputeObjandGradient_GENRAMPC(SOPFLOW sopflow, Vec X,
                                                     PetscScalar *obj,
                                                     Vec Grad) {
  PetscErrorCode ierr;
  PetscFunctionBegin;

  ierr = SOPFLOWComputeObjective_GENRAMPC(sopflow, X, obj);
  CHKERRQ(ierr);
  ierr = SOPFLOWComputeGradient_GENRAMPC(sopflow, X, Grad);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWModelSetNumVariablesandConstraints_GENRAMPC(
    SOPFLOW sopflow, PetscInt *nxi, PetscInt *ngi, PetscInt *nconeqcoup,
    PetscInt *nconineqcoup) {
  PetscInt s, ngenON;
  SCOPFLOW scopflow;
  OPFLOW opflow;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  for (s = 0; s < sopflow->ns; s++) {
    scopflow = sopflow->scopflows[s];
    if (scopflow->ismultiperiod) {
      opflow = scopflow->tcopflows[0]->opflows[0];
    } else {
      opflow = scopflow->opflows[0];
    }
    ierr = PSGetNumActiveGenerators(opflow->ps, &ngenON, NULL);
    CHKERRQ(ierr);

    ierr = SCOPFLOWGetNumVariablesandConstraints(scopflow, &nxi[s], &ngi[s]);
    CHKERRQ(ierr);
    ;
    if (sopflow->iscoupling)
      nconineqcoup[s] = (s == 0) ? 0 : ngenON;
    ngi[s] += nconineqcoup[s];
  }
  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWComputeHessian_GENRAMPC(SOPFLOW sopflow, Vec X,
                                              Vec Lambda, Mat H) {
  PetscErrorCode ierr;
  PetscInt nrow;
  SCOPFLOW scopflow;
  PetscScalar *x, *xi, *lambda, *lambdai;
  PetscInt i;
  PetscInt roffset;
  PetscInt nvals;
  const PetscInt *cols;
  const PetscScalar *vals;
  PetscInt j, k;
  PetscInt row, col;
  PetscScalar val;

  PetscFunctionBegin;

  ierr = MatZeroEntries(H);
  CHKERRQ(ierr);

  ierr = VecGetArray(X, &x);
  CHKERRQ(ierr);
  ierr = VecGetArray(Lambda, &lambda);
  CHKERRQ(ierr);
  for (i = 0; i < sopflow->ns; i++) {
    scopflow = sopflow->scopflows[i];
    scopflow->obj_factor = sopflow->obj_factor;

    roffset = sopflow->xstarti[i];

    xi = x + roffset;
    lambdai = lambda + sopflow->gstarti[i];

    ierr = VecPlaceArray(scopflow->X, xi);
    CHKERRQ(ierr);
    ierr = VecPlaceArray(scopflow->Lambda, lambdai);
    CHKERRQ(ierr);

    ierr = (*scopflow->modelops.computehessian)(
        scopflow, scopflow->X, scopflow->Lambda, scopflow->Hes);
    CHKERRQ(ierr);

    ierr = VecResetArray(scopflow->X);
    CHKERRQ(ierr);
    ierr = VecResetArray(scopflow->Lambda);
    CHKERRQ(ierr);

    /* Copy over values */
    ierr = MatGetSize(scopflow->Hes, &nrow, &nrow);
    CHKERRQ(ierr);
    for (j = 0; j < nrow; j++) {
      ierr = MatGetRow(scopflow->Hes, j, &nvals, &cols, &vals);
      CHKERRQ(ierr);

      row = roffset + j;
      for (k = 0; k < nvals; k++) {
        col = roffset + cols[k];
        val = vals[k];
        ierr = MatSetValues(H, 1, &row, 1, &col, &val, INSERT_VALUES);
        CHKERRQ(ierr);
      }
      ierr = MatRestoreRow(scopflow->Hes, j, &nvals, &cols, &vals);
      CHKERRQ(ierr);
    }
  }

  ierr = MatAssemblyBegin(H, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = MatAssemblyEnd(H, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode SOPFLOWModelCreate_GENRAMPC(SOPFLOW sopflow) {
  GENRAMPC genrampc;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = PetscCalloc1(1, &genrampc);
  CHKERRQ(ierr);

  sopflow->model = genrampc;

  /* Inherit Ops */
  sopflow->modelops.destroy = SOPFLOWModelDestroy_GENRAMPC;
  sopflow->modelops.setnumvariablesandconstraints =
      SOPFLOWModelSetNumVariablesandConstraints_GENRAMPC;
  sopflow->modelops.setvariablebounds = SOPFLOWSetVariableBounds_GENRAMPC;
  sopflow->modelops.setconstraintbounds = SOPFLOWSetConstraintBounds_GENRAMPC;
  sopflow->modelops.setvariableandconstraintbounds =
      SOPFLOWSetVariableandConstraintBounds_GENRAMPC;
  sopflow->modelops.setinitialguess = SOPFLOWSetInitialGuess_GENRAMPC;
  sopflow->modelops.computeconstraints = SOPFLOWComputeConstraints_GENRAMPC;
  sopflow->modelops.computejacobian = SOPFLOWComputeJacobian_GENRAMPC;
  sopflow->modelops.computehessian = SOPFLOWComputeHessian_GENRAMPC;
  sopflow->modelops.computeobjandgradient =
      SOPFLOWComputeObjandGradient_GENRAMPC;
  sopflow->modelops.computeobjective = SOPFLOWComputeObjective_GENRAMPC;
  sopflow->modelops.computegradient = SOPFLOWComputeGradient_GENRAMPC;

  PetscFunctionReturn(0);
}
