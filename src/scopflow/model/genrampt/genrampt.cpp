
#include "genrampt.h"
#include <private/opflowimpl.h>
#include <private/scopflowimpl.h>
#include <private/tcopflowimpl.h>

PetscErrorCode SCOPFLOWModelDestroy_GENRAMPT(SCOPFLOW scopflow) {
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = PetscFree(scopflow->model);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSetVariableBounds_GENRAMPT(SCOPFLOW scopflow, Vec Xl,
                                                  Vec Xu) {
  PetscErrorCode ierr;
  TCOPFLOW tcopflow;
  PetscScalar *xl, *xu, *xli, *xui;
  PetscInt i;
  PetscFunctionBegin;

  ierr = VecGetArray(Xl, &xl);
  CHKERRQ(ierr);
  ierr = VecGetArray(Xu, &xu);
  CHKERRQ(ierr);

  for (i = 0; i < scopflow->nc; i++) {
    tcopflow = scopflow->tcopflows[i];

    /* Set bounds on variables */
    xli = xl + scopflow->xstarti[i];
    xui = xu + scopflow->xstarti[i];

    ierr = VecPlaceArray(tcopflow->Xl, xli);
    CHKERRQ(ierr);
    ierr = VecPlaceArray(tcopflow->Xu, xui);
    CHKERRQ(ierr);

    /* Set bounds */
    ierr = (*tcopflow->modelops.setvariablebounds)(tcopflow, tcopflow->Xl,
                                                   tcopflow->Xu);
    CHKERRQ(ierr);

    ierr = VecResetArray(tcopflow->Xl);
    CHKERRQ(ierr);
    ierr = VecResetArray(tcopflow->Xu);
    CHKERRQ(ierr);
  }

  ierr = VecRestoreArray(Xl, &xl);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(Xu, &xu);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSetConstraintBounds_GENRAMPT(SCOPFLOW scopflow, Vec Gl,
                                                    Vec Gu) {
  PetscInt i, j, k, ctr;
  PetscErrorCode ierr;
  PetscScalar *gl, *gu, *gli, *gui;
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

  tcopflow0 = scopflow->tcopflows[0];
  opflow0 = tcopflow0->opflows[0];
  for (i = 0; i < scopflow->nc; i++) {
    tcopflow = scopflow->tcopflows[i];

    /* Set bounds on constraints */
    gli = gl + scopflow->gstarti[i];
    gui = gu + scopflow->gstarti[i];

    ierr = VecPlaceArray(tcopflow->Gl, gli);
    CHKERRQ(ierr);
    ierr = VecPlaceArray(tcopflow->Gu, gui);
    CHKERRQ(ierr);

    ierr = (*tcopflow->modelops.setconstraintbounds)(tcopflow, tcopflow->Gl,
                                                     tcopflow->Gu);
    CHKERRQ(ierr);

    ierr = VecResetArray(tcopflow->Gl);
    CHKERRQ(ierr);
    ierr = VecResetArray(tcopflow->Gu);
    CHKERRQ(ierr);

    gli += tcopflow->Ncon;
    gui += tcopflow->Ncon;

    if (scopflow->nconineqcoup[i]) {
      ctr = 0;
      opflow = tcopflow->opflows[0];
      ps = opflow->ps;
      ps0 = opflow0->ps;
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

          /* Ramp constraints */
          if (scopflow->mode == 0) {
            /* Only ref. bus responsible for make-up power for contingencies */
            if (bus->ide == REF_BUS) {
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

PetscErrorCode
SCOPFLOWSetVariableandConstraintBounds_GENRAMPT(SCOPFLOW scopflow, Vec Xl,
                                                Vec Xu, Vec Gl, Vec Gu) {
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = SCOPFLOWSetVariableBounds_GENRAMPT(scopflow, Xl, Xu);
  CHKERRQ(ierr);
  ierr = SCOPFLOWSetConstraintBounds_GENRAMPT(scopflow, Gl, Gu);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSetInitialGuess_GENRAMPT(SCOPFLOW scopflow, Vec X) {
  PetscErrorCode ierr;
  TCOPFLOW tcopflow;
  PetscScalar *x, *xi, *xl, *xu, *xli, *xui;
  PetscInt i;
  PetscFunctionBegin;

  ierr = VecGetArray(X, &x);
  CHKERRQ(ierr);
  ierr = VecGetArray(scopflow->Xl, &xl);
  CHKERRQ(ierr);
  ierr = VecGetArray(scopflow->Xu, &xu);
  CHKERRQ(ierr);

  for (i = 0; i < scopflow->nc; i++) {
    tcopflow = scopflow->tcopflows[i];
    /* Set initial guess and bounds on variables */
    xi = x + scopflow->xstarti[i];
    xli = xl + scopflow->xstarti[i];
    xui = xu + scopflow->xstarti[i];

    ierr = VecPlaceArray(tcopflow->X, xi);
    CHKERRQ(ierr);
    ierr = VecPlaceArray(tcopflow->Xl, xli);
    CHKERRQ(ierr);
    ierr = VecPlaceArray(tcopflow->Xu, xui);
    CHKERRQ(ierr);

    /* Set initial guess */
    ierr = (*tcopflow->modelops.setinitialguess)(tcopflow, tcopflow->X);
    CHKERRQ(ierr);

    ierr = VecResetArray(tcopflow->X);
    CHKERRQ(ierr);
    ierr = VecResetArray(tcopflow->Xl);
    CHKERRQ(ierr);
    ierr = VecResetArray(tcopflow->Xu);
    CHKERRQ(ierr);
  }

  ierr = VecRestoreArray(X, &x);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(scopflow->Xl, &xl);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(scopflow->Xu, &xu);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWComputeJacobian_GENRAMPT(SCOPFLOW scopflow, Vec X,
                                                Mat J) {
  PetscErrorCode ierr;
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

  tcopflow0 = scopflow->tcopflows[0];
  opflow0 = tcopflow0->opflows[0];
  for (i = 0; i < scopflow->nc; i++) {
    tcopflow = scopflow->tcopflows[i];
    opflow = tcopflow->opflows[0];

    roffset = scopflow->gstarti[i];
    coffset = scopflow->xstarti[i];

    xi = x + scopflow->xstarti[i];

    ierr = VecPlaceArray(tcopflow->X, xi);
    CHKERRQ(ierr);

    ierr = (*tcopflow->modelops.computejacobian)(tcopflow, tcopflow->X,
                                                 tcopflow->Jac);
    CHKERRQ(ierr);

    ierr = MatGetSize(tcopflow->Jac, &nrow, &ncol);
    CHKERRQ(ierr);
    /* Copy over locations and values to Jac */
    for (j = 0; j < nrow; j++) {
      ierr = MatGetRow(tcopflow->Jac, j, &nvals, &cols, &vals);
      CHKERRQ(ierr);
      for (k = 0; k < nvals; k++) {
        row = roffset + j;
        col = coffset + cols[k];
        val = vals[k];
        ierr = MatSetValues(J, 1, &row, 1, &col, &val, INSERT_VALUES);
      }
      ierr = MatRestoreRow(tcopflow->Jac, j, &nvals, &cols, &vals);
      CHKERRQ(ierr);
    }

    ierr = VecResetArray(tcopflow->X);
    CHKERRQ(ierr);

    roffset += tcopflow->Ncon;

    if (scopflow->nconineqcoup[i]) {
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

          x0loc = scopflow->xstarti[0] + loc0;
          xiloc = scopflow->xstarti[i] + loc;
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

PetscErrorCode SCOPFLOWComputeConstraints_GENRAMPT(SCOPFLOW scopflow, Vec X,
                                                   Vec G) {
  PetscErrorCode ierr;
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

  tcopflow0 = scopflow->tcopflows[0];
  opflow0 = tcopflow0->opflows[0];
  x0 = x + scopflow->xstarti[0];

  for (i = 0; i < scopflow->nc; i++) {
    xi = x + scopflow->xstarti[i];
    gi = g + scopflow->gstarti[i];

    tcopflow = scopflow->tcopflows[i];
    opflow = tcopflow->opflows[0];

    ierr = VecPlaceArray(tcopflow->X, xi);
    CHKERRQ(ierr);
    ierr = VecPlaceArray(tcopflow->G, gi);
    CHKERRQ(ierr);

    ierr = (*tcopflow->modelops.computeconstraints)(tcopflow, tcopflow->X,
                                                    tcopflow->G);
    CHKERRQ(ierr);

    ierr = VecResetArray(tcopflow->X);
    CHKERRQ(ierr);
    ierr = VecResetArray(tcopflow->G);
    CHKERRQ(ierr);

    gi += tcopflow->Ncon;

    if (scopflow->nconineqcoup[i]) {
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

PetscErrorCode SCOPFLOWComputeTotalObjective_GENRAMPT(SCOPFLOW scopflow, Vec X,
                                                      PetscScalar *obj) {
  PetscErrorCode ierr;
  TCOPFLOW tcopflow;
  PetscInt i;
  PetscScalar *xi;
  PetscScalar tcopflowobj = 0.0;
  PetscScalar *x;

  PetscFunctionBegin;
  *obj = 0.0;

  ierr = VecGetArray(X, &x);
  CHKERRQ(ierr);
  for (i = 0; i < scopflow->nc; i++) {
    tcopflowobj = 0.0;
    xi = x + scopflow->xstarti[i];
    tcopflow = scopflow->tcopflows[i];
    ierr = VecPlaceArray(tcopflow->X, xi);
    CHKERRQ(ierr);
    ierr = (*tcopflow->modelops.computeobjective)(tcopflow, tcopflow->X,
                                                  &tcopflowobj);
    CHKERRQ(ierr);
    *obj += tcopflowobj;
    ierr = VecResetArray(tcopflow->X);
    CHKERRQ(ierr);
  }
  ierr = VecRestoreArray(X, &x);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWComputeBaseObjective_GENRAMPT(SCOPFLOW scopflow, Vec X,
                                                     PetscScalar *obj) {
  PetscErrorCode ierr;
  TCOPFLOW tcopflow;
  PetscScalar *xi;
  PetscScalar tcopflowobj = 0.0;
  PetscScalar *x;

  PetscFunctionBegin;
  *obj = 0.0;

  ierr = VecGetArray(X, &x);
  CHKERRQ(ierr);
  xi = x + scopflow->xstarti[0];
  tcopflow = scopflow->tcopflows[0];
  ierr = VecPlaceArray(tcopflow->X, xi);
  CHKERRQ(ierr);
  ierr = (*tcopflow->modelops.computeobjective)(tcopflow, tcopflow->X,
                                                &tcopflowobj);
  CHKERRQ(ierr);
  *obj = tcopflowobj;
  ierr = VecResetArray(tcopflow->X);
  CHKERRQ(ierr);

  ierr = VecRestoreArray(X, &x);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWComputeGradient_GENRAMPT(SCOPFLOW scopflow, Vec X,
                                                Vec Grad) {
  PetscErrorCode ierr;
  TCOPFLOW tcopflow;
  PetscInt i;
  PetscScalar *x, *xi, *grad, *gradi;

  ierr = VecGetArray(X, &x);
  CHKERRQ(ierr);
  ierr = VecGetArray(Grad, &grad);
  CHKERRQ(ierr);

  for (i = 0; i < scopflow->nc; i++) {
    tcopflow = scopflow->tcopflows[i];
    xi = x + scopflow->xstarti[i];
    gradi = grad + scopflow->xstarti[i];

    ierr = VecPlaceArray(tcopflow->X, xi);
    CHKERRQ(ierr);
    ierr = VecPlaceArray(tcopflow->gradobj, gradi);
    CHKERRQ(ierr);
    ierr = VecSet(tcopflow->gradobj, 0.0);
    CHKERRQ(ierr);

    ierr = (*tcopflow->modelops.computegradient)(tcopflow, tcopflow->X,
                                                 tcopflow->gradobj);
    CHKERRQ(ierr);

    ierr = VecResetArray(tcopflow->X);
    CHKERRQ(ierr);
    ierr = VecResetArray(tcopflow->gradobj);
    CHKERRQ(ierr);
  }

  ierr = VecRestoreArray(X, &x);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(Grad, &grad);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWModelSetNumVariablesandConstraints_GENRAMPT(
    SCOPFLOW scopflow, PetscInt *nxi, PetscInt *ngi, PetscInt *nconeqcoup,
    PetscInt *nconineqcoup) {
  (void)nconeqcoup;
  PetscInt c, ngenON;
  TCOPFLOW tcopflow;
  OPFLOW opflow;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  for (c = 0; c < scopflow->nc; c++) {
    tcopflow = scopflow->tcopflows[c];
    opflow = tcopflow->opflows[0];
    ierr = PSGetNumActiveGenerators(opflow->ps, &ngenON, NULL);
    CHKERRQ(ierr);

    nxi[c] = tcopflow->Nx;
    if (scopflow->iscoupling)
      nconineqcoup[c] = (c == 0) ? 0 : ngenON;
    ngi[c] = tcopflow->Ncon + nconineqcoup[c];
  }
  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWComputeHessian_GENRAMPT(SCOPFLOW scopflow, Vec X,
                                               Vec Lambda, Mat H) {
  PetscErrorCode ierr;
  PetscInt nrow;
  TCOPFLOW tcopflow;
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
  for (i = 0; i < scopflow->nc; i++) {
    tcopflow = scopflow->tcopflows[i];
    tcopflow->obj_factor = scopflow->obj_factor;

    roffset = scopflow->xstarti[i];

    xi = x + roffset;
    lambdai = lambda + scopflow->gstarti[i];

    ierr = VecPlaceArray(tcopflow->X, xi);
    CHKERRQ(ierr);
    ierr = VecPlaceArray(tcopflow->Lambda, lambdai);
    CHKERRQ(ierr);

    ierr = (*tcopflow->modelops.computehessian)(
        tcopflow, tcopflow->X, tcopflow->Lambda, tcopflow->Hes);
    CHKERRQ(ierr);

    ierr = VecResetArray(tcopflow->X);
    CHKERRQ(ierr);
    ierr = VecResetArray(tcopflow->Lambda);
    CHKERRQ(ierr);

    /* Copy over values */
    ierr = MatGetSize(tcopflow->Hes, &nrow, &nrow);
    CHKERRQ(ierr);
    for (j = 0; j < nrow; j++) {
      ierr = MatGetRow(tcopflow->Hes, j, &nvals, &cols, &vals);
      CHKERRQ(ierr);
      row = roffset + j;
      for (k = 0; k < nvals; k++) {
        col = roffset + cols[k];
        val = vals[k];
        ierr = MatSetValues(H, 1, &row, 1, &col, &val, INSERT_VALUES);
        CHKERRQ(ierr);
      }
      ierr = MatRestoreRow(tcopflow->Hes, j, &nvals, &cols, &vals);
      CHKERRQ(ierr);
    }
  }

  ierr = MatAssemblyBegin(H, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = MatAssemblyEnd(H, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWModelCreate_GENRAMPT(SCOPFLOW scopflow) {
  GENRAMPT genrampt;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = PetscCalloc1(1, &genrampt);
  CHKERRQ(ierr);

  scopflow->model = genrampt;

  /* Inherit Ops */
  scopflow->modelops.destroy = SCOPFLOWModelDestroy_GENRAMPT;
  scopflow->modelops.setnumvariablesandconstraints =
      SCOPFLOWModelSetNumVariablesandConstraints_GENRAMPT;
  scopflow->modelops.setvariablebounds = SCOPFLOWSetVariableBounds_GENRAMPT;
  scopflow->modelops.setconstraintbounds = SCOPFLOWSetConstraintBounds_GENRAMPT;
  scopflow->modelops.setvariableandconstraintbounds =
      SCOPFLOWSetVariableandConstraintBounds_GENRAMPT;
  scopflow->modelops.setinitialguess = SCOPFLOWSetInitialGuess_GENRAMPT;
  scopflow->modelops.computeconstraints = SCOPFLOWComputeConstraints_GENRAMPT;
  scopflow->modelops.computejacobian = SCOPFLOWComputeJacobian_GENRAMPT;
  scopflow->modelops.computehessian = SCOPFLOWComputeHessian_GENRAMPT;
  scopflow->modelops.computebaseobjective =
      SCOPFLOWComputeBaseObjective_GENRAMPT;
  scopflow->modelops.computetotalobjective =
      SCOPFLOWComputeTotalObjective_GENRAMPT;
  scopflow->modelops.computegradient = SCOPFLOWComputeGradient_GENRAMPT;

  PetscFunctionReturn(0);
}
