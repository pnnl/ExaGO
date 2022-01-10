#include "genramp.h"
#include <exago_config.h>
#include <private/opflowimpl.h>
#include <private/scopflowimpl.h>

PetscErrorCode SCOPFLOWModelDestroy_GENRAMP(SCOPFLOW scopflow) {
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = PetscFree(scopflow->model);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSetVariableBounds_GENRAMP(SCOPFLOW scopflow, Vec Xl,
                                                 Vec Xu) {
  PetscErrorCode ierr;
  OPFLOW opflow;
  PetscScalar *xl, *xu, *xli, *xui;
  PetscInt i;
  PetscFunctionBegin;

  ierr = VecGetArray(Xl, &xl);
  CHKERRQ(ierr);
  ierr = VecGetArray(Xu, &xu);
  CHKERRQ(ierr);

  for (i = 0; i < scopflow->nc; i++) {
    opflow = scopflow->opflows[i];

    /* Set bounds on variables */
    xli = xl + scopflow->xstarti[i];
    xui = xu + scopflow->xstarti[i];

    ierr = VecPlaceArray(opflow->Xl, xli);
    CHKERRQ(ierr);
    ierr = VecPlaceArray(opflow->Xu, xui);
    CHKERRQ(ierr);

    /* Set bounds */
    ierr = OPFLOWComputeVariableBounds(opflow, opflow->Xl, opflow->Xu);
    CHKERRQ(ierr);

    ierr = VecResetArray(opflow->Xl);
    CHKERRQ(ierr);
    ierr = VecResetArray(opflow->Xu);
    CHKERRQ(ierr);
  }

  ierr = VecRestoreArray(Xl, &xl);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(Xu, &xu);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSetConstraintBounds_GENRAMP(SCOPFLOW scopflow, Vec Gl,
                                                   Vec Gu) {
  PetscInt i, j, k, ctr;
  PetscErrorCode ierr;
  PetscScalar *gl, *gu, *gli, *gui;
  OPFLOW opflow, opflow0;
  PS ps, ps0;
  PSBUS bus, bus0;
  PSGEN gen, gen0;

  PetscFunctionBegin;

  ierr = VecGetArray(Gl, &gl);
  CHKERRQ(ierr);
  ierr = VecGetArray(Gu, &gu);
  CHKERRQ(ierr);

  opflow0 = scopflow->opflows[0];
  for (i = 0; i < scopflow->nc; i++) {
    opflow = scopflow->opflows[i];

    /* Set bounds on constraints */
    gli = gl + scopflow->gstarti[i];
    gui = gu + scopflow->gstarti[i];

    ierr = VecPlaceArray(opflow->Gl, gli);
    CHKERRQ(ierr);
    ierr = VecPlaceArray(opflow->Gu, gui);
    CHKERRQ(ierr);

    ierr =
        (*opflow->modelops.setconstraintbounds)(opflow, opflow->Gl, opflow->Gu);
    CHKERRQ(ierr);

    ierr = VecResetArray(opflow->Gl);
    CHKERRQ(ierr);
    ierr = VecResetArray(opflow->Gu);
    CHKERRQ(ierr);

    if (scopflow->nconineqcoup[i]) {
      ctr = 0;
      ps = opflow->ps;
      ps0 = opflow0->ps;
      /* Bounds on inequality coupling constraints */
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
              gli[opflow->ncon + ctr] = -10000;
              gui[opflow->ncon + ctr] = 10000;
            } else {
              gli[opflow->ncon + ctr] = 0.0;
              gui[opflow->ncon + ctr] = 0.0;
            }
          } else {
            gli[opflow->ncon + ctr] = -gen->ramp_rate_30min;
            gui[opflow->ncon + ctr] = gen->ramp_rate_30min;
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

PetscErrorCode SCOPFLOWSetVariableandConstraintBounds_GENRAMP(SCOPFLOW scopflow,
                                                              Vec Xl, Vec Xu,
                                                              Vec Gl, Vec Gu) {
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = SCOPFLOWSetVariableBounds_GENRAMP(scopflow, Xl, Xu);
  CHKERRQ(ierr);
  ierr = SCOPFLOWSetConstraintBounds_GENRAMP(scopflow, Gl, Gu);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWSetInitialGuess_GENRAMP(SCOPFLOW scopflow, Vec X) {
  PetscErrorCode ierr;
  OPFLOW opflow;
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
    opflow = scopflow->opflows[i];
    /* Set initial guess and bounds on variables */
    xi = x + scopflow->xstarti[i];
    xli = xl + scopflow->xstarti[i];
    xui = xu + scopflow->xstarti[i];

    ierr = VecPlaceArray(opflow->X, xi);
    CHKERRQ(ierr);
    ierr = VecPlaceArray(opflow->Xl, xli);
    CHKERRQ(ierr);
    ierr = VecPlaceArray(opflow->Xu, xui);
    CHKERRQ(ierr);

    /* Set initial guess */
    ierr = OPFLOWSetInitialGuess(opflow, opflow->X);
    CHKERRQ(ierr);

    ierr = VecResetArray(opflow->X);
    CHKERRQ(ierr);
    ierr = VecResetArray(opflow->Xl);
    CHKERRQ(ierr);
    ierr = VecResetArray(opflow->Xu);
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

PetscErrorCode SCOPFLOWComputeJacobian_GENRAMP(SCOPFLOW scopflow, Vec X,
                                               Mat J) {
  PetscErrorCode ierr;
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
  PetscInt row, col, gloc;
  PetscScalar val;
  PetscInt flps = 0;

  PetscFunctionBegin;

  ierr = VecGetArray(X, &x);
  CHKERRQ(ierr);

  opflow0 = scopflow->opflows[0];
  for (i = 0; i < scopflow->nc; i++) {
    opflow = scopflow->opflows[i];

    roffset = scopflow->gstarti[i];
    coffset = scopflow->xstarti[i];

    xi = x + scopflow->xstarti[i];
    ierr = VecPlaceArray(opflow->X, xi);
    CHKERRQ(ierr);

    ierr = (*opflow->modelops.computeequalityconstraintjacobian)(
        opflow, opflow->X, opflow->Jac_Ge);
    CHKERRQ(ierr);

    ierr = MatGetSize(opflow->Jac_Ge, &nrow, &ncol);
    CHKERRQ(ierr);
    /* Copy over locations and values to Jac */
    for (j = 0; j < nrow; j++) {
      ierr = MatGetRow(opflow->Jac_Ge, j, &nvals, &cols, &vals);
      CHKERRQ(ierr);
      for (k = 0; k < nvals; k++) {
        row = roffset + j;
        col = coffset + cols[k];
        val = vals[k];
        ierr = MatSetValues(J, 1, &row, 1, &col, &val, INSERT_VALUES);
      }
      ierr = MatRestoreRow(opflow->Jac_Ge, j, &nvals, &cols, &vals);
      CHKERRQ(ierr);
    }

    if (opflow->has_gensetpoint) {
      ps = opflow->ps;
      ps0 = opflow0->ps;
      for (j = 0; j < ps->nbus; j++) {
        bus = &ps->bus[j];
        bus0 = &ps0->bus[j];

        for (k = 0; k < bus->ngen; k++) {
          ierr = PSBUSGetGen(bus, k, &gen);
          CHKERRQ(ierr);
          ierr = PSBUSGetGen(bus0, k, &gen0);
          CHKERRQ(ierr);

          if (!gen->status)
            continue;

          loc0 = gen0->startxpowloc;
          gloc = gen->starteqloc + 1;

          x0loc = scopflow->xstarti[0] + loc0;
          row = roffset + gloc;
          col = x0loc;
          val = -1.;
          ierr = MatSetValues(J, 1, &row, 1, &col, &val, INSERT_VALUES);
          CHKERRQ(ierr);
        }
      }
    }

    roffset += opflow->nconeq;
    if (opflow->Nconineq) {
      /* Inequality constrained Jacobian */
      ierr = (*opflow->modelops.computeinequalityconstraintjacobian)(
          opflow, opflow->X, opflow->Jac_Gi);
      CHKERRQ(ierr);

      ierr = MatGetSize(opflow->Jac_Gi, &nrow, &ncol);
      CHKERRQ(ierr);
      /* Copy over locations to triplet format */
      for (j = 0; j < nrow; j++) {
        ierr = MatGetRow(opflow->Jac_Gi, j, &nvals, &cols, &vals);
        CHKERRQ(ierr);
        for (k = 0; k < nvals; k++) {
          row = roffset + j;
          col = coffset + cols[k];
          val = vals[k];
          ierr = MatSetValues(J, 1, &row, 1, &col, &val, INSERT_VALUES);
        }
        ierr = MatRestoreRow(opflow->Jac_Gi, j, &nvals, &cols, &vals);
        CHKERRQ(ierr);
      }

      roffset += opflow->nconineq;
    }

    ierr = VecResetArray(opflow->X);
    CHKERRQ(ierr);

    if (scopflow->nconineqcoup[i]) {
      ps = opflow->ps;
      ps0 = opflow0->ps;
      for (j = 0; j < ps->nbus; j++) {
        bus = &ps->bus[j];
        bus0 = &ps0->bus[j];

        for (k = 0; k < bus->ngen; k++) {
          ierr = PSBUSGetGen(bus, k, &gen);
          CHKERRQ(ierr);
          ierr = PSBUSGetGen(bus0, k, &gen0);
          CHKERRQ(ierr);

          if (!gen->status) {
            if (gen0->status)
              loc0 = gen0->startxpowloc;
            continue;
          } else {
            loc = gen->startxpowloc;
            if (!gen0->status)
              continue;
            loc0 = gen0->startxpowloc;
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

PetscErrorCode SCOPFLOWComputeConstraints_GENRAMP(SCOPFLOW scopflow, Vec X,
                                                  Vec G) {
  PetscErrorCode ierr;
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

  opflow0 = scopflow->opflows[0];
  ps0 = opflow0->ps;
  x0 = x + scopflow->xstarti[0];

  for (i = 0; i < scopflow->nc; i++) {
    xi = x + scopflow->xstarti[i];
    gi = g + scopflow->gstarti[i];

    opflow = scopflow->opflows[i];

    if (opflow->has_gensetpoint) {
      ps = opflow->ps;
      for (j = 0; j < ps->nbus; j++) {
        bus = &ps->bus[j];
        bus0 = &ps0->bus[j];

        for (k = 0; k < bus->ngen; k++) {
          ierr = PSBUSGetGen(bus, k, &gen);
          CHKERRQ(ierr);
          ierr = PSBUSGetGen(bus0, k, &gen0);
          CHKERRQ(ierr);
          if (!gen->status)
            continue;
          /* Update the generator set-point */
          gen->pgs = x0[gen0->startxpowloc];
        }
      }
    }

    ierr = VecPlaceArray(opflow->X, xi);
    CHKERRQ(ierr);

    /* Equality constraints */
    ierr = VecPlaceArray(opflow->Ge, gi);
    CHKERRQ(ierr);
    ierr = (*opflow->modelops.computeequalityconstraints)(opflow, opflow->X,
                                                          opflow->Ge);
    CHKERRQ(ierr);
    ierr = VecResetArray(opflow->Ge);
    CHKERRQ(ierr);
    gi = gi + opflow->nconeq;

    if (opflow->Nconineq) {
      /* Inequality constraints */
      ierr = VecPlaceArray(opflow->Gi, gi);
      CHKERRQ(ierr);
      ierr = (*opflow->modelops.computeinequalityconstraints)(opflow, opflow->X,
                                                              opflow->Gi);
      CHKERRQ(ierr);
      ierr = VecResetArray(opflow->Gi);
      CHKERRQ(ierr);
      gi = gi + opflow->nconineq;
    }

    if (scopflow->nconineqcoup[i]) {
      ctr = 0;
      ps = opflow->ps;
      ps0 = opflow0->ps;
      for (j = 0; j < ps->nbus; j++) {
        bus = &ps->bus[j];
        bus0 = &ps0->bus[j];

        for (k = 0; k < bus->ngen; k++) {
          ierr = PSBUSGetGen(bus, k, &gen);
          CHKERRQ(ierr);
          ierr = PSBUSGetGen(bus0, k, &gen0);
          CHKERRQ(ierr);

          if (!gen->status) {
            if (gen0->status)
              loc0 = gen0->startxpowloc;
            continue;
          } else {
            loc = gen->startxpowloc;
            if (!gen0->status)
              continue;
            loc0 = gen0->startxpowloc;
          }

          gi[ctr] = xi[loc] - x0[loc0]; /* PG(i) - PG(0) */
          ctr++;
        }
      }
    }

    ierr = VecResetArray(opflow->X);
    CHKERRQ(ierr);
  }

  ierr = VecRestoreArray(X, &x);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(G, &g);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWComputeTotalObjective_GENRAMP(SCOPFLOW scopflow, Vec X,
                                                     PetscScalar *obj) {
  PetscErrorCode ierr;
  OPFLOW opflow;
  PetscInt i;
  PetscScalar *xi;
  PetscScalar opflowobj = 0.0;
  PetscScalar *x;

  PetscFunctionBegin;
  *obj = 0.0;

  ierr = VecGetArray(X, &x);
  CHKERRQ(ierr);
  for (i = 0; i < scopflow->nc; i++) {
    opflowobj = 0.0;
    xi = x + scopflow->xstarti[i];
    opflow = scopflow->opflows[i];
    ierr = VecPlaceArray(opflow->X, xi);
    CHKERRQ(ierr);
    ierr = (*opflow->modelops.computeobjective)(opflow, opflow->X, &opflowobj);
    CHKERRQ(ierr);
    *obj += opflowobj;
    ierr = VecResetArray(opflow->X);
    CHKERRQ(ierr);
  }
  ierr = VecRestoreArray(X, &x);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWComputeBaseObjective_GENRAMP(SCOPFLOW scopflow, Vec X,
                                                    PetscScalar *obj) {
  PetscErrorCode ierr;
  OPFLOW opflow;
  PetscScalar *xi;
  PetscScalar opflowobj = 0.0;
  PetscScalar *x;

  PetscFunctionBegin;
  *obj = 0.0;

  ierr = VecGetArray(X, &x);
  CHKERRQ(ierr);

  xi = x + scopflow->xstarti[0];
  opflow = scopflow->opflows[0];
  ierr = VecPlaceArray(opflow->X, xi);
  CHKERRQ(ierr);
  ierr = (*opflow->modelops.computeobjective)(opflow, opflow->X, &opflowobj);
  CHKERRQ(ierr);
  *obj = opflowobj;
  ierr = VecResetArray(opflow->X);
  CHKERRQ(ierr);

  ierr = VecRestoreArray(X, &x);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWComputeGradient_GENRAMP(SCOPFLOW scopflow, Vec X,
                                               Vec Grad) {
  PetscErrorCode ierr;
  OPFLOW opflow;
  PetscInt i;
  PetscScalar *x, *xi, *grad, *gradi;

  ierr = VecGetArray(X, &x);
  CHKERRQ(ierr);
  ierr = VecGetArray(Grad, &grad);
  CHKERRQ(ierr);

  for (i = 0; i < scopflow->nc; i++) {
    opflow = scopflow->opflows[i];
    xi = x + scopflow->xstarti[i];
    gradi = grad + scopflow->xstarti[i];

    ierr = VecPlaceArray(opflow->X, xi);
    CHKERRQ(ierr);
    ierr = VecPlaceArray(opflow->gradobj, gradi);
    CHKERRQ(ierr);
    ierr = VecSet(opflow->gradobj, 0.0);
    CHKERRQ(ierr);

    ierr =
        (*opflow->modelops.computegradient)(opflow, opflow->X, opflow->gradobj);
    CHKERRQ(ierr);

    ierr = VecResetArray(opflow->X);
    CHKERRQ(ierr);
    ierr = VecResetArray(opflow->gradobj);
    CHKERRQ(ierr);
  }

  ierr = VecRestoreArray(X, &x);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(Grad, &grad);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWModelSetNumVariablesandConstraints_GENRAMP(
    SCOPFLOW scopflow, PetscInt *nxi, PetscInt *ngi, PetscInt *nconeqcoup,
    PetscInt *nconineqcoup) {
  PetscInt i, ngenON;
  OPFLOW opflow;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  for (i = 0; i < scopflow->nc; i++) {
    opflow = scopflow->opflows[i];
    ierr = PSGetNumActiveGenerators(opflow->ps, &ngenON, NULL);
    CHKERRQ(ierr);
    nxi[i] = opflow->nx;
    if (scopflow->iscoupling) {
      if (opflow->has_gensetpoint)
        nconeqcoup[i] = 0; //(i == 0)?0:ngenON;
      else
        nconineqcoup[i] = (i == 0) ? 0 : ngenON;
    }
    ngi[i] = opflow->ncon + nconeqcoup[i] + nconineqcoup[i];
  }
  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWComputeHessian_GENRAMP(SCOPFLOW scopflow, Vec X,
                                              Vec Lambda, Mat H) {
  PetscErrorCode ierr;
  PetscInt nrow;
  OPFLOW opflow;
  PetscScalar *x, *xi, *lambda, *lameqi, *lamineqi;
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
    opflow = scopflow->opflows[i];
    opflow->obj_factor = scopflow->obj_factor;

    roffset = scopflow->xstarti[i];

    xi = x + roffset;
    lameqi = lambda + scopflow->gstarti[i];

    ierr = VecPlaceArray(opflow->X, xi);
    CHKERRQ(ierr);
    ierr = VecPlaceArray(opflow->Lambdae, lameqi);
    CHKERRQ(ierr);
    if (opflow->Nconineq) {
      lamineqi = lameqi + opflow->nconeq;
      ierr = VecPlaceArray(opflow->Lambdai, lamineqi);
      CHKERRQ(ierr);
    }

    ierr = (*opflow->modelops.computehessian)(
        opflow, opflow->X, opflow->Lambdae, opflow->Lambdai, opflow->Hes);
    CHKERRQ(ierr);

    ierr = VecResetArray(opflow->X);
    CHKERRQ(ierr);
    ierr = VecResetArray(opflow->Lambdae);
    CHKERRQ(ierr);
    if (opflow->Nconineq) {
      ierr = VecResetArray(opflow->Lambdai);
      CHKERRQ(ierr);
    }

    /* Copy over values */
    ierr = MatGetSize(opflow->Hes, &nrow, &nrow);
    CHKERRQ(ierr);
    for (j = 0; j < nrow; j++) {
      ierr = MatGetRow(opflow->Hes, j, &nvals, &cols, &vals);
      CHKERRQ(ierr);
      row = roffset + j;
      for (k = 0; k < nvals; k++) {
        col = roffset + cols[k];
        val = vals[k];
        ierr = MatSetValues(H, 1, &row, 1, &col, &val, INSERT_VALUES);
        CHKERRQ(ierr);
      }
      ierr = MatRestoreRow(opflow->Hes, j, &nvals, &cols, &vals);
      CHKERRQ(ierr);
    }
  }

  ierr = MatAssemblyBegin(H, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = MatAssemblyEnd(H, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode SCOPFLOWModelCreate_GENRAMP(SCOPFLOW scopflow) {
  GENRAMP genramp;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = PetscCalloc1(1, &genramp);
  CHKERRQ(ierr);

  scopflow->model = genramp;

  /* Inherit Ops */
  scopflow->modelops.destroy = SCOPFLOWModelDestroy_GENRAMP;
  scopflow->modelops.setnumvariablesandconstraints =
      SCOPFLOWModelSetNumVariablesandConstraints_GENRAMP;
  scopflow->modelops.setvariablebounds = SCOPFLOWSetVariableBounds_GENRAMP;
  scopflow->modelops.setconstraintbounds = SCOPFLOWSetConstraintBounds_GENRAMP;
  scopflow->modelops.setvariableandconstraintbounds =
      SCOPFLOWSetVariableandConstraintBounds_GENRAMP;
  scopflow->modelops.setinitialguess = SCOPFLOWSetInitialGuess_GENRAMP;
  scopflow->modelops.computeconstraints = SCOPFLOWComputeConstraints_GENRAMP;
  scopflow->modelops.computejacobian = SCOPFLOWComputeJacobian_GENRAMP;
  scopflow->modelops.computehessian = SCOPFLOWComputeHessian_GENRAMP;
  scopflow->modelops.computebaseobjective =
      SCOPFLOWComputeBaseObjective_GENRAMP;
  scopflow->modelops.computetotalobjective =
      SCOPFLOWComputeTotalObjective_GENRAMP;

  scopflow->modelops.computegradient = SCOPFLOWComputeGradient_GENRAMP;

  PetscFunctionReturn(0);
}
