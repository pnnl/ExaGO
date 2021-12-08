
#include "genramp.h"
#include <private/opflowimpl.h>
#include <private/tcopflowimpl.h>

PetscErrorCode TCOPFLOWModelDestroy_GENRAMP(TCOPFLOW tcopflow) {
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = PetscFree(tcopflow->model);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode TCOPFLOWSetVariableBounds_GENRAMP(TCOPFLOW tcopflow, Vec Xl,
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

  for (i = 0; i < tcopflow->Nt; i++) {
    opflow = tcopflow->opflows[i];

    /* Set bounds on variables */
    xli = xl + tcopflow->xstarti[i];
    xui = xu + tcopflow->xstarti[i];

    ierr = VecPlaceArray(opflow->Xl, xli);
    CHKERRQ(ierr);
    ierr = VecPlaceArray(opflow->Xu, xui);
    CHKERRQ(ierr);

    /* Set bounds */
    ierr =
        (*opflow->modelops.setvariablebounds)(opflow, opflow->Xl, opflow->Xu);
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

PetscErrorCode TCOPFLOWSetConstraintBounds_GENRAMP(TCOPFLOW tcopflow, Vec Gl,
                                                   Vec Gu) {
  PetscInt i, j, k, ctr;
  PetscErrorCode ierr;
  PetscScalar *gl, *gu, *gli, *gui;
  OPFLOW opflow, opflowpre;
  PS ps, pspre;
  PSBUS bus, buspre;
  PSGEN gen, genpre;

  PetscFunctionBegin;

  ierr = VecGetArray(tcopflow->Gl, &gl);
  CHKERRQ(ierr);
  ierr = VecGetArray(tcopflow->Gu, &gu);
  CHKERRQ(ierr);

  opflowpre = tcopflow->opflows[0];
  for (i = 0; i < tcopflow->Nt; i++) {
    if (i > 0)
      opflowpre = tcopflow->opflows[i - 1];
    opflow = tcopflow->opflows[i];

    /* Set bounds on constraints */
    gli = gl + tcopflow->gstarti[i];
    gui = gu + tcopflow->gstarti[i];

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

    if (tcopflow->nconineqcoup[i]) {
      ctr = 0;
      ps = opflow->ps;
      pspre = opflowpre->ps;
      /* Bounds on coupling constraints */
      for (j = 0; j < ps->nbus; j++) {
        bus = &ps->bus[j];
        buspre = &pspre->bus[j];

        for (k = 0; k < bus->ngen; k++) {
          ierr = PSBUSGetGen(bus, k, &gen);
          CHKERRQ(ierr);
          ierr = PSBUSGetGen(buspre, k, &genpre);
          CHKERRQ(ierr);
          if (!gen->status || !genpre->status)
            continue;
          /* Ramp constraints */
          gli[opflow->ncon + ctr] = -gen->ramp_rate_min * tcopflow->dT;
          gui[opflow->ncon + ctr] = gen->ramp_rate_min * tcopflow->dT;
          ctr++;
        }
      }
    }
  }

  ierr = VecRestoreArray(tcopflow->Gl, &gl);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(tcopflow->Gu, &gu);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode TCOPFLOWSetVariableandConstraintBounds_GENRAMP(TCOPFLOW tcopflow,
                                                              Vec Xl, Vec Xu,
                                                              Vec Gl, Vec Gu) {
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = TCOPFLOWSetVariableBounds_GENRAMP(tcopflow, Xl, Xu);
  CHKERRQ(ierr);
  ierr = TCOPFLOWSetConstraintBounds_GENRAMP(tcopflow, Gl, Gu);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode TCOPFLOWSetInitialGuess_GENRAMP(TCOPFLOW tcopflow, Vec X) {
  PetscErrorCode ierr;
  OPFLOW opflow;
  PetscScalar *x, *xi, *xl, *xu, *xli, *xui;
  PetscInt i;
  PetscFunctionBegin;

  ierr = VecGetArray(X, &x);
  CHKERRQ(ierr);
  ierr = VecGetArray(tcopflow->Xl, &xl);
  CHKERRQ(ierr);
  ierr = VecGetArray(tcopflow->Xu, &xu);
  CHKERRQ(ierr);

  for (i = 0; i < tcopflow->Nt; i++) {
    opflow = tcopflow->opflows[i];
    /* Set initial guess and bounds on variables */
    xi = x + tcopflow->xstarti[i];
    xli = xl + tcopflow->xstarti[i];
    xui = xu + tcopflow->xstarti[i];

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
  ierr = VecRestoreArray(tcopflow->Xl, &xl);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(tcopflow->Xu, &xu);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode TCOPFLOWComputeJacobian_GENRAMP(TCOPFLOW tcopflow, Vec X,
                                               Mat J) {
  PetscErrorCode ierr;
  OPFLOW opflow, opflowtpre;
  PetscInt roffset, coffset;
  PetscInt nrow, ncol;
  PetscScalar *xi, *x;
  PetscInt i, j, k, loc, loctpre, xtpreloc, xiloc;
  PS ps, pstpre;
  PSBUS bus, bustpre;
  PSGEN gen, gentpre;
  PetscInt nvals;
  const PetscInt *cols;
  const PetscScalar *vals;
  PetscInt row, col;
  PetscScalar val;

  PetscFunctionBegin;

  ierr = VecGetArray(X, &x);
  CHKERRQ(ierr);

  for (i = 0; i < tcopflow->Nt; i++) {
    opflow = tcopflow->opflows[i];

    roffset = tcopflow->gstarti[i];
    coffset = tcopflow->xstarti[i];

    xi = x + tcopflow->xstarti[i];
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

    if (tcopflow->nconineqcoup[i]) {
      ps = opflow->ps;
      pstpre = opflowtpre->ps;
      for (j = 0; j < ps->nbus; j++) {
        bus = &ps->bus[j];
        bustpre = &pstpre->bus[j];
        ierr = PSBUSGetVariableLocation(bus, &loc);
        CHKERRQ(ierr);
        ierr = PSBUSGetVariableLocation(bustpre, &loctpre);
        CHKERRQ(ierr);
        for (k = 0; k < bus->ngen; k++) {
          ierr = PSBUSGetGen(bus, k, &gen);
          CHKERRQ(ierr);
          ierr = PSBUSGetGen(bustpre, k, &gentpre);
          CHKERRQ(ierr);

          if (!gen->status) {
            if (gentpre->status)
              loctpre += 2;
            continue;
          } else {
            loc += 2;
            if (!gentpre->status)
              continue;
            loctpre += 2;
          }

          xtpreloc = tcopflow->xstarti[i - 1] + loctpre;
          xiloc = tcopflow->xstarti[i] + loc;
          row = roffset;
          col = xtpreloc;
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

    opflowtpre = opflow;
  }

  ierr = VecRestoreArray(X, &x);
  CHKERRQ(ierr);

  ierr = MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode TCOPFLOWComputeConstraints_GENRAMP(TCOPFLOW tcopflow, Vec X,
                                                  Vec G) {
  PetscErrorCode ierr;
  OPFLOW opflowtpre, opflow;
  PetscInt i, j, k, loc, loctpre, ctr;
  PetscScalar *xtpre, *x, *xi, *g, *gi;
  PS ps, pstpre;
  PSBUS bus, bustpre;
  PSGEN gen, gentpre;

  PetscFunctionBegin;

  ierr = VecGetArray(X, &x);
  CHKERRQ(ierr);
  ierr = VecGetArray(G, &g);
  CHKERRQ(ierr);

  for (i = 0; i < tcopflow->Nt; i++) {
    xi = x + tcopflow->xstarti[i];
    gi = g + tcopflow->gstarti[i];

    opflow = tcopflow->opflows[i];

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

    if (tcopflow->nconineqcoup[i]) {
      ctr = 0;
      ps = opflow->ps;
      pstpre = opflowtpre->ps;
      for (j = 0; j < ps->nbus; j++) {
        bus = &ps->bus[j];
        bustpre = &pstpre->bus[j];
        ierr = PSBUSGetVariableLocation(bus, &loc);
        CHKERRQ(ierr);
        ierr = PSBUSGetVariableLocation(bustpre, &loctpre);
        CHKERRQ(ierr);

        for (k = 0; k < bus->ngen; k++) {
          ierr = PSBUSGetGen(bus, k, &gen);
          CHKERRQ(ierr);
          ierr = PSBUSGetGen(bustpre, k, &gentpre);
          CHKERRQ(ierr);

          if (!gen->status) {
            if (gentpre->status)
              loctpre += 2;
            continue;
          } else {
            loc += 2;
            if (!gentpre->status)
              continue;
            loctpre += 2;
          }

          gi[ctr] = xi[loc] - xtpre[loctpre]; /* PG(t) - PG(t-dT) */
          ctr++;
        }
      }
    }

    ierr = VecResetArray(opflow->X);
    CHKERRQ(ierr);
    opflowtpre = opflow;
    xtpre = xi;
  }

  ierr = VecRestoreArray(X, &x);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(G, &g);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode TCOPFLOWComputeObjective_GENRAMP(TCOPFLOW tcopflow, Vec X,
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
  for (i = 0; i < tcopflow->Nt; i++) {
    opflowobj = 0.0;
    xi = x + tcopflow->xstarti[i];
    opflow = tcopflow->opflows[i];
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

PetscErrorCode TCOPFLOWComputeGradient_GENRAMP(TCOPFLOW tcopflow, Vec X,
                                               Vec Grad) {
  PetscErrorCode ierr;
  OPFLOW opflow;
  PetscInt i;
  PetscScalar *x, *xi, *grad, *gradi;

  ierr = VecGetArray(X, &x);
  CHKERRQ(ierr);
  ierr = VecGetArray(Grad, &grad);
  CHKERRQ(ierr);

  for (i = 0; i < tcopflow->Nt; i++) {
    opflow = tcopflow->opflows[i];
    xi = x + tcopflow->xstarti[i];
    gradi = grad + tcopflow->xstarti[i];

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

PetscErrorCode TCOPFLOWComputeObjandGradient_GENRAMP(TCOPFLOW tcopflow, Vec X,
                                                     PetscScalar *obj,
                                                     Vec Grad) {
  PetscErrorCode ierr;
  PetscFunctionBegin;

  ierr = TCOPFLOWComputeObjective_GENRAMP(tcopflow, X, obj);
  CHKERRQ(ierr);
  ierr = TCOPFLOWComputeGradient_GENRAMP(tcopflow, X, Grad);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode TCOPFLOWModelSetNumVariablesandConstraints_GENRAMP(
    TCOPFLOW tcopflow, PetscInt *nxi, PetscInt *ngi, PetscInt *nconineqcoup) {
  PetscInt i, ngenON;
  OPFLOW opflow;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  for (i = 0; i < tcopflow->Nt; i++) {
    opflow = tcopflow->opflows[i];
    ierr = PSGetNumActiveGenerators(opflow->ps, &ngenON, NULL);
    CHKERRQ(ierr);
    nxi[i] = opflow->nx;
    if (tcopflow->iscoupling)
      nconineqcoup[i] = (i == 0) ? 0 : ngenON;
    ngi[i] = opflow->ncon + nconineqcoup[i];
  }
  PetscFunctionReturn(0);
}

PetscErrorCode TCOPFLOWComputeHessian_GENRAMP(TCOPFLOW tcopflow, Vec X,
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
  for (i = 0; i < tcopflow->Nt; i++) {
    opflow = tcopflow->opflows[i];
    opflow->obj_factor = tcopflow->obj_factor;

    roffset = tcopflow->xstarti[i];

    xi = x + roffset;
    lameqi = lambda + tcopflow->gstarti[i];

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

PetscErrorCode TCOPFLOWModelCreate_GENRAMP(TCOPFLOW tcopflow) {
  GENRAMP genramp;
  PetscErrorCode ierr;

  PetscFunctionBegin;

  ierr = PetscCalloc1(1, &genramp);
  CHKERRQ(ierr);

  tcopflow->model = genramp;

  /* Inherit Ops */
  tcopflow->modelops.destroy = TCOPFLOWModelDestroy_GENRAMP;
  tcopflow->modelops.setnumvariablesandconstraints =
      TCOPFLOWModelSetNumVariablesandConstraints_GENRAMP;
  tcopflow->modelops.setvariablebounds = TCOPFLOWSetVariableBounds_GENRAMP;
  tcopflow->modelops.setconstraintbounds = TCOPFLOWSetConstraintBounds_GENRAMP;
  tcopflow->modelops.setvariableandconstraintbounds =
      TCOPFLOWSetVariableandConstraintBounds_GENRAMP;
  tcopflow->modelops.setinitialguess = TCOPFLOWSetInitialGuess_GENRAMP;
  tcopflow->modelops.computeconstraints = TCOPFLOWComputeConstraints_GENRAMP;
  tcopflow->modelops.computejacobian = TCOPFLOWComputeJacobian_GENRAMP;
  tcopflow->modelops.computehessian = TCOPFLOWComputeHessian_GENRAMP;
  tcopflow->modelops.computeobjandgradient =
      TCOPFLOWComputeObjandGradient_GENRAMP;
  tcopflow->modelops.computeobjective = TCOPFLOWComputeObjective_GENRAMP;
  tcopflow->modelops.computegradient = TCOPFLOWComputeGradient_GENRAMP;

  PetscFunctionReturn(0);
}
