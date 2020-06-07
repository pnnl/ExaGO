#pragma once

namespace exago{

  namespace tests{

    int TestOPFLOWComputeObjective(OPFLOW opflow, Vec X, PetscScalar obj_ref)
    {
      PetscErrorCode ierr;
      int            return_val=0;
      double         obj_val;
      
      ierr = OPFLOWComputeObjective(opflow,X,&obj_val);CHKERRQ(ierr);
      if(abs((obj_val - obj_ref)/obj_ref) > 1e-6) return_val = 1;
      return return_val;
    }
    
    int TestOPFLOWComputeGradient(OPFLOW opflow, Vec X, Vec grad_ref)
    {
      PetscErrorCode ierr;
      int            return_val=0;
      Vec            grad;
      PetscReal      vecnorm;
      
      ierr = VecDuplicate(grad_ref,&grad);CHKERRQ(ierr);
      ierr = OPFLOWComputeGradient(opflow,X,grad);CHKERRQ(ierr);
      
      ierr = VecAXPY(grad,-1.0,grad_ref);CHKERRQ(ierr);
      ierr = VecNorm(grad,NORM_2,&vecnorm);CHKERRQ(ierr);
      if(vecnorm > 1e-6) return_val = 1;
      ierr = VecDestroy(&grad);CHKERRQ(ierr);
      return return_val;
    }

    int TestOPFLOWComputeVariableBounds(OPFLOW opflow,Vec Xlref,Vec Xuref)
    {
      PetscErrorCode ierr;
      int            return_val=0;
      Vec            Xl,Xu;
      PetscReal      dxlnorm,dxunorm;
      
      ierr = VecDuplicate(Xlref,&Xl);CHKERRQ(ierr);
      ierr = VecDuplicate(Xuref,&Xu);CHKERRQ(ierr);

      ierr = OPFLOWComputeVariableBounds(opflow,Xl,Xu);CHKERRQ(ierr);
      
      ierr = VecAXPY(Xl,-1.0,Xlref);CHKERRQ(ierr);
      ierr = VecAXPY(Xu,-1.0,Xuref);CHKERRQ(ierr);

      ierr = VecNorm(Xl,NORM_2,&dxlnorm);CHKERRQ(ierr);
      ierr = VecNorm(Xu,NORM_2,&dxunorm);CHKERRQ(ierr);

      if(dxlnorm > 1e-6 || dxunorm > 1e-6) return_val = 1;
      ierr = VecDestroy(&Xl);CHKERRQ(ierr);
      ierr = VecDestroy(&Xu);CHKERRQ(ierr);
      return return_val;
    }

    int TestOPFLOWComputeConstraints(OPFLOW opflow,Vec X,Vec Gref)
    {
      PetscErrorCode ierr;
      int            return_val=0;
      Vec            G;
      PetscReal      dgnorm;
      
      ierr = VecDuplicate(Gref,&G);CHKERRQ(ierr);
      ierr = OPFLOWComputeConstraints(opflow,X,G);CHKERRQ(ierr);
      
      ierr = VecAXPY(G,-1.0,Gref);CHKERRQ(ierr);
      ierr = VecNorm(G,NORM_2,&dgnorm);CHKERRQ(ierr);
      if(dgnorm > 1e-4) return_val = 1;
      ierr = VecDestroy(&G);CHKERRQ(ierr);
      return return_val;
    }

    int TestOPFLOWComputeConstraintBounds(OPFLOW opflow,Vec Glref,Vec Guref)
    {
      PetscErrorCode ierr;
      int            return_val=0;
      Vec            Gl,Gu;
      PetscReal      dglnorm,dgunorm;
      
      ierr = VecDuplicate(Glref,&Gl);CHKERRQ(ierr);
      ierr = VecDuplicate(Guref,&Gu);CHKERRQ(ierr);

      ierr = OPFLOWComputeConstraintBounds(opflow,Gl,Gu);CHKERRQ(ierr);
      
      ierr = VecAXPY(Gl,-1.0,Glref);CHKERRQ(ierr);
      ierr = VecAXPY(Gu,-1.0,Guref);CHKERRQ(ierr);

      ierr = VecNorm(Gl,NORM_2,&dglnorm);CHKERRQ(ierr);
      ierr = VecNorm(Gu,NORM_2,&dgunorm);CHKERRQ(ierr);

      if(dglnorm > 1e-6 || dgunorm > 1e-6) return_val = 1;
      ierr = VecDestroy(&Gl);CHKERRQ(ierr);
      ierr = VecDestroy(&Gu);CHKERRQ(ierr);
      return return_val;
    }
    
  } // namespace tests
  
} // namespace exago
