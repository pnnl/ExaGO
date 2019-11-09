# To use forward sensitivity, add option -DFWDSA to CFLAGS
# override CFLAGS += -DFWDSA -O2 -Iinclude
# For debugging
# override CFLAGS += -Iinclude -DPFLOW_DISPLAY_RESULTS -DDEBUGPS

CFLAGS += -Iinclude 
CFLAGS += -I${IPOPT_BUILD_DIR}/include/coin # For IPOPT
CFLAGS += -I${PIPS_DIR}/PIPS-NLP # For PIPS-NLP
FFLAGS           =
CPPFLAGS         =
FPPFLAGS         =
## Turn this ON for using IPOPT
CFLAGS_IPOPT     = #-DSCOPFLOW_HAVE_IPOPT
# Turn this ON for using PIPS
CFLAGS_PIPS      = #-DSCOPFLOW_HAVE_PIPS

CFLAGS += ${CFLAGS_IPOPT} ${CFLAGS_PIPS}

OS := $(shell uname)
ifeq ($(OS),Darwin)
  OTHER_LIB =
  LIB_EXT = dylib
  LDFLAGS = -dynamiclib
else
  OTHER_LIB = -lrt
  LIB_EXT = so
  LDFLAGS = -shared
endif

ALL:

include $(PETSC_DIR)/lib/petsc/conf/variables
include $(PETSC_DIR)/lib/petsc/conf/rules

#******************************
#	Genral use Objects
#******************************

PS_SRC_OBJECTS = src/ps/ps.o src/ps/psreaddata.o src/ps/psislanding.o src/utils/comm.o src/utils/utils.o

#******************************
#	PFLOW Specific Make
#******************************
PFLOW_SRC_OBJECTS = src/pflow/pflow.o ${PS_SRC_OBJECTS}

#******** Option 1 **********
PFLOW_APP_OBJECTS = applications/pflow-main.o
OBJECTS_PFLOW = $(PFLOW_APP_OBJECTS)
PFLOW: $(OBJECTS_PFLOW) libpflow chkopts
	 -$(CLINKER) -o PFLOW $(OBJECTS_PFLOW) ${PETSC_SNES_LIB} -L${SCOPFLOW_DIR} -lpflow
	$(RM) $(OBJECTS_PFLOW)


PFLOW_APP2_OBJECTS = applications/pflow-main2.o
#******** Option 2 **********
OBJECTS_PFLOW2 = $(PFLOW_APP2_OBJECTS)
PFLOW2: $(OBJECTS_PFLOW2) libpflow chkopts
	 -$(CLINKER) -o PFLOW2 $(OBJECTS_PFLOW2) ${PETSC_SNES_LIB} -L${SCOPFLOW_DIR} -lpflow
	$(RM) $(OBJECTS_PFLOW2)


#******************************
#	OPFLOW Specific Make
#******************************
OPFLOW_SRC_OBJECTS = src/opflow/interface/opflow.o src/opflow/interface/opflowregi.o src/opflow/formulation/power-bal-polar/pbpol.o src/opflow/formulation/power-bal-cartesian/pbcar.o src/opflow/formulation/current-bal-cartesian/ibcar.o src/opflow/solver/ipopt/opflow-ipopt.o src/opflow/solver/tao/opflow-tao.o ${PS_SRC_OBJECTS}

OPFLOW_APP_OBJECTS = applications/opflow-main.o
OBJECTS_OPFLOW = $(OPFLOW_APP_OBJECTS)
OPFLOW: $(OBJECTS_OPFLOW) libopflow chkopts
	 -$(CLINKER) -o OPFLOW $(OBJECTS_OPFLOW) ${PETSC_TAO_LIB} -L${SCOPFLOW_DIR} -lopflow
	$(RM) $(OBJECTS_OPFLOW)


#******************************
#	SCOPFLOW Specific Make
#******************************
SCOPFLOW_APP_OBJECTS = applications/scopflow-main.o
#******** Option 2 **********
OBJECTS_SCOPFLOW2 = $(SCOPFLOW_APP_OBJECTS) 
SCOPFLOW_IPOPT: $(OBJECTS_SCOPFLOW2) libscopflowipopt chkopts
	 -$(CLINKER) -o SCOPFLOW_IPOPT $(OBJECTS_SCOPFLOW2) -L${SCOPFLOW_DIR} -lscopflowipopt ${PETSC_LIB}
	$(RM) $(OBJECTS_SCOPFLOW2)

# SCOPFLOW with PIPS
OBJECTS_SCOPFLOW3 = $(SCOPFLOW_APP_OBJECTS) 
SCOPFLOW_PIPS: $(OBJECTS_SCOPFLOW3) libscopflowpips chkopts
	 -$(CLINKER) -o SCOPFLOW_PIPS $(OBJECTS_SCOPFLOW3) -L${SCOPFLOW_DIR} -lscopflowpips ${PETSC_LIB}
	$(RM) $(OBJECTS_SCOPFLOW3)

#***************************
#	Make Library Commands
#***************************
libpflow:$(PFLOW_SRC_OBJECTS) chkopts
	 -$(CLINKER) $(LDFLAGS) -o libpflow.$(LIB_EXT) $(PFLOW_SRC_OBJECTS) $(PETSC_TS_LIB)

libopflow:$(OPFLOW_SRC_OBJECTS) chkopts
	 -$(CLINKER) $(LDFLAGS) -o libopflow.$(LIB_EXT) $(OPFLOW_SRC_OBJECTS)  -L${IPOPT_BUILD_DIR}/lib -lipopt $(PETSC_TAO_LIB)


SCOPFLOW_IPOPT_SRC_OBJECTS = src/scopflow/scopflow-ipopt.o src/scopflow/scopflow-ipopt-constraints.o src/scopflow/scopflow-ipopt-objective.o src/scopflow/scopflow-ipopt-hessian.o src/scopflow/scopflow-ipopt-pipsfunctions.o ${OPFLOW_IPOPT_SRC_OBJECTS}
libscopflowipopt:$(SCOPFLOW_IPOPT_SRC_OBJECTS) chkopts
	 -$(CLINKER) $(LDFLAGS) -o libscopflowipopt.$(LIB_EXT) $(SCOPFLOW_IPOPT_SRC_OBJECTS) -L${IPOPT_BUILD_DIR}/lib -lipopt $(PETSC_TAO_LIB)

SCOPFLOW_PIPS_SRC_OBJECTS = src/scopflow/scopflow-pips.o src/scopflow/scopflow-pips-constraints.o src/scopflow/scopflow-pips-hessian.o src/scopflow/scopflow-pips-objective.o ${OPFLOW_SRC_OBJECTS}
libscopflowpips:$(SCOPFLOW_PIPS_SRC_OBJECTS) chkopts
	 -$(CLINKER) $(LDFLAGS) -o libscopflowpips.$(LIB_EXT) $(SCOPFLOW_PIPS_SRC_OBJECTS) -L${PIPS_DIR}/build/PIPS-NLP -lparpipsnlp $(PETSC_TAO_LIB)

#******************************
#	Remove .o Command
#******************************
cleanobj:
	rm -rf $(OBJECTS_PFLOW) $(OBJECTS_PFLOW2) $(PFLOW_SRC_OBJECTS) $(OBJECTS_OPFLOW) $(OPFLOW_SRC_OBJECTS) $(OBJECTS_SCOPFLOW2) $(SCOPFLOW_IPOPT_SRC_OBJECTS) $(SCOPFLOW_PIPS_SRC_OBJECTS) *.dylib *.dSYM PFLOW PFLOW2 OPFLOW SCOPFLOW_IPOPT SCOPFLOW_PIPS
