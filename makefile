# To use forward sensitivity, add option -DFWDSA to CFLAGS
# override CFLAGS += -DFWDSA -O2 -Iinclude
# For debugging
# override CFLAGS += -Iinclude -DPFLOW_DISPLAY_RESULTS -DDEBUGPS

CFLAGS += -Iinclude -I${IPOPT_BUILD_DIR}/include/coin
FFLAGS           =
CPPFLAGS         =
FPPFLAGS         =
CFLAGS_IPOPT     = #-DPSAPPS_HAVE_IPOPT

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

DYNGENMODEL_OBJECTS = src/dyn/dyngenmodels/dyngenmodels.o src/dyn/dyngenmodels/dyngenrou.o src/dyn/dyngenmodels/dynpvd1.o src/dyn/dyngenmodels/dyncv.o
DYNEXCMODEL_OBJECTS = src/dyn/dynexcmodels/dynexcmodels.o src/dyn/dynexcmodels/dynieeet1.o src/dyn/dynexcmodels/dynexst1.o src/dyn/dynexcmodels/dynsexs.o
DYNTURBGOVMODEL_OBJECTS = src/dyn/dynturbgovmodels/dynturbgovmodels.o src/dyn/dynturbgovmodels/dyntgov1.o
DYNSTABMODEL_OBJECTS = src/dyn/dynstabmodels/dynstabmodels.o src/dyn/dynstabmodels/dynstab1.o
DYNLOADMODEL_OBJECTS = src/dyn/dynloadmodels/dynloadmodels.o src/dyn/dynloadmodels/dynzip.o src/dyn/dynloadmodels/dyncompload.o
DYNEVENT_OBJECTS = src/dyn/dynevents.o src/dyn/dynfaultevents.o src/dyn/dynlineswevents.o src/dyn/dyngentripevents.o

DYN_SRC_OBJECTS = src/ps/ps.o src/utils/comm.o src/utils/utils.o src/pflow/pflow.o src/dyn/dyn.o ${DYNGENMODEL_OBJECTS} ${DYNEXCMODEL_OBJECTS} ${DYNTURBGOVMODEL_OBJECTS} ${DYNSTABMODEL_OBJECTS} ${DYNLOADMODEL_OBJECTS} ${DYNEVENT_OBJECTS}

DYN_APP_OBJECTS = applications/dyn-main.o
OBJECTS_DYN = $(DYN_APP_OBJECTS)
DYN: $(OBJECTS_DYN) libdyn chkopts
	 -$(CLINKER) -o DYN $(OBJECTS_DYN) ${PETSC_TS_LIB} -L${PSAPPS_DIR} -ldyn
	$(RM) $(OBJECTS_DYN)

PFLOW_SRC_OBJECTS = src/ps/ps.o src/utils/comm.o src/utils/utils.o src/pflow/pflow.o ${DYNGENMODEL_OBJECTS} ${DYNEXCMODEL_OBJECTS} ${DYNTURBGOVMODEL_OBJECTS} ${DYNSTABMODEL_OBJECTS} ${DYNLOADMODEL_OBJECTS}

PFLOW_APP_OBJECTS = applications/pflow-main.o
OBJECTS_PFLOW = $(PFLOW_APP_OBJECTS)
PFLOW: $(OBJECTS_PFLOW) libpflow chkopts
	 -$(CLINKER) -o PFLOW $(OBJECTS_PFLOW) ${PETSC_SNES_LIB} -L${PSAPPS_DIR} -lpflow
	$(RM) $(OBJECTS_PFLOW)

PFLOW_APP2_OBJECTS = applications/pflow-main2.o
OBJECTS_PFLOW2 = $(PFLOW_APP2_OBJECTS)
PFLOW2: $(OBJECTS_PFLOW2) libpflow chkopts
	 -$(CLINKER) -o PFLOW2 $(OBJECTS_PFLOW2) ${PETSC_SNES_LIB} -L${PSAPPS_DIR} -lpflow
	$(RM) $(OBJECTS_PFLOW2)

OPFLOW_SRC_OBJECTS = src/ps/ps.o src/utils/comm.o src/utils/utils.o src/opflow/opflow.o src/opflow/econstraints.o src/opflow/iconstraints.o ${DYNGENMODEL_OBJECTS} ${DYNEXCMODEL_OBJECTS} ${DYNTURBGOVMODEL_OBJECTS} ${DYNSTABMODEL_OBJECTS} ${DYNLOADMODEL_OBJECTS}

OPFLOW_APP_OBJECTS = applications/opflow-main.o
OBJECTS_OPFLOW = $(OPFLOW_APP_OBJECTS)
OPFLOW: $(OBJECTS_OPFLOW) libopflow chkopts
	 -$(CLINKER) -o OPFLOW $(OBJECTS_OPFLOW) ${PETSC_TAO_LIB} -L${PSAPPS_DIR} -lopflow
	$(RM) $(OBJECTS_OPFLOW)

OPFLOW_IPOPT_SRC_OBJECTS = src/ps/ps.o src/utils/comm.o src/utils/utils.o src/opflow/opflow-ipopt.o ${DYNGENMODEL_OBJECTS} ${DYNEXCMODEL_OBJECTS} ${DYNTURBGOVMODEL_OBJECTS} ${DYNSTABMODEL_OBJECTS} ${DYNLOADMODEL_OBJECTS}


OBJECTS_OPFLOW2 = $(OPFLOW_APP_OBJECTS)
OPFLOW_IPOPT: $(OBJECTS_OPFLOW2) libopflowipopt chkopts
	 -$(CLINKER) -o OPFLOW_IPOPT $(OBJECTS_OPFLOW2) -L${PSAPPS_DIR} -lopflowipopt
	$(RM) $(OBJECTS_OPFLOW2)


libdyn:$(DYN_SRC_OBJECTS) chkopts
	 -$(CLINKER) $(LDFLAGS) -o libdyn.$(LIB_EXT) $(DYN_SRC_OBJECTS) $(PETSC_TS_LIB)

CFLAGS += ${CFLAGS_IPOPT}
libpflow:$(PFLOW_SRC_OBJECTS) chkopts
	 -$(CLINKER) $(LDFLAGS) -o libpflow.$(LIB_EXT) $(PFLOW_SRC_OBJECTS) $(PETSC_TS_LIB)

libopflow:$(OPFLOW_SRC_OBJECTS) chkopts
	 -$(CLINKER) $(LDFLAGS) -o libopflow.$(LIB_EXT) $(OPFLOW_SRC_OBJECTS) $(PETSC_TAO_LIB)

libopflowipopt:$(OPFLOW_IPOPT_SRC_OBJECTS) chkopts
	 -$(CLINKER) $(LDFLAGS) -o libopflowipopt.$(LIB_EXT) $(OPFLOW_IPOPT_SRC_OBJECTS) -L${IPOPT_BUILD_DIR}/lib -lipopt $(PETSC_TAO_LIB)

cleanobj:
	rm -rf $(OBJECTS_PFLOW) $(OBJECTS_PFLOW2) $(PFLOW_SRC_OBJECTS) $(OBJECTS_OPFLOW) $(OPFLOW_SRC_OBJECTS) $(OBJECTS_DYN) $(DYN_SRC_OBJECTS) *.dylib *.dSYM PFLOW PFLOW2 DYN OPFLOW
