
BUILD_WITH_IPOPT=${WITH_IPOPT}
BUILD_WITH_PIPS=${WITH_PIPS}
BUILD_WITH_HIOP=${WITH_HIOP}

CFLAGS += -Iinclude 
CXXFLAGS += -Iinclude
FFLAGS           =
CPPFLAGS         =
FPPFLAGS         =

ifeq ($(BUILD_WITH_PIPS),1)
  CFLAGS_PIPS = -DEXAGO_HAVE_PIPS
  PIPS_LIB    = -lparpipsnlp
  CFLAGS += -I${PIPS_DIR}/PIPS-NLP # For PIPS-NLP
  CXXFLAGS += -I${PIPS_DIR}/PIPS-NLP # For PIPS-NLP

  CFLAGS_IPOPT = -DEXAGO_HAVE_IPOPT
  IPOPT_LIB    = -lipopt
  CFLAGS += -I${IPOPT_BUILD_DIR}/include/coin
  CXXFLAGS += -I${IPOPT_BUILD_DIR}/include/coin

endif

ifeq ($(BUILD_WITH_IPOPT),1)
  CFLAGS_IPOPT = -DEXAGO_HAVE_IPOPT
  IPOPT_LIB    = -lipopt
  CFLAGS += -I${IPOPT_BUILD_DIR}/include/coin
  CXXFLAGS += -I${IPOPT_BUILD_DIR}/include/coin
endif

ifeq ($(BUILD_WITH_HIOP),1)
  CFLAGS_HIOP = -DEXAGO_HAVE_HIOP
  HIOP_LIB    = -lhiop
  CFLAGS += -I${HIOP_DIR}/include
  CXXFLAGS += -I${HIOP_DIR}/include #-I${HIOP_SRC_DIR}/LinAlg -I${HIOP_SRC_DIR}/Optimization
endif

CFLAGS += ${CFLAGS_IPOPT} ${CFLAGS_PIPS} ${CFLAGS_HIOP}
CXXFLAGS += ${CFLAGS_IPOPT} ${CFLAGS_PIPS} ${CFLAGS_HIOP}


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
#	General use Objects
#******************************

PS_SRC_OBJECTS = src/ps/ps.o src/ps/psreaddata.o src/ps/psislanding.o src/ps/psoutput.o src/utils/comm.o src/utils/utils.o

#******************************
#	PFLOW Specific Make
#******************************
PFLOW_SRC_OBJECTS = src/pflow/pflow.o

#******** Option 1 **********
PFLOW_APP_OBJECTS = applications/pflow-main.o
OBJECTS_PFLOW = $(PFLOW_APP_OBJECTS)

PFLOW: $(OBJECTS_PFLOW) libpflow chkopts
	 -$(CLINKER) -o PFLOW $(OBJECTS_PFLOW) ${PETSC_SNES_LIB} -L${EXAGO_DIR} -lpflow
	$(RM) $(OBJECTS_PFLOW)

PFLOW_PROTOAPP_OBJECTS = applications/pflow-proto.o
OBJECTS_PFLOWPROTO = $(PFLOW_PROTOAPP_OBJECTS)
PFLOWPROTO: $(OBJECTS_PFLOWPROTO) libpflow chkopts
	 -$(CLINKER) -o PFLOWPROTO $(OBJECTS_PFLOWPROTO) ${PETSC_SNES_LIB} -L${EXAGO_DIR} -lpflow
	$(RM) $(OBJECTS_PFLOWPROTO)

#******** Option 2 **********
PFLOW_APP2_OBJECTS = applications/pflow-main2.o
OBJECTS_PFLOW2 = $(PFLOW_APP2_OBJECTS)

PFLOW2: $(OBJECTS_PFLOW2) libpflow chkopts
	 -$(CLINKER) -o PFLOW2 $(OBJECTS_PFLOW2) ${PETSC_SNES_LIB} -L${EXAGO_DIR} -lpflow
	$(RM) $(OBJECTS_PFLOW2)


#******************************
#	OPFLOW Specific Make
#******************************
OPFLOW_INTERFACE_OBJECTS = src/opflow/interface/opflow.o src/opflow/interface/opflowregi.o src/opflow/interface/opflowoutput.o
OPFLOW_FORMULATION_OBJECTS = src/opflow/formulation/power-bal-polar/pbpol.o src/opflow/formulation/power-bal-cartesian/pbcar.o src/opflow/formulation/current-bal-cartesian/ibcar.o src/opflow/formulation/current-bal-cartesian/ibcar2.o

OPFLOW_SOLVER_OBJECTS = src/opflow/solver/ipopt/opflow-ipopt.o src/opflow/solver/tao/opflow-tao.o src/opflow/solver/hiop/opflow-hiop.o

OPFLOW_SRC_OBJECTS = ${PFLOW_SRC_OBJECTS} ${OPFLOW_INTERFACE_OBJECTS} ${OPFLOW_FORMULATION_OBJECTS} ${OPFLOW_SOLVER_OBJECTS} ${PS_SRC_OBJECTS}

OPFLOW_APP_OBJECTS = applications/opflow-main.o
OBJECTS_OPFLOW = $(OPFLOW_APP_OBJECTS)

OPFLOW: $(OBJECTS_OPFLOW) libopflow chkopts
	 -$(CLINKER) -o OPFLOW $(OBJECTS_OPFLOW) ${PETSC_TAO_LIB} -L${EXAGO_DIR} -lopflow
	$(RM) $(OBJECTS_OPFLOW)

OPFLOW_PROTOAPP_OBJECTS = applications/opflow-proto.o
OBJECTS_OPFLOWP = $(OPFLOW_PROTOAPP_OBJECTS)
OPFLOWPROTO: $(OBJECTS_OPFLOWP) libopflow chkopts
	 -$(CLINKER) -o OPFLOWPROTO $(OBJECTS_OPFLOWP) ${PETSC_TAO_LIB} -L${EXAGO_DIR} -lopflow
	$(RM) $(OBJECTS_OPFLOWP)


OPFLOW_PROTOAPP2_OBJECTS = applications/opflow-proto2.o
OBJECTS_OPFLOWP2 = $(OPFLOW_PROTOAPP2_OBJECTS)
OPFLOWPROTO2: $(OBJECTS_OPFLOWP2) libopflow chkopts
	 -$(CLINKER) -o OPFLOWPROTO2 $(OBJECTS_OPFLOWP2) ${PETSC_TAO_LIB} -L${EXAGO_DIR} -lopflow
	$(RM) $(OBJECTS_OPFLOWP2)

#******************************
#	SCOPFLOW Specific Make
#******************************
EXAGO_INTERFACE_OBJECTS = src/scopflow/interface/scopflow.o src/scopflow/interface/scopflowregi.o
EXAGO_SOLVER_OBJECTS = src/scopflow/solver/ipopt/scopflow-ipopt.o src/scopflow/solver/pips/scopflow-pips.o

EXAGO_SRC_OBJECTS = ${EXAGO_INTERFACE_OBJECTS} ${EXAGO_SOLVER_OBJECTS} ${OPFLOW_SRC_OBJECTS}

EXAGO_APP_OBJECTS = applications/scopflow-main.o

OBJECTS_SCOPFLOW = $(EXAGO_APP_OBJECTS) 
SCOPFLOW: $(OBJECTS_SCOPFLOW) libscopflow chkopts
	 -$(CLINKER) -o SCOPFLOW $(OBJECTS_SCOPFLOW) -L${EXAGO_DIR} -lscopflow ${PETSC_LIB}
	$(RM) $(OBJECTS_SCOPFLOW)

#***************************
#	Make Library Commands
#***************************
libpflow:${PS_SRC_OBJECTS} $(PFLOW_SRC_OBJECTS) chkopts
	 -$(CLINKER) $(LDFLAGS) -o libpflow.$(LIB_EXT) ${PS_SRC_OBJECTS} $(PFLOW_SRC_OBJECTS) $(PETSC_TS_LIB)

libopflow:$(OPFLOW_SRC_OBJECTS) chkopts
	 -$(CLINKER) $(LDFLAGS) -o libopflow.$(LIB_EXT) $(OPFLOW_SRC_OBJECTS) -L${IPOPT_BUILD_DIR}/lib ${IPOPT_LIB} -L${HIOP_DIR}/lib ${HIOP_LIB} $(PETSC_TAO_LIB)

libscopflow:${EXAGO_SRC_OBJECTS} chkopts
	 -$(CLINKER) $(LDFLAGS) -o libscopflow.$(LIB_EXT) $(EXAGO_SRC_OBJECTS) -L${PIPS_DIR}/build_pips/PIPS-NLP ${PIPS_LIB} -L${IPOPT_BUILD_DIR}/lib ${IPOPT_LIB} ${PETSC_TAO_LIB}
#******************************
#	Remove .o Command
#******************************
cleanobj:
	rm -rf $(OBJECTS_PFLOW) $(OBJECTS_PFLOW2) $(PFLOW_SRC_OBJECTS) $(OBJECTS_OPFLOW) $(OPFLOW_SRC_OBJECTS) $(EXAGO_SRC_OBJECTS) ${OBJECTS_SCOPFLOW} *.dylib *.so* *.dSYM PFLOW PFLOW2 PFLOWPROTO OPFLOW SCOPFLOW
