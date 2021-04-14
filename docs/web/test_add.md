### Creating a new test
Create a test for a new application (referred to here as "myapp") requires
several steps. First the new application is available as a library
and all runtime parameters can be specified as runtime arguments. In addition,
the new application has already been incorporated into the build system using
the previous applications (such as pflow, opflow, sopflow, etc.) as templates.

The first step in creating tests is to create a directory myapp in the
tests/functionality directory to hold the CMakeLists.txt and C files for the new
application. Two C files are needed for the test. The first is myapp.c that is
used to actually run the application. Files from other applications such as
pflow, opflow, etc. can be used as templates. Second, a file myappselfcheck.c
and its corresponding header file myappselfcheck.h need to be created in the
directory.

The myappselfcheck.h file defines a C-struct that defines the outputs from the
calculation that are going to be checked as well as some of the attributes of
the calculations themselves (solver, network, etc). The results that are checked
are typically the number of iterations and the value of the objective function
which are obtained from some calculation that is considered to be the "truth".
The myappselfcheck.c file defines a collection of results from different tests
that are stored as an array of the structs defined in the myappselfcheck.h file.
Again, files from other applications can be used as a template. The selfcheck
function defined in myappselfcheck.c is called from the main routine in myapp.c.

### Adding test to CMake build
Adding a new test to the CMake build requires modifications to several files.
The first is to add a new CMakeLists.txt file to the testing directory in
tests/functionality/myapp. This is fairly minimal and has the form
```
#[[ Functionality tests for ExaGO MYAPP applications ]]

add_executable(testMyappFunctionality myapp.c myappselfcheck.c)
target_include_directories(testMyappFunctionality PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

#[[ Configure MYAPP functionality tests ]]
target_link_libraries(testMyappFunctionality ${EXAGO_APP_LIBS})
if(EXAGO_INSTALL_TESTS)
  install(TARGETS testMyappFunctionality RUNTIME DESTINATION bin/tests/functionality)
endif()
```
The <code>target_include_directories</code> line is included so that the notation
```
#include <myapp.h>
```
using the angle brackets works instead of having to use double quotes.

The CMakeLists.txt file in the tests/functionality directory needs to be
modified so that the library <code>EXAGO_MYAPP_LIBS</code> is included in the
<code>EXAGO_APP_LIBS</code> variable.
```
set(EXAGO_APP_LIBS
  ${EXAGO_MYAPP_LIBS}
  ${EXAGO_SOPFLOW_LIBS}
  ${EXAGO_SCOPFLOW_LIBS}
  ${EXAGO_TCOPFLOW_LIBS}
  ExaGO::OPFLOW
  ${EXAGO_PFLOW_LIBS}
  ExaGO::UTILS
  ExaGO::PS
  ${EXAGO_MATH_LIB}
  PETSC::SNES
)
```
The position of <code>EXAGO_MYAPP_LIBS</code> in the list may vary depending on
what dependencies the new application has on other applications and vice-versa.
The subdirectory containing the new tests must also be added to to the
CMakeLists.txt file using the command
```
add_subdirectory(myapp)
```

Finally, the CMakeLists.txt file in the top level ExaGO directory must by
modified to actually specify which tests of myapp are to be run. All tests are
included in the
```
if(RUN_EXAGO_TESTS)
```
block, which is located at the bottom of the file. These tests assume that it is
possible to completely specify all parameters needed to run the application
using runtime arguments. These include network configuration files, contingency
files and scenario files. Note that some network files have already been
included in the repository and are available to any application
```
    datafiles/case9/case9mod.m
    datafiles/case118.m
    datafiles/case_ACTIVSg200.m
```
Existing applications can be used as prototypes for adding more to the top level
CMakeLists.txt file. A relatively simple example is the SOPFLOW application test
```
  message(STATUS "Configuring SOPFLOW functionality tests")
  foreach(ns RANGE 1 3)
    exago_add_test(
      NAME FUNCTIONALITY_TEST_SOPFLOW_200bus_IPOPT_${ns}scen
      DEPENDS IPOPT
      NETFILES datafiles/case_ACTIVSg200.m
      COMMAND ${RUNCMD} $<TARGET_FILE:testSopflowFunctionality>
      -scenfile datafiles/TAMU200_scenarios/scenarios_200bus.csv
      -sopflow_solver IPOPT
      -sopflow_enable_multicontingency 0
      -sopflow_Ns ${ns}
      -sopflow_tolerance ${tolerance}
      )
  endforeach(ns)
```
This tests SOPFLOW on 1 to 3 scenarios using <code>case_ACTIVSg200.m</code>
network configuration file and the corresponding
<code>scenarios_200bus.csv</code> scenario file. The <code>DEPENDS IPOPT</code>
line means that this test requires that ExaGO is configured with the IPOPT
solver. The
```
foreach<ns RANGE 1 3)
endforeach(ns)
```
specifies that the number of scenarios varies from 1 to 3 (the
<code>scenarios_200bus.csv</code> file must contain at least 3 scenarios). The
<code>testSopflowFunctionality</code> variable has already been defined in the
CMakeLists.txt file in the test/functionality/sopflow directory.
