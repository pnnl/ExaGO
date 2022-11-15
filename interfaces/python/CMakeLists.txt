# New true python bindings pybind11_add_module(exago exago_python.cpp)
# target_link_libraries(exago PRIVATE ExaGO::UTILS ExaGO::OPFLOW)

# Install Python bindings to same installation path: install(TARGETS exago
# DESTINATION ${EXAGO_PYTHON_SITELIB})

# New true python bindings
pybind11_add_module(exago exago_python.cpp)

# Link against relevant exago libraries
target_link_libraries(
  exago PRIVATE ExaGO::UTILS ExaGO::SOPFLOW ExaGO::SCOPFLOW ExaGO::OPFLOW
)

# Find mpi4py c headers
include(ExaGOFindPythonModule)

find_python_module(mpi4py REQUIRED)

# Get the include directories needed
if(mpi4py_FOUND)
  execute_process(
    COMMAND "${Python_EXECUTABLE}" "-c"
            "import mpi4py as m; print(m.__version__); print(m.get_include());"
    RESULT_VARIABLE _mpi4py_SEARCH_SUCCESS
    OUTPUT_VARIABLE _mpi4py_VALUES
    ERROR_VARIABLE _mpi4py_ERROR_VALUE
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  # Find pytest c headers
  find_python_module(pytest REQUIRED)

  # Get the include directories needed
  if(pytest_FOUND)
    execute_process(
      COMMAND "${Python_EXECUTABLE}" "-c"
              "import pytest as p; print(p.__version__);"
      RESULT_VARIABLE _pytest_SEARCH_SUCCESS
      OUTPUT_VARIABLE _pytest_VALUES
      ERROR_VARIABLE _pytest_ERROR_VALUE
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    # Convert the process output into a list
    string(REGEX REPLACE ";" "\\\\;" _mpi4py_VALUES ${_mpi4py_VALUES})
    string(REGEX REPLACE "\n" ";" _mpi4py_VALUES ${_mpi4py_VALUES})
    list(GET _mpi4py_VALUES 0 mpi4py_VERSION)
    list(GET _mpi4py_VALUES 1 mpi4py_INCLUDE_DIRS)

    # Make sure all directory separators are '/'
    string(REGEX REPLACE "\\\\" "/" mpi4py_INCLUDE_DIRS ${mpi4py_INCLUDE_DIRS})

    # Get the major and minor version numbers
    string(REGEX REPLACE "\\." ";" _mpi4py_VERSION_LIST ${mpi4py_VERSION})
    list(GET _mpi4py_VERSION_LIST 0 mpi4py_VERSION_MAJOR)
    list(GET _mpi4py_VERSION_LIST 1 mpi4py_VERSION_MINOR)
    list(GET _mpi4py_VERSION_LIST 2 mpi4py_VERSION_PATCH)
    string(REGEX MATCH "[0-9]*" mpi4py_VERSION_PATCH ${mpi4py_VERSION_PATCH})
    math(
      EXPR
      mpi4py_VERSION_DECIMAL
      "(${mpi4py_VERSION_MAJOR} * 10000) + (${mpi4py_VERSION_MINOR} * 100) + ${mpi4py_VERSION_PATCH}"
    )
  endif()

  target_include_directories(exago SYSTEM PRIVATE ${mpi4py_INCLUDE_DIRS})

  string(REGEX REPLACE ";" "\\\\;" _pytest_VALUES ${_pytest_VALUES})
  string(REGEX REPLACE "\n" ";" _pytest_VALUES ${_pytest_VALUES})
  list(GET _pytest_VALUES 0 pytest_VERSION)
  list(GET _pytest_VALUES 0 pytest_INCLUDE_DIRS)

  # Make sure all directory separators are '/'
  string(REGEX REPLACE "\\\\" "/" pytest_INCLUDE_DIRS ${pytest_INCLUDE_DIRS})

  # Get the major and minor version numbers
  string(REGEX REPLACE "\\." ";" _pytest_VERSION_LIST ${pytest_VERSION})
  list(GET _pytest_VERSION_LIST 0 pytest_VERSION_MAJOR)
  list(GET _pytest_VERSION_LIST 1 pytest_VERSION_MINOR)
  list(GET _pytest_VERSION_LIST 2 pytest_VERSION_PATCH)
  string(REGEX MATCH "[0-9]*" pytest_VERSION_PATCH ${pytest_VERSION_PATCH})
  math(
    EXPR
    pytest_VERSION_DECIMAL
    "(${pytest_VERSION_MAJOR} * 10000) + (${pytest_VERSION_MINOR} * 100) + ${pytest_VERSION_PATCH}"
  )
endif()

target_include_directories(exago SYSTEM PRIVATE ${pytest_INCLUDE_DIRS})

# Install Python bindings to same installation path:
install(TARGETS exago DESTINATION ${EXAGO_PYTHON_SITELIB})
