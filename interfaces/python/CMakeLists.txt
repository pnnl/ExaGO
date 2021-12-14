# New true python bindings
pybind11_add_module(exago opflow_python.cpp)
target_link_libraries(exago PRIVATE ExaGO::UTILS ExaGO::OPFLOW)

# Install Python bindings to same installation path:
install(TARGETS exago DESTINATION ${EXAGO_PYTHON_SITELIB})
