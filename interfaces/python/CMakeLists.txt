add_subdirectory(typedefs)

configure_file(exago_config.py.in ${PROJECT_BINARY_DIR}/exago_config.py)

set(PYTHON_INTERFACE ${PROJECT_BINARY_DIR}/exago_config.py __init__.py pflow.py
                     opflow.py
)

install(FILES ${PYTHON_INTERFACE} DESTINATION ${EXAGO_PYTHON_SITELIB}/exago)
