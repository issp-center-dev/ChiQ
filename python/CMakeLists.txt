if (${CMAKE_PROJECT_NAME} STREQUAL "Project")
    message(FATAL_ERROR "cmake should be executed not for 'python' subdirectory, but for the top directory of ChiQ.")
endif (${CMAKE_PROJECT_NAME} STREQUAL "Project")

if(CHIQ_WHEEL_BUILD)
  # Wheel: pure-Python files come from scikit-build-core wheel.packages; the CLI
  # from pyproject entry points. The extension is installed by src/CMakeLists.txt
  # to DESTINATION chiq.
else()
  file(GLOB scripts LIST_DIRECTORIES false scripts/*.py)
  foreach(script ${scripts})
    get_filename_component(command_name "${script}" NAME_WE)
    set(deprecated_wrapper
      "${CMAKE_CURRENT_BINARY_DIR}/deprecated-scripts/${command_name}.py")
    file(GENERATE OUTPUT "${deprecated_wrapper}" CONTENT
"#!/usr/bin/env python3
from chiq.cli._deprecated import ${command_name}_py

if __name__ == \"__main__\":
    ${command_name}_py()
")
    install(PROGRAMS "${deprecated_wrapper}"
      DESTINATION ${CMAKE_INSTALL_BINDIR})
    install(PROGRAMS "${script}"
      DESTINATION ${CMAKE_INSTALL_BINDIR}
      RENAME "${command_name}")
  endforeach()
  set(_chiq_install_excludes
    PATTERN "__pycache__" EXCLUDE
    PATTERN "*.pyc" EXCLUDE
    PATTERN "*.so" EXCLUDE
    PATTERN "*.pyd" EXCLUDE
    PATTERN "*.dylib" EXCLUDE)
  install(DIRECTORY package/chiq DESTINATION ${CMAKE_INSTALL_LIBDIR}/bse-python
    ${_chiq_install_excludes})
  # Remove only the exact stale symlink before installing the real package.
  install(CODE "
    set(_obsolete_bse \"\$ENV{DESTDIR}${CMAKE_INSTALL_FULL_LIBDIR}/bse-python/bse\")
    if(IS_SYMLINK \"\${_obsolete_bse}\")
      file(REMOVE \"\${_obsolete_bse}\")
    endif()
  ")
  install(DIRECTORY package/bse DESTINATION ${CMAKE_INSTALL_LIBDIR}/bse-python
    ${_chiq_install_excludes})
  install(DIRECTORY package/bse_solver DESTINATION ${CMAKE_INSTALL_LIBDIR}/bse-python
    ${_chiq_install_excludes})
endif()
