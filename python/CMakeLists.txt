if (${CMAKE_PROJECT_NAME} STREQUAL "Project")
    message(FATAL_ERROR "cmake should be executed not for 'python' subdirectory, but for the top directory of ChiQ.")
endif (${CMAKE_PROJECT_NAME} STREQUAL "Project")

if(CHIQ_WHEEL_BUILD)
  # Wheel: pure-Python files come from scikit-build-core wheel.packages; the CLI
  # from pyproject entry points. The extension is installed by src/CMakeLists.txt
  # to DESTINATION chiq.
else()
  file(GLOB scripts LIST_DIRECTORIES false scripts/*.py)
  install(FILES ${scripts}
   PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ
   DESTINATION ${CMAKE_INSTALL_BINDIR})
  install(DIRECTORY package/chiq DESTINATION ${CMAKE_INSTALL_LIBDIR}/bse-python)
  # Remove any stale symlink before installing the real backward-compat packages.
  install(CODE "file(REMOVE \"\$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/bse-python/bse\")")
  install(DIRECTORY package/bse        DESTINATION ${CMAKE_INSTALL_LIBDIR}/bse-python)
  install(DIRECTORY package/bse_solver DESTINATION ${CMAKE_INSTALL_LIBDIR}/bse-python)
endif()
