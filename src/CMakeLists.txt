# include guard
if (${CMAKE_PROJECT_NAME} STREQUAL "Project")
    message(FATAL_ERROR "cmake should be executed not for 'src' subdirectory, but for the top directory of ChiQ.")
endif (${CMAKE_PROJECT_NAME} STREQUAL "Project")

# include_directories(./)
# set(SOURCES main.cpp set_initial.cpp)
# add_executable(BSE ${SOURCES})
# target_link_libraries(BSE ${CMAKE_THREAD_LIBS_INIT})
#
# install(TARGETS BSE DESTINATION bin)

pybind11_add_module(_bse_solver bse_solver_pybind.cpp)

if(CHIQ_WHEEL_BUILD)
  install(TARGETS _bse_solver LIBRARY DESTINATION chiq)
else()
  install(SCRIPT ${CHIQ_LEGACY_CLEANUP_SCRIPT})
  install(TARGETS _bse_solver
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/bse-python/chiq)
endif()
