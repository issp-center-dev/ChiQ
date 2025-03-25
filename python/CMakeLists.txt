if (${CMAKE_PROJECT_NAME} STREQUAL "Project")
    message(FATAL_ERROR "cmake should be executed not for 'python' subdirectory, but for the top directory of ChiQ.")
endif (${CMAKE_PROJECT_NAME} STREQUAL "Project")

file(GLOB scripts LIST_DIRECTORIES false scripts/*.py)
install(FILES ${scripts}
 PERMISSIONS # 755
 OWNER_EXECUTE OWNER_WRITE OWNER_READ
 GROUP_EXECUTE GROUP_READ
 WORLD_EXECUTE WORLD_READ
 DESTINATION ${CMAKE_INSTALL_BINDIR})
install(DIRECTORY package/chiq DESTINATION ${CMAKE_INSTALL_LIBDIR}/bse-python)
# make bse package for backward compatibility
install(CODE "execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/bse-python/chiq ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/bse-python/bse)")
