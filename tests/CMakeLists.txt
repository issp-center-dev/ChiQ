# Copy all *.h5 in database/
file(GLOB_RECURSE DATA_FILES "${CMAKE_SOURCE_DIR}/database/*.h5")
foreach(datafile ${DATA_FILES})
    configure_file(${datafile} ${CMAKE_BINARY_DIR}/tests/ COPYONLY)
endforeach()

add_subdirectory(c++)
# add_subdirectory(python)

# for python
file(COPY ${CMAKE_SOURCE_DIR}/tests/python/ DESTINATION ${CMAKE_BINARY_DIR}/tests/python)
configure_file(python/pytest.ini ${CMAKE_BINARY_DIR}/tests/python/pytest.ini)
configure_file(python/pytest.ini ${CMAKE_BINARY_DIR}/pytest.ini)
