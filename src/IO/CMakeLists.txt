set(LIBIO "Cubism${DIR_IO}" CACHE INTERNAL "Cubism library archives")
add_library(${LIBIO} STATIC HDFDriver.cpp)

# hdf5
if (HDF5_FOUND)
    target_include_directories(${LIBIO} PRIVATE ${HDF5_INCLUDE_DIRS})
    target_link_libraries(${LIBIO} ${HDF5_LIBRARIES})
endif ()

install(TARGETS ${LIBIO} ARCHIVE DESTINATION ${LIB_DIR})
