add_library(${DIR_IO} STATIC HDFDriver.cpp)

# hdf5
if (HDF5_FOUND)
    target_link_libraries(${DIR_IO} ${HDF5_LIBRARIES})
endif ()

install(TARGETS ${DIR_IO} ARCHIVE DESTINATION ${LIB_DIR})
