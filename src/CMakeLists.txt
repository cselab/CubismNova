# CubismNova library components

# utilities
if (CUBISM_UTIL)
    set(DIR_UTIL Util)
    add_subdirectory(${DIR_UTIL})
endif ()

# IO
if (CUBISM_IO)
    set(DIR_IO IO)
    add_subdirectory(${DIR_IO})
endif ()
