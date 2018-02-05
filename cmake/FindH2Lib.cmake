# - Try to find the h2lib library
# Define a path using H2Lib_ROOT
# Once done this will define
#
#  H2Lib_FOUND - system has h2lib
#  H2Lib_INCLUDE_DIRS - the h2lib include directory
#  H2Lib_LIBRARIES - Link these to use h2lib
#  H2Lib_LINKER_FLAGS - add these to your CMAKE_EXE_LINKER_FLAGS
#

find_path(H2Lib_INCLUDE_DIRS h2matrix.h
            PATH_SUFFIXES Library
            PATHS ${H2Lib_ROOT} 
            )
find_library(H2Lib_LIBRARIES libh2.a
            PATHS ${H2Lib_ROOT} 
            )

if (H2Lib_INCLUDE_DIRS AND H2Lib_LIBRARIES)
    message("-- Found the H2Lib library: ${H2Lib_LIBRARIES}")
    set(H2Lib_FOUND on)
else()
    set(H2Lib_FOUND off)
endif()

if (H2Lib_FOUND)
    find_package(Cairo)
    if (Chairo_FOUND)
        list(APPEND H2Lib_LIBRARIES ${CAIRO_LIBRARIES})
        list(APPEND H2Lib_INCLUDES ${CAIRO_INCLUDE_DIRS})
    endif()

    find_package(GLUT)
    if (GLUT_FOUND)
        list(APPEND H2Lib_LIBRARIES ${GLUT_LIBRARIES})
        list(APPEND H2Lib_INCLUDES ${GLUT_INCLUDE_DIR})
    endif()

    find_package(OpenGL)
    if (OpenGL_FOUND)
        list(APPEND H2Lib_LIBRARIES ${OPENGL_LIBRARIES})
        list(APPEND H2Lib_INCLUDES ${OPENGL_INCLUDE_DIR})
    endif()

    find_package(LAPACK REQUIRED)
    list(APPEND H2Lib_LIBRARIES ${LAPACK_LIBRARIES})
    set(H2Lib_LINKER_FLAGS ${LAPACK_LINKER_FLAGS})
endif()

# Hide advanced variables from CMake GUIs
MARK_AS_ADVANCED(H2Lib_LIBRARIES H2Lib_INCLUDE_DIRS)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(H2Lib DEFAULT_MSG H2Lib_LIBRARIES H2Lib_INCLUDE_DIRS H2Lib_FOUND)
