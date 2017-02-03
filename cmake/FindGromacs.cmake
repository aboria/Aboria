#
# FindGromacs
#
# This module finds the Gromacs library
#
# GROMACS_LIBRARIES -  libraries to link against
find_library( GROMACS_LIBRARIES 
            NAMES gromacs
            DOC "Gromacs libraries"
  )

# Check for required components
include( FindPackageHandleStandardArgs )
find_package_handle_standard_args( Gromacs
    REQUIRED_VARS GROMACS_LIBRARIES
  )

mark_as_advanced(GROMACS_LIBRARIES)
