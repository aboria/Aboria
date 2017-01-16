################################################################################
# Copyright (c) 2010 Daniel Pfeifer                                            #
################################################################################

find_program(FOP_EXECUTABLE fop)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FOP DEFAULT_MSG FOP_EXECUTABLE)
