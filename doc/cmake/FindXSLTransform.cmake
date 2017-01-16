
# Find xsltproc to transform XML documents via XSLT
find_program(XSLTPROC_EXECUTABLE xsltproc DOC "xsltproc transforms XML via XSLT")
if(NOT XSLTPROC_EXECUTABLE)
  message(SEND_ERROR "xsltproc could not be found!")
endif(NOT XSLTPROC_EXECUTABLE)

set(XSLTPROC_FLAGS "--xinclude" CACHE STRING 
  "Flags to pass to xsltproc to transform XML documents")

find_package(ParseArguments REQUIRED)

# Transforms the source XML file by applying the given XSL stylesheet.
#
#   xsl_transform(output input [input2 input3 ...]
#                 STYLESHEET stylesheet
#                 [CATALOG catalog]
#                 [DIRECTORY mainfile]
#                 [PARAMETERS param1=value1 param2=value2 ...]
#                 [[MAKE_ALL_TARGET | MAKE_TARGET] target]
#                 [COMMENT comment])
#
# This macro builds a custom command that transforms an XML file
# (input) via the given XSL stylesheet. The output will either be a
# single file (the default) or a directory (if the DIRECTION argument
# is specified). The STYLESHEET stylesheet must be a valid XSL
# stylesheet. Any extra input files will be used as additional
# dependencies for the target. For example, these extra input files
# might refer to other XML files that are included by the input file
# through XInclude.
#
# When the XSL transform output is going to a directory, the mainfile
# argument provides the name of a file that will be generated within
# the output directory. This file will be used for dependency tracking.
# 
# XML catalogs can be used to remap parts of URIs within the
# stylesheet to other (typically local) entities. To provide an XML
# catalog file, specify the name of the XML catalog file via the
# CATALOG argument. It will be provided to the XSL transform.
# 
# The PARAMETERS argument is followed by param=value pairs that set
# additional parameters to the XSL stylesheet. The parameter names
# that can be used correspond to the <xsl:param> elements within the
# stylesheet.
# 
# To associate a target name with the result of the XSL
# transformation, use the MAKE_TARGET or MAKE_ALL_TARGET option and
# provide the name of the target. The MAKE_ALL_TARGET option only
# differs from MAKE_TARGET in that MAKE_ALL_TARGET will make the
# resulting target a part of the default build.
#
# If a COMMENT argument is provided, it will be used as the comment
# CMake provides when running this XSL transformation. Otherwise, the
# comment will be "Generating "output" via XSL transformation...".
macro(xsl_transform OUTPUT INPUT)
  parse_arguments(THIS_XSL
    "STYLESHEET;CATALOG;MAKE_ALL_TARGET;MAKE_TARGET;PARAMETERS;DIRECTORY;COMMENT"
    ""
    ${ARGN}
    )
  
  # TODO: Is this the best way to handle catalogs? The alternative is
  # that we could provide explicit remappings to the xsl_transform
  # macro, and it could generate a temporary XML catalog file.
  if (THIS_XSL_CATALOG)
    if(WIN32)
      set(THIS_XSL_CATALOG "")
    else(WIN32)
      set(THIS_XSL_CATALOG "XML_CATALOG_FILES=${THIS_XSL_CATALOG}")
    endif(WIN32)
  endif ()

  # Translate XSL parameters into a form that xsltproc can use.
  set(THIS_XSL_EXTRA_FLAGS)
  foreach(PARAM ${THIS_XSL_PARAMETERS})
    string(REGEX REPLACE "([^=]*)=([^;]*)" "\\1;\\2"
      XSL_PARAM_LIST ${PARAM})
    list(GET XSL_PARAM_LIST 0 XSL_PARAM_NAME)
    list(GET XSL_PARAM_LIST 1 XSL_PARAM_VALUE)
    if(${XSL_PARAM_VALUE} MATCHES "^\"([^\"]*)\"$")
      set(XSL_PARAM_VALUE ${CMAKE_MATCH_1})
    endif()

    list(APPEND THIS_XSL_EXTRA_FLAGS 
      --stringparam ${XSL_PARAM_NAME} ${XSL_PARAM_VALUE})
  endforeach(PARAM)

  # If the user didn't provide a comment for this transformation,
  # create a default one.
  if(NOT THIS_XSL_COMMENT)
    set(THIS_XSL_COMMENT "Generating ${OUTPUT} via XSL transformation...")
  endif()

  # Figure out the actual output file that we tell CMake about
  # (THIS_XSL_OUTPUT_FILE) and the output file or directory that we
  # tell xsltproc about (THIS_XSL_OUTPUT).
  if (THIS_XSL_DIRECTORY)
    set(THIS_XSL_OUTPUT_FILE ${OUTPUT}/${THIS_XSL_DIRECTORY})
    set(THIS_XSL_OUTPUT      ${OUTPUT}/)
  else()
    set(THIS_XSL_OUTPUT_FILE ${OUTPUT})
    set(THIS_XSL_OUTPUT      ${OUTPUT})
  endif()

  if(NOT THIS_XSL_STYLESHEET)
    message(SEND_ERROR 
      "xsl_transform macro invoked without a STYLESHEET argument")
  else()
  
    # Run the XSLT processor to do the XML transformation.
    add_custom_command(OUTPUT ${THIS_XSL_OUTPUT_FILE}
      COMMAND ${THIS_XSL_CATALOG} ${XSLTPROC_EXECUTABLE} ${XSLTPROC_FLAGS} 
              ${THIS_XSL_EXTRA_FLAGS} -o ${THIS_XSL_OUTPUT} 
              --path ${CMAKE_CURRENT_BINARY_DIR}
              ${THIS_XSL_STYLESHEET} ${INPUT}
      COMMENT ${THIS_XSL_COMMENT}
      DEPENDS ${INPUT} ${THIS_XSL_DEFAULT_ARGS}
      VERBATIM )
    set_source_files_properties(${THIS_XSL_OUTPUT_FILE}
      PROPERTIES GENERATED TRUE)
    # Create a custom target to refer to the result of this
    # transformation.
    if (THIS_XSL_MAKE_ALL_TARGET)
      add_custom_target(${THIS_XSL_MAKE_ALL_TARGET} ALL
        DEPENDS ${THIS_XSL_OUTPUT_FILE})
    elseif(THIS_XSL_MAKE_TARGET)
      add_custom_target(${THIS_XSL_MAKE_TARGET}
        DEPENDS ${THIS_XSL_OUTPUT_FILE})
      set_target_properties(${THIS_XSL_MAKE_TARGET}
        PROPERTIES
        EXCLUDE_FROM_ALL ON)
    endif()
  endif()
endmacro(xsl_transform)
