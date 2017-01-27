find_package(ParseArguments REQUIRED)
find_package(XSLTransform REQUIRED)
find_package(BoostBook REQUIRED)
find_package(FOP)



# Transform Quickbook into BoostBook XML
macro(add_documentation INPUT)
  parse_arguments(DOCUMENTATION "IMAGES;PARAMETERS" "" ${ARGN})

  # If INPUT is not a full path, it's in the current source directory.
  get_filename_component(INPUT_PATH ${INPUT} PATH)
  if(INPUT_PATH STREQUAL "")
    set(INPUT_PATH "${CMAKE_CURRENT_SOURCE_DIR}/${INPUT}")
  else(INPUT_PATH STREQUAL "")
    set(INPUT_PATH ${INPUT})
  endif(INPUT_PATH STREQUAL "")

  set(QBK_FILE ${CMAKE_CURRENT_BINARY_DIR}/${THIS_PROJECT_NAME}.qbk)
  set(DBK_FILE ${CMAKE_CURRENT_BINARY_DIR}/${THIS_PROJECT_NAME}.docbook)
  set(FOP_FILE ${CMAKE_CURRENT_BINARY_DIR}/${THIS_PROJECT_NAME}.fo)
  set(PDF_FILE ${CMAKE_CURRENT_BINARY_DIR}/${THIS_PROJECT_NAME}.pdf)

  # copy to destination directory because quickbook screws up xinclude paths 
  # when the output is not in the source directory
  add_custom_command(OUTPUT ${QBK_FILE}
    COMMAND ${CMAKE_COMMAND} -E copy ${INPUT_PATH} ${QBK_FILE}
    DEPENDS ${INPUT_PATH})

  # copy all dependencies that are not built
  set(DEPENDENCIES)
  foreach(file ${DOCUMENTATION_DEFAULT_ARGS})
    set(srcfile ${CMAKE_CURRENT_SOURCE_DIR}/${file})
    set(binfile ${CMAKE_CURRENT_BINARY_DIR}/${file})
    if(EXISTS ${srcfile})
      add_custom_command(OUTPUT ${binfile}
        COMMAND ${CMAKE_COMMAND} -E copy ${srcfile} ${binfile}
        DEPENDS ${srcfile})
    endif(EXISTS ${srcfile})

    set(DEPENDENCIES ${DEPENDENCIES} ${binfile})
  endforeach(file ${DOCUMENTATION_DEFAULT_ARGS})

  # copy images into different locations... quick and dirty...
  file(COPY images DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  file(COPY images DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/html)
  file(COPY css DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/html)
  file(COPY images DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/html-website)
  file(COPY css DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/html-website)

  quickbook_to_docbook(${DBK_FILE} ${QBK_FILE} ${DEPENDENCIES})

  xsl_transform(${CMAKE_CURRENT_BINARY_DIR}/html-website ${DBK_FILE}
    STYLESHEET ${CMAKE_CURRENT_BINARY_DIR}/html-website.xsl
    CATALOG ${BOOSTBOOK_CATALOG}
    PARAMETERS ${DOCUMENTATION_PARAMETERS}
    DIRECTORY HTML-WEBSITE.manifest
    COMMENT "Generating HTML documentation for ${THIS_PROJECT_NAME}."
    MAKE_TARGET ${THIS_PROJECT_NAME}-html-website)

 xsl_transform(${CMAKE_CURRENT_BINARY_DIR}/html ${DBK_FILE}
     #STYLESHEET ${BOOSTBOOK_XSL_DIR}/html.xsl
    STYLESHEET ${CMAKE_CURRENT_BINARY_DIR}/html.xsl
    CATALOG ${BOOSTBOOK_CATALOG}
    PARAMETERS ${DOCUMENTATION_PARAMETERS}
    DIRECTORY HTML.manifest
    COMMENT "Generating HTML documentation for ${THIS_PROJECT_NAME}."
    MAKE_TARGET ${THIS_PROJECT_NAME}-html)

  xsl_transform(${CMAKE_CURRENT_BINARY_DIR}/man ${DBK_FILE}
    STYLESHEET ${BOOSTBOOK_XSL_DIR}/manpages.xsl
    CATALOG ${BOOSTBOOK_CATALOG}
    DIRECTORY man.manifest
    COMMENT "Generating man pages for ${THIS_PROJECT_NAME}."
    MAKE_TARGET ${THIS_PROJECT_NAME}-man)

  if(FOP_FOUND)
    xsl_transform(${FOP_FILE} ${DBK_FILE}
      STYLESHEET ${BOOSTBOOK_XSL_DIR}/fo.xsl
      CATALOG ${BOOSTBOOK_CATALOG}
      PARAMETERS img.src.path=${CMAKE_CURRENT_BINARY_DIR}/images/
      MAKE_TARGET ${THIS_PROJECT_NAME}-fo)

    add_custom_command(OUTPUT ${PDF_FILE}
      COMMAND ${FOP_EXECUTABLE} ${FOP_FILE} ${PDF_FILE} 2>/dev/null
      DEPENDS ${FOP_FILE}
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      )
    add_custom_target(${THIS_PROJECT_NAME}-pdf DEPENDS ${PDF_FILE})
    set_target_properties(${THIS_PROJECT_NAME}-pdf
      PROPERTIES EXCLUDE_FROM_ALL ON)
  else(FOP_FOUND)
    message(STATUS "Install Apache FOP to build PDF documentation.")
  endif(FOP_FOUND)

endmacro(add_documentation INPUT)
