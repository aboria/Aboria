
include(FindPackageHandleStandardArgs)

# Find the QuickBook executable
find_program(QUICKBOOK_EXECUTABLE quickbook
  PATHS $ENV{BOOST_ROOT}/dist/bin)
find_package_handle_standard_args(QUICKBOOK DEFAULT_MSG QUICKBOOK_EXECUTABLE)

# Find the DocBook DTD (version 4.2)
find_path(DOCBOOK_DTD_DIR docbookx.dtd
  PATHS /usr/share/xml/docbook/schema/dtd/4.2
  DOC "Path to the DocBook DTD")
find_package_handle_standard_args(DOCBOOK_DTD DEFAULT_MSG DOCBOOK_DTD_DIR)

# Find the DocBook XSL stylesheets
find_path(DOCBOOK_XSL_DIR html/html.xsl
  PATHS /usr/share/xml/docbook/stylesheet/docbook-xsl 
  DOC "Path to the DocBook XSL stylesheets")
find_package_handle_standard_args(DOCBOOK_XSL DEFAULT_MSG DOCBOOK_XSL_DIR)

#Find the BoostBook DTD
find_path(BOOSTBOOK_DTD_DIR boostbook.dtd
  PATHS /usr/share/boostbook/dtd $ENV{BOOST_ROOT}/dist/share/boostbook/dtd
  DOC "Path to the BoostBook DTD")
find_package_handle_standard_args(BOOSTBOOK_DTD DEFAULT_MSG BOOSTBOOK_DTD_DIR)

#Find the BoostBook XSL stylesheets
find_path(BOOSTBOOK_XSL_DIR docbook.xsl
  PATHS /usr/share/boostbook/xsl $ENV{BOOST_ROOT}/dist/share/boostbook/xsl
  DOC "Path to the BoostBook XSL stylesheets")
find_package_handle_standard_args(BOOSTBOOK_XSL DEFAULT_MSG BOOSTBOOK_XSL_DIR)

set(BOOSTBOOK_CATALOG ${CMAKE_BINARY_DIR}/boostbook_catalog.xml)
file(WRITE ${BOOSTBOOK_CATALOG}
  "<?xml version=\"1.0\"?>\n"
  "<!DOCTYPE catalog\n"
  "  PUBLIC \"-//OASIS/DTD Entity Resolution XML Catalog V1.0//EN\"\n"
  "  \"http://www.oasis-open.org/committees/entity/release/1.0/catalog.dtd\">\n"
  "<catalog xmlns=\"urn:oasis:names:tc:entity:xmlns:xml:catalog\">\n"
  "  <rewriteURI"
    " uriStartString=\"http://www.oasis-open.org/docbook/xml/4.2/\""
    " rewritePrefix=\"file://${DOCBOOK_DTD_DIR}/\""
    "/>\n"
  "  <rewriteURI"
    " uriStartString=\"http://docbook.sourceforge.net/release/xsl/current/\""
    " rewritePrefix=\"file://${DOCBOOK_XSL_DIR}/\""
    "/>\n"
  "  <rewriteURI"
    " uriStartString=\"http://www.boost.org/tools/boostbook/dtd/\""
    " rewritePrefix=\"file://${BOOSTBOOK_DTD_DIR}/\""
    "/>\n"
  "  <rewriteURI"
    " uriStartString=\"http://www.boost.org/tools/boostbook/xsl/\""
    " rewritePrefix=\"file://${BOOSTBOOK_XSL_DIR}/\""
    "/>\n"
  "</catalog>\n"
  )

# Transform Quickbook into BoostBook XML
macro(quickbook_to_boostbook OUTPUT INPUT)

  # If INPUT is not a full path, it's in the current source directory.
  get_filename_component(INPUT_PATH ${INPUT} PATH)
  if(INPUT_PATH STREQUAL "")
    set(INPUT_PATH "${CMAKE_CURRENT_SOURCE_DIR}/${INPUT}")
  else(INPUT_PATH STREQUAL "")
    set(INPUT_PATH ${INPUT})
  endif(INPUT_PATH STREQUAL "")

  add_custom_command(OUTPUT ${OUTPUT}
    COMMAND ${QUICKBOOK_EXECUTABLE} "--output-file=${OUTPUT}" ${INPUT_PATH}
    COMMENT "Generating BoostBook XML from ${INPUT}."
    DEPENDS ${INPUT} ${ARGN})

endmacro(quickbook_to_boostbook OUTPUT INPUT)

# Transform BoostBook XML into DocBook
macro(boostbook_to_docbook OUTPUT INPUT)
  xsl_transform(${OUTPUT} ${INPUT}
    STYLESHEET ${BOOSTBOOK_XSL_DIR}/docbook.xsl
    CATALOG ${BOOSTBOOK_CATALOG}
    COMMENT "Generating DocBook from ${INPUT}."
    DEPENDS ${INPUT} ${ARGN}
    )
endmacro(boostbook_to_docbook OUTPUT INPUT)

# Transform Quickbook into BoostBook XML, then BoostBook XML into DocBook
macro(quickbook_to_docbook OUTPUT INPUT)
  get_filename_component(XML_FILE ${INPUT} NAME_WE)
  set(XML_FILE ${XML_FILE}.xml)
  quickbook_to_boostbook(${XML_FILE} ${INPUT} ${ARGN})
  boostbook_to_docbook(${OUTPUT} ${XML_FILE} ${ARGN})
endmacro(quickbook_to_docbook OUTPUT INPUT)
