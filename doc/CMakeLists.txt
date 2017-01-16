cmake_minimum_required(VERSION 2.6)

# Building of Documentation is tested on Linux only
if(CMAKE_SYSTEM_NAME STREQUAL Linux)

# adapted from http://github.com/purpleKarrot/CToolBox.git
  set(THIS_PROJECT_NAME aboria)
  set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})
  find_package(Reference)
  find_package(Documentation)

  set(ABORIA_HEADERS
    ../src/Aboria.h
    ../src/Particles.h
    ../src/Variable.h
    ../src/DataVector.h
    ../src/BucketSearchParallel.h
    ../src/Functions.h
    ../src/Log.h
    ../src/Symbolic.h
    ../src/Traits.h
    ../src/Get.h
    ../src/Assemble.h
    ../src/Operators.h
    ../src/Vector.h
    )

  add_reference(libaboria-ref.xml    ${ABORIA_HEADERS}   BOOSTBOOK boost.doxygen.reftitle="Aboria")
  
    set(QBK_SOURCES
        aboria.qbk
        installation.qbk
        neighbourhood.qbk
        container.qbk
        introduction.qbk
        reference.qbk
        examples.qbk
        linear_algebra.qbk
        symbolic.qbk
        )

    # PARAMETERS : http://docbook.sourceforge.net/release/xsl/current/doc/
  add_documentation(
  #  quickstart.cpp
    ${QBK_SOURCES}
    libaboria-ref.xml
  PARAMETERS
    #boost.root=./
    html.stylesheet=css/boostbook.css
    boost.image.src=images/aboria2.jpg
    boost.image.alt="Aboria"
    boost.image.w=200
    boost.image.h=62
    chunk.section.depth=4
    #<xsl:param>chunk.first.sections=1
    toc.section.depth=3
    toc.max.depth=2
    generate.section.toc.level=4
    nav.layout=none
    #admon.graphics.path=images/admonition/
    #navig.graphics.path=images/navigation/
    #callout.graphics.path=images/callout/
    #boost.doxygen.strip_from_path=.
    #img.src.path=../images/useless_directory/
  IMAGES
    images/aboria2.png
    images/next.png
    images/prev.png
    images/up.png
    images/down.png
    images/home.png
    css/boostbook.css
  )

endif(CMAKE_SYSTEM_NAME STREQUAL Linux)

