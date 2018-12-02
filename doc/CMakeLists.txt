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
    ../src/CellListOrdered.h
    ../src/CellList.h
    ../src/NanoFlannAdaptor.h
    ../src/Kdtree.h
    ../src/OctTree.h
    ../src/NeighbourSearchBase.h
    ../src/Operators.h
    ../src/Chebyshev.h
    ../src/Kernels.h
    ../src/Functions.h
    ../src/Log.h
    ../src/Utils.h
    ../src/Symbolic.h
    ../src/Transform.h
    ../src/LatticeIterator.h
    ../src/SpatialUtil.h
    ../src/Traits.h
    ../src/Get.h
    ../src/Evaluate.h
    ../src/Search.h
    ../src/Vector.h
    ../src/FastMultipoleMethod.h
    ../src/H2Lib.h
    )

  add_reference(libaboria.xml ${ABORIA_HEADERS}
      DOXYGEN
      MACRO_EXPANSION=YES
      EXPAND_ONLY_PREDEF=YES
      EXTRACT_ALL=YES
      EXTRACT_PRIVATE=YES
      EXTRACT_STATIC=YES
      TAB_SIZE=4
      PREDEFINED="BOOST_COMPUTE_DOXYGEN_INVOKED\" \\ \"HAVE_EIGEN\" \\ \" HAVE_THRUST\" \\ \"HAVE_H2LIB\" \\ \"HAVE_VTK\" \\ \"CUDA_HOST_DEVICE=\" \\ \"ABORIA_HOST_DEVICE_IGNORE_WARN="
      EXPAND_AS_DEFINED=UFUNC
      EXPAND_AS_DEFINED="COMPOUND_ASSIGN\" \"COMPARISON\" \"OPERATOR\" \"UNARY_OPERATOR\" \"ABORIA_UNARY_FUNCTION\" \"ABORIA_BINARY_FUNCTION"
      BOOSTBOOK
      boost.doxygen.reftitle="libaboria")


    set(QBK_SOURCES
        aboria.qbk
        abstract.qbk
        introduction.qbk
        installation.qbk
        reference.qbk
        examples.qbk
        benchmarks.qbk
        symbolic.qbk
        )
    set(CPP_EXAMPLES
        ../tests/particle_container.h
        ../tests/operators.h
        ../tests/neighbours.h
        ../tests/data_structures.h
        ../tests/symbolic.h
        ../tests/md.h
        ../tests/example_container.h
        ../tests/bd.h
        ../tests/dem.h
        ../tests/fput.h
        ../tests/sph.h
        ../tests/rbf_interpolation.h
        ../tests/rbf_pde.h
        ../tests/speed_test.h
        ../tests/doc_getting_started.h
        ../tests/id_search.h
        ../tests/parallel.h
        )


    set (HTML_XSL_PRE "<?xml version='1.0'?>
<xsl:stylesheet xmlns:xsl=\"http://www.w3.org/1999/XSL/Transform\" version=\"1.0\">
<xsl:import href=\"${BOOSTBOOK_XSL_DIR}/html.xsl\"/>
<xsl:template name=\"user.head.content\">
<script type=\"text/x-mathjax-config\">
MathJax.Hub.Config({
  tex2jax: {
    inlineMath: [['$','$'], ['\\\\(','\\\\)']],
    processEscapes: true
    },
  TeX: {
    extensions: [\"AMSmath.js\", \"AMSsymbols.js\"]
    }
});
</script>
<script type=\"text/javascript\"
  src=\"https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS_CHTML\">
</script>
")
    set (HTML_XSL_GOOGLE_ANALYTICS "
<script>
  (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
  (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
  m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
  })(window,document,'script','https://www.google-analytics.com/analytics.js','ga');

  ga('create', 'UA-48149829-4', 'auto');
  ga('send', 'pageview');
</script>")

    set (HTML_XSL_POST "</xsl:template>
</xsl:stylesheet>")

   file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/html.xsl "${HTML_XSL_PRE}${HTML_XSL_POST}")
   file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/html-website.xsl "${HTML_XSL_PRE}${HTML_XSL_GOOGLE_ANALYTICS}${HTML_XSL_POST}")

    # PARAMETERS : http://docbook.sourceforge.net/release/xsl/current/doc/
  add_documentation(
  #  quickstart.cpp
    ${QBK_SOURCES}
    ${CPP_EXAMPLES}
    libaboria.xml
  PARAMETERS
    #boost.root=./
    html.stylesheet=css/boostbook.css
    boost.image.src=images/aboria2.jpg
    boost.image.alt="Aboria"
    boost.image.w=200
    boost.image.h=62
    boost.mathjax=1
    chunk.section.depth=1
    #<xsl:param>chunk.first.sections=1
    #toc.section.depth=3
    #toc.max.depth=2
    #generate.section.toc.level=4
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
    images/note.png
    images/neighbour/cell_lists.svg
    images/neighbour/octtree.svg
    images/fput/fput999.svg
    css/boostbook.css
  )

endif(CMAKE_SYSTEM_NAME STREQUAL Linux)

