function(aboria_cxx_test target source)
    get_filename_component(CPP_FILE_NAME ${source} NAME)
if (Aboria_USE_THRUST)
    string(REGEX REPLACE "h$|hpp$" "cu" CPP_FILE_NAME ${CPP_FILE_NAME})
else()
    string(REGEX REPLACE "h$|hpp$" "cpp" CPP_FILE_NAME ${CPP_FILE_NAME})
endif()
    message(${CPP_FILE_NAME})
    set(CPP_FULL_NAME "${CMAKE_CURRENT_BINARY_DIR}/${CPP_FILE_NAME}")
    add_custom_command(
        OUTPUT "${CPP_FULL_NAME}"
        COMMAND ${CXXTESTGEN} --runner=ErrorPrinter --output "${CPP_FULL_NAME}" "${source}"
        DEPENDS "${source}"
    )
if (Aboria_USE_THRUST)
    cuda_add_executable(${target} ${CPP_FULL_NAME})
else()
    add_executable(${target} ${CPP_FULL_NAME})
endif()
    set_target_properties(${target} PROPERTIES COMPILE_FLAGS "-Wno-effc++")
    add_test(${target} ${RUNTIME_OUTPUT_DIRECTORY}/${target})
endfunction(abori_cxx_test)

set(test_package 
    operators
    constructors 
    variables
    particle_container
    neighbours
    #geometry
    diffusion_around_spheres
    assemble
    symbolic
    rbf_pde
    rbf_interpolation
    metafunctions
    dem
    sph
    bd
    md
    )

set(test_package
    neighbours
    )

foreach(test ${test_package})
    aboria_cxx_test(${test} ${CMAKE_CURRENT_SOURCE_DIR}/${test}.h)
    #cxxtest_add_test(${test} ${test}.cpp ${CMAKE_CURRENT_SOURCE_DIR}/${test}.h)
    target_link_libraries(${test} ${VTK_LIBRARIES} ${Boost_LIBRARIES})
endforeach(test)
