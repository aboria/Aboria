function(aboria_cxx_test target source ext)
    get_filename_component(CPP_FILE_NAME ${source} NAME)
    string(REGEX REPLACE "h$|hpp$" ${ext} CPP_FILE_NAME ${CPP_FILE_NAME})
    message(${CPP_FILE_NAME})
    set(CPP_FULL_NAME "${CMAKE_CURRENT_BINARY_DIR}/${CPP_FILE_NAME}")
    add_custom_command(
        OUTPUT "${CPP_FULL_NAME}"
        COMMAND ${CXXTESTGEN} --runner=ErrorPrinter --output "${CPP_FULL_NAME}" "${source}"
        DEPENDS "${source}"
    )
if (${ext} MATCHES cu)
    cuda_add_executable(${target} ${CPP_FULL_NAME})
else()
    add_executable(${target} ${CPP_FULL_NAME})
endif()
    set_target_properties(${target} PROPERTIES COMPILE_FLAGS "-Wno-effc++")
    add_test(${target} ${RUNTIME_OUTPUT_DIRECTORY}/${target})
endfunction(aboria_cxx_test)

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
#    rbf_interpolation
    neighbours
    )

foreach(test ${test_package})
    aboria_cxx_test(${test} ${CMAKE_CURRENT_SOURCE_DIR}/${test}.h cpp)
    #cxxtest_add_test(${test} ${test}.cpp ${CMAKE_CURRENT_SOURCE_DIR}/${test}.h)
    target_link_libraries(${test} ${VTK_LIBRARIES} ${Boost_LIBRARIES})
    if (Aboria_THRUST_TEST_DEVICE_VECTOR)
        aboria_cxx_test(${test}_cuda ${CMAKE_CURRENT_SOURCE_DIR}/${test}.h cu)
        target_link_libraries(${test}_cuda ${VTK_LIBRARIES} ${Boost_LIBRARIES} ${CUDA_LIBRARIES})
    endif()
endforeach(test)
