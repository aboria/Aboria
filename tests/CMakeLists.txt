macro(aboria_cxx_test target target_file)
    set(_test_files "")
    set(list_var "${ARGN}")
    foreach(filename IN LISTS list_var)
        set(_test_files ${_test_files} ${CMAKE_CURRENT_SOURCE_DIR}/${filename})
    endforeach()
    message(${_test_files})
    get_filename_component(CPP_FILE_NAME ${target_file} NAME)
    set(CPP_FULL_NAME "${CMAKE_CURRENT_BINARY_DIR}/${CPP_FILE_NAME}")
    add_custom_command(
        OUTPUT "${CPP_FULL_NAME}"
        COMMAND ${CXXTESTGEN} --runner=ErrorPrinter --output "${CPP_FULL_NAME}" ${_test_files}
        DEPENDS ${_test_files}
        VERBATIM
    )
if (${target_file} MATCHES .cu$)
    cuda_add_executable(${target} ${CPP_FULL_NAME})
else()
    add_executable(${target} ${CPP_FULL_NAME})
endif()
    set_target_properties(${target} PROPERTIES COMPILE_FLAGS "-Wno-effc++")
endmacro(aboria_cxx_test)

set(test_files
    operators.h
    constructors.h 
    variables.h
    particle_container.h
    neighbours.h
    diffusion_around_spheres.h
    symbolic.h
    rbf_pde.h
    rbf_interpolation.h
    metafunctions.h
    dem.h
    sph.h
    bd.h
    md.h
    linear_algebra_compare.h
    speed_test.h
    )

set(test_suites
    NeighboursTest
    BDTest
    SPHTest
    MDTest
    RbfInterpolationTest 
    RbfPdeTest 
    DiffusionAroundSpheres 
    ParticleContainerTest 
    SymbolicTest 
    VariablesTest
    ConstructorsTest
    OperatorsTest 
    SpeedTest
    )


#set(test_files metafunctions.h)
#set(test_suites )
set(test_files speed_test.h)
set(test_suites SpeedTest)

set(SpeedTest 
    #test_aboria
    #test_eigen
    test_vector_addition
    )

set(SymbolicTest
    test_default
    )

set(VariablesTest
    test_std_vector
    )

set(OperatorsTest
    test_Eigen
    )

set(ConstructorsTest
    test_std_vector
    )

set(ConstructorsTestCuda
    test_thrust_vector
    )

set(NeighboursTest
    test_std_vector_bucket_search_serial
    test_std_vector_bucket_search_parallel
    )

set(ParticleContainerTest 
    test_std_vector_bucket_search_serial
    test_std_vector_bucket_search_parallel
    )

set(BDTest
    test_bucket_search_serial
    test_bucket_search_parallel
    )

set(SPHTest
    test_bucket_search_serial
    test_bucket_search_parallel
    )
 
set(MDTest
    test_bucket_search_serial
    test_bucket_search_parallel
    )
 
set(RbfInterpolationTest 
    test_bucket_search_serial
    test_bucket_search_parallel
    )
    
set(RbfPdeTest 
    test_bucket_search_serial
    test_bucket_search_parallel
    )

set(DiffusionAroundSpheres 
    test_bucket_search_serial
    test_bucket_search_parallel
    )

aboria_cxx_test(tests tests.cpp ${test_files})
target_link_libraries(tests ${VTK_LIBRARIES} ${Boost_LIBRARIES} ${CUDA_LIBRARIES} ${GPERFTOOLS_LIBRARIES})
if (Aboria_THRUST_TEST_DEVICE_VECTOR)
    aboria_cxx_test(tests_cuda tests_cuda.cu ${test_files})
    target_link_libraries(tests_cuda ${VTK_LIBRARIES} ${Boost_LIBRARIES} ${CUDA_LIBRARIES} )
endif()

foreach(test_suite ${test_suites})
    foreach(test_name ${${test_suite}})
        add_test(NAME ${test_suite}:__${test_name} COMMAND ${RUNTIME_OUTPUT_DIRECTORY}/tests ${test_suite} ${test_name})
    endforeach()
    foreach(test_name ${${test_suite}Cuda})
        if (Aboria_THRUST_TEST_DEVICE_VECTOR)
            add_test(NAME ${test_suite}:__${test_name} COMMAND ${RUNTIME_OUTPUT_DIRECTORY}/tests ${test_suite} ${test_name})
        endif()
    endforeach()
endforeach()

