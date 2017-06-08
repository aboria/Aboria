macro(aboria_cxx_test target target_file)
    set(_test_files "")
    set(list_var "${ARGN}")
    foreach(filename IN LISTS list_var)
        set(_test_files ${_test_files} ${CMAKE_CURRENT_SOURCE_DIR}/${filename})
    endforeach()
    message("\t${_test_files}")
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

option(Aboria_RUN_TESTS "run main test suite with ctest" ON)
if (Aboria_RUN_TESTS)
    set(all_test_suites
        MetafunctionsTest 
        NeighboursTest
        BDTest
        SPHTest
        MDTest
        RbfInterpolationTest
        RbfPdeTest
        RbfStokesTest
        DiffusionAroundSpheres
        ParticleContainerTest
        SymbolicTest
        VariablesTest
        ConstructorsTest
        OperatorsTest
        SpatialUtilsTest
        ChebyshevTest
        FMMTest
        H2Test
        )
    foreach(test_suite ${all_test_suites})
        option(Aboria_RUN_TEST_${test_suite} "run ${test_suite} test suite with ctest" ON)
        if (Aboria_RUN_TEST_${test_suite})
            set(test_suites ${test_suites} ${test_suite})
        endif()
    endforeach()
endif()

option(Aboria_RUN_BENCHMARKS "run benchmarks during testing" OFF)
if (Aboria_RUN_BENCHMARKS)
    find_package(Gromacs)
    if (GROMACS_FOUND)
        add_definitions(-DHAVE_GROMACS)
        list(APPEND Aboria_LIBRARIES "${GROMACS_LIBRARIES}")
        include_directories(${GROMACS_INCLUDE_DIRS})
    endif()
    set(test_suites ${test_suites} SpeedTest)
endif()

set(SpeedTestFile speed_test.h)
set(SpeedTest 
    test_vector_addition
    test_daxpy
    test_finite_difference
    test_multiquadric
    test_multiquadric_scaling
    test_linear_spring
    )

set(MetafunctionsTestFile metafunctions.h) 
set(MetafunctionsTest
    test_aliasing
    test_expr_matching
    test_get_labels
    test_thrust
    ) 

set(SpatialUtilsTestFile spatial_utils.h)
set(SpatialUtilsTest
    test_bucket_indicies
    )

set(SymbolicTestFile symbolic.h)
set(SymbolicTest
    test_default
    )

set(VariablesTestFile variables.h)
set(VariablesTest
    test_std_vector
    )

set(OperatorsTestFile operators.h)
set(OperatorsTest
    test_Eigen
    test_Eigen_block
    test_documentation
    )

set(ConstructorsTestFile constructors.h)
set(ConstructorsTest
    test_std_vector
    )

set(ConstructorsTestCuda
    test_thrust_vector
    )

set(NeighboursTestFile neighbours.h)
set(NeighboursTest
    test_std_vector_bucket_search_serial
    test_std_vector_bucket_search_parallel
    test_std_vector_nanoflann_adaptor
    test_std_vector_octtree
    test_documentation
    )

set(ParticleContainerTestFile particle_container.h)
set(ParticleContainerTest 
    test_std_vector_bucket_search_serial
    test_std_vector_bucket_search_parallel
    test_documentation
    )

set(BDTestFile bd.h)
set(BDTest
    test_bucket_search_serial
    test_bucket_search_parallel
    test_nanoflann_adaptor
    test_octtree
    )

set(SPHTestFile sph.h)
set(SPHTest
    test_bucket_search_serial
    test_bucket_search_parallel
    test_nanoflann_adaptor
    )
 
set(MDTestFile md.h)
set(MDTest
    test_bucket_search_serial
    test_bucket_search_parallel
    test_nanoflann_adaptor
    test_octtree
    )
 
set(RbfInterpolationTestFile rbf_interpolation.h)
set(RbfInterpolationTest 
    test_bucket_search_serial
    test_bucket_search_parallel
    )
    
set(RbfPdeTestFile rbf_pde.h)
set(RbfPdeTest 
    test_bucket_search_serial
    test_bucket_search_parallel
    )

set(RbfStokesTestFile rbf_stokes.h)
set(RbfStokesTest 
    test_bucket_search_serial
    test_bucket_search_parallel
    )

set(DiffusionAroundSpheresFile diffusion_around_spheres.h)
set(DiffusionAroundSpheres 
    test_bucket_search_serial
    test_bucket_search_parallel
    test_nanoflann_adaptor
    )

set(ChebyshevTestFile chebyshev.h)
set(ChebyshevTest
    test_Rn_calculation
    test_fast_methods_bucket_search_serial
    test_fast_methods_bucket_search_parallel
    test_fast_methods_kd_tree
    test_fast_methods_octtree
    test_fmm_operators
    test_chebyshev_polynomial_calculation
    )


set(FMMTestFile fmm.h)
set(FMMTest
    test_fast_methods_bucket_search_serial
    test_fast_methods_bucket_search_parallel
    test_fast_methods_kd_tree
    test_fast_methods_octtree
    test_fmm_operators
    )

set(H2TestFile h2.h)
set(H2Test
    test_fast_methods_bucket_search_serial
    test_fast_methods_bucket_search_parallel
    test_fast_methods_kd_tree
    test_fast_methods_octtree
    )






message("-- Adding test suites:")
foreach(test_suite ${test_suites})

    aboria_cxx_test(${test_suite} ${test_suite}.cpp ${${test_suite}File})
    target_link_libraries(${test_suite} ${VTK_LIBRARIES} ${Boost_LIBRARIES} ${CUDA_LIBRARIES} ${GPERFTOOLS_LIBRARIES} ${Aboria_LIBRARIES})
    if (Aboria_THRUST_TEST_DEVICE_VECTOR)
        #aboria_cxx_test(tests_cuda tests_cuda.cu ${test_files})
        aboria_cxx_test(${test_suite}_cuda ${test_suite}.cu ${${test_suite}File})
        target_link_libraries(${test_suite}_cuda ${VTK_LIBRARIES} ${Boost_LIBRARIES} ${CUDA_LIBRARIES} ${Aboria_LIBRARIES})
    endif()

    foreach(test_name ${${test_suite}})
        add_test(NAME ${test_suite}:__${test_name} COMMAND ${RUNTIME_OUTPUT_DIRECTORY}/${test_suite} ${test_suite} ${test_name})
    endforeach()
    foreach(test_name ${${test_suite}Cuda})
        if (Aboria_THRUST_TEST_DEVICE_VECTOR)
            add_test(NAME ${test_suite}:__${test_name} COMMAND ${RUNTIME_OUTPUT_DIRECTORY}/${test_suite}_cuda ${test_suite} ${test_name})
        endif()
    endforeach()
endforeach()

