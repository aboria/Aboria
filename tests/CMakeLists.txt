macro(aboria_cxx_test target target_file use_cuda_function)
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
if (${use_cuda_function})
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
        IDSearchTest
        ParticleContainerTest
        SymbolicTest
        VariablesTest
        ConstructorsTest
        OperatorsTest
        UtilsTest
        IteratorsTest
        ChebyshevTest
        FMMTest
        DocGettingStartedTest
        )
    set(notdebug_test_suites
        BDTest
        SPHTest
        MDTest
        MDLevel1Test
        DiffusionAroundSpheres
        RbfInterpolationTest
        RbfPdeTest
        #RbfStokesTest
        H2Test
        )

    foreach(test_suite ${all_test_suites})
        option(Aboria_RUN_TEST_${test_suite} "run ${test_suite} test suite with ctest" ON)
        if (Aboria_RUN_TEST_${test_suite})
            set(test_suites ${test_suites} ${test_suite})
        endif()
    endforeach()

    foreach(test_suite ${notdebug_test_suites})
        if (CMAKE_BUILD_TYPE STREQUAL "Debug")
            option(Aboria_RUN_TEST_${test_suite} "run ${test_suite} test suite with ctest" OFF)
        else()
            option(Aboria_RUN_TEST_${test_suite} "run ${test_suite} test suite with ctest" ON)
        endif()
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
    set(benchmark_test_suites
            SpeedTest
            BenchmarkFMM
            )
    foreach(test_suite ${benchmark_test_suites})
        option(Aboria_RUN_TEST_${test_suite} "run ${test_suite} test suite with ctest" ON)
        if (Aboria_RUN_TEST_${test_suite})
            set(test_suites ${test_suites} ${test_suite})
        endif()
    endforeach()
endif()

set(DocGettingStartedTestFile doc_getting_started.h)
set(DocGettingStartedTest
    test_getting_started
    )


set(BenchmarkFMMFile benchmark_fmm.h)
set(BenchmarkFMM
    test_multiquadric
    )

set(SpeedTestFile speed_test.h)
set(SpeedTest 
    test_vector_addition
    test_daxpy
    test_finite_difference
    test_multiquadric
    test_multiquadric_scaling
    test_linear_spring
    )

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

set(UtilsTestFile utils.h)
set(UtilsTest
    test_bucket_indicies
    test_point_to_bucket_indicies
    test_low_rank
    )

set(IteratorsTestFile iterators.h)
set(IteratorsTest
    test_lattice_within_distance
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
    test_dense_operator
    test_sparse_operator
    test_block_operator
    test_documentation
    )

set(ConstructorsTestFile constructors.h)

set(ConstructorsTest
    test_std_vector
    )
if (Aboria_USE_THRUST)
    list(APPEND ConstructorsTest
    test_thrust_vector
    )
endif()


set(NeighboursTestFile neighbours.h)
set(NeighboursTest
    test_std_vector_bucket_search_serial
    test_std_vector_bucket_search_serial_fast_bucketsearch
    test_std_vector_bucket_search_parallel
    test_std_vector_bucket_search_parallel_fast_bucketsearch
    test_std_vector_nanoflann_adaptor
    test_std_vector_octtree
    test_documentation
    )
if (Aboria_USE_THRUST)
    list(APPEND NeighboursTest
    test_thrust_vector_bucket_search_parallel
    test_thrust_vector_octtree
    )
endif()

set(IDSearchTestFile id_search.h)
set(IDSearchTest
    test_std_vector_bucket_search_serial
    test_std_vector_bucket_search_parallel
    test_std_vector_nanoflann_adaptor
    test_std_vector_octtree
    test_documentation
    )
if (Aboria_USE_THRUST)
    list(APPEND IDSearchTest
        #test_thrust_vector_bucket_search_serial
    test_thrust_vector_bucket_search_parallel
    test_thrust_vector_octtree
    )
endif()

set(ParticleContainerTestFile particle_container.h)
set(ParticleContainerTest 
    test_std_vector_bucket_search_serial
    test_std_vector_bucket_search_parallel
    test_documentation
    )
if (Aboria_USE_THRUST)
    list(APPEND ParticleContainerTest 
    test_thrust_vector_bucket_search_parallel
    )
endif()



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

set(MDLevel1TestFile md_level1.h)
set(MDLevel1Test "")
if (Aboria_USE_THRUST)
    list(APPEND MDLevel1Test
    test_std_vector_bucket_search_serial
    test_std_vector_bucket_search_parallel
    test_std_vector_octtree
    test_thrust_vector_bucket_search_parallel
    test_thrust_vector_octtree
    )
endif()


 
set(RbfInterpolationTestFile rbf_interpolation.h)
set(RbfInterpolationTest 
    test_bucket_search_serial
    test_bucket_search_parallel
    test_kdtree
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
    test_fmm_matrix_operators
    )
if (Aboria_USE_THRUST)
    list(APPEND H2Test
    test_fast_methods_bucket_search_serial_thrust
    test_fast_methods_bucket_search_parallel_thrust
    test_fast_methods_kd_tree_thrust
    test_fast_methods_octtree_thrust
    )
endif()


message("-- Adding test suites:")
foreach(test_suite ${test_suites})
    if (Aboria_THRUST_TEST_NVCC)
        #aboria_cxx_test(tests_cuda tests_cuda.cu ${test_files})
        aboria_cxx_test(${test_suite}_nvcc ${test_suite}.cu ON ${${test_suite}File} )
        target_link_libraries(${test_suite}_nvcc ${VTK_LIBRARIES} ${Boost_LIBRARIES} ${CUDA_LIBRARIES} ${Aboria_LIBRARIES})
    endif()

    if (Aboria_THRUST_TEST_CUDA_CLANG)
        aboria_cxx_test(${test_suite}_cuda_clang ${test_suite}.cpp OFF ${${test_suite}File})
        target_compile_options(${test_suite}_cuda_clang PRIVATE "--cuda-gpu-arch=sm_35" "-xcuda" "-pthread")
        target_link_libraries(${test_suite}_cuda_clang ${VTK_LIBRARIES} ${Boost_LIBRARIES} ${CUDA_LIBRARIES} ${Aboria_LIBRARIES} cudart_static dt rt)
    endif()


    aboria_cxx_test(${test_suite} ${test_suite}.cpp OFF ${${test_suite}File})
    target_link_libraries(${test_suite} ${VTK_LIBRARIES} ${Boost_LIBRARIES} ${CUDA_LIBRARIES} ${GPERFTOOLS_LIBRARIES} ${Aboria_LIBRARIES})
    
    foreach(test_name ${${test_suite}})
        add_test(NAME ${test_suite}:__${test_name} COMMAND ${RUNTIME_OUTPUT_DIRECTORY}/${test_suite} ${test_suite} ${test_name})
    endforeach()

    if (Aboria_USE_OPENMP)
        foreach(test_name ${${test_suite}})
            add_test(NAME ${test_suite}__OpenMP:__${test_name} 
                     COMMAND ${RUNTIME_OUTPUT_DIRECTORY}/${test_suite} ${test_suite} ${test_name})
            set_tests_properties(${test_suite}:__${test_name} 
                PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1")
            set_tests_properties(${test_suite}__OpenMP:__${test_name} 
                PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=3")
        endforeach()
    endif()
    foreach(test_name ${${test_suite}})
        if (Aboria_THRUST_TEST_NVCC)
            add_test(NAME ${test_suite}__NVCC:__${test_name} COMMAND ${RUNTIME_OUTPUT_DIRECTORY}/${test_suite}_nvcc ${test_suite} ${test_name})
        endif()
    endforeach()
endforeach()

