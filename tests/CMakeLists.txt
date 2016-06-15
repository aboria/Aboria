
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
    metafunctions
    )

#set(test_package
#    diffusion_around_spheres
#    metafunctions
#    )

foreach(test ${test_package})
    cxx_test(${test} ${CMAKE_CURRENT_SOURCE_DIR}/${test}.h)
    #cxxtest_add_test(${test} ${test}.cpp ${CMAKE_CURRENT_SOURCE_DIR}/${test}.h)
    target_link_libraries(${test} ${VTK_LIBRARIES} ${Boost_LIBRARIES})
endforeach(test)
