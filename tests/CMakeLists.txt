
set(test_package 
    variables
    particle_container
    constructors 
    geometry
    neighbours
    diffusion_around_spheres
    symbolic
    )

foreach(test ${test_package})
    cxx_test(${test} ${CMAKE_CURRENT_SOURCE_DIR}/${test}.h)
    #cxxtest_add_test(${test} ${test}.cpp ${CMAKE_CURRENT_SOURCE_DIR}/${test}.h)
    target_link_libraries(${test} ${VTK_LIBRARIES} ${Boost_LIBRARIES})
endforeach(test)
