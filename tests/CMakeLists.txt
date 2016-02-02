
set(test_package 
    constructors 
    particle_container
    neighbours
    geometry
    variables
    diffusion_around_spheres
    symbolic
    )

foreach(test ${test_package})
    cxxtest_add_test(${test} ${test}.cpp ${CMAKE_CURRENT_SOURCE_DIR}/${test}.h)
    target_link_libraries(${test} ${VTK_LIBRARIES} ${Boost_LIBRARIES})
endforeach(test)
