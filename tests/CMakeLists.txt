
set(test_package 
    constructors 
    particle_container
    neighbours
    geometry
    variables
    diffusion_around_spheres
    )

foreach(test ${test_package})
    cxx_test(${test} ${CMAKE_CURRENT_SOURCE_DIR}/${test}.h)
    target_link_libraries(${test} ${VTK_LIBRARIES} ${Boost_LIBRARIES})
endforeach(test)
