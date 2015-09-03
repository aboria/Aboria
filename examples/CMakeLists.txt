
set(examples
    bd.cpp 
    dem.cpp 
    sph.cpp 
    )

foreach(example ${examples})
    get_filename_component(example_name ${example} NAME_WE)
    add_executable(${example_name} ${example}) 
    target_link_libraries(${example_name} ${VTK_LIBRARIES} ${Boost_LIBRARIES})
endforeach(test)
