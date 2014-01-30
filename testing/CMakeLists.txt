
LINK_DIRECTORIES(${Aboria_BINARY_DIR}/lib)
include_directories(${Aboria_INCLUDE_DIRECTORIES})
add_executable (ParticleForEach particle_for_each.cpp)
target_link_libraries (ParticleForEach
                       Aboria
                       ${Boost_LIBRARIES}
                       )
