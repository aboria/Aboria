file(GLOB Aboria_SOURCES "*.cpp")
file(GLOB Aboria_HEADERS "*.h")
include_directories(${Aboria_INCLUDE_DIRECTORIES})
message(STATUS ${Tyche_INCLUDE_DIRECTORIES})

add_library(Aboria SHARED ${Aboria_SOURCES})
TARGET_LINK_LIBRARIES(Aboria ${VTK_LIBRARIES} ${Boost_LIBRARIES} ${Boost_TIMER_LIBRARY})

install (TARGETS Aboria DESTINATION lib)
install (FILES ${Aboria_HEADERS} DESTINATION include)
install (DIRECTORY ${Aboria_SOURCE_DIR}/eigen3/ DESTINATION include)
